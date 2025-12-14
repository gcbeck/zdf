/**
 * *****************************************************************************
 * \file zdf.h
 * \author Graham Beck
 * \brief ZDF: Contructs the ZDF filter coefficients for derivatives of any order and any 
 *                     'timescale' mu. Compensates for minimal-filter delay by adding specific 
 *                      (mu,kappa) filter combinations together. 
 *                      Implementation of research by Mboup et al.
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#pragma once

#include <bit>
#include <climits>
#include <utility>

#include "cx_math.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mkl_vml.h"

#include "constants.h"
#include "proto.h"
#include "types.h"


namespace zdf
{
    using wkx_t = unsigned short;
    namespace wkx {
        static const wkx_t kOnex   = 0; 
        static const wkx_t kTaux    = 1; 
        static const wkx_t kTauKx  = 2; 
        static const wkx_t kTauMx = 3;

        static const wkx_t N         = kTauMx + 1;
    } // namespace wkx

    template <zdf_t T>
    class ZDF
    {
      public:
        static constexpr MKL_INT N = zdfix::decode<zdfix::kN>(T);                                                                    // Length of filter
        static constexpr unsigned short nD = std::popcount(zdfix::decode<zdfix::kD>(T));                           // Number of derivatives
        static constexpr unsigned short nM = zdfix::decode<zdfix::kU>(T) - zdfix::decode<zdfix::kM>(T); // Number of mu timescales
        static constexpr unsigned short d0 = std::countr_zero(zdfix::decode<zdfix::kD>(T));                       // Smallest-order derivative
        static constexpr MKL_INT kF = nD*nM;
        
        /**
        * @brief The delay for a (minimal) filter of derivative order n and parameters kappa & mu. It is this delay that is 
        *              compensated by introducing non-minimality (q>0) at the cost of filter performance
        */
        static constexpr float delay(const unsigned short n, const unsigned short kappa, const unsigned short mu) {
            return static_cast<float>(kappa + n + 1) / (mu + kappa + 2*(n+1));
        }
        /**
        * @brief A recommended minimum length for a filter that will return second derivatives given a signal-to-noise ratio
        *              and a required 'similarity'
        * 
        * @details The signal-to-noise ratio may be estimated as the RMS of the (notional) signal divided by the
        *                  standard deviation of the noise
        *                  The similarity is a measure of how well the phase of the calculated second derivative matches its theoretical form.
        *                  As the filter length decreases a phase shift is induced which in turn decreases similarity. Increased noise also 
        *                  hurts similarity, requiring a longer filter length. 
        */
        static constexpr unsigned short d2N(const float& snr, const float& similarity) {
            const float ilsnr = 1/cx::log(snr);
            return static_cast<unsigned short>(cx::pow(2.0f, D2COEFS[0] + D2COEFS[1]*ilsnr + D2COEFS[2]*cx::log(1-similarity)*ilsnr));
        }
        static constexpr unsigned short kDelay = static_cast<unsigned short>(N * delay(d0, zdfix::decode<zdfix::kK>(T), zdfix::decode<zdfix::kM>(T)));                                               
       /**
        * @brief Contructs a Zero Delay Filter from a binary file containing the N signal values used for initialization
        */
        template<size_t M>
        ZDF(const std::array<char, M>& from)
            : _proto(from)
            , _hx(0)
        {
            _proto.template get<proto::kFCore>(_X);
            _proto.close();

            std::memset(_fir, 0, N*kF*sizeof(float));

            float working[N*wkx::N];
            std::fill(working+N*wkx::kOnex, working+N*(wkx::kOnex+1), 1);
            for (unsigned short ix = 0; ix < N; ++ix ) {working[N*wkx::kTaux+ix] = N-ix-1; }
            const float n = static_cast<float>(N);
            vsDivI(N, working+N*wkx::kTaux, SINGLESTEP, &n, NOSTEP, working+N*wkx::kTaux, SINGLESTEP);

            constexpr unsigned short q = zdfix::decode<zdfix::kQ>(T);
            float lmbd[q+1];

            float derivs[nD]; derivs[0] = d0;
            const unsigned short de = sizeof(zdf_t)*CHAR_BIT - std::countl_zero(zdfix::decode<zdfix::kD>(T));
            for (unsigned short dx = d0+1, ix = 1; dx < de; ++dx) {
                if ((zdfix::decode<zdfix::kD>(T) >> dx) & 1) {
                    derivs[ix++] = dx;
                }
            }

            for (unsigned short nj = 0; nj < nD; ++nj) {
                const unsigned short nx = derivs[nj];
                for (unsigned short mj = 0; mj < nM; ++mj) {
                    const unsigned short mx = zdfix::decode<zdfix::kM>(T) + mj;
                    const unsigned short fx = mj+nM*nj;

                    std::memset(lmbd, 0, q*sizeof(float)); lmbd[q] = 1.0f;
                    float* const firx = _fir+N*fx;
                    // Get the minimal filter
                    h(nx, zdfix::decode<zdfix::kK>(T), mx, lmbd[q], working, reinterpret_cast<float(&)[N]>(*firx));
                    // Get the normalizing term for the minimal filter
                    if (nx > 0) { 
                        _filtered[fx] = sdot(&N, firx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP) / N;
                        vsSubI(N, firx, SINGLESTEP, _filtered+fx, NOSTEP, firx, SINGLESTEP);
                        vsAbs(N, firx, working+N*wkx::kTauMx);
                        _filtered[fx] = sdot(&N, working+N*wkx::kTauMx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP);
                    } else {
                        _filtered[fx] = sdot(&N, firx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP);
                    }

                    if constexpr (q > 0) {
                        // Solve for the Q-hull coefficients lambda
                        constexpr MKL_INT nq = q+1;
                        float gram[nq*nq];
                        const float c = 1.0f / nCr(mx+zdfix::decode<zdfix::kK>(T)+2*(nx+q)+1, q);
                        for (unsigned short ix = 0; ix <= q; ++ix) {
                            for (unsigned short jx = 0; jx <= q; ++jx) {
                                gram[ix+nq*jx] = c*nCr(mx+nx+ix+jx, mx+nx+jx);
                                gram[ix+nq*jx] *= nCr(zdfix::decode<zdfix::kK>(T)+nx+2*q-ix-jx, zdfix::decode<zdfix::kK>(T)+nx+q-jx);
                            }
                        }
                        MKL_INT ipiv[nq]; MKL_INT outcome;
                        sgesv(&nq, &SINGLESTEP, gram, &nq, ipiv, lmbd, &nq, &outcome);

                        // Add together the non-minimal components
                        std::memset(firx, 0, N*sizeof(float));
                        for (unsigned short ix = 0; ix <= q; ++ix) {
                            h(nx, zdfix::decode<zdfix::kK>(T)+q-ix, mx+ix, lmbd[ix], working, reinterpret_cast<float(&)[N]>(*firx));
                        }
                        // Calculate the normalization
                        if (nx > 0) { 
                            const float Z = sdot(&N, firx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP) / N;
                            vsSubI(N, firx, SINGLESTEP, &Z, NOSTEP, firx, SINGLESTEP);
                            vsAbs(N, firx, working+N*wkx::kTauMx);
                            _filtered[fx] /= sdot(&N, working+N*wkx::kTauMx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP);
                        } else {
                            _filtered[fx] /= sdot(&N, firx, &SINGLESTEP, working+N*wkx::kOnex, &SINGLESTEP);
                        }
                    } else {
                        _filtered[fx] = 1/_filtered[fx];
                    }
                    vsMulI(N, firx, SINGLESTEP, _filtered+fx, NOSTEP, firx, SINGLESTEP);
                }
            }
            // Apply the filter to the initialization data
            sgemv(&TRANSPOSED, &N, &kF, &ONEf, _fir, &N, _X, &SINGLESTEP, &ZEROf, _filtered, &SINGLESTEP);
        }

        /**
        * @brief Perform the filtering for a new signal value, returning all nD derivatives at their nM mu timescales
        */
        const float(&update(const float& x))[kF] {
            _X[_hx++] = x;
            sgemv(&TRANSPOSED, &_hx, &kF, &ONEf, _fir+N-_hx, &N, _X, &SINGLESTEP, &ZEROf, _filtered, &SINGLESTEP);
            if (_hx < N) {
                const MKL_INT tx = N-_hx;
                sgemv(&TRANSPOSED, &tx, &kF, &ONEf, _fir, &N, _X+_hx, &SINGLESTEP, &ONEf, _filtered, &SINGLESTEP);
            }
            _hx %= N;
            return _filtered;
        }
        /**
        * @brief Convenience/demonstration function that performs all online updates from the rows in file 'from'. 
        */
        template<size_t P, size_t M>
        void update(const std::array<char, M>& from) { 
            _proto.template open<M, proto::kFIn>(from);
            float x; float cache[P*kF]; float* cachex = cache;
            while(_proto.template get<proto::kFIn>(x)) {
                std::memcpy(cachex, update(x), kF*sizeof(float));
                cachex += kF;
            }
            _proto.close();
            const size_t nUpdates = (cachex - cache) / kF;
            _proto.template open<M, proto::kFOut>(from);
            for (size_t ix = 0; ix < nUpdates; ++ix) {
                _proto.template set<proto::kFOut>(reinterpret_cast<float(&)[kF]>(cache[ix*kF]));
            }
            _proto.close();
        }

        /**
        * @brief Persist the updated observable series to file. 
        */
        template<size_t M>
        void write(const std::array<char, M>& to) {
            if (_hx > 0) { rectify<1, N>(_X, _hx, _hx); _hx = 0; }
            _proto.open(to);
            _proto.template set<proto::kFCore>(_X);
            _proto.close();
        }

      private:
        // Generate the minimal filter
        static void h(const unsigned short n, const unsigned short kappa, const unsigned short mu, const float& lambda, float (&wkg)[N*wkx::N], float (&out)[N]) {
            const auto gamma = (lambda*((mu+n+1)*nCr(kappa+mu+2*n+1, kappa+n))) / N;
            float w =1;
            for (unsigned short ix = 1; ix <= n; ++ix) { w *= (kappa + ix); }
            for (unsigned short ix = 0; ix <= n; ++ix) {
                vsPowx(N, wkg+N*wkx::kTaux, kappa+ix, wkg+N*wkx::kTauKx);
                vsSubI(N, &ONEf, NOSTEP, wkg+N*wkx::kTaux, SINGLESTEP, wkg+N*wkx::kTauMx, SINGLESTEP);
                vsPowx(N, wkg+N*wkx::kTauMx, mu+n-ix, wkg+N*wkx::kTauMx);
                vsMul(N, wkg+N*wkx::kTauKx, wkg+N*wkx::kTauMx, wkg+N*wkx::kTauKx);
                const float c = (ix%2 ? -w : w) * nCr(n, ix) * gamma;
                saxpy(&N, &c, wkg+N*wkx::kTauKx, &SINGLESTEP, out, &SINGLESTEP);
                w *= static_cast<float>(mu+n-ix) / (kappa + ix+1); 
            }
        }

        Proto<T> _proto;
        float _X[N];
        float _fir[N*kF];
        float _filtered[kF];
        MKL_INT _hx;
    };

} // namespace zdf