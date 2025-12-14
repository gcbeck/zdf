/**
 * *****************************************************************************
 * \file util.h
 * \author Graham Beck
 * \brief ZDF: Useful & rather random utilities. 
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#pragma once

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdint>

constexpr uint32_t nCr(const unsigned short n, unsigned short r) {
    uint32_t x = 1;
    r = n - r > r ? n - r : r;
    for (unsigned short ix = 1; ix <= n - r; ++ix) { x *= r + ix; x /= ix; }
    return x;
}

template<size_t T>
constexpr unsigned short nchar() {
    unsigned short l = 0;
    for (auto n = T; n; l++, n /= 10);
    return l;
}

/**
* @brief Compile-time concatenation of multiple raw char arrays
*              into a std::array<char, .> 
*/
template<unsigned ...L>
constexpr auto join(const char (&...strings)[L]) {
    constexpr unsigned short N = (... + L) - sizeof...(L);
    std::array<char, N + 1> joined = {};
    joined[N] = '\0';

    auto it = joined.begin();
    (void)((it = std::copy_n(strings, L-1, it), 0), ...);
    return joined;
}

/**
* @brief Compile-time concatenation of a std::array<char, .> prefix, 
*               a integer T and a const char[.] postfix, useful in particular 
*               as a filepath for bppr files
*/
template<size_t T, size_t L, unsigned short U>
constexpr auto join(const std::array<char, L>& prefix, const char (&postfix)[U]) {
    constexpr unsigned short M = nchar<T>();
    constexpr unsigned short N = L + U + M - 1;
    std::array<char, N> joined = {};
    joined[N-1] = '\0';
    auto it = joined.begin();
    it = std::copy_n(prefix.begin(), L-1, it);
    it += M;
    for (auto n = T; n; n /= 10) {*--it = "0123456789"[(n % 10)]; }
    std::copy_n(postfix, U-1, it+M);
    return joined;
}

/**
* @brief Realigns circular/ring buffers that are col-major matrices
*              with U rows and L columns, given current head and tail indices. 
*
* @details Each segment entering or leaving the buffers has U rows and 
*                  runtime-dependent columns. The tail and head indices index 
*                  those columns. This function pulls all data back so that the 
*                  resulting tail index would be zero.  
*/
template<unsigned short U, unsigned short L>
void rectify(float (&circ)[U*L], const MKL_UINT& hx, const MKL_UINT& tx) {
    if (hx > tx) {
        std::memmove(circ, circ+tx*U, (hx > tx)*U*sizeof(float));
    } else if (2*tx >= L + hx) {
        std::memmove(circ+U*(L-tx), circ, hx*U*sizeof(float));
        std::memcpy(circ, circ+U*tx, U*(L-tx)*sizeof(float));
    } else {
        float* buffer = new float[hx*U];
        std::memcpy(buffer, circ, hx*U*sizeof(float));
        std::memmove(circ, circ+U*tx, U*(L-tx)*sizeof(float));
        std::memcpy(circ+U*(L-tx), buffer, hx*U*sizeof(float));
        delete[] buffer;
    }
}
