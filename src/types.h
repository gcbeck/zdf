/**
 * *****************************************************************************
 * \file types.h
 * \author Graham Beck
 * \brief ZDF: Primarily sets up the en- & decoding of the ZDF type
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#pragma once

#include <cassert>
#include <cstddef>

namespace zdf
{
    using zdf_t = size_t;
    using zdfix_t = unsigned short;

    namespace zdfix
    {
        static const zdfix_t kN = 0x0000;
        static const zdf_t kNMask = 0x000000000000FFFF; // N has max 65535
        static const zdfix_t kNShift = 0;

        static const zdfix_t kD = 0x0001;
        static const zdf_t kDMask = 0x00000000FFFF0000; // D has max 65535
        static const zdfix_t kDShift = 16;

        static const zdfix_t kQ = 0x0002;
        static const zdf_t kQMask = 0x0000000F00000000; // Q has max 16
        static const zdfix_t kQShift = 32;

        static const zdfix_t kK = 0x0004;
        static const zdf_t kKMask = 0x000000F000000000; // K has max 16
        static const zdfix_t kKShift = 36;

        static const zdfix_t kM = 0x0008;
        static const zdf_t kMMask = 0x00000F0000000000; // M has max 16
        static const zdfix_t kMShift = 40;

        static const zdfix_t kU = 0x0010;
        static const zdf_t kUMask = 0x0000F00000000000; // U has max 16
        static const zdfix_t kUShift = 44;

        // static const zdfix_t kNext = 0x0020;
        // static const zdfix_t kNextShift = 48;

        static const zdf_t kInvalid = 0;

        template <zdfix_t T>
        static constexpr zdf_t decode(const zdf_t& type) { return kInvalid; }

        template <>
        constexpr zdf_t decode<kN>(const zdf_t& type) { return (type & kNMask) >> kNShift; }
        template <>
        constexpr zdf_t decode<kD>(const zdf_t& type) { return (type & kDMask) >> kDShift; }
        template <>
        constexpr zdf_t decode<kQ>(const zdf_t& type) { return (type & kQMask) >> kQShift; }
        template <>
        constexpr zdf_t decode<kK>(const zdf_t& type) { return (type & kKMask) >> kKShift; }
        template <>
        constexpr zdf_t decode<kM>(const zdf_t& type) { return (type & kMMask) >> kMShift; }
        template <>
        constexpr zdf_t decode<kU>(const zdf_t& type) { return (type & kUMask) >> kUShift; }

        template <zdfix_t T>
        static constexpr zdf_t encode(const zdf_t& type) { return kInvalid; }

        template <>
        constexpr zdf_t encode<kN>(const zdf_t& type)
        {
            assert(type <= (kNMask >> kNShift));
            return type << kNShift;
        }
        template <>
        constexpr zdf_t encode<kD>(const zdf_t& type)
        {
            assert(type <= (kDMask >> kDShift));
            return type << kDShift;
        }
        template<zdfix_t T, unsigned short nD, std::size_t... iD>
        constexpr zdf_t encode(unsigned short (&&derivs)[nD], std::index_sequence<iD...>)
        {
            return encode<kD>((... | (1 << derivs[iD])));
        }
        template <>
        constexpr zdf_t encode<kQ>(const zdf_t& type)
        {
            assert(type <= (kQMask >> kQShift));
            return type << kQShift;
        }
        template <>
        constexpr zdf_t encode<kK>(const zdf_t& type)
        {
            assert(type <= (kKMask >> kKShift));
            return type << kKShift;
        }
        template <>
        constexpr zdf_t encode<kM>(const zdf_t& type)
        {
            assert(type <= (kMMask >> kMShift));
            return type << kMShift;
        }
        template <>
        constexpr zdf_t encode<kU>(const zdf_t& type)
        {
            assert(type <= (kUMask >> kUShift));
            return type << kUShift;
        }

        static constexpr zdf_t encode(unsigned short&& N, unsigned short &&D, unsigned short&& Q, unsigned short&& K, unsigned short&& M, unsigned short&& U)
        {
            return encode<kN>(N) | encode<kD>(D) | encode<kQ>(Q) | encode<kK>(K) | encode<kM>(M) | encode<kU>(U);
        }

        template <unsigned short nD>
        static constexpr zdf_t encode(unsigned short&& N, unsigned short (&&derivs)[nD], unsigned short&& Q, unsigned short&& K, unsigned short&& M, unsigned short&& U)
        {
            return encode<kN>(N) | encode<kD>(std::move(derivs), std::make_index_sequence<nD>{}) | encode<kQ>(Q) | encode<kK>(K) | encode<kM>(M) | encode<kU>(U);
        }

    } // namespace zdfix

} // namespace zdf