/**
 * *****************************************************************************
 * \file proto.h
 * \author Graham Beck
 * \brief ZDF: Protocol for reading/writing the observable time series and its derivatives 
 *                     to binary file. 
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#pragma once

#include <fstream>
#include <stdexcept>

#include "constants.h"
#include "types.h"
#include "util.h"


namespace zdf
{
    using protix_t = unsigned short; // Enumeration of the file type
    namespace proto {
        static constexpr protix_t kFInvalid  = 0x0000;
        static constexpr protix_t kFCore      = 0x0001;
        static constexpr protix_t kFIn           = 0x0002;
        static constexpr protix_t kFOut        = 0x0004;
    } // namespace proto

    template <zdf_t T>
    class Proto{
      public:
      /**
        * @brief The Proto class establishes the persistence protocol, reading and writing to file
        *              for initialization / warmstarting of state. The constructor takes as argument the 
        *              directory where zdf files are persisted. 
        * 
        * @note The expected .zdf file is opened for reading by construction and must be closed 
        *               elsewhere once initialization has been completed. 
        */
        template<size_t M>
        explicit Proto(const std::array<char, M>& from)
            : _from(join<T>(from, suffix<proto::kFCore>()).data(), std::ios::in | std::ios::binary)
        {
            if (!_from.is_open()) {
                throw std::runtime_error(join<T>(from, suffix<proto::kFCore>()).data());
            }
        }

        void close() { _from.close(); }

      /**
        * @brief Opens a binary file for reading or writing, depending on file purpose
        * 
        * @details When open(.) is called explicitly on the .zdft file (signaled by proto::kFCore)
        *                  the file mode is set to write for persistence of the most recent observable series. 
        *                  On the other hand when called on the input data file (signaled by 
        *                  proto::kFIn) then reading is assumed. 
        */
        template<size_t M, protix_t U=proto::kFCore>
        void open(const std::array<char, M>& to) { 
            _from.open(join<T>(to, suffix<U>()).data(), mode<U>());
            if (!_from.is_open()) {
                throw std::runtime_error(join<T>(to, suffix<U>()).data());
            }
        }

      /**
        * @brief Reads template-specific data into the arguments submitted
        */
        template<protix_t U, typename... Ts>
        bool get(Ts&... args) {
            return io<U>::get(*this, args...);
        }

      /**
        * @brief Writes the protocol-observant data to the file; otherwise similar to get(...)
        */
        template<protix_t U, typename... Ts>
        void set(Ts&... args) {
            return io<U>::set(*this, args...);
        }

      private:
        template<protix_t U>
        const char(&suffix())[std::size(zdf::SUFFIX)] { return zdf::SUFFIX; }
        template<> const char(&suffix<proto::kFIn>())[std::size(zdf::INPUT)] { return zdf::INPUT; }
        template<> const char(&suffix<proto::kFOut>())[std::size(zdf::OUTPUT)] { return zdf::OUTPUT; }

        template<protix_t U>
        const auto mode() { return std::ios::out | std::ios::binary; }
        template<> const auto mode<proto::kFIn>() { return std::ios::in | std::ios::binary; }

        template<protix_t U> struct io;

        template<> struct io<proto::kFCore> {
            static bool get(Proto<T>& p, float(&X)[zdfix::decode<zdfix::kN>(T)]) {
                p._from.read(reinterpret_cast<char*>(X), sizeof(X));
                return true;
            }
            static void set(Proto<T>& p, const float(&X)[zdfix::decode<zdfix::kN>(T)]) {
                p._from.write(reinterpret_cast<const char*>(X), sizeof(X));
            }
        };

        template<> struct io<proto::kFIn> {
            static bool get(Proto<T>& p, float& X) {
                p._from.read(reinterpret_cast<char*>(&X), sizeof(float));
                if (p._from) { return true; }
                return false;
            }
        };

        template<> struct io<proto::kFOut> {
            static void set(Proto<T>& p, const float(&X)[std::popcount(zdfix::decode<zdfix::kD>(T))*(zdfix::decode<zdfix::kU>(T) - zdfix::decode<zdfix::kM>(T))]) {
                p._from.write(reinterpret_cast<const char*>(X), sizeof(X));
            }
        };

        std::fstream _from;
    };

} // namespace zdf