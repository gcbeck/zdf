/**
 * *****************************************************************************
 * \file constants.h
 * \author Graham Beck
 * \brief ZDF: Constants, largely for facilitating BLAS/LAPACK operations. 
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#pragma once

#include <filesystem>


namespace zdf
{
    static constexpr MKL_INT SINGLESTEP    = 1;
    static constexpr MKL_INT NOSTEP            = 0;
    static constexpr float ONEf                         = 1; 
    static constexpr float NEGATIVEONEf       = -1; 
    static constexpr float ZEROf                        = 0; 
    static constexpr char TRANSPOSED           = 'T'; 
    static constexpr char UNTRANSPOSED     = 'N'; 
    static constexpr char LEFTSIDE                   = 'L'; 
    static constexpr char RIGHTSIDE                = 'R'; 
   
    static constexpr char PATHSEP[2]               = {std::filesystem::path::preferred_separator, '\0'};
    static constexpr char SUFFIX[]                     = ".zdft"; 
    static constexpr char INPUT[]                      = ".zdfi"; 
    static constexpr char OUTPUT[]                  = ".zdfo"; 

    static constexpr float D2COEFS[]                = {6.0, 0.75, -3.5};

} // namespace zdf

// Satisfy the extern error messages in cx::err namespace
namespace cx
{
  namespace err
  {
    namespace
    {
      const char* abs_runtime_error = "Abs Runtime error";
      const char* sqrt_domain_error = "Sqrt Domain error";
      const char* exp_runtime_error = "Exp Runtime error";
      const char* floor_runtime_error= "Floor Runtime error";
      const char* ceil_runtime_error= "Ceil Runtime error";
      const char* fmod_domain_error= "Fmod Runtime error";
      const char* remainder_domain_error= "Remainder Domain error";
      const char* log_domain_error= "Log Domain error";
      const char* pow_runtime_error= "Pow Runtime error";
    }
  } // namespace err

} // namespace cx