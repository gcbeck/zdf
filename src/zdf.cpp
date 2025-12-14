/**
 * *****************************************************************************
 * \file zdf.cpp
 * \author Graham Beck
 * \brief ZDF: Provides filter information; constructs and convolves a filter over input data. 
 * \version 0.1
 * \date 2025-12-01
 *
 * \copyright Copyright (c) 2025
 * *****************************************************************************
 */
#include "zdf.h"

#include <concepts>
#include <iostream>
#include <sstream> 
#include <string>
#include <unistd.h>

#include "mkl_service.h"
#include "mkl_vml_defines.h"
#include "mkl_vml_functions.h"

constexpr auto REPO = join(".", zdf::PATHSEP, "dat", zdf::PATHSEP);
constexpr auto T = zdf::zdfix::encode(512, {0,1,2}, 2, 0, 1, 2);

constexpr char OPTS[] = "d:e:in:t:w";
constexpr char OPTSEP = ',';

template <typename U>
concept Numeric = std::is_arithmetic_v<U>;

template<Numeric U>
U next(std::stringstream& ss) {
    std::string token;
    if (std::getline(ss, token, OPTSEP)) { 
        if constexpr (std::is_integral<U>::value) { return std::stoi(token); } else { return std::stof(token); }
    } else { 
        throw std::runtime_error("Unparseable Parameter"); 
    }
}

int main(int argc, char *argv[])
{
    vmlSetMode(VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_DEFAULT);  

     int opt;
    bool write = false;

    while ((opt = getopt(argc, argv, OPTS)) != -1) {
        switch (opt) {
          case 'd': {
            std::stringstream ss(optarg);
            std::string token; unsigned short denc = 0;
            while (std::getline(ss, token, OPTSEP)) { 
                denc |= 1 << std::stoi(token);
            }
            std::cout << "Encoded Derivatives: " << denc << std::endl;
            return 0;
          }
          case 'e': {
            std::stringstream ss(optarg);
            std::cout << "Encoding: " << zdf::zdfix::encode(next<unsigned short>(ss), 
                                                                                               next<unsigned short>(ss), 
                                                                                               next<unsigned short>(ss), 
                                                                                               next<unsigned short>(ss), 
                                                                                               next<unsigned short>(ss), 
                                                                                               next<unsigned short>(ss)) << std::endl;
            return 0;
          }
          case 'i': {
            std::cout << "Minimal-Filter <" << T << "> Delay: " << zdf::ZDF<T>::kDelay << std::endl;
            return 0;
          }
          case 'n': {
            std::stringstream ss(optarg);
            std::cout << "Minimum Second-Derivative-Based Filter Length: " << zdf::ZDF<T>::d2N(next<float>(ss), next<float>(ss)) << std::endl;
            return 0;
          }
          case 't':
            mkl_set_num_threads(std::stoi(optarg));
            break;
          case 'w':
            write = true;
            break;
        }
    }

    zdf::ZDF<T> zdf(REPO);

    zdf.update<256>(REPO);

    if (write) { zdf.write(REPO); }

    return 0;
}