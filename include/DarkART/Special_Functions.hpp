#ifndef __DarkART_Special_Functions_hpp_
#define __DarkART_Special_Functions_hpp_

#include <complex>

namespace DarkART
{

// 1. Gaunt coefficients
extern double Gaunt_Coefficient(int j1, int j2, int j3, int m1, int m2, int m3);

// 2. Spherical Bessel function j_L
extern double Spherical_Bessel_jL(int L, double x);

// 3. Coulomb wave
extern double Coulomb_Wave(int L, double eta, double rho);

// 4. Hypergeometric functions
extern std::complex<double> Hypergeometric_2F1(double a, double b, double c, double z);

}	// namespace DarkART

#endif