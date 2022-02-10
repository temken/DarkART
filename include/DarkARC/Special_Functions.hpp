#ifndef __Special_Functions_hpp_
#define __Special_Functions_hpp_

#include <complex>
#include <vector>

namespace DarkARC
{

// 1. Spherical Harmonics
extern std::complex<double> Spherical_Harmonics(int l, int m, double theta, double phi);

extern std::complex<double> VSH_Y_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Y(int l, int m, double theta, double phi);

extern std::complex<double> VSH_Psi_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Psi(int l, int m, double theta, double phi);

extern double Gaunt_Coefficient(int j1, int j2, int j3, int m1, int m2, int m3);

// 2. Spherical Bessel function j_L
extern double Spherical_Bessel_jL(int L, double x);

// 3. Coulomb wave
extern double Coulomb_Wave(int L, double eta, double rho);

// 4. Legendre polynomial roots and weights for Gauss-Legendre Quadrature
extern std::vector<std::vector<double>> Compute_Gauss_Legendre_Roots_and_Weights(unsigned int n, double x_min = -1.0, double x_max = 1.0);

unsigned int Locate_Closest_Location(const std::vector<double>& list, double target);

}	// namespace DarkARC

#endif