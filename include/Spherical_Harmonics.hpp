#ifndef __Spherical_Harmonics_hpp_
#define __Spherical_Harmonics_hpp_

#include <complex>
#include <vector>

namespace DarkARC
{

extern std::complex<double> Spherical_Harmonics(int l, int m, double theta, double phi);

extern std::complex<double> VSH_Y_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Y(int l, int m, double theta, double phi);

extern std::complex<double> VSH_Psi_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Psi(int l, int m, double theta, double phi);

}	// namespace DarkARC

#endif