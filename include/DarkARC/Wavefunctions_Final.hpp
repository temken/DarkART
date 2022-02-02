#ifndef __Wavefunctions_Final_hpp_
#define __Wavefunctions_Final_hpp_

#include <string>
#include <vector>

namespace DarkARC
{

// 1. Final state: Positive energy continuum solution of Schroedinger equation with hydrogenic potential
extern double Radial_Wavefunction_Final(double k_final, unsigned l_prime, double Z_eff, double r);

}	// namespace DarkARC

#endif