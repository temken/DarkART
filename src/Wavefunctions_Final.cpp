#include "DarkARC/Wavefunctions_Final.hpp"

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Special_Functions.hpp"
#include "version.hpp"

namespace DarkARC
{
using namespace libphysica::natural_units;

// 1. Final state: Positive energy continuum solution of Schroedinger equation with hydrogenic potential
double Radial_Wavefunction_Final(double k_final, unsigned l_final, double Z_eff, double r)
{
	double eta = -Z_eff / k_final / Bohr_Radius;
	double rho = k_final * r;
	return 4.0 * M_PI / rho * Coulomb_Wave(l_final, eta, rho);
}

}	// namespace DarkARC