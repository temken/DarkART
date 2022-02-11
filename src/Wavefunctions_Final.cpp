#include "DarkART/Wavefunctions_Final.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

#include "DarkART/Special_Functions.hpp"
#include "DarkART/version.hpp"

namespace DarkART
{
using namespace libphysica::natural_units;

// 1. Base class for final state wave function

Final_Electron_State* Final_Electron_State::Clone() const
{
	return new Final_Electron_State(*this);
}

// 2. Positive energy continuum solution of Schroedinger equation with hydrogenic potential
double Radial_Wavefunction_Hydrogenic(double k_final, unsigned l_final, double Z_eff, double r)
{
	double eta = -Z_eff / k_final / Bohr_Radius;
	double rho = k_final * r;
	return 4.0 * M_PI / rho * Coulomb_Wave(l_final, eta, rho);
}

Final_Electron_State_Hydrogenic::Final_Electron_State_Hydrogenic(double Z_eff)
: Z_effective(Z_eff)
{
}

double Final_Electron_State_Hydrogenic::Radial_Wavefunction(double r, double k_final, unsigned int l_final)
{
	return Radial_Wavefunction_Hydrogenic(k_final, l_final, Z_effective, r);
}

Final_Electron_State_Hydrogenic* Final_Electron_State_Hydrogenic::Clone() const
{
	return new Final_Electron_State_Hydrogenic(*this);
}

}	// namespace DarkART