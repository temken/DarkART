#include "DarkARC/Wavefunctions_Final.hpp"

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Special_Functions.hpp"
#include "version.hpp"

namespace DarkARC
{
using namespace libphysica::natural_units;

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

void Final_Electron_State_Hydrogenic::Fit_Zeff(int n, double binding_energy)
{
	double au	= 27.211386245988 * eV;
	Z_effective = sqrt(-2.0 * binding_energy / au) * n;
}

double Final_Electron_State_Hydrogenic::Radial_Wavefunction(double r, double k_final, unsigned int l_final)
{
	return Radial_Wavefunction_Hydrogenic(k_final, l_final, Z_effective, r);
}

}	// namespace DarkARC