#include "DarkART/Wavefunctions_Final.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

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

// 2. Positive energy continuum solution of Schroedinger equation with hydrogenic potential (i.e. constant Z_eff)
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

// 3. Positive energy continuum solution of Schroedinger equation for a given potential Z_eff(r)

Final_Electron_State_Schroedinger::Final_Electron_State_Schroedinger(Initial_Electron_State& ini_state, double Z_eff)
: initial_state(ini_state), r_min(0.0), r_max(50.0 * Bohr_Radius)
{
	std::vector<double> r_list = libphysica::Linear_Space(r_min, r_max, 5);
	std::vector<double> Z_eff_list(r_list.size(), Z_eff);
	Z_effective_interpolation = libphysica::Interpolation(r_list, Z_eff_list);
}

double Final_Electron_State_Schroedinger::Z_effective(double r)
{
	return Z_effective_interpolation(r);
}

void Final_Electron_State_Schroedinger::Determine_Z_effective()
{
}

double Final_Electron_State_Schroedinger::Radial_Wavefunction(double r, double k_final, unsigned int l_final)
{
	return 0.0;
}

Final_Electron_State_Schroedinger* Final_Electron_State_Schroedinger::Clone() const
{
	return new Final_Electron_State_Schroedinger(*this);
}
}	// namespace DarkART