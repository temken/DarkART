#include "DarkART/Wavefunctions_Final.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkART/Schroedinger_Equation.hpp"
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
: initial_state(ini_state), r_min(1e-4 * Bohr_Radius), r_max(51.0 * Bohr_Radius)
{
	std::vector<double> r_list = libphysica::Log_Space(r_min, r_max, 5);
	std::vector<double> Z_eff_list(r_list.size(), Z_eff);
	Z_effective_interpolation = libphysica::Interpolation(r_list, Z_eff_list);
}

double Final_Electron_State_Schroedinger::Z_effective(double r)
{
	if(r <= r_min)
		return initial_state.Z;
	else if(r >= r_max)
		return 1.0;
	else
		return Z_effective_interpolation(r);
}

std::vector<Initial_Electron_State> Final_Electron_State_Schroedinger::Import_All_States(const std::string& element)
{
	// Import all wave functions of the element.
	std::vector<Initial_Electron_State> all_states = {};
	std::vector<std::string> l_orbital_names	   = {"s", "p", "d", "f", "g"};
	for(int n = 1; n < 7; n++)
	{
		bool shell_data_found = false;
		for(int l = 0; l < n; l++)
		{
			std::string path = TOP_LEVEL_DIR "data/" + initial_state.element_name + "_" + std::to_string(n) + l_orbital_names[l] + ".txt";
			if(libphysica::File_Exists(path))
			{
				shell_data_found = true;
				all_states.push_back(Initial_Electron_State(initial_state.element_name, n, l));
			}
			else
				break;
		}
		if(!shell_data_found)
			break;
	}
	return all_states;
}

void Final_Electron_State_Schroedinger::Determine_Z_effective()
{
	// 1. Tabulate Z_effective
	auto all_states				   = Import_All_States(initial_state.element_name);
	auto r_list					   = libphysica::Log_Space(r_min, r_max, 50);
	std::vector<double> Z_eff_list = {};
	for(auto& r : r_list)
	{
		double Z_eff = initial_state.Z;
		for(auto& state : all_states)
			Z_eff -= ((state.n == initial_state.n && state.l == initial_state.l) ? (4.0 * state.l + 1.0) : 2.0 * (2.0 * state.l + 1)) * state.Radial_Integral(r);
		// std::cout << r / Bohr_Radius << "\t" << Z_eff << std::endl;
		Z_eff_list.push_back(Z_eff);
	}
	// 2. Interpolate Z_effective
	Z_effective_interpolation = libphysica::Interpolation(r_list, Z_eff_list);
}

void Final_Electron_State_Schroedinger::Determine_Z_effective_2()
{
	// 1. Tabulate Z_effective
	auto r_list					   = libphysica::Log_Space(r_min, r_max, 50);
	std::vector<double> Z_eff_list = {};
	bool Z_eff_is_one			   = false;
	for(auto r : r_list)
	{
		int l		  = initial_state.l;
		double EB	  = initial_state.binding_energy;
		double R	  = initial_state.Radial_Wavefunction(r);
		double dRdr	  = initial_state.Radial_Wavefunction_Derivative(r);
		double d2Rdr2 = initial_state.Radial_Wavefunction_Derivative_2(r);

		r  = r / Bohr_Radius;
		EB = -EB / Rydberg;
		R *= std::pow(Bohr_Radius, 1.5);
		dRdr *= std::pow(Bohr_Radius, 2.5);
		d2Rdr2 *= std::pow(Bohr_Radius, 3.5);
		double Z_eff = Z_eff_is_one ? 1.0 : 1.0 + (-r / 2.0 * (-EB - l * (l + 1.0) / r / r + (2.0 * dRdr + r * d2Rdr2) / r / R) - 1.0) * (libphysica::StepFunction(5.0 - r) + libphysica::StepFunction(r - 5.0) * std::exp(-(r - 5.0) * (r - 5.0)));
		Z_eff_list.push_back(Z_eff);
		if(Z_eff_is_one == false && libphysica::Relative_Difference(Z_eff, 1.0) < 1e-6)
			Z_eff_is_one = true;
	}
	// 2. Interpolate Z_effective
	Z_effective_interpolation = libphysica::Interpolation(r_list, Z_eff_list);
}

void Final_Electron_State_Schroedinger::Solve_Schroedinger_Equation(double k_final, unsigned int l_final)
{
	Schroedinger_Solver schroedinger_solution(k_final, l_final);
	radial_wavefunction = schroedinger_solution.Solve_Radial_Equation();
}

double Final_Electron_State_Schroedinger::Radial_Wavefunction(double r, double k_final, unsigned int l_final)
{
	return radial_wavefunction(r);
}

Final_Electron_State_Schroedinger* Final_Electron_State_Schroedinger::Clone() const
{
	return new Final_Electron_State_Schroedinger(*this);
}
}	// namespace DarkART