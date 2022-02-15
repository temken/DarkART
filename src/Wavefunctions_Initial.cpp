#include "DarkART/Wavefunctions_Initial.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

#include <boost/math/special_functions/factorials.hpp>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkART/Special_Functions.hpp"
#include "DarkART/version.hpp"

namespace DarkART
{
using namespace std::complex_literals;
using namespace libphysica::natural_units;
using boost::math::factorial;

double a0 = Bohr_Radius;
double au = 27.211386245988 * eV;

// 1. Initial state: Roothaan-Hartree-Fock Ground-State Atomic Wave Functions

void Initial_Electron_State::Import_RHF_Coefficients()
{
	std::string filepath = TOP_LEVEL_DIR "data/" + Orbital_Name() + ".txt";
	if(libphysica::File_Exists(filepath) == false)
	{
		std::cerr << "Error in Initial_Electron_State::Import_RHF_Coefficients(): Coefficient table for " << Orbital_Name() << " does not exist." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	std::ifstream f;
	f.open(filepath);
	if(f.is_open())
	{
		f >> binding_energy;
		binding_energy *= au;
		double C, Z;
		unsigned int nin;
		while(f >> nin >> Z >> C)
		{
			n_lj.push_back(nin);
			Z_lj.push_back(Z);
			C_nlj.push_back(C);
		}
		f.close();
	}
	Z_eff = sqrt(-2.0 * binding_energy / au) * n;
}

void Initial_Electron_State::Check_Normalization()
{
	double norm = Normalization();
	if(std::fabs(1.0 - norm) > 1.0e-4)
	{
		std::cout << "Error in Initial_Electron_State(): Normalization of " << element_name << " = " << norm << " != 1.0" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

std::vector<std::string> element_names = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
Initial_Electron_State::Initial_Electron_State()
: element_name("none"), n(0), l(0), binding_energy(0.0), Z_eff(0.0)
{
}

Initial_Electron_State::Initial_Electron_State(int z, int N, int L)
: element_name(element_names[z - 1]), Z(z), n(N), l(L)
{
	Import_RHF_Coefficients();
	Check_Normalization();
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, int N, int L)
: element_name(element), n(N), l(L)
{
	auto const it = std::find(element_names.begin(), element_names.end(), element);
	Z			  = std::distance(element_names.begin(), it) + 1;
	Import_RHF_Coefficients();
	Check_Normalization();
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, std::string shell_name)
: element_name(element)
{
	auto const it = std::find(element_names.begin(), element_names.end(), element);
	Z			  = std::distance(element_names.begin(), it) + 1;
	n			  = shell_name[0] - '0';
	for(l = 0; l < l_orbital_names.size(); l++)
		if(shell_name[1] == l_orbital_names[l][0])
			break;
	Import_RHF_Coefficients();
	Check_Normalization();
}

std::string Initial_Electron_State::Orbital_Name() const
{
	return element_name + "_" + std::to_string(n) + l_orbital_names[l];
}

double Initial_Electron_State::Radial_Wavefunction(double r) const
{
	double R_nl = 0.0;
	for(unsigned int j = 0; j < C_nlj.size(); j++)
		R_nl += C_nlj[j] * std::pow(2.0 * Z_lj[j], n_lj[j] + 0.5) / sqrt(factorial<double>(2.0 * n_lj[j])) * std::pow(r / a0, n_lj[j] - 1.0) * std::exp(-Z_lj[j] * r / a0);

	return std::pow(a0, -1.5) * R_nl;
}

double Initial_Electron_State::Radial_Wavefunction_Derivative(double r) const
{
	double dR_dr = 0.0;
	for(unsigned int j = 0; j < C_nlj.size(); j++)
		dR_dr += C_nlj[j] * std::pow(2.0 * Z_lj[j], n_lj[j] + 0.5) / sqrt(factorial<double>(2.0 * n_lj[j])) * ((n_lj[j] - 1.0) / a0 * std::pow(r / a0, n_lj[j] - 2.0) - Z_lj[j] / a0 * std::pow(r / a0, n_lj[j] - 1.0)) * std::exp(-Z_lj[j] * r / a0);

	return std::pow(a0, -1.5) * dR_dr;
}

double Initial_Electron_State::Radial_Wavefunction_Derivative_2(double r) const
{
	double dR2_dr2 = 0.0;
	for(unsigned int j = 0; j < C_nlj.size(); j++)
		dR2_dr2 += C_nlj[j] * std::pow(2.0 * Z_lj[j], n_lj[j] + 0.5) / sqrt(factorial<double>(2.0 * n_lj[j])) * std::exp(-Z_lj[j] * r / a0) * (std::pow(Z_lj[j] / a0, 2.0) * std::pow(r / a0, n_lj[j] - 1.0) - 2.0 * Z_lj[j] * (n_lj[j] - 1.0) / a0 / a0 * std::pow(r / a0, n_lj[j] - 2.0) + (n_lj[j] - 1.0) * (n_lj[j] - 2.0) / a0 / a0 * std::pow(r / a0, n_lj[j] - 3.0));

	return std::pow(a0, -1.5) * dR2_dr2;
}

double Initial_Electron_State::Normalization() const
{
<<<<<<< HEAD
	std::function<double(double)> integrand = [this](double r) {
		double R = Radial_Wavefunction(r);
		return r * r * R * R;
	};
	// Integrate with Gauss Legendre
	return libphysica::Integrate_Gauss_Legendre(integrand, 0.0, 50.0 * Bohr_Radius, 1000);
=======
	double r_max = 50.0 * Bohr_Radius;
	return Radial_Integral(r_max);
>>>>>>> c55d261 (Small extension to initial wavefunction class:)
}

double Initial_Electron_State::Radial_Integral(double r) const
{
	std::function<double(double)> integrand = [this](double rprime) {
		double R = Radial_Wavefunction(rprime);
		return rprime * rprime * R * R;
	};
	return libphysica::Integrate_Gauss_Legendre(integrand, 0.0, r, 1000);
}

void Initial_Electron_State::Print_Summary(unsigned int mpi_rank) const
{
	if(mpi_rank == 0)
	{
		std::cout << SEPARATOR
				  << Orbital_Name() << " - Summary" << std::endl
				  << std::endl
				  << "Binding energy [eV]:\t" << In_Units(binding_energy, eV) << std::endl
				  << "Z_effective:\t\t" << Z_eff << std::endl
				  << std::endl
				  << "n_lj\tZ_lj\tC_nlj" << std::endl;
		for(unsigned int i = 0; i < C_nlj.size(); i++)
			std::cout << n_lj[i] << "\t" << Z_lj[i] << "\t" << C_nlj[i] << std::endl;
	}
}

}	// namespace DarkART