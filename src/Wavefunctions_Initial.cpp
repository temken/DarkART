#include "DarkARC/Wavefunctions_Initial.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Special_Functions.hpp"
#include "version.hpp"

namespace DarkARC
{
using namespace std::complex_literals;
using namespace libphysica::natural_units;
using namespace boost::math::quadrature;
using boost::math::factorial;

double a0 = Bohr_Radius;
double au = 27.211386245988 * eV;

// 1. Initial state: Roothaan-Hartree-Fock Ground-State Atomic Wave Functions

void Initial_Electron_State::Import_RHF_Coefficients()
{
	std::string filepath = TOP_LEVEL_DIR "data/" + Orbital_Name() + ".txt";
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

	double norm = Normalization();
	if(std::fabs(1.0 - norm) > 0.01)
	{
		std::cout << "Error in Initial_Electron_State(): Normalization of " << element_name << " = " << norm << " != 1.0" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, int N, int L)
: element_name(element), n(N), l(L)
{
	Import_RHF_Coefficients();
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, std::string shell_name)
: element_name(element)
{
	n = shell_name[0] - '0';
	for(l = 0; l < l_orbital_names.size(); l++)
		if(shell_name[1] == l_orbital_names[l][0])
			break;
	Import_RHF_Coefficients();
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

double Initial_Electron_State::Normalization() const
{
	std::function<double(double)> integrand = [this](double r) {
		double R = Radial_Wavefunction(r);
		return r * r * R * R;
	};
	// Integrate stepwise
	double stepsize	 = Bohr_Radius;
	double integral	 = 0.0;
	double epsilon_1 = 1.0, epsilon_2 = 1.0;
	double tolerance = 1.0e-6;
	for(unsigned int i = 0; epsilon_1 > tolerance || epsilon_2 > tolerance; i++)
	{
		epsilon_2				= epsilon_1;
		double new_contribution = gauss_kronrod<double, 31>::integrate(integrand, i * stepsize, (i + 1) * stepsize, 5, 1e-9);
		integral += new_contribution;
		epsilon_1 = std::fabs(new_contribution / integral);
	}
	return integral;
}

double Initial_Electron_State::Radial_Integral(double r) const
{
	std::function<double(double)> integrand = [this](double rprime) {
		double R = Radial_Wavefunction(rprime);
		return rprime * rprime * R * R;
	};
	return gauss_kronrod<double, 31>::integrate(integrand, 0.0, r, 5, 1e-9);
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

}	// namespace DarkARC