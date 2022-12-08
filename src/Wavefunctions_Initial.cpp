#include "DarkART/Wavefunctions_Initial.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkART/Special_Functions.hpp"
#include "DarkART/version.hpp"

namespace DarkART
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

Initial_Electron_State::Initial_Electron_State()
: element_name("none"), n(0), l(0), binding_energy(0.0), Z_eff(0.0)
{
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, int N, int L)
: element_name(element), n(N), l(L)
{
	Import_RHF_Coefficients();
	Check_Normalization();
}

Initial_Electron_State::Initial_Electron_State(const std::string& element, std::string shell_name)
: element_name(element)
{
	n = shell_name[0] - '0';
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

double Initial_Electron_State::Normalization() const
{
	std::function<double(double)> integrand = [this](double r) {
		double R = Radial_Wavefunction(r);
		return r * r * R * R;
	};
	// Integrate with Gauss Legendre
	return libphysica::Integrate_Gauss_Legendre(integrand, 0.0, 50.0 * Bohr_Radius, 1000);
}

double Initial_Electron_State::Radial_Integral(double r) const
{
	std::function<double(double)> integrand = [this](double rprime) {
		double R = Radial_Wavefunction(rprime);
		return rprime * rprime * R * R;
	};
	return gauss_kronrod<double, 31>::integrate(integrand, 0.0, r, 5, 1e-9);
}

std::complex<double> Initial_Electron_State::Radial_Wavefunction_Momentum(double p) const
{
	using namespace std::complex_literals;
	std::complex<double> R_nl = 0.0;
	for(unsigned int j = 0; j < C_nlj.size(); j++)
	{
		// Hypergeometric function 2F1(a,b,c,z)
		double a								= 0.5 * (2.0 + l + n_lj[j]);
		double b								= 0.5 * (3.0 + l + n_lj[j]);
		double c								= 1.5 + l;
		double z								= -std::pow(p * Bohr_Radius / Z_lj[j], 2.0);
		std::complex<double> hypergeometric_2F1 = Hypergeometric_2F1(a, b, c, z);

		R_nl += std::pow(2.0 * M_PI * Bohr_Radius / Z_lj[j], 1.5) * C_nlj[j] * std::pow(2.0, -l + n_lj[j]) * factorial<double>(1.0 + n_lj[j] + l) / sqrt(factorial<double>(2.0 * n_lj[j])) * std::pow(1i * p * Bohr_Radius / Z_lj[j], l) * hypergeometric_2F1 / tgamma(1.5 + l);
	}
	return R_nl;
}

double Initial_Electron_State::Normalization_Momentum() const
{
	std::function<double(double)> integrand = [this](double k) {
		std::complex<double> R = Radial_Wavefunction_Momentum(k);
		return k * k * std::norm(R);
	};
	// Integrate with Gauss Legendre
	return std::pow(2.0 * M_PI, -3) * libphysica::Integrate_Gauss_Legendre(integrand, 0.0, 3000.0 * keV, 1000);
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