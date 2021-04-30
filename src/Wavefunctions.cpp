#include "Wavefunctions.hpp"

#include <cmath>
#include <complex>
#include <fstream>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/factorials.hpp>

// #include <arb_hypgeom.h>
#include <gsl/gsl_sf_coulomb.h>

#include "arb_hypgeom.h"

#include "libphysica/Natural_Units.hpp"

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
Initial_Electron_State::Initial_Electron_State(const std::string& element, int N, int L)
: element_name(element), n(N), l(L)
{
	// Import RHF coefficients from file
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

// 2. Final electron state wavefunction
double Coulomb_Wave_ARB(int L, double eta, double rho)
{
	slong prec;
	prec = 53;
	arb_t F, l, eta_2, rho_2;
	arb_init(F);
	arb_init(l);
	arb_init(eta_2);
	arb_init(rho_2);
	arb_set_d(l, L);
	arb_set_d(eta_2, eta);
	arb_set_d(rho_2, rho);
	arb_hypgeom_coulomb(F, NULL, l, eta_2, rho_2, prec);
	arb_clear(F);
	arb_clear(eta_2);
	arb_clear(rho_2);
	arb_clear(l);
	double re = arf_get_d(arb_midref(F), ARF_RND_NEAR);
	return re;
}

double Coulomb_Wave_GSL(int L, double eta, double rho, int& status)
{
	double fc_array[1];
	double F_exponent[1];
	status = gsl_sf_coulomb_wave_F_array(L, 0, eta, rho, fc_array, F_exponent);
	return fc_array[0];
}

double Coulomb_Wave(int L, double eta, double rho)
{
	int status;
	double cw = Coulomb_Wave_GSL(L, eta, rho, status);
	if(status != 0)
		cw = Coulomb_Wave_ARB(L, eta, rho);
	return cw;
}

double Radial_Wavefunction_Final(double k_final, unsigned l_final, double Z_eff, double r)
{
	double eta = -Z_eff / k_final / a0;
	double rho = k_final * r;
	return 4.0 * M_PI / rho * Coulomb_Wave(l_final, eta, rho);
}

}	// namespace DarkARC