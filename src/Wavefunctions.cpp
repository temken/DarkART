#include "Wavefunctions.hpp"

#include <cmath>
#include <complex>
#include <fstream>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>

// Headers from libphysica
#include "Natural_Units.hpp"

#include "version.hpp"

namespace DarkARC
{
using namespace std::complex_literals;
using namespace libphysica::natural_units;
using namespace boost::math::quadrature;
using boost::math::factorial;
using boost::math::hypergeometric_1F1;

double a0 = Bohr_Radius;
double au = 27.211386245988 * eV;

// 1. Initial state: Roothaan-Hartree-Fock Ground-State Atomic Wave Functions
Initial_Electron_State::Initial_Electron_State(const std::string& element, unsigned int N, unsigned int L)
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
	Z_eff = sqrt(-binding_energy / (13.6 * eV)) * n;
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

double Radial_Wavefunction_Final(double k_final, unsigned l_final, double Z_eff, double r)
{
	// 1. Compute the hypergeometric function 1F1(l'+ 1 + i Z / k' /a0, 2 l' + 1, 2 i k' r)
	// 1.1 Compute |Gamma(l' + 1 + i Zeff / k' / a0)|
	double b		 = Z_eff / k_final / a0;
	double gamma_sqr = M_PI * b / sinh(M_PI * b);
	for(int k = 1; k <= l_final; k++)
		gamma_sqr *= (k * k + b * b);
	double gamma = sqrt(gamma_sqr);
	// 1.2 Use 13.4.1 of https://dlmf.nist.gov/13.4#E1 to compute 1F1 exp(-ikr) = |1F1|
	std::function<double(double)> integrand = [k_final, l_final, Z_eff, r](double t) {
		std::complex<double> z				   = 1.0 * l_final + 1.0i * Z_eff / k_final / a0;
		std::complex<double> complex_integrand = std::pow(t, z) * std::pow(1.0 - t, std::conj(z)) * std::exp(2.0i * k_final * r * t - 1.0i * k_final * r);
		return complex_integrand.real();
	};
	double integral				   = gauss_kronrod<double, 31>::integrate(integrand, 0.0, 1.0, 5, 1e-9);
	double hypergeometric_1F1_norm = tgamma(2.0 * l_final + 2.0) / gamma_sqr * integral;
	double R_final				   = 4.0 * M_PI * std::pow(2.0 * k_final * r, l_final) * exp(M_PI * Z_eff / 2.0 / k_final / a0) / factorial<double>(2.0 * l_final + 1.0) * gamma * hypergeometric_1F1_norm;
	return R_final;
}

}	// namespace DarkARC