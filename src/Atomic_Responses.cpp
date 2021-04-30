#include "Atomic_Responses.hpp"

#include <complex>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>

namespace DarkARC
{
using namespace std::complex_literals;
using namespace libphysica::natural_units;
using namespace boost::math::quadrature;

double Radial_Integral(unsigned int integral_index, double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int l_final, unsigned int L)
{
	std::function<double(double)> integrand = [&bound_electron, L, q, k_final, l_final](double r) {
		return r * r * bound_electron.Radial_Wavefunction(r) * Radial_Wavefunction_Final(k_final, l_final, bound_electron.Z_eff, r) * gsl_sf_bessel_jl(L, q * r);
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

double Gaunt_Coefficient(int j1, int j2, int j3, int m1, int m2, int m3)
{
	return sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0)) / sqrt(4.0 * M_PI) * gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0) * gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3);
}

std::complex<double> Scalar_Atomic_Formfactor(double q, const Initial_Electron_State& bound_electron, int m, double k_final, int l_final, int m_final)
{
	std::complex<double> f_12 = 0.0;
	for(int L = std::fabs(bound_electron.l - l_final); L <= bound_electron.l + l_final; L++)
		if(m - m_final == 0 && (bound_electron.l + l_final + L) % 2 == 0)
		{
			double radial_integral = Radial_Integral(1, k_final, q, bound_electron, l_final, L);
			f_12 += std::pow(1.0i, L) * radial_integral * pow(-1.0, m_final) * sqrt(2.0 * L + 1.0) * Gaunt_Coefficient(bound_electron.l, l_final, L, m, -m_final, 0);
		}
	return sqrt(4.0 * M_PI) * f_12;
}

std::vector<std::complex<double>> Vectorial_Atomic_Formfactor(double q, const Initial_Electron_State& bound_electron, int m, double k_final, unsigned int l_final, int m_final)
{
	return {};
}

// double Transition_Response_Function(response,element,n,l,m,kPrime,lPrime,mPrime,q):
double Transition_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, int m, unsigned int l_final, int m_final, unsigned int response)
{
	double W_12;
	if(response == 1)
	{
		std::complex<double> f_12 = Scalar_Atomic_Formfactor(q, bound_electron, m, k_final, l_final, m_final);
		W_12					  = std::norm(f_12);
	}
	else
		W_12 = 0.0;
	return W_12;
}

double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response)
{
	double convergence_level = 0.1;
	double prefactor		 = 4.0 * std::pow(k_final / 2.0 / M_PI, 3.0);
	double response_function = 0.0;
	std::vector<double> terms;
	bool function_converged = false;
	int l_final				= 0;
	for(l_final = 0; l_final < 100; l_final++)
	{
		double new_term = 0.0;
		for(int m = -bound_electron.l; m <= bound_electron.l; m++)
			for(int m_final = -l_final; m_final <= l_final; m_final++)
				new_term += prefactor * Transition_Response_Function(k_final, q, bound_electron, m, l_final, m_final, response);
		terms.push_back(new_term);
		response_function += new_term;

		// Test for convergence
		if(terms.size() >= 2 * (bound_electron.l + 1))
		{
			std::vector<double> aux(terms.end() - 2 * (bound_electron.l + 1), terms.end());
			double mean = libphysica::Arithmetic_Mean(aux);
			if(mean < convergence_level * response_function / terms.size())
				break;
		}
		// std::cout << l_final << "\t" << response_function << "\t" << new_term << std::endl;
	}
	std::cout << l_final << std::endl;
	return response_function;
}

}	// namespace DarkARC
