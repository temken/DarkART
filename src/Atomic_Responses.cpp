#include "Atomic_Responses.hpp"

#include <complex>

// Headers from libphysica
#include "Statistics.hpp"

#include <gsl/gsl_sf_coupling.h>

namespace DarkARC
{
using namespace std::complex_literals;

double Radial_Integral(unsigned int i, double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int l_final, unsigned int L)
{
	return 0.0;
}

double Gaunt_Coefficient(int j1, int j2, int j3, int m1, int m2, int m3)
{
	return sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0)) / sqrt(4.0 * M_PI) * gsl_sf_coupling_3j(j1, j2, j3, 0, 0, 0) * gsl_sf_coupling_3j(j1, j2, j3, m1, m2, m3);
}

std::complex<double> Scalar_Atomic_Formfactor(double q, const Initial_Electron_State& bound_electron, int m, double k_final, unsigned int l_final, int m_final)
{
	std::complex<double> f_12 = 0.0;
	for(unsigned int L = std::fabs(bound_electron.l - l_final); L < bound_electron.l + l_final; L++)
	{
		double radial_integral = Radial_Integral(1, k_final, q, bound_electron, l_final, L);
		f_12 += std::pow(1.0i, L) * radial_integral * pow(-1.0, m_final) * sqrt(2.0 * L + 1.0) * Gaunt_Coefficient(bound_electron.l, l_final, L, m, m_final, 0);
	}
	return sqrt(4.0 * M_PI) * f_12;
}

std::vector<std::complex<double>> Vectorial_Atomic_Formfactor(double q, const Initial_Electron_State& bound_electron, int m, double k_final, unsigned int l_final, int m_final)
{
	return {};
}

double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response)
{
	double convergence_level = 0.1;
	double response_function = 0.0;
	std::vector<double> terms;
	bool function_converged = false;
	unsigned int l_prime	= 0;
	while(!function_converged)
		for(int m = -bound_electron.l; m <= bound_electron.l; m++)
			for(int m_prime = -l_prime; m_prime <= l_prime; m_prime++)
			{
				double new_term = 1.0;
				terms.push_back(new_term);
				response_function += new_term;
				// Test for convergence
				if(terms.size() >= 2 * (bound_electron.l + 1))
				{
					std::vector<double> aux(terms.end() - 2 * (bound_electron.l + 1), terms.end());
					double mean = libphysica::Arithmetic_Mean(aux);
					if(mean < convergence_level * response_function / terms.size())
						function_converged = true;
				}
			}
	return response_function;
}

}	// namespace DarkARC
