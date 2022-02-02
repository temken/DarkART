#include "DarkARC/Atomic_Responses.hpp"

#include <complex>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

#include "DarkARC/Special_Functions.hpp"
#include "DarkARC/Wavefunctions_Final.hpp"

namespace DarkARC
{
using namespace std::complex_literals;
using namespace libphysica::natural_units;

double Radial_Integral(unsigned int integral_index, double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int l_final, unsigned int L)
{
	std::function<double(double)> integrand = [integral_index, &bound_electron, L, q, k_final, l_final](double r) {
		switch(integral_index)
		{
			case 1:
				return r * r * bound_electron.Radial_Wavefunction(r) * Radial_Wavefunction_Final(k_final, l_final, bound_electron.Z_eff, r) * Spherical_Bessel_jL(L, q * r);
			case 2:
				return r * r * bound_electron.Radial_Wavefunction_Derivative(r) * Radial_Wavefunction_Final(k_final, l_final, bound_electron.Z_eff, r) * Spherical_Bessel_jL(L, q * r);
			case 3:
				return r * bound_electron.Radial_Wavefunction(r) * Radial_Wavefunction_Final(k_final, l_final, bound_electron.Z_eff, r) * Spherical_Bessel_jL(L, q * r);
			default:
				std::cerr << "Error in Radial_Integral(): Integral I_" << integral_index << " not defined." << std::endl;
				std::exit(EXIT_FAILURE);
		}
	};
	// Integrate stepwise
	double stepsize	 = 0.5 * Bohr_Radius;
	double integral	 = 0.0;
	double epsilon_1 = 1.0, epsilon_2 = 1.0;
	double tolerance = 1.0e-6;
	unsigned int i;
	for(i = 0; epsilon_1 > tolerance || epsilon_2 > tolerance; i++)
	{
		epsilon_2				= epsilon_1;
		double new_contribution = boost::math::quadrature::gauss_kronrod<double, 41>::integrate(integrand, i * stepsize, (i + 1) * stepsize, 5, 1e-9);

		integral += new_contribution;
		epsilon_1 = std::fabs(new_contribution / integral);
	}
	return integral;
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

std::complex<double> Vectorial_Atomic_Formfactor(int component, double q, const Initial_Electron_State& bound_electron, int m, double k_final, int l_final, int m_final)
{
	std::complex<double> f_12 = 0.0;
	for(int l_hat : {bound_electron.l - 1, bound_electron.l + 1})
		for(int L = std::fabs(l_hat - l_final); L <= l_hat + l_final; L++)
		{
			double radial_integral_2 = Radial_Integral(2, k_final, q, bound_electron, l_final, L);
			double radial_integral_3 = Radial_Integral(3, k_final, q, bound_electron, l_final, L);
			for(int m_hat = m - 1; m_hat < m + 2; m_hat++)
				f_12 += std::pow(1.0i, L) * (VSH_Y_Component(component, bound_electron.l, m, l_hat, m_hat) * radial_integral_2 + VSH_Psi_Component(component, bound_electron.l, m, l_hat, m_hat) * radial_integral_3) * std::pow(-1.0, m_final) * sqrt(4.0 * M_PI * (2.0 * L + 1.0)) * Gaunt_Coefficient(l_hat, l_final, L, m_hat, -m_final, 0);
		}
	return 1.0i / mElectron * f_12;
}

// double Transition_Response_Function(response,element,n,l,m,kPrime,lPrime,mPrime,q):
double Transition_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, int m, unsigned int l_final, int m_final, unsigned int response)
{
	double W_12;
	std::complex<double> f, f_1, f_2, f_3;
	switch(response)
	{
		case 1:
			f	 = Scalar_Atomic_Formfactor(q, bound_electron, m, k_final, l_final, m_final);
			W_12 = std::norm(f);
			break;
		case 2:
			f	 = Scalar_Atomic_Formfactor(q, bound_electron, m, k_final, l_final, m_final);
			f_3	 = Vectorial_Atomic_Formfactor(2, q, bound_electron, m, k_final, l_final, m_final);
			W_12 = (q / mElectron * f * std::conj(f_3)).real();
			break;
		case 3:
			f_1	 = Vectorial_Atomic_Formfactor(0, q, bound_electron, m, k_final, l_final, m_final);
			f_2	 = Vectorial_Atomic_Formfactor(1, q, bound_electron, m, k_final, l_final, m_final);
			f_3	 = Vectorial_Atomic_Formfactor(2, q, bound_electron, m, k_final, l_final, m_final);
			W_12 = std::norm(f_1) + std::norm(f_2) + std::norm(f_3);
			break;
		case 4:
			f_3	 = Vectorial_Atomic_Formfactor(2, q, bound_electron, m, k_final, l_final, m_final);
			W_12 = q * q / mElectron / mElectron * std::norm(f_3);
			break;
		default:
			std::cerr << "Error in Transition_Response_Function(): Response " << response << " out of bound." << std::endl;
			std::exit(EXIT_FAILURE);
	}
	return W_12;
}

double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response, int& l_convergence)
{
	double convergence_level = 0.01;
	double prefactor		 = 4.0 * std::pow(k_final / 2.0 / M_PI, 3.0);
	double response_function = 0.0;
	std::vector<double> terms;
	for(int l_final = 0; l_final < 1000; l_final++)
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
			{
				l_convergence = l_final;
				break;
			}
		}
	}
	// One correction
	if(response == 1)
	{
		double radial_integral_1   = Radial_Integral(1, k_final, q, bound_electron, bound_electron.l, 0);
		double radial_integral_one = Radial_Integral(1, k_final, 0.0, bound_electron, bound_electron.l, 0);
		double correction		   = 4.0 * std::pow(k_final, 3) / std::pow(2 * M_PI, 3) * (2 * bound_electron.l + 1) * (radial_integral_one * radial_integral_one - 2 * radial_integral_one * radial_integral_1);
		if(correction < 0.0)
			response_function += correction;
	}
	return response_function;
}

double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response)
{
	int l_convergence;
	return Atomic_Response_Function(k_final, q, bound_electron, response, l_convergence);
}

}	// namespace DarkARC
