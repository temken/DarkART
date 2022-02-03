#include "DarkARC/Radial_Integrator.hpp"

#include <cmath>
#include <functional>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Special_Functions.hpp"

namespace DarkARC
{

using namespace libphysica::natural_units;

double Radial_Integrator::Radial_Integral_Adaptive(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L)
{
	std::function<double(double)> integrand = [this, integral_index, L, q, k_final, l_final](double r) {
		switch(integral_index)
		{
			case 1:
				return r * r * initial_state.Radial_Wavefunction(r) * final_state->Radial_Wavefunction(r, k_final, l_final) * Spherical_Bessel_jL(L, q * r);
			case 2:
				return r * r * initial_state.Radial_Wavefunction_Derivative(r) * final_state->Radial_Wavefunction(r, k_final, l_final) * Spherical_Bessel_jL(L, q * r);
			case 3:
				return r * initial_state.Radial_Wavefunction(r) * final_state->Radial_Wavefunction(r, k_final, l_final) * Spherical_Bessel_jL(L, q * r);
			default:
				std::cerr << "Radial_Integrator::Radial_Integral_Adaptive(): Integral I_" << integral_index << " not defined." << std::endl;
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

double Radial_Integrator::Radial_Integral_Table(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L)
{
	return 0.0;
}

Radial_Integrator::Radial_Integrator(const Initial_Electron_State& ini_state, const Final_Electron_State& fin_state)
: initial_state(ini_state), using_function_tabulation(false)
{
	final_state = fin_state.Clone();
}

void Radial_Integrator::Tabulate_Functions()
{
	using_function_tabulation = true;
	// tabulate...
}

double Radial_Integrator::Radial_Integral(int integral_index, double k_final, double q, unsigned int l_final, unsigned int L)
{
	if(using_function_tabulation)
		return Radial_Integral_Table(integral_index, k_final, q, l_final, L);
	else
		return Radial_Integral_Adaptive(integral_index, k_final, q, l_final, L);
}

}	// namespace DarkARC