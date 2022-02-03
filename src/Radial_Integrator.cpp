#include "DarkARC/Radial_Integrator.hpp"

#include <cmath>
#include <functional>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

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

void Radial_Integrator::Tabulate_Functions(unsigned int r_points, const std::vector<double>& k_list, const std::vector<double>& q_list)
{
	using_function_tabulation = true;

	k_grid = k_list;
	q_grid = q_list;

	r_max		= 100.0 * Bohr_Radius;
	r_grid		= libphysica::Linear_Space(r_points, 0.0, r_max);
	delta_r		= r_grid[1] - r_grid[0];
	l_final_max = 150;
	L_max		= initial_state.l + 1 + l_final_max;

	std::cout << "hmm" << std::endl;
	initial_radial_wavefunction_list = std::vector<double>(r_grid.size(), 0.0);
	final_radial_wavefunction_list	 = std::vector<std::vector<std::vector<double>>>(k_grid.size(), std::vector<std::vector<double>>(l_final_max, std::vector<double>(r_grid.size(), 0.0)));
	bessel_function_list			 = std::vector<std::vector<std::vector<double>>>(q_grid.size(), std::vector<std::vector<double>>(L_max, std::vector<double>(r_grid.size(), 0.0)));

	std::cout << "hmm" << std::endl;
	// 1.  Tabulate initial radial wavefunctions
	std::cout << "Tabulate initial radial wavefunctions..." << std::endl;
	for(unsigned int ri = 0; ri < r_grid.size(); ri++)
		initial_radial_wavefunction_list[ri] = initial_state.Radial_Wavefunction(r_grid[ri]);

	// 2. Tabulate final radial wavefunctions
	std::cout << "Tabulate final radial wavefunctions..." << std::endl;
	for(unsigned int ki = 0; ki < k_grid.size(); ki++)
		for(unsigned int li = 0; li <= l_final_max; li++)
			for(unsigned int ri = 0; ri < r_grid.size(); ri++)
				final_radial_wavefunction_list[ki][li][ri] = final_state->Radial_Wavefunction(r_grid[ri], k_grid[ki], li);

	// 3. Tabulate spherical Bessel functions
	std::cout << "Tabulate bessel functions..." << std::endl;
	for(unsigned int qi = 0; qi < q_grid.size(); qi++)
		for(unsigned int L = 0; L <= L_max; L++)
			for(unsigned int ri = 0; ri < r_grid.size(); ri++)
				bessel_function_list[qi][L][ri] = Spherical_Bessel_jL(L, q_grid[qi] * r_grid[ri]);
}

double Radial_Integrator::Radial_Integral(int integral_index, double k_final, double q, unsigned int l_final, unsigned int L)
{
	if(using_function_tabulation)
		return Radial_Integral_Table(integral_index, k_final, q, l_final, L);
	else
		return Radial_Integral_Adaptive(integral_index, k_final, q, l_final, L);
}

}	// namespace DarkARC