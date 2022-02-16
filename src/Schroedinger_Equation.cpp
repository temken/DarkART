#include "DarkART/Schroedinger_Equation.hpp"

#include <boost/numeric/odeint.hpp>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

namespace DarkART
{
using namespace libphysica::natural_units;
using namespace boost::numeric::odeint;

Schroedinger_Solver::Schroedinger_Solver(double k, unsigned int l)
: k_final(k), l_final(l)
{
	Z_effective = libphysica::Interpolation(
		libphysica::Log_Space(1e-4 * Bohr_Radius, 51.0 * Bohr_Radius, 5),
		std::vector<double>(5, 1.0));
};

libphysica::Interpolation Schroedinger_Solver::Solve_Radial_Equation()
{
	std::vector<state_type> x_vec;
	std::vector<double> times;
	state_type x(2);
	x[0]					 = 1.0;	  // start at x=1.0, p=0.0
	x[1]					 = 0.0;
	double t_ini			 = 0.0;
	double t_final			 = 40.0;
	double initial_step_size = 0.1;

	size_t steps = integrate(*this,
							 x, t_ini, t_final, initial_step_size,
							 push_back_state_and_time(x_vec, times));

	for(size_t i = 0; i <= steps; i++)
		std::cout << times[i] << '\t' << x_vec[i][0] << '\t' << x_vec[i][1] << std::endl;
	std::vector<double> R_list;
	for(auto& x : x_vec)
		R_list.push_back(x[0]);
	return libphysica::Interpolation(times, R_list);
}

void Schroedinger_Solver::operator()(const state_type& x, state_type& dxdt, const double /* t */)
{
	dxdt[0] = x[1];
	dxdt[1] = -x[0] - k_final * x[1];
}

}	// namespace DarkART