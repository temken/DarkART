#include "DarkART/Schroedinger_Equation.hpp"

#include <boost/numeric/odeint.hpp>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

namespace DarkART
{
using namespace libphysica::natural_units;

struct push_back_state_and_time
{
	std::vector<state_type>& m_states;
	std::vector<double>& m_times;

	push_back_state_and_time(std::vector<state_type>& states, std::vector<double>& times)
	: m_states(states), m_times(times) {}

	void operator()(const state_type& x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}
};

Schroedinger_Solver::Schroedinger_Solver(double k, unsigned int l)
: k_final(k), l_final(l)
{
	Z_effective = libphysica::Interpolation(
		libphysica::Linear_Space(0.0, 51.0 * Bohr_Radius, 5),
		std::vector<double>(5, 1.0));
};

libphysica::Interpolation Schroedinger_Solver::Solve_Radial_Equation()
{
	std::vector<state_type> x_vec;
	std::vector<double> times;
	state_type x(2);
	x[0]					 = 0.0;
	x[1]					 = 0.01;
	double t_ini			 = 1e-4 * Bohr_Radius;
	double t_final			 = 10.0 * Bohr_Radius;
	double initial_step_size = 0.01 * Bohr_Radius;

	boost::numeric::odeint::integrate(*this,
									  x, t_ini, t_final, initial_step_size,
									  push_back_state_and_time(x_vec, times));
	std::vector<double> R_list;
	for(auto& x : x_vec)
		R_list.push_back(x[0]);
	return libphysica::Interpolation(times, R_list);
}

void Schroedinger_Solver::operator()(const state_type& x, state_type& dxdt, const double t)
{
	double mu	 = mElectron;
	double E	 = k_final * k_final / 2.0 / mu;
	double e	 = Elementary_Charge;
	double Z_eff = 1.0;
	dxdt[0]		 = x[1];
	dxdt[1]		 = -2.0 * x[1] / t + (-2.0 * mu * Z_eff * e * e / t + l_final * (l_final + 1.0) / t / t - 2.0 * mu * E) * x[0];
}

}	// namespace DarkART