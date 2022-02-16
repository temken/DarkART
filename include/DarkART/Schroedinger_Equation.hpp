#ifndef __Schroedinger_Equation_hpp_
#define __Schroedinger_Equation_hpp_

#include <vector>

#include "libphysica/Numerics.hpp"

namespace DarkART
{

typedef std::vector<double> state_type;
class Schroedinger_Solver
{

	double k_final;
	unsigned int l_final;

	libphysica::Interpolation Z_effective;

  public:
	Schroedinger_Solver(double k, unsigned int l);

	libphysica::Interpolation Solve_Radial_Equation();

	void operator()(const state_type& x, state_type& dxdt, const double /* t */);
};

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

}	// namespace DarkART

#endif