#include "DarkART/Schroedinger_Equation.hpp"

#include <fstream>

#include <boost/numeric/odeint.hpp>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

namespace DarkART
{
using namespace libphysica::natural_units;
using namespace std::complex_literals;

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

Schroedinger_Solver::Schroedinger_Solver(double kf, unsigned int l)
: k_final(kf), l_final(l)
{
	double k	= In_Units(k_final, 1.0 / Bohr_Radius);
	r_min		= 1.0e-4;
	r_max		= 50.0 / k;
	Z_effective = libphysica::Interpolation(
		libphysica::Linear_Space(0.0, r_max, 5),
		std::vector<double>(5, 1.0));
};

std::vector<std::complex<double>> Schroedinger_Solver::Solve_Radial_Equation_f()
{
	double k = In_Units(k_final, 1.0 / Bohr_Radius);

	double dr = -1e-4;	 // * Bohr_Radius;

	std::vector<state_type> x_vec;
	state_type x(2);
	x[0] = 1.0;
	x[1] = 1.0i / r_max / k;

	typedef boost::numeric::odeint::runge_kutta4<state_type> stepper_type;
	boost::numeric::odeint::integrate_const(stepper_type(),
											*this,
											x, r_max, r_min, dr,
											push_back_state_and_time(x_vec, r_grid));
	std::reverse(std::begin(x_vec), std::end(x_vec));
	std::reverse(std::begin(r_grid), std::end(r_grid));

	std::vector<std::complex<double>> f_list;
	for(auto& x : x_vec)
		f_list.push_back(x[0]);
	return f_list;
}

libphysica::Interpolation Schroedinger_Solver::Solve_Radial_Equation()
{
	double k									 = In_Units(k_final, 1.0 / Bohr_Radius);
	std::vector<std::complex<double>> f_solution = Solve_Radial_Equation_f();
	double phase								 = std::arg(f_solution[0]);
	std::vector<double> rRout_list, rRdiv_list;
	for(int i = 0; i < f_solution.size(); i++)
	{
		double r			   = r_grid[i];
		std::complex<double> f = f_solution[i];

		double rRout = 4.0 * M_PI / k * (f * std::exp(1.0i * k * (r - r_min) - phase)).imag();
		double rRdiv = 4.0 * M_PI / k * (f * std::exp(1.0i * k * (r - r_min) - phase)).real();
		rRout_list.push_back(rRout);
		rRdiv_list.push_back(rRdiv);
	}
	double r_divergence = r_min;
	if(l_final != 0)
	{
		auto rRdiv = libphysica::Interpolation(r_grid, rRdiv_list);
		auto func  = [k, &rRdiv](double r) {
			 return 1e4 - std::pow(k / 4.0 / M_PI * rRdiv(r), 2.0);
		};
		r_divergence = libphysica::Find_Root(func, 2.0 * r_min, r_max, 1e-10);
	}
	std::cout << "r_divergence = " << r_divergence << std::endl;
	double r_cut = r_min;
	if(r_divergence > 2.0 * r_min && r_divergence < r_max)
	{
		auto rRout = libphysica::Interpolation(r_grid, rRout_list);
		auto func  = [&rRout](double r) {
			 double rR = rRout(r);
			 return std::log(rR * rR);
		};
		r_cut = libphysica::Find_Minimum(func, 2.0 * r_min, r_divergence);
	}
	std::cout << "r_cut = " << r_cut << std::endl;
	// Tabulate R_out
	std::vector<double> Rout_list;
	for(int i = 0; i < r_grid.size(); i++)
	{
		double r = r_grid[i];
		if(r > r_cut)
			Rout_list.push_back(r * rRout_list[i]);
		else
			Rout_list.push_back(0.0);
	}

	return libphysica::Interpolation(r_grid, Rout_list);
}

void Schroedinger_Solver::operator()(const state_type& x, state_type& dxdr, const double r)
{
	double mu	 = 1.0;
	double e	 = 1.0;
	double Z_eff = 1.0;
	double k	 = In_Units(k_final, 1.0 / Bohr_Radius);

	dxdr[0] = x[1];
	dxdr[1] = -2.0i * k * x[1] - 2.0 * e * e * Z_eff * mu / r * x[0] + l_final * (l_final + 1.0) / r / r * x[0];
}

}	// namespace DarkART