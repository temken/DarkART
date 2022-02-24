#ifndef __Schroedinger_Equation_hpp_
#define __Schroedinger_Equation_hpp_

#include <complex>
#include <vector>

#include "libphysica/Numerics.hpp"

namespace DarkART
{

typedef std::vector<std::complex<double>> state_type;
class Schroedinger_Solver
{
	double r_min, r_max;
	std::vector<double> r_grid;
	double k_final;
	unsigned int l_final;

	libphysica::Interpolation Z_effective;

	std::vector<std::complex<double>> Solve_Radial_Equation_f();

  public:
	Schroedinger_Solver(double kf, unsigned int l);

	libphysica::Interpolation Solve_Radial_Equation();

	void operator()(const state_type& x, state_type& dxdr, const double r);
};

}	// namespace DarkART

#endif