#ifndef __Radial_Integrator_hpp__
#define __Radial_Integrator_hpp__

#include "DarkARC/Wavefunctions_Final.hpp"
#include "DarkARC/Wavefunctions_Initial.hpp"

namespace DarkARC
{

class Radial_Integrator
{
  protected:
	Initial_Electron_State initial_state;
	Final_Electron_State* final_state = nullptr;

	// 1. Method: Adaptive integration using boost.

	// 2. Method: Integration using tabulated functions.
	bool using_function_tabulation;
	std::vector<double> k_grid, q_grid;
	double r_max;
	std::vector<double> initial_radial_wavefunction_list, initial_radial_wavefunction_derivative_list;

	unsigned int l_final_max = 150, L_max;
	std::vector<std::vector<std::vector<double>>> final_radial_wavefunction_list;
	std::vector<std::vector<std::vector<double>>> bessel_function_list;

	std::vector<std::vector<double>> r_values_and_weights;

	double Radial_Integral_Table(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);

  public:
	double Radial_Integral_Adaptive(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);
	Radial_Integrator(const Initial_Electron_State& ini_state, const Final_Electron_State& fin_state);

	void Tabulate_Functions(unsigned int r_points, const std::vector<double>& k_list, const std::vector<double>& q_list, int threads);

	double Radial_Integral(int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);
};

}	// namespace DarkARC

#endif