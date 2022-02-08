#ifndef __Radial_Integrator_hpp__
#define __Radial_Integrator_hpp__

#include "DarkARC/Wavefunctions_Final.hpp"
#include "DarkARC/Wavefunctions_Initial.hpp"

namespace DarkARC
{

class Radial_Integrator
{
  protected:
	// 1. Method: Adaptive integration using boost.
	double Radial_Integral_Adaptive(unsigned int integral_index, double k_final, double q, int l_final, int L);

	// 2. Method: Integration using tabulated functions.
	bool using_function_tabulation;
	std::vector<double> k_grid, q_grid;
	std::vector<int> l_final_max, L_max;

	unsigned int r_points;
	std::vector<std::vector<double>> r_values_and_weights;
	double r_max;

	unsigned int l_final_max_max = 400;
	std::vector<double> initial_radial_wavefunction_list, initial_radial_wavefunction_derivative_list;
	std::vector<std::vector<std::vector<double>>> final_radial_wavefunction_list;
	std::vector<std::vector<std::vector<double>>> bessel_function_list;
	void Tabulate_Final_Wavefunction(int l_final, int ki);
	void Tabulate_Bessel_Function(int Lmax, int qi);

	double Radial_Integral_Table(unsigned int integral_index, double k_final, double q, int l_final, int L);

  public:
	Initial_Electron_State initial_state;
	Final_Electron_State* final_state = nullptr;

	Radial_Integrator(const Initial_Electron_State& ini_state, const Final_Electron_State& fin_state);

	void Use_Tabulated_Functions(unsigned int rpoints, const std::vector<double>& k_list, const std::vector<double>& q_list);

	double Radial_Integral(int integral_index, double k_final, double q, int l_final, int L);
};

}	// namespace DarkARC

#endif