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
	Final_Electron_State* final_state;

	bool using_function_tabulation;

	double Radial_Integral_Adaptive(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);
	double Radial_Integral_Table(unsigned int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);

  public:
	Radial_Integrator(const Initial_Electron_State& ini_state, const Final_Electron_State& fin_state);

	void Tabulate_Functions();

	double Radial_Integral(int integral_index, double k_final, double q, unsigned int l_final, unsigned int L);
};

}	// namespace DarkARC

#endif