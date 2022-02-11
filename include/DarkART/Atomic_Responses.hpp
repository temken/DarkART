#ifndef __Atomic_Responses_hpp_
#define __Atomic_Responses_hpp_

#include "Radial_Integrator.hpp"
#include "Wavefunctions_Final.hpp"
#include "Wavefunctions_Initial.hpp"

namespace DarkART
{
extern double Atomic_Response_Function(unsigned int response, double k_final, double q, Radial_Integrator& radial_integrator, int& l_convergence);
extern double Atomic_Response_Function(unsigned int response, double k_final, double q, Radial_Integrator& radial_integrator);

extern double Atomic_Response_Function(unsigned int response, double k_final, double q, const Initial_Electron_State& bound_electron, const Final_Electron_State& final_state);
extern double Atomic_Response_Function(unsigned int response, double k_final, double q, const Initial_Electron_State& bound_electron, const Final_Electron_State& final_state, int& l_convergence);
}	// namespace DarkART

#endif