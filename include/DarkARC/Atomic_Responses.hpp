#ifndef __Atomic_Responses_hpp_
#define __Atomic_Responses_hpp_

#include "Wavefunctions_Initial.hpp"

namespace DarkARC
{

extern double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response);
extern double Atomic_Response_Function(double k_final, double q, const Initial_Electron_State& bound_electron, unsigned int response, int& l_convergence);

}	// namespace DarkARC

#endif