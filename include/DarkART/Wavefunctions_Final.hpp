#ifndef __Wavefunctions_Final_hpp_
#define __Wavefunctions_Final_hpp_

#include <string>
#include <vector>

#include "libphysica/Numerics.hpp"

#include "DarkART/Wavefunctions_Initial.hpp"

namespace DarkART
{
// 1. Base class for final state wave function
class Final_Electron_State
{
  protected:
  public:
	Final_Electron_State() {};

	virtual double Radial_Wavefunction(double r, double k_final, unsigned int l_final) { return 0.0; };

	virtual Final_Electron_State* Clone() const;
};

// 2. Positive energy continuum solution of Schroedinger equation with hydrogenic potential (i.e. constant Z_eff)
extern double Radial_Wavefunction_Hydrogenic(double k_final, unsigned l_prime, double Z_eff, double r);

class Final_Electron_State_Hydrogenic : public Final_Electron_State
{
  protected:
  public:
	double Z_effective;
	Final_Electron_State_Hydrogenic(double Z_eff = 1.0);

	virtual double Radial_Wavefunction(double r, double k_final, unsigned int l_final) override;

	virtual Final_Electron_State_Hydrogenic* Clone() const override;
};

// 3. Positive energy continuum solution of Schroedinger equation for a given potential Z_eff(r)
class Final_Electron_State_Schroedinger : Final_Electron_State
{
  protected:
	Initial_Electron_State initial_state;
	double r_min, r_max;
	libphysica::Interpolation Z_effective_interpolation;

  public:
	Final_Electron_State_Schroedinger(Initial_Electron_State& ini_state, double Z_eff = 1.0);

	double Z_effective(double r);

	void Determine_Z_effective();

	void Solve_Schroedinger_Equation();

	virtual double Radial_Wavefunction(double r, double k_final, unsigned int l_final) override;

	virtual Final_Electron_State_Schroedinger* Clone() const override;
};

}	// namespace DarkART

#endif