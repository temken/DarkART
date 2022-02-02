#ifndef __Wavefunctions_Final_hpp_
#define __Wavefunctions_Final_hpp_

#include <string>
#include <vector>

namespace DarkARC
{
// 1. Base class for final state wave function
class Final_Electron_State
{
  protected:
  public:
	Final_Electron_State() {};

	virtual double Radial_Wavefunction(double r, double k_final, unsigned int l_final) { return 0.0; };
};

// 2. Positive energy continuum solution of Schroedinger equation with hydrogenic potential
extern double Radial_Wavefunction_Hydrogenic(double k_final, unsigned l_prime, double Z_eff, double r);

class Final_Electron_State_Hydrogenic : public Final_Electron_State
{
  protected:
	double Z_effective;

  public:
	Final_Electron_State_Hydrogenic(double Z_eff = 1.0);

	void Fit_Zeff(int n, double binding_energy);

	virtual double Radial_Wavefunction(double r, double k_final, unsigned int l_final) override;
};

}	// namespace DarkARC

#endif