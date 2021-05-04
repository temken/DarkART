#ifndef __Wavefunctions_hpp_
#define __Wavefunctions_hpp_

#include <string>
#include <vector>

namespace DarkARC
{

// 1. Initial state: Roothaan-Hartree-Fock Ground-State Atomic Wave Functions
class Initial_Electron_State
{
  private:
	std::string element_name;

	// RHF coefficients
	std::vector<double> C_nlj, Z_lj, n_lj;
	std::vector<std::string> l_orbital_names = {"s", "p", "d", "f", "g"};

  public:
	int n, l;
	double binding_energy, Z_eff;
	Initial_Electron_State(const std::string& element, int N, int L);

	std::string Orbital_Name() const;
	double Radial_Wavefunction(double r) const;
	double Radial_Wavefunction_Derivative(double r) const;
	double Normalization() const;

	void Print_Summary(unsigned int mpi_rank = 0) const;
};

// 2. Final state: Positive energy continuum solution of Schroedinger equation with hydrogenic potential
extern double Coulomb_Wave(int L, double eta, double rho);
extern double Radial_Wavefunction_Final(double k_final, unsigned l_prime, double Z_eff, double r);

}	// namespace DarkARC

#endif