#ifndef __Wavefunctions_Initial_hpp_
#define __Wavefunctions_Initial_hpp_

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

	void Import_RHF_Coefficients();

  public:
	int n, l;
	double binding_energy, Z_eff;
	Initial_Electron_State(const std::string& element, int N, int L);
	Initial_Electron_State(const std::string& element, std::string shell_name);

	std::string Orbital_Name() const;
	double Radial_Wavefunction(double r) const;
	double Radial_Wavefunction_Derivative(double r) const;
	double Normalization() const;
	double Radial_Integral(double r) const;

	void Print_Summary(unsigned int mpi_rank = 0) const;
};

}	// namespace DarkARC

#endif