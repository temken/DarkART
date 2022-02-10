#include "gtest/gtest.h"

#include "DarkARC/Wavefunctions_Final.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestWavefunctions, TestFinalElectronState)
{
	double k_final = 10.0 * eV;
	int l_final	   = 3;
	double Z_eff   = 1.0;
	double r	   = Bohr_Radius;

	// ACT
	Final_Electron_State_Hydrogenic state(Z_eff);

	// ASSERT
	ASSERT_FLOAT_EQ(state.Radial_Wavefunction(r, k_final, l_final), 0.74927193);
}

TEST(TestWavefunctions, TestClone)
{
	double k_final = 10.0 * eV;
	int l_final	   = 3;
	double Z_eff   = 3.8;
	double r	   = 1.9 * Bohr_Radius;

	// ACT
	Final_Electron_State_Hydrogenic state_1(Z_eff);
	Final_Electron_State* state_2 = state_1.Clone();

	// ASSERT
	ASSERT_FLOAT_EQ(state_1.Radial_Wavefunction(r, k_final, l_final), state_2->Radial_Wavefunction(r, k_final, l_final));
}

TEST(TestWavefunctions, TestRadialWavefunctionHydrogenic)
{
	// ARRANGE
	double k_final = 10.0 * eV;
	int l_final	   = 3;
	double Z_eff   = 1.0;
	double r	   = Bohr_Radius;
	// ACT & ASSERT
	ASSERT_FLOAT_EQ(Radial_Wavefunction_Hydrogenic(k_final, l_final, Z_eff, r), 0.74927193);
}
