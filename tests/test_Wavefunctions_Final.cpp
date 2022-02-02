#include "gtest/gtest.h"

#include "DarkARC/Wavefunctions_Final.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestWavefunctions, TestRadialWavefunctionFinal)
{
	// ARRANGE
	double k_final = 10.0 * eV;
	int l_final	   = 3;
	double Z_eff   = 1.0;
	double r	   = Bohr_Radius;
	// ACT & ASSERT
	ASSERT_FLOAT_EQ(Radial_Wavefunction_Final(k_final, l_final, Z_eff, r), 0.74927193);
}
