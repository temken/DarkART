#include "gtest/gtest.h"

#include <cmath>
#include <complex>

#include "DarkART/Special_Functions.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace DarkART;
using namespace libphysica::natural_units;
using namespace std::complex_literals;

TEST(TestSpecialFunctions, TestGauntCoefficients)
{
	// ARRANGE
	int l_1 = 1;
	int l_2 = 2;
	int l_3 = 3;
	int m_1 = 0;
	int m_2 = -1;
	int m_3 = +1;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, m_1, m_2, m_3), -sqrt(6.0 / 35.0 / M_PI));
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, 0, 1, 1), 0.0);
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, 0, 0, 0), 3.0 / 2.0 * sqrt(3.0 / 35.0 / M_PI));
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(1, 1, 1, 0, 0, 0), 0.0);
}

TEST(TestSpecialFunctions, TestSphericalBessel)
{
	// ARRANGE
	double x = 1.5;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(0, x), sin(x) / x);
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(1, x), sin(x) / x / x - cos(x) / x);
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(10, x), 3.99344e-9);
}

TEST(TestSpecialFunctions, TestCoulombWave)
{
	// ARRANGE
	int L		  = 10;
	double eta	  = -0.5;
	double rho	  = 20.0;
	double result = 1.0587974;
	// ACT
	// ASSERT
	EXPECT_FLOAT_EQ(Coulomb_Wave(L, eta, rho), result);
	EXPECT_FLOAT_EQ(Coulomb_Wave(L, eta, 0.0), 0.0);
	EXPECT_FLOAT_EQ(Coulomb_Wave(0, 0.0, rho), sin(rho));
}