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

TEST(TestSpecialFunctions, TestHypergeometricFunction)
{
	// ARRANGE
	double a = 2.0;
	double b = 3.0;
	double z = 0.4;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(1, 1, 2, -z).real(), log(1.0 + z) / z);
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(a, b, b, z).real(), std::pow(1.0 - z, -a));
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(0.5, 0.5, 1.5, z * z).real(), asin(z) / z);
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(1.2, 32.2, 14.2, 1.2).real(), -16290.21576);
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(1.2, 32.2, 14.2, 1.2).imag(), 11835.53454);
}

TEST(TestSpecialFunctions, TestHypergeometricFunction2)
{
	// ARRANGE
	double a0	  = Bohr_Radius;
	double Z	  = 4.3;
	double k	  = 0.5 * keV;
	double a	  = 3.0;
	double b	  = 3.5;
	double c	  = 2.5;
	double z	  = -k * k * a0 * a0 / Z / Z;
	double result = -std::pow(Z, 6) / 5.0 * (k * k * a0 * a0 - 5.0 * Z * Z) / std::pow(k * k * a0 * a0 + Z * Z, 4);
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(Hypergeometric_2F1(a, b, c, z).real(), result);
}