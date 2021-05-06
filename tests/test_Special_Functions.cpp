#include "gtest/gtest.h"

#include "Special_Functions.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;
using namespace std::complex_literals;

TEST(TestSpecialFunction, TestSphericalHarmonics)
{
	// ARRANGE
	int l						= 5;
	int m						= 3;
	double theta				= 25 * deg;
	double phi					= 13 * deg;
	std::complex<double> result = -0.129726 - 1.0i * 0.10505;
	double tol					= 1.0e-6;
	// ACT & ASSERT
	EXPECT_NEAR(Spherical_Harmonics(1, 0, theta, phi).real(), 0.5 * sqrt(3.0 / M_PI) * cos(theta), tol);

	EXPECT_NEAR(Spherical_Harmonics(l, m, theta, phi).real(), result.real(), tol);
	EXPECT_NEAR(Spherical_Harmonics(l, m, theta, phi).imag(), result.imag(), tol);

	EXPECT_NEAR(Spherical_Harmonics(l, m, 0.0, 0.0).imag(), 0.0, tol);
	EXPECT_NEAR(Spherical_Harmonics(l, m, 0.0, 0.0).imag(), 0.0, tol);

	EXPECT_NEAR(Spherical_Harmonics(l, 2 * l, theta, phi).imag(), 0.0, tol);
	EXPECT_NEAR(Spherical_Harmonics(l, 2 * l, theta, phi).imag(), 0.0, tol);
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsY)
{
	// ARRANGE
	// ACT
	//ASSERT
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsPsi)
{
	// ARRANGE
	// ACT
	//ASSERT
}

TEST(TestSpecialFunction, TestGauntCoefficients)
{
	// ARRANGE
	// ACT
	//ASSERT
}

TEST(TestSpecialFunction, TestSphericalBessel)
{
	// ARRANGE
	double x = 1.5;
	//ACT & ASSERT
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(0, x), sin(x) / x);
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(1, x), sin(x) / x / x - cos(x) / x);
	EXPECT_FLOAT_EQ(Spherical_Bessel_jL(10, x), 3.99344e-9);
}

TEST(TestSpecialFunction, TestCoulombWave)
{
	// ARRANGE
	int L		  = 10;
	double eta	  = -0.5;
	double rho	  = 20.0;
	double result = 1.0587974;
	// ACT
	//ASSERT
	EXPECT_FLOAT_EQ(Coulomb_Wave(L, eta, rho), result);
}