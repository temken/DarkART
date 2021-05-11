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
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Y	= Vector_Spherical_Harmonics_Y(l, m, theta, phi);
	auto Ym = Vector_Spherical_Harmonics_Y(l, -m, theta, phi);
	//ASSERT
	for(int i = 0; i < 3; i++)
	{
		EXPECT_NEAR(Ym[i].real(), std::pow(-1.0, m) * std::conj(Y[i]).real(), 1.0e-6);
		EXPECT_NEAR(Ym[i].imag(), std::pow(-1.0, m) * std::conj(Y[i]).imag(), 1.0e-6);
	}
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsPsi)
{
	// ARRANGE
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Psi  = Vector_Spherical_Harmonics_Psi(l, m, theta, phi);
	auto Psim = Vector_Spherical_Harmonics_Psi(l, -m, theta, phi);
	//ASSERT
	for(int i = 0; i < 3; i++)
	{
		EXPECT_NEAR(Psim[i].real(), std::pow(-1.0, m) * std::conj(Psi[i]).real(), 1.0e-6);
		EXPECT_NEAR(Psim[i].imag(), std::pow(-1.0, m) * std::conj(Psi[i]).imag(), 1.0e-6);
	}
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsOrthogonality)
{
	// ARRANGE
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Y	 = Vector_Spherical_Harmonics_Y(l, m, theta, phi);
	auto Psi = Vector_Spherical_Harmonics_Psi(l, m, theta, phi);
	// ASSERT
	EXPECT_NEAR((Y[0] * Psi[0] + Y[1] * Psi[1] + Y[2] * Psi[2]).real(), 0.0, 1e-16);
	EXPECT_NEAR((Y[0] * Psi[0] + Y[1] * Psi[1] + Y[2] * Psi[2]).imag(), 0.0, 1e-16);
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsYFailure)
{
	// ACT & ASSERT
	ASSERT_DEATH({ auto comp = VSH_Y_Component(5, 5, 3, 5, 3); }, "");
}

TEST(TestSpecialFunction, TestVectorialSphericalHarmonicsPsiFailure)
{
	// ACT & ASSERT
	ASSERT_DEATH({ auto comp = VSH_Psi_Component(5, 5, 3, 5, 3); }, "");
}

TEST(TestSpecialFunction, TestGauntCoefficients)
{
	// ARRANGE
	int l_1 = 1;
	int l_2 = 2;
	int l_3 = 3;
	int m_1 = 0;
	int m_2 = -1;
	int m_3 = +1;
	//ACT & ASSERT
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, m_1, m_2, m_3), -sqrt(6.0 / 35.0 / M_PI));
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, 0, 1, 1), 0.0);
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(l_1, l_2, l_3, 0, 0, 0), 3.0 / 2.0 * sqrt(3.0 / 35.0 / M_PI));
	EXPECT_FLOAT_EQ(Gaunt_Coefficient(1, 1, 1, 0, 0, 0), 0.0);
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
	EXPECT_FLOAT_EQ(Coulomb_Wave(L, eta, 0.0), 0.0);
	EXPECT_FLOAT_EQ(Coulomb_Wave(0, 0.0, rho), sin(rho));
}