#include "gtest/gtest.h"

#include "DarkART/Wavefunctions_Initial.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

using namespace DarkART;
using namespace libphysica::natural_units;

double au = 27.211386245988 * eV;

TEST(TestWavefunctions, TestConstructor)
{
	// ARRANGE
	std::string element = "Xe";
	int n				= 5;
	int l				= 0;
	// ACT
	Initial_Electron_State Xe5p(element, n, l);
	// ASSERT
	EXPECT_EQ(Xe5p.Orbital_Name(), "Xe_5s");
	EXPECT_EQ(Xe5p.n, 5);
	EXPECT_EQ(Xe5p.l, 0);
	EXPECT_FLOAT_EQ(Xe5p.binding_energy, -25.698624 * eV);
	EXPECT_FLOAT_EQ(Xe5p.Z_eff, sqrt(-2.0 * Xe5p.binding_energy / au) * n);
}

TEST(TestWavefunctions, TestRadialWavefunctionDerivative)
{
	// ACT
	double eps = 1.0e-6 * Bohr_Radius;
	double r   = 0.8 * Bohr_Radius;

	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 1),
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Xe", 4, 2),
		Initial_Electron_State("Xe", 4, 1),
		Initial_Electron_State("Xe", 4, 0),
		Initial_Electron_State("Xe", 3, 2),
		Initial_Electron_State("Xe", 3, 1),
		Initial_Electron_State("Xe", 3, 0),
		Initial_Electron_State("Xe", 2, 1),
		Initial_Electron_State("Xe", 2, 0),
		Initial_Electron_State("Xe", 1, 0),
		Initial_Electron_State("Ar", 3, 1),
		Initial_Electron_State("Ar", 3, 0),
		Initial_Electron_State("Ar", 2, 1),
		Initial_Electron_State("Ar", 2, 0),
		Initial_Electron_State("Ar", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		if(electron.Radial_Wavefunction_Derivative(r) > 0.0)
			EXPECT_GT(electron.Radial_Wavefunction(r + eps), electron.Radial_Wavefunction(r));
		else
			EXPECT_LT(electron.Radial_Wavefunction(r + eps), electron.Radial_Wavefunction(r));
}

TEST(TestWavefunctions, TestNormalizationArgon)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Ar", 3, 1),
		Initial_Electron_State("Ar", 3, 0),
		Initial_Electron_State("Ar", 2, 1),
		Initial_Electron_State("Ar", 2, 0),
		Initial_Electron_State("Ar", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestNormalizationGermanium)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Ge", 4, 1),
		Initial_Electron_State("Ge", 4, 0),
		Initial_Electron_State("Ge", 3, 2),
		Initial_Electron_State("Ge", 3, 1),
		Initial_Electron_State("Ge", 3, 0),
		Initial_Electron_State("Ge", 2, 1),
		Initial_Electron_State("Ge", 2, 0),
		Initial_Electron_State("Ge", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestNormalizationIodine)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("I", 5, 1),
		Initial_Electron_State("I", 5, 0),
		Initial_Electron_State("I", 4, 2),
		Initial_Electron_State("I", 4, 1),
		Initial_Electron_State("I", 4, 0),
		Initial_Electron_State("I", 3, 2),
		Initial_Electron_State("I", 3, 1),
		Initial_Electron_State("I", 3, 0),
		Initial_Electron_State("I", 2, 1),
		Initial_Electron_State("I", 2, 0),
		Initial_Electron_State("I", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestNormalizationSodium)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Na", 3, 0),
		Initial_Electron_State("Na", 2, 1),
		Initial_Electron_State("Na", 2, 0),
		Initial_Electron_State("Na", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestNormalizationSilicon)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Si", 3, 1),
		Initial_Electron_State("Si", 3, 0),
		Initial_Electron_State("Si", 2, 1),
		Initial_Electron_State("Si", 2, 0),
		Initial_Electron_State("Si", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestNormalizationXenon)
{
	// ARRANGE
	double tol									  = 1.0e-5;
	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 1),
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Xe", 4, 2),
		Initial_Electron_State("Xe", 4, 1),
		Initial_Electron_State("Xe", 4, 0),
		Initial_Electron_State("Xe", 3, 2),
		Initial_Electron_State("Xe", 3, 1),
		Initial_Electron_State("Xe", 3, 0),
		Initial_Electron_State("Xe", 2, 1),
		Initial_Electron_State("Xe", 2, 0),
		Initial_Electron_State("Xe", 1, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		EXPECT_NEAR(electron.Normalization(), 1.0, tol);
}

TEST(TestWavefunctions, TestPrintSummary)
{
	// ARRANGE
	std::string element = "Xe";
	int n				= 5;
	int l				= 0;
	// ACT
	Initial_Electron_State Xe5p(element, n, l);
	// ASSERT
	Xe5p.Print_Summary();
}