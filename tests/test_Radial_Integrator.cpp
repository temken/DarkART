#include "gtest/gtest.h"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Radial_Integrator.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestRadialIntegrator, TestIntegrator)
{
	// ARRANGE
	double tol	   = 1e-6;
	double k_final = keV;
	double q	   = 0.1 * keV;
	int l_final	   = 5;
	int L		   = 4;

	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Ar", 3, 0),
	};
	Final_Electron_State_Hydrogenic final_state(1.0);

	// ACT & ASSERT
	for(auto& electron : electrons)
		for(int integral = 1; integral < 3; integral++)
		{
			Radial_Integrator integrator(electron, final_state);
			double integral_adaptive = integrator.Radial_Integral(integral, k_final, q, l_final, L);

			integrator.Use_Tabulated_Functions(5000, {k_final}, {q});
			double integral_tabulated = integrator.Radial_Integral(integral, k_final, q, l_final, L);

			EXPECT_NEAR(integral_adaptive, integral_tabulated, std::fabs(tol * integral_adaptive));
		}
}

TEST(TestRadialIntegrator, TestSetNewStates)
{
	// ARRANGE
	double tol	   = 1e-6;
	double k_final = keV;
	double q	   = 0.1 * keV;
	int l_final	   = 5;
	int L		   = 4;

	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Ar", 3, 0),
	};
	Final_Electron_State_Hydrogenic final_state(1.0);
	// ACT
	Radial_Integrator integrator;

	// ASSERT
	for(auto& electron : electrons)
		for(int integral = 1; integral < 3; integral++)
		{
			integrator.Set_New_States(electron, final_state);
			EXPECT_EQ(integrator.initial_state.Orbital_Name(), electron.Orbital_Name());

			double integral_adaptive = integrator.Radial_Integral(integral, k_final, q, l_final, L);
			integrator.Use_Tabulated_Functions(5000, {k_final}, {q});
			double integral_tabulated = integrator.Radial_Integral(integral, k_final, q, l_final, L);
			EXPECT_NEAR(integral_adaptive, integral_tabulated, std::fabs(tol * integral_adaptive));
		}
}