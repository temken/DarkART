#include "gtest/gtest.h"

#include "DarkARC/Atomic_Responses.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Wavefunctions.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestAtomicResponse, TestDipoleApproximation)
{
	// ARRANGE
	double k_final = keV;
	double q_1	   = 0.1 * keV;
	double q_2	   = 0.2 * keV;

	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Ar", 3, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
	{
		double W_1 = Atomic_Response_Function(k_final, q_1, electron, 1);
		double W_2 = Atomic_Response_Function(k_final, q_2, electron, 1);
		double tol = 1.0e-2 * std::fabs(W_2);
		EXPECT_NEAR(W_2, q_2 / q_1 * q_2 / q_1 * W_1, tol);
	}
}

TEST(TestAtomicResponse, TestResponses)
{
	// ARRANGE
	double k_final = keV;
	double q	   = 0.1 * keV;

	std::vector<Initial_Electron_State> electrons = {
		Initial_Electron_State("Xe", 5, 0),
		Initial_Electron_State("Ar", 3, 0),
	};
	// ACT & ASSERT
	for(auto& electron : electrons)
		for(int response = 1; response < 5; response++)
			if(response == 2)
				EXPECT_NE(Atomic_Response_Function(k_final, q, electron, response), 0.0);
			else
				EXPECT_GT(Atomic_Response_Function(k_final, q, electron, response), 0.0);
}