#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "DarkARC/Configuration.hpp"
#include "version.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestConfiguration, TestConfiguration)
{
	std::string ID				 = "test_cfg";
	std::string run_modus		 = "Tabulation";
	std::string element			 = "Ar";
	int atomic_shell_list_length = 5;
	int atomic_responses_length	 = 4;
	bool overwrite				 = true;
	int threads					 = 4;
	int k_gridpoints			 = 250;
	int q_gridpoints			 = 150;
	double kmin					 = 0.1 * keV;
	double k_max				 = 100 * keV;
	double q_min				 = 1 * keV;
	double q_max				 = 1000.0 * keV;

	double k_prime = 10.0 * keV;
	double q	   = 10 * keV;

	bool tabulate_radial_functions = true;
	int r_gridpoints			   = 5000;

	// ARRANGE
	std::string filename = "test.cfg";
	// ACT
	Configuration config(filename);

	// ASSERT
	EXPECT_EQ(config.ID, ID);
	EXPECT_EQ(config.run_modus, run_modus);
	EXPECT_EQ(config.element, element);
	EXPECT_EQ(config.atomic_shell_list.size(), atomic_shell_list_length);
	EXPECT_EQ(config.atomic_responses.size(), atomic_responses_length);
	EXPECT_EQ(config.overwrite_old_tables, overwrite);
	EXPECT_EQ(config.threads, threads);
	EXPECT_EQ(config.k_gridpoints, k_gridpoints);
	EXPECT_EQ(config.q_gridpoints, q_gridpoints);
	EXPECT_DOUBLE_EQ(config.k_min, kmin);
	EXPECT_DOUBLE_EQ(config.k_max, k_max);
	EXPECT_DOUBLE_EQ(config.q_min, q_min);
	EXPECT_DOUBLE_EQ(config.q_max, q_max);
	EXPECT_DOUBLE_EQ(config.k_prime, k_prime);
	EXPECT_DOUBLE_EQ(config.q, q);
	EXPECT_EQ(config.tabulate_radial_functions, tabulate_radial_functions);
	EXPECT_EQ(config.r_gridpoints, r_gridpoints);
}

TEST(TestConfiguration, TestPrintSummary)
{
	// ARRANGE
	std::string filename = "test.cfg";

	Configuration config(filename);
	// ACT & ASSERT
	config.Print_Summary();
}