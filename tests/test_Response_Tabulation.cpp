#include "gtest/gtest.h"

#include "Response_Tabulation.hpp"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "Wavefunctions.hpp"
#include "version.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

TEST(TestTabulator, TestTabulator)
{
	// ARRANGE
	int gridsize = 3;
	double k_min = 0.1 * keV;
	double k_max = 1.0 * keV;
	double q_min = 1.0 * keV;
	double q_max = 2.0 * keV;
	auto k_grid	 = libphysica::Log_Space(k_min, k_max, gridsize);
	auto q_grid	 = libphysica::Log_Space(q_min, q_max, gridsize);
	Initial_Electron_State Xenon_5s("Xe", 5, 0);
	std::string filepath_1 = TOP_LEVEL_DIR "tests/Xe_5s_1_Table.txt";
	std::string filepath_2 = TOP_LEVEL_DIR "tests/Xe_5s_1_List.txt";
	// ACT
	Response_Tabulator tabulator(k_min, k_max, q_min, q_max);
	tabulator.Resize_Grid(gridsize, gridsize);
	tabulator.Tabulate(1, Xenon_5s, 1);
	tabulator.Export_Tables(TOP_LEVEL_DIR "tests/");
	// ASSERT
	EXPECT_TRUE(libphysica::File_Exists(filepath_1));
	// std::vector<std::vector<double>> table = libphysica::Import_Table(filepath_1, {}, 3);
	// EXPECT_EQ(table.size(), gridsize);
	// for(int i = 0; i < gridsize; i++)
	// 	EXPECT_EQ(table[i].size(), gridsize);
	// for(int i = 0; i < gridsize; i++)
	// 	for(int j = 0; j < gridsize; j++)
	// 		EXPECT_GT(table[i][j], 0.0);

	EXPECT_TRUE(libphysica::File_Exists(filepath_2));
	// std::vector<std::vector<double>> list = libphysica::Import_Table(filepath_2, {}, 2);
	// EXPECT_EQ(list.size(), 9);
	// for(int i = 0; i < gridsize; i++)
	// 	EXPECT_EQ(list[i].size(), gridsize);
	// for(int i = 0; i < gridsize * gridsize; i++)
	// {
	// 	EXPECT_NEAR(list[i][0] * keV, k_grid[i / 3], 1.0e-6);
	// 	EXPECT_NEAR(list[i][1] * keV, q_grid[i % 3], 1.0e-6);
	// 	EXPECT_GT(list[i][2], 0.0);
	// }
}