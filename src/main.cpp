#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkARC/Atomic_Responses.hpp"
#include "DarkARC/Response_Tabulation.hpp"
#include "DarkARC/Special_Functions.hpp"
#include "DarkARC/Wavefunctions.hpp"
#include "version.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

int main()
{
	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << LOGO << std::endl;
	////////////////////////////////////////////////////////////////////////

	// Initial_Electron_State Xenon_5p("Xe", 5, 0);
	// double k_final = 100.0 * keV;
	// // double q	   = 100.0 * keV;
	// double q_min = 1.0 * keV;
	// double q_max = 1000 * keV;
	// auto q_grid	 = libphysica::Log_Space(q_min, q_max, 100);

	// std::ofstream f;
	// f.open("test_5.txt");
	// for(auto& q : q_grid)
	// 	f << Atomic_Response_Function(k_final, q, Xenon_5p, 1) << std::endl;
	// f.close();
	// std::cout << Atomic_Response_Function(k_final, q, Xenon_5p, 1) << std::endl;

	// // Input
	// int num_threads								  = 4;
	// std::vector<int> responses					  = {1, 2, 3, 4};
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
		// Initial_Electron_State("Ar", 3, 1),
		// Initial_Electron_State("Ar", 3, 0),
		// Initial_Electron_State("Ar", 2, 1),
		// Initial_Electron_State("Ar", 2, 0),
		// Initial_Electron_State("Ar", 1, 0),
	};
	double r	= 1.0 * Bohr_Radius;
	auto r_list = libphysica::Log_Space(1e-3 * Bohr_Radius, 1e2 * Bohr_Radius, 100);
	for(auto r : r_list)
	{
		double sum = 0.0;
		for(auto& electron : electrons)
			sum += 2.0 * (2.0 * electron.l + 1) * electron.Radial_Integral(r);
		double Z_eff = 54.0 - sum + electrons.back().Radial_Integral(r);
		std::cout << r / Bohr_Radius << "\t" << Z_eff << std::endl;
	}
	// double k_min = 0.1 * keV;
	// double k_max = 100.0 * keV;
	// double q_min = 1.0 * keV;
	// double q_max = 1000 * keV;
	// int k_points = 100;
	// int q_points = 200;

	// auto k_grid = libphysica::Log_Space(k_min, k_max, k_points);
	// auto q_grid = libphysica::Log_Space(q_min, q_max, q_points);

	// Response_Tabulator tabulator(k_min, k_max, q_min, q_max);
	// tabulator.Resize_Grid(k_points, q_points);
	// int counter		  = 1;
	// int num_responses = responses.size() * electrons.size();
	// for(auto& electron : electrons)
	// 	for(auto& response : responses)
	// 	{
	// 		std::cout << counter++ << "/" << num_responses << ")" << std::endl;
	// 		if(libphysica::File_Exists(TOP_LEVEL_DIR "results/" + electron.Orbital_Name() + "_" + std::to_string(response) + "_Table.txt"))
	// 			std::cout << "\tResponse " << response << " of " << electron.Orbital_Name() << " was already tabulated.\n\tTo re-calculate this response, remove the corresponding files from the /results/ folder." << std::endl;
	// 		else
	// 		{
	// 			tabulator.Tabulate(response, electron, num_threads);
	// 			tabulator.Export_Tables(TOP_LEVEL_DIR "results/");
	// 		}
	// 	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
	if(durationTotal > 60.0)
		std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
	else
		std::cout << "]" << std::endl;

	return 0;
}