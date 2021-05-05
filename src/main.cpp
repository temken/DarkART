#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "Atomic_Responses.hpp"
#include "Response_Tabulation.hpp"
#include "Special_Functions.hpp"
#include "Wavefunctions.hpp"
#include "version.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

int main()
{
	//Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	std::cout << "[Started on " << ctime_start << "]" << std::endl;
	std::cout << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
			  << LOGO << std::endl;
	////////////////////////////////////////////////////////////////////////

	Initial_Electron_State Xenon_5p("Xe", 5, 1);
	int response = 1;

	double q_min = 1.0 * keV;
	double q_max = 1000 * keV;
	double k_min = 0.1 * keV;
	double k_max = 100.0 * keV;

	Response_Tabulator tabulator(q_min, q_max, k_min, k_max);
	tabulator.Resize_Grid(100);
	tabulator.Tabulate(response, Xenon_5p);
	tabulator.Export_Tables(TOP_LEVEL_DIR "results/");

	// for(auto& q : qs)
	// {
	// 	double W = Atomic_Response_Function(k, q, Xenon_5p, index);
	// 	// std::cout << q / keV << "\t" << W << std::endl;
	// }

	////////////////////////////////////////////////////////////////////////
	//Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << std::round(1000. * durationTotal) / 1000. << "s";
	if(durationTotal > 60.0)
		std::cout << " (" << floor(durationTotal / 3600.0) << ":" << floor(fmod(durationTotal / 60.0, 60.0)) << ":" << floor(fmod(durationTotal, 60.0)) << ")]." << std::endl;
	else
		std::cout << "]" << std::endl;

	return 0;
}