#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkARC/Atomic_Responses.hpp"
#include "DarkARC/Configuration.hpp"
#include "DarkARC/Radial_Integrator.hpp"
#include "DarkARC/Response_Tabulation.hpp"
#include "DarkARC/Special_Functions.hpp"
#include "DarkARC/Wavefunctions_Final.hpp"
#include "DarkARC/Wavefunctions_Initial.hpp"
#include "version.hpp"

using namespace DarkARC;
using namespace libphysica::natural_units;

int main(int argc, char* argv[])
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

	DarkARC::Configuration cfg(argv[1]);
	cfg.Print_Summary();

	if(cfg.run_modus == "Tabulation")
	{
		auto k_grid = libphysica::Log_Space(cfg.k_min, cfg.k_max, cfg.k_gridpoints);
		auto q_grid = libphysica::Log_Space(cfg.q_min, cfg.q_max, cfg.q_gridpoints);

		Response_Tabulator tabulator(cfg.k_min, cfg.k_max, cfg.q_min, cfg.q_max);
		tabulator.Resize_Grid(cfg.k_gridpoints, cfg.q_gridpoints);

		int counter		  = 1;
		int num_responses = cfg.atomic_responses.size() * cfg.atomic_shell_list.size();
		for(auto& response : cfg.atomic_responses)
			for(auto& atomic_shell_name : cfg.atomic_shell_list)
			{
				Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
				Final_Electron_State* final_state = new Final_Electron_State_Hydrogenic(initial_state.Z_eff);
				std::cout << counter++ << "/" << num_responses << ")" << std::endl;
				if(!cfg.overwrite_old_tables && libphysica::File_Exists(cfg.results_path + initial_state.Orbital_Name() + "_" + std::to_string(response) + "_Table.txt"))
					std::cout << "\tResponse " << response << " of " << initial_state.Orbital_Name() << " was already tabulated.\n\tTo re-calculate this response, remove the corresponding files from the /results/ folder." << std::endl;
				else
				{
					tabulator.Tabulate(response, initial_state, *final_state, cfg.threads);
					tabulator.Export_Tables(cfg.results_path);
				}
			}
	}
	else if(cfg.run_modus == "Evaluation")
	{
		std::cout << "Evaluate atomic responses for k' = " << cfg.k_prime / keV << " and q = " << cfg.q / keV << " keV" << std::endl;
		for(auto& response : cfg.atomic_responses)
		{
			std::cout << std::endl;
			for(auto& atomic_shell_name : cfg.atomic_shell_list)
			{
				Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
				Final_Electron_State_Hydrogenic final_state(initial_state.Z_eff);

				int l_convergence;
				double W = Atomic_Response_Function(response, cfg.k_prime, cfg.q, initial_state, final_state, l_convergence);
				std::cout << "\t" << initial_state.Orbital_Name() << "\tW_" << response << "(k',q) = " << W << "\t(maximum l' = " << l_convergence << ")" << std::endl;
			}
		}
	}
	else
	{
		Initial_Electron_State Xe_5p("Xe", 5, 1);
		Final_Electron_State_Hydrogenic hydrogenic_state(Xe_5p.Z_eff);
		Radial_Integrator integrator(Xe_5p, hydrogenic_state);

		auto k_grid = libphysica::Log_Space(0.1 * keV, 100 * keV, 100);	  // cfg.k_gridpoints);
		auto q_grid = libphysica::Log_Space(keV, 1000 * keV, 100);		  // cfg.q_gridpoints);

		unsigned int l_final = 2;
		unsigned int L		 = 3;
		integrator.Use_Tabulated_Functions(10000, k_grid, q_grid);
		// std::cout << integrator.Radial_Integral_Adaptive(1, k_grid[40], q_grid[40], l_final, L) << std::endl;
		std::cout << integrator.Radial_Integral(1, k_grid[40], q_grid[40], l_final, L) << std::endl;
		// for(int l_final = 0; l_final <= 150; l_final++)
		// {
		// 	std::cout << std::endl
		// 			  << l_final << "\t" << integrator.Radial_Integral_Adaptive(1, k_grid[0], q_grid[0], l_final, l_final) << "\t" << std::flush;
		// 	std::cout << integrator.Radial_Integral(1, k_grid[0], q_grid[0], l_final, l_final) << "\t" << std::flush;
		// 	std::cout << "\t" << integrator.Radial_Integral_Adaptive(1, k_grid[1], q_grid[1], l_final, l_final) << "\t" << std::flush;
		// 	std::cout << integrator.Radial_Integral(1, k_grid[1], q_grid[1], l_final, l_final) << std::endl;
		// }
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]" << std::endl;

	return 0;
}