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

		for(auto& atomic_shell_name : cfg.atomic_shell_list)
		{
			Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
			Final_Electron_State_Hydrogenic final_state(initial_state.Z_eff);

			// Check if tables exist already
			if(!cfg.overwrite_old_tables)
			{
				bool all_tables_exist = true;
				for(auto& response : cfg.atomic_responses)
				{
					bool table_exists = libphysica::File_Exists(cfg.results_path + initial_state.Orbital_Name() + "_" + std::to_string(response) + "_Table.txt");
					all_tables_exist *= table_exists;
					if(table_exists)
						counter++;
				}
				if(all_tables_exist)
				{
					std::cout << "\nResponses for " << initial_state.Orbital_Name() << " were already tabulated." << std::endl;
					continue;
				}
			}

			Radial_Integrator radial_integrator(initial_state, final_state);
			if(cfg.tabulate_radial_functions)
				radial_integrator.Use_Tabulated_Functions(cfg.r_gridpoints, k_grid, q_grid);

			for(auto& response : cfg.atomic_responses)
			{
				std::cout << std::endl
						  << counter++ << "/" << num_responses << ")" << std::endl;
				tabulator.Tabulate(response, radial_integrator, cfg.threads);
				tabulator.Export_Tables(cfg.results_path);
			}
		}
	}
	else if(cfg.run_modus == "Evaluation")
	{
		std::cout << "Evaluate atomic responses for k' = " << cfg.k_prime / keV << " and q = " << cfg.q / keV << " keV" << std::endl;
		for(auto& atomic_shell_name : cfg.atomic_shell_list)
		{
			std::cout << std::endl;

			Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
			Final_Electron_State_Hydrogenic final_state(initial_state.Z_eff);
			Radial_Integrator radial_integrator(initial_state, final_state);

			if(cfg.tabulate_radial_functions)
				radial_integrator.Use_Tabulated_Functions(cfg.r_gridpoints, {cfg.k_prime}, {cfg.q});

			for(auto& response : cfg.atomic_responses)
			{
				int l_convergence;
				double W = Atomic_Response_Function(response, cfg.k_prime, cfg.q, radial_integrator, l_convergence);
				std::cout << "\t" << initial_state.Orbital_Name() << "\tW_" << response << "(k',q) = " << W << "\t(maximum l' = " << l_convergence << ")" << std::endl;
			}
		}
	}
	else
	{
	}

	////////////////////////////////////////////////////////////////////////
	// Final terminal output
	auto time_end		 = std::chrono::system_clock::now();
	double durationTotal = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
	std::cout << "\n[Finished in " << libphysica::Time_Display(durationTotal) << "]" << std::endl;

	return 0;
}