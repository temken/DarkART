#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "DarkART/Atomic_Responses.hpp"
#include "DarkART/Configuration.hpp"
#include "DarkART/Radial_Integrator.hpp"
#include "DarkART/Response_Tabulation.hpp"
#include "DarkART/Special_Functions.hpp"
#include "DarkART/Wavefunctions_Final.hpp"
#include "DarkART/Wavefunctions_Initial.hpp"
#include "DarkART/version.hpp"

using namespace DarkART;
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
	if(argc < 2)
	{
		std::cerr << "\nError: No configuration file provided." << std::endl
				  << "\tCorrect usage:\t" << argv[0] << " <configuration file>" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	DarkART::Configuration cfg(argv[1]);
	cfg.Print_Summary();

	if(cfg.run_modus == "Tabulation")
	{

		std::cout << "\nTabulate atomic responses for " << cfg.element << std::endl
				  << std::endl;

		Response_Tabulator tabulator(cfg.k_min, cfg.k_max, cfg.q_min, cfg.q_max);
		tabulator.Resize_Grid(cfg.k_gridpoints, cfg.q_gridpoints);

		Radial_Integrator radial_integrator;
		if(cfg.tabulate_radial_functions)
		{
			auto k_grid = libphysica::Log_Space(cfg.k_min, cfg.k_max, cfg.k_gridpoints);
			auto q_grid = libphysica::Log_Space(cfg.q_min, cfg.q_max, cfg.q_gridpoints);
			radial_integrator.Use_Tabulated_Functions(cfg.r_gridpoints, k_grid, q_grid);
		}

		int counter		  = 1;
		int num_responses = cfg.atomic_responses.size() * cfg.atomic_shell_list.size();
		for(auto& atomic_shell_name : cfg.atomic_shell_list)
		{
			Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
			Final_Electron_State_Hydrogenic final_state(initial_state.Z_eff);

			if(!cfg.overwrite_old_tables)
			{
				bool atomic_shell_finished = true;
				for(auto& response : cfg.atomic_responses)
					atomic_shell_finished *= libphysica::File_Exists(cfg.results_path + initial_state.Orbital_Name() + "_" + std::to_string(response) + "_Table.txt");
				if(atomic_shell_finished)
				{
					counter += cfg.atomic_responses.size();
					std::cout << "\t" << initial_state.Orbital_Name() << ": All responses were already tabulated." << std::endl;
					continue;
				}
			}

			radial_integrator.Set_New_States(initial_state, final_state);

			for(auto& response : cfg.atomic_responses)
			{
				counter++;
				if(!cfg.overwrite_old_tables && libphysica::File_Exists(cfg.results_path + initial_state.Orbital_Name() + "_" + std::to_string(response) + "_Table.txt"))
					std::cout << "\t" << initial_state.Orbital_Name() << ": Response " << response << " was already tabulated." << std::endl;
				else
				{
					std::cout << counter - 1 << " / " << num_responses << ":\t";
					tabulator.Tabulate(response, radial_integrator, cfg.threads);
					tabulator.Export_Tables(cfg.results_path);
					std::cout << std::endl
							  << SEPARATOR << std::endl;
				}
			}
		}
	}
	else if(cfg.run_modus == "Evaluation")
	{
		std::cout << "Evaluate atomic responses for k' = " << cfg.k_prime / keV << " and q = " << cfg.q / keV << " keV\n"
				  << std::endl;
		for(auto& atomic_shell_name : cfg.atomic_shell_list)
		{
			Initial_Electron_State initial_state(cfg.element, atomic_shell_name);
			Final_Electron_State_Hydrogenic final_state(initial_state.Z_eff);
			Radial_Integrator radial_integrator(initial_state, final_state);

			if(cfg.tabulate_radial_functions)
				radial_integrator.Use_Tabulated_Functions(cfg.r_gridpoints, {cfg.k_prime}, {cfg.q});

			for(auto& response : cfg.atomic_responses)
			{
				int l_convergence;
				double W = Atomic_Response_Function(response, cfg.k_prime, cfg.q, radial_integrator, l_convergence);
				std::cout << "\t" << initial_state.Orbital_Name() << "\tW_" << response << "(k',q) = " << W << "\t(l' ≤ " << l_convergence << ")" << std::endl;
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