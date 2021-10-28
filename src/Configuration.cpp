#include "DarkARC/Configuration.hpp"

#include <fstream>
#include <iostream>
#include <sys/stat.h>	 //required to create a folder
#include <sys/types.h>	 // required for stat.h

#include "libphysica/Natural_Units.hpp"

#include "version.hpp"

namespace DarkARC
{

using namespace libconfig;
using namespace libphysica::natural_units;

void Configuration::Initialize_Parameters()
{
	try
	{
		run_modus = config.lookup("run_modus").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): No 'run_modus' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		element = config.lookup("element").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): No 'element' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		int atomic_shell_count = config.lookup("atomic_shells").getLength();
		for(int j = 0; j < atomic_shell_count; j++)
			atomic_shell_list.push_back(config.lookup("atomic_shells")[j]);
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'atomic_shells' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		int atomic_response_count = config.lookup("atomic_responses").getLength();
		for(int j = 0; j < atomic_response_count; j++)
			atomic_responses.push_back(config.lookup("atomic_responses")[j]);
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "No 'atomic_shells' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		threads = config.lookup("threads");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): No 'threads' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		overwrite_old_tables = config.lookup("overwrite_old_tables");
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): No 'overwrite_old_tables' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		k_min		 = config.lookup("k_min");
		k_max		 = config.lookup("k_max");
		k_gridpoints = config.lookup("k_points");

		k_prime = config.lookup("k_prime");

		k_min *= keV;
		k_max *= keV;
		k_prime *= keV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): Faulty k' grid setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		q_min		 = config.lookup("q_min");
		q_max		 = config.lookup("q_max");
		q_gridpoints = config.lookup("q_points");

		q = config.lookup("q");

		q_min *= keV;
		q_max *= keV;
		q *= keV;
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in Configuration::Initialize_Parameters(): Faulty q grid setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
: libphysica::Configuration(cfg_filename, MPI_rank)
{
	Initialize_Parameters();
}

void Configuration::Print_Summary(int MPI_rank)
{
	if(MPI_rank == 0)
	{
		std::cout << SEPARATOR << std::endl
				  << "DarkARC++ configuration summary" << std::endl
				  << "ID:\t" << ID << std::endl
				  << "\nElement:\t" << element << std::endl
				  << "Atomic shells:\t";
		for(unsigned int i = 0; i < atomic_shell_list.size(); i++)
			std::cout << atomic_shell_list[i] << ((i < atomic_shell_list.size() - 1) ? ", " : "\n");
		std::cout << "\nRun modus:\t" << run_modus << std::endl;
		if(run_modus != "Evaluation")
			std::cout << "Threads:\t" << threads << std::endl
					  << "k' grid [keV]:\t" << k_min / keV << " - " << k_max / keV << "\t(" << k_gridpoints << " points)" << std::endl
					  << "q grid [keV]:\t" << q_min / keV << " - " << q_max / keV << "\t(" << q_gridpoints << " points)" << std::endl;
		if(run_modus != "Tabulation")
			std::cout << "k' [keV]:\t" << k_prime / keV << std::endl
					  << "q [keV]:\t" << q / keV << std::endl;
		std::cout << SEPARATOR << std::endl;
	}
}

}	// namespace DarkARC