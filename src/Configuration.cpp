#include "DarkARC/Configuration.hpp"

#include <fstream>
#include <iostream>
#include <sys/stat.h>	 //required to create a folder
#include <sys/types.h>	 // required for stat.h

#include "version.hpp"

namespace DarkARC
{

using namespace libconfig;

void Configuration::Read_Config_File()
{
	try
	{
		config.readFile(cfg_file.c_str());
	}
	catch(const FileIOException& fioex)
	{
		std::cerr << "Error in DarkARC::Configuration::Read_Config_File(): I/O error while reading configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	catch(const ParseException& pex)
	{
		std::cerr << "Error in DarkARC::Configuration::Read_Config_File(): Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Configuration::Initialize_Result_Folder(int MPI_rank)
{
	try
	{
		ID = config.lookup("ID").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in DarkARC::Configuration::Initialize_Result_Folder(): No 'ID' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	results_path = TOP_LEVEL_DIR "results/" + ID + "/";
	Create_Result_Folder(MPI_rank);
	Copy_Config_File(MPI_rank);
}

void Configuration::Create_Result_Folder(int MPI_rank)
{
	if(MPI_rank == 0)
	{
		// 1. Create the /results/ folder if necessary
		std::string results_folder = TOP_LEVEL_DIR "results";
		mode_t nMode			   = 0733;	 // UNIX style permissions
		int nError_1			   = 0;
#if defined(_WIN32)
		nError_1 = _mkdir(results_folder.c_str());	 // can be used on Windows
#else
		nError_1 = mkdir(results_folder.c_str(), nMode);   // can be used on non-Windows
#endif

		// 2. Create a /result/<ID>/ folder for result files.
		int nError = 0;
#if defined(_WIN32)
		nError = _mkdir(results_path.c_str());	 // can be used on Windows
#else
		nError	 = mkdir(results_path.c_str(), nMode);	   // can be used on non-Windows
#endif
		if(nError != 0)
		{
			std::cerr << "\nWarning in Configuration::Create_Result_Folder(int): The folder exists already, data will be overwritten." << std::endl
					  << std::endl;
		}
	}
}

void Configuration::Copy_Config_File(int MPI_rank)
{
	if(MPI_rank == 0)
	{
		std::ifstream inFile;
		std::ofstream outFile;
		inFile.open(cfg_file);
		outFile.open(TOP_LEVEL_DIR "results/" + ID + "/" + ID + ".cfg");
		outFile << "// " << PROJECT_NAME << "-v" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl;
		outFile << inFile.rdbuf();
		inFile.close();
		outFile.close();
	}
}

void Configuration::Initialize_Parameters()
{
	// try
	// {
	// 	constraints_certainty = config.lookup("constraints_certainty");
	// }
	// catch(const SettingNotFoundException& nfex)
	// {
	// 	std::cerr << "Error in Configuration::Initialize_Parameters(): No 'constraints_certainty' setting in configuration file." << std::endl;
	// 	std::exit(EXIT_FAILURE);
	// }
}

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
: cfg_file(cfg_filename), results_path("./")
{

	// 1. Read the cfg file.
	Read_Config_File();

	// 2. Find the run ID, create a folder and copy the cfg file.
	Initialize_Result_Folder(MPI_rank);

	// 3. Set parameters
	Initialize_Parameters();
}

void Configuration::Print_Summary(int MPI_rank)
{
	if(MPI_rank == 0)
	{
	}
}

}	// namespace DarkARC