#ifndef __Configuration_hpp__
#define __Configuration_hpp__

#include <libconfig.h++>
#include <vector>

namespace DarkARC
{

class Configuration
{
  protected:
	libconfig::Config config;
	std::string cfg_file;

	void Read_Config_File();

	void Initialize_Result_Folder(int MPI_rank = 0);
	void Create_Result_Folder(int MPI_rank = 0);
	void Copy_Config_File(int MPI_rank = 0);

	void Initialize_Parameters();

  public:
	std::string ID, results_path, run_modus, element;
	std::vector<std::string> atomic_shells;
	std::vector<int> responses;
	double k_min, k_max, q_min, q_max, k_prime, q;
	int k_gridpoints, q_gridpoints, threads;

	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	void Print_Summary(int MPI_rank = 0);
};

}	// namespace DarkARC

#endif
