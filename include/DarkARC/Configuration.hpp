#ifndef __Configuration_hpp__
#define __Configuration_hpp__

#include <vector>

#include "libphysica/Utilities.hpp"

namespace DarkARC
{

class Configuration : public libphysica::Configuration
{
  protected:
	virtual void Initialize_Parameters() override;

  public:
	std::string run_modus, element;
	std::vector<std::string> atomic_shell_list;
	std::vector<int> atomic_responses;
	bool overwrite_old_tables;
	double k_min, k_max, q_min, q_max, k_prime, q;
	int k_gridpoints, q_gridpoints, threads;

	bool tabulate_radial_functions;
	int r_gridpoints;

	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	virtual void Print_Summary(int MPI_rank = 0) override;
};

}	// namespace DarkARC

#endif
