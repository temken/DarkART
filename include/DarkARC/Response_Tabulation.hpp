#ifndef __Response_Tabulation_hpp_
#define __Response_Tabulation_hpp_

#include "Wavefunctions.hpp"

namespace DarkARC
{

class Response_Tabulator
{
  private:
	double k_min, k_max, q_min, q_max;
	std::vector<double> k_grid;
	std::vector<double> q_grid;
	std::vector<std::vector<double>> response_table;
	std::string electron_orbital;
	int tabulated_response;

	void Initialize_Lists(int k_points, int q_points);

  public:
	Response_Tabulator(double kmin, double kmax, double qmin, double qmax);

	void Resize_Grid(int k_points, int q_points = 0);

	void Tabulate(int response, const Initial_Electron_State& bound_electron, int threads = 1);

	void Export_Tables(const std::string& path);
};

}	// namespace DarkARC

#endif