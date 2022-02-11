#ifndef __Response_Tabulation_hpp_
#define __Response_Tabulation_hpp_

#include "Radial_Integrator.hpp"
#include "Wavefunctions_Final.hpp"
#include "Wavefunctions_Initial.hpp"

namespace DarkART
{

class Response_Tabulator
{
  private:
	double k_min, k_max, q_min, q_max;
	std::vector<double> k_grid;
	std::vector<double> q_grid;
	std::vector<std::vector<double>> response_table;
	std::vector<std::vector<int>> l_prime_table;
	std::string electron_orbital;
	int tabulated_response;

	void Initialize_Lists(int k_points, int q_points);

  public:
	Response_Tabulator(double kmin, double kmax, double qmin, double qmax);

	void Resize_Grid(int k_points, int q_points = 0);

	void Tabulate(int response, Radial_Integrator& radial_integrator, int threads = 1);
	void Tabulate(int response, const Initial_Electron_State& bound_electron, Final_Electron_State& final_state, bool use_tables = false, int threads = 1);

	void Export_Tables(const std::string& path);
};

}	// namespace DarkART

#endif