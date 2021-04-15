#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

// Headers from libphysica
#include "Natural_Units.hpp"

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
	Xenon_5p.Print_Summary();
	std::cout << Xenon_5p.Radial_Wavefunction(Bohr_Radius) << std::endl;
	std::cout << Xenon_5p.Radial_Wavefunction_Derivative(Bohr_Radius) << std::endl;
	std::cout << Xenon_5p.Normalization() << std::endl;

	double k_final = keV;
	int l_final	   = 1;
	double Z_eff   = 1.0;
	double r	   = Bohr_Radius;
	std::cout << Radial_Wavefunction_Final(k_final, l_final, Z_eff, r) << std::endl;
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