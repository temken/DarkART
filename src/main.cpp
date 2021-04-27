#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>

#include "libphysica/Natural_Units.hpp"

#include "Atomic_Responses.hpp"
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
	// Xenon_5p.Print_Summary();
	// std::cout << Xenon_5p.Radial_Wavefunction(Bohr_Radius) << std::endl;
	// std::cout << Xenon_5p.Radial_Wavefunction_Derivative(Bohr_Radius) << std::endl;
	// std::cout << Xenon_5p.Normalization() << std::endl;

	// double k_final = keV;
	// double q	   = keV;
	// int l_final	   = 1;
	// double Z_eff   = Xenon_5p.Z_eff;
	// double r	   = Bohr_Radius;
	// std::cout << Radial_Wavefunction_Final(k_final, l_final, Z_eff, r) << std::endl;

	double l = 10.0, Z = 10.0, k = 100 * keV, r = 100 * Bohr_Radius;

	std::cout << Radial_Wavefunction_Final(k, l, Z, r) << std::endl;

	// std::complex<double> a = l + 1.0 + 1.0i * Z / k / Bohr_Radius;
	// std::complex<double> b = 2.0 * l + 2.0;
	// std::complex<double> z = 2.0i * k * r;
	// std::cout << a << "\t" << b << "\t" << z << std::endl;
	// std::cout << Hypergeometric_1F1_asymptotic(a, b, z) << std::endl;
	// std::cout << Hypergeometric_1F1_series(a, b, z) << std::endl;

	// if(status)
	// 	std::cout << "error" << std::endl;
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