#include "Spherical_Harmonics.hpp"

#include <iostream>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace DarkARC
{
using namespace std::complex_literals;

std::complex<double> Spherical_Harmonics(int l, int m, double theta, double phi)
{
	return boost::math::spherical_harmonic(l, m, theta, phi);
}

std::complex<double> VSH_Y_Component(int component, int l, int m, int l_hat, int m_hat)
{
	switch(component)
	{
		case 0:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 && m_hat == m + 1)
				return -1.0 / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 && m_hat == m - 1)
				return 1.0 / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m + 1)
				return 1.0 / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m - 1)
				return -1.0 / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 1:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 && m_hat == m + 1)
				return 1.0i / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 && m_hat == m - 1)
				return 1.0i / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 2:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m))
				return 0.0;
			else if(l_hat == l + 1)
				return 1.0 * std::sqrt(1.0 * (l - m + 1) * (l + m + 1) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1)
				return 1.0 * std::sqrt(1.0 * (l - m) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		default:
			std::cerr << "Error in VSH_Y_Component(): Component " << component << " out of bound." << std::endl;
			std::exit(EXIT_FAILURE);
	}
	return 0.0;
}

std::vector<std::complex<double>> Vector_Spherical_Harmonics_Y(int l, int m, double theta, double phi)
{
	std::vector<std::complex<double>> Y(3, 0.0);
	for(unsigned int i = 0; i < 3; i++)
		for(int l_hat = l - 1; l_hat < l + 2; l_hat += 2)
			for(int m_hat = m - 1; m_hat < m + 2; m_hat++)
				if(std::abs(m_hat) <= l_hat)
					Y[i] += VSH_Y_Component(i, l, m, l_hat, m_hat) * Spherical_Harmonics(l_hat, m_hat, theta, phi);
	return Y;
}

std::complex<double> VSH_Psi_Component(int component, int l, int m, int l_hat, int m_hat)
{
	switch(component)
	{
		case 0:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 and m_hat == m + 1)
				return l / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 and m_hat == m - 1)
				return -l / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m + 1)
				return (l + 1) / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m - 1)
				return -(l + 1) / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 1:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 and m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * l * l * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 and m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * l * l * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + 1) * (l + 1) * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + 1) * (l + 1) * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 2:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m))
				return 0.0;
			else if(l_hat == l + 1)
				return -l * std::sqrt(1.0 * (l - m + 1) * (l + m + 1) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1)
				return (1 + l) * std::sqrt(1.0 * (l - m) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		default:
			std::cerr << "Error in VSH_Psi_Component(): Component " << component << " out of bound." << std::endl;
			std::exit(EXIT_FAILURE);
	}
	return 0.0;
}

std::vector<std::complex<double>> Vector_Spherical_Harmonics_Psi(int l, int m, double theta, double phi)
{
	std::vector<std::complex<double>> Psi(3, 0.0);
	for(unsigned int i = 0; i < 3; i++)
		for(int l_hat = l - 1; l_hat < l + 2; l_hat += 2)
			for(int m_hat = m - 1; m_hat < m + 2; m_hat++)
				if(std::abs(m_hat) <= l_hat)
					Psi[i] += VSH_Psi_Component(i, l, m, l_hat, m_hat) * Spherical_Harmonics(l_hat, m_hat, theta, phi);
	return Psi;
}

}	// namespace DarkARC