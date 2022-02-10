#include "DarkARC/Special_Functions.hpp"

#include <algorithm>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_coupling.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "arb_hypgeom.h"

namespace DarkARC
{
using namespace std::complex_literals;

// 1. Spherical Harmonics
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

double Gaunt_Coefficient(int j1, int j2, int j3, int m1, int m2, int m3)
{
	return sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0)) / sqrt(4.0 * M_PI) * gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0) * gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3);
}

// 2. Spherical Bessel function j_L
double Spherical_Bessel_jL_arb(int L, double x)
{
	double prefactor = sqrt(M_PI / 2.0 / x);
	double result;
	slong prec;
	arb_t J_arb, x_arb, n_arb;
	arb_init(J_arb);
	arb_init(x_arb);
	arb_init(n_arb);
	arb_set_d(x_arb, x);
	arb_set_d(n_arb, 1.0 * L + 0.5);
	for(prec = 80;; prec *= 2)
	{
		arb_hypgeom_bessel_j(J_arb, n_arb, x_arb, prec);
		if(arb_rel_accuracy_bits(J_arb) >= 53)
		{
			result = arf_get_d(arb_midref(J_arb), ARF_RND_NEAR);
			break;
		}
		else if(prec > 10000)
		{
			result = NAN;
			break;
		}
	}
	arb_clear(J_arb);
	arb_clear(n_arb);
	arb_clear(x_arb);
	return prefactor * result;
}

double Spherical_Bessel_jL(int L, double x)
{
	if(L < 10)
		return gsl_sf_bessel_jl(L, x);
	else
		return Spherical_Bessel_jL_arb(L, x);
}

// 3. Coulomb wave
double Coulomb_Wave_ARB(int L, double eta, double rho)
{
	double result;
	slong prec;
	arb_t F, l, eta_2, rho_2;
	arb_init(F);
	arb_init(l);
	arb_init(eta_2);
	arb_init(rho_2);
	arb_set_d(l, L);
	arb_set_d(eta_2, eta);
	arb_set_d(rho_2, rho);
	for(prec = 80;; prec *= 2)
	{
		arb_hypgeom_coulomb(F, NULL, l, eta_2, rho_2, prec);
		if(arb_rel_accuracy_bits(F) >= 53)
		{
			result = arf_get_d(arb_midref(F), ARF_RND_NEAR);
			break;
		}
		else if(prec > 10000)
		{
			result = NAN;
			break;
		}
	}
	arb_clear(F);
	arb_clear(eta_2);
	arb_clear(rho_2);
	arb_clear(l);
	return result;
}

double Coulomb_Wave_GSL(int L, double eta, double rho, int& status)
{
	double fc_array[1];
	double F_exponent[1];
	status = gsl_sf_coulomb_wave_F_array(L, 0, eta, rho, fc_array, F_exponent);
	return fc_array[0];
}

double Coulomb_Wave(int L, double eta, double rho)
{
	int status;
	double cw = Coulomb_Wave_GSL(L, eta, rho, status);
	if(status != 0 || std::isnan(cw))
		cw = Coulomb_Wave_ARB(L, eta, rho);
	return cw;
}

std::vector<std::vector<double>> Compute_Gauss_Legendre_Roots_and_Weights(unsigned int n, double x_min, double x_max)
{
	std::vector<std::vector<double>> roots_and_weights(n, std::vector<double>(2, 0.0));

	double eps			= 1.0e-14;
	int m				= (n + 1) / 2;
	double x_middle		= 0.5 * (x_max + x_min);
	double x_half_width = 0.5 * (x_max - x_min);

	for(int i = 0; i < m; i++)
	{
		double pp;
		double z = cos(M_PI * (i + 0.75) / (n + 0.5));
		while(true)
		{
			double p1 = 1.0;
			double p2 = 0.0;
			for(int j = 0; j < n; j++)
			{
				double p3 = p2;
				p2		  = p1;
				p1		  = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
			}
			pp		  = n * (z * p1 - p2) / (z * z - 1.0);
			double z1 = z;
			z		  = z1 - p1 / pp;
			if(std::fabs(z - z1) <= eps)
				break;
		}
		roots_and_weights[i][0]			= x_middle - x_half_width * z;
		roots_and_weights[n - i - 1][0] = x_middle + x_half_width * z;
		roots_and_weights[i][1]			= 2.0 * x_half_width / ((1.0 - z * z) * pp * pp);
		roots_and_weights[n - i - 1][1] = roots_and_weights[i][1];
	}
	return roots_and_weights;
}

unsigned int Locate_Closest_Location(const std::vector<double>& list, double target)
{
	auto const it = std::upper_bound(list.begin(), list.end(), target);
	if(it == list.end())
		return list.size() - 1;
	else
	{
		int index = std::distance(list.begin(), it);
		if(index == 0)
			return 0;
		else if(index == list.size())
			return list.size() - 1;
		else
		{
			double diff1 = std::fabs(list[index - 1] - target);
			double diff2 = std::fabs(list[index] - target);
			if(diff1 < diff2)
				return index - 1;
			else
				return index;
		}
	}
}

}	// namespace DarkARC