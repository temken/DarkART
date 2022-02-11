#include "DarkART/Special_Functions.hpp"

#include <algorithm>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_coupling.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "arb_hypgeom.h"

namespace DarkART
{
// using namespace std::complex_literals;

// 1. Gaunt coefficients
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

}	// namespace DarkART