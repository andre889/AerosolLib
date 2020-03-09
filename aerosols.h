#pragma once
#include <math.h>
#include <cmath>
#include "constants.h"
#include <iostream>
namespace aerosols
{

	template <class T>
	T getMFP(const T& temp, const T& P)
	{
		//based on Willeke's (1976) Temperature Dependence of Particle Slip in a Gaseous Medium. Relations valid from 200k to 1000k
		//https://doi.org/10.1016/0021-8502(76)90024-0

#ifdef DEBUG_
		if (temp > 1000 || temp < 200)
		{
			std::cout << "CRITICAL ERROR: TEMP OUT OF BOUNDS FOR THIS MFP CALCULATION Willeke(1976)/n";
			exit(2);
		}
#endif
		double mfp0 = 67.3e-9;
		double T0 = 296.15;
		double p0 = 101325.0;

		return mfp0* (temp / T0)* (p0 / P)* (1.0 + 110.4 / T0) / (1.0 + 110.4 / temp);
	}

	//calculate slip correction fror particle size and mfp
	template <class T>
	T getCc(const T& dp, const T& mfp)
	{	
#define _KIM_ //sellected Model
#ifdef _KIM_ //from Kim(2005) PSL with uncertainty if interested. http://dx.doi.org/10.6028/jres.110.005.
#define a 1.165
#define b 0.486
#define c -0.997
#endif 
#ifdef _BASE_ //unkown auther, given by Chris Hogan
#define a 1.257
#define b 0.4
#define c -1.1
#endif
#ifdef _Rader_ //Rader 1990 oil drops
#define a 1.209
#define b 0.441
#define c -0.779
#endif
#ifdef _Hutchins_ //Hutchins 1995 PSL
#define a 1.231
#define b 0.469
#define c -1.178
#endif

		T kn = 2.0 * mfp / dp;
		T Cc = 1.0 + kn * (a + b * exp(c / kn));
		return Cc;
	}

	// Calculates drag coeficient based on Morse Alexander (1972) 
	// for 0 < Re < 50,000
	// https://doi.org/10.1017/S0022112072001806
	template <class T>
	T getCD(const T& Re)
	{
		T Cd = 0.0;

		if (Re < 0.1)
		{
			Cd = 24.0 / Re;
		}
		else if (Re < 1.0)
		{
			Cd = 22.73 / Re + 0.0903 / (Re * Re) + 3.69;
		}
		else if (Re < 10.0)
		{
			Cd = 29.1667 / Re - 3.8889 / (Re * Re) + 1.222;
		}
		else if (Re < 100.0)
		{
			Cd = 46.5 / Re - 116.67 / (Re * Re) + 0.6167;
		}
		else if (Re < 1000.0)
		{
			Cd = 98.33 / Re - 2778.0 / (Re * Re) + 0.3644;
		}
		else if (Re < 5000.0)
		{
			Cd = 148.62 / Re - 47500.0 / (Re * Re) + 0.357;
		}
		else if (Re < 10000.0)
		{
			Cd = -490.546 / Re + 578700.0 / (Re * Re) + 0.46;
		}
		else if (Re < 50000.0)
		{
			Cd = -1662.5 / Re - 5416700.0 / (Re * Re) + 0.5191;
		}
		else
		{
			Cd = 0.4;
		}

		return Cd;
	}

	//calculate terminal velocity assuming rhoP >> rho with slip
	template <class T>
	T termVelCalc(const T& rhoP, const T& dp, const T& mu, const T& g, const T& rho,const T& mfp)
	{
		double w = 0.3;

		T Cc = getCc(dp, mfp);

		T Vterm = rhoP * dp * dp * g * Cc / 18.0 / mu;

		T Vterm_new = 0.0;
		//Cc = 1 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
		T error = 10.0;

		while (error > 1e-6)
		{
			error = 0;

			T Re = rho * dp * Vterm / mu;

			T Cd_star = getCD(Re);
			Cc = getCc(dp, mfp);
			T Cd = Cd_star / Cc;

			Vterm_new = sqrt(dp*4.0*g*rhoP/(3.0*rho*Cd));

			error = abs(Vterm_new - Vterm) / Vterm;

			Vterm = Vterm_new*(w) + Vterm * (1.0 - w);
		}
	
			return Vterm;

	}

	//calculate diameter assuming rhoP >> rho with slip
	template <class T>
	T diatermVelCalc(const T& rhoP, const T& Vt, const T& mu, const T& g, const T& rho, const T& mfp)
	{
		double w = 0.1;
		T Cc = 1.0;

		//guess
		T dp = sqrt(18.0*mu*Vt/rhoP/g);
		std::cout << dp << std::endl;
		T dp_new = 0.0;
		//Cc = 1 + kn * (1.257 + 0.4 * exp(-1.1 / kn));
		T error = 10.0;

		while (error > 1e-6)
		{
			error = 0;

			T Re = rho * dp * Vt / mu;

			T Cd_star = getCD(Re);
			Cc = getCc(dp, mfp);
			T Cd = Cd_star / Cc;

			dp_new = (3.0*rho*Cd*Vt*Vt)/ (4.0 * g * rhoP );

			error = abs(dp_new - dp) / dp;

			dp = dp_new*(w) + dp*(1.0 - w);
			//std::cout << dp << std::endl;
		}

		return dp;

	}


	


	//get friction factor with slip
	template <class T>
	T getff(const T& dp, const T& mu, const T& mfp)
	{
		T Cc = getCc(dp, mfp);
		T ff = 3.0 * constants::PI * mu * dp / Cc;
		return ff;

	}

	//calculate diffusion coeff baised on stokes einstien assuming f = 3pidpmu/cc
	template <class T>
	T getDiff(const T& temp, const T& dp, const T& mu, const T& mfp)
	{
		T ff = getff(dp, mu, mfp);

		T D = constants::K * temp / ff;

		return D;
	}
	
};
