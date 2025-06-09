#pragma once

#include "stdlib.h"
#include <complex>
#include <cmath>

using namespace std;

// constants
extern const double pi, Pi, E;

// flags on infinite, semi-infinite problems!
extern const double full_flag; 

extern const int annual_flag;

// excitation frequency
extern const double omega, f_m;
extern const complex<double> FF;

// heat conductivity and capacity, heat transfer coefficient
extern const double K_0, Cp_0, alpha_0;

// Kornecter delta function
extern const double d[3][3];

// imaged matrix
extern const double M[3];

// define 16 - point Gaussian integrals

extern const int ngp;
extern const double gp[16];
extern const double w[16];
