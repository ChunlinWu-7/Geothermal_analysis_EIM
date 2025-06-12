
#include "eyemat.h"

// constants
const double pi = 2.0 * asin(1.0);
const double Pi = 2.0 * asin(1.0);
const double E = exp(1.0);

/*
	-1: semi-infinite with infinite-conductivity
	0: full-space solution
	1: semi-infinite with zero-conductivity
*/
const double full_flag = -1.0;

/*
	0: daily flag;
	1: annual flag;
*/
const int annual_flag = 1;

// thermal constants for annual temperature purposes:
double K_0 = 3.62;\
double Cp_0 = 2.787E6;\
double alpha_0 = K_0 / Cp_0;\

// annual change:
const double omega = 1.0 / (364.0 * 24 * 3600.0) * 2.0 * pi;

const double sdd = 6.0;

// excitation frequency: 
double f_m = sqrt(omega / (2.0 * alpha_0));

complex<double> FF = f_m * (1.0 + 1.0i);

// Kornecter delta function
const double d[3][3] = { 1,0,0,0,1,0,0,0,1 };

// imaged matrix
const double M[3] = { 1, 1, -1 };

// define 16-point Gauss quadrature
const int ngp = 16;\
const double gp[16] = { -0.09501251,0.09501251,-0.281603551,0.281603551,-0.458016778,0.458016778,-0.617876244,0.617876244,-0.755404408,0.755404408,\
-0.865631202,0.865631202,-0.944575023,0.944575023,-0.989400935,0.989400935 };\
const double w[16] = { 0.18945061,0.18945061,0.182603415,0.182603415,0.169156519,0.169156519,0.149595989,0.149595989,0.124628971,\
0.124628971,0.095158512,0.095158512,0.062253524,0.062253524,0.027152459,0.027152459 };