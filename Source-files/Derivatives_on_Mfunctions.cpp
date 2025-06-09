
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

template <typename R>
promote<R> DIE::M0_real(const R& r)
{
	if (r > a) {
		return (pi * (cosh(a * f_m) * (-(cos(f_m * r) * sin(a * f_m))
			+ (-2 * a * f_m * cos(a * f_m) + sin(a * f_m)) * sin(f_m * r))
			+ (2 * a * f_m * cos(f_m * r) * sin(a * f_m) + cos(a * f_m) * (cos(f_m * r) + sin(f_m * r))) * sinh(a * f_m))
			) / (pow(E, f_m * r) * pow(f_m, 3) * r);
	}
	else {
		return (
			pi * cosh(f_m * r) * (-((1 + 2 * a * f_m) * cos(a * f_m)) + sin(a * f_m)) * sin(f_m * r) 
			+ pi * cos(f_m * r) * (cos(a * f_m) + (1 + 2 * a * f_m) * sin(a * f_m)) * sinh(f_m * r)
			) /
			(pow(E, a * f_m) * pow(f_m, 3) * r);

	}
}

template <typename R>
promote<R> DIE::M0_imag(const R& r)
{
	if (r > a) {
		return (pi * (cosh(a * f_m) * (2 * a * f_m * cos(a * f_m) * cos(f_m * r) - sin(a * f_m) * (cos(f_m * r) + sin(f_m * r)))
			+ (2 * a * f_m * sin(a * f_m) * sin(f_m * r) + cos(a * f_m) * (-cos(f_m * r) + sin(f_m * r))) * sinh(a * f_m))
			) / (pow(E, f_m * r) * pow(f_m, 3) * r);
	}
	else {
		return (pi * (2 * pow(E, a * f_m) * f_m * r - cosh(f_m * r) * (cos(a * f_m) + (1 + 2 * a * f_m) * sin(a * f_m)) * sin(f_m * r) +
			cos(f_m * r) * (-((1 + 2 * a * f_m) * cos(a * f_m)) + sin(a * f_m)) * sinh(f_m * r))) / (pow(E, a * f_m) * pow(f_m, 3) * r);
	}
}

template <typename R>
promote<R> DIE::M1_real(const R& r)
{
	if (r > a) {
		return (pi * (cosh(a * f_m) * (sin(a * f_m) * (-(f_m * (2 * pow(a, 2) * f_m + r) * cos(f_m * r)) + (3 + f_m * r) * sin(f_m * r)) - a * f_m * cos(a * f_m) * (3 * cos(f_m * r) + (3 + 2 * f_m * r) * sin(f_m * r))) +
			(a * f_m * sin(a * f_m) * ((3 + 2 * f_m * r) * cos(f_m * r) - 3 * sin(f_m * r)) + cos(a * f_m) * ((3 + f_m * r) * cos(f_m * r) + f_m * (2 * pow(a, 2) * f_m + r) * sin(f_m * r))) * sinh(a * f_m))) / (pow(E, f_m * r) * pow(f_m, 4) * r);
	}
	else {
		return (pi * (-2 * pow(E, a * f_m) * f_m * r + cosh(f_m * r) * (-(f_m * r * cos(f_m * r) * (cos(a * f_m)
			+ (1 + 2 * a * f_m) * sin(a * f_m))) + (-(a * f_m * (3 + 2 * a * f_m) * cos(a * f_m))
				+ 3 * (1 + a * f_m) * sin(a * f_m)) * sin(f_m * r)) +
			(cos(f_m * r) * (3 * (1 + a * f_m) * cos(a * f_m) + a * f_m * (3 + 2 * a * f_m) * sin(a * f_m))
				+ f_m * r * ((1 + 2 * a * f_m) * cos(a * f_m) - sin(a * f_m)) * sin(f_m * r)) * sinh(f_m * r))
			) / (pow(E, a * f_m) * pow(f_m, 4) * r);

	}
}

template <typename R>
promote<R> DIE::M1_imag(const R& r)
{
	if (r > a) {
		return (pi * cosh(a * f_m) * (a * f_m * cos(a * f_m) * ((3 + 2 * f_m * r) * cos(f_m * r) - 3 * sin(f_m * r)) + sin(a * f_m) * (-((3 + f_m * r) * cos(f_m * r)) - f_m * (2 * pow(a, 2) * f_m + r) * sin(f_m * r))) +
			pi * (cos(a * f_m) * (-(f_m * (2 * pow(a, 2) * f_m + r) * cos(f_m * r)) + (3 + f_m * r) * sin(f_m * r)) + a * f_m * sin(a * f_m) * (3 * cos(f_m * r) + (3 + 2 * f_m * r) * sin(f_m * r))) * sinh(a * f_m)) / (pow(E, f_m * r) * pow(f_m, 4) * r);
	}
	else {
		return (pi * (2 * pow(E, a * f_m) * f_m * r + cosh(f_m * r) * (cos(a * f_m) * (f_m * (1 + 2 * a * f_m) * r * cos(f_m * r)
			- 3 * (1 + a * f_m) * sin(f_m * r)) - f_m * sin(a * f_m) * (r * cos(f_m * r) + a * (3 + 2 * a * f_m) * sin(f_m * r))) +
			(cos(f_m * r) * (-(a * f_m * (3 + 2 * a * f_m) * cos(a * f_m))
				+ 3 * (1 + a * f_m) * sin(a * f_m)) + f_m * r * (cos(a * f_m)
					+ (1 + 2 * a * f_m) * sin(a * f_m)) * sin(f_m * r)) * sinh(f_m * r))
			) / (pow(E, a * f_m) * pow(f_m, 4) * r);

	}
}

template <typename R>
promote<R> DIE::M2_real(const R& r)
{
	if (r > a) {
		return (pi * (-(cosh(a * f_m) * (sin(a * f_m) * ((-6 + pow(f_m, 2) * (pow(r, 2) + pow(a, 2) * (5 + 4 * f_m * r))) * cos(f_m * r) - (6 + 5 * pow(a, 2) * pow(f_m, 2) + f_m * r * (6 + f_m * r)) * sin(f_m * r)) +
			2 * a * f_m * cos(a * f_m) * (3 * (2 + f_m * r) * cos(f_m * r) + f_m * (3 * r + f_m * (pow(a, 2) + pow(r, 2))) * sin(f_m * r)))) +
			(2 * a * f_m * sin(a * f_m) * (f_m * (3 * r + f_m * (pow(a, 2) + pow(r, 2))) * cos(f_m * r) - 3 * (2 + f_m * r) * sin(f_m * r)) +
				cos(a * f_m) * ((6 + 5 * pow(a, 2) * pow(f_m, 2) + f_m * r * (6 + f_m * r)) * cos(f_m * r) + (-6 + pow(f_m, 2) * (pow(r, 2) + pow(a, 2) * (5 + 4 * f_m * r))) * sin(f_m * r))) * sinh(a * f_m))) / (pow(E, f_m * r) * pow(f_m, 5) * r);
	}
	else {
		return (pi * (-6 * pow(E, a * f_m) * f_m * r + cosh(f_m * r) * (2 * f_m * r * cos(f_m * r) * (-3 * (1 + a * f_m) * cos(a * f_m) - a * f_m * (3 + 2 * a * f_m) * sin(a * f_m)) +
			(-((-6 + pow(f_m, 2) * (pow(a, 2) * (5 + 2 * a * f_m) + (1 + 2 * a * f_m) * pow(r, 2))) * cos(a * f_m)) + (6 + f_m * (a * (12 + 5 * a * f_m) + f_m * pow(r, 2))) * sin(a * f_m)) * sin(f_m * r)) +
			(cos(f_m * r) * ((6 + f_m * (a * (12 + 5 * a * f_m) + f_m * pow(r, 2))) * cos(a * f_m) + (-6 + pow(f_m, 2) * (pow(a, 2) * (5 + 2 * a * f_m) + (1 + 2 * a * f_m) * pow(r, 2))) * sin(a * f_m)) +
				2 * f_m * r * (a * f_m * (3 + 2 * a * f_m) * cos(a * f_m) - 3 * (1 + a * f_m) * sin(a * f_m)) * sin(f_m * r)) * sinh(f_m * r))) / (pow(E, a * f_m) * pow(f_m, 5) * r);
	}
}

template <typename R>
promote<R> DIE::M2_imag(const R& r)
{
	if (r > a) {
		return (pi * (cosh(a * f_m) * (2 * a * f_m * cos(a * f_m) * (f_m * (pow(a, 2) * f_m + r * (3 + f_m * r)) * cos(f_m * r) - 3 * (2 + f_m * r) * sin(f_m * r)) -
			sin(a * f_m) * ((6 + 5 * pow(a, 2) * pow(f_m, 2) + f_m * r * (6 + f_m * r)) * cos(f_m * r) + (-6 + pow(f_m, 2) * (pow(r, 2) + pow(a, 2) * (5 + 4 * f_m * r))) * sin(f_m * r))) +
			(cos(a * f_m) * (-((-6 + pow(f_m, 2) * (pow(r, 2) + pow(a, 2) * (5 + 4 * f_m * r))) * cos(f_m * r)) + (6 + 5 * pow(a, 2) * pow(f_m, 2) + f_m * r * (6 + f_m * r)) * sin(f_m * r)) +
				2 * a * f_m * sin(a * f_m) * (3 * (2 + f_m * r) * cos(f_m * r) + f_m * (3 * r + f_m * (pow(a, 2) + pow(r, 2))) * sin(f_m * r))) * sinh(a * f_m))) / (pow(E, f_m * r) * pow(f_m, 5) * r);
	}
	else {
		return -((pi * (cosh(f_m * r) * (2 * f_m * r * cos(f_m * r) * (-(a * f_m * (3 + 2 * a * f_m) * cos(a * f_m)) + 3 * (1 + a * f_m) * sin(a * f_m)) +
			((6 + f_m * (a * (12 + 5 * a * f_m) + f_m * pow(r, 2))) * cos(a * f_m) + (-6 + pow(f_m, 2) * (pow(a, 2) * (5 + 2 * a * f_m) + (1 + 2 * a * f_m) * pow(r, 2))) * sin(a * f_m)) * sin(f_m * r)) +
			(cos(f_m * r) * ((-6 + pow(f_m, 2) * (pow(a, 2) * (5 + 2 * a * f_m) + (1 + 2 * a * f_m) * pow(r, 2))) * cos(a * f_m) - (6 + f_m * (a * (12 + 5 * a * f_m) + f_m * pow(r, 2))) * sin(a * f_m)) +
				2 * f_m * r * (-3 * (1 + a * f_m) * cos(a * f_m) - a * f_m * (3 + 2 * a * f_m) * sin(a * f_m)) * sin(f_m * r)) * sinh(f_m * r))) / (pow(E, a * f_m) * pow(f_m, 5) * r));

	}
}

// before calculation: must assign value of beta
void DIE::Store_derivs(int nsolve, double* x, double* x_p)
{
	if (nsolve == 4) {

		// uniform case: only require M_0 and second order
		M0_r = new complex<double>[3];
		M0_r_ori = new complex<double>[3];

		// imaged terms
		auto var = make_ftuple<float_50, 2>(dist);
		auto const result_real = M0_real(get<0>(var));
		auto const result_imag = M0_imag(get<0>(var));

		for (int h = 0; h < 3; h++) {
			M0_r[h] = double(result_real.derivative(h)) + 1.0i * double(result_imag.derivative(h));
		}

		/*
			Differ cases with original terms!
		*/
		if (dist_flag != 0) {
			auto var_ori = make_ftuple<float_50, 2>(dist_ori);
			auto const result_real_ori = M0_real(get<0>(var_ori));
			auto const result_imag_ori = M0_imag(get<0>(var_ori));

			for (int h = 0; h < 3; h++) {
				M0_r_ori[h] = double(result_real_ori.derivative(h)) + 1.0i * double(result_imag_ori.derivative(h));
			}
		}

	}

	else if (nsolve == 16) {
		// linear case: (i) require M0 to third order; (ii) require M1 up to fourth order
		M0_r = new complex<double>[4];  M1_r = new complex<double>[5];
		M0_r_ori = new complex<double>[4];  M1_r_ori = new complex<double>[5];

		// imaged terms!
		auto var = make_ftuple<float_50, 4>(dist);

		auto const result_real = M0_real(get<0>(var));
		auto const result_imag = M0_imag(get<0>(var));

		for (int h = 0; h < 4; h++) {
			M0_r[h] = double(result_real.derivative(h)) + 1.0i * double(result_imag.derivative(h));
		}

		auto const result_real1 = M1_real(get<0>(var));
		auto const result_imag1 = M1_imag(get<0>(var));

		for (int h = 0; h < 5; h++) {
			M1_r[h] = double(result_real1.derivative(h)) + 1.0i * double(result_imag1.derivative(h));
		}

		/*
			Differ cases with original terms!
		*/
		if (dist_flag != 0) {
			auto var_ori = make_ftuple<float_50, 4>(dist_ori);
			auto const result_real_ori = M0_real(get<0>(var_ori));
			auto const result_imag_ori = M0_imag(get<0>(var_ori));

			for (int h = 0; h < 4; h++) {
				M0_r_ori[h] = double(result_real_ori.derivative(h)) + 1.0i * double(result_imag_ori.derivative(h));
			}

			auto const result_real1_ori = M1_real(get<0>(var_ori));
			auto const result_imag1_ori = M1_imag(get<0>(var_ori));

			for (int h = 0; h < 5; h++) {
				M1_r_ori[h] = double(result_real1_ori.derivative(h)) + 1.0i * double(result_imag1_ori.derivative(h));
			}

		}

	}
	else if (nsolve == 52) {
		// quadratic case: (i) require M0 to fourth order; (ii) require M1 up to fifth order
		// (iii) require M2 to sixth order
		M0_r = new complex<double>[5];  M1_r = new complex<double>[6];
		M2_r = new complex<double>[7];

		M0_r_ori = new complex<double>[5];  M1_r_ori = new complex<double>[6];
		M2_r_ori = new complex<double>[7];

		// imaged terms!
		auto var = make_ftuple<float_50, 6>(dist);
		auto const result_real = M0_real(get<0>(var));
		auto const result_imag = M0_imag(get<0>(var));

		for (int h = 0; h < 5; h++) {
			M0_r[h] = double(result_real.derivative(h)) + 1.0i * double(result_imag.derivative(h));
		}

		auto const result_real1 = M1_real(get<0>(var));
		auto const result_imag1 = M1_imag(get<0>(var));

		for (int h = 0; h < 6; h++) {
			M1_r[h] = double(result_real1.derivative(h)) + 1.0i * double(result_imag1.derivative(h));
		}

		auto const result_real2 = M2_real(get<0>(var));
		auto const result_imag2 = M2_imag(get<0>(var));

		for (int h = 0; h < 7; h++) {
			M2_r[h] = double(result_real2.derivative(h)) + 1.0i * double(result_imag2.derivative(h));
		}

		if (dist_flag != 0) {
			auto var_ori = make_ftuple<float_50, 6>(dist_ori);
			auto const result_real_ori = M0_real(get<0>(var_ori));
			auto const result_imag_ori = M0_imag(get<0>(var_ori));
			for (int h = 0; h < 5; h++) {
				M0_r_ori[h] = double(result_real_ori.derivative(h)) + 1.0i * double(result_imag_ori.derivative(h));
			}

			auto const result_real1_ori = M1_real(get<0>(var_ori));
			auto const result_imag1_ori = M1_imag(get<0>(var_ori));

			for (int h = 0; h < 6; h++) {
				M1_r_ori[h] = double(result_real1_ori.derivative(h)) + 1.0i * double(result_imag1_ori.derivative(h));
			}

			auto const result_real2_ori = M2_real(get<0>(var_ori));
			auto const result_imag2_ori = M2_imag(get<0>(var_ori));

			for (int h = 0; h < 7; h++) {
				M2_r_ori[h] = double(result_real2_ori.derivative(h)) + 1.0i * double(result_imag2_ori.derivative(h));
			}

		}


	}
}