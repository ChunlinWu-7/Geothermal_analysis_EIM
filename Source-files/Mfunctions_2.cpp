
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

/*
	M2 functions: original
*/

complex<double> DIE::M2_ori()
{
	if (dist_flag != 0) {
		return M2_r_ori[0];
	}
	else {
		return 4.0 * pi / pow(beta, 4.0) * (
			6.0 + exp(1.0i * a * beta) * (
				-6.0 + a * beta * (
					6.0i + a * beta * (3.0 - 1.0i * a * beta)
					)
				)
			);
	}
}

complex<double> DIE::M2_ori_1(int i)
{
	if (dist_flag != 0) {
		return M2_r_ori[1] * psi_ori_1(i);
	}
	else {
		return 0.0 + 0.0i;
	}
}

complex<double> DIE::M2_ori_2(int i, int j)
{
	if (dist_flag != 0) {
		return M2_r_ori[1] * psi_ori_2(i, j) + M2_r_ori[2] * psi_ori_1(i) * psi_ori_1(j);
	}
	else {
		return d[i][j] * 4.0 / 3.0 * pi * a * a * exp(1.0i * a * beta) * (1.0 + 1.0i * a * beta);
	}
}

complex<double> DIE::M2_ori_3(int i, int j, int k)
{
	if (dist_flag != 0) {
		return M2_r_ori[1] * psi_ori_3(i, j, k)
			+ M2_r_ori[2] * (
				psi_ori_1(k) * psi_ori_2(i, j)
				+ psi_ori_1(i) * psi_ori_2(j, k)
				+ psi_ori_1(j) * psi_ori_2(i, k)
				)
			+ M2_r_ori[3] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k);
	}
	else {
		return 0.0 + 0.0i;
	}
}

complex<double> DIE::M2_ori_4(int i, int j, int k, int l)
{
	if (dist_flag != 0) {
		return M2_r_ori[1] * psi_ori_4(i, j, k, l)
			+ M2_r_ori[2] * (
				psi_ori_2(k, l) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, l)
				+ psi_ori_2(i, l) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, l)
				+ psi_ori_2(j, l) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, l)
				+ psi_ori_3(i, j, k) * psi_ori_1(l)
				)
			+ M2_r_ori[3] * (
				psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, l)
				+ psi_ori_1(k) * psi_ori_2(i, j) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_2(j, k) * psi_ori_1(l)
				+ psi_ori_1(j) * psi_ori_2(i, k) * psi_ori_1(l)
				)
			+ M2_r_ori[4] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l);
	}
	else {
		return (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k]) * (
			-4.0i / 15.0 * exp(1.0i * a * beta) * pi *
			(-1.0i + a * beta) * (
				2.0 + a * beta * (-4.0i + a * beta)
				)
			);
	}
}

complex<double> DIE::M2_ori_5(int i, int j, int k, int l, int m)
{
	if (dist_flag != 0) {
		return M2_r_ori[1] * psi_ori_5(i, j, k, l, m)
			+ M2_r_ori[2] * psi_ori_4(i, j, k, l) * psi_ori_1(m)

			+ M2_r_ori[2] * (psi_ori_4(i, j, k, m) * psi_ori_1(l) + psi_ori_3(i, j, k) * psi_ori_2(l, m))
			+ M2_r_ori[3] * psi_ori_1(m) * (psi_ori_3(i, j, k) * psi_ori_1(l))


			+ (M2_r_ori[4] * psi_ori_1(l) * psi_ori_1(m) + M2_r_ori[3] * psi_ori_2(l, m)) * (
				psi_ori_1(k) * psi_ori_2(i, j)
				+ psi_ori_1(i) * psi_ori_2(j, k)
				+ psi_ori_1(j) * psi_ori_2(i, k)
				)

			+ M2_r_ori[3] * psi_ori_1(l) * (
				psi_ori_2(k, m) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, m)
				+ psi_ori_2(i, m) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, m)
				+ psi_ori_2(j, m) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, m)
				)

			+ M2_r_ori[3] * psi_ori_1(m) * (
				psi_ori_2(k, l) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, l)
				+ psi_ori_2(i, l) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, l)
				+ psi_ori_2(j, l) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, l)
				)

			+ M2_r_ori[2] * (
				psi_ori_3(k, l, m) * psi_ori_2(i, j) + psi_ori_2(k, l) * psi_ori_3(i, j, m)
				+ psi_ori_2(k, m) * psi_ori_3(i, j, l) + psi_ori_1(k) * psi_ori_4(i, j, l, m)
				+ psi_ori_3(i, l, m) * psi_ori_2(j, k) + psi_ori_2(i, l) * psi_ori_3(j, k, m)
				+ psi_ori_2(i, m) * psi_ori_3(j, k, l) + psi_ori_1(i) * psi_ori_4(j, k, l, m)
				+ psi_ori_3(j, l, m) * psi_ori_2(i, k) + psi_ori_2(j, l) * psi_ori_3(i, k, m)
				+ psi_ori_2(j, m) * psi_ori_3(i, k, l) + psi_ori_1(j) * psi_ori_4(i, k, l, m)
				)

			+ M2_r_ori[4] * psi_ori_1(m) * (
				psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, l)
				)

			+ M2_r_ori[3] * (
				psi_ori_3(i, l, m) * psi_ori_1(j) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_2(j, m) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_2(k, m)
				+ psi_ori_2(i, m) * psi_ori_2(j, l) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_3(j, l, m) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_2(k, m)
				+ psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_3(k, l, m)
				)

			+ M2_r_ori[4] * (
				psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_1(k) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, m) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_2(l, m)
				)

			+ M2_r_ori[5] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l) * psi_ori_1(m);
	}
	else {
		return 0.0 + 0.0i;
	}
}

complex<double> DIE::M2_ori_6(int i, int j, int k, int l, int m, int n)
{
	if (dist_flag != 0) {
		complex<double> result = M2_r_ori[1] * psi_ori_6(i, j, k, l, m, n) + M2_r_ori[2] * psi_ori_1(n) * psi_ori_5(i, j, k, l, m)
			+ M2_r_ori[2] * psi_ori_5(i, j, k, l, n) * psi_ori_1(m) + M2_r_ori[2] * psi_ori_4(i, j, k, l) * psi_ori_2(m, n)
			+ M2_r_ori[3] * psi_ori_4(i, j, k, l) * psi_ori_1(m) * psi_ori_1(n);

		result += M2_r_ori[3] * psi_ori_1(n) * (psi_ori_4(i, j, k, m) * psi_ori_1(l) + psi_ori_3(i, j, k) * psi_ori_2(l, m))
			+ M2_r_ori[2] * (
				psi_ori_5(i, j, k, m, n) * psi_ori_1(l) + psi_ori_4(i, j, k, m) * psi_ori_2(l, n)
				+ psi_ori_4(i, j, k, n) * psi_ori_2(l, m) + psi_ori_3(i, j, k) * psi_ori_3(l, m, n)
				)
			+ M2_r_ori[4] * psi_ori_1(n) * psi_ori_1(m) * psi_ori_3(i, j, k) * psi_ori_1(l)
			+ M2_r_ori[3] * (
				psi_ori_2(m, n) * psi_ori_3(i, j, k) * psi_ori_1(l)
				+ psi_ori_1(m) * psi_ori_4(i, j, k, n) * psi_ori_1(l)
				+ psi_ori_1(m) * psi_ori_3(i, j, k) * psi_ori_2(l, n)
				);


		result += (M2_r_ori[5] * psi_ori_1(l) * psi_ori_1(m) * psi_ori_1(n) + M2_r_ori[4] * psi_ori_2(l, n) * psi_ori_1(m)
			+ M2_r_ori[4] * psi_ori_1(l) * psi_ori_2(m, n) + M2_r_ori[4] * psi_ori_1(n) * psi_ori_2(l, m) + M2_r_ori[3] * psi_ori_3(l, m, n))
			* (
				psi_ori_1(k) * psi_ori_2(i, j)
				+ psi_ori_1(i) * psi_ori_2(j, k)
				+ psi_ori_1(j) * psi_ori_2(i, k)
				)
			+ (M2_r_ori[4] * psi_ori_1(l) * psi_ori_1(m) + M2_r_ori[3] * psi_ori_2(l, m)) * (
				psi_ori_2(k, n) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, n)
				+ psi_ori_2(i, n) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, n)
				+ psi_ori_2(j, n) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, n)
				);

		result += (M2_r_ori[4] * psi_ori_1(l) * psi_ori_1(n) + M2_r_ori[3] * psi_ori_2(l, n)) * (
			psi_ori_2(k, m) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, m)
			+ psi_ori_2(i, m) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, m)
			+ psi_ori_2(j, m) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, m)
			) + M2_r_ori[3] * psi_ori_1(l) * (
				psi_ori_3(k, m, n) * psi_ori_2(i, j) + psi_ori_2(k, n) * psi_ori_3(i, j, m)
				+ psi_ori_2(k, m) * psi_ori_3(i, j, n) + psi_ori_1(k) * psi_ori_4(i, j, m, n)
				+ psi_ori_3(i, m, n) * psi_ori_2(j, k) + psi_ori_2(i, n) * psi_ori_3(j, k, m)
				+ psi_ori_2(i, m) * psi_ori_3(j, k, n) + psi_ori_1(i) * psi_ori_4(j, k, m, n)
				+ psi_ori_3(j, m, n) * psi_ori_2(i, k) + psi_ori_2(j, n) * psi_ori_3(i, k, m)
				+ psi_ori_2(j, m) * psi_ori_3(i, k, n) + psi_ori_1(j) * psi_ori_4(i, k, m, n)
				);

		result += (M2_r_ori[4] * psi_ori_1(n) * psi_ori_1(m) + M2_r_ori[3] * psi_ori_2(m, n)) * (
			psi_ori_2(k, l) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, l)
			+ psi_ori_2(i, l) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, l)
			+ psi_ori_2(j, l) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, l)
			) + M2_r_ori[3] * psi_ori_1(m) * (
				psi_ori_3(k, l, n) * psi_ori_2(i, j) + psi_ori_2(k, n) * psi_ori_3(i, j, l)
				+ psi_ori_2(k, l) * psi_ori_3(i, j, n) + psi_ori_1(k) * psi_ori_4(i, j, l, n)
				+ psi_ori_3(i, l, n) * psi_ori_2(j, k) + psi_ori_2(i, n) * psi_ori_3(j, k, l)
				+ psi_ori_2(i, l) * psi_ori_3(j, k, n) + psi_ori_1(i) * psi_ori_4(j, k, l, n)
				+ psi_ori_3(j, l, n) * psi_ori_2(i, k) + psi_ori_2(j, n) * psi_ori_3(i, k, l)
				+ psi_ori_2(j, l) * psi_ori_3(i, k, n) + psi_ori_1(j) * psi_ori_4(i, k, l, n)
				);

		result += M2_r_ori[3] * psi_ori_1(n) * (
			psi_ori_3(k, l, m) * psi_ori_2(i, j) + psi_ori_2(k, l) * psi_ori_3(i, j, m)
			+ psi_ori_2(k, m) * psi_ori_3(i, j, l) + psi_ori_1(k) * psi_ori_4(i, j, l, m)
			+ psi_ori_3(i, l, m) * psi_ori_2(j, k) + psi_ori_2(i, l) * psi_ori_3(j, k, m)
			+ psi_ori_2(i, m) * psi_ori_3(j, k, l) + psi_ori_1(i) * psi_ori_4(j, k, l, m)
			+ psi_ori_3(j, l, m) * psi_ori_2(i, k) + psi_ori_2(j, l) * psi_ori_3(i, k, m)
			+ psi_ori_2(j, m) * psi_ori_3(i, k, l) + psi_ori_1(j) * psi_ori_4(i, k, l, m)
			) + M2_r_ori[2] * (
				psi_ori_4(k, l, m, n) * psi_ori_2(i, j) + psi_ori_3(k, l, n) * psi_ori_3(i, j, m)
				+ psi_ori_3(k, l, m) * psi_ori_3(i, j, n) + psi_ori_2(k, l) * psi_ori_4(i, j, m, n)
				+ psi_ori_3(k, m, n) * psi_ori_3(i, j, l) + psi_ori_2(k, n) * psi_ori_4(i, j, l, m)
				+ psi_ori_2(k, m) * psi_ori_4(i, j, l, n) + psi_ori_1(k) * psi_ori_5(i, j, l, m, n)
				+ psi_ori_4(i, l, m, n) * psi_ori_2(j, k) + psi_ori_3(i, l, n) * psi_ori_3(j, k, m)
				+ psi_ori_3(i, l, m) * psi_ori_3(j, k, n) + psi_ori_2(i, l) * psi_ori_4(j, k, m, n)
				+ psi_ori_3(i, m, n) * psi_ori_3(j, k, l) + psi_ori_2(i, n) * psi_ori_4(j, k, l, m)
				+ psi_ori_2(i, m) * psi_ori_4(j, k, l, n) + psi_ori_1(i) * psi_ori_5(j, k, l, m, n)
				+ psi_ori_4(j, l, m, n) * psi_ori_2(i, k) + psi_ori_3(j, l, n) * psi_ori_3(i, k, m)
				+ psi_ori_3(j, l, m) * psi_ori_3(i, k, n) + psi_ori_2(j, l) * psi_ori_4(i, k, m, n)
				+ psi_ori_3(j, m, n) * psi_ori_3(i, k, l) + psi_ori_2(j, n) * psi_ori_4(i, k, l, m)
				+ psi_ori_2(j, m) * psi_ori_4(i, k, l, n) + psi_ori_1(j) * psi_ori_5(i, k, l, m, n)
				);

		result += (M2_r_ori[5] * psi_ori_1(n) * psi_ori_1(m) + M2_r_ori[4] * psi_ori_2(m, n)) * (
			psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_1(k)
			+ psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_1(k)
			+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, l)
			) + M2_r_ori[4] * psi_ori_1(m) * (
				psi_ori_3(i, l, n) * psi_ori_1(j) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_2(j, n) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_2(k, n)
				+ psi_ori_2(i, n) * psi_ori_2(j, l) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_3(j, l, n) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_2(k, n)
				+ psi_ori_2(i, n) * psi_ori_1(j) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_2(j, n) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_3(k, l, n)
				);

		result += M2_r_ori[4] * psi_ori_1(n) * (
			psi_ori_3(i, l, m) * psi_ori_1(j) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_2(j, m) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_2(k, m)
			+ psi_ori_2(i, m) * psi_ori_2(j, l) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_3(j, l, m) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_2(k, m)
			+ psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_3(k, l, m)
			) + M2_r_ori[3] * (
				psi_ori_4(i, l, m, n) * psi_ori_1(j) * psi_ori_1(k) + psi_ori_3(i, l, n) * psi_ori_2(j, m) * psi_ori_1(k) + psi_ori_3(i, l, n) * psi_ori_1(j) * psi_ori_2(k, m)
				+ psi_ori_3(i, l, m) * psi_ori_2(j, n) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_3(j, m, n) * psi_ori_1(k) + psi_ori_2(i, l) * psi_ori_2(j, n) * psi_ori_2(k, m)
				+ psi_ori_3(i, l, m) * psi_ori_1(j) * psi_ori_2(k, n) + psi_ori_2(i, l) * psi_ori_2(j, m) * psi_ori_2(k, n) + psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_3(k, m, n)

				+ psi_ori_3(i, m, n) * psi_ori_2(j, l) * psi_ori_1(k) + psi_ori_2(i, n) * psi_ori_3(j, l, m) * psi_ori_1(k) + psi_ori_2(i, n) * psi_ori_2(j, l) * psi_ori_2(k, m)
				+ psi_ori_2(i, m) * psi_ori_3(j, l, n) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_4(j, l, m, n) * psi_ori_1(k) + psi_ori_1(i) * psi_ori_3(j, l, n) * psi_ori_2(k, m)
				+ psi_ori_2(i, m) * psi_ori_2(j, l) * psi_ori_2(k, n) + psi_ori_1(i) * psi_ori_3(j, l, m) * psi_ori_2(k, n) + psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_3(k, m, n)

				+ psi_ori_3(i, m, n) * psi_ori_1(j) * psi_ori_2(k, l) + psi_ori_2(i, n) * psi_ori_2(j, m) * psi_ori_2(k, l) + psi_ori_2(i, n) * psi_ori_1(j) * psi_ori_3(k, l, m)
				+ psi_ori_2(i, m) * psi_ori_2(j, n) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_3(j, m, n) * psi_ori_2(k, l) + psi_ori_1(i) * psi_ori_2(j, n) * psi_ori_3(k, l, m)
				+ psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_3(k, l, n) + psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_3(k, l, n) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_4(k, l, m, n)
				);

		result += M2_r_ori[5] * psi_ori_1(n) * (
			psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l)
			+ psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_1(k) * psi_ori_1(l)
			+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, m) * psi_ori_1(l)
			+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_2(l, m)
			) + M2_r_ori[4] * (
				psi_ori_3(i, m, n) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l) + psi_ori_2(i, m) * psi_ori_2(j, n) * psi_ori_1(k) * psi_ori_1(l)
				+ psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_2(k, n) * psi_ori_1(l) + psi_ori_2(i, m) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_2(l, n)

				+ psi_ori_2(i, n) * psi_ori_2(j, m) * psi_ori_1(k) * psi_ori_1(l) + psi_ori_1(i) * psi_ori_3(j, m, n) * psi_ori_1(k) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_2(k, n) * psi_ori_1(l) + psi_ori_1(i) * psi_ori_2(j, m) * psi_ori_1(k) * psi_ori_2(l, n)

				+ psi_ori_2(i, n) * psi_ori_1(j) * psi_ori_2(k, m) * psi_ori_1(l) + psi_ori_1(i) * psi_ori_2(j, n) * psi_ori_2(k, m) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_3(k, m, n) * psi_ori_1(l) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, m) * psi_ori_2(l, n)

				+ psi_ori_2(i, n) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_2(l, m) + psi_ori_1(i) * psi_ori_2(j, n) * psi_ori_1(k) * psi_ori_2(l, m)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, n) * psi_ori_2(l, m) + psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_3(l, m, n)
				);

		result += M2_r_ori[6] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l) * psi_ori_1(m) * psi_ori_1(n)
			+ M2_r_ori[5] * (
				psi_ori_2(i, n) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l) * psi_ori_1(m)
				+ psi_ori_1(i) * psi_ori_2(j, n) * psi_ori_1(k) * psi_ori_1(l) * psi_ori_1(m)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, n) * psi_ori_1(l) * psi_ori_1(m)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_2(l, n) * psi_ori_1(m)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l) * psi_ori_2(m, n)
				);

		return result;
	}
	else {
		complex<double> result = 4.0 / 105.0 * exp(1.0i * a * beta) * pi
			* pow(beta, 2.0) * (12.0 + a * beta * (
				-12.0i + a * beta * (9.0 + 1.0i * a * beta)
				)
				);

		return result * (
			(d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k]) * d[m][n]
			+ (d[j][m] * d[k][l] + d[k][m] * d[j][l] + d[l][m] * d[j][k]) * d[i][n]
			+ (d[i][m] * d[k][l] + d[i][k] * d[l][m] + d[i][l] * d[k][m]) * d[j][n]
			+ (d[i][j] * d[l][m] + d[i][m] * d[j][l] + d[i][l] * d[j][m]) * d[k][n]
			+ (d[i][j] * d[k][m] + d[i][k] * d[j][m] + d[i][m] * d[j][k]) * d[l][n]
			);
	}
}

/*
	M2 functions: imaged
*/
complex<double> DIE::M2()
{
	return M2_r[0];
}

complex<double> DIE::M2_1(int i)
{
	return M2_r[1] * psi_1(i);
}

complex<double> DIE::M2_2(int i, int j)
{
	return M2_r[1] * psi_2(i, j) + M2_r[2] * psi_1(i) * psi_1(j);
}

complex<double> DIE::M2_3(int i, int j, int k)
{
	return M2_r[1] * psi_3(i, j, k)
		+ M2_r[2] * (
			psi_1(k) * psi_2(i, j)
			+ psi_1(i) * psi_2(j, k)
			+ psi_1(j) * psi_2(i, k)
			)
		+ M2_r[3] * psi_1(i) * psi_1(j) * psi_1(k);
}

complex<double> DIE::M2_4(int i, int j, int k, int l)
{
	return M2_r[1] * psi_4(i, j, k, l)
		+ M2_r[2] * (
			psi_2(k, l) * psi_2(i, j) + psi_1(k) * psi_3(i, j, l)
			+ psi_2(i, l) * psi_2(j, k) + psi_1(i) * psi_3(j, k, l)
			+ psi_2(j, l) * psi_2(i, k) + psi_1(j) * psi_3(i, k, l)
			+ psi_3(i, j, k) * psi_1(l)
			)
		+ M2_r[3] * (
			psi_2(i, l) * psi_1(j) * psi_1(k)
			+ psi_1(i) * psi_2(j, l) * psi_1(k)
			+ psi_1(i) * psi_1(j) * psi_2(k, l)
			+ psi_1(k) * psi_2(i, j) * psi_1(l)
			+ psi_1(i) * psi_2(j, k) * psi_1(l)
			+ psi_1(j) * psi_2(i, k) * psi_1(l)
			)
		+ M2_r[4] * psi_1(i) * psi_1(j) * psi_1(k) * psi_1(l);
}

complex<double> DIE::M2_5(int i, int j, int k, int l, int m)
{
	return M2_r[1] * psi_5(i, j, k, l, m)
		+ M2_r[2] * psi_4(i, j, k, l) * psi_1(m)

		+ M2_r[2] * (psi_4(i, j, k, m) * psi_1(l) + psi_3(i, j, k) * psi_2(l, m))
		+ M2_r[3] * psi_1(m) * (psi_3(i, j, k) * psi_1(l))


		+ (M2_r[4] * psi_1(l) * psi_1(m) + M2_r[3] * psi_2(l, m)) * (
			psi_1(k) * psi_2(i, j)
			+ psi_1(i) * psi_2(j, k)
			+ psi_1(j) * psi_2(i, k)
			)

		+ M2_r[3] * psi_1(l) * (
			psi_2(k, m) * psi_2(i, j) + psi_1(k) * psi_3(i, j, m)
			+ psi_2(i, m) * psi_2(j, k) + psi_1(i) * psi_3(j, k, m)
			+ psi_2(j, m) * psi_2(i, k) + psi_1(j) * psi_3(i, k, m)
			)

		+ M2_r[3] * psi_1(m) * (
			psi_2(k, l) * psi_2(i, j) + psi_1(k) * psi_3(i, j, l)
			+ psi_2(i, l) * psi_2(j, k) + psi_1(i) * psi_3(j, k, l)
			+ psi_2(j, l) * psi_2(i, k) + psi_1(j) * psi_3(i, k, l)
			)

		+ M2_r[2] * (
			psi_3(k, l, m) * psi_2(i, j) + psi_2(k, l) * psi_3(i, j, m)
			+ psi_2(k, m) * psi_3(i, j, l) + psi_1(k) * psi_4(i, j, l, m)
			+ psi_3(i, l, m) * psi_2(j, k) + psi_2(i, l) * psi_3(j, k, m)
			+ psi_2(i, m) * psi_3(j, k, l) + psi_1(i) * psi_4(j, k, l, m)
			+ psi_3(j, l, m) * psi_2(i, k) + psi_2(j, l) * psi_3(i, k, m)
			+ psi_2(j, m) * psi_3(i, k, l) + psi_1(j) * psi_4(i, k, l, m)
			)

		+ M2_r[4] * psi_1(m) * (
			psi_2(i, l) * psi_1(j) * psi_1(k)
			+ psi_1(i) * psi_2(j, l) * psi_1(k)
			+ psi_1(i) * psi_1(j) * psi_2(k, l)
			)

		+ M2_r[3] * (
			psi_3(i, l, m) * psi_1(j) * psi_1(k) + psi_2(i, l) * psi_2(j, m) * psi_1(k) + psi_2(i, l) * psi_1(j) * psi_2(k, m)
			+ psi_2(i, m) * psi_2(j, l) * psi_1(k) + psi_1(i) * psi_3(j, l, m) * psi_1(k) + psi_1(i) * psi_2(j, l) * psi_2(k, m)
			+ psi_2(i, m) * psi_1(j) * psi_2(k, l) + psi_1(i) * psi_2(j, m) * psi_2(k, l) + psi_1(i) * psi_1(j) * psi_3(k, l, m)
			)

		+ M2_r[4] * (
			psi_2(i, m) * psi_1(j) * psi_1(k) * psi_1(l)
			+ psi_1(i) * psi_2(j, m) * psi_1(k) * psi_1(l)
			+ psi_1(i) * psi_1(j) * psi_2(k, m) * psi_1(l)
			+ psi_1(i) * psi_1(j) * psi_1(k) * psi_2(l, m)
			)

		+ M2_r[5] * psi_1(i) * psi_1(j) * psi_1(k) * psi_1(l) * psi_1(m);
}

complex<double> DIE::M2_6(int i, int j, int k, int l, int m, int n)
{
		complex<double> result = M2_r[1] * psi_6(i, j, k, l, m, n) + M2_r[2] * psi_1(n) * psi_5(i, j, k, l, m)
			+ M2_r[2] * psi_5(i, j, k, l, n) * psi_1(m) + M2_r[2] * psi_4(i, j, k, l) * psi_2(m, n)
			+ M2_r[3] * psi_4(i, j, k, l) * psi_1(m) * psi_1(n);

		result += M2_r[3] * psi_1(n) * (psi_4(i, j, k, m) * psi_1(l) + psi_3(i, j, k) * psi_2(l, m))
			+ M2_r[2] * (
				psi_5(i, j, k, m, n) * psi_1(l) + psi_4(i, j, k, m) * psi_2(l, n)
				+ psi_4(i, j, k, n) * psi_2(l, m) + psi_3(i, j, k) * psi_3(l, m, n)
				)
			+ M2_r[4] * psi_1(n) * psi_1(m) * psi_3(i, j, k) * psi_1(l)
			+ M2_r[3] * (
				psi_2(m, n) * psi_3(i, j, k) * psi_1(l)
				+ psi_1(m) * psi_4(i, j, k, n) * psi_1(l)
				+ psi_1(m) * psi_3(i, j, k) * psi_2(l, n)
				);


		result += (M2_r[5] * psi_1(l) * psi_1(m) * psi_1(n) + M2_r[4] * psi_2(l, n) * psi_1(m)
			+ M2_r[4] * psi_1(l) * psi_2(m, n) + M2_r[4] * psi_1(n) * psi_2(l, m) + M2_r[3] * psi_3(l, m, n))
			* (
				psi_1(k) * psi_2(i, j)
				+ psi_1(i) * psi_2(j, k)
				+ psi_1(j) * psi_2(i, k)
				)
			+ (M2_r[4] * psi_1(l) * psi_1(m) + M2_r[3] * psi_2(l, m)) * (
				psi_2(k, n) * psi_2(i, j) + psi_1(k) * psi_3(i, j, n)
				+ psi_2(i, n) * psi_2(j, k) + psi_1(i) * psi_3(j, k, n)
				+ psi_2(j, n) * psi_2(i, k) + psi_1(j) * psi_3(i, k, n)
				);

		result += (M2_r[4] * psi_1(l) * psi_1(n) + M2_r[3] * psi_2(l, n)) * (
			psi_2(k, m) * psi_2(i, j) + psi_1(k) * psi_3(i, j, m)
			+ psi_2(i, m) * psi_2(j, k) + psi_1(i) * psi_3(j, k, m)
			+ psi_2(j, m) * psi_2(i, k) + psi_1(j) * psi_3(i, k, m)
			) + M2_r[3] * psi_1(l) * (
				psi_3(k, m, n) * psi_2(i, j) + psi_2(k, n) * psi_3(i, j, m)
				+ psi_2(k, m) * psi_3(i, j, n) + psi_1(k) * psi_4(i, j, m, n)
				+ psi_3(i, m, n) * psi_2(j, k) + psi_2(i, n) * psi_3(j, k, m)
				+ psi_2(i, m) * psi_3(j, k, n) + psi_1(i) * psi_4(j, k, m, n)
				+ psi_3(j, m, n) * psi_2(i, k) + psi_2(j, n) * psi_3(i, k, m)
				+ psi_2(j, m) * psi_3(i, k, n) + psi_1(j) * psi_4(i, k, m, n)
				);

		result += (M2_r[4] * psi_1(n) * psi_1(m) + M2_r[3] * psi_2(m, n)) * (
			psi_2(k, l) * psi_2(i, j) + psi_1(k) * psi_3(i, j, l)
			+ psi_2(i, l) * psi_2(j, k) + psi_1(i) * psi_3(j, k, l)
			+ psi_2(j, l) * psi_2(i, k) + psi_1(j) * psi_3(i, k, l)
			) + M2_r[3] * psi_1(m) * (
				psi_3(k, l, n) * psi_2(i, j) + psi_2(k, n) * psi_3(i, j, l)
				+ psi_2(k, l) * psi_3(i, j, n) + psi_1(k) * psi_4(i, j, l, n)
				+ psi_3(i, l, n) * psi_2(j, k) + psi_2(i, n) * psi_3(j, k, l)
				+ psi_2(i, l) * psi_3(j, k, n) + psi_1(i) * psi_4(j, k, l, n)
				+ psi_3(j, l, n) * psi_2(i, k) + psi_2(j, n) * psi_3(i, k, l)
				+ psi_2(j, l) * psi_3(i, k, n) + psi_1(j) * psi_4(i, k, l, n)
				);

		result += M2_r[3] * psi_1(n) * (
			psi_3(k, l, m) * psi_2(i, j) + psi_2(k, l) * psi_3(i, j, m)
			+ psi_2(k, m) * psi_3(i, j, l) + psi_1(k) * psi_4(i, j, l, m)
			+ psi_3(i, l, m) * psi_2(j, k) + psi_2(i, l) * psi_3(j, k, m)
			+ psi_2(i, m) * psi_3(j, k, l) + psi_1(i) * psi_4(j, k, l, m)
			+ psi_3(j, l, m) * psi_2(i, k) + psi_2(j, l) * psi_3(i, k, m)
			+ psi_2(j, m) * psi_3(i, k, l) + psi_1(j) * psi_4(i, k, l, m)
			) + M2_r[2] * (
				psi_4(k, l, m, n) * psi_2(i, j) + psi_3(k, l, n) * psi_3(i, j, m)
				+ psi_3(k, l, m) * psi_3(i, j, n) + psi_2(k, l) * psi_4(i, j, m, n)
				+ psi_3(k, m, n) * psi_3(i, j, l) + psi_2(k, n) * psi_4(i, j, l, m)
				+ psi_2(k, m) * psi_4(i, j, l, n) + psi_1(k) * psi_5(i, j, l, m, n)
				+ psi_4(i, l, m, n) * psi_2(j, k) + psi_3(i, l, n) * psi_3(j, k, m)
				+ psi_3(i, l, m) * psi_3(j, k, n) + psi_2(i, l) * psi_4(j, k, m, n)
				+ psi_3(i, m, n) * psi_3(j, k, l) + psi_2(i, n) * psi_4(j, k, l, m)
				+ psi_2(i, m) * psi_4(j, k, l, n) + psi_1(i) * psi_5(j, k, l, m, n)
				+ psi_4(j, l, m, n) * psi_2(i, k) + psi_3(j, l, n) * psi_3(i, k, m)
				+ psi_3(j, l, m) * psi_3(i, k, n) + psi_2(j, l) * psi_4(i, k, m, n)
				+ psi_3(j, m, n) * psi_3(i, k, l) + psi_2(j, n) * psi_4(i, k, l, m)
				+ psi_2(j, m) * psi_4(i, k, l, n) + psi_1(j) * psi_5(i, k, l, m, n)
				);

		result += (M2_r[5] * psi_1(n) * psi_1(m) + M2_r[4] * psi_2(m, n)) * (
			psi_2(i, l) * psi_1(j) * psi_1(k)
			+ psi_1(i) * psi_2(j, l) * psi_1(k)
			+ psi_1(i) * psi_1(j) * psi_2(k, l)
			) + M2_r[4] * psi_1(m) * (
				psi_3(i, l, n) * psi_1(j) * psi_1(k) + psi_2(i, l) * psi_2(j, n) * psi_1(k) + psi_2(i, l) * psi_1(j) * psi_2(k, n)
				+ psi_2(i, n) * psi_2(j, l) * psi_1(k) + psi_1(i) * psi_3(j, l, n) * psi_1(k) + psi_1(i) * psi_2(j, l) * psi_2(k, n)
				+ psi_2(i, n) * psi_1(j) * psi_2(k, l) + psi_1(i) * psi_2(j, n) * psi_2(k, l) + psi_1(i) * psi_1(j) * psi_3(k, l, n)
				);

		result += M2_r[4] * psi_1(n) * (
			psi_3(i, l, m) * psi_1(j) * psi_1(k) + psi_2(i, l) * psi_2(j, m) * psi_1(k) + psi_2(i, l) * psi_1(j) * psi_2(k, m)
			+ psi_2(i, m) * psi_2(j, l) * psi_1(k) + psi_1(i) * psi_3(j, l, m) * psi_1(k) + psi_1(i) * psi_2(j, l) * psi_2(k, m)
			+ psi_2(i, m) * psi_1(j) * psi_2(k, l) + psi_1(i) * psi_2(j, m) * psi_2(k, l) + psi_1(i) * psi_1(j) * psi_3(k, l, m)
			) + M2_r[3] * (
				psi_4(i, l, m, n) * psi_1(j) * psi_1(k) + psi_3(i, l, n) * psi_2(j, m) * psi_1(k) + psi_3(i, l, n) * psi_1(j) * psi_2(k, m)
				+ psi_3(i, l, m) * psi_2(j, n) * psi_1(k) + psi_2(i, l) * psi_3(j, m, n) * psi_1(k) + psi_2(i, l) * psi_2(j, n) * psi_2(k, m)
				+ psi_3(i, l, m) * psi_1(j) * psi_2(k, n) + psi_2(i, l) * psi_2(j, m) * psi_2(k, n) + psi_2(i, l) * psi_1(j) * psi_3(k, m, n)

				+ psi_3(i, m, n) * psi_2(j, l) * psi_1(k) + psi_2(i, n) * psi_3(j, l, m) * psi_1(k) + psi_2(i, n) * psi_2(j, l) * psi_2(k, m)
				+ psi_2(i, m) * psi_3(j, l, n) * psi_1(k) + psi_1(i) * psi_4(j, l, m, n) * psi_1(k) + psi_1(i) * psi_3(j, l, n) * psi_2(k, m)
				+ psi_2(i, m) * psi_2(j, l) * psi_2(k, n) + psi_1(i) * psi_3(j, l, m) * psi_2(k, n) + psi_1(i) * psi_2(j, l) * psi_3(k, m, n)

				+ psi_3(i, m, n) * psi_1(j) * psi_2(k, l) + psi_2(i, n) * psi_2(j, m) * psi_2(k, l) + psi_2(i, n) * psi_1(j) * psi_3(k, l, m)
				+ psi_2(i, m) * psi_2(j, n) * psi_2(k, l) + psi_1(i) * psi_3(j, m, n) * psi_2(k, l) + psi_1(i) * psi_2(j, n) * psi_3(k, l, m)
				+ psi_2(i, m) * psi_1(j) * psi_3(k, l, n) + psi_1(i) * psi_2(j, m) * psi_3(k, l, n) + psi_1(i) * psi_1(j) * psi_4(k, l, m, n)
				);

		result += M2_r[5] * psi_1(n) * (
			psi_2(i, m) * psi_1(j) * psi_1(k) * psi_1(l)
			+ psi_1(i) * psi_2(j, m) * psi_1(k) * psi_1(l)
			+ psi_1(i) * psi_1(j) * psi_2(k, m) * psi_1(l)
			+ psi_1(i) * psi_1(j) * psi_1(k) * psi_2(l, m)
			) + M2_r[4] * (
				psi_3(i, m, n) * psi_1(j) * psi_1(k) * psi_1(l) + psi_2(i, m) * psi_2(j, n) * psi_1(k) * psi_1(l)
				+ psi_2(i, m) * psi_1(j) * psi_2(k, n) * psi_1(l) + psi_2(i, m) * psi_1(j) * psi_1(k) * psi_2(l, n)

				+ psi_2(i, n) * psi_2(j, m) * psi_1(k) * psi_1(l) + psi_1(i) * psi_3(j, m, n) * psi_1(k) * psi_1(l)
				+ psi_1(i) * psi_2(j, m) * psi_2(k, n) * psi_1(l) + psi_1(i) * psi_2(j, m) * psi_1(k) * psi_2(l, n)

				+ psi_2(i, n) * psi_1(j) * psi_2(k, m) * psi_1(l) + psi_1(i) * psi_2(j, n) * psi_2(k, m) * psi_1(l)
				+ psi_1(i) * psi_1(j) * psi_3(k, m, n) * psi_1(l) + psi_1(i) * psi_1(j) * psi_2(k, m) * psi_2(l, n)

				+ psi_2(i, n) * psi_1(j) * psi_1(k) * psi_2(l, m) + psi_1(i) * psi_2(j, n) * psi_1(k) * psi_2(l, m)
				+ psi_1(i) * psi_1(j) * psi_2(k, n) * psi_2(l, m) + psi_1(i) * psi_1(j) * psi_1(k) * psi_3(l, m, n)
				);

		result += M2_r[6] * psi_1(i) * psi_1(j) * psi_1(k) * psi_1(l) * psi_1(m) * psi_1(n)
			+ M2_r[5] * (
				psi_2(i, n) * psi_1(j) * psi_1(k) * psi_1(l) * psi_1(m)
				+ psi_1(i) * psi_2(j, n) * psi_1(k) * psi_1(l) * psi_1(m)
				+ psi_1(i) * psi_1(j) * psi_2(k, n) * psi_1(l) * psi_1(m)
				+ psi_1(i) * psi_1(j) * psi_1(k) * psi_2(l, n) * psi_1(m)
				+ psi_1(i) * psi_1(j) * psi_1(k) * psi_1(l) * psi_2(m, n)
				);

		return result;
}