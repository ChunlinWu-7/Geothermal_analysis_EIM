
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

/*
	M0 functions: imaged

	No need to discuss whether it can reach dist_flag due to
	no singularity!
*/

complex<double> DIE::M0()
{
	return M0_r[0];
}

complex<double> DIE::M0_1(int i)
{
	return M0_r[1] * psi_1(i);
}

complex<double> DIE::M0_2(int i, int j)
{
	return M0_r[1] * psi_2(i, j) + M0_r[2] * psi_1(i) * psi_1(j);
}

complex<double> DIE::M0_3(int i, int j, int k)
{
	return M0_r[1] * psi_3(i, j, k)
		+ M0_r[2] * (
			psi_1(k) * psi_2(i, j)
			+ psi_1(i) * psi_2(j, k)
			+ psi_1(j) * psi_2(i, k)
			)
		+ M0_r[3] * psi_1(i) * psi_1(j) * psi_1(k);
}

complex<double> DIE::M0_4(int i, int j, int k, int l)
{
	return M0_r[1] * psi_4(i, j, k, l)
		+ M0_r[2] * (
			psi_2(k, l) * psi_2(i, j) + psi_1(k) * psi_3(i, j, l)
			+ psi_2(i, l) * psi_2(j, k) + psi_1(i) * psi_3(j, k, l)
			+ psi_2(j, l) * psi_2(i, k) + psi_1(j) * psi_3(i, k, l)
			+ psi_3(i, j, k) * psi_1(l)
			)
		+ M0_r[3] * (
			psi_2(i, l) * psi_1(j) * psi_1(k)
			+ psi_1(i) * psi_2(j, l) * psi_1(k)
			+ psi_1(i) * psi_1(j) * psi_2(k, l)
			+ psi_1(k) * psi_2(i, j) * psi_1(l)
			+ psi_1(i) * psi_2(j, k) * psi_1(l)
			+ psi_1(j) * psi_2(i, k) * psi_1(l)
			)
		+ M0_r[4] * psi_1(i) * psi_1(j) * psi_1(k) * psi_1(l);
}

/*
	M0 functions: original
*/
complex<double> DIE::M0_ori()
{
	if (dist_flag != 0) {
		return M0_r_ori[0];
	}
	else {
		return 4.0 * pi * (-1.0 + exp(1.0i * a * beta) * (1.0 - 1.0i * a * beta)) / (beta * beta);
	}
}

complex<double> DIE::M0_ori_1(int i)
{
	if (dist_flag != 0) {
		return M0_r_ori[1] * psi_ori_1(i);
	}
	else {
		return 0.0 + 0.0i;
	}
}

complex<double> DIE::M0_ori_2(int i, int j)
{
	if (dist_flag != 0) {
		return M0_r_ori[1] * psi_ori_2(i, j) + M0_r_ori[2] * psi_ori_1(i) * psi_ori_1(j);
	}
	else {
		return d[i][j] * 4.0i / 3.0 * pi * exp(1.0i * a * beta) * (1.0i + a * beta);
	}
}

complex<double> DIE::M0_ori_3(int i, int j, int k)
{
	if (dist_flag != 0) {
		return M0_r_ori[1] * psi_ori_3(i, j, k)
			+ M0_r_ori[2] * (
				psi_ori_1(k) * psi_ori_2(i, j)
				+ psi_ori_1(i) * psi_ori_2(j, k)
				+ psi_ori_1(j) * psi_ori_2(i, k)
				)
			+ M0_r_ori[3] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k);
	}
	else {
		return 0.0 + 0.0i;
	}
}

complex<double> DIE::M0_ori_4(int i, int j, int k, int l)
{
	if (dist_flag != 0) {
		return M0_r_ori[1] * psi_ori_4(i, j, k, l)
			+ M0_r_ori[2] * (
				psi_ori_2(k, l) * psi_ori_2(i, j) + psi_ori_1(k) * psi_ori_3(i, j, l)
				+ psi_ori_2(i, l) * psi_ori_2(j, k) + psi_ori_1(i) * psi_ori_3(j, k, l)
				+ psi_ori_2(j, l) * psi_ori_2(i, k) + psi_ori_1(j) * psi_ori_3(i, k, l)
				+ psi_ori_3(i, j, k) * psi_ori_1(l)
				)
			+ M0_r_ori[3] * (
				psi_ori_2(i, l) * psi_ori_1(j) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_2(j, l) * psi_ori_1(k)
				+ psi_ori_1(i) * psi_ori_1(j) * psi_ori_2(k, l)
				+ psi_ori_1(k) * psi_ori_2(i, j) * psi_ori_1(l)
				+ psi_ori_1(i) * psi_ori_2(j, k) * psi_ori_1(l)
				+ psi_ori_1(j) * psi_ori_2(i, k) * psi_ori_1(l)
				)
			+ M0_r_ori[4] * psi_ori_1(i) * psi_ori_1(j) * psi_ori_1(k) * psi_ori_1(l);
	}
	else {
		return (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[i][l] * d[j][k]) * (
			4.0 / 15.0 * exp(1.0i * a * beta) * pi * beta * beta * (1.0 - 1.0i * a * beta)
			);
	}
}