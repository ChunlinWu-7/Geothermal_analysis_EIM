
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

/*
	Domain integrals of Helmholtz function H_0
				Imaged!
*/
complex<double> DIE::Helm_0()
{
	return M0();
}

complex<double> DIE::Helm_1(int i)
{
	return M0_1(i);
}

complex<double> DIE::Helm_2(int i, int j)
{
	return M0_2(i, j);
}

complex<double> DIE::Helm_3(int i, int j, int k)
{
	return M0_3(i, j, k);
}

complex<double> DIE::Helm_4(int i, int j, int k, int l)
{
	return M0_4(i, j, k, l);
}

/*
	Domain integrals of Helmholtz function H_ori_0
				Original!
*/
complex<double> DIE::Helm_ori_0()
{
	return M0_ori();
}

complex<double> DIE::Helm_ori_1(int i)
{
	return M0_ori_1(i);
}

complex<double> DIE::Helm_ori_2(int i, int j)
{
	return M0_ori_2(i, j);
}

complex<double> DIE::Helm_ori_3(int i, int j, int k)
{
	return M0_ori_3(i, j, k);
}

complex<double> DIE::Helm_ori_4(int i, int j, int k, int l)
{
	return M0_ori_4(i, j, k, l);
}


/*
	Domain integrals of Helmholtz function H_1
				Imaged!
*/
complex<double> DIE::Helm_1_0(int p)
{
	return 1.0i / beta * M1_1(p) + difx[p] * Helm_0();
}

complex<double> DIE::Helm_1_1(int p, int i)
{
	return 1.0i / beta * M1_2(p, i) + d[p][i] * Helm_0()
		+ difx[p] * Helm_1(i);
}

complex<double> DIE::Helm_1_2(int p, int i, int j)
{
	return 1.0i / beta * M1_3(p, i, j) + d[p][i] * Helm_1(j)
		+ d[p][j] * Helm_1(i) + difx[p] * Helm_2(i, j);
}

complex<double> DIE::Helm_1_3(int p, int i, int j, int k)
{
	return 1.0i / beta * M1_4(p, i, j, k) + d[p][i] * Helm_2(j, k)
		+ d[p][j] * Helm_2(i, k) + d[p][k] * Helm_2(i, j)
		+ difx[p] * Helm_3(i, j, k);
}

complex<double> DIE::Helm_1_4(int p, int i, int j, int k, int l)
{
	return 1.0i / beta * M1_5(p, i, j, k, l) + d[p][i] * Helm_3(j, k, l)
		+ d[p][j] * Helm_3(i, k, l) + d[p][k] * Helm_3(i, j, l)
		+ d[p][l] * Helm_3(i, j, k) + difx[p] * Helm_4(i, j, k, l);
}


/*
	Domain integrals of Helmholtz function H_1
				Original!
*/
complex<double> DIE::Helm_ori_1_0(int p)
{
	return 1.0i / beta * M1_ori_1(p) + difx_ori[p] * Helm_ori_0();
}

complex<double> DIE::Helm_ori_1_1(int p, int i)
{
	return 1.0i / beta * M1_ori_2(p, i) + d[p][i] * Helm_ori_0()
		+ difx_ori[p] * Helm_ori_1(i);
}

complex<double> DIE::Helm_ori_1_2(int p, int i, int j)
{
	return 1.0i / beta * M1_ori_3(p, i, j) + d[p][i] * Helm_ori_1(j)
		+ d[p][j] * Helm_ori_1(i) + difx_ori[p] * Helm_ori_2(i, j);
}

complex<double> DIE::Helm_ori_1_3(int p, int i, int j, int k)
{
	return 1.0i / beta * M1_ori_4(p, i, j, k) + d[p][i] * Helm_ori_2(j, k)
		+ d[p][j] * Helm_ori_2(i, k) + d[p][k] * Helm_ori_2(i, j)
		+ difx_ori[p] * Helm_ori_3(i, j, k);
}

complex<double> DIE::Helm_ori_1_4(int p, int i, int j, int k, int l)
{
	return 1.0i / beta * M1_ori_5(p, i, j, k, l) + d[p][i] * Helm_ori_3(j, k, l)
		+ d[p][j] * Helm_ori_3(i, k, l) + d[p][k] * Helm_ori_3(i, j, l)
		+ d[p][l] * Helm_ori_3(i, j, k) + difx_ori[p] * Helm_ori_4(i, j, k, l);
}


/*
	Domain integrals of Helmholtz function H_2
				imaged!
*/
complex<double> DIE::Helm_2_0(int p, int q)
{
	return -1.0 / pow(beta, 2.0) * M2_2(p, q) - 1.0i / pow(beta, 3.0) * M1_2(p, q)
		+ 1.0i / beta * d[p][q] * M1() + difx[p] * Helm_1_0(q) + difx[q] * Helm_1_0(p)
		- difx[p] * difx[q] * Helm_0();
}

complex<double> DIE::Helm_2_1(int p, int q, int i)
{
	return -1.0 / pow(beta, 2.0) * M2_3(p, q, i) - 1.0i / pow(beta, 3.0) * M1_3(p, q, i)
		+ 1.0i / beta * d[p][q] * M1_1(i)
		+ d[p][i] * Helm_1_0(q) + d[q][i] * Helm_1_0(p)
		+ difx[p] * Helm_1_1(q, i) + difx[q] * Helm_1_1(p, i)
		- (d[p][i] * difx[q] + d[q][i] * difx[p]) * Helm_0()
		- difx[p] * difx[q] * Helm_1(i);
}

complex<double> DIE::Helm_2_2(int p, int q, int i, int j)
{
	return -1.0 / pow(beta, 2.0) * M2_4(p, q, i, j) - 1.0i / pow(beta, 3.0) * M1_4(p, q, i, j)
		+ 1.0i / beta * d[p][q] * M1_2(i, j)
		+ d[p][i] * Helm_1_1(q, j) + d[q][i] * Helm_1_1(p, j)
		+ d[p][j] * Helm_1_1(q, i) + d[q][j] * Helm_1_1(p, i)
		+ difx[p] * Helm_1_2(q, i, j) + difx[q] * Helm_1_2(p, i, j)
		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_0()
		- (d[p][i] * difx[q] + d[q][i] * difx[p]) * Helm_1(j)
		- (d[p][j] * difx[q] + d[q][j] * difx[p]) * Helm_1(i)
		- difx[p] * difx[q] * Helm_2(i, j);
}

complex<double> DIE::Helm_2_3(int p, int q, int i, int j, int k)
{
	return -1.0 / pow(beta, 2.0) * M2_5(p, q, i, j, k) - 1.0i / pow(beta, 3.0) * M1_5(p, q, i, j, k)
		+ 1.0i / beta * d[p][q] * M1_3(i, j, k)
		+ d[p][i] * Helm_1_2(q, j, k) + d[q][i] * Helm_1_2(p, j, k)
		+ d[p][j] * Helm_1_2(q, i, k) + d[q][j] * Helm_1_2(p, i, k)
		+ d[p][k] * Helm_1_2(q, i, j) + d[q][k] * Helm_1_2(p, i, j)
		+ difx[p] * Helm_1_3(q, i, j, k) + difx[q] * Helm_1_3(p, i, j, k)

		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_1(k)
		- (d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Helm_1(j)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Helm_1(i)
		- (d[p][i] * difx[q] + d[q][i] * difx[p]) * Helm_2(j, k)
		- (d[p][j] * difx[q] + d[q][j] * difx[p]) * Helm_2(i, k)
		- (d[p][k] * difx[q] + d[q][k] * difx[p]) * Helm_2(i, j)
		- difx[p] * difx[q] * Helm_3(i, j, k);

}

complex<double> DIE::Helm_2_4(int p, int q, int i, int j, int k, int l)
{
	return -1.0 / pow(beta, 2.0) * M2_6(p, q, i, j, k, l) - 1.0i / pow(beta, 3.0) * M1_6(p, q, i, j, k, l)
		+ 1.0i / beta * d[p][q] * M1_4(i, j, k, l)
		+ d[p][i] * Helm_1_3(q, j, k, l) + d[q][i] * Helm_1_3(p, j, k, l)
		+ d[p][j] * Helm_1_3(q, i, k, l) + d[q][j] * Helm_1_3(p, i, k, l)
		+ d[p][k] * Helm_1_3(q, i, j, l) + d[q][k] * Helm_1_3(p, i, j, l)
		+ d[p][l] * Helm_1_3(q, i, j, k) + d[q][l] * Helm_1_3(p, i, j, k)
		+ difx[p] * Helm_1_4(q, i, j, k, l) + difx[q] * Helm_1_4(p, i, j, k, l)

		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_2(k, l)
		- (d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Helm_2(j, l)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Helm_2(i, l)
		- (d[p][i] * d[q][l] + d[q][i] * d[p][l]) * Helm_2(j, k)
		- (d[p][j] * d[q][l] + d[q][j] * d[p][l]) * Helm_2(i, k)
		- (d[p][k] * d[q][l] + d[q][k] * d[p][l]) * Helm_2(i, j)

		- (d[p][i] * difx[q] + d[q][i] * difx[p]) * Helm_3(j, k, l)
		- (d[p][j] * difx[q] + d[q][j] * difx[p]) * Helm_3(i, k, l)
		- (d[p][k] * difx[q] + d[q][k] * difx[p]) * Helm_3(i, j, l)
		- (d[p][l] * difx[q] + d[q][l] * difx[p]) * Helm_3(i, j, k)

		- difx[p] * difx[q] * Helm_4(i, j, k, l);
}

/*
	Domain integrals of Helmholtz function H_2
				original!
*/
complex<double> DIE::Helm_ori_2_0(int p, int q)
{
	return -1.0 / pow(beta, 2.0) * M2_ori_2(p, q) - 1.0i / pow(beta, 3.0) * M1_ori_2(p, q)
		+ 1.0i / beta * d[p][q] * M1_ori() + difx_ori[p] * Helm_ori_1_0(q) + difx_ori[q] * Helm_ori_1_0(p)
		- difx_ori[p] * difx_ori[q] * Helm_ori_0();
}

complex<double> DIE::Helm_ori_2_1(int p, int q, int i)
{
	return -1.0 / pow(beta, 2.0) * M2_ori_3(p, q, i) - 1.0i / pow(beta, 3.0) * M1_ori_3(p, q, i)
		+ 1.0i / beta * d[p][q] * M1_ori_1(i)
		+ d[p][i] * Helm_ori_1_0(q) + d[q][i] * Helm_ori_1_0(p)
		+ difx_ori[p] * Helm_ori_1_1(q, i) + difx_ori[q] * Helm_ori_1_1(p, i)
		- (d[p][i] * difx_ori[q] + d[q][i] * difx_ori[p]) * Helm_ori_0()
		- difx_ori[p] * difx_ori[q] * Helm_ori_1(i);
}

complex<double> DIE::Helm_ori_2_2(int p, int q, int i, int j)
{
	return -1.0 / pow(beta, 2.0) * M2_ori_4(p, q, i, j) - 1.0i / pow(beta, 3.0) * M1_ori_4(p, q, i, j)
		+ 1.0i / beta * d[p][q] * M1_ori_2(i, j)
		+ d[p][i] * Helm_ori_1_1(q, j) + d[q][i] * Helm_ori_1_1(p, j)
		+ d[p][j] * Helm_ori_1_1(q, i) + d[q][j] * Helm_ori_1_1(p, i)
		+ difx_ori[p] * Helm_ori_1_2(q, i, j) + difx_ori[q] * Helm_ori_1_2(p, i, j)
		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_ori_0()
		- (d[p][i] * difx_ori[q] + d[q][i] * difx_ori[p]) * Helm_ori_1(j)
		- (d[p][j] * difx_ori[q] + d[q][j] * difx_ori[p]) * Helm_ori_1(i)
		- difx_ori[p] * difx_ori[q] * Helm_ori_2(i, j);
}

complex<double> DIE::Helm_ori_2_3(int p, int q, int i, int j, int k)
{
	return -1.0 / pow(beta, 2.0) * M2_ori_5(p, q, i, j, k) - 1.0i / pow(beta, 3.0) * M1_ori_5(p, q, i, j, k)
		+ 1.0i / beta * d[p][q] * M1_ori_3(i, j, k)
		+ d[p][i] * Helm_ori_1_2(q, j, k) + d[q][i] * Helm_ori_1_2(p, j, k)
		+ d[p][j] * Helm_ori_1_2(q, i, k) + d[q][j] * Helm_ori_1_2(p, i, k)
		+ d[p][k] * Helm_ori_1_2(q, i, j) + d[q][k] * Helm_ori_1_2(p, i, j)
		+ difx_ori[p] * Helm_ori_1_3(q, i, j, k) + difx_ori[q] * Helm_ori_1_3(p, i, j, k)

		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_ori_1(k)
		- (d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Helm_ori_1(j)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Helm_ori_1(i)
		- (d[p][i] * difx_ori[q] + d[q][i] * difx_ori[p]) * Helm_ori_2(j, k)
		- (d[p][j] * difx_ori[q] + d[q][j] * difx_ori[p]) * Helm_ori_2(i, k)
		- (d[p][k] * difx_ori[q] + d[q][k] * difx_ori[p]) * Helm_ori_2(i, j)
		- difx_ori[p] * difx_ori[q] * Helm_ori_3(i, j, k);

}

complex<double> DIE::Helm_ori_2_4(int p, int q, int i, int j, int k, int l)
{
	return -1.0 / pow(beta, 2.0) * M2_ori_6(p, q, i, j, k, l) - 1.0i / pow(beta, 3.0) * M1_ori_6(p, q, i, j, k, l)
		+ 1.0i / beta * d[p][q] * M1_ori_4(i, j, k, l)
		+ d[p][i] * Helm_ori_1_3(q, j, k, l) + d[q][i] * Helm_ori_1_3(p, j, k, l)
		+ d[p][j] * Helm_ori_1_3(q, i, k, l) + d[q][j] * Helm_ori_1_3(p, i, k, l)
		+ d[p][k] * Helm_ori_1_3(q, i, j, l) + d[q][k] * Helm_ori_1_3(p, i, j, l)
		+ d[p][l] * Helm_ori_1_3(q, i, j, k) + d[q][l] * Helm_ori_1_3(p, i, j, k)
		+ difx_ori[p] * Helm_ori_1_4(q, i, j, k, l) + difx_ori[q] * Helm_ori_1_4(p, i, j, k, l)

		- (d[p][i] * d[q][j] + d[q][i] * d[p][j]) * Helm_ori_2(k, l)
		- (d[p][i] * d[q][k] + d[q][i] * d[p][k]) * Helm_ori_2(j, l)
		- (d[p][j] * d[q][k] + d[q][j] * d[p][k]) * Helm_ori_2(i, l)
		- (d[p][i] * d[q][l] + d[q][i] * d[p][l]) * Helm_ori_2(j, k)
		- (d[p][j] * d[q][l] + d[q][j] * d[p][l]) * Helm_ori_2(i, k)
		- (d[p][k] * d[q][l] + d[q][k] * d[p][l]) * Helm_ori_2(i, j)

		- (d[p][i] * difx_ori[q] + d[q][i] * difx_ori[p]) * Helm_ori_3(j, k, l)
		- (d[p][j] * difx_ori[q] + d[q][j] * difx_ori[p]) * Helm_ori_3(i, k, l)
		- (d[p][k] * difx_ori[q] + d[q][k] * difx_ori[p]) * Helm_ori_3(i, j, l)
		- (d[p][l] * difx_ori[q] + d[q][l] * difx_ori[p]) * Helm_ori_3(i, j, k)

		- difx_ori[p] * difx_ori[q] * Helm_ori_4(i, j, k, l);
}