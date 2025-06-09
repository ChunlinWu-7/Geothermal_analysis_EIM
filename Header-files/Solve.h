#pragma once

#include <stdio.h>
#include "eigen-3.4.0/Eigen/Dense"
#include <complex>

using namespace std;
using namespace Eigen;

class Solve_BC
{
public:

	void Start_solve(int nsolve, int npart, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT, Ref<VectorXcd> RHS, Ref<VectorXd> H_source, Ref<VectorXcd> U);

private:

	//MatrixXcd AA, XX; VectorXcd BB;

	/*
		Math: AA XX = BB, solve XX in linear algebra
	*/
	//void Solve_matrix();
};