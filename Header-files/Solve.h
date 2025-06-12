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
};