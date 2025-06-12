#pragma once
/*
	This EIM_run.h file calls all function:

	(1) call equivalent conditions for EIM;
	(2) Solve the coefficient matrix to obtain eigen-fields;
	(3) Postprocess
*/

#include <iostream>
#include <complex>
#include "eigen-3.4.0/Eigen/Dense"
#include "Equivalent_conditions.h"
#include "Solve.h"
#include "Input_parameters.h"

using namespace Eigen;
using namespace std;

// Construct EIM matrix
Reading_inputs* construct_EIM(Reading_inputs& READ)
{
	EIM_integrals EIM; 

	EIM.addInclusion(READ.nsolve, READ.num, READ.eigen_point, READ.radius, READ.eigen_mat, \
		READ.HMAT, READ.GMAT, READ.RHS, READ.index_Q, READ.index_Q_i, READ.index_Q_ij, \
		READ.index_T_i, READ.index_T_ij, READ.index_T_ijk);

	return &READ;
}

// Solve the EIM matrix:
Reading_inputs* solution_EIM(Reading_inputs& READ)
{
	Solve_BC SOL; 
	SOL.Start_solve(READ.nsolve, READ.num, READ.HMAT, READ.GMAT, READ.RHS, READ.Heat_source, READ.U);

	return &READ;
}

// Postprocess thermal fields:
Reading_inputs* Post_process_thermal(Reading_inputs& READ)
{
	EIM_integrals EIM; 

	EIM.post_eigen(READ.nsolve, READ.nump, READ.Points, READ.num, READ.eigen_point, READ.radius, READ.U, \
		READ.Heat_source, READ.index_Q, READ.index_Q_i, READ.index_Q_ij, READ.index_T_i, READ.index_T_ij, READ.index_T_ijk, \
		READ.temp, READ.flux);

	return &READ;
}