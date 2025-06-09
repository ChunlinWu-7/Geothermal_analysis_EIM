
#include <iostream>
#include <complex>
#include "eigen-3.4.0/Eigen/Dense"

#include "Equivalent_conditions.h"
#include "Solve.h"
#include "Input_parameters.h"

using namespace Eigen; 
using namespace std;


int main()
{
	Reading_inputs READ; READ.Start_read(); 

	EIM_integrals EIM; EIM.addInclusion(READ.nsolve, READ.num, READ.eigen_point, READ.radius, READ.eigen_mat, \
		READ.HMAT, READ.GMAT, READ.RHS, READ.index_Q, READ.index_Q_i, READ.index_Q_ij, READ.index_T_i, READ.index_T_ij, READ.index_T_ijk);

	Solve_BC SOL; SOL.Start_solve(READ.nsolve, READ.num, READ.HMAT, READ.GMAT, READ.RHS, READ.Heat_source, READ.U);

	EIM.post_eigen(READ.nsolve, READ.nump, READ.Points, READ.num, READ.eigen_point, READ.radius, READ.U, \
		READ.Heat_source, READ.index_Q, READ.index_Q_i, READ.index_Q_ij, READ.index_T_i, READ.index_T_ij, READ.index_T_ijk, \
		READ.temp, READ.flux);

	return 0;
}
