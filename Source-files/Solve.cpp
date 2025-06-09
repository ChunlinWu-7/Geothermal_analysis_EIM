
#include "Solve.h"
#include <iostream>
#include <fstream>

void Solve_BC::Start_solve(int nsolve, int num, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT, Ref<VectorXcd> RHS, Ref<VectorXd> H_source, Ref<VectorXcd> U)
{
	VectorXcd BB = VectorXcd::Zero(num * nsolve);

	for (int i = 0; i < num * nsolve; i++) {
		for (int j = 0; j < num * nsolve; j++) {
			BB(i) = BB(i) - H_source(j) * GMAT(i, j);
		}
		BB(i) = BB(i) + RHS(i);
	}

	/*multi-core solve*/
	Eigen::initParallel(); int core = Eigen::nbThreads();
	Eigen::setNbThreads(core);

	VectorXcd XX = HMAT.fullPivLu().solve(BB);

	ofstream myfile; myfile.open("eigen-field.txt");

	for (int i = 0; i < num * nsolve; i++) {
		U(i) = XX(i);
		myfile << U(i).real() << " " << U(i).imag() << endl;
	}

	myfile.close();

	cout << "residue = " << '\t' << (HMAT * U - BB).norm() << endl;



}

