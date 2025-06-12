
#include "Equivalent_conditions.h"
#include <fstream>

void EIM_integrals::addInclusion(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT\
	, Ref<MatrixXcd> GMAT, Ref<VectorXcd> RHS, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk)
{
	add_equiv_0(nsolve, num, eigen_point, radius, eigen_mat, HMAT, GMAT, index_B\
		, index_B_i, index_B_ij, index_E_i, index_E_ij, index_E_ijk);

	if (nsolve != 4) {
		add_equiv_1(nsolve, num, eigen_point, radius, eigen_mat, HMAT, GMAT, index_B\
			, index_B_i, index_B_ij, index_E_i, index_E_ij, index_E_ijk);
		if (nsolve != 16) {
			add_equiv_2(nsolve, num, eigen_point, radius, eigen_mat, HMAT, GMAT, index_B\
				, index_B_i, index_B_ij, index_E_i, index_E_ij, index_E_ijk);
		}
	}

	// call functions to create RHS for temperature gradient
	add_temp_grad(annual_flag, nsolve, num, eigen_point, eigen_mat, RHS, index_B, index_B_i, index_B_ij, index_E_i, index_E_ij, index_E_ijk);

}

void EIM_integrals::add_equiv_0(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT\
	, Ref<MatrixXcd> GMAT, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk)
{
	complex<double>** A = new complex<double>*[nsolve * num];

	complex<double>** B = new complex<double> *[nsolve * num];

	for (int i = 0; i < nsolve * num; i++) {
		A[i] = new complex<double>[nsolve * num];

		B[i] = new complex<double>[nsolve * num];
	}

	for (int i = 0; i < nsolve * num; i++) {
		for (int j = 0; j < nsolve * num; j++) {
			A[i][j] = 0.0;
 			
			B[i][j] = 0.0;
		}
	}

# pragma omp parallel shared(A,num,radius,index_E_i,index_E_ij,index_E_ijk)

	{
		int t = 0;
		int m, n, p, q, h;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for schedule(dynamic)
		// s: build point; 
		for (int s = 0; s < num; s++) {
			//h: source point
			double K_p, Cp_p; double* x = new double[3];

			x[0] = eigen_point(s, 0); x[1] = eigen_point(s, 1); x[2] = eigen_point(s, 2);
			K_p = eigen_mat(s, 0); Cp_p = eigen_mat(s, 1); 

			for (h = 0; h < num; h++) {
				if (h == s) {
					sym = 1.0;
				}
				else {
					sym = 0.0;
				}

				double* x_p = new double[3];
				x_p[0] = eigen_point(h, 0);
				x_p[1] = eigen_point(h, 1);
				x_p[2] = eigen_point(h, 2);

				// heat flux
				complex<double>* M_i, ** Mp_i, *** Mpq_i; complex<double>** J_ij, *** Jp_ij, **** Jpq_ij;

				// temperature 
				complex<double> M_tensor, * Mp, ** Mpq; complex<double>* J_i, ** Jp_i, *** Jpq_i;
			
				inil_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);

				DIE SPE; SPE.D10_heat(nsolve, x, x_p, radius[h], M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij);
				SPE.D00_heat(nsolve, x, x_p, radius[h], M_tensor, Mp, Mpq, J_i, Jp_i, Jpq_i);

				// ETG related equivalent conditions
				for (int i = 0; i < 3; i++) {

					A[index_E_i[3 * s + i]][index_B[h]] = -(K_0 - K_p) * M_i[i];

					B[index_E_i[3 * s + i]][index_B[h]] = -(K_0 - K_p) * M_i[i];

					if (nsolve != 4) {
						for (int p = 0; p < 3; p++) {
							A[index_E_i[3 * s + i]][index_B_i[3 * h + p]] = -(K_0 - K_p) * Mp_i[p][i];
						}

						if (nsolve != 16) {
							for (int p = 0; p < 3; p++) {
								for (int q = 0; q < 3; q++) {
									A[index_E_i[3 * s + i]][index_B_ij[3 * h + p][q]] = -(K_0 - K_p) * Mpq_i[p][q][i];
								}
							}
						}
					}

					for (int j = 0; j < 3; j++) {

						complex<double> D_ij = -(K_0 - K_p) * J_ij[j][i];

						A[index_E_i[3 * s + i]][index_E_i[3 * h + j]] = D_ij + sym * K_0 * d[i][j];

						if (nsolve != 4) {
							for (p = 0; p < 3; p++) {

								A[index_E_i[3 * s + i]][index_E_ij[3 * h + j][p]] = -(K_0 - K_p) * Jp_ij[p][j][i];
							}
							if (nsolve != 16) {
								for (p = 0; p < 3; p++) {
									for (q = 0; q < 3; q++) {

										A[index_E_i[3 * s + i]][index_E_ijk[3 * h + j][p][q]] = -(K_0 - K_p) * Jpq_ij[p][q][j][i];

									}
								}
							}
						}

					}
				}

				// eigen heat source related equivalent conditions
				A[index_B[s]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_tensor + sym;

				B[index_B[s]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_tensor;

				for (int i = 0; i < 3; i++) {
					A[index_B[s]][index_E_i[3 * h + i]] = 1.0i * omega * (Cp_0 - Cp_p) * J_i[i];
				}

				if (nsolve != 4) {
					for (int p = 0; p < 3; p++) {
						A[index_B[s]][index_B_i[3 * h + p]] = 1.0i * omega * (Cp_0 - Cp_p) * Mp[p];

						for (int i = 0; i < 3; i++) {
							A[index_B[s]][index_E_ij[3 * h + i][p]] = 1.0i * omega * (Cp_0 - Cp_p) * Jp_i[p][i];
						}
					}

					if (nsolve != 16) {
						for (int p = 0; p < 3; p++) {
							for (int q = 0; q < 3; q++) {
								A[index_B[s]][index_B_ij[3 * h + p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Mpq[p][q];

								for (int i = 0; i < 3; i++) {
									A[index_B[s]][index_E_ijk[3 * h + i][p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Jpq_i[p][q][i];
								}

							}
						}
					}
				}


				destroy_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);
				delete[] x_p;
			}
			//cout << "s =" << s << " ";
			delete[] x;
		}

	}
	
	for (int i = 0; i < nsolve * num; i++)
	{
		for (int j = 0; j < nsolve * num; j++)
		{
			HMAT(i, j) = A[i][j];
			
			GMAT(i, j) = B[i][j];
			//cout << B[i][j] << " ";
		}
		//cout << endl;
	}

	for (int i = 0; i < nsolve * num; i++) {
		delete A[i];
		delete B[i];
	}
	delete[] A; delete[] B;
}

void EIM_integrals::add_equiv_1(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT\
	, Ref<MatrixXcd> GMAT, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk)
{
	complex<double>** A = new complex<double>*[nsolve * num];

	complex<double>** B = new complex<double> *[nsolve * num];

	for (int i = 0; i < nsolve * num; i++) {
		A[i] = new complex<double>[nsolve * num];

		B[i] = new complex<double>[nsolve * num];
	}

	for (int i = 0; i < nsolve * num; i++) {
		for (int j = 0; j < nsolve * num; j++) {
			A[i][j] = 0.0;

			B[i][j] = 0.0;
		}
	}

# pragma omp parallel shared(A,num,radius,index_E_i,index_E_ij,index_E_ijk)

	{
		int t = 0;
		int m, n, p, q, h;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for schedule(dynamic)
		// s: build point; 
		for (int s = 0; s < num; s++) {
			//h: source point
			double K_p, Cp_p; double* x = new double[3];

			x[0] = eigen_point(s, 0); x[1] = eigen_point(s, 1); x[2] = eigen_point(s, 2);
			K_p = eigen_mat(s, 0); Cp_p = eigen_mat(s, 1);

			for (h = 0; h < num; h++) {
				if (h == s) {
					sym = 1.0;
				}
				else {
					sym = 0.0;
				}

				double* x_p = new double[3];
				x_p[0] = eigen_point(h, 0);
				x_p[1] = eigen_point(h, 1);
				x_p[2] = eigen_point(h, 2);

				// heat flux
				complex<double>** M_ij, *** Mp_ij, **** Mpq_ij; complex<double>*** J_ijk, **** Jp_ijk, ***** Jpq_ijk;

				// temperature 
				complex<double> *M_i, ** Mp_i, *** Mpq_i; complex<double>** J_ij, *** Jp_ij, **** Jpq_ij;

				inil_1(M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk, M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij);
				
				DIE SPE; SPE.D20_heat(nsolve, x, x_p, radius[h], M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk);\
				SPE.D10_heat(nsolve, x, x_p, radius[h], M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij);


				// ETG - related equivalence conditions
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						A[index_E_ij[3 * s + i][j]][index_B[h]] = -(K_0 - K_p) * M_ij[i][j];

						B[index_E_ij[3 * s + i][j]][index_B[h]] = -(K_0 - K_p) * M_ij[i][j];
						if (nsolve != 4) {
							for (int p = 0; p < 3; p++) {
								A[index_E_ij[3 * s + i][j]][index_B_i[3 * h + p]] = -(K_0 - K_p) * Mp_ij[p][i][j];
							}

							if (nsolve != 16) {
								for (int p = 0; p < 3; p++) {
									for (int q = 0; q < 3; q++) {
										A[index_E_ij[3 * s + i][j]][index_B_ij[3 * h + p][q]] = -(K_0 - K_p) * Mpq_ij[p][q][i][j];
									}
								}
							}

						}
					}

					for (int j = 0; j < 3; j++) {

						for (int k = 0; k < 3; k++) {

							A[index_E_ij[3 * s + i][j]][index_E_i[3 * h + k]] = -(K_0 - K_p) * J_ijk[k][i][j];

							if (nsolve != 4) {
								for (p = 0; p < 3; p++) {

									A[index_E_ij[3 * s + i][j]][index_E_ij[3 * h + k][p]] = -(K_0 - K_p) * Jp_ijk[p][k][i][j] + K_0 * sym * d[i][k] * d[j][p];
								}
								if (nsolve != 16) {
									for (p = 0; p < 3; p++) {
										for (q = 0; q < 3; q++) {

											A[index_E_ij[3 * s + i][j]][index_E_ijk[3 * h + k][p][q]] = -(K_0 - K_p) * Jpq_ijk[p][q][k][i][j];

										}
									}
								}
							}

						}
					}
				}

				// eigen-heat-source related equivalence conditions
				for (int i = 0; i < 3; i++) {
					A[index_B_i[3 * s + i]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_i[i];

					B[index_B_i[3 * s + i]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_i[i];
					for (int j = 0; j < 3; j++) {
						A[index_B_i[3 * s + i]][index_E_i[3 * h + j]] = 1.0i * omega * (Cp_0 - Cp_p) * J_ij[j][i];
					}
				}

				if (nsolve != 4) {

					for (int p = 0; p < 3; p++) {
						for (int i = 0; i < 3; i++) {
							A[index_B_i[3 * s + i]][index_B_i[3 * h + p]] = 1.0i * omega * (Cp_0 - Cp_p) * Mp_i[p][i] + sym * d[i][p];
							for (int j = 0; j < 3; j++) {
								A[index_B_i[3 * s + i]][index_E_ij[3 * h + j][p]] = 1.0i * omega * (Cp_0 - Cp_p) * Jp_ij[p][j][i];
							}
						}
					}

					if (nsolve != 16) {

						for (int p = 0; p < 3; p++) {
							for (int q = 0; q < 3; q++) {
								for (int i = 0; i < 3; i++) {
									A[index_B_i[3 * s + i]][index_B_ij[3 * h + p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Mpq_i[p][q][i];
									for (int j = 0; j < 3; j++) {
										A[index_B_i[3 * s + i]][index_E_ijk[3 * h + j][p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Jpq_ij[p][q][j][i];
									}
								}
							}
						}

					}

				}

				destroy_1(M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk, M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij);
				delete[] x_p;
			}
			//cout << "s =" << s << " ";
			delete[] x;
		}

	}

	for (int i = 0; i < nsolve * num; i++)
	{
		for (int j = 0; j < nsolve * num; j++)
		{
			HMAT(i, j) += A[i][j];

			GMAT(i, j) += B[i][j];
			//cout << HMAT(i, j) << " ";
		}
		//cout << endl;
	}

	for (int i = 0; i < nsolve * num; i++) {
		delete A[i];
		delete B[i];
	}
	delete[] A; delete[] B;
}

void EIM_integrals::add_equiv_2(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT\
	, Ref<MatrixXcd> GMAT, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk)
{
	complex<double>** A = new complex<double>*[nsolve * num];

	complex<double>** B = new complex<double> *[nsolve * num];

	for (int i = 0; i < nsolve * num; i++) {
		A[i] = new complex<double>[nsolve * num];

		B[i] = new complex<double>[nsolve * num];
	}

	for (int i = 0; i < nsolve * num; i++) {
		for (int j = 0; j < nsolve * num; j++) {
			A[i][j] = 0.0;

			B[i][j] = 0.0;
		}
	}

# pragma omp parallel shared(A,num,radius,index_E_i,index_E_ij,index_E_ijk)

	{
		int t = 0;
		int m, n, p, q, h;
		double sym = 0.0;

		////Assembly Matrix Tensor///////

#	pragma omp for schedule(dynamic)
		// s: build point; 
		for (int s = 0; s < num; s++) {
			//h: source point
			double K_p, Cp_p; double* x = new double[3];

			x[0] = eigen_point(s, 0); x[1] = eigen_point(s, 1); x[2] = eigen_point(s, 2);
			K_p = eigen_mat(s, 0); Cp_p = eigen_mat(s, 1);

			for (h = 0; h < num; h++) {
				if (h == s) {
					sym = 1.0;
				}
				else {
					sym = 0.0;
				}

				double* x_p = new double[3];
				x_p[0] = eigen_point(h, 0);
				x_p[1] = eigen_point(h, 1);
				x_p[2] = eigen_point(h, 2);

				// temperature 
				complex<double>*** M_ijk, **** Mp_ijk, ***** Mpq_ijk; complex<double>**** J_ijkl, ***** Jp_ijkl, ****** Jpq_ijkl;

				// heat flux
				complex<double>** M_ij, *** Mp_ij, **** Mpq_ij; complex<double>*** J_ijk, **** Jp_ijk, ***** Jpq_ijk;

				inil_2(M_ijk, Mp_ijk, Mpq_ijk, J_ijkl, Jp_ijkl, Jpq_ijkl, M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk);

				DIE SPE; SPE.D30_heat(nsolve, x, x_p, radius[h], M_ijk, Mp_ijk, Mpq_ijk, J_ijkl, Jp_ijkl, Jpq_ijkl);
				SPE.D20_heat(nsolve, x, x_p, radius[h], M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk);

				// ETG - related equivalence conditions
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 3; k++) {
							A[index_E_ijk[3 * s + i][j][k]][index_B[h]] = -0.5 * (K_0 - K_p) * M_ijk[i][j][k];

							B[index_E_ijk[3 * s + i][j][k]][index_B[h]] = -0.5 * (K_0 - K_p) * M_ijk[i][j][k];
							if (nsolve != 4) {
								for (int p = 0; p < 3; p++) {
									A[index_E_ijk[3 * s + i][j][k]][index_B_i[3 * h + p]] = -0.5 * (K_0 - K_p) * Mp_ijk[p][i][j][k];
								}

								if (nsolve != 16) {
									for (int p = 0; p < 3; p++) {
										for (int q = 0; q < 3; q++) {
											A[index_E_ijk[3 * s + i][j][k]][index_B_ij[3 * h + p][q]] = -0.5 * (K_0 - K_p) * Mpq_ijk[p][q][i][j][k];
										}
									}
								}

							}
						}
					}

					for (int j = 0; j < 3; j++) {

						for (int k = 0; k < 3; k++) {
							for (int l = 0; l < 3; l++) {
								A[index_E_ijk[3 * s + i][j][k]][index_E_i[3 * h + l]] = -0.5 * (K_0 - K_p) * J_ijkl[l][i][j][k];

								if (nsolve != 4) {
									for (p = 0; p < 3; p++) {

										A[index_E_ijk[3 * s + i][j][k]][index_E_ij[3 * h + l][p]] = -0.5 * (K_0 - K_p) * Jp_ijkl[p][l][i][j][k];
									}
									if (nsolve != 16) {
										for (p = 0; p < 3; p++) {
											for (q = 0; q < 3; q++) {
												A[index_E_ijk[3 * s + i][j][k]][index_E_ijk[3 * h + l][p][q]] = (-0.5 * (K_0 - K_p) * Jpq_ijkl[p][q][l][i][j][k] + K_0 * sym * d[i][l] * d[j][p] * d[k][q]);
											}
										}
									}
								}
							}
						}
					}
				}

				// eigen-heat-source equivalence conditions
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						A[index_B_ij[3 * s + i][j]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_ij[i][j];

						B[index_B_ij[3 * s + i][j]][index_B[h]] = 1.0i * omega * (Cp_0 - Cp_p) * M_ij[i][j];
						for (int k = 0; k < 3; k++) {
							A[index_B_ij[3 * s + i][j]][index_E_i[3 * h + k]] = 1.0i * omega * (Cp_0 - Cp_p) * J_ijk[k][i][j];
						}
					}
				}

				if (nsolve != 4) {

					for (int p = 0; p < 3; p++) {
						for (int i = 0; i < 3; i++) {
							for (int j = 0; j < 3; j++) {
								A[index_B_ij[3 * s + i][j]][index_B_i[3 * h + p]] = 1.0i * omega * (Cp_0 - Cp_p) * Mp_ij[p][i][j];
								for (int k = 0; k < 3; k++) {
									A[index_B_ij[3 * s + i][j]][index_E_ij[3 * h + k][p]] = 1.0i * omega * (Cp_0 - Cp_p) * Jp_ijk[p][k][i][j];
								}
							}
						}
					}

					if (nsolve != 16) {

						for (int p = 0; p < 3; p++) {
							for (int q = 0; q < 3; q++) {
								for (int i = 0; i < 3; i++) {
									for (int j = 0; j < 3; j++) {
										A[index_B_ij[3 * s + i][j]][index_B_ij[3 * h + p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Mpq_ij[p][q][i][j] 
											+ sym * (d[i][p] * d[j][q] + d[i][q] * d[j][p]);
										for (int k = 0; k < 3; k++) {
											A[index_B_ij[3 * s + i][j]][index_E_ijk[3 * h + k][p][q]] = 1.0i * omega * (Cp_0 - Cp_p) * Jpq_ijk[p][q][k][i][j];
										}
									}
								}
							}
						}

					}

				}

				destroy_2(M_ijk, Mp_ijk, Mpq_ijk, J_ijkl, Jp_ijkl, Jpq_ijkl, M_ij, Mp_ij, Mpq_ij, J_ijk, Jp_ijk, Jpq_ijk);
				delete[] x_p;
			}
			cout << "s =" << s << " ";
			delete[] x;
		}

	}

	for (int i = 0; i < nsolve * num; i++)
	{
		for (int j = 0; j < nsolve * num; j++)
		{
			HMAT(i, j) += A[i][j];
			GMAT(i, j) += B[i][j];
		}
	}

	for (int i = 0; i < nsolve * num; i++) {
		delete A[i]; delete B[i];
	}
	delete[] A; delete[] B;
}

// Postprocess, displacement and stress
void EIM_integrals::post_eigen(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, \
	Ref<VectorXd> Heat_source, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>* temp, complex<double>** flux)
{
	/*
		Postprocess of temperature and heat flux
	*/

	post_eigen_temp(nsolve, num_post, Points, num, eigen_point, radius, U, Heat_source, index_B, index_B_i, index_B_ij\
		, index_E_i, index_E_ij, index_E_ijk, temp);

	post_eigen_flux(nsolve, num_post, Points, num, eigen_point, radius, U, Heat_source, index_B, index_B_i, index_B_ij\
		, index_E_i, index_E_ij, index_E_ijk, flux);

	/*
		post_temp.txt: postprocess temperature, including disturbances of tanks;
		post_flux.txt: postprocess heat flux;
		post_temp_ori.txt: temperature profile (the unperturbed solution)
	*/

	ofstream myfile_T, myfile_F, myfile_T_ori;

	myfile_T.open("post_temp.txt"); myfile_F.open("post_flux.txt"); myfile_T_ori.open("post_temp_ori.txt");

	/*
		sd refers to division: (1) see several months, sd = 6.0; (2) see daily variation: sd = 364.0;
	*/

	double sd = sdd;

	for (int i = 0; i < num_post; i++) {

		for (int j = 0; j < sd; j++) {

			myfile_T << (temp[i] * exp(-1.0i * (double(j) - sd / 2.0) * 2.0 * pi / sd)).real() + 0.0\
				+ (13.89 * exp((1.0 - 1.0i) * f_m * Points(i, 2)) \
					* exp(-1.0i * (double(j) - sd / 2.0) * 2.0 * pi / sd)).real() + 12.78 << '\t';

			myfile_T_ori << (13.89 * exp((1.0 - 1.0i) * f_m * Points(i, 2)) \
					* exp(-1.0i * (double(j) - sd / 2.0) * 2.0 * pi / sd)).real() + 12.78 << '\t';

			for (int s = 2; s < 3; s++) {
				myfile_F << (flux[i][s] * exp(-1.0i * (double(j) - sd / 2.0) * 2.0 * pi / sd)).real() +
					(-K_0 * 13.89 * d[s][2] * exp((1.0 - 1.0i) * Points(i, 2) * f_m) * (1.0 - 1.0i) \
						* f_m * exp(-1.0i * (double(j) - sd / 2.0) * 2.0 * pi / sd)).real()\
					<< '\t';
			}

		}

		myfile_T << endl; myfile_F << endl; myfile_T_ori << endl;
	}

	myfile_T.close(); myfile_F.close(); myfile_T_ori.close();
}

void EIM_integrals::post_eigen_temp(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, \
	Ref<VectorXd> Heat_source, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>* temp)
{
	// retrieve
	complex<double>* Eigen = new complex<double>[num * nsolve];

	for (int i = 0; i < num * nsolve; i++) {
		Eigen[i] = Heat_source(i) + U(i);
		
	}

#pragma omp parallel shared(Points, eigen_point, radius)

	{
		int s, h;
		double* x = new double[3]; double* x_p = new double[3];

		// For flux equivalence:time integral
		complex<double>* M_i, ** Mp_i, *** Mpq_i, ** J_ij, *** Jp_ij, **** Jpq_ij;
		// For temperature equivalence: without time integral: Non-use in postprocess for gradient!
		complex<double> M_tensor, * Mp, ** Mpq, * J_i, ** Jp_i, *** Jpq_i;

		inil_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);

#pragma omp for schedule(dynamic)

		for (s = 0; s < num_post; s++) {

			x[0] = Points(s, 0); x[1] = Points(s, 1); x[2] = Points(s, 2);

			for (int h = 0; h < num; h++) {
				x_p[0] = eigen_point(h, 0); x_p[1] = eigen_point(h, 1); x_p[2] = eigen_point(h, 2);
				DIE SPE;
				SPE.D00_heat(nsolve, x, x_p, radius[h], M_tensor, Mp, Mpq, J_i, Jp_i, Jpq_i);

				// by ETG

				for (int j = 0; j < 3; j++) {
					temp[s] = temp[s] + J_i[j] * Eigen[index_E_i[3 * h + j]];
					if (nsolve != 4) {
						for (int p = 0; p < 3; p++) {
							temp[s] = temp[s] + Jp_i[p][j] * Eigen[index_E_ij[3 * h + j][p]];
							if (nsolve != 16) {
								for (int q = 0; q < 3; q++) {
									temp[s] = temp[s] + Jpq_i[p][q][j] * Eigen[index_E_ijk[3 * h + j][p][q]];
								}
							}
						}
					}
				}

				// by eigen-heat-source
				temp[s] = temp[s] + M_tensor * Eigen[index_B[h]];
				if (nsolve != 4) {
					for (int p = 0; p < 3; p++) {
						temp[s] = temp[s] + Mp[p] * Eigen[index_B_i[3 * h + p]];
						if (nsolve != 16) {
							for (int q = 0; q < 3; q++) {
								temp[s] = temp[s] + Mpq[p][q] * Eigen[index_B_ij[3 * h + p][q]];
							}
						}
					}
				}


			}

		}

		destroy_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);
		delete[] x; delete[] x_p; 
	}
}

void EIM_integrals::post_eigen_flux(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, \
	Ref<VectorXd> Heat_source, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>** flux)
{
	complex<double> * Eigen = new complex<double>[num * nsolve];

	for (int i = 0; i < num * nsolve; i++) {
		Eigen[i] = U(i) + Heat_source(i);
	}

#pragma omp parallel shared(Points, eigen_point, radius)

	{
		// judge whether in the particles!
		double** ep_s = new double* [num_post];
		// judge in or outside
		int* test = new int[num_post];
		for (int i = 0; i < num_post; i++) {
			test[i] = -1;
			ep_s[i] = new double[3];			
		}

		for (int s = 0; s < num_post; s++) {
			for (int h = 0; h < num; h++) {
				double dist = 0.0;
				for (int i = 0; i < 3; i++) {
					dist += pow(Points(s, i) - eigen_point(h, i), 2.0);
				}
				dist = sqrt(dist);
				if (dist <= 1.0 * radius(h)) {
					test[s] = h;
				}
			}

		}

		double* x = new double[3]; double* x_p = new double[3];

		// For flux equivalence:time integral
		complex<double>* M_i, ** Mp_i, *** Mpq_i, ** J_ij, *** Jp_ij, **** Jpq_ij;
		// For temperature equivalence: without time integral: Non-use in postprocess for gradient!
		complex<double> M_tensor, * Mp, ** Mpq, * J_i, ** Jp_i, *** Jpq_i;

		inil_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);

#pragma omp for schedule(dynamic)

		for (int s = 0; s < num_post; s++) {

			complex<double>* ETG = new complex<double>[3];
			for (int i = 0; i < 3; i++) {
				ETG[i] = 0.0;
			}

			x[0] = Points(s, 0); x[1] = Points(s, 1); x[2] = Points(s, 2);

			if (test[s] != -1) {

				for (int i = 0; i < 3; i++) {

					ETG[i] = Eigen[index_E_i[3 * test[s] + i]];

					if (nsolve != 4) {
						for (int k = 0; k < 3; k++) {
							ETG[i] = ETG[i] + Eigen[index_E_ij[3 * test[s] + i][k]] * (x[k] - eigen_point(test[s], k));
						}

						if (nsolve != 16) {
							for (int k = 0; k < 3; k++) {
								for (int l = 0; l < 3; l++) {
									ETG[i] = ETG[i] + Eigen[index_E_ijk[3 * test[s] + i][k][l]] * (x[k] - eigen_point(test[s], k)) * (x[l] - eigen_point(test[s], l));
								}
							}

						}
					}

				}
			}

			for (int h = 0; h < num; h++) {
				x_p[0] = eigen_point(h, 0); x_p[1] = eigen_point(h, 1); x_p[2] = eigen_point(h, 2);
				DIE SPE;
				SPE.D10_heat(nsolve, x, x_p, radius[h], M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij);

				// by ETG
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						flux[s][i] = flux[s][i] - K_0 * J_ij[j][i] * Eigen[index_E_i[3 * h + j]];
						if (nsolve != 4) {
							for (int p = 0; p < 3; p++) {
								flux[s][i] = flux[s][i] - K_0 * Jp_ij[p][j][i] * Eigen[index_E_ij[3 * h + j][p]];
								if (nsolve != 16) {
									for (int q = 0; q < 3; q++) {
										flux[s][i] = flux[s][i] - K_0 * Jpq_ij[p][q][j][i] * Eigen[index_E_ijk[3 * h + j][p][q]];
									}
								}
							}
						}
					}
				}

				// by eigen-heat-source
				for (int i = 0; i < 3; i++) {					
					flux[s][i] = flux[s][i] - K_0 * M_i[i] * Eigen[index_B[h]];
					if (nsolve != 4) {
						for (int p = 0; p < 3; p++) {
							flux[s][i] = flux[s][i] - K_0 * Mp_i[p][i] * Eigen[index_B_i[3 * h + p]];
							if (nsolve != 16) {
								for (int q = 0; q < 3; q++) {
									flux[s][i] = flux[s][i] - K_0 * Mpq_i[p][q][i] * Eigen[index_B_ij[3 * h + p][q]];
								}
							}
						}
					}
				}

			}

			// Eigen part
			for (int l = 0; l < 3; l++) {
				flux[s][l] = flux[s][l] - (-K_0) * ETG[l];
			}
			delete[] ETG;
		}

		destroy_0(M_i, Mp_i, Mpq_i, J_ij, Jp_ij, Jpq_ij, Mp, Mpq, J_i, Jp_i, Jpq_i);
		delete[] x; delete[] x_p; delete[] test;
	}

}

// add temperature gradient and temperature difference
void EIM_integrals::add_temp_grad(int flag, int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<MatrixXd> eigen_mat, Ref<VectorXcd> RHS, int* index_B, int* index_B_i, int** index_B_ij, \
	int* index_E_i, int** index_E_ij, int*** index_E_ijk)
{
	// assign values to RHS vector
	double parameter = 0.0;
	if (flag == 1) {
		parameter = 13.89;
	}
	else {
		parameter = 7.78 * 0.5;
	}

	// ETG-related problem:
	for (int s = 0; s < num; s++) {
		double z = eigen_point(s, 2);
		double K_s = eigen_mat(s, 0); double Cp_s = eigen_mat(s, 1);
		
		// uniform EHS
		RHS(index_B[s]) = exp((1.0 - 1.0i) * z * f_m) * \
			(-1.0i * omega) * (Cp_0 - Cp_s) * parameter;

		// uniform ETG
		for (int i = 0; i < 3; i++) {
			RHS(index_E_i[3 * s + i]) = d[i][2] * exp((1.0 - 1.0i) * z * f_m) * (1.0 - 1.0i) * f_m\
				* (K_0 - K_s) * parameter;
		}

		if (nsolve != 4) {

			for (int p = 0; p < 3; p++) {
				RHS(index_B_i[3 * s + p]) = d[p][2] * exp((1.0 - 1.0i) * z * f_m) * (1.0 - 1.0i) * f_m\
					* (-1.0i * omega) * (Cp_0 - Cp_s) * parameter;
			}

			for (int p = 0; p < 3; p++) {
				for (int i = 0; i < 3; i++) {
					RHS(index_E_ij[3 * s + i][p]) = d[i][2] * d[p][2] * exp((1.0 - 1.0i) * z * f_m) * -2.0i * f_m * f_m\
						* (K_0 - K_s) * parameter;
				}
			}

			if (nsolve != 16) {

				for (int p = 0; p < 3; p++) {
					for (int q = 0; q < 3; q++) {
						RHS(index_B_ij[3 * s + p][q]) = 0.5 * d[p][2] * d[q][2] * exp((1.0 - 1.0i) * z * f_m) * -2.0i * f_m * f_m\
							* (-1.0i * omega) * (Cp_0 - Cp_s) * parameter;
					}
				}

				for (int p = 0; p < 3; p++) {
					for (int q = 0; q < 3; q++) {
						for (int i = 0; i < 3; i++) {
							RHS(index_E_ijk[3 * s + i][p][q]) = 0.5 * d[i][2] * d[p][2] * d[q][2] * exp((1.0 - 1.0i) * z * f_m) * -2.0i * (1.0 - 1.0i) * f_m * f_m * f_m\
								* (K_0 - K_s) * parameter;
						}
					}
				}

			}

		}

	}



}
