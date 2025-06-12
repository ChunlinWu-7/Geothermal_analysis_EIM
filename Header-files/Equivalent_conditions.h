#pragma once
#include "Domain_integrals.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eyemat.h"
#include <complex>

using namespace std;
using namespace Eigen;

class EIM_integrals {
public:
	/*
		addInclusion: arrange coefficient matrices for equivalent conditions:

		(1) add_equiv_0: uniform order; (2) add_equiv_1: linear order; (3) add_equiv_2: quadratic order
	*/

	void addInclusion(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT\
		,Ref<VectorXcd> RHS, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk);

	void add_equiv_0(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT \
		, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk);

	void add_equiv_1(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT\
		, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk);

	void add_equiv_2(int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<MatrixXd> eigen_mat, Ref<MatrixXcd> HMAT, Ref<MatrixXcd> GMAT\
		, int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk);

	// Postprocess, temperature and heat flux:
	void post_eigen(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, Ref<VectorXd> Heat_source,\
		int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>* temp, complex<double>** flux);

	void post_eigen_temp(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, Ref<VectorXd> Heat_source, \
		int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>* temp);

	void post_eigen_flux(int nsolve, int num_post, Ref<MatrixXd> Points, int num, Ref<MatrixXd> eigen_point, Ref<VectorXd> radius, Ref<VectorXcd> U, Ref<VectorXd> Heat_source, \
		int* index_B, int* index_B_i, int** index_B_ij, int* index_E_i, int** index_E_ij, int*** index_E_ijk, complex<double>** flux);

private:

	// temperature gradients by central park data annual:
	void add_temp_grad(int flag, int nsolve, int num, Ref<MatrixXd> eigen_point, Ref<MatrixXd> eigen_mat, Ref<VectorXcd> RHS, int* index_B, int* index_B_i, int** index_B_ij, \
		int* index_E_i, int** index_E_ij, int*** index_E_ijk);

	void inil_0(complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij, \
		complex<double>*& Mp, complex<double>**& Mpq, complex<double>*& J_i, complex<double>**& Jp_i, complex<double>***& Jpq_i);

	void destroy_0(complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij, \
		complex<double>*& Mp, complex<double>**& Mpq, complex<double>*& J_i, complex<double>**& Jp_i, complex<double>***& Jpq_i);

	void inil_1(complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk, \
		complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij);
	void destroy_1(complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk, \
		complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij);

	void inil_2(complex<double>***& M_ijk, complex<double>****& Mp_ijk, complex<double>*****& Mpq_ijk, complex<double>****& J_ijkl, complex<double>*****& Jp_ijkl, complex<double>******& Jpq_ijkl, \
		complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk);
	void destroy_2(complex<double>***& M_ijk, complex<double>****& Mp_ijk, complex<double>*****& Mpq_ijk, complex<double>****& J_ijkl, complex<double>*****& Jp_ijkl, complex<double>******& Jpq_ijkl, \
		complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk);
};