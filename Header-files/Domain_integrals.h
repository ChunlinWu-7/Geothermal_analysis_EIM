#pragma once

#include <stdio.h>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;
using namespace boost::math::differentiation;
using namespace boost::math;
using float_50 = boost::multiprecision::cpp_bin_float_100;

class DIE {

public:

	/*
		Notice that: 
			(i) D00_heat_source can be directly used for temperature calculation;
			(ii) D10_heat can be directly used for temperature gradient calculation;
	*/

	// construct Eshelby's tensor for six equivalent conditions
	void D00_heat(int nsolve, double* x, double* x_p, double radius, complex<double>& Ms, complex<double>* Mp, complex<double>** Mpq, \
		complex<double>* J_i, complex<double>** Jp_i, complex<double>*** Jpq_i);

	void D10_heat(int nsolve, double* x, double* x_p, double radius, complex<double>* M_i, complex<double>** Mp_i, complex<double>*** Mpq_i, \
		complex<double>** J_ij, complex<double>*** Jp_ij, complex<double>**** Jpq_ij);

	void D20_heat(int nsolve, double* x, double* x_p, double radius, complex<double>** M_ij, complex<double>*** Mp_ij, complex<double>**** Mpq_ij, \
		complex<double>*** J_ijk, complex<double>**** Jp_ijk, complex<double>***** Jpq_ijk);

	void D30_heat(int nsolve, double* x, double* x_p, double radius, complex<double>*** M_ijk, complex<double>**** Mp_ijk, complex<double>***** Mpq_ijk, \
		complex<double>**** J_ijkl, complex<double>***** Jp_ijkl, complex<double>****** Jpq_ijkl);

private:

	/*
		(i) For semi-infinite case: source and field are always in the same phase;
		(ii) For full-space case: full_flag = 0: full-space or full_flag = -1 semi-infinite case; 
		(iii) z_flag : check whether it is zero; r_flag: interior or exterior
		(iv) a: radius
		(v): beta is the frequency \sqrt{\omega / (2.0 * \alpha)} * (1.0 + 1.0i)
		(vi) xi is the frequency \sqrt{\omega / (2.0 * \alpha)}
	*/

	int dist_flag, r_flag; double a; complex<double> beta;

	// distance and dist groups
	double dist, dist_1, dist_3, dist_5, dist_7, dist_9, dist_11;
	double dist_ori, dist_ori_1, dist_ori_3, dist_ori_5, dist_ori_7, dist_ori_9, dist_ori_11;
	double* difx, * difx_ori;

	// start initialization of the problem!
	void initial_distance(int nsolve, double* x, double* x_p);

	// arrays to store for partial differentiation w.r.t x_{i}
	complex<double>* M0_r_ori, * M1_r_ori, * M2_r_ori;
	complex<double>* M0_r, * M1_r, * M2_r;

	void Store_derivs(int nsolve, double* x, double* x_p);

	/*
		Mfunctions: M0_ori and M0
	*/
	complex<double> M0(); complex<double> M0_1(int i);
	complex<double> M0_2(int i, int j); complex<double> M0_3(int i, int j, int k);
	complex<double> M0_4(int i, int j, int k, int l);

	complex<double> M0_ori(); complex<double> M0_ori_1(int i);
	complex<double> M0_ori_2(int i, int j); complex<double> M0_ori_3(int i, int j, int k);
	complex<double> M0_ori_4(int i, int j, int k, int l);

	/*
		Mfunctions: M1_ori and M1
	*/
	complex<double> M1(); complex<double> M1_1(int i);
	complex<double> M1_2(int i, int j); complex<double> M1_3(int i, int j, int k);
	complex<double> M1_4(int i, int j, int k, int l); complex<double> M1_5(int i, int j, int k, int l, int m);
	complex<double> M1_6(int i, int j, int k, int l, int m, int n);

	complex<double> M1_ori(); complex<double> M1_ori_1(int i);
	complex<double> M1_ori_2(int i, int j); complex<double> M1_ori_3(int i, int j, int k);
	complex<double> M1_ori_4(int i, int j, int k, int l); complex<double> M1_ori_5(int i, int j, int k, int l, int m);
	complex<double> M1_ori_6(int i, int j, int k, int l, int m, int n);

	/*
		Mfunctions: M2_ori and M2
	*/
	complex<double> M2(); complex<double> M2_1(int i);
	complex<double> M2_2(int i, int j); complex<double> M2_3(int i, int j, int k);
	complex<double> M2_4(int i, int j, int k, int l); complex<double> M2_5(int i, int j, int k, int l, int m);
	complex<double> M2_6(int i, int j, int k, int l, int m, int n);

	complex<double> M2_ori(); complex<double> M2_ori_1(int i);
	complex<double> M2_ori_2(int i, int j); complex<double> M2_ori_3(int i, int j, int k);
	complex<double> M2_ori_4(int i, int j, int k, int l); complex<double> M2_ori_5(int i, int j, int k, int l, int m);
	complex<double> M2_ori_6(int i, int j, int k, int l, int m, int n);

	/*
		Real and Imaginary part of M0
		Only when z_flag != 0

		There is no need differentiating ori and imag as r is a parameter!
		The following functions use f_m directly from eyemat!
	*/
	template <typename R> promote<R> M0_real(const R& r);
	template <typename R> promote<R> M0_imag(const R& r);
	/*
		Real and Imaginary part of M1
	*/
	template <typename R> promote<R> M1_real(const R& r);
	template <typename R> promote<R> M1_imag(const R& r);
	/*
		Real and Imaginary part of M2
	*/
	template <typename R> promote<R> M2_real(const R& r);
	template <typename R> promote<R> M2_imag(const R& r);


	/*
		Derivatives of \psi and \phi
		: including ori and imag
	*/
	double psi(); double psi_1(int i);
	double psi_2(int i, int j); double psi_3(int i, int j, int k);
	double psi_4(int i, int j, int k, int l); double psi_5(int i, int j, int k, int l, int m);
	double psi_6(int i, int j, int k, int l, int m, int n);

	double psi_ori(); double psi_ori_1(int i);
	double psi_ori_2(int i, int j); double psi_ori_3(int i, int j, int k);
	double psi_ori_4(int i, int j, int k, int l); double psi_ori_5(int i, int j, int k, int l, int m);
	double psi_ori_6(int i, int j, int k, int l, int m, int n);

	// Since psi is dependent on phi()
	double phi(); double phi_1(int i);
	double phi_2(int i, int j); double phi_3(int i, int j, int k);
	double phi_4(int i, int j, int k, int l); double phi_5(int i, int j, int k, int l, int m);

	double phi_ori(); double phi_ori_1(int i);
	double phi_ori_2(int i, int j); double phi_ori_3(int i, int j, int k);
	double phi_ori_4(int i, int j, int k, int l); double phi_ori_5(int i, int j, int k, int l, int m);

	/*
		Domain integrals of Helmholtz's potential functions: imaged
	*/
	complex<double> Helm_0(); complex<double> Helm_1(int i);
	complex<double> Helm_2(int i, int j); complex<double> Helm_3(int i, int j, int k);
	complex<double> Helm_4(int i, int j, int k, int l);

	/*
		Domain integrals of Helmholtz's potential functions: original
	*/
	complex<double> Helm_ori_0(); complex<double> Helm_ori_1(int i);
	complex<double> Helm_ori_2(int i, int j); complex<double> Helm_ori_3(int i, int j, int k);
	complex<double> Helm_ori_4(int i, int j, int k, int l);

	/*
		Domain integrals of Helmholtz's potential function I:imaged
	*/
	complex<double> Helm_1_0(int p); complex<double> Helm_1_1(int p, int i);
	complex<double> Helm_1_2(int p, int i, int j); complex<double> Helm_1_3(int p, int i, int j, int k);
	complex<double> Helm_1_4(int p, int i, int j, int k, int l);

	/*
		Domain integrals of Helmholtz's potential function I: original
	*/
	complex<double> Helm_ori_1_0(int p); complex<double> Helm_ori_1_1(int p, int i);
	complex<double> Helm_ori_1_2(int p, int i, int j); complex<double> Helm_ori_1_3(int p, int i, int j, int k);
	complex<double> Helm_ori_1_4(int p, int i, int j, int k, int l);

	/*
		Domain integrals of Helmholtz's potential function II: imaged
	*/
	complex<double> Helm_2_0(int p, int q); complex<double> Helm_2_1(int p, int q, int i);
	complex<double> Helm_2_2(int p, int q, int i, int j); complex<double> Helm_2_3(int p, int q, int i, int j, int k);
	complex<double> Helm_2_4(int p, int q, int i, int j, int k, int l);

	/*
		Domain integrals of Helmholtz's potential function II: original
	*/
	complex<double> Helm_ori_2_0(int p, int q); complex<double> Helm_ori_2_1(int p, int q, int i);
	complex<double> Helm_ori_2_2(int p, int q, int i, int j); complex<double> Helm_ori_2_3(int p, int q, int i, int j, int k);
	complex<double> Helm_ori_2_4(int p, int q, int i, int j, int k, int l);

	/*
		Assign domain integrals to Integrated Green's function!

		Initialization and destroy of arrays
	*/
	// temperature
	complex<double> G_0, * Gp_0, ** Gpq_0; complex<double>* G_1, ** Gp_1, *** Gpq_1;

	// first order derivative
	complex<double>** G_2, *** Gp_2, **** Gpq_2;

	// second order derivative
	complex<double>*** G_3, **** Gp_3, ***** Gpq_3;

	// third order derivative
	complex<double>**** G_4, ***** Gp_4, ****** Gpq_4;

	void destroy_0(int nsolve); void destroy_1(int nsolve); void destroy_2(); void destroy_3();
	void inil_0(int nsolve); void inil_1(int nsolve); void inil_2(); void inil_3();

	// temperature:
	void Int_Gfunc_temp(int nsolve, double* x, double* x_p);

	void Int_Gfunc_flux(int nsolve, double* x, double* x_p);

	void Int_Gfunc_flux_grad1(int nsolve, double* x, double* x_p);

	void Int_Gfunc_flux_grad2(int nsolve, double* x, double* x_p);

	void clear_all();

};


