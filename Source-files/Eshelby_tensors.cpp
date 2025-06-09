#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

void DIE::D00_heat(int nsolve, double* x, double* x_p, double radius, complex<double>& Ms, complex<double>* Mp, complex<double>** Mpq, \
	complex<double>* J_i, complex<double>** Jp_i, complex<double>*** Jpq_i)
{
	// initialization 
	inil_0(nsolve);
	a = radius;
	
	// for temperature
	Int_Gfunc_temp(nsolve, x, x_p);

	// for eigen-temeprature-gradient
	for (int i = 0; i < 3; i++) {
		J_i[i] = G_1[i];
		if (nsolve != 4) {
			for (int p = 0; p < 3; p++) {
				Jp_i[p][i] = Gp_1[p][i];
				if (nsolve != 16) {
					for (int q = 0; q < 3; q++) {
						Jpq_i[p][q][i] = Gpq_1[p][q][i];
					}
				}
			}
		}
	}

	// for eigen-heat-source
	Ms = G_0;
	if (nsolve != 4) {
		for (int p = 0; p < 3; p++) {
			Mp[p] = Gp_0[p];
			if (nsolve != 16) {
				for (int q = 0; q < 3; q++) {
					Mpq[p][q] = Gpq_0[p][q];
				}
			}
		}
	}

	// set them to nil
	destroy_0(nsolve);
	clear_all();
}

void DIE::D10_heat(int nsolve, double* x, double* x_p, double radius, complex<double>* M_i, complex<double>** Mp_i, complex<double>*** Mpq_i, \
	complex<double>** J_ij, complex<double>*** Jp_ij, complex<double>**** Jpq_ij)
{
	// initialization 
	inil_1(nsolve);
	a = radius;

	// for temperature
	Int_Gfunc_flux(nsolve, x, x_p);

	// for eigen-temeprature-gradient
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			J_ij[i][j] = G_2[i][j];
			if (nsolve != 4) {
				for (int p = 0; p < 3; p++) {
					Jp_ij[p][i][j] = Gp_2[p][i][j];
					if (nsolve != 16) {
						for (int q = 0; q < 3; q++) {
							Jpq_ij[p][q][i][j] = Gpq_2[p][q][i][j];
						}
					}
				}
			}

		}
	}

	// for eigen-heat-source
	for (int i = 0; i < 3; i++) {
		M_i[i] = G_1[i];
		if (nsolve != 4) {
			for (int p = 0; p < 3; p++) {
				Mp_i[p][i] = Gp_1[p][i];
				if (nsolve != 16) {
					for (int q = 0; q < 3; q++) {
						Mpq_i[p][q][i] = Gpq_1[p][q][i];
					}
				}
			}
		}
	}

	// set them to nil
	destroy_1(nsolve);
	clear_all();
}

void DIE::D20_heat(int nsolve, double* x, double* x_p, double radius, complex<double>** M_ij, complex<double>*** Mp_ij, complex<double>**** Mpq_ij, \
	complex<double>*** J_ijk, complex<double>**** Jp_ijk, complex<double>***** Jpq_ijk)
{
	// initialization 
	inil_2();
	a = radius;

	// for temperature
	Int_Gfunc_flux_grad1(nsolve, x, x_p);

	// for eigen-temeprature-gradient
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {

				J_ijk[i][j][k] = G_3[i][j][k];
				if (nsolve != 4) {
					for (int p = 0; p < 3; p++) {
						Jp_ijk[p][i][j][k] = Gp_3[p][i][j][k];
						if (nsolve != 16) {
							for (int q = 0; q < 3; q++) {
								Jpq_ijk[p][q][i][j][k] = Gpq_3[p][q][i][j][k];
							}
						}
					}
				}

			}
		}
	}

	// for eigen-heat-source
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			M_ij[i][j] = G_2[i][j];
			if (nsolve != 4) {
				for (int p = 0; p < 3; p++) {
					Mp_ij[p][i][j] = Gp_2[p][i][j];
					if (nsolve != 16) {
						for (int q = 0; q < 3; q++) {
							Mpq_ij[p][q][i][j] = Gpq_2[p][q][i][j];
						}
					}
				}
			}

		}
	}

	// set them to nil
	destroy_2();
	clear_all();
}

void DIE::D30_heat(int nsolve, double* x, double* x_p, double radius, complex<double>*** M_ijk, complex<double>**** Mp_ijk, complex<double>***** Mpq_ijk, \
	complex<double>**** J_ijkl, complex<double>***** Jp_ijkl, complex<double>****** Jpq_ijkl)
{
	// initialization 
	inil_3();
	a = radius;

	// for temperature
	Int_Gfunc_flux_grad2(nsolve, x, x_p);

	// for eigen-temeprature-gradient
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {

					J_ijkl[i][j][k][l] = G_4[i][j][k][l];
					if (nsolve != 4) {
						for (int p = 0; p < 3; p++) {
							Jp_ijkl[p][i][j][k][l] = Gp_4[p][i][j][k][l];
							if (nsolve != 16) {
								for (int q = 0; q < 3; q++) {
									Jpq_ijkl[p][q][i][j][k][l] = Gpq_4[p][q][i][j][k][l];
								}
							}
						}
					}

				}
			}
		}
	}

	// for eigen-heat-source
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {

				M_ijk[i][j][k] = G_3[i][j][k];
				if (nsolve != 4) {
					for (int p = 0; p < 3; p++) {
						Mp_ijk[p][i][j][k] = Gp_3[p][i][j][k];
						if (nsolve != 16) {
							for (int q = 0; q < 3; q++) {
								Mpq_ijk[p][q][i][j][k] = Gpq_3[p][q][i][j][k];
							}
						}
					}
				}

			}
		}
	}

	// set them to nil
	destroy_3();
	clear_all();
}