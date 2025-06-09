
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

// temperature:
void DIE::Int_Gfunc_temp(int nsolve, double* x, double* x_p)
{
	/*
		Must assign values before operating!
	*/
	beta = f_m * (1.0 + 1.0i);
	initial_distance(nsolve, x, x_p);
	Store_derivs(nsolve, x, x_p);

	// assign values to Green's function with integrals: G_{,i'} for ETG
	for (int i = 0; i < 3; i++) {

		G_1[i] = -0.25 / (pi) * (
			Helm_ori_1(i) + full_flag * M[i] * Helm_1(i)
			);

		if (nsolve != 4) {

			for (int p = 0; p < 3; p++) {

				Gp_1[p][i] = -0.25 / (pi) * (
					Helm_ori_1_1(p, i) + full_flag * M[i] * M[p] * Helm_1_1(p, i)
					);

				if (nsolve != 16) {

					for (int q = 0; q < 3; q++) {

						Gpq_1[p][q][i] = -0.25 / (pi) * (
							Helm_ori_2_1(p, q, i) + full_flag * M[i] * M[p] * M[q] * Helm_2_1(p, q, i)
							);

					}

				}

			}
		}

	}

	// assign values to Green's function with integrals: G for heat source

	G_0 = 0.25 / (pi * K_0) * (
		Helm_ori_0() + full_flag * Helm_0()
		);

	if (nsolve != 4) {

		for (int p = 0; p < 3; p++) {

			Gp_0[p] = 0.25 / (pi * K_0) * (
				Helm_ori_1_0(p) + full_flag * M[p] * Helm_1_0(p)
				);

			if (nsolve != 16) {

				for (int q = 0; q < 3; q++) {

					Gpq_0[p][q] = 0.25 / (pi * K_0) * (
						Helm_ori_2_0(p, q) + full_flag * M[p] * M[q] * Helm_2_0(p, q)
						);

				}

			}

		}
	}
}

void DIE::Int_Gfunc_flux(int nsolve, double* x, double* x_p)
{
	/*
		Must assign values before operating!
	*/
	beta = f_m * (1.0 + 1.0i);
	initial_distance(nsolve, x, x_p);
	Store_derivs(nsolve, x, x_p);

	// assign values to Green's function with integrals: G_{,i'} for ETG
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			G_2[i][j] = -0.25 / (pi) * (
				Helm_ori_2(i, j) + full_flag * M[i] * Helm_2(i, j)
				);

			if (nsolve != 4) {

				for (int p = 0; p < 3; p++) {

					Gp_2[p][i][j] = -0.25 / (pi) * (
						Helm_ori_1_2(p, i, j) + full_flag * M[i] * M[p] * Helm_1_2(p, i, j)
						);

					if (nsolve != 16) {

						for (int q = 0; q < 3; q++) {

							Gpq_2[p][q][i][j] = -0.25 / (pi) * (
								Helm_ori_2_2(p, q, i, j) + full_flag * M[i] * M[p] * M[q] * Helm_2_2(p, q, i, j)
								);

						}

					}

				}
			}

		}
	}

	// assign values to Green's function with integrals: G for heat source

	for (int i = 0; i < 3; i++) {
		G_1[i] = 0.25 / (pi * K_0) * (
			Helm_ori_1(i) + full_flag * Helm_1(i)
			);

		if (nsolve != 4) {

			for (int p = 0; p < 3; p++) {

				Gp_1[p][i] = 0.25 / (pi * K_0) * (
					Helm_ori_1_1(p, i) + full_flag * M[p] * Helm_1_1(p, i)
					);

				if (nsolve != 16) {

					for (int q = 0; q < 3; q++) {

						Gpq_1[p][q][i] = 0.25 / (pi * K_0) * (
							Helm_ori_2_1(p, q, i) + full_flag * M[p] * M[q] * Helm_2_1(p, q, i)
							);

					}

				}

			}
		}
	}
}

void DIE::Int_Gfunc_flux_grad1(int nsolve, double* x, double* x_p)
{
	/*
		Must assign values before operating!
	*/
	beta = f_m * (1.0 + 1.0i);
	initial_distance(nsolve, x, x_p);
	Store_derivs(nsolve, x, x_p);

	// assign values to Green's function with integrals: G_{,i'} for ETG
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {

				G_3[i][j][k] = -0.25 / (pi) * (
					Helm_ori_3(i, j, k) + full_flag * M[i] * Helm_3(i, j, k)
					);

				if (nsolve != 4) {

					for (int p = 0; p < 3; p++) {

						Gp_3[p][i][j][k] = -0.25 / (pi) * (
							Helm_ori_1_3(p, i, j, k) + full_flag * M[i] * M[p] * Helm_1_3(p, i, j, k)
							);

						if (nsolve != 16) {

							for (int q = 0; q < 3; q++) {

								Gpq_3[p][q][i][j][k] = -0.25 / (pi) * (
									Helm_ori_2_3(p, q, i, j, k) + full_flag * M[i] * M[p] * M[q] * Helm_2_3(p, q, i, j, k)
									);

							}

						}

					}
				}

			}
		}
	}

	// assign values to Green's function with integrals: G for heat source

	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			G_2[i][j] = 0.25 / (pi * K_0) * (
				Helm_ori_2(i, j) + full_flag * Helm_2(i, j)
				);

			if (nsolve != 4) {

				for (int p = 0; p < 3; p++) {

					Gp_2[p][i][j] = 0.25 / (pi * K_0) * (
						Helm_ori_1_2(p, i, j) + full_flag * M[p] * Helm_1_2(p, i, j)
						);

					if (nsolve != 16) {

						for (int q = 0; q < 3; q++) {

							Gpq_2[p][q][i][j] = 0.25 / (pi * K_0) * (
								Helm_ori_2_2(p, q, i, j) + full_flag * M[p] * M[q] * Helm_2_2(p, q, i, j)
								);

						}

					}

				}
			}
		}
	}
}

void DIE::Int_Gfunc_flux_grad2(int nsolve, double* x, double* x_p)
{
	/*
		Must assign values before operating!
	*/
	beta = f_m * (1.0 + 1.0i);
	initial_distance(nsolve, x, x_p);
	Store_derivs(nsolve, x, x_p);

	// assign values to Green's function with integrals: G_{,i'} for ETG
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {

				for (int m = 0; m < 3; m++) {

					G_4[i][j][k][m] = -0.25 / (pi) * (
						Helm_ori_4(i, j, k, m) + full_flag * M[i] * Helm_4(i, j, k, m)
						);

					if (nsolve != 4) {

						for (int p = 0; p < 3; p++) {

							Gp_4[p][i][j][k][m] = -0.25 / (pi) * (
								Helm_ori_1_4(p, i, j, k, m) + full_flag * M[i] * M[p] * Helm_1_4(p, i, j, k, m)
								);

							if (nsolve != 16) {

								for (int q = 0; q < 3; q++) {

									Gpq_4[p][q][i][j][k][m] = -0.25 / (pi) * (
										Helm_ori_2_4(p, q, i, j, k, m) + full_flag * M[i] * M[p] * M[q] * Helm_2_4(p, q, i, j, k, m)
										);

								}

							}

						}
					}

				}
			}
		}
	}

	// assign values to Green's function with integrals: G for heat source

	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {

				G_3[i][j][k] = 0.25 / (pi * K_0) * (
					Helm_ori_3(i, j, k) + full_flag * Helm_3(i, j, k)
					);

				if (nsolve != 4) {

					for (int p = 0; p < 3; p++) {

						Gp_3[p][i][j][k] = 0.25 / (pi * K_0) * (
							Helm_ori_1_3(p, i, j, k) + full_flag * M[p] * Helm_1_3(p, i, j, k)
							);

						if (nsolve != 16) {

							for (int q = 0; q < 3; q++) {

								Gpq_3[p][q][i][j][k] = 0.25 / (pi * K_0) * (
									Helm_ori_2_3(p, q, i, j, k) + full_flag * M[p] * M[q] * Helm_2_3(p, q, i, j, k)
									);

							}

						}

					}
				}
			}
		}
	}
}

void DIE::clear_all()
{
	delete[] difx; delete[] difx_ori;
}