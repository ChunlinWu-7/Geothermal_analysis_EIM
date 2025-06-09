
#include <complex>
#include "eyemat.h"
#include "Domain_integrals.h"

using namespace std;

/*
	Inil_0 and destroy_0 serves for temperature and equivalent temperature conditions!
*/

void DIE::inil_0(int nsolve)
{
	Gp_0 = new complex<double>[3]; Gpq_0 = new complex<double>* [3];

	G_1 = new complex<double>[3]; Gp_1 = new complex<double> *[3];
	Gpq_1 = new complex<double> **[3];

	for (int i = 0; i < 3; i++) {
		Gpq_0[i] = new complex<double>[3];
		Gp_1[i] = new complex<double>[3];

		Gpq_1[i] = new complex<double> *[3];

		for (int j = 0; j < 3; j++) {
			Gpq_1[i][j] = new complex<double>[3];
		}
	}
}

void DIE::destroy_0(int nsolve)
{
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {
			delete Gpq_1[i][j];
		}

		delete Gpq_0[i];

		delete Gp_1[i]; delete Gpq_1[i];
	}

	delete[] Gp_0; delete[] Gpq_0; 
	delete[] G_1; delete[] Gp_1; delete[] Gpq_1;
}

/*
	Inil_0 and destroy_0 serves for heat flux
	(i) equivalent heat source gradient;
	(ii) or equivalent heat flux condition
*/

void DIE::inil_1(int nsolve)
{
	G_1 = new complex<double>[3]; Gp_1 = new complex<double> *[3]; Gpq_1 = new complex<double> **[3];

	G_2 = new complex<double>*[3]; Gp_2 = new complex<double> **[3]; Gpq_2 = new complex<double> ***[3];
	
	for (int i = 0; i < 3; i++) {

		Gp_1[i] = new complex<double>[3]; Gpq_1[i] = new complex<double> *[3];

		G_2[i] = new complex<double>[3]; Gp_2[i] = new complex<double> *[3]; Gpq_2[i] = new complex<double> **[3];

		for (int j = 0; j < 3; j++) {
			Gpq_1[i][j] = new complex<double>[3];

			Gp_2[i][j] = new complex<double>[3]; Gpq_2[i][j] = new complex<double> *[3];

			for (int k = 0; k < 3; k++) {
				Gpq_2[i][j][k] = new complex<double>[3];
			}
		}
	}

}

void DIE::destroy_1(int nsolve)
{
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {

				delete Gpq_2[i][j][k];
			}

			delete Gpq_1[i][j];

			delete Gp_2[i][j]; delete Gpq_2[i][j];
		}
		delete Gp_1[i]; delete Gpq_1[i];

		delete G_2[i]; delete Gp_2[i]; delete Gpq_2[i];
	}

	delete[]G_1; delete[]Gp_1; delete[] Gpq_1;

	delete[]G_2; delete[]Gp_2; delete[]Gpq_2;
}

/*
	Inil_2 and destroy_2 serves
	(i) equivalent heat source gradient 2;
	(ii) equivalent heat flux gradient;
*/

void DIE::inil_2()
{
	G_2 = new complex<double>*[3]; Gp_2 = new complex<double> **[3]; Gpq_2 = new complex<double> ***[3];

	G_3 = new complex<double> **[3]; Gp_3 = new complex<double> ***[3]; Gpq_3 = new complex<double> ****[3];

	for (int i = 0; i < 3; i++) {

		G_2[i] = new complex<double>[3]; Gp_2[i] = new complex<double> *[3]; Gpq_2[i] = new complex<double> **[3];

		G_3[i] = new complex<double> *[3]; Gp_3[i] = new complex<double> **[3]; Gpq_3[i] = new complex<double> ***[3];

		for (int j = 0; j < 3; j++) {

			Gp_2[i][j] = new complex<double>[3]; Gpq_2[i][j] = new complex<double> *[3];

			G_3[i][j] = new complex<double>[3]; Gp_3[i][j] = new complex<double> *[3]; Gpq_3[i][j] = new complex<double>** [3];

			for (int k = 0; k < 3; k++) {
				Gpq_2[i][j][k] = new complex<double>[3];

				Gp_3[i][j][k] = new complex<double>[3]; Gpq_3[i][j][k] = new complex<double> *[3];

				for (int l = 0; l < 3; l++) {
					Gpq_3[i][j][k][l] = new complex<double>[3];
				}
			}
		}
	}
}

void DIE::destroy_2()
{
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {
				
				for (int l = 0; l < 3; l++) {
					delete Gpq_3[i][j][k][l];
				}

				delete Gpq_2[i][j][k];

				delete Gp_3[i][j][k]; delete Gpq_3[i][j][k];
			}

			delete Gp_2[i][j]; delete Gpq_2[i][j];

			delete G_3[i][j]; delete Gp_3[i][j]; delete Gpq_3[i][j];
		}

		delete G_2[i]; delete Gp_2[i]; delete Gpq_2[i];

		delete G_3[i]; delete Gp_3[i]; delete Gpq_3[i];
	}

	delete[] G_2; delete[] Gp_2; delete[] Gpq_2;

	delete[] G_3; delete[] Gp_3; delete[] Gpq_3;
}

/*
	Inil_3 and destroy_3 serves only for:
	(i) equivalent heat flux gradient 2;
*/

void DIE::inil_3()
{
	G_3 = new complex<double> **[3]; Gp_3 = new complex<double> ***[3]; Gpq_3 = new complex<double> ****[3];

	G_4 = new complex<double> ***[3]; Gp_4 = new complex<double> ****[3]; Gpq_4 = new complex<double> *****[3];

	for (int i = 0; i < 3; i++) {

		G_3[i] = new complex<double> *[3]; Gp_3[i] = new complex<double> **[3]; Gpq_3[i] = new complex<double> ***[3];

		G_4[i] = new complex<double> **[3]; Gp_4[i] = new complex<double> ***[3]; Gpq_4[i] = new complex<double> ****[3];

		for (int j = 0; j < 3; j++) {

			G_3[i][j] = new complex<double>[3]; Gp_3[i][j] = new complex<double> *[3]; Gpq_3[i][j] = new complex<double>**[3];

			G_4[i][j] = new complex<double> *[3]; Gp_4[i][j] = new complex<double> **[3]; Gpq_4[i][j] = new complex<double> ***[3];

			for (int k = 0; k < 3; k++) {

				Gp_3[i][j][k] = new complex<double>[3]; Gpq_3[i][j][k] = new complex<double> *[3];

				G_4[i][j][k] = new complex<double> [3]; Gp_4[i][j][k] = new complex<double> *[3]; Gpq_4[i][j][k] = new complex<double> **[3];

				for (int l = 0; l < 3; l++) {
					Gpq_3[i][j][k][l] = new complex<double>[3];

					Gp_4[i][j][k][l] = new complex<double>[3]; Gpq_4[i][j][k][l] = new complex<double> *[3];

					for (int m = 0; m < 3; m++) {
						Gpq_4[i][j][k][l][m] = new complex<double>[3]; 
					}

				}
			}
		}
	}
}

void DIE::destroy_3()
{
	for (int i = 0; i < 3; i++) {

		for (int j = 0; j < 3; j++) {

			for (int k = 0; k < 3; k++) {

				for (int l = 0; l < 3; l++) {
					
					for (int m = 0; m < 3; m++) {
						delete Gpq_4[i][j][k][l][m];
					}
					delete Gpq_3[i][j][k][l];

					delete Gp_4[i][j][k][l]; delete Gpq_4[i][j][k][l];
				}

				delete Gp_3[i][j][k]; delete Gpq_3[i][j][k];

				delete G_4[i][j][k]; delete Gp_4[i][j][k]; delete Gpq_4[i][j][k];

			}

			delete G_3[i][j]; delete Gp_3[i][j]; delete Gpq_3[i][j];

			delete G_4[i][j]; delete Gp_4[i][j]; delete Gpq_4[i][j];

		}

		delete G_3[i]; delete Gp_3[i]; delete Gpq_3[i];

		delete G_4[i]; delete Gp_4[i]; delete Gpq_4[i];

	}

	delete[] G_3; delete[] Gp_3; delete[] Gpq_3;

	delete[] G_4; delete[] Gp_4; delete[] Gpq_4;

}
