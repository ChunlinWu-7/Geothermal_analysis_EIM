
#include "Equivalent_conditions.h"

void EIM_integrals::inil_0(complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij, \
	complex<double>*& Mp, complex<double>**& Mpq, complex<double>*& J_i, complex<double>**& Jp_i, complex<double>***& Jpq_i)
{
	M_i = new complex<double>[3]; Mp_i = new complex<double>* [3]; Mpq_i = new complex<double>** [3];
	J_ij = new complex<double>* [3]; Jp_ij = new complex<double>** [3]; Jpq_ij = new complex<double>*** [3];

	Mp = new complex<double>[3];  Mpq = new complex<double>* [3]; J_i = new complex<double>[3];
	Jp_i = new complex<double>* [3]; Jpq_i = new complex<double>** [3];

	for (int i = 0; i < 3; i++) {
		Mp_i[i] = new complex<double>[3]; Mpq_i[i] = new complex<double>* [3];
		J_ij[i] = new complex<double>[3]; Jp_ij[i] = new complex<double>* [3]; Jpq_ij[i] = new complex<double>** [3];

		Mpq[i] = new complex<double>[3]; Jp_i[i] = new complex<double>[3]; Jpq_i[i] = new complex<double>* [3];

		for (int j = 0; j < 3; j++) {

			Mpq_i[i][j] = new complex<double>[3]; Jp_ij[i][j] = new complex<double>[3]; Jpq_ij[i][j] = new complex<double>* [3];
			Jpq_i[i][j] = new complex<double>[3];

			for (int k = 0; k < 3; k++) {
				Jpq_ij[i][j][k] = new complex<double>[3];
			}
		}

	}

}

void EIM_integrals::destroy_0(complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij, \
	complex<double>*& Mp, complex<double>**& Mpq, complex<double>*& J_i, complex<double>**& Jp_i, complex<double>***& Jpq_i)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				delete Jpq_ij[i][j][k];
			}
			delete Mpq_i[i][j];  delete Jp_ij[i][j];
			delete Jpq_ij[i][j]; delete Jpq_i[i][j];
		}

		delete Mp_i[i]; delete Mpq_i[i];
		delete J_ij[i]; delete Jp_ij[i]; delete Jpq_ij[i];

		delete Mpq[i]; delete Jp_i[i]; delete Jpq_i[i];
	}
	delete[] M_i; delete[] Mp_i; delete[] Mpq_i;
	delete[] J_ij; delete[] Jp_ij; delete[] Jpq_ij;

	delete[] Mp; delete[] Mpq; delete[] J_i;
	delete[] Jp_i; delete[] Jpq_i;
}

void EIM_integrals::inil_1(complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk, \
	complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij)
{
	M_ij = new complex<double>* [3]; Mp_ij = new complex<double>** [3]; Mpq_ij = new complex<double>*** [3];
	J_ijk = new complex<double>** [3]; Jp_ijk = new complex<double>*** [3]; Jpq_ijk = new complex<double>**** [3];

	M_i = new complex<double>[3]; Mp_i = new complex<double>* [3]; Mpq_i = new complex<double>** [3];
	J_ij = new complex<double>* [3]; Jp_ij = new complex<double>** [3]; Jpq_ij = new complex<double>*** [3];

	for (int i = 0; i < 3; i++) {
		M_ij[i] = new complex<double>[3]; Mp_ij[i] = new complex<double>* [3]; Mpq_ij[i] = new complex<double>** [3];
		J_ijk[i] = new complex<double>* [3]; Jp_ijk[i] = new complex<double>** [3]; Jpq_ijk[i] = new complex<double>*** [3];

		Mp_i[i] = new complex<double>[3]; Mpq_i[i] = new complex<double>* [3];
		J_ij[i] = new complex<double>[3]; Jp_ij[i] = new complex<double>* [3]; Jpq_ij[i] = new complex<double>** [3];

		for (int j = 0; j < 3; j++) {
			Mp_ij[i][j] = new complex<double>[3]; Mpq_ij[i][j] = new complex<double>* [3];
			J_ijk[i][j] = new complex<double>[3]; Jp_ijk[i][j] = new complex<double>* [3]; Jpq_ijk[i][j] = new complex<double>** [3];

			Mpq_i[i][j] = new complex<double>[3];
			Jp_ij[i][j] = new complex<double>[3]; Jpq_ij[i][j] = new complex<double>* [3];

			for (int k = 0; k < 3; k++) {
				Mpq_ij[i][j][k] = new complex<double>[3];
				Jp_ijk[i][j][k] = new complex<double>[3]; Jpq_ijk[i][j][k] = new complex<double>* [3];

				Jpq_ij[i][j][k] = new complex<double>[3];

				for (int l = 0; l < 3; l++) {
					Jpq_ijk[i][j][k][l] = new complex<double>[3];
				}
			}

		}
	}

}

void EIM_integrals::destroy_1(complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk, \
	complex<double>*& M_i, complex<double>**& Mp_i, complex<double>***& Mpq_i, complex<double>**& J_ij, complex<double>***& Jp_ij, complex<double>****& Jpq_ij)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					delete Jpq_ijk[i][j][k][l];
				}
				delete Mpq_ij[i][j][k];
				delete Jp_ijk[i][j][k]; delete Jpq_ijk[i][j][k];

				delete Jpq_ij[i][j][k];
			}
			delete Mp_ij[i][j]; delete Mpq_ij[i][j];
			delete J_ijk[i][j]; delete Jp_ijk[i][j]; delete Jpq_ijk[i][j];

			delete Mpq_i[i][j];
			delete Jp_ij[i][j]; Jpq_ij[i][j];
		}
		delete M_ij[i]; delete Mp_ij[i]; delete Mpq_ij[i];
		delete J_ijk[i]; delete Jp_ijk[i]; delete Jpq_ijk[i];

		delete Mp_i[i]; delete Mpq_i[i];
		delete J_ij[i]; delete Jp_ij[i]; delete Jpq_ij[i];
	}
	delete[] M_ij; delete[] Mp_ij; delete[] Mpq_ij;
	delete[] J_ijk; delete[] Jp_ijk; delete[] Jpq_ijk;

	delete[] M_i; delete[] Mp_i; delete[] Mpq_i;
	delete[] J_ij; delete[] Jp_ij; delete[] Jpq_ij;
}

void  EIM_integrals::inil_2(complex<double>***& M_ijk, complex<double>****& Mp_ijk, complex<double>*****& Mpq_ijk, complex<double>****& J_ijkl, complex<double>*****& Jp_ijkl, complex<double>******& Jpq_ijkl, \
	complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk)
{
	M_ijk = new complex<double>** [3]; Mp_ijk = new complex<double>*** [3]; Mpq_ijk = new complex<double>**** [3];
	J_ijkl = new complex<double>*** [3]; Jp_ijkl = new complex<double>**** [3]; Jpq_ijkl = new complex<double>***** [3];

	M_ij = new complex<double>* [3]; Mp_ij = new complex<double>** [3]; Mpq_ij = new complex<double>*** [3];
	J_ijk = new complex<double>** [3]; Jp_ijk = new complex<double>*** [3]; Jpq_ijk = new complex<double>**** [3];

	for (int i = 0; i < 3; i++) {

		M_ijk[i] = new complex<double>* [3]; Mp_ijk[i] = new complex<double>** [3]; Mpq_ijk[i] = new complex<double>*** [3];
		J_ijkl[i] = new complex<double>** [3]; Jp_ijkl[i] = new complex<double>*** [3]; Jpq_ijkl[i] = new complex<double>**** [3];

		M_ij[i] = new complex<double>[3]; Mp_ij[i] = new complex<double>* [3]; Mpq_ij[i] = new complex<double>** [3];
		J_ijk[i] = new complex<double>* [3]; Jp_ijk[i] = new complex<double>** [3]; Jpq_ijk[i] = new complex<double>*** [3];

		for (int j = 0; j < 3; j++) {
			M_ijk[i][j] = new complex<double>[3]; Mp_ijk[i][j] = new complex<double>* [3]; Mpq_ijk[i][j] = new complex<double>** [3];
			J_ijkl[i][j] = new complex<double>* [3]; Jp_ijkl[i][j] = new complex<double>** [3]; Jpq_ijkl[i][j] = new complex<double>*** [3];

			Mp_ij[i][j] = new complex<double>[3]; Mpq_ij[i][j] = new complex<double>* [3];
			J_ijk[i][j] = new complex<double>[3]; Jp_ijk[i][j] = new complex<double>* [3]; Jpq_ijk[i][j] = new complex<double>** [3];

			for (int k = 0; k < 3; k++) {

				Mp_ijk[i][j][k] = new complex<double>[3]; Mpq_ijk[i][j][k] = new complex<double>* [3];
				J_ijkl[i][j][k] = new complex<double>[3]; Jp_ijkl[i][j][k] = new complex<double>* [3]; Jpq_ijkl[i][j][k] = new complex<double>** [3];

				Mpq_ij[i][j][k] = new complex<double>[3];
				Jp_ijk[i][j][k] = new complex<double>[3]; Jpq_ijk[i][j][k] = new complex<double>* [3];

				for (int l = 0; l < 3; l++) {
					Mpq_ijk[i][j][k][l] = new complex<double>[3];
					Jp_ijkl[i][j][k][l] = new complex<double>[3]; Jpq_ijkl[i][j][k][l] = new complex<double>* [3];

					Jpq_ijk[i][j][k][l] = new complex<double>[3];

					for (int m = 0; m < 3; m++) {
						Jpq_ijkl[i][j][k][l][m] = new complex<double>[3];
					}

				}

			}


		}


	}


}

void  EIM_integrals::destroy_2(complex<double>***& M_ijk, complex<double>****& Mp_ijk, complex<double>*****& Mpq_ijk, complex<double>****& J_ijkl, complex<double>*****& Jp_ijkl, complex<double>******& Jpq_ijkl, \
	complex<double>**& M_ij, complex<double>***& Mp_ij, complex<double>****& Mpq_ij, complex<double>***& J_ijk, complex<double>****& Jp_ijk, complex<double>*****& Jpq_ijk)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						delete Jpq_ijkl[i][j][k][l][m];
					}
					delete Mpq_ijk[i][j][k][l];
					delete Jp_ijkl[i][j][k][l]; delete Jpq_ijkl[i][j][k][l];

					delete Jpq_ijk[i][j][k][l];
				}
				delete Mp_ijk[i][j][k]; delete Mpq_ijk[i][j][k];
				delete J_ijkl[i][j][k]; delete Jp_ijkl[i][j][k]; delete Jpq_ijkl[i][j][k];

				delete Mpq_ij[i][j][k];
				delete Jp_ijk[i][j][k]; delete Jpq_ijk[i][j][k];
			}
			delete M_ijk[i][j]; delete Mp_ijk[i][j]; delete Mpq_ijk[i][j];
			delete J_ijkl[i][j]; delete Jp_ijkl[i][j]; delete Jpq_ijkl[i][j];

			delete Mp_ij[i][j]; delete Mpq_ij[i][j];
			delete J_ijk[i][j]; delete Jp_ijk[i][j]; delete Jpq_ijk[i][j];
		}
		delete M_ijk[i]; delete Mp_ijk[i]; delete Mpq_ijk[i];
		delete J_ijkl[i]; delete Jp_ijkl[i]; delete Jpq_ijkl[i];

		delete M_ij[i]; delete Mp_ij[i]; delete Mpq_ij[i];
		delete J_ijk[i]; delete Jp_ijk[i]; delete Jpq_ijk[i];

	}

	delete[] M_ijk; delete[] Mp_ijk; delete[] Mpq_ijk;
	delete[] J_ijkl; delete[] Jp_ijkl; delete[] Jpq_ijkl;

	delete[] M_ij; delete[] Mp_ij; delete[] Mpq_ij;
	delete[] J_ijk; delete[] Jp_ijk; delete[] Jpq_ijk;
}