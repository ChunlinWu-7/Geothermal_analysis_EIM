#pragma once

#include "eigen-3.4.0/Eigen/Dense"
#include <complex>

using namespace std;
using namespace Eigen;

class Reading_inputs
{
public:
	// boundary information
	int NN, nump;

	// Matrix collection of BEM info
	MatrixXcd HMAT; MatrixXd Points, NODES;
	VectorXcd U;

	// initial heat source:
	MatrixXcd GMAT; VectorXd Heat_source;

	// temperature gradient:
	VectorXcd RHS;

	// output definition, temperature and heat flux
	complex<double> * temp;  complex<double>** flux;
	complex<double> * temp_nodes; complex<double>** flux_nodes;

	// particle information
	int num, nsolve;

	// point, radius and properties
	MatrixXd eigen_point, eigen_mat; VectorXd radius;

	int* index_T_i, ** index_T_ij, *** index_T_ijk;

	int* index_Q, * index_Q_i, ** index_Q_ij;
	void Start_read();
	
	// postprocess of temperature
	void temperature_work();

private:
	// particle information
	void Start_input_particle();
	void particle_index();
	void Start_input_post();
	void Start_input_post_BC();
	void Start_input_heatsource();

	// generate particles
	void generate_particles();

};