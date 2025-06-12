
#include <iostream>
#include <complex>
#include "eigen-3.4.0/Eigen/Dense"

#include "Equivalent_conditions.h"
#include "Solve.h"
#include "Input_parameters.h"
#include "EIM_run.h"

int main()
{
	Reading_inputs READ; READ.Start_read(); 

	construct_EIM(READ);

	solution_EIM(READ);

	Post_process_thermal(READ);

	return 0;
}
