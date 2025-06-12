
#include "Input_parameters.h"
#include "eyemat.h"
#include <fstream>
#include <iostream>
#include <complex>

using namespace std;

void Reading_inputs::Start_read()
{
    read_soil();

	Start_input_particle();
	
	Start_input_post();
}

void Reading_inputs::read_soil()
{
    ifstream myfile; myfile.open("soil_properties.txt");

    int line_num = 0; string line; stringstream ss;
    while (!myfile.eof()) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0 && line[0] != '\t') {

            if (line_num == 0) {
                ss << line;
                ss >> K_0 >> Cp_0;
                line_num++;
                cout << K_0 << '\t' << Cp_0 << endl;
                // calculate properties:
                alpha_0 = K_0 / Cp_0;
                f_m = sqrt(omega / (2.0 * alpha_0));
                FF = f_m * (1.0 + 1.0i);

            }
        }
        else if (line[0] == '#') {
            cout << line << endl;
        }

    }
}

void Reading_inputs::Start_input_particle()
{
    ifstream myfile; myfile.open("particle.txt");

    int line_num = 0; string line; stringstream ss;
    while (!myfile.eof()) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0 && line[0] != '\t') {

            if (line_num == 0) {                
                num = stoi(line); line_num++;
                eigen_point.resize(num, 3); radius.resize(num); eigen_mat.resize(num, 3);
                cout << "Number of water tanks = " << " " << num << endl;

            }
            else if (line_num >= 1 && line_num < num + 1) {
                ss << line;
                ss >> eigen_point(line_num - 1, 0) >> eigen_point(line_num - 1, 1) >> eigen_point(line_num - 1, 2) \
                    >> radius(line_num - 1) >> eigen_mat(line_num - 1, 0) >> eigen_mat(line_num - 1, 1) >> eigen_mat(line_num - 1, 2);
                ss.clear();
                line_num++;
            }
            else if (line_num == num + 1) {
                if (line == "UNI") {
                    cout << "UNIFORM" << " " << "order of eigen-fields" << endl;
                    nsolve = 4;
                }
                else if (line == "LIN") {
                    cout << "LINEAR" << " " << "order of eigen-fields" << endl;
                    nsolve = 4 + 12;
                }
                else if (line == "QUA") {
                    cout << "QUADRATIC" << " " << "order of eigen-fields" << endl;
                    nsolve = 4 + 12 + 36;
                }
                else {
                    cout << "Warning: Wrong accuracy orders!" << endl;
                    system("pause");
                }
                particle_index();
            }

        }
        else if (line[0] == '#') {
            cout << line << endl;
        }

    }

    HMAT.resize(nsolve * num, nsolve * num); U.resize(nsolve * num);
    GMAT.resize(nsolve * num, nsolve * num); Heat_source.resize(nsolve * num);

    RHS.resize(nsolve * num);

    Start_input_heatsource();

    myfile.close();
}

void Reading_inputs::particle_index()
{
    /*
       The sequence follow the order:
       (1) T_i and Q
       (2) T_{ip} and Q_{p}
       (3) T_{ipq} and Q_{pq}
   */

    index_T_i = new int[3 * num]; index_T_ij = new int* [3 * num]; index_T_ijk = new int** [3 * num];

    index_Q = new int[num]; index_Q_i = new int[3 * num]; index_Q_ij = new int* [3 * num];

    for (int i = 0; i < 3 * num; i++) {
        index_T_ij[i] = new int[3];
        index_Q_ij[i] = new int[3];
    }

    for (int i = 0; i < 3 * num; i++)
    {
        index_T_ijk[i] = new int* [3];
        for (int j = 0; j < 3; j++)
        {
            index_T_ijk[i][j] = new int[3];
        }
    }

    int id = 0;

    for (int h = 0; h < num; h++) {

        // Uniform
        for (int i = 0; i < 3; i++) {
            index_T_i[3 * h + i] = id;
            id++;
        }

        index_Q[h] = id; id++;

        if (nsolve != 4) {
            // Linear
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    index_T_ij[3 * h + i][j] = id;
                    id++;
                }
            }

            for (int i = 0; i < 3; i++) {
                index_Q_i[3 * h + i] = id;
                id++;
            }

            // Quadratic
            if (nsolve != 16) {
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {

                            index_T_ijk[3 * h + i][j][k] = id;
                            id++;

                        }
                    }
                }

                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        index_Q_ij[3 * h + i][j] = id;
                        id++;
                    }
                }

            }


        }



    }

}

void Reading_inputs::Start_input_heatsource()
{
    for (int i = 0; i < num * nsolve; i++) {
        Heat_source(i) = 0.0;
    }

    for (int i = 0; i < num; i++) {
        Heat_source(index_Q[i]) = eigen_mat(i, 2);
    }
}

void Reading_inputs::Start_input_post()
{
    ifstream myfile; string line; stringstream ss;

    int line_num = 0;

    myfile.open("post_info.txt", ios::in);

    while (!myfile.eof()) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0 && line[0] != '\t') {

            if (line_num == 0) {
                nump = stoi(line);

                // initialization of temperature and flux in postprocess
                temp = new complex<double>[nump]; flux = new complex<double>*[nump];

                for (int i = 0; i < nump; i++) {
                    temp[i] = 0.0 + 0.0i;
                    flux[i] = new complex<double>[3];
                    for (int j = 0; j < 3; j++) {
                        flux[i][j] = 0.0 + 0.0i;
                    }
                }

                Points.resize(nump, 3);

                line_num++;
            }
            else if (line_num < nump + 1) {
                ss << line; ss >> Points(line_num - 1, 0) >> Points(line_num - 1, 1) >> Points(line_num - 1, 2);
                ss.clear();

                line_num++;
            }
        }
        else if (line[0] == '#') {
            // print the comments
            cout << line << endl;
        }

    }

    myfile.close();
}