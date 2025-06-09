
#include "Input_parameters.h"
#include "eyemat.h"
#include <fstream>
#include <iostream>
#include <complex>

using namespace std;

void Reading_inputs::Start_read()
{
   // generate_particles();

	Start_input_particle();
	
	Start_input_post();
}

void Reading_inputs::temperature_work()
{
    // read temp[s] on Gauss integral points:

    double g_rho, g_theta, g_phi;

    complex<double> energy_saved = 0.0 + 0.0i;

    for (int s = 0; s < num; s++) {

        for (int i = 0; i < ngp; i++) {

            g_rho = 0.5 * (1.0 + gp[i]) * radius(s);
            for (int j = 0; j < ngp; j++) {
                g_theta = 0.5 * (1.0 + gp[j]) * 2.0 * pi;
                for (int k = 0; k < ngp; k++) {
                    g_phi = 0.5 * (1.0 + gp[k]) * 1.0 * pi;

                    energy_saved += (0.5 * radius(s)) * (0.5 * 2.0 * pi) * (0.5 * 1.0 * pi) * w[i] * w[j] * w[k] \
                        * g_rho * g_rho * sin(g_phi) * \
                        (
                            temp[s * ngp * ngp * ngp + i * ngp * ngp + j * ngp + k] * eigen_mat(s, 1) * 1339.59 / 3600.0
                            );
                }
            }
        }

    }

    // need to multiply with time difference:
    ofstream myfile; myfile.open("Energy_stored.txt");

    for (int j = 0; j < 364; j++) {
       // energy about the time at zero days: Jan 1st;
        myfile << (energy_saved * exp(-1.0i * omega * (double(j) - 364.0 / 2.0) * 364.0 / 364.0)
            - energy_saved * exp(-1.0i * omega * (0.0 - 364.0 / 2.0) * 364.0 / 364.0)).real();

        myfile << " ";\
        myfile << (20.0 * 4.0 * pi / 3.0 * pow(0.5, 3.0) * 1.0i / omega * (\
            exp(-1.0i * omega * (double(j) - 364.0 / 2.0) * 364.0 / 364.0)\
            - exp(-1.0i * omega * (0.0 - 364.0 / 2.0) * 364.0 / 364.0)\
            )).real();

        myfile << endl;

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
                    cout << "UNIFORM" << '\t' << "order of eigenstain" << endl;
                    nsolve = 4;
                }
                else if (line == "LIN") {
                    cout << "LINEAR" << '\t' << "order of eigenstain" << endl;
                    nsolve = 4 + 12;
                }
                else if (line == "QUA") {
                    cout << "QUADRATIC" << '\t' << "order of eigenstain" << endl;
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

    //cout << "input Qv" << endl;\
    double div;\
    cin >> div;

    for (int i = 0; i < num; i++) {
       // Heat_source(index_Q[i]) = eigen_mat(i, 2) * pow(0.5, 3.0) / pow(radius(0), 3.0) * 12.5 / 12.5 * (10.0 / double(num));
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
                nump = 801;

                // initialization of temperature and flux in postprocess
                temp = new complex<double> [nump]; flux = new complex<double>* [nump];

                for (int i = 0; i < nump; i++) {
                    temp[i] = 0.0 + 0.0i;
                    flux[i] = new complex<double> [3];
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

   // along the axis symmetric:
    for (int i = 0; i < nump; i++) {\
        //Points(i, 2) = eigen_point(0, 2); Points(i, 1) = eigen_point(0, 1); \
        //Points(i, 2) = -25.0 / (nump - 1) * i;
        //Points(i, 0) = eigen_point(0,0) - 15.0 * radius(0) + 30.0 * radius(0) / 5000.0 * i; 
        //Points(i, 0) = 0.0 - 25.0 * radius(0) + 50.0 * radius(0) / 5000.0 * i; \

        Points(i, 0) = eigen_point(0, 0); Points(i, 1) = eigen_point(0, 1);

        Points(i, 2) = 0.0 - 0.1 * 0.25 * i;
    }
    //cout << Points << endl;
   // center of the tank:
    //nump = 1; \
        Points(0, 0) = 0.0; eigen_point(0, 0); \
    Points(0, 2) = eigen_point(0, 2);

    // contour plot:
    /*
    nump = 201 * 201;
    Points.resize(nump, 3);
    temp = new complex<double>[nump];
    for (int i = 0; i < nump; i++) {
        temp[i] = 0.0;
    }

    for (int i = 0; i < 201; i++) {
        for (int j = 0; j < 201; j++) {
            Points(i * 201 + j, 0) = 0.0;
            Points(i * 201 + j, 1) = eigen_point(0, 1) + 2.0 * radius(0) - 4.0 * radius(0) / (200.0) * i;
            Points(i * 201 + j, 2) = eigen_point(0, 2) + 2.0 * radius(0) - 4.0 * radius(0) / (200.0) * j;
        }
    }
    */
    
    // volume integral:
    /*
    nump = num * ngp * ngp * ngp;
    Points.resize(nump, 3);
    temp = new complex<double>[nump];
    for (int i = 0; i < nump; i++) {
        temp[i] = 0.0;
    }

    double g_theta, g_phi, g_rho;

    for (int s = 0; s < num; s++) {

        for (int i = 0; i < ngp; i++) {

            g_rho = 0.5 * (1.0 + gp[i]) * radius(s);

            for (int j = 0; j < ngp; j++) {

                g_theta = 0.5 * (1.0 + gp[j]) * 2.0 * pi;

                for (int k = 0; k < ngp; k++) {
                    g_phi = 0.5 * (1.0 + gp[k]) * pi;

                    Points(s * ngp * ngp * ngp + i * ngp * ngp + j * ngp + k, 0) = eigen_point(s, 0) + g_rho * cos(g_theta) * sin(g_phi);
                    Points(s * ngp * ngp * ngp + i * ngp * ngp + j * ngp + k, 1) = eigen_point(s, 1) + g_rho * sin(g_theta) * sin(g_phi);
                    Points(s * ngp * ngp * ngp + i * ngp * ngp + j * ngp + k, 2) = eigen_point(s, 2) + g_rho * cos(g_phi);

                }
            }
        }

    }
    */


    /*
    double rrr = radius(0);

    for (int i = 0; i < 200; i++) {
        
        Points(i, 0) = 0.0; Points(i, 1) = 0.0;

        Points(i, 2) = (eigen_point(0, 2) + 4.0 * rrr) / (200.0) * (i + 1.0);

    }

    for (int i = 200; i < 800 + 200; i++) {

        Points(i, 0) = 0.0; Points(i, 1) = 0.0;

        Points(i, 2) = (eigen_point(0, 2) + 4.0 * rrr) - 8.0 * rrr / 800.0 * (i - 200);

    }

    for (int i = 800 + 200; i < nump; i++) {

        Points(i, 0) = 0.0; Points(i, 1) = 0.0;

        Points(i, 2) = (eigen_point(0, 2) - 4.0 * rrr) - (eigen_point(0, 2) - 4.0 * rrr + 20.0) / (double(nump) - 1000.0) * (i - 1000.0);
    }
    
    */
   // for (int i = 0; i < nump; i++) {\
        Points(i, 2) = eigen_point(0, 2) - 4.0 * radius(0) + 8.0 * radius(0) / (nump - 1) * i;\
    }



    ofstream myfile_s; myfile_s.open("post_info_origin.txt");

    for (int i = 0; i < nump; i++) {
        myfile_s << (Points(i, 0)) / radius(0) << endl;
    }

    myfile_s.close();
}

void Reading_inputs::Start_input_post_BC()
{
    ifstream myfile; string line; stringstream ss;

    int line_num = 0;

    myfile.open("post_nodes.txt", ios::in);

    while (!myfile.eof()) {
        getline(myfile, line);

        if (line[0] != '#' && line.length() != 0 && line[0] != '\t') {

            if (line_num == 0) {
                NN = stoi(line);

                // initialization of temperature and flux in postprocess
                temp_nodes = new complex<double>[NN]; 

                for (int i = 0; i < NN; i++) {
                    temp_nodes[i] = 0.0 + 0.0i;
                }

                NODES.resize(NN, 3);

                line_num++;
            }
            else if (line_num < NN + 1) {
                ss << line; ss >> NODES(line_num - 1, 0) >> NODES(line_num - 1, 1) >> NODES(line_num - 1, 2);
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

void Reading_inputs::generate_particles()
{
    ofstream myfile; myfile.open("particle.txt");

    double volume_total = 4.0 * pi / 3.0 * pow(1.75, 3.0) * 2.0;

    num = 7; double a = pow(volume_total / double(num) / (4.0 * pi) * 3.0, 1.0 / 3.0);

    double CD; 

    cout << "input CD" << endl;

    cin >> CD;

    myfile << num << endl;

    for (int i = 0; i < num; i++) {

        // annual
        myfile << -double(num - 1) * 0.5 * double(CD) * a + i * double(CD) * a << " " << \
            0 << " " << -15.0 << " " << a << " " \
            << 15.8308 << " " << 3.03533 << " " << 33.5923 << endl;


        // daily
        //myfile << -double(num - 1) * 0.5 * 2.0 * a + i * 2.0 * a << " " << \
            0 << " " << -10.0 * a << " " << a << " " \
            << 0.00592531 * 1000.0 << " " << 5.39823 * 1000.0 << " " << 7.79681 * 1000.0 << endl;
    }

    myfile << "QUA" << endl;

    myfile.close();
}