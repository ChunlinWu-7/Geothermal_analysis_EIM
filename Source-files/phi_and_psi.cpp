
#include "Domain_integrals.h"
#include <cmath>
#include "eyemat.h"

void DIE::initial_distance(int nsolve, double* x, double* x_p)
{
    difx = new double[3]; difx_ori = new double[3];

    dist = 0.0; dist_ori = 0.0;

    for (int i = 0; i < 3; i++) {

        dist += (x[i] - x_p[i] * M[i]) * (x[i] - x_p[i] * M[i]);
        dist_ori += (x[i] - x_p[i]) * (x[i] - x_p[i]);
       
        difx[i] = x[i] - x_p[i] * M[i];
        difx_ori[i] = (x[i] - x_p[i]);
      
    }
    dist = sqrt(dist);
    dist_ori = sqrt(dist_ori);

    if (dist_ori <= 1E-6) {
        dist_flag = 0;
    }
    else {
        dist_flag = 1;
    }

    dist_1 = 1.0 / dist; dist_3 = pow(dist_1, 3.0); dist_5 = pow(dist_1, 5.0); 
    dist_7 = pow(dist_1, 7.0); dist_9 = pow(dist_1, 9.0); dist_11 = pow(dist_1, 11.0);

    dist_ori_1 = 1.0 / dist_ori; dist_ori_3 = pow(dist_ori_1, 3.0); dist_ori_5 = pow(dist_ori_1, 5.0); 
    dist_ori_7 = pow(dist_ori_1, 7.0); dist_ori_9 = pow(dist_ori_1, 9.0); dist_ori_11 = pow(dist_ori_1, 11.0);
}

/*
    phi(): imaged
*/
double DIE::phi()
{
    return dist_1;
}

double DIE::phi_1(int i)
{
    return -difx[i] * dist_3;
}

double DIE::phi_2(int i, int j)
{
    return -d[i][j] * dist_3 + 3.0 * dist_5 * difx[i] * difx[j];
}

double DIE::phi_3(int i, int j, int k)
{
    return 3.0 * dist_5 * (
        d[i][j] * difx[k] + d[i][k] * difx[j] + d[j][k] * difx[i]
        )
        - 15.0 * dist_7 * difx[i] * difx[j] * difx[k];
}

double DIE::phi_4(int i, int j, int k, int l)
{
    return 3.0 * dist_5 * (
        d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[j][k] * d[i][l]
        )
        - 15.0 * dist_7 * difx[l] * (
            d[i][j] * difx[k] + d[i][k] * difx[j] + d[j][k] * difx[i]
            )
        - 15.0 * dist_7 * (
            d[i][l] * difx[j] * difx[k]
            + d[j][l] * difx[i] * difx[k]
            + d[k][l] * difx[i] * difx[j]
            )
        + 105.0 * dist_9 * difx[i] * difx[j] * difx[k] * difx[l];
}

double DIE::phi_5(int i, int j, int k, int l, int m)
{

    return -15.0 * dist_7 * (
        difx[m] * (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[j][k] * d[i][l])
        + difx[i] * (d[j][k] * d[l][m] + d[j][l] * d[k][m] + d[k][l] * d[j][m])
        + difx[j] * (d[i][k] * d[l][m] + d[i][l] * d[k][m] + d[k][l] * d[i][m])
        + difx[k] * (d[i][j] * d[l][m] + d[i][l] * d[j][m] + d[j][l] * d[i][m])
        + difx[l] * (d[i][k] * d[j][m] + d[j][k] * d[i][m] + d[i][j] * d[k][m])
        )

        + 105.0 * dist_9 * (
            difx[i] * difx[j] * (d[k][l] * difx[m] + d[l][m] * difx[k] + d[k][m] * difx[l])
            + difx[i] * difx[k] * (d[j][l] * difx[m] + d[j][m] * difx[l])
            + difx[j] * difx[k] * (d[i][m] * difx[l] + d[i][l] * difx[m])
            + difx[l] * difx[m] * (d[i][j] * difx[k] + d[i][k] * difx[j] + d[j][k] * difx[i])
            )

        - 945.0 * dist_11 * difx[i] * difx[j] * difx[k] * difx[l] * difx[m];

}

/*
    Psi is totally dependent on phi(): imaged!
*/
double DIE::psi()
{
    return dist;
}

double DIE::psi_1(int i)
{
    return difx[i] * dist_1;
}

double DIE::psi_2(int i, int j)
{
    return d[i][j] * phi() + difx[i] * phi_1(j);
}

double DIE::psi_3(int i, int j, int k)
{
    return d[i][j] * phi_1(k) + d[i][k] * phi_1(j) + difx[i] * phi_2(j, k);
}

double DIE::psi_4(int i, int j, int k, int l)
{
    return d[i][j] * phi_2(k, l) + d[i][k] * phi_2(j, l) + d[i][l] * phi_2(j, k) + difx[i] * phi_3(j, k, l);
}

double DIE::psi_5(int i, int j, int k, int l, int m)
{
    return difx[i] * phi_4(j, k, l, m) + d[i][j] * phi_3(k, l, m)
        + d[i][k] * phi_3(j, l, m) + d[i][l] * phi_3(j, k, m) + d[i][m] * phi_3(j, k, l);
}

double DIE::psi_6(int i, int j, int k, int l, int m, int n)
{
    return difx[i] * phi_5(j, k, l, m, n) + d[i][j] * phi_4(k, l, m, n)
        + d[i][k] * phi_4(j, l, m, n) + d[i][l] * phi_4(j, k, m, n)
        + d[i][m] * phi_4(j, k, l, n) + d[i][n] * phi_4(j, k, l, m);
}

/*
    phi: original
*/

double DIE::phi_ori()
{
    return dist_ori_1;
}

double DIE::phi_ori_1(int i)
{
    return -difx_ori[i] * dist_ori_3;
}

double DIE::phi_ori_2(int i, int j)
{
    return -d[i][j] * dist_ori_3 + 3.0 * dist_ori_5 * difx_ori[i] * difx_ori[j];
}

double DIE::phi_ori_3(int i, int j, int k)
{
    return 3.0 * dist_ori_5 * (
        d[i][j] * difx_ori[k] + d[i][k] * difx_ori[j] + d[j][k] * difx_ori[i]
        )
        - 15.0 * dist_ori_7 * difx_ori[i] * difx_ori[j] * difx_ori[k];
}

// original form fourth derivatives of phi
double DIE::phi_ori_4(int i, int j, int k, int l)
{
    return 3.0 * dist_ori_5 * (
        d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[j][k] * d[i][l]
        )
        - 15.0 * dist_ori_7 * difx_ori[l] * (
            d[i][j] * difx_ori[k] + d[i][k] * difx_ori[j] + d[j][k] * difx_ori[i]
            )
        - 15.0 * dist_ori_7 * (
            d[i][l] * difx_ori[j] * difx_ori[k]
            + d[j][l] * difx_ori[i] * difx_ori[k]
            + d[k][l] * difx_ori[i] * difx_ori[j]
            )
        + 105.0 * dist_ori_9 * difx_ori[i] * difx_ori[j] * difx_ori[k] * difx_ori[l];
}

// original form fifth derivatives of phi
double DIE::phi_ori_5(int i, int j, int k, int l, int m)
{
    return -15.0 * dist_ori_7 * (
        difx_ori[m] * (d[i][j] * d[k][l] + d[i][k] * d[j][l] + d[j][k] * d[i][l])
        + difx_ori[i] * (d[j][k] * d[l][m] + d[j][l] * d[k][m] + d[k][l] * d[j][m])
        + difx_ori[j] * (d[i][k] * d[l][m] + d[i][l] * d[k][m] + d[k][l] * d[i][m])
        + difx_ori[k] * (d[i][j] * d[l][m] + d[i][l] * d[j][m] + d[j][l] * d[i][m])
        + difx_ori[l] * (d[i][k] * d[j][m] + d[j][k] * d[i][m] + d[i][j] * d[k][m])
        )

        + 105.0 * dist_ori_9 * (
            difx_ori[i] * difx_ori[j] * (d[k][l] * difx_ori[m] + d[l][m] * difx_ori[k] + d[k][m] * difx_ori[l])
            + difx_ori[i] * difx_ori[k] * (d[j][l] * difx_ori[m] + d[j][m] * difx_ori[l])
            + difx_ori[j] * difx_ori[k] * (d[i][m] * difx_ori[l] + d[i][l] * difx_ori[m])
            + difx_ori[l] * difx_ori[m] * (d[i][j] * difx_ori[k] + d[i][k] * difx_ori[j] + d[j][k] * difx_ori[i])
            )

        - 945.0 * dist_ori_11 * difx_ori[i] * difx_ori[j] * difx_ori[k] * difx_ori[l] * difx_ori[m];


}

/*
    Psi is totally dependent on phi(): original!
*/
double DIE::psi_ori()
{
    return dist_ori;
}

double DIE::psi_ori_1(int i)
{
    return difx_ori[i] * dist_ori_1;
}

double DIE::psi_ori_2(int i, int j)
{
    return d[i][j] * phi_ori() + difx_ori[i] * phi_ori_1(j);
}

double DIE::psi_ori_3(int i, int j, int k)
{
    return d[i][j] * phi_ori_1(k) + d[i][k] * phi_ori_1(j) + difx_ori[i] * phi_ori_2(j, k);
}

double DIE::psi_ori_4(int i, int j, int k, int l)
{
    return d[i][j] * phi_ori_2(k, l) + d[i][k] * phi_ori_2(j, l) + d[i][l] * phi_ori_2(j, k) + difx_ori[i] * phi_ori_3(j, k, l);
}

double DIE::psi_ori_5(int i, int j, int k, int l, int m)
{
    return difx_ori[i] * phi_ori_4(j, k, l, m) + d[i][j] * phi_ori_3(k, l, m)
        + d[i][k] * phi_ori_3(j, l, m) + d[i][l] * phi_ori_3(j, k, m) + d[i][m] * phi_ori_3(j, k, l);
}

double DIE::psi_ori_6(int i, int j, int k, int l, int m, int n)
{
    return difx_ori[i] * phi_ori_5(j, k, l, m, n) + d[i][j] * phi_ori_4(k, l, m, n)
        + d[i][k] * phi_ori_4(j, l, m, n) + d[i][l] * phi_ori_4(j, k, m, n)
        + d[i][m] * phi_ori_4(j, k, l, n) + d[i][n] * phi_ori_4(j, k, l, m);
}
