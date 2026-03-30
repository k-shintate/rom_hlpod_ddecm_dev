
#include "elemmat.h"

double BBFE_elemmat_mag_mat_curl(const double *curlNi, const double *curlNj, double k) {
    double dot = curlNi[0]*curlNj[0] + curlNi[1]*curlNj[1] + curlNi[2]*curlNj[2];
    return k * dot;
}

double BBFE_elemmat_mag_mat_mass(const double *Ni, const double *Nj, double a) {
    double dot = Ni[0]*Nj[0] + Ni[1]*Nj[1] + Ni[2]*Nj[2];
    return a * dot;
}
