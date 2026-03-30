
#include "set_reduced_matvec.h"

void set_element_reduced_rhs(
    double*			vec,
    double**		pod_modes,
    const int		index,
    const double	integ_val,
    const int		num_modes)
{
    for(int k = 0; k < num_modes; k++){
        double val = integ_val * pod_modes[index][k];
        vec[k] += val;
    }
}

void set_element_reduced_rhs_Dbc(
    BBFE_BC*     	bc,
    double*			vec,
    double**		pod_modes,
    const int		index_i,
    const int		index_j,
    const double	integ_val,
    const int		num_modes)
{
    if( bc->D_bc_exists[index_j]) {
        for(int k1 = 0; k1 < num_modes; k1++){
            double val = pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
            vec[k1] += - val;
        }
    }

}

void set_element_reduced_mat(
    double**		mat,
    double**		pod_modes,
    const int		index_i,
    const int		index_j,
    const double	integ_val,
    const int		num_modes)
{
    for(int k1 = 0; k1 < num_modes; k1++){
        for(int k2 = 0; k2 < num_modes; k2++){
            double val = pod_modes[index_i][k1] * integ_val * pod_modes[index_j][k2];
            mat[k1][k2] += val;
        }
    }
}