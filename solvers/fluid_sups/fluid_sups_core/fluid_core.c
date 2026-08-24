#include "fluid_core.h"

#include <stdio.h>


/* ==========================================================================
 * Global variables
 * ========================================================================== */

const char* CODENAME = "fluid >";


/*
 * These filenames are private to fluid_core.c.
 */
static const char* INPUT_FILENAME_NODE  = "node.dat";
static const char* INPUT_FILENAME_GRAPH = "graph.dat";
static const char* INPUT_FILENAME_ELEM  = "elem.dat";


/* ==========================================================================
 * Basic fluid utilities
 * ========================================================================== */

const char* BBFE_fluid_get_directory_name(
        int         argc,
        char*       argv[],
        const char* codename)
{
    const char* dir_name;

    if(argc < 2){
        dir_name = ".";
    }
    else{
        dir_name = argv[1];
    }

    printf(
        "%s Main directory: %s\n",
        codename,
        dir_name);

    return dir_name;
}


void BBFE_fluid_pre(
        BBFE_DATA*    fe,
        BBFE_BASIS*   basis,
        int           argc,
        char*         argv[],
        const char*   directory,
        int           num_integ_points_each_axis)
{
    /*
     * argc / argv are retained for API compatibility.
     */
    (void)argc;
    (void)argv;

    BB_calc_void();

    const int n_axis =
        num_integ_points_each_axis;

    const char* filename;


    filename =
        monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            INPUT_FILENAME_NODE);

    BBFE_sys_read_node(
        fe,
        filename,
        directory);


    filename =
        monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            INPUT_FILENAME_ELEM);

    BBFE_sys_read_elem(
        fe,
        filename,
        directory,
        n_axis * n_axis * n_axis);


    BBFE_sys_memory_allocation_integ(
        basis,
        n_axis * n_axis * n_axis,
        3);


    BBFE_sys_memory_allocation_shapefunc(
        basis,
        fe->local_num_nodes,
        1,
        n_axis * n_axis * n_axis);


    BBFE_fluid_set_basis(
        basis,
        fe->local_num_nodes,
        n_axis);
}


void BBFE_fluid_set_basis(
        BBFE_BASIS*   basis,
        int           local_num_nodes,
        int           num_integ_points_each_axis)
{
    switch(local_num_nodes){

    case 4:

        basis->num_integ_points =
            BBFE_std_integ_tet_set_arbitrary_points(
                num_integ_points_each_axis,
                basis->integ_point,
                basis->integ_weight);


        for(int i = 0;
            i < basis->num_integ_points;
            ++i){

            BBFE_std_shapefunc_tet1st_get_val(
                basis->integ_point[i],
                basis->N[i]);


            BBFE_std_shapefunc_tet1st_get_derivative(
                basis->integ_point[i],
                basis->dN_dxi[i],
                basis->dN_det[i],
                basis->dN_dze[i]);
        }


        printf(
            "%s Element type: "
            "1st-order tetrahedron.\n",
            CODENAME);

        break;


    case 8:

        basis->num_integ_points =
            BBFE_std_integ_hex_set_arbitrary_points(
                num_integ_points_each_axis,
                basis->integ_point,
                basis->integ_weight);


        for(int i = 0;
            i < basis->num_integ_points;
            ++i){

            BBFE_std_shapefunc_hex1st_get_val(
                basis->integ_point[i],
                basis->N[i]);


            BBFE_std_shapefunc_hex1st_get_derivative(
                basis->integ_point[i],
                basis->dN_dxi[i],
                basis->dN_det[i],
                basis->dN_dze[i]);
        }


        printf(
            "%s Element type: "
            "1st-order hexahedron.\n",
            CODENAME);

        break;


    default:

        fprintf(
            stderr,
            "%s ERROR: unsupported local_num_nodes=%d "
            "in BBFE_fluid_set_basis().\n",
            CODENAME,
            local_num_nodes);

        exit(EXIT_FAILURE);
    }


    printf(
        "%s The number of integration points: %d\n",
        CODENAME,
        basis->num_integ_points);
}


void BBFE_fluid_sups_renew_velocity(
        double**  v,
        double*   ans_vec,
        const int total_num_nodes)
{
    for(int i = 0;
        i < total_num_nodes;
        ++i){

        for(int d = 0; d < 3; ++d){

            v[i][d] =
                ans_vec[4*i + d];
        }
    }
}


void BBFE_fluid_sups_renew_pressure(
        double*   p,
        double*   ans_vec,
        const int total_num_nodes)
{
    for(int i = 0;
        i < total_num_nodes;
        ++i){

        p[i] =
            ans_vec[4*i + 3];
    }
}


void BBFE_fluid_sups_update_vec(
        double**  v,
        double*   p,
        double*   ans_vec,
        const int total_num_nodes)
{
    for(int i = 0;
        i < total_num_nodes;
        ++i){

        ans_vec[4*i + 0] = v[i][0];
        ans_vec[4*i + 1] = v[i][1];
        ans_vec[4*i + 2] = v[i][2];
        ans_vec[4*i + 3] = p[i];
    }
}


void BBFE_fluid_finalize(
        BBFE_DATA*   fe,
        BBFE_BASIS*  basis)
{
    BBFE_sys_memory_free_integ(
        basis,
        3);


    BBFE_sys_memory_free_shapefunc(
        basis);


    BBFE_sys_memory_free_node(
        fe,
        3);


    BBFE_sys_memory_free_elem(
        fe,
        basis->num_integ_points,
        3);
}


/* ==========================================================================
 * Shared FOM / NR / ST-FOM API
 * ========================================================================== */


/*
 * Generic zero initial condition.
 */
void BBFE_fluid_initialize_velocity_pressure(
        double**    v,
        double*     p,
        const int   total_num_nodes)
{
    if(v == NULL || p == NULL){
        return;
    }


    for(int i = 0;
        i < total_num_nodes;
        ++i){

        v[i][0] = 0.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;

        p[i] = 0.0;
    }
}


/*
 * Karman-vortex initial condition.
 *
 * Existing core_FOM.c uses:
 *
 *   u_x = 1
 *   u_y = 0
 *   u_z = 0
 *   p   = 0
 *
 * v_old is initialized consistently with v.
 */
void BBFE_fluid_initialize_velocity_pressure_karman_vortex(
        double**    v,
        double**    v_old,
        double*     p,
        const int   total_num_nodes)
{
    if(v == NULL || p == NULL){
        return;
    }


    for(int i = 0;
        i < total_num_nodes;
        ++i){

        v[i][0] = 1.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;


        if(v_old != NULL){

            v_old[i][0] = v[i][0];
            v_old[i][1] = v[i][1];
            v_old[i][2] = v[i][2];
        }


        p[i] = 0.0;
    }
}


/*
 * Pack velocity + pressure into a 4-DOF vector:
 *
 *   [ux0, uy0, uz0, p0,
 *    ux1, uy1, uz1, p1,
 *    ...]
 */
void BBFE_fluid_sups_add_velocity_pressure(
        double**    v,
        double*     p,
        double*     sol_vec,
        const int   total_num_nodes)
{
    if(v == NULL ||
       p == NULL ||
       sol_vec == NULL){

        return;
    }


    for(int i = 0;
        i < total_num_nodes;
        ++i){

        sol_vec[4*i + 0] = v[i][0];
        sol_vec[4*i + 1] = v[i][1];
        sol_vec[4*i + 2] = v[i][2];
        sol_vec[4*i + 3] = p[i];
    }
}


/*
 * Apply prescribed Dirichlet values to a flattened vector.
 */
void BBFE_fluid_add_Dbc(
        double*         ansvec,
        BBFE_BC*        bc,
        const int       total_num_nodes,
        const int       dof)
{
    if(ansvec == NULL ||
       bc == NULL ||
       bc->D_bc_exists == NULL ||
       bc->imposed_D_val == NULL){

        return;
    }


    if(dof <= 0){
        return;
    }


    for(int i = 0;
        i < total_num_nodes;
        ++i){

        for(int d = 0;
            d < dof;
            ++d){

            /*
             * Avoid reading beyond the BC block.
             */
            if(d >= bc->block_size){
                continue;
            }


            const int bc_index =
                bc->block_size * i + d;

            const int vec_index =
                dof * i + d;


            if(bc->D_bc_exists[bc_index]){

                ansvec[vec_index] =
                    bc->imposed_D_val[bc_index];
            }
        }
    }
}


/*
 * Newton update:
 *
 *   v <- v + delta_v
 *   p <- p + delta_p
 */
void BBFE_fluid_update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes)
{
    if(v == NULL ||
       delta_v == NULL ||
       p == NULL ||
       delta_p == NULL){

        return;
    }


    for(int i = 0;
        i < total_num_nodes;
        ++i){

        for(int d = 0;
            d < 3;
            ++d){

            v[i][d] +=
                delta_v[i][d];
        }


        p[i] +=
            delta_p[i];
    }
}


/*
 * Local squared L2 norm.
 *
 * sqrt() is intentionally NOT applied here.
 * This allows MPI SUM reduction before taking the square root.
 */
double BBFE_fluid_calc_internal_norm_2d(
        double**    in,
        const int   num1,
        const int   num2)
{
    if(in == NULL){
        return 0.0;
    }


    double norm2 = 0.0;


    for(int i = 0;
        i < num1;
        ++i){

        for(int j = 0;
            j < num2;
            ++j){

            const double x =
                in[i][j];

            norm2 +=
                x * x;
        }
    }


    return norm2;
}


double BBFE_fluid_calc_internal_norm_1d(
        double*     in,
        const int   num1,
        const int   num2)
{
    /*
     * Retained for compatibility with the historical API.
     */
    (void)num2;


    if(in == NULL){
        return 0.0;
    }


    double norm2 = 0.0;


    for(int i = 0;
        i < num1;
        ++i){

        const double x =
            in[i];

        norm2 +=
            x * x;
    }


    return norm2;
}


/*
 * Copy 2-D array.
 *
 * Historical argument order:
 *
 *   in  -> source
 *   out -> destination
 */
void BBFE_fluid_vec_copy_2d(
        double**    in,
        double**    out,
        const int   num1,
        const int   num2)
{
    if(in == NULL ||
       out == NULL){

        return;
    }


    for(int i = 0;
        i < num1;
        ++i){

        for(int j = 0;
            j < num2;
            ++j){

            out[i][j] =
                in[i][j];
        }
    }
}


/*
 * Copy 1-D array.
 *
 * dst <- src
 */
void BBFE_fluid_vec_copy_1d(
        double*       dst,
        const double* src,
        int           size)
{
    if(dst == NULL ||
       src == NULL ||
       size <= 0){

        return;
    }


    for(int i = 0;
        i < size;
        ++i){

        dst[i] =
            src[i];
    }
}


/*
 * Residual layout:
 *
 *   rvec[4*i + 0] : x momentum
 *   rvec[4*i + 1] : y momentum
 *   rvec[4*i + 2] : z momentum
 *   rvec[4*i + 3] : continuity / mass
 *
 * Return LOCAL squared norms.
 *
 * The caller can perform MPI_SUM and sqrt afterwards.
 */
double BBFE_fluid_calc_internal_norm_total(
        const double* rvec,
        int           n_internal_vertex)
{
    if(rvec == NULL ||
       n_internal_vertex <= 0){

        return 0.0;
    }


    double norm2 = 0.0;


    for(int i = 0;
        i < n_internal_vertex;
        ++i){

        for(int d = 0;
            d < 4;
            ++d){

            const double x =
                rvec[4*i + d];

            norm2 +=
                x * x;
        }
    }


    return norm2;
}


double BBFE_fluid_calc_internal_norm_momentum(
        const double* rvec,
        int           n_internal_vertex)
{
    if(rvec == NULL ||
       n_internal_vertex <= 0){

        return 0.0;
    }


    double norm2 = 0.0;


    for(int i = 0;
        i < n_internal_vertex;
        ++i){

        for(int d = 0;
            d < 3;
            ++d){

            const double x =
                rvec[4*i + d];

            norm2 +=
                x * x;
        }
    }


    return norm2;
}


double BBFE_fluid_calc_internal_norm_mass(
        const double* rvec,
        int           n_internal_vertex)
{
    if(rvec == NULL ||
       n_internal_vertex <= 0){

        return 0.0;
    }


    double norm2 = 0.0;


    for(int i = 0;
        i < n_internal_vertex;
        ++i){

        const double x =
            rvec[4*i + 3];

        norm2 +=
            x * x;
    }


    return norm2;
}