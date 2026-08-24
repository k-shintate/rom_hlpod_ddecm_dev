#include "core_FOM_mixed.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FOM_MIXED_BUFFER_SIZE 10000

static const char* FOM_MIXED_OUTPUT_VTK    = "result_%06d.vtk";
static const char* FOM_MIXED_OUTPUT_TEMP   = "temparature_%06d.dat";
static const char* FOM_MIXED_OUTPUT_SOURCE = "source_%06d.dat";

/* ------------------------------------------------------------------------- */
/* Prism reference element                                                   */
/* ------------------------------------------------------------------------- */

/*
 * Map tensor-product Gauss points (u,v,w) in [-1,1]^3 to the reference
 * triangular prism
 *
 *     0 <= xi, 0 <= eta, xi + eta <= 1,  -1 <= zeta <= 1
 *
 * by the Duffy transform
 *
 *     xi   = (1 + u)/2
 *     eta  = (1 - u)(1 + v)/4
 *     zeta = w
 *
 * whose Jacobian is (1-u)/8.
 *
 * This fixes the old experimental implementation that overwrote every prism
 * quadrature weight with 1/6: that constant does not preserve the tensor
 * Gauss rule and gives an incorrect integral for n_axis > 1.
 */
int FOM_mixed_integ_prism_set_arbitrary_points(
        int      num_points_in_each_axis,
        double** integ_point,
        double*  integ_weight)
{
    const int n = num_points_in_each_axis;
    const int num_integ_points = BBFE_std_integ_hex_set_arbitrary_points(
            n, integ_point, integ_weight);

    for(int g = 0; g < num_integ_points; g++) {
        const double u = integ_point[g][0];
        const double v = integ_point[g][1];
        const double w = integ_point[g][2];

        const double weight_cube = integ_weight[g];
        const double jac_duffy = (1.0 - u) / 8.0;

        integ_point[g][0] = 0.5 * (1.0 + u);
        integ_point[g][1] = 0.25 * (1.0 - u) * (1.0 + v);
        integ_point[g][2] = w;

        integ_weight[g] = weight_cube * jac_duffy;
    }

    return num_integ_points;
}


void FOM_mixed_shapefunc_prism1st_get_val(
        const double xi[3],
        double*      N)
{
    const double r = xi[0];
    const double s = xi[1];
    const double z = xi[2];
    const double c = 0.5;

    N[0] = c * (1.0 - r - s) * (1.0 - z);
    N[1] = c * r               * (1.0 - z);
    N[2] = c * s               * (1.0 - z);
    N[3] = c * (1.0 - r - s) * (1.0 + z);
    N[4] = c * r               * (1.0 + z);
    N[5] = c * s               * (1.0 + z);
}


void FOM_mixed_shapefunc_prism1st_get_derivative(
        const double xi[3],
        double*      dN_dxi,
        double*      dN_deta,
        double*      dN_dzeta)
{
    const double r = xi[0];
    const double s = xi[1];
    const double z = xi[2];
    const double c = 0.5;

    dN_dxi[0] = -c * (1.0 - z);
    dN_dxi[1] =  c * (1.0 - z);
    dN_dxi[2] =  0.0;
    dN_dxi[3] = -c * (1.0 + z);
    dN_dxi[4] =  c * (1.0 + z);
    dN_dxi[5] =  0.0;

    dN_deta[0] = -c * (1.0 - z);
    dN_deta[1] =  0.0;
    dN_deta[2] =  c * (1.0 - z);
    dN_deta[3] = -c * (1.0 + z);
    dN_deta[4] =  0.0;
    dN_deta[5] =  c * (1.0 + z);

    dN_dzeta[0] = -c * (1.0 - r - s);
    dN_dzeta[1] = -c * r;
    dN_dzeta[2] = -c * s;
    dN_dzeta[3] =  c * (1.0 - r - s);
    dN_dzeta[4] =  c * r;
    dN_dzeta[5] =  c * s;
}


void BBFE_convdiff_set_basis_prism(
        BBFE_BASIS* basis,
        int         num_integ_points_each_axis)
{
    basis->num_integ_points = FOM_mixed_integ_prism_set_arbitrary_points(
            num_integ_points_each_axis,
            basis->integ_point,
            basis->integ_weight);

    for(int p = 0; p < basis->num_integ_points; p++) {
        FOM_mixed_shapefunc_prism1st_get_val(
                basis->integ_point[p],
                basis->N[p]);

        FOM_mixed_shapefunc_prism1st_get_derivative(
                basis->integ_point[p],
                basis->dN_dxi[p],
                basis->dN_det[p],
                basis->dN_dze[p]);
    }

    printf("%s Element type: 1st-order prism.\n", CODENAME);
    printf("%s Prism integration points: %d\n",
            CODENAME, basis->num_integ_points);
}

/* ------------------------------------------------------------------------- */
/* Mixed input                                                               */
/* ------------------------------------------------------------------------- */

static void FOM_mixed_check_shared_nodes(
        const BBFE_DATA* fe_tet,
        const BBFE_DATA* fe_pri)
{
    if(fe_tet->total_num_nodes != fe_pri->total_num_nodes) {
        fprintf(stderr,
                "%s ERROR: tetra/prism node counts differ: tet=%d pri=%d\n",
                CODENAME,
                fe_tet->total_num_nodes,
                fe_pri->total_num_nodes);
        exit(EXIT_FAILURE);
    }

    if(fe_tet->local_num_nodes <= 0 || fe_pri->local_num_nodes <= 0) {
        fprintf(stderr,
                "%s ERROR: invalid element node count: tet=%d pri=%d\n",
                CODENAME,
                fe_tet->local_num_nodes,
                fe_pri->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    if(fe_tet->local_num_nodes != 4) {
        fprintf(stderr,
                "%s ERROR: elem_tet.dat must contain 4-node tetrahedra; got %d nodes/element.\n",
                CODENAME, fe_tet->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    if(fe_pri->local_num_nodes != 6) {
        fprintf(stderr,
                "%s ERROR: elem_prism.dat must contain 6-node prisms; got %d nodes/element.\n",
                CODENAME, fe_pri->local_num_nodes);
        exit(EXIT_FAILURE);
    }
}


void BBFE_convdiff_pre_mixed(
        BBFE_DATA*    fe_tet,
        BBFE_BASIS*   basis_tet,
        BBFE_DATA*    fe_pri,
        BBFE_BASIS*   basis_pri,
        BBFE_BC*      bc,
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        const char*   directory,
        int           num_integ_points_each_axis)
{
    BB_calc_void();

    const int n_axis = num_integ_points_each_axis;
    const int num_ip_alloc = n_axis * n_axis * n_axis;
    const char* filename;

    /* shared node numbering; each FE object owns its own coordinate array */
    filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            INPUT_FILENAME_NODE);

    BBFE_sys_read_node(fe_tet, filename, directory);
    BBFE_sys_read_node(fe_pri, filename, directory);

    /* tetrahedra */
    /*
    filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            FOM_MIXED_ELEM_TET_FILE);

    BBFE_sys_read_elem(
            fe_tet,
            filename,
            directory,
            num_ip_alloc);
    */
        
    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, FOM_MIXED_ELEM_TET_FILE);

    read_connectivity_graph_surf(
                    fe_tet,
                    fname,
                    directory,
                    num_ip_alloc);

    /* prisms */
    /*
    filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            FOM_MIXED_ELEM_PRI_FILE);

    BBFE_sys_read_elem(
            fe_pri,
            filename,
            directory,
            num_ip_alloc);
    */
        
    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, FOM_MIXED_ELEM_PRI_FILE);

    read_connectivity_graph_surf(
                    fe_pri,
                    fname,
                    directory,
                    num_ip_alloc);

    FOM_mixed_check_shared_nodes(fe_tet, fe_pri);

    /* tetra integration and shape functions */
    BBFE_sys_memory_allocation_integ(
            basis_tet,
            num_ip_alloc,
            3);

    BBFE_sys_memory_allocation_shapefunc(
            basis_tet,
            fe_tet->local_num_nodes,
            1,
            num_ip_alloc);

    BBFE_convdiff_set_basis(
            basis_tet,
            fe_tet->local_num_nodes,
            n_axis);

    /* prism integration and shape functions */
    BBFE_sys_memory_allocation_integ(
            basis_pri,
            num_ip_alloc,
            3);

    BBFE_sys_memory_allocation_shapefunc(
            basis_pri,
            fe_pri->local_num_nodes,
            1,
            num_ip_alloc);

    BBFE_convdiff_set_basis_prism(
            basis_pri,
            n_axis);

    /* one BC array because the node numbering is common */
    filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            INPUT_FILENAME_D_BC);

    BBFE_sys_read_Dirichlet_bc(
            bc,
            filename,
            directory,
            fe_tet->total_num_nodes,
            BLOCK_SIZE);

    /* matrix object used during time stepping */
    monolis_initialize(monolis);

    /* communication tables produced by the graph partitioning preprocessing */
    monolis_com_initialize_by_parted_files(
            monolis_com,
            monolis_mpi_get_global_comm(),
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            INPUT_FILENAME_GRAPH);

    printf("%s Mixed mesh: nodes=%d, tetra=%d, prism=%d\n",
            CODENAME,
            fe_tet->total_num_nodes,
            fe_tet->total_num_elems,
            fe_pri->total_num_elems);
}

/* ------------------------------------------------------------------------- */
/* Unified nodal graph                                                       */
/* ------------------------------------------------------------------------- */

static void FOM_mixed_read_nodal_graph(
        const char* directory,
        const char* filename,
        int*        num_nodes_out,
        int**       index_out,
        int**       item_out)
{
    FILE* fp = NULL;
    int num_nodes = 0;
    long long num_items_ll = 0;

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    if(fscanf(fp, "%d", &num_nodes) != 1 || num_nodes <= 0) {
        fprintf(stderr, "%s ERROR: invalid graph header in %s.\n",
                CODENAME, filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    for(int row = 0; row < num_nodes; row++) {
        int node_id = -1;
        int degree = -1;

        if(fscanf(fp, "%d %d", &node_id, &degree) != 2 || degree < 0) {
            fprintf(stderr, "%s ERROR: invalid graph row %d in %s.\n",
                    CODENAME, row, filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        for(int j = 0; j < degree; j++) {
            int dummy;
            if(fscanf(fp, "%d", &dummy) != 1) {
                fprintf(stderr, "%s ERROR: truncated adjacency at row %d in %s.\n",
                        CODENAME, row, filename);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
        }

        (void)node_id;
        num_items_ll += degree;
    }
    fclose(fp);

    if(num_items_ll > 2147483647LL) {
        fprintf(stderr, "%s ERROR: graph adjacency is too large for int CSR.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    const int num_items = (int)num_items_ll;
    int* index = NULL;
    int* item = NULL;

    index = BB_std_calloc_1d_int(index, num_nodes + 1);
    item  = BB_std_calloc_1d_int(item, num_items);

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    fscanf(fp, "%d", &num_nodes);

    int cursor = 0;
    index[0] = 0;

    for(int row = 0; row < num_nodes; row++) {
        int node_id = -1;
        int degree = -1;
        fscanf(fp, "%d %d", &node_id, &degree);

        for(int j = 0; j < degree; j++) {
            int adj = -1;
            fscanf(fp, "%d", &adj);

            if(adj < 0 || adj >= num_nodes) {
                fprintf(stderr,
                        "%s ERROR: graph adjacency out of range: row=%d adj=%d nnode=%d.\n",
                        CODENAME, row, adj, num_nodes);
                fclose(fp);
                BB_std_free_1d_int(index, num_nodes + 1);
                BB_std_free_1d_int(item, num_items);
                exit(EXIT_FAILURE);
            }

            item[cursor++] = adj;
        }

        index[row + 1] = cursor;
        (void)node_id;
    }

    fclose(fp);

    *num_nodes_out = num_nodes;
    *index_out = index;
    *item_out = item;
}


void BBFE_convdiff_set_nonzero_pattern_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_ref,
        const char* directory)
{
    const char* filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            FOM_MIXED_GRAPH_FILE);

    int num_nodes = 0;
    int* index = NULL;
    int* item = NULL;

    FOM_mixed_read_nodal_graph(
            directory,
            filename,
            &num_nodes,
            &index,
            &item);

    const int num_items = index[num_nodes];

    if(num_nodes != fe_ref->total_num_nodes) {
        fprintf(stderr,
                "%s ERROR: graph/node mismatch: graph=%d node.dat=%d (%s).\n",
                CODENAME,
                num_nodes,
                fe_ref->total_num_nodes,
                filename);
        BB_std_free_1d_int(index, num_nodes + 1);
        BB_std_free_1d_int(item, num_items);
        exit(EXIT_FAILURE);
    }

    monolis_get_nonzero_pattern_by_nodal_graph_R(
            monolis,
            num_nodes,
            BLOCK_SIZE,
            index,
            item);

    printf("%s Unified tetra+prism nodal graph: nodes=%d directed_adjacencies=%d\n",
            CODENAME, num_nodes, num_items);

    BB_std_free_1d_int(index, num_nodes + 1);
    BB_std_free_1d_int(item, num_items);
}

/* ------------------------------------------------------------------------- */
/* Mixed assembly / solve                                                    */
/* ------------------------------------------------------------------------- */

void set_element_mat_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        VALUES*     vals)
{
    set_element_mat(monolis, fe_tet, basis_tet, vals);
    set_element_mat(monolis, fe_pri, basis_pri, vals);
}


void set_element_vec_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        VALUES*     vals,
        double      t)
{
    set_element_vec(monolis, fe_tet, basis_tet, vals, t);
    set_element_vec(monolis, fe_pri, basis_pri, vals, t);
}


void solver_fom_mixed(
        FE_SYSTEM_MIXED* sys,
        double           t,
        int              step)
{
    printf("\n%s ----------------- mixed FOM step %d ----------------\n",
            CODENAME, step);

    monolis_copy_mat_value_R(&(sys->monolis0), &(sys->monolis));
    monolis_clear_mat_value_rhs_R(&(sys->monolis));

    if(manusol_rank_benchmark_is_active()) {
        /* Preserve the algebraic manufactured benchmark behavior. */
        manusol_set_theo_sol(&(sys->fe_tet), sys->vals.theo_sol, t);

        monolis_matvec_product_R(
                &(sys->monolis),
                &(sys->monolis_com),
                sys->vals.theo_sol,
                sys->monolis.mat.R.B);
    }
    else {
        set_element_vec_mixed(
                &(sys->monolis),
                &(sys->fe_tet),
                &(sys->basis_tet),
                &(sys->fe_pri),
                &(sys->basis_pri),
                &(sys->vals),
                t);

        /* node.dat is common, so either FE object is valid for nodal fields */
        manusol_set_theo_sol(&(sys->fe_tet), sys->vals.theo_sol, t);
    }

    BBFE_manusol_set_bc_scalar(
            &(sys->fe_tet),
            &(sys->bc),
            sys->vals.theo_sol,
            t);

    BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            sys->fe_tet.total_num_nodes,
            BLOCK_SIZE,
            &(sys->bc),
            sys->monolis.mat.R.B);

    BBFE_sys_monowrap_solve(
            &(sys->monolis),
            &(sys->monolis_com),
            sys->vals.T,
            MONOLIS_ITER_CG,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);
}

/* ------------------------------------------------------------------------- */
/* Mixed VTK / L2 output                                                     */
/* ------------------------------------------------------------------------- */

static void FOM_mixed_vtk_write_cells(
        FILE* fp,
        int   num_cells,
        int   num_points_in_cell,
        int** connectivity)
{
    for(int e = 0; e < num_cells; e++) {
        fprintf(fp, "%d", num_points_in_cell);
        for(int i = 0; i < num_points_in_cell; i++) {
            fprintf(fp, " %d", connectivity[e][i]);
        }
        fprintf(fp, "\n");
    }
}


static void FOM_mixed_vtk_write_cell_types(
        FILE* fp,
        int   num_cells,
        int   vtk_type)
{
    for(int e = 0; e < num_cells; e++) {
        fprintf(fp, "%d\n", vtk_type);
    }
}


static void FOM_mixed_vtk_write_shape(
        FILE*      fp,
        BBFE_DATA* fe_tet,
        BBFE_DATA* fe_pri)
{
    const int ncell = fe_tet->total_num_elems + fe_pri->total_num_elems;
    const int cell_storage =
        fe_tet->total_num_elems * (fe_tet->local_num_nodes + 1)
        + fe_pri->total_num_elems * (fe_pri->local_num_nodes + 1);

    BB_vtk_write_header(fp);
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    BB_vtk_write_points_3d(fp, fe_tet->total_num_nodes, fe_tet->x);

    fprintf(fp, "CELLS %d %d\n", ncell, cell_storage);
    FOM_mixed_vtk_write_cells(
            fp,
            fe_tet->total_num_elems,
            fe_tet->local_num_nodes,
            fe_tet->conn);
    FOM_mixed_vtk_write_cells(
            fp,
            fe_pri->total_num_elems,
            fe_pri->local_num_nodes,
            fe_pri->conn);

    fprintf(fp, "CELL_TYPES %d\n", ncell);
    FOM_mixed_vtk_write_cell_types(fp, fe_tet->total_num_elems, TYPE_VTK_TETRA);
    FOM_mixed_vtk_write_cell_types(fp, fe_pri->total_num_elems, TYPE_VTK_WEDGE);
}


static void FOM_mixed_l2_accumulate_local(
        BBFE_DATA*    fe,
        BBFE_BASIS*   basis,
        MONOLIS_COM*  monolis_com,
        double        t,
        const double* comp_vec,
        double        (*func)(double, double, double, double),
        double*       abs_error,
        double*       abs_theo)
{
    const int np = basis->num_integ_points;
    const int nl = fe->local_num_nodes;

    double* val_ip_error = NULL;
    double* val_ip_theo = NULL;
    double* Jacobian_ip = NULL;
    bool* is_internal_elem = NULL;
    double** local_x = NULL;
    double* local_val = NULL;

    val_ip_error = BB_std_calloc_1d_double(val_ip_error, np);
    val_ip_theo  = BB_std_calloc_1d_double(val_ip_theo, np);
    Jacobian_ip  = BB_std_calloc_1d_double(Jacobian_ip, np);
    is_internal_elem = BB_std_calloc_1d_bool(
            is_internal_elem, fe->total_num_elems);
    local_x = BB_std_calloc_2d_double(local_x, nl, 3);
    local_val = BB_std_calloc_1d_double(local_val, nl);

    monolis_get_bool_list_of_internal_simple_mesh(
            monolis_com,
            fe->total_num_nodes,
            fe->total_num_elems,
            nl,
            fe->conn,
            is_internal_elem);

    for(int e = 0; e < fe->total_num_elems; e++) {
        if(!is_internal_elem[e]) {
            continue;
        }

        for(int i = 0; i < nl; i++) {
            const int node = fe->conn[e][i];
            local_val[i] = comp_vec[node];
            local_x[i][0] = fe->x[node][0];
            local_x[i][1] = fe->x[node][1];
            local_x[i][2] = fe->x[node][2];
        }

        BBFE_elemmat_set_Jacobian_array(
                Jacobian_ip,
                np,
                e,
                fe);

        for(int p = 0; p < np; p++) {
            const double value = BBFE_std_mapping_scalar(
                    nl, local_val, basis->N[p]);

            double x_ip[3];
            BBFE_std_mapping_vector3d(
                    x_ip, nl, local_x, basis->N[p]);

            const double exact = func(x_ip[0], x_ip[1], x_ip[2], t);
            const double diff = value - exact;

            val_ip_error[p] = diff * diff;
            val_ip_theo[p]  = exact * exact;
        }

        *abs_error += BBFE_std_integ_calc(
                np,
                val_ip_error,
                basis->integ_weight,
                Jacobian_ip);

        *abs_theo += BBFE_std_integ_calc(
                np,
                val_ip_theo,
                basis->integ_weight,
                Jacobian_ip);
    }

    BB_std_free_1d_double(val_ip_error, np);
    BB_std_free_1d_double(val_ip_theo, np);
    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);
    BB_std_free_2d_double(local_x, nl, 3);
    BB_std_free_1d_double(local_val, nl);
}


double BBFE_convdiff_equivval_relative_L2_error_scalar_mixed(
        BBFE_DATA*    fe_tet,
        BBFE_BASIS*   basis_tet,
        BBFE_DATA*    fe_pri,
        BBFE_BASIS*   basis_pri,
        MONOLIS_COM*  monolis_com,
        double        t,
        const double* comp_vec,
        double        (*func)(double, double, double, double))
{
    double abs_error = 0.0;
    double abs_theo = 0.0;

    FOM_mixed_l2_accumulate_local(
            fe_tet,
            basis_tet,
            monolis_com,
            t,
            comp_vec,
            func,
            &abs_error,
            &abs_theo);

    FOM_mixed_l2_accumulate_local(
            fe_pri,
            basis_pri,
            monolis_com,
            t,
            comp_vec,
            func,
            &abs_error,
            &abs_theo);

    monolis_allreduce_R(
            1,
            &abs_error,
            MONOLIS_MPI_SUM,
            monolis_com->comm);

    monolis_allreduce_R(
            1,
            &abs_theo,
            MONOLIS_MPI_SUM,
            monolis_com->comm);

    if(abs_theo <= 0.0) {
        return sqrt(abs_error);
    }

    return sqrt(abs_error / abs_theo);
}


void output_result_file_vtk_mixed(
        BBFE_DATA*  fe_tet,
        BBFE_DATA*  fe_pri,
        VALUES*     vals,
        const char* filename,
        const char* directory,
        double      t)
{
    FILE* fp = NULL;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    FOM_mixed_vtk_write_shape(fp, fe_tet, fe_pri);

    fprintf(fp, "POINT_DATA %d\n", fe_tet->total_num_nodes);
    BB_vtk_write_point_vals_scalar(
            fp, vals->T, fe_tet->total_num_nodes, "temperature");

    BBFE_manusol_calc_nodal_error_scalar(
            fe_tet, vals->error, vals->theo_sol, vals->T);
    BB_vtk_write_point_vals_scalar(
            fp, vals->error, fe_tet->total_num_nodes, "abs_error");
    BB_vtk_write_point_vals_scalar(
            fp, vals->theo_sol, fe_tet->total_num_nodes, "theoretical");

    double* source = NULL;
    source = BB_std_calloc_1d_double(source, fe_tet->total_num_nodes);
    manusol_set_source(fe_tet, source, t);
    BB_vtk_write_point_vals_scalar(
            fp, source, fe_tet->total_num_nodes, "source");
    BB_std_free_1d_double(source, fe_tet->total_num_nodes);

    fclose(fp);
}


void output_files_mixed(
        FE_SYSTEM_MIXED* sys,
        int              file_num,
        double           t)
{
    const char* filename;
    char fname_vtk[FOM_MIXED_BUFFER_SIZE];
    char fname_temp[FOM_MIXED_BUFFER_SIZE];
    char fname_source[FOM_MIXED_BUFFER_SIZE];

    snprintf(fname_vtk, sizeof(fname_vtk), FOM_MIXED_OUTPUT_VTK, file_num);
    snprintf(fname_temp, sizeof(fname_temp), FOM_MIXED_OUTPUT_TEMP, file_num);
    snprintf(fname_source, sizeof(fname_source), FOM_MIXED_OUTPUT_SOURCE, file_num);

    filename = monolis_get_global_output_file_name(
            MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

    output_result_file_vtk_mixed(
            &(sys->fe_tet),
            &(sys->fe_pri),
            &(sys->vals),
            filename,
            sys->cond.directory,
            t);

    filename = monolis_get_global_output_file_name(
            MONOLIS_DEFAULT_TOP_DIR, "./", fname_temp);

    BBFE_write_ascii_nodal_vals_scalar(
            &(sys->fe_tet),
            sys->vals.T,
            filename,
            sys->cond.directory);

    double* source = NULL;
    source = BB_std_calloc_1d_double(source, sys->fe_tet.total_num_nodes);
    manusol_set_source(&(sys->fe_tet), source, t);

    filename = monolis_get_global_output_file_name(
            MONOLIS_DEFAULT_TOP_DIR, "./", fname_source);

    BBFE_write_ascii_nodal_vals_scalar(
            &(sys->fe_tet),
            source,
            filename,
            sys->cond.directory);

    BB_std_free_1d_double(source, sys->fe_tet.total_num_nodes);

    const double l2_error = BBFE_convdiff_equivval_relative_L2_error_scalar_mixed(
            &(sys->fe_tet),
            &(sys->basis_tet),
            &(sys->fe_pri),
            &(sys->basis_pri),
            &(sys->monolis_com),
            t,
            sys->vals.T,
            manusol_get_sol);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf("%s Mixed tetra+prism relative L2 error: %e\n", CODENAME, l2_error);

        FILE* fp = NULL;
        fp = BBFE_sys_write_add_fopen(fp, "l2_error.txt", sys->cond.directory);
        fprintf(fp, "%e %e\n", t, l2_error);
        fclose(fp);
    }
}

/* ------------------------------------------------------------------------- */
/* Cleanup                                                                   */
/* ------------------------------------------------------------------------- */

void BBFE_convdiff_finalize_mixed(
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        BBFE_BC*    bc)
{
    const int total_num_nodes = fe_tet->total_num_nodes;
    const int num_ip_tet = basis_tet->num_integ_points;
    const int num_ip_pri = basis_pri->num_integ_points;

    BBFE_sys_memory_free_elem(fe_tet, num_ip_tet, 3);
    BBFE_sys_memory_free_node(fe_tet, 3);
    BBFE_sys_memory_free_integ(basis_tet, 3);
    BBFE_sys_memory_free_shapefunc(basis_tet);

    BBFE_sys_memory_free_elem(fe_pri, num_ip_pri, 3);
    BBFE_sys_memory_free_node(fe_pri, 3);
    BBFE_sys_memory_free_integ(basis_pri, 3);
    BBFE_sys_memory_free_shapefunc(basis_pri);

    BBFE_sys_memory_free_Dirichlet_bc(
            bc,
            total_num_nodes,
            BLOCK_SIZE);
}
