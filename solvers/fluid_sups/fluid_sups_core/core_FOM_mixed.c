#include "core_FOM_mixed.h"
#include "core_NR.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FLUID_MIXED_IO_BUFFER 4096

static const char* FLUID_MIXED_NODE_FILE = "node.dat";

/* -------------------------------------------------------------------------- */
/* PRI6 reference element                                                     */
/* -------------------------------------------------------------------------- */

/*
 * Duffy map from [-1,1]^3 to
 *   0 <= r, 0 <= s, r+s <= 1, -1 <= z <= 1.
 *
 *   r = (1+u)/2
 *   s = (1-u)(1+v)/4
 *   z = w
 *
 * det(d(r,s,z)/d(u,v,w)) = (1-u)/8.
 *
 * Reusing BBFE's tensor-product Gauss points keeps the allocation convention
 * n_axis^3 used by the existing fluid code, while giving the correct prism
 * quadrature weights.
 */
int BBFE_fluid_mixed_integ_prism_set_arbitrary_points(
        int      num_points_in_each_axis,
        double** integ_point,
        double*  integ_weight)
{
    const int n = num_points_in_each_axis;
    const int np = BBFE_std_integ_hex_set_arbitrary_points(
        n, integ_point, integ_weight);

    for(int g = 0; g < np; ++g) {
        const double u = integ_point[g][0];
        const double v = integ_point[g][1];
        const double w = integ_point[g][2];
        const double weight_cube = integ_weight[g];

        integ_point[g][0] = 0.5 * (1.0 + u);
        integ_point[g][1] = 0.25 * (1.0 - u) * (1.0 + v);
        integ_point[g][2] = w;
        integ_weight[g] = weight_cube * (1.0 - u) / 8.0;
    }

    return np;
}

void BBFE_fluid_mixed_shapefunc_prism1st_get_val(
        const double xi[3],
        double*      N)
{
    const double r = xi[0];
    const double s = xi[1];
    const double z = xi[2];
    const double c = 0.5;

    N[0] = c * (1.0-r-s) * (1.0-z);
    N[1] = c * r           * (1.0-z);
    N[2] = c * s           * (1.0-z);
    N[3] = c * (1.0-r-s) * (1.0+z);
    N[4] = c * r           * (1.0+z);
    N[5] = c * s           * (1.0+z);
}

void BBFE_fluid_mixed_shapefunc_prism1st_get_derivative(
        const double xi[3],
        double*      dN_dxi,
        double*      dN_deta,
        double*      dN_dzeta)
{
    const double r = xi[0];
    const double s = xi[1];
    const double z = xi[2];
    const double c = 0.5;

    dN_dxi[0] = -c*(1.0-z);
    dN_dxi[1] =  c*(1.0-z);
    dN_dxi[2] =  0.0;
    dN_dxi[3] = -c*(1.0+z);
    dN_dxi[4] =  c*(1.0+z);
    dN_dxi[5] =  0.0;

    dN_deta[0] = -c*(1.0-z);
    dN_deta[1] =  0.0;
    dN_deta[2] =  c*(1.0-z);
    dN_deta[3] = -c*(1.0+z);
    dN_deta[4] =  0.0;
    dN_deta[5] =  c*(1.0+z);

    dN_dzeta[0] = -c*(1.0-r-s);
    dN_dzeta[1] = -c*r;
    dN_dzeta[2] = -c*s;
    dN_dzeta[3] =  c*(1.0-r-s);
    dN_dzeta[4] =  c*r;
    dN_dzeta[5] =  c*s;
}

void BBFE_fluid_set_basis_prism(
        BBFE_BASIS* basis,
        int         num_integ_points_each_axis)
{
    basis->num_integ_points =
        BBFE_fluid_mixed_integ_prism_set_arbitrary_points(
            num_integ_points_each_axis,
            basis->integ_point,
            basis->integ_weight);

    for(int p = 0; p < basis->num_integ_points; ++p) {
        BBFE_fluid_mixed_shapefunc_prism1st_get_val(
            basis->integ_point[p], basis->N[p]);

        BBFE_fluid_mixed_shapefunc_prism1st_get_derivative(
            basis->integ_point[p],
            basis->dN_dxi[p],
            basis->dN_det[p],
            basis->dN_dze[p]);
    }

    printf("%s Element type: 1st-order prism (PRI6).\n", CODENAME);
    printf("%s Prism integration points: %d\n",
           CODENAME, basis->num_integ_points);
}

/* -------------------------------------------------------------------------- */
/* Mixed input                                                                */
/* -------------------------------------------------------------------------- */

static int fluid_mixed_read_partitioned_elements(
        BBFE_DATA*  fe,
        const char* filename,
        const char* directory,
        int         num_integ_points)
{
    FILE* fp = NULL;
    const int rank = monolis_mpi_get_global_my_rank();

    /* Always initialize the element-family state before reading. */
    fe->total_num_elems = 0;
    fe->local_num_nodes = 0;

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    if(fp == NULL) {
        fprintf(stderr,
            "%s rank=%d ERROR: cannot open element file %s\n",
            CODENAME, rank, filename);
        return -1;
    }

    /* File format:
     *   Nelem
     *   elem_id NEN node0 node1 ...
     *   elem_id NEN node0 node1 ...
     */
    if(fscanf(fp, "%d", &(fe->total_num_elems)) != 1) {
        fprintf(stderr,
            "%s rank=%d ERROR: cannot read element count from %s\n",
            CODENAME, rank, filename);
        fclose(fp);
        return -1;
    }

    printf("%s rank=%d: %s Num. elements: %d\n",
           CODENAME, rank, filename, fe->total_num_elems);

    if(fe->total_num_elems < 0) {
        fprintf(stderr,
            "%s rank=%d ERROR: negative element count %d in %s\n",
            CODENAME, rank, fe->total_num_elems, filename);
        fclose(fp);
        return -1;
    }

    /* An MPI partition is allowed to contain none of this element family. */
    if(fe->total_num_elems == 0) {
        fe->local_num_nodes = 0;
        fclose(fp);
        printf("%s rank=%d: %s is empty -> skip this element family.\n",
               CODENAME, rank, filename);
        return 0;
    }

    int elem_id = -1;
    int nen = 0;

    /* Read NEN from the first active element before allocation. */
    if(fscanf(fp, "%d %d", &elem_id, &nen) != 2 || nen <= 0) {
        fprintf(stderr,
            "%s rank=%d ERROR: invalid first element header in %s\n",
            CODENAME, rank, filename);
        fclose(fp);
        return -1;
    }

    fe->local_num_nodes = nen;

    BBFE_sys_memory_allocation_elem(fe, num_integ_points, 3);

    for(int i = 0; i < fe->local_num_nodes; ++i) {
        if(fscanf(fp, "%d", &(fe->conn[0][i])) != 1) {
            fprintf(stderr,
                "%s rank=%d ERROR: truncated connectivity in %s at e=0 i=%d\n",
                CODENAME, rank, filename, i);
            fclose(fp);
            return -1;
        }
    }

    for(int e = 1; e < fe->total_num_elems; ++e) {
        int current_id = -1;
        int current_nen = 0;

        if(fscanf(fp, "%d %d", &current_id, &current_nen) != 2) {
            fprintf(stderr,
                "%s rank=%d ERROR: invalid element header in %s at e=%d\n",
                CODENAME, rank, filename, e);
            fclose(fp);
            return -1;
        }

        /* One partitioned element file must contain one element family only. */
        if(current_nen != fe->local_num_nodes) {
            fprintf(stderr,
                "%s rank=%d ERROR: NEN changes inside %s at e=%d: first=%d current=%d\n",
                CODENAME, rank, filename, e,
                fe->local_num_nodes, current_nen);
            fclose(fp);
            return -1;
        }

        for(int i = 0; i < fe->local_num_nodes; ++i) {
            if(fscanf(fp, "%d", &(fe->conn[e][i])) != 1) {
                fprintf(stderr,
                    "%s rank=%d ERROR: truncated connectivity in %s at e=%d i=%d\n",
                    CODENAME, rank, filename, e, i);
                fclose(fp);
                return -1;
            }
        }
    }

    fclose(fp);

    printf("%s rank=%d: loaded %s elements=%d NEN=%d\n",
           CODENAME, rank, filename,
           fe->total_num_elems, fe->local_num_nodes);

    return 0;
}

static void fluid_mixed_check_mesh_pair(
        const BBFE_DATA* fe_tet,
        const BBFE_DATA* fe_pri)
{
    if(fe_tet->total_num_nodes != fe_pri->total_num_nodes) {
        fprintf(stderr,
            "%s ERROR: tetra/prism node counts differ: tet=%d prism=%d\n",
            CODENAME, fe_tet->total_num_nodes, fe_pri->total_num_nodes);
        exit(EXIT_FAILURE);
    }

    /* Empty element families are valid after graph partitioning. */
    if(fe_tet->total_num_elems > 0 && fe_tet->local_num_nodes != 4) {
        fprintf(stderr,
            "%s ERROR: active %s must be TET4; elems=%d local_num_nodes=%d\n",
            CODENAME, FLUID_MIXED_ELEM_TET_FILE,
            fe_tet->total_num_elems, fe_tet->local_num_nodes);
        exit(EXIT_FAILURE);
    }
    if(fe_pri->total_num_elems > 0 && fe_pri->local_num_nodes != 6) {
        fprintf(stderr,
            "%s ERROR: active %s must be PRI6; elems=%d local_num_nodes=%d\n",
            CODENAME, FLUID_MIXED_ELEM_PRI_FILE,
            fe_pri->total_num_elems, fe_pri->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    /* For an empty family, both 0 and the nominal NEN are accepted. */
    if(fe_tet->total_num_elems == 0 &&
       fe_tet->local_num_nodes != 0 && fe_tet->local_num_nodes != 4) {
        fprintf(stderr,
            "%s ERROR: empty TET family has invalid local_num_nodes=%d\n",
            CODENAME, fe_tet->local_num_nodes);
        exit(EXIT_FAILURE);
    }
    if(fe_pri->total_num_elems == 0 &&
       fe_pri->local_num_nodes != 0 && fe_pri->local_num_nodes != 6) {
        fprintf(stderr,
            "%s ERROR: empty PRI family has invalid local_num_nodes=%d\n",
            CODENAME, fe_pri->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    /* Both objects read the same node.dat. */
    for(int i = 0; i < fe_tet->total_num_nodes; ++i) {
        for(int d = 0; d < 3; ++d) {
            const double a = fe_tet->x[i][d];
            const double b = fe_pri->x[i][d];
            const double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
            if(fabs(a-b) > 1.0e-12 * scale) {
                fprintf(stderr,
                    "%s ERROR: tetra/prism node.dat mismatch at node=%d dir=%d: %.17e %.17e\n",
                    CODENAME, i, d, a, b);
                exit(EXIT_FAILURE);
            }
        }
    }
}

void BBFE_fluid_pre_mixed(
        BBFE_DATA*   fe_tet,
        BBFE_BASIS*  basis_tet,
        BBFE_DATA*   fe_pri,
        BBFE_BASIS*  basis_pri,
        int          argc,
        char*        argv[],
        const char*  directory,
        int          num_integ_points_each_axis)
{
    (void)argc;
    (void)argv;

    BB_calc_void();

    const int n_axis = num_integ_points_each_axis;
    const int np_alloc = n_axis*n_axis*n_axis;
    const char* filename;

    filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_MIXED_NODE_FILE);

    BBFE_sys_read_node(fe_tet, filename, directory);
    BBFE_sys_read_node(fe_pri, filename, directory);

    filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_MIXED_ELEM_TET_FILE);
    if(fluid_mixed_read_partitioned_elements(
            fe_tet, filename, directory, np_alloc) != 0) {
        exit(EXIT_FAILURE);
    }

    filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_MIXED_ELEM_PRI_FILE);
    if(fluid_mixed_read_partitioned_elements(
            fe_pri, filename, directory, np_alloc) != 0) {
        exit(EXIT_FAILURE);
    }

    fluid_mixed_check_mesh_pair(fe_tet, fe_pri);

    if(BBFE_fluid_mixed_has_elements(fe_tet)) {
        BBFE_sys_memory_allocation_integ(basis_tet, np_alloc, 3);
        BBFE_sys_memory_allocation_shapefunc(
            basis_tet, fe_tet->local_num_nodes, 1, np_alloc);
        BBFE_fluid_set_basis(
            basis_tet, fe_tet->local_num_nodes, n_axis);
    } else {
        memset(basis_tet, 0, sizeof(*basis_tet));
        printf("%s TET4: no local elements on this rank -> skip basis/geometry assembly.\n",
               CODENAME);
    }

    if(BBFE_fluid_mixed_has_elements(fe_pri)) {
        BBFE_sys_memory_allocation_integ(basis_pri, np_alloc, 3);
        BBFE_sys_memory_allocation_shapefunc(
            basis_pri, fe_pri->local_num_nodes, 1, np_alloc);
        BBFE_fluid_set_basis_prism(basis_pri, n_axis);
    } else {
        memset(basis_pri, 0, sizeof(*basis_pri));
        printf("%s PRI6: no local elements on this rank -> skip basis/geometry assembly.\n",
               CODENAME);
    }

    printf("%s Mixed volume mesh: nodes=%d tetra=%d prism=%d\n",
        CODENAME,
        fe_tet->total_num_nodes,
        fe_tet->total_num_elems,
        fe_pri->total_num_elems);
}


/* -------------------------------------------------------------------------- */
/* Robust partitioned TRI3 surface input                                      */
/* -------------------------------------------------------------------------- */

void BBFE_fluid_pre_surface_mixed(
        BBFE_DATA*   surf,
        BBFE_BASIS*  basis,
        const char*  directory,
        const char*  filename_base,
        int          num_integ_points_each_axis)
{
    if(surf == NULL || basis == NULL || directory == NULL ||
       filename_base == NULL || num_integ_points_each_axis <= 0) {
        fprintf(stderr, "%s ERROR: invalid mixed surface input.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    memset(surf, 0, sizeof(*surf));
    memset(basis, 0, sizeof(*basis));

    const int n_axis = num_integ_points_each_axis;
    const int np_elem = n_axis*n_axis;
    const int np_alloc = n_axis*n_axis*n_axis;

    const char* filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        filename_base);

    if(fluid_mixed_read_partitioned_elements(
            surf, filename, directory, np_elem) != 0) {
        fprintf(stderr, "%s ERROR: failed to read surface file %s.\n",
                CODENAME, filename);
        exit(EXIT_FAILURE);
    }

    if(surf->total_num_elems == 0) {
        surf->local_num_nodes = 0;
        printf("%s rank=%d: %s is empty -> skip TRI3 basis.\n",
               CODENAME, monolis_mpi_get_global_my_rank(), filename_base);
        return;
    }

    if(surf->local_num_nodes != 3) {
        fprintf(stderr,
            "%s ERROR: active %s must contain TRI3; NEN=%d.\n",
            CODENAME, filename_base, surf->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    BBFE_sys_memory_allocation_integ(basis, np_alloc, 2);
    BBFE_sys_memory_allocation_shapefunc(
        basis, surf->local_num_nodes, 1, np_alloc);

    BBFE_convdiff_set_basis_surface(
        basis, surf->local_num_nodes, n_axis);
}

void BBFE_fluid_finalize_surface_mixed(
        BBFE_DATA*  surf,
        BBFE_BASIS* basis)
{
    if(surf == NULL || basis == NULL) return;

    if(surf->total_num_elems > 0 && surf->local_num_nodes > 0) {
        const int np = basis->num_integ_points;
        BBFE_sys_memory_free_integ(basis, 2);
        BBFE_sys_memory_free_shapefunc(basis);
        BBFE_sys_memory_free_elem(surf, np, 3);
    }

    memset(surf, 0, sizeof(*surf));
    memset(basis, 0, sizeof(*basis));
}

/* -------------------------------------------------------------------------- */
/* Unified nodal graph / Monolis matrix                                       */
/* -------------------------------------------------------------------------- */

static void fluid_mixed_read_nodal_graph(
        const char* directory,
        const char* filename,
        int*        num_nodes_out,
        int**       index_out,
        int**       item_out)
{
    FILE* fp = NULL;
    int nnode = 0;
    long long nitem_ll = 0;

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    if(fp == NULL || fscanf(fp, "%d", &nnode) != 1 || nnode <= 0) {
        fprintf(stderr, "%s ERROR: cannot read graph header: %s\n",
                CODENAME, filename);
        if(fp) fclose(fp);
        exit(EXIT_FAILURE);
    }

    for(int row = 0; row < nnode; ++row) {
        int node_id = -1;
        int degree = -1;
        if(fscanf(fp, "%d %d", &node_id, &degree) != 2 || degree < 0) {
            fprintf(stderr, "%s ERROR: invalid graph row %d in %s\n",
                    CODENAME, row, filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        if(node_id != row) {
            fprintf(stderr,
                "%s ERROR: graph rows must be ordered 0..N-1; row=%d id=%d\n",
                CODENAME, row, node_id);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        for(int j = 0; j < degree; ++j) {
            int dummy = -1;
            if(fscanf(fp, "%d", &dummy) != 1) {
                fprintf(stderr, "%s ERROR: truncated graph row %d in %s\n",
                        CODENAME, row, filename);
                fclose(fp);
                exit(EXIT_FAILURE);
            }
        }
        nitem_ll += degree;
    }
    fclose(fp);

    if(nitem_ll > 2147483647LL) {
        fprintf(stderr, "%s ERROR: graph adjacency exceeds 32-bit int.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    const int nitem = (int)nitem_ll;
    int* index = (int*)calloc((size_t)nnode + 1u, sizeof(int));
    int* item  = (int*)calloc((size_t)(nitem > 0 ? nitem : 1), sizeof(int));
    if(index == NULL || item == NULL) {
        fprintf(stderr, "%s ERROR: graph allocation failed.\n", CODENAME);
        free(index);
        free(item);
        exit(EXIT_FAILURE);
    }

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    if(fp == NULL || fscanf(fp, "%d", &nnode) != 1) {
        fprintf(stderr, "%s ERROR: cannot reopen graph: %s\n", CODENAME, filename);
        free(index); free(item);
        if(fp) fclose(fp);
        exit(EXIT_FAILURE);
    }

    int cursor = 0;
    index[0] = 0;
    for(int row = 0; row < nnode; ++row) {
        int node_id = -1;
        int degree = -1;
        if(fscanf(fp, "%d %d", &node_id, &degree) != 2) {
            fprintf(stderr, "%s ERROR: graph second-pass read failed row=%d\n",
                    CODENAME, row);
            fclose(fp); free(index); free(item); exit(EXIT_FAILURE);
        }
        for(int j = 0; j < degree; ++j) {
            int adj = -1;
            if(fscanf(fp, "%d", &adj) != 1 || adj < 0 || adj >= nnode) {
                fprintf(stderr,
                    "%s ERROR: invalid graph adjacency row=%d value=%d\n",
                    CODENAME, row, adj);
                fclose(fp); free(index); free(item); exit(EXIT_FAILURE);
            }
            item[cursor++] = adj;
        }
        index[row+1] = cursor;
    }
    fclose(fp);

    *num_nodes_out = nnode;
    *index_out = index;
    *item_out = item;
}

void BBFE_fluid_init_monomat_mixed(
        MONOLIS*     monolis,
        MONOLIS_COM* monolis_com,
        BBFE_DATA*   fe_ref,
        const char*  directory)
{
    monolis_initialize(monolis);

    /*
     * The communication table (.recv/.send) is generated from graph.dat.
     * node.dat contains coordinates only and must not be used as the
     * communication-table base name.
     */
    monolis_com_initialize_by_parted_files(
        monolis_com,
        monolis_mpi_get_global_comm(),
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_MIXED_GRAPH_FILE);

    const char* filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        FLUID_MIXED_GRAPH_FILE);

    int nnode = 0;
    int* index = NULL;
    int* item = NULL;
    fluid_mixed_read_nodal_graph(
        directory, filename, &nnode, &index, &item);

    if(nnode != fe_ref->total_num_nodes) {
        fprintf(stderr,
            "%s ERROR: graph/node mismatch: graph=%d node.dat=%d file=%s\n",
            CODENAME, nnode, fe_ref->total_num_nodes, filename);
        free(index); free(item);
        exit(EXIT_FAILURE);
    }

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        monolis,
        nnode,
        FLUID_MIXED_BLOCK_SIZE,
        index,
        item);

    printf("%s Unified TET4+PRI6 graph: nodes=%d directed_adjacencies=%d block=%d\n",
        CODENAME, nnode, index[nnode], FLUID_MIXED_BLOCK_SIZE);

    free(index);
    free(item);
}


/* -------------------------------------------------------------------------- */
/* MPI overlap synchronization / PRI6 diagnostics                             */
/* -------------------------------------------------------------------------- */

void BBFE_fluid_mixed_sync_solution_vector(
        MONOLIS_COM* mono_com,
        int          total_num_nodes,
        double*      vec4)
{
    if(mono_com == NULL || vec4 == NULL || total_num_nodes <= 0) return;

    /* Owner/internal node values overwrite overlap copies. */
    monolis_mpi_update_R(
        mono_com,
        total_num_nodes,
        FLUID_MIXED_BLOCK_SIZE,
        vec4);
}

void BBFE_fluid_mixed_sync_state(
        MONOLIS_COM* mono_com,
        VALUES*      vals,
        int          total_num_nodes)
{
    if(mono_com == NULL || vals == NULL || total_num_nodes <= 0) return;

    double* work = BB_std_calloc_1d_double(NULL, 4*total_num_nodes);
    if(work == NULL) {
        fprintf(stderr, "%s ERROR: state-sync allocation failed.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    BBFE_fluid_sups_add_velocity_pressure(
        vals->v,
        vals->p,
        work,
        total_num_nodes);

    monolis_mpi_update_R(
        mono_com,
        total_num_nodes,
        FLUID_MIXED_BLOCK_SIZE,
        work);

    BBFE_fluid_sups_renew_velocity(
        vals->v,
        work,
        total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        vals->p,
        work,
        total_num_nodes);

    BB_std_free_1d_double(work, 4*total_num_nodes);
}

void BBFE_fluid_mixed_report_prism_bc(
        const FE_SYSTEM_FLUID_MIXED* sys)
{
    if(sys == NULL) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    unsigned char* used = (unsigned char*)calloc(
        (size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));
    if(used == NULL) return;

    int local_nodes = 0;
    int local_bc[4] = {0,0,0,0};

    if(BBFE_fluid_mixed_has_elements(&(sys->fe_pri))) {
        for(int e = 0; e < sys->fe_pri.total_num_elems; ++e) {
            for(int a = 0; a < sys->fe_pri.local_num_nodes; ++a) {
                const int i = sys->fe_pri.conn[e][a];
                if(i >= 0 && i < nnode) used[i] = 1;
            }
        }

        for(int i = 0; i < nnode; ++i) {
            if(!used[i]) continue;
            ++local_nodes;
            for(int d = 0; d < 4; ++d) {
                if(sys->bc_NR.D_bc_exists != NULL &&
                   sys->bc_NR.D_bc_exists[4*i+d]) {
                    ++local_bc[d];
                }
            }
        }
    }

    int global_nodes = local_nodes;
    int global_bc[4] = {local_bc[0],local_bc[1],local_bc[2],local_bc[3]};

    monolis_allreduce_I(1, &global_nodes, MONOLIS_MPI_SUM, sys->mono_com.comm);
    for(int d = 0; d < 4; ++d) {
        monolis_allreduce_I(1, &global_bc[d], MONOLIS_MPI_SUM, sys->mono_com.comm);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "[PRI6 BC] nodes=%d constrained ux=%d uy=%d uz=%d p=%d\n",
            global_nodes,
            global_bc[0], global_bc[1], global_bc[2], global_bc[3]);
    }

    free(used);
}

void BBFE_fluid_mixed_report_prism_components(
        const FE_SYSTEM_FLUID_MIXED* sys,
        int step,
        int newton_iteration)
{
    if(sys == NULL) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    unsigned char* used = (unsigned char*)calloc(
        (size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));
    if(used == NULL) return;

    double local_x_l1[4]  = {0,0,0,0};
    double local_x_max[4] = {0,0,0,0};
    double local_v_l1[4]  = {0,0,0,0};
    double local_v_max[4] = {0,0,0,0};
    int local_internal_nodes = 0;

    if(BBFE_fluid_mixed_has_elements(&(sys->fe_pri))) {
        for(int e = 0; e < sys->fe_pri.total_num_elems; ++e) {
            for(int a = 0; a < sys->fe_pri.local_num_nodes; ++a) {
                const int i = sys->fe_pri.conn[e][a];
                if(i >= 0 && i < nnode) used[i] = 1;
            }
        }

        /* Count owner/internal values only to avoid overlap double counting. */
        const int nint = sys->mono_com.n_internal_vertex;
        for(int i = 0; i < nint; ++i) {
            if(!used[i]) continue;
            ++local_internal_nodes;

            const double dx[4] = {
                sys->monolis.mat.R.X[4*i+0],
                sys->monolis.mat.R.X[4*i+1],
                sys->monolis.mat.R.X[4*i+2],
                sys->monolis.mat.R.X[4*i+3]
            };
            const double vv[4] = {
                sys->vals.v[i][0],
                sys->vals.v[i][1],
                sys->vals.v[i][2],
                sys->vals.p[i]
            };

            for(int d = 0; d < 4; ++d) {
                const double ax = fabs(dx[d]);
                const double av = fabs(vv[d]);
                local_x_l1[d] += ax;
                local_v_l1[d] += av;
                if(ax > local_x_max[d]) local_x_max[d] = ax;
                if(av > local_v_max[d]) local_v_max[d] = av;
            }
        }
    }

    int global_nodes = local_internal_nodes;
    double gx_l1[4], gx_max[4], gv_l1[4], gv_max[4];
    for(int d=0; d<4; ++d) {
        gx_l1[d]=local_x_l1[d]; gx_max[d]=local_x_max[d];
        gv_l1[d]=local_v_l1[d]; gv_max[d]=local_v_max[d];
    }

    monolis_allreduce_I(1, &global_nodes, MONOLIS_MPI_SUM, sys->mono_com.comm);
    for(int d=0; d<4; ++d) {
        monolis_allreduce_R(1, &gx_l1[d],  MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &gx_max[d], MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &gv_l1[d],  MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &gv_max[d], MONOLIS_MPI_MAX, sys->mono_com.comm);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "[PRI6 components step=%d NR=%d] internal_nodes=%d\n"
            "  dX ux: L1=%e max=%e   state ux: L1=%e max=%e\n"
            "  dX uy: L1=%e max=%e   state uy: L1=%e max=%e\n"
            "  dX uz: L1=%e max=%e   state uz: L1=%e max=%e\n"
            "  dX  p: L1=%e max=%e   state  p: L1=%e max=%e\n",
            step, newton_iteration, global_nodes,
            gx_l1[0],gx_max[0],gv_l1[0],gv_max[0],
            gx_l1[1],gx_max[1],gv_l1[1],gv_max[1],
            gx_l1[2],gx_max[2],gv_l1[2],gv_max[2],
            gx_l1[3],gx_max[3],gv_l1[3],gv_max[3]);
    }

    free(used);
}

void BBFE_fluid_mixed_report_prism_partition(
        const BBFE_DATA*   fe_pri,
        const MONOLIS_COM* mono_com)
{
    if(fe_pri == NULL || mono_com == NULL) return;

    const int nnode = fe_pri->total_num_nodes;
    unsigned char* used = (unsigned char*)calloc(
        (size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));
    if(used == NULL) return;

    int local_unique = 0;
    int local_internal = 0;
    int local_overlap = 0;

    if(BBFE_fluid_mixed_has_elements(fe_pri)) {
        for(int e = 0; e < fe_pri->total_num_elems; ++e) {
            for(int a = 0; a < fe_pri->local_num_nodes; ++a) {
                const int i = fe_pri->conn[e][a];
                if(i >= 0 && i < nnode) used[i] = 1;
            }
        }

        for(int i = 0; i < nnode; ++i) {
            if(!used[i]) continue;
            ++local_unique;
            if(i < mono_com->n_internal_vertex) ++local_internal;
            else ++local_overlap;
        }
    }

    printf(
        "%s rank=%d PRI6 partition: elems=%d unique_nodes=%d internal=%d overlap=%d\n",
        CODENAME,
        monolis_mpi_get_global_my_rank(),
        fe_pri->total_num_elems,
        local_unique,
        local_internal,
        local_overlap);

    free(used);
}

void BBFE_fluid_mixed_report_prism_update(
        const FE_SYSTEM_FLUID_MIXED* sys,
        int                          step,
        int                          newton_iteration)
{
    if(sys == NULL) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    unsigned char* used = (unsigned char*)calloc(
        (size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));
    if(used == NULL) return;

    int local_nodes = 0;
    int local_internal = 0;
    int local_overlap = 0;
    int local_changed = 0;
    double local_max_dv = 0.0;
    double local_max_dp = 0.0;

    if(BBFE_fluid_mixed_has_elements(&(sys->fe_pri))) {
        for(int e = 0; e < sys->fe_pri.total_num_elems; ++e) {
            for(int a = 0; a < sys->fe_pri.local_num_nodes; ++a) {
                const int i = sys->fe_pri.conn[e][a];
                if(i >= 0 && i < nnode) used[i] = 1;
            }
        }

        for(int i = 0; i < nnode; ++i) {
            if(!used[i]) continue;

            ++local_nodes;
            if(i < sys->mono_com.n_internal_vertex) ++local_internal;
            else ++local_overlap;

            double dv = 0.0;
            for(int d = 0; d < 3; ++d) {
                const double a = fabs(sys->vals.delta_v[i][d]);
                if(a > dv) dv = a;
            }

            const double dp = fabs(sys->vals.delta_p[i]);
            if(dv > local_max_dv) local_max_dv = dv;
            if(dp > local_max_dp) local_max_dp = dp;
            if(dv > 1.0e-14 || dp > 1.0e-14) ++local_changed;
        }
    }

    free(used);

    int global_nodes = local_nodes;
    int global_internal = local_internal;
    int global_overlap = local_overlap;
    int global_changed = local_changed;
    double global_max_dv = local_max_dv;
    double global_max_dp = local_max_dp;

    monolis_allreduce_I(1, &global_nodes,    MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_I(1, &global_internal, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_I(1, &global_overlap,  MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_I(1, &global_changed,  MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_R(1, &global_max_dv,   MONOLIS_MPI_MAX, sys->mono_com.comm);
    monolis_allreduce_R(1, &global_max_dp,   MONOLIS_MPI_MAX, sys->mono_com.comm);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "[PRI6 update step=%d NR=%d] nodes=%d internal=%d overlap=%d "
            "changed=%d max|delta_v|=%e max|delta_p|=%e\n",
            step,
            newton_iteration,
            global_nodes,
            global_internal,
            global_overlap,
            global_changed,
            global_max_dv,
            global_max_dp);
    }
}

void BBFE_fluid_mixed_check_positive_jacobian(
        const BBFE_DATA*  fe,
        const BBFE_BASIS* basis,
        const char*       element_name,
        double            tolerance)
{
    int local_bad = 0;
    double local_min = DBL_MAX;

    if(fe != NULL && basis != NULL && BBFE_fluid_mixed_has_elements(fe)) {
        for(int e = 0; e < fe->total_num_elems; ++e) {
            for(int p = 0; p < basis->num_integ_points; ++p) {
                const double J = fe->geo[e][p].Jacobian;
                if(J < local_min) local_min = J;
                if(!isfinite(J) || J <= tolerance) ++local_bad;
            }
        }
    }

    int global_bad = local_bad;
    monolis_allreduce_I(
        1, &global_bad, MONOLIS_MPI_SUM, monolis_mpi_get_global_comm());

    double neg_min = (local_min < DBL_MAX) ? -local_min : -DBL_MAX;
    monolis_allreduce_R(
        1, &neg_min, MONOLIS_MPI_MAX, monolis_mpi_get_global_comm());

    if(monolis_mpi_get_global_my_rank() == 0) {
        const double global_min = (neg_min > -DBL_MAX) ? -neg_min : 0.0;
        printf("%s %s Jacobian: min=%e bad_qp=%d\n",
               CODENAME, element_name, global_min, global_bad);
    }

    if(global_bad > 0) {
        fprintf(stderr, "%s ERROR: non-positive Jacobian in %s.\n",
                CODENAME, element_name);
        exit(EXIT_FAILURE);
    }
}

/* -------------------------------------------------------------------------- */
/* Mixed SUPS / Newton assembly                                               */
/* -------------------------------------------------------------------------- */

void set_element_mat_NR_Tezuer_mixed(
        MONOLIS*    monolis,
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri,
        VALUES*     vals)
{
    /* Each MPI rank may contain only one of the two element families. */
    if(BBFE_fluid_mixed_has_elements(fe_tet)) {
        set_element_mat_NR_Tezuer(monolis, fe_tet, basis_tet, vals);
    }
    if(BBFE_fluid_mixed_has_elements(fe_pri)) {
        set_element_mat_NR_Tezuer(monolis, fe_pri, basis_pri, vals);
    }
}

void solver_fom_NR_mixed(
        FE_SYSTEM_FLUID_MIXED* sys,
        double                 t,
        int                    step)
{
    (void)t;

    printf("\n%s ----------------- mixed Time step %d ----------------\n",
           CODENAME, step);

    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double tiny = 1.0e-30;

    /* Preserve the current FOM behavior.  Increase this if true Newton
       iterations are desired. */
    const int max_iter_NR = 1;

    const int nnode = sys->fe_tet.total_num_nodes;

    for(int it = 0; it < max_iter_NR; ++it) {
        printf("\n%s ---- step %d : NR %d ----\n", CODENAME, step, it);

        monolis_clear_mat_value_R(&(sys->monolis));

        set_element_mat_NR_Tezuer_mixed(
            &(sys->monolis),
            &(sys->fe_tet), &(sys->basis_tet),
            &(sys->fe_pri), &(sys->basis_pri),
            &(sys->vals));

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            FLUID_MIXED_BLOCK_SIZE,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        BBFE_sys_monowrap_solve(
            &(sys->monolis),
            &(sys->mono_com),
            sys->monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        BBFE_fluid_sups_renew_velocity(
            sys->vals.delta_v,
            sys->monolis.mat.R.X,
            nnode);

        BBFE_fluid_sups_renew_pressure(
            sys->vals.delta_p,
            sys->monolis.mat.R.X,
            nnode);

        update_velocity_pressure_NR(
            sys->vals.v,
            sys->vals.delta_v,
            sys->vals.p,
            sys->vals.delta_p,
            nnode);

        double norm_v = calc_internal_norm_2d(
            sys->vals.v,
            sys->mono_com.n_internal_vertex,
            3);
        double norm_delta_v = calc_internal_norm_2d(
            sys->vals.delta_v,
            sys->mono_com.n_internal_vertex,
            3);

        double norm_p = 0.0;
        double norm_delta_p = 0.0;
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;

        for(int i = 0; i < sys->mono_com.n_internal_vertex; ++i) {
            const double pv = sys->vals.p[i];
            const double dp = sys->vals.delta_p[i];
            norm_p += pv*pv;
            norm_delta_p += dp*dp;

            for(int d = 0; d < 3; ++d) {
                const double a = fabs(sys->vals.delta_v[i][d]);
                if(a > linf_delta_v_local) linf_delta_v_local = a;
            }
            if(fabs(dp) > linf_delta_p_local) linf_delta_p_local = fabs(dp);
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys->mono_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys->mono_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys->mono_com.comm);

        const double nrm_v  = sqrt(norm_v);
        const double nrm_dv = sqrt(norm_delta_v);
        const double nrm_p  = sqrt(norm_p);
        const double nrm_dp = sqrt(norm_delta_p);
        const double denom_v = fmax(nrm_v, tiny);
        const double denom_p = fmax(nrm_p, tiny);

        const int conv_v =
            (nrm_dv <= abs_tol_v) ||
            (nrm_dv/denom_v <= rel_tol_v) ||
            (linf_delta_v <= abs_tol_v);
        const int conv_p =
            (nrm_dp <= abs_tol_p) ||
            (nrm_dp/denom_p <= rel_tol_p) ||
            (linf_delta_p <= abs_tol_p);

        printf("[mixed NR %2d] ||dv||2=%.3e ||v||2=%.3e rel=%.3e Linf(dv)=%.3e\n",
            it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
        printf("[mixed NR %2d] ||dp||2=%.3e ||p||2=%.3e rel=%.3e Linf(dp)=%.3e\n",
            it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p);

        if(conv_v && conv_p) {
            double max_du = 0.0;
            for(int i = 0; i < nnode; ++i) {
                for(int d = 0; d < 3; ++d) {
                    const double du = fabs(sys->vals.v[i][d] - sys->vals.v_old[i][d]);
                    if(du > max_du) max_du = du;
                }
            }
            printf("[mixed step %d] max|v^{n+1}-v^{n}| = %.6e\n",
                   step, max_du);

            ROM_BB_vec_copy_2d(
                sys->vals.v,
                sys->vals.v_old,
                nnode,
                3);
            break;
        }
    }
}

/* -------------------------------------------------------------------------- */
/* Mixed VTK output                                                           */
/* -------------------------------------------------------------------------- */

void output_result_file_vtk_mixed(
        const BBFE_DATA* fe_tet,
        const BBFE_DATA* fe_pri,
        const VALUES*    vals,
        const char*      filename,
        const char*      directory,
        double           t)
{
    (void)t;

    FILE* fp = NULL;
    fp = BBFE_sys_write_fopen(fp, filename, directory);
    if(fp == NULL) {
        fprintf(stderr, "%s ERROR: cannot open mixed VTK output: %s\n",
                CODENAME, filename);
        return;
    }

    const int nnode = fe_tet->total_num_nodes;
    const int ntet  = fe_tet->total_num_elems;
    const int npri  = fe_pri->total_num_elems;
    const int ncells = ntet + npri;
    const long long cell_list_size =
        (long long)ntet * 5LL + (long long)npri * 7LL;

    unsigned char* is_prism_node = (unsigned char*)calloc(
        (size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));

    if(is_prism_node != NULL) {
        for(int e = 0; e < npri; ++e) {
            for(int a = 0; a < 6; ++a) {
                const int i = fe_pri->conn[e][a];
                if(i >= 0 && i < nnode) is_prism_node[i] = 1;
            }
        }
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "TET4+PRI6 fluid result\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d double\n", nnode);
    for(int i = 0; i < nnode; ++i) {
        fprintf(fp, "%.15e %.15e %.15e\n",
            fe_tet->x[i][0], fe_tet->x[i][1], fe_tet->x[i][2]);
    }

    fprintf(fp, "CELLS %d %lld\n", ncells, cell_list_size);
    for(int e = 0; e < ntet; ++e) {
        fprintf(fp, "4 %d %d %d %d\n",
            fe_tet->conn[e][0], fe_tet->conn[e][1],
            fe_tet->conn[e][2], fe_tet->conn[e][3]);
    }
    for(int e = 0; e < npri; ++e) {
        fprintf(fp, "6 %d %d %d %d %d %d\n",
            fe_pri->conn[e][0], fe_pri->conn[e][1], fe_pri->conn[e][2],
            fe_pri->conn[e][3], fe_pri->conn[e][4], fe_pri->conn[e][5]);
    }

    fprintf(fp, "CELL_TYPES %d\n", ncells);
    for(int e = 0; e < ntet; ++e) fprintf(fp, "10\n"); /* VTK_TETRA */
    for(int e = 0; e < npri; ++e) fprintf(fp, "13\n"); /* VTK_WEDGE */

    fprintf(fp, "POINT_DATA %d\n", nnode);
    fprintf(fp, "VECTORS Velocity double\n");
    for(int i = 0; i < nnode; ++i) {
        fprintf(fp, "%.15e %.15e %.15e\n",
            vals->v[i][0], vals->v[i][1], vals->v[i][2]);
    }

    fprintf(fp, "SCALARS Pressure double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(int i = 0; i < nnode; ++i) {
        fprintf(fp, "%.15e\n", vals->p[i]);
    }

    fprintf(fp, "SCALARS IsPrismNode int 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(int i = 0; i < nnode; ++i) {
        fprintf(fp, "%d\n", is_prism_node ? (int)is_prism_node[i] : 0);
    }

    if(vals->y_plus != NULL) {
        fprintf(fp, "SCALARS yPlus double 1\n");
        fprintf(fp, "LOOKUP_TABLE default\n");
        for(int i = 0; i < nnode; ++i) {
            fprintf(fp, "%.15e\n", vals->y_plus[i]);
        }
    }

    fprintf(fp, "CELL_DATA %d\n", ncells);
    fprintf(fp, "SCALARS ElementFamily int 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(int e = 0; e < ntet; ++e) fprintf(fp, "4\n");
    for(int e = 0; e < npri; ++e) fprintf(fp, "6\n");

    free(is_prism_node);
    fclose(fp);
}

void output_files_mixed(
        FE_SYSTEM_FLUID_MIXED* sys,
        int                    file_num,
        double                 t)
{
    /* Keep overlap-node values consistent before writing partition VTK files. */
    BBFE_fluid_mixed_sync_state(
        &(sys->mono_com), &(sys->vals), sys->fe_tet.total_num_nodes);

    const int myrank = monolis_mpi_get_global_my_rank();
    char local_name[FLUID_MIXED_IO_BUFFER];
    snprintf(local_name, sizeof(local_name),
             "result_mixed_%06d.%d.vtk", file_num, myrank);

    const char* filename = monolis_get_global_output_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        "./",
        local_name);

    output_result_file_vtk_mixed(
        &(sys->fe_tet),
        &(sys->fe_pri),
        &(sys->vals),
        filename,
        sys->cond.directory,
        t);
}

void BBFE_fluid_finalize_mixed(
        BBFE_DATA*  fe_tet,
        BBFE_BASIS* basis_tet,
        BBFE_DATA*  fe_pri,
        BBFE_BASIS* basis_pri)
{
    /* basis arrays are allocated only for active element families. */
    if(BBFE_fluid_mixed_has_elements(fe_tet)) {
        BBFE_fluid_finalize(fe_tet, basis_tet);
    } else {
        /* node.dat is still read for the empty family. */
        BBFE_sys_memory_free_node(fe_tet, 3);
    }

    if(BBFE_fluid_mixed_has_elements(fe_pri)) {
        BBFE_fluid_finalize(fe_pri, basis_pri);
    } else {
        BBFE_sys_memory_free_node(fe_pri, 3);
    }
}
