#include "core_FOM_ST_spaceblock.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int st_block(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    int alpha)
{
    return ST_fom_block_dof(layout, slab, alpha, 0);
}

static int st_resolve_owned_nodes(ST_DIFFUSION_SPACEBLOCK* solver)
{
    const int size = monolis_mpi_get_global_comm_size();
    const int local_nodes = solver->fom->fe.total_num_nodes;
    const int comm_owned = solver->fom->monolis_com.n_internal_vertex;

    return ST_fom_resolve_owned_space_points(
        local_nodes,
        size,
        comm_owned,
        &solver->num_owned_space_nodes);
}

static int st_solve_mode_from_environment(void)
{
    const char* value = getenv("ST_FOM_SOLVE_MODE");

    if(value == NULL || value[0] == '\0' || strcmp(value, "causal") == 0){
        return ST_FOM_SOLVE_CAUSAL;
    }
    if(strcmp(value, "window") == 0 || strcmp(value, "all-at-once") == 0){
        return ST_FOM_SOLVE_WINDOW;
    }

    fprintf(
        stderr,
        "WARNING: unknown ST_FOM_SOLVE_MODE='%s'; using causal.\n",
        value);
    return ST_FOM_SOLVE_CAUSAL;
}


static int st_env_flag(const char* name, int default_value)
{
    const char* value = getenv(name);
    if(value == NULL || value[0] == '\0'){
        return default_value;
    }
    if(strcmp(value, "0") == 0 || strcmp(value, "false") == 0 ||
       strcmp(value, "off") == 0 || strcmp(value, "no") == 0)
    {
        return 0;
    }
    return 1;
}

static int st_env_positive_int(const char* name, int default_value)
{
    const char* value = getenv(name);
    char* end = NULL;
    long parsed;

    if(value == NULL || value[0] == '\0'){
        return default_value;
    }

    parsed = strtol(value, &end, 10);
    if(end == value || *end != '\0' || parsed <= 0 || parsed > INT_MAX){
        fprintf(
            stderr,
            "WARNING: invalid %s='%s'; using %d.\n",
            name,
            value,
            default_value);
        return default_value;
    }
    return (int)parsed;
}

static void st_vtk_write_scalar(
    FILE* fp,
    const char* name,
    const double* value,
    int n)
{
    fprintf(fp, "SCALARS %s double 1\n", name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(int i = 0; i < n; i++){
        fprintf(fp, "%.17e\n", value[i]);
    }
}

static int st_open_vtk_collection(ST_DIFFUSION_SPACEBLOCK* solver)
{
    const int rank = monolis_mpi_get_global_my_rank();
    const char* filename = getenv("ST_FOM_VTK_PVD");

    solver->write_vtk = st_env_flag("ST_FOM_WRITE_VTK", 1);
    solver->vtk_interval = st_env_positive_int(
        "ST_FOM_VTK_INTERVAL",
        solver->fom->vals.output_interval > 0
            ? solver->fom->vals.output_interval : 1);

    if(!solver->write_vtk){
        return ST_FOM_SUCCESS;
    }

    if(filename == NULL || filename[0] == '\0'){
        filename = "st_results.pvd";
    }
    snprintf(
        solver->vtk_pvd_filename,
        sizeof(solver->vtk_pvd_filename),
        "%s",
        filename);

    if(rank == 0){
        solver->vtk_pvd_fp = fopen(solver->vtk_pvd_filename, "w");
        if(solver->vtk_pvd_fp == NULL){
            fprintf(
                stderr,
                "ERROR: cannot open VTK collection '%s'.\n",
                solver->vtk_pvd_filename);
            return ST_FOM_ERR_IO;
        }
        fprintf(solver->vtk_pvd_fp, "<?xml version=\"1.0\"?>\n");
        fprintf(
            solver->vtk_pvd_fp,
            "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(solver->vtk_pvd_fp, "  <Collection>\n");
        fflush(solver->vtk_pvd_fp);
    }

    if(rank == 0){
        printf(
            "ST-FOM VTK output enabled: interval=%d collection=%s\n",
            solver->vtk_interval,
            solver->vtk_pvd_filename);
        fflush(stdout);
    }
    return ST_FOM_SUCCESS;
}

static void st_close_vtk_collection(ST_DIFFUSION_SPACEBLOCK* solver)
{
    if(solver == NULL || solver->vtk_pvd_fp == NULL){
        return;
    }
    fprintf(solver->vtk_pvd_fp, "  </Collection>\n");
    fprintf(solver->vtk_pvd_fp, "</VTKFile>\n");
    fclose(solver->vtk_pvd_fp);
    solver->vtk_pvd_fp = NULL;
}

static int st_write_vtk_piece(
    ST_DIFFUSION_SPACEBLOCK* solver,
    const double* numerical,
    int global_step,
    double time)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    const int rank = monolis_mpi_get_global_my_rank();
    const int nnode = fe->total_num_nodes;
    const int nelem = fe->total_num_elems;
    const int nl = fe->local_num_nodes;
    int vtk_cell_type;
    char filename[4096];
    FILE* fp = NULL;
    double* exact = NULL;
    double* abs_error = NULL;
    double* source = NULL;
    double* dirichlet = NULL;
    double* rank_value = NULL;

    if(numerical == NULL || nnode <= 0 || nelem <= 0){
        return ST_FOM_ERR_ARGUMENT;
    }

    if(nl == 8){
        vtk_cell_type = 12; /* VTK_HEXAHEDRON */
    }
    else if(nl == 4){
        vtk_cell_type = 10; /* VTK_TETRA */
    }
    else{
        fprintf(
            stderr,
            "ERROR [rank %d]: unsupported VTK element with %d nodes.\n",
            rank,
            nl);
        return ST_FOM_ERR_UNSUPPORTED;
    }

    exact = (double*)calloc((size_t)nnode, sizeof(double));
    abs_error = (double*)calloc((size_t)nnode, sizeof(double));
    source = (double*)calloc((size_t)nnode, sizeof(double));
    dirichlet = (double*)calloc((size_t)nnode, sizeof(double));
    rank_value = (double*)calloc((size_t)nnode, sizeof(double));
    if(exact == NULL || abs_error == NULL || source == NULL ||
       dirichlet == NULL || rank_value == NULL)
    {
        free(exact); free(abs_error); free(source); free(dirichlet); free(rank_value);
        return ST_FOM_ERR_ALLOCATION;
    }

    for(int node = 0; node < nnode; node++){
        double a;
        double k;
        double v[3];

        exact[node] = manusol_get_sol(
            fe->x[node][0], fe->x[node][1], fe->x[node][2], time);
        abs_error[node] = fabs(numerical[node] - exact[node]);
        a = manusol_get_mass_coef(fe->x[node]);
        k = manusol_get_diff_coef(fe->x[node]);
        manusol_get_conv_vel(v, fe->x[node]);
        source[node] = manusol_get_source(fe->x[node], time, a, v, k);
        dirichlet[node] = fom->bc.D_bc_exists[node] ? 1.0 : 0.0;
        rank_value[node] = (double)rank;
    }

    snprintf(
        filename,
        sizeof(filename),
        "st_result_%06d_rank%04d.vtk",
        global_step,
        rank);
    fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "ERROR [rank %d]: cannot open '%s'.\n", rank, filename);
        free(exact); free(abs_error); free(source); free(dirichlet); free(rank_value);
        return ST_FOM_ERR_IO;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "ST-FOM diffusion step %d time %.17g rank %d\n", global_step, time, rank);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", nnode);
    for(int node = 0; node < nnode; node++){
        fprintf(
            fp,
            "%.17e %.17e %.17e\n",
            fe->x[node][0],
            fe->x[node][1],
            fe->x[node][2]);
    }

    fprintf(fp, "CELLS %d %d\n", nelem, nelem * (nl + 1));
    for(int e = 0; e < nelem; e++){
        fprintf(fp, "%d", nl);
        for(int i = 0; i < nl; i++){
            fprintf(fp, " %d", fe->conn[e][i]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "CELL_TYPES %d\n", nelem);
    for(int e = 0; e < nelem; e++){
        fprintf(fp, "%d\n", vtk_cell_type);
    }

    fprintf(fp, "POINT_DATA %d\n", nnode);
    st_vtk_write_scalar(fp, "temperature", numerical, nnode);
    st_vtk_write_scalar(fp, "theoretical", exact, nnode);
    st_vtk_write_scalar(fp, "abs_error", abs_error, nnode);
    st_vtk_write_scalar(fp, "source", source, nnode);
    st_vtk_write_scalar(fp, "dirichlet_flag", dirichlet, nnode);
    st_vtk_write_scalar(fp, "mpi_rank", rank_value, nnode);

    fclose(fp);
    free(exact); free(abs_error); free(source); free(dirichlet); free(rank_value);
    return ST_FOM_SUCCESS;
}

static int st_write_vtk_snapshot(
    ST_DIFFUSION_SPACEBLOCK* solver,
    const double* numerical,
    int global_step,
    double time,
    int force_write)
{
    const int rank = monolis_mpi_get_global_my_rank();
    const int nrank = monolis_mpi_get_global_comm_size();
    int status;

    if(!solver->write_vtk){
        return ST_FOM_SUCCESS;
    }
    if(!force_write && global_step % solver->vtk_interval != 0){
        return ST_FOM_SUCCESS;
    }

    status = st_write_vtk_piece(solver, numerical, global_step, time);
    if(status != ST_FOM_SUCCESS){
        return status;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0 && solver->vtk_pvd_fp != NULL){
        for(int part = 0; part < nrank; part++){
            fprintf(
                solver->vtk_pvd_fp,
                "    <DataSet timestep=\"%.17g\" group=\"\" part=\"%d\" "
                "file=\"st_result_%06d_rank%04d.vtk\"/>\n",
                time,
                part,
                global_step,
                part);
        }
        fflush(solver->vtk_pvd_fp);
    }

    return ST_FOM_SUCCESS;
}

static int st_initialize_window_bc(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    const int nnode = fom->fe.total_num_nodes;
    int count = 0;

    solver->bc_window.total_num_nodes = nnode;
    solver->bc_window.block_size = solver->window_dof;

    BBFE_sys_memory_allocation_Dirichlet_bc(
        &solver->bc_window,
        nnode,
        solver->window_dof);
    solver->window_bc_initialized = 1;

    for(int node = 0; node < nnode; node++){
        const int constrained = fom->bc.D_bc_exists[node] ? 1 : 0;
        for(int slab = 0; slab < solver->layout.num_slabs; slab++){
            for(int alpha = 0; alpha < solver->layout.time.n_dof; alpha++){
                const int block = st_block(&solver->layout, slab, alpha);
                const int index = node * solver->window_dof + block;
                solver->bc_window.D_bc_exists[index] = constrained != 0;
                solver->bc_window.imposed_D_val[index] = 0.0;
                if(constrained){
                    count++;
                }
            }
        }
    }

    solver->bc_window.num_D_bcs = count;
    return ST_FOM_SUCCESS;
}

static int st_initialize_slab_bc(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    const int nnode = fom->fe.total_num_nodes;
    const int nt = solver->layout.time.n_dof;
    int count = 0;

    solver->bc_slab.total_num_nodes = nnode;
    solver->bc_slab.block_size = nt;
    BBFE_sys_memory_allocation_Dirichlet_bc(&solver->bc_slab, nnode, nt);
    solver->slab_bc_initialized = 1;

    for(int node = 0; node < nnode; node++){
        const int constrained = fom->bc.D_bc_exists[node] ? 1 : 0;
        for(int alpha = 0; alpha < nt; alpha++){
            const int index = node * nt + alpha;
            solver->bc_slab.D_bc_exists[index] = constrained != 0;
            solver->bc_slab.imposed_D_val[index] = 0.0;
            if(constrained){
                count++;
            }
        }
    }

    solver->bc_slab.num_D_bcs = count;
    return ST_FOM_SUCCESS;
}

static void st_set_window_bc_values(
    ST_DIFFUSION_SPACEBLOCK* solver,
    double time_window_start)
{
    FE_SYSTEM* fom = solver->fom;
    const double dt = fom->vals.dt;

    for(int node = 0; node < fom->fe.total_num_nodes; node++){
        if(!fom->bc.D_bc_exists[node]){
            continue;
        }

        for(int slab = 0; slab < solver->layout.num_slabs; slab++){
            const double slab_start =
                time_window_start + (double)slab * dt;

            for(int alpha = 0; alpha < solver->layout.time.n_dof; alpha++){
                const double tau = solver->layout.time.dof_tau[alpha];
                const double time = slab_start + tau * dt;
                const int block = st_block(&solver->layout, slab, alpha);
                const int index = node * solver->window_dof + block;

                solver->bc_window.imposed_D_val[index] = manusol_get_sol(
                    fom->fe.x[node][0],
                    fom->fe.x[node][1],
                    fom->fe.x[node][2],
                    time);
            }
        }
    }
}

static void st_set_slab_bc_values(
    ST_DIFFUSION_SPACEBLOCK* solver,
    double time_slab_start)
{
    FE_SYSTEM* fom = solver->fom;
    const double dt = fom->vals.dt;
    const int nt = solver->layout.time.n_dof;

    for(int node = 0; node < fom->fe.total_num_nodes; node++){
        if(!fom->bc.D_bc_exists[node]){
            continue;
        }

        for(int alpha = 0; alpha < nt; alpha++){
            const double time = time_slab_start
                + solver->layout.time.dof_tau[alpha] * dt;
            solver->bc_slab.imposed_D_val[node * nt + alpha] =
                manusol_get_sol(
                    fom->fe.x[node][0],
                    fom->fe.x[node][1],
                    fom->fe.x[node][2],
                    time);
        }
    }
}

static int st_read_partition_graph(ST_DIFFUSION_SPACEBLOCK* solver)
{
    const int rank = monolis_mpi_get_global_my_rank();
    const int size = monolis_mpi_get_global_comm_size();
    const char* override_name = getenv("ST_FOM_GRAPH_FILE");
    FILE* fp = NULL;
    int nnode = 0;
    long long nitem_ll = 0;

    if(override_name != NULL && override_name[0] != '\0'){
        snprintf(
            solver->graph_filename,
            sizeof(solver->graph_filename),
            "%s",
            override_name);
    }
    else if(size > 1){
        snprintf(
            solver->graph_filename,
            sizeof(solver->graph_filename),
            "parted.0/graph.dat.%d",
            rank);
    }
    else{
        snprintf(
            solver->graph_filename,
            sizeof(solver->graph_filename),
            "%s",
            "graph.dat");
    }

    fp = fopen(solver->graph_filename, "r");
    if(fp == NULL){
        return ST_FOM_ERR_IO;
    }

    if(fscanf(fp, "%d", &nnode) != 1 || nnode <= 0){
        fclose(fp);
        return ST_FOM_ERR_IO;
    }

    if(nnode != solver->fom->fe.total_num_nodes){
        fprintf(
            stderr,
            "ERROR [rank %d]: graph/node mismatch: graph=%d fe=%d file=%s\n",
            rank,
            nnode,
            solver->fom->fe.total_num_nodes,
            solver->graph_filename);
        fclose(fp);
        return ST_FOM_ERR_DIMENSION;
    }

    for(int row = 0; row < nnode; row++){
        int node_id = -1;
        int degree = -1;

        if(fscanf(fp, "%d %d", &node_id, &degree) != 2 ||
           node_id != row || degree < 0)
        {
            fclose(fp);
            return ST_FOM_ERR_IO;
        }

        for(int j = 0; j < degree; j++){
            int dummy = -1;
            if(fscanf(fp, "%d", &dummy) != 1){
                fclose(fp);
                return ST_FOM_ERR_IO;
            }
        }
        nitem_ll += (long long)degree;
    }
    fclose(fp);
    fp = NULL;

    if(nitem_ll < 0 || nitem_ll > (long long)INT_MAX){
        return ST_FOM_ERR_DIMENSION;
    }

    solver->spatial_graph.num_nodes = nnode;
    solver->spatial_graph.num_items = (int)nitem_ll;
    solver->spatial_graph.index =
        (int*)calloc((size_t)nnode + 1u, sizeof(int));
    solver->spatial_graph.item =
        (int*)calloc((size_t)(nitem_ll > 0 ? nitem_ll : 1), sizeof(int));

    if(solver->spatial_graph.index == NULL || solver->spatial_graph.item == NULL){
        return ST_FOM_ERR_ALLOCATION;
    }

    fp = fopen(solver->graph_filename, "r");
    if(fp == NULL || fscanf(fp, "%d", &nnode) != 1){
        if(fp != NULL) fclose(fp);
        return ST_FOM_ERR_IO;
    }

    {
        int cursor = 0;
        int self_rows = 0;
        solver->spatial_graph.index[0] = 0;

        for(int row = 0; row < nnode; row++){
            int node_id = -1;
            int degree = -1;
            int has_self = 0;

            if(fscanf(fp, "%d %d", &node_id, &degree) != 2 || node_id != row){
                fclose(fp);
                return ST_FOM_ERR_IO;
            }

            for(int j = 0; j < degree; j++){
                int adj = -1;
                if(fscanf(fp, "%d", &adj) != 1 || adj < 0 || adj >= nnode){
                    fclose(fp);
                    return ST_FOM_ERR_IO;
                }
                if(adj == row){
                    has_self = 1;
                }
                solver->spatial_graph.item[cursor++] = adj;
            }
            self_rows += has_self;
            solver->spatial_graph.index[row + 1] = cursor;
        }
        fclose(fp);

        printf(
            "ST-FOM rank %d/%d: using GEDATSU partition graph '%s': "
            "nodes=%d directed_adjacencies=%d self_rows=%d\n",
            rank,
            size,
            solver->graph_filename,
            solver->spatial_graph.num_nodes,
            solver->spatial_graph.num_items,
            self_rows);
        fflush(stdout);
    }

    return ST_FOM_SUCCESS;
}

static int st_initialize_spatial_graph(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    const int rank = monolis_mpi_get_global_my_rank();
    const char* source = getenv("ST_FOM_GRAPH_SOURCE");
    int status;

    status = ST_fom_spatial_graph_initialize(&solver->spatial_graph);
    if(status != ST_FOM_SUCCESS){
        return status;
    }
    solver->graph_initialized = 1;

    if(source == NULL || source[0] == '\0' || strcmp(source, "partition") == 0){
        status = st_read_partition_graph(solver);
        if(status == ST_FOM_SUCCESS){
            return ST_FOM_SUCCESS;
        }

        fprintf(
            stderr,
            "WARNING [rank %d]: could not use GEDATSU partition graph "
            "(status=%d); falling back to FE-connectivity graph.\n",
            rank,
            status);
        ST_fom_spatial_graph_finalize(&solver->spatial_graph);
        ST_fom_spatial_graph_initialize(&solver->spatial_graph);
    }
    else if(strcmp(source, "rebuild") != 0){
        fprintf(
            stderr,
            "ERROR: ST_FOM_GRAPH_SOURCE must be partition or rebuild; got '%s'.\n",
            source);
        return ST_FOM_ERR_ARGUMENT;
    }

    status = ST_fom_spatial_graph_build_from_connectivity(
        &solver->spatial_graph,
        fom->fe.total_num_nodes,
        fom->fe.total_num_elems,
        fom->fe.local_num_nodes,
        fom->fe.conn);
    if(status != ST_FOM_SUCCESS){
        return status;
    }

    status = ST_fom_spatial_graph_validate_connectivity(
        &solver->spatial_graph,
        fom->fe.total_num_nodes,
        fom->fe.total_num_elems,
        fom->fe.local_num_nodes,
        fom->fe.conn);
    if(status != ST_FOM_SUCCESS){
        return status;
    }

    snprintf(
        solver->graph_filename,
        sizeof(solver->graph_filename),
        "%s",
        "FE-connectivity-rebuilt");
    return ST_FOM_SUCCESS;
}

static int st_initialize_matrix_pattern(
    ST_DIFFUSION_SPACEBLOCK* solver,
    MONOLIS* matrix,
    int block_size,
    const char* label)
{
    FE_SYSTEM* fom = solver->fom;
    const int rank = monolis_mpi_get_global_my_rank();
    const int size = monolis_mpi_get_global_comm_size();
    const int local_nodes = fom->fe.total_num_nodes;
    const int owned = solver->num_owned_space_nodes;

    if(matrix == NULL || block_size <= 0 ||
       solver->spatial_graph.num_nodes != local_nodes)
    {
        return ST_FOM_ERR_DIMENSION;
    }

    monolis_initialize(matrix);
    monolis_get_nonzero_pattern_by_nodal_graph_R(
        matrix,
        solver->spatial_graph.num_nodes,
        block_size,
        solver->spatial_graph.index,
        solver->spatial_graph.item);

    matrix->mat.N = owned;
    matrix->mat.NP = local_nodes;

    if(size > 1 && fom->monolis_com.n_internal_vertex != owned){
        fprintf(
            stderr,
            "ERROR [rank %d]: ST ownership/communicator mismatch: "
            "matrix N=%d communicator n_internal_vertex=%d.\n",
            rank,
            owned,
            fom->monolis_com.n_internal_vertex);
        return ST_FOM_ERR_DIMENSION;
    }

    printf(
        "ST-FOM rank %d/%d: %s pattern block=%d "
        "matrix(N,NP)=(%d,%d) graph=%s\n",
        rank,
        size,
        label,
        block_size,
        matrix->mat.N,
        matrix->mat.NP,
        solver->graph_filename);
    fflush(stdout);

    return ST_FOM_SUCCESS;
}

static int st_allocate_vectors(ST_DIFFUSION_SPACEBLOCK* solver)
{
    const size_t window_length = ST_fom_window_vector_length(&solver->layout);
    const size_t spatial_length = (size_t)solver->fom->fe.total_num_nodes;
    const size_t slab_length = spatial_length * (size_t)solver->slab_dof;

    solver->window_solution =
        (double*)calloc(window_length, sizeof(double));
    solver->slab_solution =
        (double*)calloc(slab_length, sizeof(double));
    solver->previous_trace =
        (double*)calloc(spatial_length, sizeof(double));
    solver->trace_work =
        (double*)calloc(spatial_length, sizeof(double));

    if(solver->window_solution == NULL || solver->slab_solution == NULL ||
       solver->previous_trace == NULL || solver->trace_work == NULL)
    {
        return ST_FOM_ERR_ALLOCATION;
    }
    return ST_FOM_SUCCESS;
}

static void st_set_initial_trace_exact(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;

    for(int node = 0; node < fom->fe.total_num_nodes; node++){
        solver->previous_trace[node] = manusol_get_sol(
            fom->fe.x[node][0],
            fom->fe.x[node][1],
            fom->fe.x[node][2],
            0.0);
    }

    monolis_mpi_update_R(
        &fom->monolis_com,
        fom->fe.total_num_nodes,
        1,
        solver->previous_trace);
}

static int st_open_accuracy_file(
    ST_DIFFUSION_SPACEBLOCK* solver,
    const char* filename)
{
    const int rank = monolis_mpi_get_global_my_rank();

    if(filename == NULL || filename[0] == '\0'){
        filename = "st_accuracy.csv";
    }

    snprintf(
        solver->accuracy_filename,
        sizeof(solver->accuracy_filename),
        "%s",
        filename);

    if(rank == 0){
        solver->accuracy_fp = fopen(solver->accuracy_filename, "w");
        if(solver->accuracy_fp == NULL){
            fprintf(stderr, "ERROR: cannot open %s\n", solver->accuracy_filename);
            return ST_FOM_ERR_IO;
        }

        fprintf(
            solver->accuracy_fp,
            "step,time,dt,degree,slabs_per_window,relative_l2,absolute_l2,linf\n");
        fflush(solver->accuracy_fp);
    }
    return ST_FOM_SUCCESS;
}

int ST_diffusion_spaceblock_initialize(
    ST_DIFFUSION_SPACEBLOCK* solver,
    FE_SYSTEM* fom,
    int time_degree,
    int slabs_per_window,
    const char* accuracy_filename)
{
    int status;

    if(solver == NULL || fom == NULL || slabs_per_window <= 0 ||
       fom->fe.total_num_nodes <= 0 || fom->fe.total_num_elems <= 0 ||
       fom->fe.local_num_nodes <= 0 || !isfinite(fom->vals.dt) ||
       fom->vals.dt <= 0.0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    memset(solver, 0, sizeof(*solver));
    solver->fom = fom;
    solver->solve_mode = st_solve_mode_from_environment();

    status = st_resolve_owned_nodes(solver);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    status = ST_fom_layout_initialize(
        &solver->layout,
        fom->fe.total_num_nodes,
        solver->num_owned_space_nodes,
        1,
        slabs_per_window,
        time_degree);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    solver->window_dof = ST_fom_window_dof_per_space_point(&solver->layout);
    solver->slab_dof = solver->layout.time.n_dof;
    if(solver->window_dof <= 0 || solver->slab_dof <= 0){
        status = ST_FOM_ERR_DIMENSION;
        goto fail;
    }

    status = st_initialize_spatial_graph(solver);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    if(solver->solve_mode == ST_FOM_SOLVE_CAUSAL){
        status = st_initialize_slab_bc(solver);
        if(status != ST_FOM_SUCCESS){
            goto fail;
        }
        status = st_initialize_matrix_pattern(
            solver,
            &solver->monolis_slab,
            solver->slab_dof,
            "causal-slab");
        if(status != ST_FOM_SUCCESS){
            goto fail;
        }
        solver->slab_matrix_initialized = 1;
    }
    else{
        status = st_initialize_window_bc(solver);
        if(status != ST_FOM_SUCCESS){
            goto fail;
        }
        status = st_initialize_matrix_pattern(
            solver,
            &solver->monolis_window,
            solver->window_dof,
            "all-at-once-window");
        if(status != ST_FOM_SUCCESS){
            goto fail;
        }
        solver->window_matrix_initialized = 1;
    }

    status = st_allocate_vectors(solver);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    st_set_initial_trace_exact(solver);

    status = st_open_accuracy_file(solver, accuracy_filename);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    status = st_open_vtk_collection(solver);
    if(status != ST_FOM_SUCCESS){
        goto fail;
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-FOM solve mode: %s (window block=%d, active linear block=%d)\n",
            solver->solve_mode == ST_FOM_SOLVE_CAUSAL
                ? "causal block-forward-substitution"
                : "all-at-once window",
            solver->window_dof,
            solver->solve_mode == ST_FOM_SOLVE_CAUSAL
                ? solver->slab_dof : solver->window_dof);
        if(solver->solve_mode == ST_FOM_SOLVE_CAUSAL){
            printf(
                "ST-FOM note: for this linear dG system, causal slab solves "
                "are the exact block forward substitution of the same "
                "lower-triangular all-at-once window equations.\n");
        }
        fflush(stdout);
    }

    return ST_FOM_SUCCESS;

fail:
    ST_diffusion_spaceblock_finalize(solver);
    return status;
}

void ST_diffusion_spaceblock_finalize(ST_DIFFUSION_SPACEBLOCK* solver)
{
    if(solver == NULL){
        return;
    }

    if(solver->accuracy_fp != NULL){
        fclose(solver->accuracy_fp);
        solver->accuracy_fp = NULL;
    }
    st_close_vtk_collection(solver);

    free(solver->window_solution);
    free(solver->slab_solution);
    free(solver->previous_trace);
    free(solver->trace_work);

    if(solver->window_bc_initialized && solver->fom != NULL){
        BBFE_sys_memory_free_Dirichlet_bc(
            &solver->bc_window,
            solver->fom->fe.total_num_nodes,
            solver->window_dof);
    }
    if(solver->slab_bc_initialized && solver->fom != NULL){
        BBFE_sys_memory_free_Dirichlet_bc(
            &solver->bc_slab,
            solver->fom->fe.total_num_nodes,
            solver->slab_dof);
    }

    if(solver->window_matrix_initialized){
        monolis_finalize(&solver->monolis_window);
    }
    if(solver->slab_matrix_initialized){
        monolis_finalize(&solver->monolis_slab);
    }
    if(solver->graph_initialized){
        ST_fom_spatial_graph_finalize(&solver->spatial_graph);
    }

    memset(solver, 0, sizeof(*solver));
}

static int st_assemble_window_matrix(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* mass_ip = (double*)calloc((size_t)np, sizeof(double));
    double* diff_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));
    double* k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(mass_ip == NULL || diff_ip == NULL || jacobian_ip == NULL ||
       local_x == NULL || x_ip == NULL || a_ip == NULL || k_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL){
            goto allocation_error;
        }
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL){
            goto allocation_error;
        }
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int row_node = fe->conn[e][i];

            for(int j = 0; j < nl; j++){
                const int col_node = fe->conn[e][j];
                double mass_value;
                double diffusion_value;

                for(int p = 0; p < np; p++){
                    mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                        basis->N[p][i],
                        basis->N[p][j],
                        a_ip[p]);
                    diff_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        k_ip[p]);
                }

                mass_value = BBFE_std_integ_calc(
                    np, mass_ip, basis->integ_weight, jacobian_ip);
                diffusion_value = BBFE_std_integ_calc(
                    np, diff_ip, basis->integ_weight, jacobian_ip);

                for(int slab = 0; slab < solver->layout.num_slabs; slab++){
                    for(int alpha = 0; alpha < te->n_dof; alpha++){
                        const int row_block = st_block(
                            &solver->layout, slab, alpha);

                        for(int beta = 0; beta < te->n_dof; beta++){
                            const int col_block = st_block(
                                &solver->layout, slab, beta);
                            const double value =
                                ((te->Dt[alpha][beta]
                                  + te->E_self[alpha][beta]) / dt)
                                    * mass_value
                                + te->Mt[alpha][beta] * diffusion_value;

                            if(fabs(value) > 0.0){
                                monolis_add_scalar_to_sparse_matrix_R(
                                    &solver->monolis_window,
                                    row_node,
                                    col_node,
                                    row_block,
                                    col_block,
                                    value);
                            }

                            if(slab > 0 &&
                               fabs(te->E_prev[alpha][beta]) > 0.0)
                            {
                                const int prev_block = st_block(
                                    &solver->layout, slab - 1, beta);
                                const double coupling =
                                    -te->E_prev[alpha][beta]
                                    * mass_value / dt;

                                monolis_add_scalar_to_sparse_matrix_R(
                                    &solver->monolis_window,
                                    row_node,
                                    col_node,
                                    row_block,
                                    prev_block,
                                    coupling);
                            }
                        }
                    }
                }
            }
        }
    }

    for(int p = 0; p < np; p++) free(x_ip[p]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(mass_ip); free(diff_ip); free(jacobian_ip); free(local_x);
    free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(mass_ip); free(diff_ip); free(jacobian_ip); free(local_x);
    free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static int st_add_previous_trace_rhs(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* mass_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));

    if(mass_ip == NULL || jacobian_ip == NULL || local_x == NULL ||
       x_ip == NULL || a_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL) goto allocation_error;
    }

    /*
     * Causal inflow term, assembled with exactly the same spatial mass
     * bilinear form as the left-hand-side matrix:
     *
     *   b_prev(alpha,i) += e_minus(alpha)/dt
     *                      * sum_j M_ij * u_prev(j).
     *
     * This deliberately avoids the interpolated vec_mass path so that the
     * history coupling is algebraically identical to the assembled M block.
     */
    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int row_node = fe->conn[e][i];

            for(int j = 0; j < nl; j++){
                const int col_node = fe->conn[e][j];
                double mass_ij;

                for(int p = 0; p < np; p++){
                    mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                        basis->N[p][i], basis->N[p][j], a_ip[p]);
                }
                mass_ij = BBFE_std_integ_calc(
                    np, mass_ip, basis->integ_weight, jacobian_ip);

                for(int alpha = 0; alpha < te->n_dof; alpha++){
                    const double temporal_coef = te->e_minus[alpha] / dt;
                    if(fabs(temporal_coef) <= 0.0) continue;

                    solver->monolis_window.mat.R.B[
                        (size_t)row_node * (size_t)solver->window_dof
                        + (size_t)st_block(&solver->layout, 0, alpha)] +=
                        temporal_coef * mass_ij
                        * solver->previous_trace[col_node];
                }
            }
        }
    }

    for(int p = 0; p < np; p++) free(x_ip[p]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(mass_ip); free(jacobian_ip); free(local_x); free(x_ip); free(a_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(mass_ip); free(jacobian_ip); free(local_x); free(x_ip); free(a_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static int st_add_source_rhs(
    ST_DIFFUSION_SPACEBLOCK* solver,
    double time_window_start)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* value_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double** v_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));
    double* k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(value_ip == NULL || jacobian_ip == NULL || local_x == NULL ||
       x_ip == NULL || v_ip == NULL || a_ip == NULL || k_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        v_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL || v_ip[p] == NULL) goto allocation_error;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            manusol_get_conv_vel(v_ip[p], x_ip[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int slab = 0; slab < solver->layout.num_slabs; slab++){
            const double t_old = time_window_start + (double)slab * dt;

            for(int alpha = 0; alpha < te->n_dof; alpha++){
                for(int i = 0; i < nl; i++){
                    double integral;

                    for(int p = 0; p < np; p++){
                        double value = 0.0;

                        for(int q = 0; q < te->num_time_ip; q++){
                            const double time = t_old + te->tau[q] * dt;
                            const double f = manusol_get_source(
                                x_ip[p],
                                time,
                                a_ip[p],
                                v_ip[p],
                                k_ip[p]);

                            value += te->weight[q] * te->psi[alpha][q]
                                * BBFE_elemmat_convdiff_vec_source(
                                    basis->N[p][i], f);
                        }
                        value_ip[p] = value;
                    }

                    integral = BBFE_std_integ_calc(
                        np, value_ip, basis->integ_weight, jacobian_ip);

                    solver->monolis_window.mat.R.B[
                        (size_t)fe->conn[e][i] * (size_t)solver->window_dof
                        + (size_t)st_block(
                            &solver->layout, slab, alpha)] += integral;
                }
            }
        }
    }

    for(int p = 0; p < np; p++){ free(x_ip[p]); free(v_ip[p]); }
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(value_ip); free(jacobian_ip); free(local_x); free(x_ip); free(v_ip);
    free(a_ip); free(k_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(v_ip != NULL){ for(int p = 0; p < np; p++) free(v_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(value_ip); free(jacobian_ip); free(local_x); free(x_ip); free(v_ip);
    free(a_ip); free(k_ip);
    return ST_FOM_ERR_ALLOCATION;
}



static int st_report_expected_slab_diagonal(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const int owned = solver->num_owned_space_nodes;
    const int nt = te->n_dof;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    const double dt = fom->vals.dt;
    double* diag = NULL;
    double* mass_ip = NULL;
    double* diff_ip = NULL;
    double* jacobian_ip = NULL;
    double** local_x = NULL;
    double** x_ip = NULL;
    double* a_ip = NULL;
    double* k_ip = NULL;
    double min_abs = DBL_MAX;
    double max_abs = 0.0;
    int exact_zero = 0;

    diag = (double*)calloc((size_t)owned * (size_t)nt, sizeof(double));
    mass_ip = (double*)calloc((size_t)np, sizeof(double));
    diff_ip = (double*)calloc((size_t)np, sizeof(double));
    jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    local_x = (double**)calloc((size_t)nl, sizeof(double*));
    x_ip = (double**)calloc((size_t)np, sizeof(double*));
    a_ip = (double*)calloc((size_t)np, sizeof(double));
    k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(diag == NULL || mass_ip == NULL || diff_ip == NULL ||
       jacobian_ip == NULL || local_x == NULL || x_ip == NULL ||
       a_ip == NULL || k_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL) goto allocation_error;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int node = fe->conn[e][i];
            double mass_value;
            double diffusion_value;
            if(node < 0 || node >= owned) continue;

            for(int p = 0; p < np; p++){
                mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                    basis->N[p][i], basis->N[p][i], a_ip[p]);
                diff_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
                    fe->geo[e][p].grad_N[i],
                    fe->geo[e][p].grad_N[i],
                    k_ip[p]);
            }
            mass_value = BBFE_std_integ_calc(
                np, mass_ip, basis->integ_weight, jacobian_ip);
            diffusion_value = BBFE_std_integ_calc(
                np, diff_ip, basis->integ_weight, jacobian_ip);

            for(int alpha = 0; alpha < nt; alpha++){
                diag[(size_t)node * (size_t)nt + (size_t)alpha] +=
                    ((te->Dt[alpha][alpha] + te->E_self[alpha][alpha]) / dt)
                        * mass_value
                    + te->Mt[alpha][alpha] * diffusion_value;
            }
        }
    }

    for(int node = 0; node < owned; node++){
        for(int alpha = 0; alpha < nt; alpha++){
            const double a = fabs(diag[(size_t)node * (size_t)nt + (size_t)alpha]);
            if(a == 0.0) exact_zero++;
            if(a < min_abs) min_abs = a;
            if(a > max_abs) max_abs = a;
        }
    }

    printf(
        "ST-FOM rank %d expected dG slab scalar diagonal: "
        "min|a_ii|=%.15e max|a_ii|=%.15e exact_zero=%d/%d\n",
        monolis_mpi_get_global_my_rank(),
        min_abs,
        max_abs,
        exact_zero,
        owned * nt);
    fflush(stdout);

    for(int p = 0; p < np; p++) free(x_ip[p]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(diag); free(mass_ip); free(diff_ip); free(jacobian_ip);
    free(local_x); free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(diag); free(mass_ip); free(diff_ip); free(jacobian_ip);
    free(local_x); free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static int st_assemble_slab_matrix(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* mass_ip = (double*)calloc((size_t)np, sizeof(double));
    double* diff_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));
    double* k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(mass_ip == NULL || diff_ip == NULL || jacobian_ip == NULL ||
       local_x == NULL || x_ip == NULL || a_ip == NULL || k_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL) goto allocation_error;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int row_node = fe->conn[e][i];
            for(int j = 0; j < nl; j++){
                const int col_node = fe->conn[e][j];
                double mass_value;
                double diffusion_value;

                for(int p = 0; p < np; p++){
                    mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                        basis->N[p][i], basis->N[p][j], a_ip[p]);
                    diff_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
                        fe->geo[e][p].grad_N[i],
                        fe->geo[e][p].grad_N[j],
                        k_ip[p]);
                }

                mass_value = BBFE_std_integ_calc(
                    np, mass_ip, basis->integ_weight, jacobian_ip);
                diffusion_value = BBFE_std_integ_calc(
                    np, diff_ip, basis->integ_weight, jacobian_ip);

                for(int alpha = 0; alpha < te->n_dof; alpha++){
                    for(int beta = 0; beta < te->n_dof; beta++){
                        const double value =
                            ((te->Dt[alpha][beta]
                              + te->E_self[alpha][beta]) / dt)
                                * mass_value
                            + te->Mt[alpha][beta] * diffusion_value;

                        if(fabs(value) > 0.0){
                            monolis_add_scalar_to_sparse_matrix_R(
                                &solver->monolis_slab,
                                row_node,
                                col_node,
                                alpha,
                                beta,
                                value);
                        }
                    }
                }
            }
        }
    }

    for(int p = 0; p < np; p++) free(x_ip[p]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(mass_ip); free(diff_ip); free(jacobian_ip); free(local_x);
    free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(mass_ip); free(diff_ip); free(jacobian_ip); free(local_x);
    free(x_ip); free(a_ip); free(k_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static int st_add_previous_trace_rhs_slab(ST_DIFFUSION_SPACEBLOCK* solver)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* mass_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));

    if(mass_ip == NULL || jacobian_ip == NULL || local_x == NULL ||
       x_ip == NULL || a_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL) goto allocation_error;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
        }

        for(int i = 0; i < nl; i++){
            const int row_node = fe->conn[e][i];

            for(int j = 0; j < nl; j++){
                const int col_node = fe->conn[e][j];
                double mass_ij;

                for(int p = 0; p < np; p++){
                    mass_ip[p] = BBFE_elemmat_convdiff_mat_mass(
                        basis->N[p][i], basis->N[p][j], a_ip[p]);
                }
                mass_ij = BBFE_std_integ_calc(
                    np, mass_ip, basis->integ_weight, jacobian_ip);

                for(int alpha = 0; alpha < te->n_dof; alpha++){
                    const double temporal_coef = te->e_minus[alpha] / dt;
                    if(fabs(temporal_coef) <= 0.0) continue;

                    solver->monolis_slab.mat.R.B[
                        (size_t)row_node * (size_t)solver->slab_dof
                        + (size_t)alpha] +=
                        temporal_coef * mass_ij
                        * solver->previous_trace[col_node];
                }
            }
        }
    }

    for(int p = 0; p < np; p++) free(x_ip[p]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(mass_ip); free(jacobian_ip); free(local_x); free(x_ip); free(a_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(mass_ip); free(jacobian_ip); free(local_x); free(x_ip); free(a_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static int st_add_source_rhs_slab(
    ST_DIFFUSION_SPACEBLOCK* solver,
    double time_slab_start)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const ST_FOM_TIME_ELEMENT* te = &solver->layout.time;
    const double dt = fom->vals.dt;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    double* value_ip = (double*)calloc((size_t)np, sizeof(double));
    double* jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    double** local_x = (double**)calloc((size_t)nl, sizeof(double*));
    double** x_ip = (double**)calloc((size_t)np, sizeof(double*));
    double** v_ip = (double**)calloc((size_t)np, sizeof(double*));
    double* a_ip = (double*)calloc((size_t)np, sizeof(double));
    double* k_ip = (double*)calloc((size_t)np, sizeof(double));

    if(value_ip == NULL || jacobian_ip == NULL || local_x == NULL ||
       x_ip == NULL || v_ip == NULL || a_ip == NULL || k_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL) goto allocation_error;
    }
    for(int p = 0; p < np; p++){
        x_ip[p] = (double*)calloc(3u, sizeof(double));
        v_ip[p] = (double*)calloc(3u, sizeof(double));
        if(x_ip[p] == NULL || v_ip[p] == NULL) goto allocation_error;
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for(int p = 0; p < np; p++){
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
            manusol_get_conv_vel(v_ip[p], x_ip[p]);
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }

        for(int alpha = 0; alpha < te->n_dof; alpha++){
            for(int i = 0; i < nl; i++){
                for(int p = 0; p < np; p++){
                    double value = 0.0;
                    for(int q = 0; q < te->num_time_ip; q++){
                        const double time =
                            time_slab_start + te->tau[q] * dt;
                        const double f = manusol_get_source(
                            x_ip[p], time, a_ip[p], v_ip[p], k_ip[p]);
                        value += te->weight[q] * te->psi[alpha][q]
                            * BBFE_elemmat_convdiff_vec_source(
                                basis->N[p][i], f);
                    }
                    value_ip[p] = value;
                }

                solver->monolis_slab.mat.R.B[
                    (size_t)fe->conn[e][i] * (size_t)solver->slab_dof
                    + (size_t)alpha] += BBFE_std_integ_calc(
                        np, value_ip, basis->integ_weight, jacobian_ip);
            }
        }
    }

    for(int p = 0; p < np; p++){ free(x_ip[p]); free(v_ip[p]); }
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(value_ip); free(jacobian_ip); free(local_x); free(x_ip); free(v_ip);
    free(a_ip); free(k_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int p = 0; p < np; p++) free(x_ip[p]); }
    if(v_ip != NULL){ for(int p = 0; p < np; p++) free(v_ip[p]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(value_ip); free(jacobian_ip); free(local_x); free(x_ip); free(v_ip);
    free(a_ip); free(k_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static void st_copy_slab_solution_to_window(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int slab)
{
    const int nnode = solver->fom->fe.total_num_nodes;
    const int nt = solver->layout.time.n_dof;

    for(int node = 0; node < nnode; node++){
        for(int alpha = 0; alpha < nt; alpha++){
            solver->window_solution[
                ST_fom_vector_index(&solver->layout, node, slab, alpha, 0)] =
                solver->slab_solution[(size_t)node * (size_t)nt + (size_t)alpha];
        }
    }
}

static void st_extract_slab_solution_right_trace(
    ST_DIFFUSION_SPACEBLOCK* solver)
{
    const int nnode = solver->fom->fe.total_num_nodes;
    const int nt = solver->layout.time.n_dof;

    for(int node = 0; node < nnode; node++){
        double value = 0.0;
        for(int alpha = 0; alpha < nt; alpha++){
            value += solver->layout.time.e_plus[alpha]
                * solver->slab_solution[
                    (size_t)node * (size_t)nt + (size_t)alpha];
        }
        solver->previous_trace[node] = value;
    }
}

static int st_report_slab_accuracy(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int global_step,
    int slab,
    double time)
{
    FE_SYSTEM* fom = solver->fom;
    BBFE_DATA* fe = &fom->fe;
    BBFE_BASIS* basis = &fom->basis;
    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;
    double local_diff2 = 0.0;
    double local_ref2 = 0.0;
    double local_linf = 0.0;
    double global_diff2 = 0.0;
    double global_ref2 = 0.0;
    double global_linf = 0.0;

    double* diff2_ip = NULL;
    double* ref2_ip = NULL;
    double* jacobian_ip = NULL;
    double* local_value = NULL;
    double** local_x = NULL;
    double** x_ip = NULL;

    if(ST_fom_extract_slab_right_trace(
        &solver->layout,
        slab,
        solver->window_solution,
        solver->trace_work) != ST_FOM_SUCCESS)
    {
        return ST_FOM_ERR_DIMENSION;
    }

    /* Linf is evaluated on unique/owned spatial nodes. */
    for(int node = 0; node < solver->num_owned_space_nodes; node++){
        const double exact = manusol_get_sol(
            fe->x[node][0],
            fe->x[node][1],
            fe->x[node][2],
            time);
        const double abs_diff = fabs(solver->trace_work[node] - exact);
        if(abs_diff > local_linf){
            local_linf = abs_diff;
        }
    }

    /*
     * True FE spatial L2 norm: interpolate T_h to integration points and
     * integrate (T_h - T_exact)^2 over each local partition element.
     */
    diff2_ip = (double*)calloc((size_t)np, sizeof(double));
    ref2_ip = (double*)calloc((size_t)np, sizeof(double));
    jacobian_ip = (double*)calloc((size_t)np, sizeof(double));
    local_value = (double*)calloc((size_t)nl, sizeof(double));
    local_x = (double**)calloc((size_t)nl, sizeof(double*));
    x_ip = (double**)calloc((size_t)np, sizeof(double*));

    if(diff2_ip == NULL || ref2_ip == NULL || jacobian_ip == NULL ||
       local_value == NULL || local_x == NULL || x_ip == NULL)
    {
        goto allocation_error;
    }

    for(int i = 0; i < nl; i++){
        local_x[i] = (double*)calloc(3u, sizeof(double));
        if(local_x[i] == NULL){
            goto allocation_error;
        }
    }
    for(int q = 0; q < np; q++){
        x_ip[q] = (double*)calloc(3u, sizeof(double));
        if(x_ip[q] == NULL){
            goto allocation_error;
        }
    }

    for(int e = 0; e < fe->total_num_elems; e++){
        BBFE_elemmat_set_Jacobian_array(jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
        BBFE_elemmat_set_local_array_scalar(
            local_value, fe, solver->trace_work, e);

        for(int q = 0; q < np; q++){
            const double numerical = BBFE_std_mapping_scalar(
                nl, local_value, basis->N[q]);
            double exact;
            double diff;

            BBFE_std_mapping_vector3d(x_ip[q], nl, local_x, basis->N[q]);
            exact = manusol_get_sol(
                x_ip[q][0], x_ip[q][1], x_ip[q][2], time);
            diff = numerical - exact;

            diff2_ip[q] = diff * diff;
            ref2_ip[q] = exact * exact;
        }

        local_diff2 += BBFE_std_integ_calc(
            np, diff2_ip, basis->integ_weight, jacobian_ip);
        local_ref2 += BBFE_std_integ_calc(
            np, ref2_ip, basis->integ_weight, jacobian_ip);
    }

    MPI_Allreduce(
        &local_diff2, &global_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_ref2, &global_ref2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_linf, &global_linf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(monolis_mpi_get_global_my_rank() == 0){
        const double abs_l2 = sqrt(fmax(0.0, global_diff2));
        const double rel_l2 = global_ref2 > DBL_MIN
            ? sqrt(fmax(0.0, global_diff2) / global_ref2)
            : abs_l2;

        printf(
            "ST ACCURACY step=%d time=%.12e rel_L2=%.15e "
            "abs_L2=%.15e Linf=%.15e\n",
            global_step,
            time,
            rel_l2,
            abs_l2,
            global_linf);

        if(solver->accuracy_fp != NULL){
            fprintf(
                solver->accuracy_fp,
                "%d,%.17g,%.17g,%d,%d,%.17g,%.17g,%.17g\n",
                global_step,
                time,
                fom->vals.dt,
                solver->layout.time.degree,
                solver->layout.num_slabs,
                rel_l2,
                abs_l2,
                global_linf);
            fflush(solver->accuracy_fp);
        }
    }

    for(int q = 0; q < np; q++) free(x_ip[q]);
    for(int i = 0; i < nl; i++) free(local_x[i]);
    free(diff2_ip); free(ref2_ip); free(jacobian_ip);
    free(local_value); free(local_x); free(x_ip);
    return ST_FOM_SUCCESS;

allocation_error:
    if(x_ip != NULL){ for(int q = 0; q < np; q++) free(x_ip[q]); }
    if(local_x != NULL){ for(int i = 0; i < nl; i++) free(local_x[i]); }
    free(diff2_ip); free(ref2_ip); free(jacobian_ip);
    free(local_value); free(local_x); free(x_ip);
    return ST_FOM_ERR_ALLOCATION;
}

static double st_owned_l2_norm(
    const ST_DIFFUSION_SPACEBLOCK* solver,
    const double* values,
    int dof)
{
    double local_sum = 0.0;
    double global_sum = 0.0;

    if(solver == NULL || values == NULL || dof <= 0){
        return 0.0;
    }

    for(int node = 0; node < solver->num_owned_space_nodes; node++){
        for(int k = 0; k < dof; k++){
            const double v = values[(size_t)node * (size_t)dof + (size_t)k];
            local_sum += v * v;
        }
    }

    MPI_Allreduce(
        &local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(fmax(0.0, global_sum));
}


/* -------------------------------------------------------------------------
 * Discrete manufactured-solution residual diagnostic.
 *
 * This checks the exact manufactured nodal field against the ACTUAL assembled
 * slab matrix/RHS.  Before Dirichlet elimination we report only unconstrained
 * owned rows.  After Dirichlet elimination we report all owned rows.
 *
 * If pre-BC residual is large, the source/history/operator assembly is
 * inconsistent.  If pre-BC is small but post-BC becomes large, the
 * Dirichlet elimination/path is the problem.
 * ------------------------------------------------------------------------- */
static double st_discrete_exact_residual_relative(
    ST_DIFFUSION_SPACEBLOCK* solver,
    double time_slab_start,
    int include_dirichlet_rows)
{
    FE_SYSTEM* fom;
    const int nt = solver != NULL ? solver->layout.time.n_dof : 0;
    const int nnode = solver != NULL && solver->fom != NULL
        ? solver->fom->fe.total_num_nodes : 0;
    const size_t n = (size_t)nnode * (size_t)nt;
    double* exact = NULL;
    double* Ax = NULL;
    double local_r2 = 0.0;
    double local_b2 = 0.0;
    double global_r2 = 0.0;
    double global_b2 = 0.0;

    if(solver == NULL || solver->fom == NULL || nt <= 0 || nnode <= 0){
        return HUGE_VAL;
    }

    fom = solver->fom;
    exact = (double*)calloc(n, sizeof(double));
    Ax = (double*)calloc(n, sizeof(double));
    if(exact == NULL || Ax == NULL){
        free(exact);
        free(Ax);
        return HUGE_VAL;
    }

    for(int node = 0; node < nnode; node++){
        for(int alpha = 0; alpha < nt; alpha++){
            const double time = time_slab_start
                + solver->layout.time.dof_tau[alpha] * fom->vals.dt;
            exact[(size_t)node * (size_t)nt + (size_t)alpha] =
                manusol_get_sol(
                    fom->fe.x[node][0],
                    fom->fe.x[node][1],
                    fom->fe.x[node][2],
                    time);
        }
    }

    monolis_mpi_update_R(&fom->monolis_com, nnode, nt, exact);
    monolis_matvec_product_R(
        &solver->monolis_slab,
        &fom->monolis_com,
        exact,
        Ax);

    for(int node = 0; node < solver->num_owned_space_nodes; node++){
        const int constrained = fom->bc.D_bc_exists[node] ? 1 : 0;
        if(!include_dirichlet_rows && constrained){
            continue;
        }
        for(int alpha = 0; alpha < nt; alpha++){
            const size_t idx = (size_t)node * (size_t)nt + (size_t)alpha;
            const double b = solver->monolis_slab.mat.R.B[idx];
            const double r = Ax[idx] - b;
            local_r2 += r * r;
            local_b2 += b * b;
        }
    }

    MPI_Allreduce(
        &local_r2, &global_r2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_b2, &global_b2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    free(exact);
    free(Ax);

    if(global_b2 <= DBL_MIN){
        return sqrt(fmax(0.0, global_r2));
    }
    return sqrt(fmax(0.0, global_r2) / global_b2);
}

static int st_window_linear_solver_type(void)
{
    const char* value = getenv("ST_FOM_LINEAR_SOLVER");

    /* Match the fluid ST-FOM default: BiCGSAFE for the nonsymmetric dG system. */
    if(value == NULL || value[0] == '\0' || strcmp(value, "bicgsafe") == 0){
        return MONOLIS_ITER_BICGSAFE;
    }
    if(strcmp(value, "bicgstab") == 0){
        return MONOLIS_ITER_BICGSTAB;
    }

    fprintf(
        stderr,
        "ERROR: ST_FOM_LINEAR_SOLVER must be bicgsafe or bicgstab; got '%s'.\n",
        value);
    return MONOLIS_ITER_BICGSAFE;
}

static const char* st_window_linear_solver_name(int solver_type)
{
    return solver_type == MONOLIS_ITER_BICGSTAB ? "BiCGSTAB" : "BiCGSAFE";
}

static int st_solve_window_all_at_once(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int window_id,
    double time_window_start)
{
    FE_SYSTEM* fom = solver->fom;
    const size_t window_length = ST_fom_window_vector_length(&solver->layout);
    int status;

    /*
     * Reassemble the ST operator for every window.  Clear BOTH A and RHS.
     * monolis_clear_mat_value_rhs_R() clears only mat.R.B; without the
     * explicit matrix clear, A accumulates window by window.
     */
    monolis_clear_mat_value_R(&solver->monolis_window);
    monolis_clear_mat_value_rhs_R(&solver->monolis_window);
    memset(solver->window_solution, 0, window_length * sizeof(double));

    status = st_assemble_window_matrix(solver);
    if(status != ST_FOM_SUCCESS) return status;

    status = st_add_previous_trace_rhs(solver);
    if(status != ST_FOM_SUCCESS) return status;

    status = st_add_source_rhs(solver, time_window_start);
    if(status != ST_FOM_SUCCESS) return status;

    st_set_window_bc_values(solver, time_window_start);
    BBFE_sys_monowrap_set_Dirichlet_bc(
        &solver->monolis_window,
        fom->fe.total_num_nodes,
        solver->window_dof,
        &solver->bc_window,
        solver->monolis_window.mat.R.B);

    {
        const int linear_solver = st_window_linear_solver_type();
        if(window_id == 0 && monolis_mpi_get_global_my_rank() == 0){
            printf(
                "ST-FOM linear solve mode=window: N=%d NP=%d NDOF=%d, "
                "solver=%s, precond=DIAG\n",
                solver->monolis_window.mat.N,
                solver->monolis_window.mat.NP,
                solver->window_dof,
                st_window_linear_solver_name(linear_solver));
            fflush(stdout);
        }

        BBFE_sys_monowrap_solve(
            &solver->monolis_window,
            &fom->monolis_com,
            solver->window_solution,
            linear_solver,
            MONOLIS_PREC_DIAG,
            fom->vals.mat_max_iter,
            fom->vals.mat_epsilon);
    }

    monolis_mpi_update_R(
        &fom->monolis_com,
        fom->fe.total_num_nodes,
        solver->window_dof,
        solver->window_solution);

    for(int slab = 0; slab < solver->layout.num_slabs; slab++){
        const int global_step =
            window_id * solver->layout.num_slabs + slab + 1;
        const double time = (double)global_step * fom->vals.dt;

        status = st_report_slab_accuracy(solver, global_step, slab, time);
        if(status != ST_FOM_SUCCESS) return status;

        status = st_write_vtk_snapshot(
            solver, solver->trace_work, global_step, time, 0);
        if(status != ST_FOM_SUCCESS) return status;
    }

    status = ST_fom_extract_last_right_trace(
        &solver->layout,
        solver->window_solution,
        solver->previous_trace);
    if(status != ST_FOM_SUCCESS) return status;

    monolis_mpi_update_R(
        &fom->monolis_com,
        fom->fe.total_num_nodes,
        1,
        solver->previous_trace);

    return ST_FOM_SUCCESS;
}

static int st_solve_window_causal(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int window_id,
    double time_window_start)
{
    FE_SYSTEM* fom = solver->fom;
    const int nt = solver->layout.time.n_dof;
    const int local_nodes = fom->fe.total_num_nodes;
    const size_t window_length = ST_fom_window_vector_length(&solver->layout);
    const size_t slab_length = (size_t)local_nodes * (size_t)nt;
    const int linear_solver = st_window_linear_solver_type();

    memset(solver->window_solution, 0, window_length * sizeof(double));

    if(window_id == 0 && monolis_mpi_get_global_my_rank() == 0){
        printf(
            "ST-FOM linear solve mode=causal: N=%d NP=%d NDOF=%d "
            "(dG time block only), solver=%s, precond=DIAG\n",
            solver->monolis_slab.mat.N,
            solver->monolis_slab.mat.NP,
            nt,
            st_window_linear_solver_name(linear_solver));
        printf(
            "ST-FOM causal equivalence: window NDOF=%d is solved by exact "
            "block forward substitution over %d lower-triangular slabs.\n",
            solver->window_dof,
            solver->layout.num_slabs);
        fflush(stdout);
    }

    for(int slab = 0; slab < solver->layout.num_slabs; slab++){
        const int global_step =
            window_id * solver->layout.num_slabs + slab + 1;
        const double time_slab_start =
            time_window_start + (double)slab * fom->vals.dt;
        const double time_right =
            time_slab_start + fom->vals.dt;
        int status;

        /*
         * The slab operator is assembled afresh at every causal step.
         * Clear matrix values as well as RHS; otherwise the same element
         * matrix is added repeatedly (A, 2A, 3A, ...), which explains why
         * step 1 is accurate and step 2 onward collapses toward zero.
         */
        monolis_clear_mat_value_R(&solver->monolis_slab);
        monolis_clear_mat_value_rhs_R(&solver->monolis_slab);
        memset(solver->slab_solution, 0, slab_length * sizeof(double));

        status = st_assemble_slab_matrix(solver);
        if(status != ST_FOM_SUCCESS) return status;

        if(window_id == 0 && slab == 0){
            status = st_report_expected_slab_diagonal(solver);
            if(status != ST_FOM_SUCCESS) return status;
        }

        status = st_add_previous_trace_rhs_slab(solver);
        if(status != ST_FOM_SUCCESS) return status;

        if(global_step <= 3){
            const double prev_norm = st_owned_l2_norm(
                solver, solver->previous_trace, 1);
            const double rhs_prev_norm = st_owned_l2_norm(
                solver, solver->monolis_slab.mat.R.B, nt);
            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "ST HISTORY step=%d before source: "
                    "||u_prev||_2=%.15e ||b_prev||_2=%.15e\n",
                    global_step, prev_norm, rhs_prev_norm);
                fflush(stdout);
            }
        }

        status = st_add_source_rhs_slab(solver, time_slab_start);
        if(status != ST_FOM_SUCCESS) return status;

        if(global_step <= 3){
            const double rhs_total_norm = st_owned_l2_norm(
                solver, solver->monolis_slab.mat.R.B, nt);
            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "ST HISTORY step=%d after source: "
                    "||b_prev+b_source||_2=%.15e\n",
                    global_step, rhs_total_norm);
                fflush(stdout);
            }
        }

        if(global_step <= 2){
            const double rel_pre = st_discrete_exact_residual_relative(
                solver, time_slab_start, 0);
            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "ST DISCRETE MMS step=%d preBC interior relative residual="
                    "%.15e\n",
                    global_step, rel_pre);
                fflush(stdout);
            }
        }

        st_set_slab_bc_values(solver, time_slab_start);
        BBFE_sys_monowrap_set_Dirichlet_bc(
            &solver->monolis_slab,
            local_nodes,
            nt,
            &solver->bc_slab,
            solver->monolis_slab.mat.R.B);

        if(global_step <= 2){
            const double rel_post = st_discrete_exact_residual_relative(
                solver, time_slab_start, 1);
            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "ST DISCRETE MMS step=%d postBC all-row relative residual="
                    "%.15e\n",
                    global_step, rel_post);
                fflush(stdout);
            }
        }

        BBFE_sys_monowrap_solve(
            &solver->monolis_slab,
            &fom->monolis_com,
            solver->slab_solution,
            linear_solver,
            MONOLIS_PREC_DIAG,
            fom->vals.mat_max_iter,
            fom->vals.mat_epsilon);

        monolis_mpi_update_R(
            &fom->monolis_com,
            local_nodes,
            nt,
            solver->slab_solution);

        st_copy_slab_solution_to_window(solver, slab);

        status = st_report_slab_accuracy(
            solver,
            global_step,
            slab,
            time_right);
        if(status != ST_FOM_SUCCESS) return status;

        /*
         * st_owned_l2_norm() contains MPI_Allreduce, so every rank MUST call it.
         * Only the printing is restricted to rank 0.  The previous v8 code put
         * the collective itself inside the rank-0 branch and deadlocked here
         * immediately after the step-1 accuracy report.
         */
        if(global_step <= 3){
            const double trace_norm = st_owned_l2_norm(
                solver, solver->trace_work, 1);
            if(monolis_mpi_get_global_my_rank() == 0){
                printf(
                    "ST HISTORY step=%d right trace: ||u_right||_2=%.15e\n",
                    global_step, trace_norm);
                fflush(stdout);
            }
        }

        status = st_write_vtk_snapshot(
            solver, solver->trace_work, global_step, time_right, 0);
        if(status != ST_FOM_SUCCESS) return status;

        /*
         * trace_work is the exact right trace that was just used for the
         * accuracy report/VTK. Carry that identical field into the next slab.
         */
        memcpy(
            solver->previous_trace,
            solver->trace_work,
            (size_t)local_nodes * sizeof(double));
        monolis_mpi_update_R(
            &fom->monolis_com,
            local_nodes,
            1,
            solver->previous_trace);
    }

    return ST_FOM_SUCCESS;
}

int ST_diffusion_spaceblock_solve_window(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int window_id,
    double time_window_start)
{
    const size_t window_length = solver != NULL
        ? ST_fom_window_vector_length(&solver->layout) : 0u;

    if(solver == NULL || solver->fom == NULL || window_id < 0 ||
       !isfinite(time_window_start) || window_length == 0u)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    if(solver->solve_mode == ST_FOM_SOLVE_CAUSAL){
        if(!solver->slab_matrix_initialized || !solver->slab_bc_initialized){
            return ST_FOM_ERR_DIMENSION;
        }
        return st_solve_window_causal(solver, window_id, time_window_start);
    }

    if(!solver->window_matrix_initialized || !solver->window_bc_initialized){
        return ST_FOM_ERR_DIMENSION;
    }
    return st_solve_window_all_at_once(solver, window_id, time_window_start);
}

int ST_diffusion_spaceblock_run(
    ST_DIFFUSION_SPACEBLOCK* solver,
    int num_windows)
{
    if(solver == NULL || solver->fom == NULL || num_windows <= 0){
        return ST_FOM_ERR_ARGUMENT;
    }

    if(monolis_mpi_get_global_my_rank() == 0){
        printf("\n============================================================\n");
        printf(" Normal space-time FOM (no ROM)\n");
        printf(" time discretization : dG(%d)\n", solver->layout.time.degree);
        printf(" slabs/window        : %d\n", solver->layout.num_slabs);
        printf(" dt                  : %.15e\n", solver->fom->vals.dt);
        printf(" number of windows   : %d\n", num_windows);
        printf(" accuracy CSV        : %s\n", solver->accuracy_filename);
        if(solver->write_vtk){
            printf(" VTK/PVD             : %s (interval=%d)\n",
                solver->vtk_pvd_filename, solver->vtk_interval);
        }
        printf("============================================================\n\n");
    }

    {
        int status = st_write_vtk_snapshot(
            solver, solver->previous_trace, 0, 0.0, 1);
        if(status != ST_FOM_SUCCESS){
            return status;
        }
    }

    for(int window = 0; window < num_windows; window++){
        const double time_start =
            (double)(window * solver->layout.num_slabs)
            * solver->fom->vals.dt;
        int status;

        if(monolis_mpi_get_global_my_rank() == 0){
            printf(
                "---- ST window %d / %d : [%.12e, %.12e] ----\n",
                window + 1,
                num_windows,
                time_start,
                time_start
                    + solver->layout.num_slabs * solver->fom->vals.dt);
        }

        status = ST_diffusion_spaceblock_solve_window(
            solver, window, time_start);
        if(status != ST_FOM_SUCCESS){
            return status;
        }
    }

    return ST_FOM_SUCCESS;
}
