/*
 * st_global_mesh_generator.c
 *
 * Generate a global windowed space-time dG graph and ST-cell connectivity
 * from an ordinary global spatial node.dat and elem.dat.
 *
 * Supported temporal elements:
 *   degree = 0 : dG(0)
 *   degree = 1 : dG(1)
 *   degree = 2 : dG(2)
 *
 * Global ST node numbering (time-major, matching the current solver):
 *
 *   G_ST(s, alpha, i) = (s * nt + alpha) * Nx + i
 *
 * where
 *   s     : slab ID in one time window
 *   alpha : temporal DOF ID
 *   i     : global spatial node ID
 *   nt    : degree + 1
 *   Nx    : number of global spatial nodes
 *
 * Input formats:
 *
 *   node.dat
 *     num_nodes
 *     x y z
 *     ...
 *
 *   spatial elem.dat
 *     num_elems num_local_nodes
 *     node_0 ... node_(num_local_nodes-1)
 *     ...
 *
 * Output files in output_dir:
 *
 *   graph.dat
 *     Symmetric ST graph. Intended for graph partitioning and also usable
 *     as a safe superset of the matrix nonzero pattern.
 *
 *   matrix_graph.dat
 *     Exact directed row graph of the current dG matrix structure.
 *
 *   elem.dat
 *     Standard fixed-width connectivity for ST cells. Each ST cell is one
 *     pair (slab, spatial element) and contains all same-slab temporal DOFs:
 *
 *       nt * num_local_nodes ST node IDs.
 *
 *   st_elem_map.dat
 *     Metadata requested for assembly:
 *
 *       ST cell ID, slab ID, global spatial element ID,
 *       global spatial node IDs.
 *
 *   st_node_map.dat
 *     global ST node ID, slab ID, temporal DOF ID, global spatial node ID.
 *
 *   spatial_graph.dat
 *     Spatial nodal graph reconstructed from elem.dat.
 *
 *   st_meta.dat
 *     Generation metadata and numbering convention.
 *
 * Compile:
 *   cc -O2 -std=c11 -Wall -Wextra -pedantic \
 *      st_global_mesh_generator.c -o st_global_mesh_generator
 *
 * Run:
 *   ./st_global_mesh_generator DEGREE NUM_WINDOW_SLABS \
 *      node.dat elem.dat OUTPUT_DIR
 *
 * Example:
 *   ./st_global_mesh_generator 1 4 ./node.dat ./elem.dat ./st_global
 */

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
    int*   data;
    size_t size;
    size_t capacity;
} INT_VECTOR;

typedef struct {
    int num_nodes;
    int num_elems;
    int num_local_nodes;
    int* conn;
    INT_VECTOR* adjacency;
} SPATIAL_MESH;

static void fail(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static void fail_path(const char* action, const char* path)
{
    fprintf(stderr, "ERROR: %s: %s\n", action, path);
    exit(EXIT_FAILURE);
}

static void* checked_calloc(size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        fail("allocation size overflow");
    }

    void* p = calloc(count, size);
    if(p == NULL) {
        fail("memory allocation failed");
    }
    return p;
}

static void* checked_realloc(void* old_ptr, size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        fail("reallocation size overflow");
    }

    void* p = realloc(old_ptr, count * size);
    if(p == NULL) {
        fail("memory reallocation failed");
    }
    return p;
}

static int compare_int(const void* lhs, const void* rhs)
{
    const int a = *(const int*)lhs;
    const int b = *(const int*)rhs;
    return (a > b) - (a < b);
}

static void int_vector_push(INT_VECTOR* v, int value)
{
    if(v->size == v->capacity) {
        size_t new_capacity = (v->capacity == 0) ? 8 : 2 * v->capacity;
        if(new_capacity < v->capacity) {
            fail("adjacency capacity overflow");
        }

        v->data = (int*)checked_realloc(
            v->data,
            new_capacity,
            sizeof(int));
        v->capacity = new_capacity;
    }

    v->data[v->size++] = value;
}

static void int_vector_sort_unique(INT_VECTOR* v)
{
    if(v->size == 0) {
        return;
    }

    qsort(v->data, v->size, sizeof(int), compare_int);

    size_t write_pos = 1;
    for(size_t read_pos = 1; read_pos < v->size; read_pos++) {
        if(v->data[read_pos] != v->data[write_pos - 1]) {
            v->data[write_pos++] = v->data[read_pos];
        }
    }
    v->size = write_pos;
}

static int64_t checked_mul_i64(int64_t a, int64_t b, const char* what)
{
    if(a < 0 || b < 0) {
        fail("negative value in checked multiplication");
    }

    if(a != 0 && b > INT64_MAX / a) {
        fprintf(stderr, "ERROR: int64 overflow while computing %s\n", what);
        exit(EXIT_FAILURE);
    }
    return a * b;
}

static int64_t st_gid(
    int slab,
    int alpha,
    int spatial_node,
    int nt,
    int num_spatial_nodes)
{
    return
        ((int64_t)slab * (int64_t)nt + (int64_t)alpha)
        * (int64_t)num_spatial_nodes
        + (int64_t)spatial_node;
}

/*
 * Structural temporal pattern of the validated nodal dG(0)-dG(2) basis.
 *
 * Same slab:
 *   ((Dt + E_self)/dt) tensor M + Mt tensor K
 *
 * For degree 0, 1, and 2 currently implemented in core_FOM_ST.c,
 * every alpha-beta pair is structurally nonzero.
 *
 * Previous slab:
 *   -(E_prev/dt) tensor M
 *
 * Only alpha = 0, beta = nt - 1 is nonzero.
 */
static int temporal_same_nonzero(int degree, int alpha, int beta)
{
    const int nt = degree + 1;
    return alpha >= 0 && alpha < nt && beta >= 0 && beta < nt;
}

static int temporal_prev_nonzero(int degree, int alpha, int beta)
{
    const int nt = degree + 1;
    return alpha == 0 && beta == nt - 1;
}

static int read_num_nodes(const char* node_path)
{
    FILE* fp = fopen(node_path, "r");
    if(fp == NULL) {
        fail_path("cannot open node file", node_path);
    }

    int num_nodes = 0;
    if(fscanf(fp, "%d", &num_nodes) != 1 || num_nodes <= 0) {
        fclose(fp);
        fail_path("invalid node.dat header", node_path);
    }

    fclose(fp);
    return num_nodes;
}

static void read_spatial_elem(
    SPATIAL_MESH* mesh,
    const char*   elem_path,
    int           num_nodes)
{
    FILE* fp = fopen(elem_path, "r");
    if(fp == NULL) {
        fail_path("cannot open spatial elem file", elem_path);
    }

    int num_elems = 0;
    int num_local_nodes = 0;

    if(fscanf(fp, "%d %d", &num_elems, &num_local_nodes) != 2 ||
       num_elems <= 0 || num_local_nodes <= 0) {
        fclose(fp);
        fail_path("invalid spatial elem.dat header", elem_path);
    }

    const size_t num_conn =
        (size_t)num_elems * (size_t)num_local_nodes;

    if(num_local_nodes != 0 &&
       num_conn / (size_t)num_local_nodes != (size_t)num_elems) {
        fclose(fp);
        fail("spatial connectivity size overflow");
    }

    int* conn = (int*)checked_calloc(num_conn, sizeof(int));

    for(int e = 0; e < num_elems; e++) {
        for(int a = 0; a < num_local_nodes; a++) {
            int node = -1;
            if(fscanf(fp, "%d", &node) != 1) {
                fclose(fp);
                free(conn);
                fail_path("failed to read spatial elem.dat", elem_path);
            }

            if(node < 0 || node >= num_nodes) {
                fprintf(
                    stderr,
                    "ERROR: spatial element node ID out of range. "
                    "elem=%d local=%d node=%d num_nodes=%d\n",
                    e,
                    a,
                    node,
                    num_nodes);
                fclose(fp);
                free(conn);
                exit(EXIT_FAILURE);
            }

            conn[(size_t)e * (size_t)num_local_nodes + (size_t)a] = node;
        }
    }

    fclose(fp);

    mesh->num_nodes = num_nodes;
    mesh->num_elems = num_elems;
    mesh->num_local_nodes = num_local_nodes;
    mesh->conn = conn;
    mesh->adjacency = NULL;
}

static void build_spatial_graph(SPATIAL_MESH* mesh)
{
    mesh->adjacency = (INT_VECTOR*)checked_calloc(
        (size_t)mesh->num_nodes,
        sizeof(INT_VECTOR));

    /*
     * P1/Q1 finite-element mass and diffusion matrices couple every pair
     * of nodes belonging to the same element. Build the node graph as the
     * union of element cliques. Self entries are included explicitly.
     */
    for(int e = 0; e < mesh->num_elems; e++) {
        const int* elem_conn =
            &mesh->conn[(size_t)e * (size_t)mesh->num_local_nodes];

        for(int a = 0; a < mesh->num_local_nodes; a++) {
            const int row_node = elem_conn[a];

            for(int b = 0; b < mesh->num_local_nodes; b++) {
                int_vector_push(
                    &mesh->adjacency[row_node],
                    elem_conn[b]);
            }
        }
    }

    for(int i = 0; i < mesh->num_nodes; i++) {
        int_vector_sort_unique(&mesh->adjacency[i]);

        if(mesh->adjacency[i].size == 0) {
            fprintf(stderr, "ERROR: isolated spatial node %d\n", i);
            exit(EXIT_FAILURE);
        }
    }
}

static void free_spatial_mesh(SPATIAL_MESH* mesh)
{
    if(mesh->adjacency != NULL) {
        for(int i = 0; i < mesh->num_nodes; i++) {
            free(mesh->adjacency[i].data);
        }
        free(mesh->adjacency);
    }

    free(mesh->conn);
    memset(mesh, 0, sizeof(*mesh));
}

static void make_directory(const char* output_dir)
{
    if(mkdir(output_dir, 0777) != 0 && errno != EEXIST) {
        fail_path("cannot create output directory", output_dir);
    }
}

static void make_path(
    char*       path,
    size_t      path_size,
    const char* directory,
    const char* filename)
{
    const int n = snprintf(path, path_size, "%s/%s", directory, filename);
    if(n < 0 || (size_t)n >= path_size) {
        fail("output path is too long");
    }
}

static FILE* open_output(const char* output_dir, const char* filename)
{
    char path[PATH_MAX];
    make_path(path, sizeof(path), output_dir, filename);

    FILE* fp = fopen(path, "w");
    if(fp == NULL) {
        fail_path("cannot open output file", path);
    }

    /* Large buffered output helps for graph files. */
    (void)setvbuf(fp, NULL, _IOFBF, 1U << 20);
    return fp;
}

static void write_spatial_graph(
    const SPATIAL_MESH* mesh,
    const char* output_dir)
{
    FILE* fp = open_output(output_dir, "spatial_graph.dat");

    fprintf(fp, "%d\n", mesh->num_nodes);

    for(int i = 0; i < mesh->num_nodes; i++) {
        const INT_VECTOR* row = &mesh->adjacency[i];
        fprintf(fp, "%d %zu", i, row->size);

        for(size_t p = 0; p < row->size; p++) {
            fprintf(fp, " %d", row->data[p]);
        }
        fputc('\n', fp);
    }

    fclose(fp);
}

static int count_same_targets(int degree, int alpha, int symmetric)
{
    const int nt = degree + 1;
    int count = 0;

    for(int beta = 0; beta < nt; beta++) {
        const int nz = symmetric
            ? (temporal_same_nonzero(degree, alpha, beta) ||
               temporal_same_nonzero(degree, beta, alpha))
            : temporal_same_nonzero(degree, alpha, beta);

        if(nz) {
            count++;
        }
    }
    return count;
}

static int count_backward_targets(int degree, int alpha)
{
    const int nt = degree + 1;
    int count = 0;

    for(int beta = 0; beta < nt; beta++) {
        if(temporal_prev_nonzero(degree, alpha, beta)) {
            count++;
        }
    }
    return count;
}

static int count_forward_targets(int degree, int alpha)
{
    const int nt = degree + 1;
    int count = 0;

    for(int beta = 0; beta < nt; beta++) {
        if(temporal_prev_nonzero(degree, beta, alpha)) {
            count++;
        }
    }
    return count;
}

static void print_spatial_target_block(
    FILE*               fp,
    const INT_VECTOR*   spatial_row,
    int                 target_slab,
    int                 target_alpha,
    int                 nt,
    int                 num_spatial_nodes)
{
    for(size_t p = 0; p < spatial_row->size; p++) {
        const int64_t col = st_gid(
            target_slab,
            target_alpha,
            spatial_row->data[p],
            nt,
            num_spatial_nodes);

        fprintf(fp, " %" PRId64, col);
    }
}

static void write_st_graph(
    const SPATIAL_MESH* mesh,
    int degree,
    int num_slabs,
    int symmetric,
    const char* output_dir,
    const char* filename)
{
    const int nt = degree + 1;

    const int64_t num_st_nodes = checked_mul_i64(
        checked_mul_i64(num_slabs, nt, "number of temporal replicas"),
        mesh->num_nodes,
        "number of ST nodes");

    FILE* fp = open_output(output_dir, filename);
    fprintf(fp, "%" PRId64 "\n", num_st_nodes);

    for(int slab = 0; slab < num_slabs; slab++) {
        for(int alpha = 0; alpha < nt; alpha++) {
            const int same_targets =
                count_same_targets(degree, alpha, symmetric);

            const int backward_targets =
                (slab > 0)
                ? count_backward_targets(degree, alpha)
                : 0;

            const int forward_targets =
                (symmetric && slab + 1 < num_slabs)
                ? count_forward_targets(degree, alpha)
                : 0;

            const int temporal_target_count =
                same_targets + backward_targets + forward_targets;

            for(int i = 0; i < mesh->num_nodes; i++) {
                const INT_VECTOR* spatial_row = &mesh->adjacency[i];

                const int64_t row = st_gid(
                    slab,
                    alpha,
                    i,
                    nt,
                    mesh->num_nodes);

                const int64_t row_size = checked_mul_i64(
                    temporal_target_count,
                    (int64_t)spatial_row->size,
                    "ST graph row size");

                fprintf(fp, "%" PRId64 " %" PRId64, row, row_size);

                /* Same-slab blocks. */
                for(int beta = 0; beta < nt; beta++) {
                    const int nz = symmetric
                        ? (temporal_same_nonzero(degree, alpha, beta) ||
                           temporal_same_nonzero(degree, beta, alpha))
                        : temporal_same_nonzero(degree, alpha, beta);

                    if(!nz) {
                        continue;
                    }

                    print_spatial_target_block(
                        fp,
                        spatial_row,
                        slab,
                        beta,
                        nt,
                        mesh->num_nodes);
                }

                /* Exact backward dG jump block. */
                if(slab > 0) {
                    for(int beta = 0; beta < nt; beta++) {
                        if(!temporal_prev_nonzero(degree, alpha, beta)) {
                            continue;
                        }

                        print_spatial_target_block(
                            fp,
                            spatial_row,
                            slab - 1,
                            beta,
                            nt,
                            mesh->num_nodes);
                    }
                }

                /* Transpose edge required by the symmetric partition graph. */
                if(symmetric && slab + 1 < num_slabs) {
                    for(int beta = 0; beta < nt; beta++) {
                        if(!temporal_prev_nonzero(degree, beta, alpha)) {
                            continue;
                        }

                        print_spatial_target_block(
                            fp,
                            spatial_row,
                            slab + 1,
                            beta,
                            nt,
                            mesh->num_nodes);
                    }
                }

                fputc('\n', fp);
            }
        }
    }

    fclose(fp);
}

static void write_st_elem_and_maps(
    const SPATIAL_MESH* mesh,
    int degree,
    int num_slabs,
    const char* output_dir)
{
    const int nt = degree + 1;

    const int64_t num_st_cells = checked_mul_i64(
        num_slabs,
        mesh->num_elems,
        "number of ST cells");

    const int64_t num_nodes_per_st_cell = checked_mul_i64(
        nt,
        mesh->num_local_nodes,
        "number of nodes per ST cell");

	FILE* elem_fp = open_output(output_dir, "st_elem.dat");
	FILE* map_fp  = open_output(output_dir, "st_elem_map.dat");

    /* Standard connectivity file: no metadata columns. */
    fprintf(
        elem_fp,
        "%" PRId64 " %" PRId64 "\n",
        num_st_cells,
        num_nodes_per_st_cell);

    /* Custom metadata file. */
    fprintf(
        map_fp,
        "%" PRId64 " %d\n",
        num_st_cells,
        mesh->num_local_nodes);

    for(int slab = 0; slab < num_slabs; slab++) {
        for(int e = 0; e < mesh->num_elems; e++) {
            const int64_t cell_id =
                (int64_t)slab * (int64_t)mesh->num_elems + (int64_t)e;

            const int* elem_conn =
                &mesh->conn[(size_t)e * (size_t)mesh->num_local_nodes];

            /*
             * ST cell connectivity used by graph/element partition tools.
             * It contains every temporal DOF on the current slab.
             */
            for(int alpha = 0; alpha < nt; alpha++) {
                for(int a = 0; a < mesh->num_local_nodes; a++) {
                    const int64_t gid = st_gid(
                        slab,
                        alpha,
                        elem_conn[a],
                        nt,
                        mesh->num_nodes);

                    if(alpha != 0 || a != 0) {
                        fputc(' ', elem_fp);
                    }
                    fprintf(elem_fp, "%" PRId64, gid);
                }
            }
            fputc('\n', elem_fp);

            /*
             * cell_id slab_id global_spatial_element_id spatial_node_ids...
             */
            /*
             * Use the same schema as the partitioned reader.
             * In the unpartitioned file local_cell_id == global_cell_id.
             *
             *   local_cell_id global_st_cell_id slab_id
             *   global_space_element_id global_space_node_ids...
             */
            fprintf(
                map_fp,
                "%" PRId64 " %" PRId64 " %d %d",
                cell_id,
                cell_id,
                slab,
                e);

            for(int a = 0; a < mesh->num_local_nodes; a++) {
                fprintf(map_fp, " %d", elem_conn[a]);
            }
            fputc('\n', map_fp);
        }
    }

    fclose(elem_fp);
    fclose(map_fp);
}

static void write_st_node_map(
    const SPATIAL_MESH* mesh,
    int degree,
    int num_slabs,
    const char* output_dir)
{
    const int nt = degree + 1;
    const int64_t num_st_nodes = checked_mul_i64(
        checked_mul_i64(num_slabs, nt, "number of temporal replicas"),
        mesh->num_nodes,
        "number of ST nodes");

    FILE* fp = open_output(output_dir, "st_node_map.dat");

    /*
     * Use the same schema as the partitioned reader.
     * In the unpartitioned file all ST nodes are internal and
     * local_st_id == global_st_id.
     *
     *   num_internal_st_nodes num_total_st_nodes
     *   local_st_id global_st_id slab_id alpha global_space_node_id
     */
    fprintf(
        fp,
        "%" PRId64 " %" PRId64 "\n",
        num_st_nodes,
        num_st_nodes);

    for(int slab = 0; slab < num_slabs; slab++) {
        for(int alpha = 0; alpha < nt; alpha++) {
            for(int i = 0; i < mesh->num_nodes; i++) {
                const int64_t gid = st_gid(
                    slab,
                    alpha,
                    i,
                    nt,
                    mesh->num_nodes);

                fprintf(
                    fp,
                    "%" PRId64 " %" PRId64 " %d %d %d\n",
                    gid,
                    gid,
                    slab,
                    alpha,
                    i);
            }
        }
    }

    fclose(fp);
}


static void write_unpartitioned_space_maps(
    const SPATIAL_MESH* mesh,
    const char* output_dir)
{
    FILE* node_fp = open_output(output_dir, "space_node_map.dat");
    FILE* elem_fp = open_output(output_dir, "space_elem_map.dat");

    /*
     * Unpartitioned spatial maps.  A partitioning postprocessor must
     * rewrite these files for each MPI rank.
     */
    fprintf(
        node_fp,
        "%d %d\n",
        mesh->num_nodes,
        mesh->num_nodes);

    for(int i = 0; i < mesh->num_nodes; i++) {
        fprintf(node_fp, "%d %d\n", i, i);
    }

    fprintf(elem_fp, "%d\n", mesh->num_elems);

    for(int e = 0; e < mesh->num_elems; e++) {
        fprintf(elem_fp, "%d %d\n", e, e);
    }

    fclose(node_fp);
    fclose(elem_fp);
}

static void write_metadata(
    const SPATIAL_MESH* mesh,
    int degree,
    int num_slabs,
    const char* node_path,
    const char* elem_path,
    const char* output_dir)
{
    const int nt = degree + 1;
    const int64_t num_st_nodes = checked_mul_i64(
        checked_mul_i64(num_slabs, nt, "number of temporal replicas"),
        mesh->num_nodes,
        "number of ST nodes");

    const int64_t num_st_cells = checked_mul_i64(
        num_slabs,
        mesh->num_elems,
        "number of ST cells");

    FILE* fp = open_output(output_dir, "st_meta.dat");

    fprintf(fp, "version 1\n");
    fprintf(fp, "index_base 0\n");
    fprintf(fp, "numbering time_major\n");
    fprintf(fp, "gid_formula ((slab*nt+alpha)*Nx)+space_id\n");
    fprintf(fp, "degree %d\n", degree);
    fprintf(fp, "num_time_dofs %d\n", nt);
    fprintf(fp, "num_window_slabs %d\n", num_slabs);
    fprintf(fp, "num_global_space_nodes %d\n", mesh->num_nodes);
    fprintf(fp, "num_global_space_elements %d\n", mesh->num_elems);
    fprintf(fp, "num_space_nodes_per_element %d\n", mesh->num_local_nodes);
    fprintf(fp, "num_global_st_nodes %" PRId64 "\n", num_st_nodes);
    fprintf(fp, "num_global_st_cells %" PRId64 "\n", num_st_cells);
    fprintf(fp, "num_st_nodes_per_cell %d\n", nt * mesh->num_local_nodes);
    fprintf(fp, "input_node_file %s\n", node_path);
    fprintf(fp, "input_spatial_elem_file %s\n", elem_path);
    fprintf(fp, "partition_graph graph.dat\n");
    fprintf(fp, "matrix_graph matrix_graph.dat\n");
    fprintf(fp, "st_connectivity st_elem.dat\n");
    fprintf(fp, "st_cell_metadata st_elem_map.dat\n");
    fprintf(fp, "st_node_metadata st_node_map.dat\n");

    fclose(fp);
}

static void check_downstream_int_range(
    const SPATIAL_MESH* mesh,
    int degree,
    int num_slabs)
{
    const int nt = degree + 1;
    const int64_t num_st_nodes = checked_mul_i64(
        checked_mul_i64(num_slabs, nt, "number of temporal replicas"),
        mesh->num_nodes,
        "number of ST nodes");

    const int64_t num_st_cells = checked_mul_i64(
        num_slabs,
        mesh->num_elems,
        "number of ST cells");

    /*
     * The current solver and many MONOLIS/GEDATSU interfaces use int IDs.
     * Stop early instead of generating files that those readers cannot use.
     */
    if(num_st_nodes > INT_MAX) {
        fprintf(
            stderr,
            "ERROR: num_global_st_nodes=%" PRId64
            " exceeds INT_MAX=%d.\n",
            num_st_nodes,
            INT_MAX);
        exit(EXIT_FAILURE);
    }

    if(num_st_cells > INT_MAX) {
        fprintf(
            stderr,
            "ERROR: num_global_st_cells=%" PRId64
            " exceeds INT_MAX=%d.\n",
            num_st_cells,
            INT_MAX);
        exit(EXIT_FAILURE);
    }
}

static void print_usage(const char* program)
{
    printf(
        "Usage:\n"
        "  %s DEGREE NUM_WINDOW_SLABS node.dat spatial_elem.dat OUTPUT_DIR\n"
        "\n"
        "Example:\n"
        "  %s 1 4 ./node.dat ./elem.dat ./st_global\n",
        program,
        program);
}

int main(int argc, char* argv[])
{
    if(argc != 6) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    const int degree = atoi(argv[1]);
    const int num_slabs = atoi(argv[2]);
    const char* node_path = argv[3];
    const char* elem_path = argv[4];
    const char* output_dir = argv[5];

    if(degree < 0 || degree > 2) {
        fail("DEGREE must be 0, 1, or 2");
    }

    if(num_slabs <= 0) {
        fail("NUM_WINDOW_SLABS must be positive");
    }

    SPATIAL_MESH mesh;
    memset(&mesh, 0, sizeof(mesh));

    mesh.num_nodes = read_num_nodes(node_path);
    read_spatial_elem(&mesh, elem_path, mesh.num_nodes);
    build_spatial_graph(&mesh);

    check_downstream_int_range(&mesh, degree, num_slabs);
    make_directory(output_dir);

    write_spatial_graph(&mesh, output_dir);

    /* Exact directed matrix pattern. */
    write_st_graph(
        &mesh,
        degree,
        num_slabs,
        0,
        output_dir,
        "matrix_graph.dat");

    /* Symmetric closure for partitioning and safe matrix allocation. */
    write_st_graph(
        &mesh,
        degree,
        num_slabs,
        1,
        output_dir,
        "graph.dat");

    write_st_elem_and_maps(&mesh, degree, num_slabs, output_dir);
    write_st_node_map(&mesh, degree, num_slabs, output_dir);
    write_unpartitioned_space_maps(&mesh, output_dir);
    write_metadata(
        &mesh,
        degree,
        num_slabs,
        node_path,
        elem_path,
        output_dir);

    const int nt = degree + 1;
    const int64_t num_st_nodes =
        (int64_t)num_slabs * nt * mesh.num_nodes;
    const int64_t num_st_cells =
        (int64_t)num_slabs * mesh.num_elems;

    printf("Generated global ST preprocessing files.\n");
    printf("  degree                    : %d\n", degree);
    printf("  temporal DOFs             : %d\n", nt);
    printf("  window slabs              : %d\n", num_slabs);
    printf("  spatial nodes             : %d\n", mesh.num_nodes);
    printf("  spatial elements          : %d\n", mesh.num_elems);
    printf("  spatial nodes/element     : %d\n", mesh.num_local_nodes);
    printf("  global ST nodes           : %" PRId64 "\n", num_st_nodes);
    printf("  global ST cells           : %" PRId64 "\n", num_st_cells);
    printf("  output directory          : %s\n", output_dir);
    printf("  partition graph           : %s/graph.dat\n", output_dir);
    printf("  exact matrix graph        : %s/matrix_graph.dat\n", output_dir);
    printf("  ST connectivity           : %s/st_elem.dat\n", output_dir);
    printf("  ST cell metadata          : %s/st_elem_map.dat\n", output_dir);
    printf("  ST node metadata          : %s/st_node_map.dat\n", output_dir);

    free_spatial_mesh(&mesh);
    return EXIT_SUCCESS;
}
