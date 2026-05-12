#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"
#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* CODENAME               = "util/vtkddm/vtkddm_hex >";
static const char* VOIDNAME               = "                        ";
static const int   BUFFER_SIZE            = 10000;
static const int   FEM_LOCAL_NUM_NODE_HEX = 8;

static const char* OUTPUT_FILENAME_VTK_LABEL  = "num_parted";

static const char* DEF_DIRECTORY            = ".";
static const char* DEF_FILENAME_NODE        = "node.dat";
static const char* DEF_FILENAME_ELEM        = "elem.dat";
static const char* DEF_FILENAME_PART_DIR    = "parted.0";
static const char* DEF_FILENAME_MPI_ID_CONN = "elem.dat.id";
static const char* DEF_FILENAME_VTK         = "vtkddm_hex.vtk";

static const char* OPTION_INFILE_NODE        = "-in";
static const char* OPTION_INFILE_ELEM        = "-ie";
static const char* OPTION_INFILE_PART_DIR    = "-id";
static const char* OPTION_INFILE_MPI_ID_CONN = "-im";
static const char* OPTION_OUTFILE_VTK        = "-ov";

static const char* OPTION_DIRECTORY = "-d";

typedef struct
{
    int         total_num_nodes;
    int         total_num_elems;
    double**    node;
    int**       conn;
} FE_DATA;

typedef struct
{
    int     num_parallel;
    int*    connectivity_id;
} MPI_PARALLEL;

typedef struct
{
    const char* infile_node;
    const char* infile_elem;

    const char* infile_part_dir;
    const char* infile_mpiid_conn;
    const char* outfile_vtk;
    const char* directory;
} SETTINGS;


/* -------------------------------------------------------------------------- */
/* Local helper functions                                                      */
/* -------------------------------------------------------------------------- */

/*
    node.dat format:
        num_nodes
        x y z
        x y z
        ...

    elem.dat format:
        num_elems num_local_nodes
        n1 n2 n3 n4 n5 n6 n7 n8
        ...
*/

static int BB_std_read_file_pointer_get_num_nodes(
    FILE*       fp,
    const int   buffer_size)
{
    (void)buffer_size;

    int num_points = 0;

    rewind(fp);

    if(fscanf(fp, "%d", &num_points) != 1) {
        fprintf(stderr,
                "%s Failed to read node file header. Expected first line: num_nodes.\n",
                CODENAME);
        exit(EXIT_FAILURE);
    }

    if(num_points <= 0) {
        fprintf(stderr,
                "%s Invalid node file header: num_nodes = %d.\n",
                CODENAME, num_points);
        exit(EXIT_FAILURE);
    }

    return num_points;
}


static int BB_std_read_file_pointer_get_num_elems(
    FILE*       fp,
    const int   num_local_nodes,
    const int   buffer_size)
{
    (void)buffer_size;

    int num_elems = 0;
    int num_vals  = 0;

    rewind(fp);

    if(fscanf(fp, "%d %d", &num_elems, &num_vals) != 2) {
        fprintf(stderr,
                "%s Failed to read element file header. Expected first line: num_elems num_local_nodes.\n",
                CODENAME);
        exit(EXIT_FAILURE);
    }

    if(num_elems <= 0 || num_vals != num_local_nodes) {
        fprintf(stderr,
                "%s Invalid element file header: num_elems = %d, num_vals = %d, expected num_vals = %d.\n",
                CODENAME, num_elems, num_vals, num_local_nodes);
        exit(EXIT_FAILURE);
    }

    return num_elems;
}


static void BB_std_read_file_pointer_get_val_2d_double(
    FILE*       fp,
    const int   row,
    const int   col,
    double**    val)
{
    int num_points = 0;

    rewind(fp);

    if(fscanf(fp, "%d", &num_points) != 1) {
        fprintf(stderr,
                "%s Failed to read node file header.\n",
                CODENAME);
        exit(EXIT_FAILURE);
    }

    if(num_points != row) {
        fprintf(stderr,
                "%s Node header mismatch: header num_nodes = %d, requested row = %d.\n",
                CODENAME, num_points, row);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            if(fscanf(fp, "%lf", &(val[i][j])) != 1) {
                fprintf(stderr,
                        "%s Failed to read node value at row %d, col %d.\n",
                        CODENAME, i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}


static void BB_std_read_file_pointer_get_val_2d_int(
    FILE*       fp,
    const int   row,
    const int   col,
    int**       val)
{
    int num_elems = 0;
    int num_vals  = 0;

    rewind(fp);

    if(fscanf(fp, "%d %d", &num_elems, &num_vals) != 2) {
        fprintf(stderr,
                "%s Failed to read element file header.\n",
                CODENAME);
        exit(EXIT_FAILURE);
    }

    if(num_elems != row || num_vals != col) {
        fprintf(stderr,
                "%s Element header mismatch: header = %d x %d, requested = %d x %d.\n",
                CODENAME, num_elems, num_vals, row, col);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < row; i++) {
        for(int j = 0; j < col; j++) {
            if(fscanf(fp, "%d", &(val[i][j])) != 1) {
                fprintf(stderr,
                        "%s Failed to read element value at row %d, col %d.\n",
                        CODENAME, i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
}


static void BB_vtk_write_elem_vals_scalar(
    FILE*       fp,
    double*     vals,
    const int   total_num_elems,
    const char* label)
{
    fprintf(fp, "SCALARS %s double 1\n", label);
    fprintf(fp, "LOOKUP_TABLE default\n");

    for(int i = 0; i < total_num_elems; i++) {
        fprintf(fp, "%.15e\n", vals[i]);
    }
}


/* -------------------------------------------------------------------------- */

int args_manager(
    SETTINGS*   set,
    int         argc,
    char*       argv[])
{
    if(argc < 2) {
        printf("%s Please input parameters.\n", CODENAME);
        printf("%s Format: \n", CODENAME);
        printf("%s    %s [the num. of parallel]\n",
                VOIDNAME, argv[0]);

        printf("%s Options: \n", CODENAME);
        printf("%s     %s [input/output directory]\n",
                VOIDNAME, OPTION_DIRECTORY);
        printf("%s     %s [input filename (nodes)]\n",
                VOIDNAME, OPTION_INFILE_NODE);
        printf("%s     %s [input filename (elements)]\n",
                VOIDNAME, OPTION_INFILE_ELEM);
        printf("%s     %s [input filename (parted directory name)]\n",
                VOIDNAME, OPTION_INFILE_PART_DIR);
        printf("%s     %s [input filename (elem id)]\n",
                VOIDNAME, OPTION_INFILE_MPI_ID_CONN);
        printf("%s     %s [output filename (vtk)]\n",
                VOIDNAME, OPTION_OUTFILE_VTK);

        printf("\n");
        exit(0);
    }

    int num_parallel = atoi(argv[1]);

    int num;

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_DIRECTORY);
    set->directory = (num == -1) ? DEF_DIRECTORY : argv[num + 1];

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_INFILE_NODE);
    set->infile_node = (num == -1) ? DEF_FILENAME_NODE : argv[num + 1];

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_INFILE_ELEM);
    set->infile_elem = (num == -1) ? DEF_FILENAME_ELEM : argv[num + 1];

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_INFILE_PART_DIR);
    set->infile_part_dir = (num == -1) ? DEF_FILENAME_PART_DIR : argv[num + 1];

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_INFILE_MPI_ID_CONN);
    set->infile_mpiid_conn = (num == -1) ? DEF_FILENAME_MPI_ID_CONN : argv[num + 1];

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_OUTFILE_VTK);
    set->outfile_vtk = (num == -1) ? DEF_FILENAME_VTK : argv[num + 1];

    printf("%s The num. of parallel              : %d\n", CODENAME, num_parallel);
    printf("%s Input/output directory            : %s\n", CODENAME, set->directory);
    printf("%s Input filename (nodes)            : %s\n", CODENAME, set->infile_node);
    printf("%s Input filename (elements)         : %s\n", CODENAME, set->infile_elem);
    printf("%s Input filename (parted directory) : %s\n", CODENAME, set->infile_part_dir);
    printf("%s Input filename (elem id)          : %s\n", CODENAME, set->infile_mpiid_conn);
    printf("%s Output filename (vtk)             : %s\n", CODENAME, set->outfile_vtk);

    return num_parallel;
}


void BB_vtk_cell_data(
    const int       cell_type,
    const int       num_nodes,
    const int       total_num_elems,
    const int       num_local_nodes,
    double**        node,
    int**           conn,
    const char*     filename,
    const char*     directory,
    double*         cell_vals)
{
    FILE* fp;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    BB_vtk_write_header(fp);
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    BB_vtk_write_points_3d(fp, num_nodes, node);
    BB_vtk_write_cells(fp, total_num_elems, num_local_nodes, conn);
    BB_vtk_write_cell_types(fp, total_num_elems, cell_type);

    fprintf(fp, "CELL_DATA %d\n", total_num_elems);
    BB_vtk_write_elem_vals_scalar(
        fp,
        cell_vals,
        total_num_elems,
        OUTPUT_FILENAME_VTK_LABEL);

    fclose(fp);
}


void read_FEM_mesh_files(
    FE_DATA*    fe,
    SETTINGS*   set)
{
    char f_name[BUFFER_SIZE];
    FILE* fp;

    /* Read node coordinates */
    snprintf(f_name, BUFFER_SIZE, "%s", set->infile_node);
    fp = BBFE_sys_read_fopen(fp, f_name, set->directory);

    fe->total_num_nodes =
        BB_std_read_file_pointer_get_num_nodes(fp, BUFFER_SIZE);

    fe->node =
        BB_std_calloc_2d_double(fe->node, fe->total_num_nodes, 3);

    BB_std_read_file_pointer_get_val_2d_double(
        fp,
        fe->total_num_nodes,
        3,
        fe->node);

    fclose(fp);

    /* Read element connectivity */
    snprintf(f_name, BUFFER_SIZE, "%s", set->infile_elem);
    fp = BBFE_sys_read_fopen(fp, f_name, set->directory);

    fe->total_num_elems =
        BB_std_read_file_pointer_get_num_elems(
            fp,
            FEM_LOCAL_NUM_NODE_HEX,
            BUFFER_SIZE);

    fe->conn =
        BB_std_calloc_2d_int(
            fe->conn,
            fe->total_num_elems,
            FEM_LOCAL_NUM_NODE_HEX);

    BB_std_read_file_pointer_get_val_2d_int(
        fp,
        fe->total_num_elems,
        FEM_LOCAL_NUM_NODE_HEX,
        fe->conn);

    fclose(fp);
}


/*
    Read:
        parted.x/elem.dat.id.rank

    Expected format:
        label
        num_points num_vals
        global_element_id
        global_element_id
        ...
*/
int read_mpi_id(
    MPI_PARALLEL*   mpi,
    const char*     filename,
    const char*     directory)
{
    FILE* fp;
    char buf[BUFFER_SIZE];
    int num_points = 0;
    int num_vals   = 0;

    fp = BBFE_sys_read_fopen(fp, filename, directory);

    BB_std_scan_line(&fp, BUFFER_SIZE, "%s", buf);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d %d", &num_points, &num_vals);

    if(num_points == 0) {
        printf("%s File \"%s\", Label \"%s\", num_points = 0.\n",
                CODENAME, filename, buf);
        exit(EXIT_FAILURE);
    }

    if(num_vals != 1) {
        printf("%s File \"%s\", Label \"%s\", num_vals != 1.\n",
                CODENAME, filename, buf);
        exit(EXIT_FAILURE);
    }

    free(mpi->connectivity_id);
    mpi->connectivity_id = NULL;

    mpi->connectivity_id =
        BB_std_calloc_1d_int(mpi->connectivity_id, num_points);

    for(int i = 0; i < num_points; i++) {
        BB_std_scan_line(
            &fp,
            BUFFER_SIZE,
            "%d",
            &(mpi->connectivity_id[i]));
    }

    fclose(fp);

    return num_points;
}


static void get_input_filename_mpi(
    char*       out,
    const int   out_size,
    const int   myrank,
    const char* filename_body)
{
    snprintf(out, out_size, "%s.%d", filename_body, myrank);
}


void build_num_parted(
    MPI_PARALLEL*   mpi,
    SETTINGS*       set,
    double*         e_global,
    const char*     directory,
    const int       total_num_elems)
{
    char filename_body[BUFFER_SIZE];
    char filename[BUFFER_SIZE];

    snprintf(
        filename_body,
        BUFFER_SIZE,
        "%s/%s",
        set->infile_part_dir,
        set->infile_mpiid_conn);

    for(int i = 0; i < mpi->num_parallel; i++) {
        get_input_filename_mpi(
            filename,
            BUFFER_SIZE,
            i,
            filename_body);

        int num_elem_par = read_mpi_id(mpi, filename, directory);

        for(int j = 0; j < num_elem_par; j++) {
            int g = mpi->connectivity_id[j];

            if(g < 0 || g >= total_num_elems) {
                fprintf(stderr,
                        "%s Invalid global element id %d in %s. total_num_elems = %d\n",
                        CODENAME, g, filename, total_num_elems);
                exit(EXIT_FAILURE);
            }

            if(e_global[g] == 0.0) {
                e_global[g] = (double)(i + 1);
            }
            else {
                e_global[g] = -1.0;
            }
        }
    }
}


int main(
    int     argc,
    char*   argv[])
{
    printf("\n");

    FE_DATA      fe;
    MPI_PARALLEL mpi;
    SETTINGS     set;

    fe.total_num_nodes = 0;
    fe.total_num_elems = 0;
    fe.node = NULL;
    fe.conn = NULL;

    mpi.num_parallel = 0;
    mpi.connectivity_id = NULL;

    BB_calc_void();

    mpi.num_parallel = args_manager(&set, argc, argv);

    if(mpi.num_parallel <= 0) {
        fprintf(stderr,
                "%s Invalid num_parallel = %d. Please pass positive number.\n",
                CODENAME, mpi.num_parallel);
        exit(EXIT_FAILURE);
    }

    read_FEM_mesh_files(&(fe), &(set));

    double* e_global = NULL;
    e_global =
        BB_std_calloc_1d_double(e_global, fe.total_num_elems);

    build_num_parted(
        &(mpi),
        &(set),
        e_global,
        set.directory,
        fe.total_num_elems);

    BB_vtk_cell_data(
        TYPE_VTK_HEXAHEDRON,
        fe.total_num_nodes,
        fe.total_num_elems,
        FEM_LOCAL_NUM_NODE_HEX,
        fe.node,
        fe.conn,
        set.outfile_vtk,
        set.directory,
        e_global);

    free(mpi.connectivity_id);
    free(e_global);

    return 0;
}