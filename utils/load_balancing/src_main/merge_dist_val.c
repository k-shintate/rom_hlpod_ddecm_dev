#include "monolis.h"

#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include "BBFE/std/integ.h"
#include "BBFE/std/shapefunc.h"
#include "BBFE/std/mapping.h"

#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "core.h"
#include "core_mpi.h"
#include "lb_dataset.h"
#include "graph.h"

static const char* OPTION_NUM_DD = "-nd";

typedef struct
{
    const char* directory;
    int num_subdomains;

} CONDITIONS;

typedef struct
{
    CONDITIONS   cond;

} FE_SYSTEM;

void merge_dist_val(
    const int   num_subdomains,
    const char* input_fname,
    const char* output_fname,
    const char* directory)
{
    enum { FIXED_DOF = 4 };

    char  path[BUFFER_SIZE];
    char  token[BUFFER_SIZE];
    FILE* fp = NULL;

    int total_num_nodes = 0;
    int ndof            = 0;

    snprintf(path, BUFFER_SIZE, "node.dat");
    fp = BBFE_sys_read_fopen(fp, path, directory);
    fscanf(fp, "%d", &total_num_nodes);
    fclose(fp);

    // 出力用ベクトル（0 初期化）
    double* vec = BB_std_calloc_1d_double(vec, total_num_nodes * FIXED_DOF);

    for (int m = 0; m < num_subdomains; ++m) {
        int n_internal = 0;
        int tmp_int    = 0;

        snprintf(path, BUFFER_SIZE, "parted.0/node.dat.n_internal.%d", m);
        fp = BBFE_sys_read_fopen(fp, path, directory);
        fscanf(fp, "%s %d", token, &tmp_int);
        fscanf(fp, "%d", &n_internal);
        fclose(fp);

        int* ids = BB_std_calloc_1d_int(ids, n_internal);

        snprintf(path, BUFFER_SIZE, "parted.0/node.dat.id.%d", m);
        fp = BBFE_sys_read_fopen(fp, path, directory);
        fscanf(fp, "%s", token);
        fscanf(fp, "%d %d", &tmp_int, &tmp_int);
        for (int i = 0; i < n_internal; ++i) {
            fscanf(fp, "%d", &ids[i]);
        }
        fclose(fp);

        snprintf(path, BUFFER_SIZE, "hot_start/%s.%d.dat", input_fname, m);
        fp = BBFE_sys_read_fopen(fp, path, directory);
        fscanf(fp, "%s", token);
        fscanf(fp, "%d %d", &tmp_int, &ndof);

        for (int i = 0; i < n_internal; ++i) {
            for (int j = 0; j < ndof; ++j) {
                double val = 0.0;
                fscanf(fp, "%lf", &val);
                vec[ids[i] * FIXED_DOF + j] = val;
            }
        }
        fclose(fp);

        BB_std_free_1d_int(ids, n_internal);
    }

    fp = BBFE_sys_write_fopen(fp, output_fname, directory);
    fprintf(fp, "initialization\n");
    fprintf(fp, "%d %d\n", total_num_nodes, ndof);
    for (int i = 0; i < total_num_nodes; ++i) {
        for (int j = 0; j < ndof; ++j) {
            fprintf(fp, "%e ", vec[i * ndof + j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void set_condition(
    int argc,
    char* argv[],
    CONDITIONS* cond)
{
    int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_DD);
    if(num == -1) {
        printf("\nargs error num_1st_subdomains");
        exit(1);
    }
    else {
        cond->num_subdomains = atoi(argv[num + 1]); //hddにおける1層目のメッシュ分割数 (POD計算領域数)
    }
    printf("\nnum_1st_dd = %d\n", cond->num_subdomains);
}

int main (
        int argc,
        char* argv[])
{
    printf("\n");
    BB_vtk_void();

    FE_SYSTEM sys;

    monolis_global_initialize();
    
    sys.cond.directory = get_directory_name(argc, argv, CODENAME);

    set_condition(argc, argv, &(sys.cond));

    merge_dist_val(sys.cond.num_subdomains, "velosity_pressure.200.000000", "hot_start.dat", sys.cond.directory);

    monolis_global_finalize();

    printf("\n");
    return 0;
}