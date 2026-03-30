
#include "read_mpi.h"

int read_total_num_D_bc(
    const int       Ddof,
    int*            mydomain,
    const char*     directory)
{
    char fname_D_bc_in[BUFFER_SIZE];
    FILE* fp;
    int sum = 0;

    int num_D_bcs;    int block_size;
    for(int m = 0; m < Ddof; m++){
        snprintf(fname_D_bc_in, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_D_BC, mydomain[m]);
        fp = BBFE_sys_read_fopen(fp, fname_D_bc_in, directory);

        fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

        sum += num_D_bcs;
        fclose(fp);
    }
    return sum;
}

int read_total_num_D_bc_p(
    const int       Ddof,
    int*            mydomain,
    const char*     label,
    const char*     directory)
{
    char fname_D_bc_in[BUFFER_SIZE];
    FILE* fp;
    int sum = 0;

    int num_D_bcs;    int block_size;
    for(int m = 0; m < Ddof; m++){
        snprintf(fname_D_bc_in, BUFFER_SIZE, "parted.0/%s.dat.%d", label, mydomain[m]);
        fp = BBFE_sys_read_fopen(fp, fname_D_bc_in, directory);

        fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

        sum += num_D_bcs;
        fclose(fp);
    }
    return sum;
}

int read_total_num_D_bc_v(
    const int Ddof,
    int* mydomain,
    const char*     directory)
{
    char fname_D_bc_in[BUFFER_SIZE];
    FILE* fp;
    int sum = 0;

    int num_D_bcs;    int block_size;
    for(int m = 0; m < Ddof; m++){
        snprintf(fname_D_bc_in, BUFFER_SIZE, "parted.0/D_bc_v.dat.%d", mydomain[m]);
        fp = BBFE_sys_read_fopen(fp, fname_D_bc_in, directory);

        fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

        sum += num_D_bcs;
        fclose(fp);
    }
    return sum;
}