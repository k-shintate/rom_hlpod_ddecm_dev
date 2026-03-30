
#include "read_std.h"

static const int BUFFER_SIZE = 10000;
static const char* OUTPUT_FILENAME_SVD          = "singularvalue.dat";


int ROM_std_hlpod_read_n_internal(
    FILE*        fp,
    const char*  fname,
    const char*  directory)
{
    char id[BUFFER_SIZE];
    int n_internal_vertex;
    int ndof;
    fp = ROM_BB_read_fopen(fp, fname, directory);
    fscanf(fp, "%s %d", id, &(ndof));
    fscanf(fp, "%d", &(n_internal_vertex));
    fclose(fp);

    return n_internal_vertex;
}

void ROM_std_hlpod_read_node_id(
    FILE*        fp,
    int*         node_id,
    const int    n_internal_vertex,
    const char*  fname,
    const char*  directory)
{
    char id[BUFFER_SIZE];
    int total_num_nodes;
    int ndof;

    fp = ROM_BB_read_fopen(fp, fname, directory);
    fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes),&(ndof));
    for(int i = 0; i < n_internal_vertex; i++) {
        fscanf(fp, "%d", &(node_id[i]));
    }
    fclose(fp);
}

void ROM_std_hlpod_read_singular_values(
    double*                     singular_value,
    const int           num_modes,
    const int           subdomain_id,
    const char*         label,
    const char*         directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE, "%s/subdomain%d/%s", label, subdomain_id, OUTPUT_FILENAME_SVD);

    fp = ROM_BB_read_fopen(fp, fname, directory);
    for(int i = 0; i < num_modes; i++){
        fscanf(fp,"%lf", &(singular_value[i]));
        //printf("%lf\n", singular_value[i]);
    }
    fclose(fp);
}


int ROM_std_hlpod_read_num_modes(
    const int		subdomain_id,
    const char*		label,
    const char*		directory)
{
    int num_modes;
    char fname[BUFFER_SIZE];
    snprintf(fname, BUFFER_SIZE,"%s/num_modes.dat.%d",label ,subdomain_id);

    FILE* fp;

    fp = ROM_BB_read_fopen(fp, fname, directory);
    fscanf(fp,"%d", &(num_modes));
    fclose(fp);

    return num_modes;
}

void ROM_std_hlpod_read_pod_modes_node(
    FILE*        fp,
    double**     S,
    const int    dof,
    const int    n_internal_vertex_in,
    const int    num_modes,
    const int    subdomain_id,
    const char*  label,
    const char*  directory)
{
    char id[BUFFER_SIZE];
    char fname[BUFFER_SIZE];
    int block_size;
    int n_internal_vertex;

    for(int j = 0; j < num_modes; j++){
        snprintf(fname, BUFFER_SIZE,"%s/subdomain%d/lb_pod_modes%d.dat", label, subdomain_id, j);

        fp = ROM_BB_read_fopen(fp, fname, directory);

        fscanf(fp,"%s", id);
        fscanf(fp,"%d %d", &n_internal_vertex, &block_size);

        if(n_internal_vertex_in != n_internal_vertex){
            printf("Error: n_internal_vertex values differ between the parted file (parted.0) = %d and the pod_modes file = %d",
                    n_internal_vertex_in, n_internal_vertex);
            exit(1);
        }

        if(dof == 1){
            for(int i = 0; i < n_internal_vertex; i++){
                fscanf(fp,"%lf", &(S[i][j]));
            }
        }
        else if(dof == 3){
            for(int i = 0; i < n_internal_vertex; i++){
                fscanf(fp,"%lf %lf %lf", &(S[i*block_size + 0][j]), &(S[i*block_size + 1][j]), &(S[i*block_size + 2][j]));
            }
        }
        else{
            printf("\n Error: dof is not 1 or 3 \n");
            exit(1);
        }
        fclose(fp);
    }
}


