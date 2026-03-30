#include "write_fe.h"

#include "write_BB.h"
#include "write_std.h"

static const char* OUTPUT_PODMODES_NUMBER     	= "pod_modes.%d";

static const int BUFFER_SIZE = 10000;

void ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
    const int 		n_internal_vertex,
    BBFE_DATA*		fe,
    double**     	pod_modes,
    const int 		num_modes,
    const char*		output_fname,
    const char*		directory)
{
    const char* fname;
    const char* pod_modes_name_vtk;

    char fname_vtk[BUFFER_SIZE];
    char pod_modes_number[BUFFER_SIZE];

    snprintf(fname_vtk, BUFFER_SIZE,"%s", output_fname);
    fname = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

    FILE* fp;
    fp = ROM_BB_write_fopen(fp, fname, directory);

    switch( fe->local_num_nodes ) {
        case 4:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
            break;

        case 8:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
            break;
    }
    
    double* pod_modes_out;
    pod_modes_out = BB_std_calloc_1d_double(pod_modes_out, fe->total_num_nodes);

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    for(int j = 0; j < num_modes; j++){		
        snprintf(pod_modes_number, BUFFER_SIZE, OUTPUT_PODMODES_NUMBER, j);
        pod_modes_name_vtk = pod_modes_number;
        for(int i = 0; i < fe->total_num_nodes; i++){
            pod_modes_out[i] = pod_modes[i][j];
        }
        BB_vtk_write_point_vals_scalar(fp, pod_modes_out, fe->total_num_nodes, pod_modes_name_vtk);
    }

    BB_std_free_1d_double(pod_modes_out, fe->total_num_nodes);

    fclose(fp);
}

void ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
    const int 		n_internal_vertex,
    BBFE_DATA*		fe,
    double**     	pod_modes,
    const int 		num_modes,
    const char*		output_fname,
    const char*		directory)
{
    const char* fname;
    const char* pod_modes_name_vtk;

    char fname_vtk[BUFFER_SIZE];
    char pod_modes_number[BUFFER_SIZE];

    snprintf(fname_vtk, BUFFER_SIZE,"%s", output_fname);
    fname = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

    FILE* fp;
    fp = ROM_BB_write_fopen(fp, fname, directory);

    switch( fe->local_num_nodes ) {
        case 4:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
            break;

        case 8:
            BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
            break;
    }
    
    double** pod_modes_out;
    pod_modes_out = BB_std_calloc_2d_double(pod_modes_out, fe->total_num_nodes, 3);

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    for(int j = 0; j < num_modes; j++){
        snprintf(pod_modes_number, BUFFER_SIZE, OUTPUT_PODMODES_NUMBER,j);
        pod_modes_name_vtk = pod_modes_number;
        for(int i = 0; i < fe->total_num_nodes; i++){
            for(int k = 0; k < 3; k++){
                pod_modes_out[i][k] = pod_modes[3*i + k][j];
            }
        }
        BB_vtk_write_point_vals_vector(fp, pod_modes_out, fe->total_num_nodes, pod_modes_name_vtk);
    }

    BB_std_free_2d_double(pod_modes_out, fe->total_num_nodes, 3);

    fclose(fp);
}


void ROM_sys_hlpod_fe_write_pod_modes(
        BBFE_DATA*		fe,
        double**		pod_modes,
        const int 		n_internal_vertex,
        const int 		num_modes,
        const int 		ndof,
        const char*     label,
        const char*		directory)
{
    double** vec;
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes*ndof, num_modes);

    for(int j = 0; j < num_modes; j++){
        for(int i = 0; i < fe->total_num_nodes; i++){
            for(int k = 0; k < ndof; k++){
                vec[i*ndof +k][j] = pod_modes[i*ndof + k][j];
            }
        }
    }

    if(ndof==1){
        ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
            fe->total_num_nodes,
            fe,
            vec,
            num_modes,
            label,
            directory);
    }
    else if(ndof==3){
        ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
            fe->total_num_nodes,
            fe,
            vec,
            num_modes,
            label,
            directory);
    }
    else{
        printf("\nError: ndof id not 1 or 3\n");
        exit(1);
    }

    BB_std_free_2d_double(vec, fe->total_num_nodes*ndof, num_modes);
}

void ROM_sys_hlpod_fe_write_local_pod_modes(
        BBFE_DATA*		fe,
        double**		pod_modes,
        const int 		n_internal_vertex,
        const int 		num_modes,
        const int		ndof,
        const char*     label,
        const char*		directory)
{
    double** vec;
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes*ndof, num_modes);

    for(int j = 0; j < num_modes; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            for(int k = 0; k < ndof; k++){
                vec[i*ndof + k][j] = pod_modes[i*ndof + k][j];
            }
        }
    }

    if(ndof==1){
        ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
            n_internal_vertex,
            fe,
            vec,
            num_modes,
            label,
            directory);
    }
    else if(ndof==3){
        ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
            n_internal_vertex,
            fe,
            vec,
            num_modes,
            label,
            directory);
    }
    else{
        printf("\nError: ndof id not 1 or 3\n");
        exit(1);
    }

    BB_std_free_2d_double(vec, fe->total_num_nodes*ndof, num_modes);
}

void ROM_sys_hlpod_fe_write_pod_modes_diag(
        BBFE_DATA*		fe,
        double**		pod_modes,
        const int 		n_internal_vertex,
        const int 		num_modes_v,
        const int 		num_modes_vi,
        const char*     label_v,
        const char*		label_p,
        const char*		directory)
{
    double** vec;
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            for(int k = 0; k < 3; k++){
                vec[i*3 +k][j] = pod_modes[i*4 + k][j];
            }
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_v,
        directory);

    BB_std_free_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            vec[i][j] = pod_modes[i*4 + 3][j + num_modes_v];
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_p,
        directory);

    BB_std_free_2d_double(vec, fe->total_num_nodes, num_modes_vi);
}

void ROM_sys_hlpod_fe_write_local_pod_modes_diag(
        BBFE_DATA*		fe,
        double**		pod_modes,
        const int 		n_internal_vertex,
        const int 		num_modes_v,
        const int 		num_modes_vi,
        const char*     label_v,
        const char*     label_p,
        const char*		directory)
{
    double** vec;
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            for(int k = 0; k < 3; k++){
                vec[i*3 + k][j] = pod_modes[i*4 + k][j];
            }
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_v,
        directory);
    
    BB_std_free_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);

    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            vec[i][j] = pod_modes[i*4 + 3][j + num_modes_v];
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_p,
        directory);

    BB_std_free_2d_double(vec, fe->total_num_nodes, num_modes_vi);
}

void ROM_sys_hlpod_fe_write_local_pod_modes_diag_id(
        BBFE_DATA*		fe,
        double**		pod_modes,
        int*			node_id_local,
        const int 		n_internal_vertex,
        const int 		num_modes_v,
        const int 		num_modes_vi,
        const char*     label_v,
        const char*     label_p,
        const char*		directory)
{
    double** vec;
    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            for(int k = 0; k < 3; k++){
                vec[node_id_local[i]*3 + k][j] = pod_modes[node_id_local[i]*4 + k][j];
            }
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_vector(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_v,
        directory);
    
    BB_std_free_2d_double(vec, fe->total_num_nodes*3, num_modes_vi);

    vec = BB_std_calloc_2d_double(vec, fe->total_num_nodes, num_modes_vi);

    for(int j = 0; j < num_modes_vi; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            vec[node_id_local[i]][j] = pod_modes[node_id_local[i]*4 + 3][j + num_modes_v];
        }
    }

    ROM_sys_hlpod_fe_write_pod_modes_vtk_scalar(
        n_internal_vertex,
        fe,
        vec,
        num_modes_vi,
        label_p,
        directory);

    BB_std_free_2d_double(vec, fe->total_num_nodes, num_modes_vi);
}