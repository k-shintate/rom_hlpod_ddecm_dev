
#include "write_mpi.h"
#include "graph.h"
#include "core_mpi.h"

void LB_write_node(
    const int ie,
    LPOD_LB*        lpod_lb,
    int* num_search,
    double** global_x_internal,
    double** global_x_overlap,
    int* overlap_id_for_sort,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];

    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.%d", INPUT_FILENAME_NODE, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"%d\n", lpod_lb->total_n_internal + num_search[lpod_lb->meta_domain_neib]);

    for(int i = 0; i < lpod_lb->total_n_internal; i++){
        int index = lpod_lb->node_id_internal_for_sort[i];   
        fprintf(fp,"%lf %lf %lf\n", 
                global_x_internal[index][0], global_x_internal[index][1], global_x_internal[index][2]);
    }

    for(int i = 0; i < ie; i++){
        int index = overlap_id_for_sort[i]- lpod_lb->total_n_internal;

        fprintf(fp,"%lf %lf %lf\n", 
                        global_x_overlap[index][0], global_x_overlap[index][1], global_x_overlap[index][2]);
    }
    fclose(fp);
}

void LB_write_node_internal(
    const int ndof, 
    const int total_n_internal,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.n_internal.%d", INPUT_FILENAME_NODE, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    fprintf(fp,"#n_internal %d\n", ndof);
    fprintf(fp,"%d", lpod_lb->total_n_internal);

    fclose(fp);
}


void LB_write_node_id(
    const int ndof, 
    const int ie,
    int* overlap_id,
    int* num_search,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.id.%d", INPUT_FILENAME_NODE, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"#id\n");

    fprintf(fp,"%d %d\n", lpod_lb->total_n_internal + num_search[lpod_lb->meta_domain_neib], ndof);


    for(int i = 0; i < lpod_lb->total_n_internal; i++){
        fprintf(fp,"%d\n",lpod_lb->node_id_internal[i]);
    }

    for(int i = 0; i < ie; i++){
        fprintf(fp,"%d\n",overlap_id[i]);
    }

    fclose(fp);
}


void LB_write_node_recv(
    int* recv_index,
    int* num_search,
    LPOD_LB*        lpod_lb,
    const char* directory)
{    
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.recv.%d", INPUT_FILENAME_NODE, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"%d %d\n", lpod_lb->meta_domain_neib, num_search[lpod_lb->meta_domain_neib]);

    for(int d = 0; d < lpod_lb->meta_domain_neib; d++){
        fprintf(fp,"%d\n",lpod_lb->meta_domain_list_neib[d]);
    }

    for(int i = 0; i < lpod_lb->meta_domain_neib + 1; i++){
        fprintf(fp,"%d\n", num_search[i]);
    }

    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(lpod_lb->recv_id[j][i] != -1){
                int index = lpod_lb->recv_id[j][i];
                int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
                fprintf(fp,"%d\n", lpod_lb->local_id_for_sort[value]);
            }
        }
    }

    fclose(fp);
}

void LB_write_node_send(
    int* recv_index,
    int* num_send,
    int** send_id,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.send.%d", INPUT_FILENAME_NODE, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"%d %d\n", lpod_lb->meta_domain_neib, num_send[lpod_lb->meta_domain_neib]);

    for(int d = 0; d < lpod_lb->meta_domain_neib; d++){
        fprintf(fp,"%d\n",lpod_lb->meta_domain_list_neib[d]);
    }

    for(int i = 0; i < lpod_lb->meta_domain_neib + 1; i++){
        fprintf(fp,"%d\n", num_send[i]);
    }

    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(send_id[j][i] != -1){
                int index = send_id[j][i];
                int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
                fprintf(fp,"%d\n", lpod_lb->local_id_for_sort[value]);
            }
        }
    }

    fclose(fp);
}

void LB_write_D_bc(
    int* global_D_bc_index,
    const int return_val,
    double* global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    char fname_out[BUFFER_SIZE];
    FILE* fp;
    int block_size = 1;

    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.%d", INPUT_FILENAME_D_BC, lpod_lb->myrank);    
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    int count = 0;

    fprintf(fp,"%d %d\n", lpod_lb->total_num_D_bc - return_val, block_size);

    for(int i = 0; i < lpod_lb->total_num_D_bc; i++){
        if(global_D_bc_index[i] != -1){
            int index = global_D_bc_index[i];
            int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
            int local_id = lpod_lb->local_id_for_sort[value];
            int block_id = 0;
        
            fprintf(fp,"%d %d %e\n", local_id, block_id, global_imposed_D_val[local_id]);
        }
    }

    fclose(fp);
}

void LB_write_D_bc_p(
    int* global_D_bc_index,
    const int return_val,
    double* global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const char*     label,
    const char* directory)
{
    char fname_out[BUFFER_SIZE];
    FILE* fp;
    int block_size = 1;

    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.dat.%d", label, lpod_lb->myrank);    
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    int count = 0;

    fprintf(fp,"%d %d\n", lpod_lb->total_num_D_bc - return_val, block_size);

    for(int i = 0; i < lpod_lb->total_num_D_bc; i++){
        if(global_D_bc_index[i] != -1){
            int index = global_D_bc_index[i];
            int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
            int local_id = lpod_lb->local_id_for_sort[value];
            int block_id = 0;
        
            fprintf(fp,"%d %d %e\n", local_id, block_id, global_imposed_D_val[i]);
        }
    }

    fclose(fp);
}

void LB_write_D_bc_v(
    int*            global_D_bc_index,
    const int       return_val,
    double**         global_imposed_D_val,
    LPOD_LB*        lpod_lb,
    const int       dof,
    const char*     directory)
{
    char fname_out[BUFFER_SIZE];
    FILE* fp;
    int block_size = dof;
    const int myrank = monolis_mpi_get_global_my_rank();

    snprintf(fname_out, BUFFER_SIZE, "merged_graph//D_bc_v.dat.%d", lpod_lb->myrank);    
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    int count = 0;

    fprintf(fp,"%d %d\n", (lpod_lb->total_num_D_bc - return_val)*dof, block_size);

    for(int i = 0; i < lpod_lb->total_num_D_bc; i++){
        if(global_D_bc_index[i] != -1){
            int index = global_D_bc_index[i];
            int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
            int local_id = lpod_lb->local_id_for_sort[value];
            
            for(int j = 0; j < dof; j++){
                fprintf(fp,"%d %d %e\n", local_id, j, global_imposed_D_val[i][j]);
            }
        }
    }
     
    fclose(fp);
}


void LB_write_elem(
    const int local_node_dof,
    const int num_internal,
    const int num_total,
    const int total_n_overlap,
    int* elem_id_bef_sort_overlap,
    int* id_overlap,
    int* sum_internal_global_id,
    int** global_conn,
    int* sum_internal_conn_id,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.%d", INPUT_FILENAME_ELEM, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"%d %d\n", num_total, local_node_dof);

    int value = 0 ;
    //内部要素とoverlap要素のidはそれぞれ分けているため、そのまま引っ張ってくることが可能
    for(int i = 0; i < num_internal; i++){
        if(sum_internal_global_id[i] != -1){
            int index = sum_internal_conn_id[i];
            for(int j = 0; j < 8; j++){
                value = ROM_BB_binarySearch(lpod_lb->local_id, global_conn[index][j], lpod_lb->total_num_nodes_lb);
                fprintf(fp,"%d ", lpod_lb->local_id_for_sort[value]);
            }
            fprintf(fp,"\n");
        }
    }

    for(int i = 0; i < total_n_overlap; i++){
        if(elem_id_bef_sort_overlap[i] != -1){
            int index = id_overlap[i];
            for(int j = 0; j < 8; j++){
                value = ROM_BB_binarySearch(lpod_lb->local_id, global_conn[index][j], lpod_lb->total_num_nodes_lb);
                fprintf(fp,"%d ", lpod_lb->local_id_for_sort[value]);
            }
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
}


void LB_write_elem_internal(
    const int ndof, 
    const int num,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.n_internal.%d", INPUT_FILENAME_ELEM, lpod_lb->myrank);
        
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    fprintf(fp,"#n_internal %d\n", ndof);
    fprintf(fp, "%d", num);
    fclose(fp);
}


void LB_write_elem_2(
    const int local_node_dof,
    const int num_total,
    const int num_internal,
    int* sum_internal_global_id,
    int** global_conn,
    int* sum_internal_conn_id,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.%d", INPUT_FILENAME_ELEM, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"%d %d\n", num_total, local_node_dof);

    int value = 0 ;
    for(int i = 0; i < num_internal; i++){
        if(sum_internal_global_id[i] != -1){
            int index = sum_internal_conn_id[i];

            for(int j = 0; j < 8; j++){
                value = ROM_BB_binarySearch(lpod_lb->local_id, global_conn[index][j], lpod_lb->total_num_nodes_lb);
                fprintf(fp,"%d ", lpod_lb->local_id_for_sort[value]);
            }
            fprintf(fp,"\n");
        }
    }
    fclose(fp);
}

void LB_write_elem_id(
    const int ndof, 
    const int total_num_elem,
    const int total_num_elem_internal,
    const int total_num_elem_overlap,
    int* sum_internal_global_id,
    int* elem_id_bef_sort_overlap,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.id.%d", INPUT_FILENAME_ELEM, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"#id\n");
    fprintf(fp,"%d %d\n", total_num_elem, ndof);

    for(int i = 0; i < total_num_elem_internal; i++){
        if(sum_internal_global_id[i] != -1){
            fprintf(fp,"%d\n",sum_internal_global_id[i]);
        }
    }

    for(int i = 0; i < total_num_elem_overlap; i++){
        if(elem_id_bef_sort_overlap[i] != -1){
            fprintf(fp,"%d\n",elem_id_bef_sort_overlap[i]);
        }
    }
    fclose(fp);
}

void LB_write_elem_id_2(
    const int ndof, 
    const int total_num_elem,
    const int total_num_elem_internal,
    int* sum_internal_global_id,
    LPOD_LB*        lpod_lb,
    const char* directory)
{
    FILE* fp;
    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph//%s.id.%d", INPUT_FILENAME_ELEM, lpod_lb->myrank);
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);

    fprintf(fp,"#id\n");
    fprintf(fp,"%d %d\n", total_num_elem, ndof);

    for(int i = 0; i < total_num_elem_internal; i++){
        if(sum_internal_global_id[i] != -1){
            fprintf(fp,"%d\n",sum_internal_global_id[i]);
        }
    }

    fclose(fp);
}