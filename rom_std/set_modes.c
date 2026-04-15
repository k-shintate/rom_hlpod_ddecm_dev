#include "set_modes.h"

static const int BUFFER_SIZE = 10000;

/*free*/
void ROM_std_hlpod_free_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		dof)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes_max);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_num_nodes*dof);
}

void ROM_std_hlpod_free_local_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_2nd_subdomains,
    const int 		num_modes,
    const int 		dof)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes*num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
}

void ROM_std_hlpod_free_local_podmodes_para(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		n_internal_vertex,
    const int 		num_2nd_subdomains,
    const int 		num_modes,
    const int 		dof)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes*num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
}

void ROM_std_hlpod_free_global_podmodes(
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes,
    const int 		dof)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes);
}

/*set*/
void ROM_std_hlpod_set_podmodes(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		num_snapshots,
    const double    rom_epsilon,
    const int 	  	dof,
    const char*  	label,
    const char*  	directory)
{
    double** S; double* V; double** D;

    const int comm = monolis_mpi_get_self_comm();
    const int total_dof = total_num_nodes * dof;
    int scalapack_comm;

    S = BB_std_calloc_2d_double(S, total_dof, num_snapshots);
    V = BB_std_calloc_1d_double(V, num_snapshots);
    D = BB_std_calloc_2d_double(D, num_snapshots, num_snapshots);
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_dof, num_modes_max);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_dof);

    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    double t1 = monolis_get_time();
    monolis_scalapack_gesvd_R(
            total_dof,
            num_snapshots, 
            hlpod_mat->snapmat, 
            S, 
            V, 
            D, 
            comm,
            scalapack_comm);
    double t2 = monolis_get_time();

    if(rom_epsilon == 1){
        hlpod_vals->num_modes = num_modes_max;
    }
    else{
        hlpod_vals->num_modes = 
            ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots,
                rom_epsilon);
    }

    for(int j = 0; j < hlpod_vals->num_modes; j++){
        for(int i = 0; i < total_dof; i++){
            hlpod_mat->pod_modes[i][j] = S[i][j];
        }
    }

    ROM_std_hlpod_write_pod_modes(S, hlpod_vals->num_modes, total_num_nodes, 0, dof, label, directory);
    ROM_std_hlpod_write_singular_values(V, num_snapshots, 0, label, directory);
    ROM_std_hlpod_write_time_svd(t2-t1, 0, label, directory);
    ROM_std_hlpod_write_num_modes(hlpod_vals->num_modes, 0, label, directory);

    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, 1);
    hlpod_mat->num_modes_internal[0] = hlpod_vals->num_modes;
    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd, 1);
    hlpod_mat->n_internal_vertex_subd[0] = total_num_nodes;

    monolis_scalapack_comm_finalize(scalapack_comm);

    BB_std_free_2d_double(S, total_dof, num_snapshots);
    BB_std_free_1d_double(V, num_snapshots);
    BB_std_free_2d_double(D, num_snapshots, num_snapshots);

}

void ROM_std_hlpod_read_podmodes(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*	    hlpod_mat,
    const int 		total_num_nodes,
    const int 		num_modes_max,
    const int 		num_snapshots,
    const double    rom_epsilon,
    const int 	  	dof,
    const char*  	label,
    const char*  	directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    const int total_dof = total_num_nodes * dof;

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_dof, num_modes_max);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_dof);

    for(int i = 0; i < total_num_nodes; i++) {
        hlpod_mat->node_id[i] = i;
    }

    int n_basis = ROM_std_hlpod_read_num_modes(
            0,
            label,
            directory);

    hlpod_vals->num_modes = n_basis;
        
    ROM_std_hlpod_read_pod_modes_node(fp, hlpod_mat->pod_modes, dof, total_num_nodes, n_basis, 0, label, directory);

    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, 1);
    hlpod_mat->num_modes_internal[0] = n_basis;
    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd, 1);
    hlpod_mat->n_internal_vertex_subd[0] = total_num_nodes;

}

void ROM_std_hlpod_set_podmodes_local(
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int       num_modes,
    const int       num_snapshots,
    const int       num_2nd_subdomains,
    const double    rom_epsilon,
    const int		dof,
    const char*     label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    double** snapmat_local;
    int n_internal_vertex;
    double** S; double* V; double** D;

    const int comm = monolis_mpi_get_self_comm();
    int scalapack_comm;
    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    int total = 0;
    int index = 0;
    int n_basis = 0;

    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);

    for(int m = 0; m < num_2nd_subdomains; m++){
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];
        snapmat_local = BB_std_calloc_2d_double(snapmat_local, n_internal_vertex*dof,num_snapshots);
        
         for(int i = 0; i < n_internal_vertex; i++){
            for(int j = 0; j < num_snapshots; j++){
                for(int k = 0; k < dof; k++){
                    snapmat_local[i*dof + k][j] = hlpod_mat->snapmat[hlpod_mat->node_id[index]*dof + k][j];
                }
            }
            index++;
        }

        S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, num_snapshots);
        V = BB_std_calloc_1d_double(V, num_snapshots);
        D = BB_std_calloc_2d_double(D, num_snapshots, num_snapshots);

        double t1 = monolis_get_time();
        monolis_scalapack_gesvd_R(
            n_internal_vertex*dof,
            num_snapshots, 
            snapmat_local, 
            S, 
            V, 
            D, 
            comm,
            scalapack_comm);
        double t2 = monolis_get_time();

        if(rom_epsilon == 1){
            n_basis = num_modes;
            
            hlpod_mat->num_modes_internal[m] = n_basis;
        }
        else{
            n_basis = ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots,
                rom_epsilon);
		if(n_basis > num_modes){
			hlpod_mat->num_modes_internal[m] = num_modes;
		}
		else{
            		hlpod_mat->num_modes_internal[m] = n_basis;
		}
        }

        ROM_std_hlpod_write_singular_values(V, num_snapshots, m, label, directory);
        ROM_std_hlpod_write_time_svd(t2-t1, m, label, directory);
        ROM_std_hlpod_write_num_modes(n_basis, m, label, directory);
        ROM_std_hlpod_write_pod_modes(S, n_basis, n_internal_vertex, m, dof, label, directory);

        BB_std_free_2d_double(S, n_internal_vertex*dof, num_snapshots);
        BB_std_free_1d_double(V, num_snapshots);
        BB_std_free_2d_double(D, num_snapshots, num_snapshots);

        BB_std_free_2d_double(snapmat_local, n_internal_vertex*dof, num_snapshots);

    }

    hlpod_vals->num_modes = total;

    monolis_scalapack_comm_finalize(scalapack_comm);
    BB_std_free_2d_double(hlpod_mat->snapmat, total_num_nodes*dof, num_snapshots);
}

void ROM_std_hlpod_read_podmodes_local(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int 		total_num_nodes,
    const int       num_modes,
    const int       num_snapshots,
    const int       num_2nd_subdomains,
    const int 		dof,
    const char*  	label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int n_internal_vertex;
    double** S;

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes * num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
    hlpod_mat->subdomain_id_in_nodes = BB_std_calloc_1d_int(hlpod_mat->subdomain_id_in_nodes, total_num_nodes*dof);

    int total = 0;
    int Total = 0;
    int index = 0;
    for(int m = 0; m < num_2nd_subdomains; m++){
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];

        int n_basis = ROM_std_hlpod_read_num_modes(
            m,
            label,
            directory);
        
        hlpod_mat->num_modes_internal[m] = n_basis;

        S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, n_basis);

        ROM_std_hlpod_read_pod_modes_node(fp, S, dof, n_internal_vertex, n_basis, m, label, directory);

        for(int i = 0; i < n_internal_vertex; i++){
            for(int j = 0; j < n_basis; j++){
                for(int k = 0; k < dof; k++){
                    hlpod_mat->pod_modes[hlpod_mat->node_id[index] * dof + k][total + j] = S[i*dof + k][j] * sqrt(n_internal_vertex * dof);
                }
            }
            hlpod_mat->subdomain_id_in_nodes[hlpod_mat->node_id[index]] = m;
            index++;
        }

        total += n_basis;
        Total += num_modes;

        BB_std_free_2d_double(S, n_internal_vertex*dof, n_basis);
    }

    hlpod_vals->num_modes = total;

}


void ROM_std_hlpod_set_podmodes_local_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int		num_2nd_subdomains,
    const double    rom_epsilon,
    const int 		dof,
    const char*     label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    double** S; double* V; double** D;

    double** snapmat_local;
    int* node_id_local;
    int n_internal_vertex;

    int comm = monolis_mpi_get_self_comm();
    int scalapack_comm;
    const int myrank = monolis_mpi_get_global_my_rank();

    int n_basis;

    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, (total_num_nodes)* dof, num_modes * num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);

    int index = 0;
    int index_row = 0;
    int index_column = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];

        snapmat_local = BB_std_calloc_2d_double(snapmat_local, n_internal_vertex * dof,num_snapshots);

        for(int j = 0; j < num_snapshots; j++){
            for(int i = 0; i < n_internal_vertex* dof; i++){
                snapmat_local[i][j] = hlpod_mat->snapmat[index_row + i][j];
            }
        }

        S = BB_std_calloc_2d_double(S, n_internal_vertex* dof, num_snapshots);
        V = BB_std_calloc_1d_double(V, num_snapshots);
        D = BB_std_calloc_2d_double(D, num_snapshots, num_snapshots);

        double t1 = monolis_get_time();
        monolis_scalapack_gesvd_R(
            n_internal_vertex* dof,
            num_snapshots, 
            snapmat_local, 
            S, 
            V, 
            D, 
            comm,
            scalapack_comm);
        double t2 = monolis_get_time();

        if(rom_epsilon == 1){
            n_basis = num_modes;
            hlpod_mat->num_modes_internal[m] = n_basis;
        }
        else{
            n_basis = ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots,
                rom_epsilon);
            
            hlpod_mat->num_modes_internal[m] = n_basis;
        }
        
        index_row += n_internal_vertex* dof;

        ROM_std_hlpod_write_singular_values(V, num_snapshots, hlpod_meta->subdomain_id[m], label, directory);		
        ROM_std_hlpod_write_time_svd(t2-t1, hlpod_meta->subdomain_id[m], label, directory);
        ROM_std_hlpod_write_num_modes(n_basis, hlpod_meta->subdomain_id[m], label, directory);
        ROM_std_hlpod_write_pod_modes(S, n_basis, n_internal_vertex, hlpod_meta->subdomain_id[m], dof, label, directory);

        BB_std_free_2d_double(S, n_internal_vertex* dof, num_snapshots);
        BB_std_free_1d_double(V, num_snapshots);
        BB_std_free_2d_double(D, num_snapshots, num_snapshots);

        BB_std_free_2d_double(snapmat_local, n_internal_vertex* dof, num_snapshots);
    }

    hlpod_vals->num_modes = index_column;

    monolis_scalapack_comm_finalize(scalapack_comm);
    BB_std_free_2d_double(hlpod_mat->snapmat, N_internal_vertex * dof, num_snapshots);
}

/*
void ROM_std_hlpod_read_podmodes_local_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int 		N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int		num_2nd_subdomains,
    const int 		dof,
    const double    rom_epsilon,
    const char*		label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** S;
    int n_internal_vertex;
    const int myrank = monolis_mpi_get_global_my_rank();

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes * num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);

        //for arbit dof ddecm
        hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);
        


    int index = 0;
    int index_row = 0;
    int index_column = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];

        int num_snapshots1 = 100;

        double* V = BB_std_calloc_1d_double(V, num_snapshots1);

        ROM_std_hlpod_read_singular_values(V, num_snapshots1, hlpod_meta->subdomain_id[m], label, directory);
        int n_basis = ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots1,
                rom_epsilon);
        //n_basis = 20;
//printf("\n\nn_basis=%d rom_epsilon = %lf\n", n_basis, rom_epsilon);
        if(n_basis > num_modes){
                n_basis = num_modes;
        }

printf("\n\nn_basis=%d rom_epsilon = %lf\n", n_basis, rom_epsilon);

//double t = monolis_get_time_global_sync();

        //exit(1);

        hlpod_mat->num_modes_internal[m] = n_basis;
        //hlpod_mat->num_modes_internal[m] = 10;

        S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, n_basis);
        ROM_std_hlpod_read_pod_modes_node(fp, S, dof, n_internal_vertex, n_basis, hlpod_meta->subdomain_id[m], label, directory);

        for(int j = 0; j < n_basis; j++){
            for(int i = 0; i < n_internal_vertex; i++){
                for(int k = 0; k < dof; k++){
                    hlpod_mat->pod_modes[index_row*dof + i*dof + k][index_column + j] = S[i*dof + k][j] * sqrt(n_internal_vertex * dof);
                    //hlpod_mat->pod_modes[index_row*dof + i*dof + k][index_column + j] = S[i*dof + k][j];
                }
            }
        }

        //for arbit dof ddecm
            for(int i = 0; i < n_internal_vertex; i++){
                        hlpod_mat->subdomain_id_in_nodes_internal[index_row + i][m] = m + 1;
        }
        

        index_row += n_internal_vertex;
        index_column += n_basis;
	//index_column += num_modes;

        BB_std_free_2d_double(S, n_internal_vertex*dof, n_basis);
    }

    double t = monolis_get_time_global_sync();

    hlpod_vals->num_modes = index_column;

}
*/


void ROM_std_hlpod_read_podmodes_local_para(
    HLPOD_VALUES*       hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*         hlpod_meta,
    const int           total_num_nodes,
    const int           N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int           num_2nd_subdomains,
    const int           dof,
    const double    rom_epsilon,
    const char*         label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** S;
    int n_internal_vertex;
    const int myrank = monolis_mpi_get_global_my_rank();

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes * num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);

        //for arbit dof ddecm
        hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);


    int index = 0;
    int index_row = 0;
    int index_column = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];

        int num_snapshots1 = 100;

        double* V = BB_std_calloc_1d_double(V, num_snapshots1);

        ROM_std_hlpod_read_singular_values(V, num_snapshots1, hlpod_meta->subdomain_id[m], label, directory);
        int n_basis = ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots1,
                rom_epsilon);
        //n_basis = 20;
//printf("\n\nn_basis=%d rom_epsilon = %lf\n", n_basis, rom_epsilon);
        if(n_basis > num_modes){
                n_basis = num_modes;
        }

    printf("\n\nn_basis=%d rom_epsilon = %lf\n", n_basis, rom_epsilon);
    //double t = monolis_get_time_global_sync();
    //exit(1);

        hlpod_mat->num_modes_internal[m] = n_basis;
	//hlpod_mat->num_modes_internal[m] = num_modes;
        //hlpod_mat->num_modes_internal[m] = 10;

        S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, n_basis);
        ROM_std_hlpod_read_pod_modes_node(fp, S, dof, n_internal_vertex, n_basis, hlpod_meta->subdomain_id[m], label, directory);

        for(int j = 0; j < n_basis; j++){
            for(int i = 0; i < n_internal_vertex; i++){
                for(int k = 0; k < dof; k++){
                    hlpod_mat->pod_modes[index_row*dof + i*dof + k][index_column + j] = S[i*dof + k][j] * sqrt(n_internal_vertex * dof);
                    //hlpod_mat->pod_modes[index_row*dof + i*dof + k][index_column + j] = S[i*dof + k][j];
                }
            }
        }

        //for arbit dof ddecm
            for(int i = 0; i < n_internal_vertex; i++){
                        hlpod_mat->subdomain_id_in_nodes_internal[index_row + i][m] = m + 1;
        }


        index_row += n_internal_vertex;
        index_column += n_basis;
        //index_column += num_modes;

        BB_std_free_2d_double(S, n_internal_vertex*dof, n_basis);
    }

    double t = monolis_get_time_global_sync();

    hlpod_vals->num_modes = index_column;

}


/*
void ROM_std_hlpod_read_podmodes_local_para(
    HLPOD_VALUES*       hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*         hlpod_meta,
    const int           total_num_nodes,
    const int           N_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const int           num_2nd_subdomains,
    const int           dof,
    const double    rom_epsilon,
    const char*         label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** S;
    int n_internal_vertex;
    const int myrank = monolis_mpi_get_global_my_rank();

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes * num_2nd_subdomains);
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);

	//for arbit dof ddecm
	hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);
	

    printf("\nrank %d, total_num_nodes = %d, n_internal_vertex = %d, num_2nd_subdomains = %d\n", myrank, total_num_nodes, N_internal_vertex, num_2nd_subdomains);


    int index = 0;
    int index_row = 0;
    int index_column = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex = hlpod_mat->n_internal_vertex_subd[m];

        int n_basis = ROM_std_hlpod_read_num_modes(
            hlpod_meta->subdomain_id[m],
            label,
            directory);

    	n_basis = 10;

        hlpod_mat->num_modes_internal[m] = n_basis;

        S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, n_basis);
        ROM_std_hlpod_read_pod_modes_node(fp, S, dof, n_internal_vertex, n_basis, hlpod_meta->subdomain_id[m], label, directory);

        for(int j = 0; j < n_basis; j++){
            for(int i = 0; i < n_internal_vertex; i++){
                for(int k = 0; k < dof; k++){
                    hlpod_mat->pod_modes[index_row*dof + i*dof + k][index_column + j] = S[i*dof + k][j] * sqrt(n_internal_vertex * dof);
                }
            }
        }

    	//for arbit dof ddecm
	    for(int i = 0; i < n_internal_vertex; i++){
			hlpod_mat->subdomain_id_in_nodes_internal[index_row + i][m] = m + 1;
        }
	

        index_row += n_internal_vertex;
        index_column += n_basis;

        BB_std_free_2d_double(S, n_internal_vertex*dof, n_basis);
    }

    double t = monolis_get_time_global_sync();

    hlpod_vals->num_modes = index_column;

}
*/


void ROM_std_hlpod_set_podmodes_global_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*   hlpod_mat,
    double**		snapmat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    const int       num_modes,
    const int       num_snapshots,
    const double    rom_epsilon,
    const int 		dof,
    const char*		label,
    const char*     directory)
{
    int n_basis;
    double** S; double* V; double** D;

    int comm = monolis_mpi_get_global_comm();
    int scalapack_comm;
    const int myrank = monolis_mpi_get_global_my_rank();
    monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    S = BB_std_calloc_2d_double(S, n_internal_vertex* dof, num_snapshots);
    V = BB_std_calloc_1d_double(V, num_snapshots);
    D = BB_std_calloc_2d_double(D, num_snapshots, num_snapshots);

    double t1 = monolis_get_time();
    monolis_scalapack_gesvd_R(
        n_internal_vertex* dof,
        num_snapshots, 
        snapmat, 
        S, 
        V, 
        D, 
        comm,
        scalapack_comm);
    double t2 = monolis_get_time();

    if(rom_epsilon == 1){
        n_basis =  num_modes;
        hlpod_vals->num_modes = num_modes;
    }
    else{
        n_basis =
            ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots,
                rom_epsilon);

        hlpod_vals->num_modes =
            ROM_BB_estimate_num_pod_modes(
                V,
                num_snapshots,
                rom_epsilon);
    }

    ROM_std_hlpod_write_singular_values(V, num_snapshots, myrank, label, directory);
    ROM_std_hlpod_write_num_modes(n_basis, myrank, label, directory);
    //ROM_std_hlpod_write_pod_modes(S, n_basis, n_internal_vertex, myrank, dof, label, directory);
    ROM_std_hlpod_write_pod_modes(snapmat, n_basis, n_internal_vertex, myrank, dof, label, directory);    
    ROM_std_hlpod_write_time_svd(t2-t1, myrank, label, directory);

    BB_std_free_2d_double(S, n_internal_vertex* dof, num_snapshots);
    BB_std_free_1d_double(V, num_snapshots);
    BB_std_free_2d_double(D, num_snapshots, num_snapshots);

}



void ROM_std_hlpod_read_podmodes_global_para(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    const int       num_modes,
    const int		dof,
    const char*		label,
    const char*     directory)
{
    FILE* fp;
    const int myrank = monolis_mpi_get_global_my_rank();

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes*dof, num_modes);

    int n_basis = ROM_std_hlpod_read_num_modes(
        myrank,
        label,
        directory);

    double** S;
    S = BB_std_calloc_2d_double(S, n_internal_vertex*dof, n_basis);

    ROM_std_hlpod_read_pod_modes_node(fp, S, dof, n_internal_vertex, n_basis, myrank, label, directory);

    for(int j = 0; j < n_basis; j++){
        for(int i = 0; i < n_internal_vertex; i++){
            for(int k = 0; k < dof; k++){
                hlpod_mat->pod_modes[i*dof + k][j] = S[i*dof + k][j];
            }
        }
    }

    BB_std_free_2d_double(S, n_internal_vertex*dof, n_basis);

    hlpod_vals->num_modes = n_basis;
}


/*set diag*/
void ROM_std_hlpod_set_podmodes_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    double** 		v,
    double** 		p,
    const int 		total_num_nodes,
    const int 		num_modes_max_1,
    const int 		num_modes_max_2,
    const int		dof_1,
    const int 		dof_2,
    const char*  	directory)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes * (dof_1+dof_2), num_modes_max_1 + num_modes_max_2);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_num_nodes);

    for(int j = 0; j < num_modes_max_1; j++){
        for(int i = 0; i < total_num_nodes; i++){
            for(int k = 0; k < 3; k++){
                hlpod_mat->pod_modes[i*4 + k][j] = v[i*3 + k][j];
            }
        }
    }

    for(int j = 0; j < num_modes_max_2; j++){
        for(int i = 0; i < total_num_nodes; i++){
            hlpod_mat->pod_modes[i * 4 + 3][j + num_modes_max_1] = p[i][j];
        }
    }
    
    for(int i = 0; i < total_num_nodes; i++) {
        hlpod_mat->node_id[i] = i;
    }

    hlpod_vals->num_modes = num_modes_max_1 + num_modes_max_2;

    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, 1);
    hlpod_mat->num_modes_internal[0] = hlpod_vals->num_modes;
    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd, 1);
    hlpod_mat->n_internal_vertex_subd[0] = total_num_nodes;
}

void ROM_std_hlpod_set_podmodes_local_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    int* 			n_internal_vertex_subd,
    int* 			node_id_local,
    const int 		num_2nd_subdomains,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, (total_num_nodes)* 4, (num_modes_max_1 + num_modes_max_2));
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd , num_2nd_subdomains);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_num_nodes);

    int index = 0;
    int index_row = 0;
    int index_column = 0;
    int index_column_v = 0;
    int index_column_p = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        for(int j = 0; j < num_modes_v[m]; j++){
            for(int i = 0; i < n_internal_vertex_subd[m]; i++){
                for(int k = 0; k < 3; k++){
                    hlpod_mat->pod_modes[node_id_local[i + index_row]*4 + k][index_column + j] = 
                    v[node_id_local[i + index_row]*3 + k][j + index_column_v];
                }
            }
        }

        index_column_v += num_modes_v[m];
        index_column += num_modes_v[m];

        for(int j = 0; j < num_modes_p[m]; j++){
            for(int i = 0; i < n_internal_vertex_subd[m]; i++){
                hlpod_mat->pod_modes[node_id_local[i + index_row]*4  + 3][index_column + j] = 
                p[node_id_local[i + index_row]][j + index_column_p];
            }
        }

        index_column_p += num_modes_p[m];
        index_column += num_modes_p[m];
        index_row += n_internal_vertex_subd[m];

        for(int i = 0; i < n_internal_vertex_subd[m]; i++){
            hlpod_mat->node_id[index] = node_id_local[index];
            index++;
        }

        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex_subd[m];
        hlpod_mat->num_modes_internal[m] = num_modes_v[m] + num_modes_p[m];
    }

    hlpod_vals->num_modes = index_column;

}

void ROM_std_hlpod_set_podmodes_local_para_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** snapmat_local;
    int* node_id_local;
    int n_internal_vertex_1stdd;

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, (total_num_nodes)* 4, (num_modes_max_1 + num_modes_max_2));
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
	//for arbit dof ddecm
	hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);
	/**/
    
    int index = 0;
    int index_row = 0;
    int index_column = 0;
    int index_column_v = 0;
    int index_column_p = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex_1stdd = hlpod_mat->n_internal_vertex_subd[m];

        for(int j = 0; j < num_modes_v[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                for(int k = 0; k < 3; k++){
                    hlpod_mat->pod_modes[index_row*4 + i*4 + k][index_column + j] = v[index_row*3 + i*3 + k][index_column_v + j];
                }
            }
        }
    
        index_column_v += num_modes_v[m];
        index_column += num_modes_v[m];

        hlpod_mat->num_modes_internal[m] = num_modes_v[m] + num_modes_p[m];

        for(int j = 0; j < num_modes_p[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                hlpod_mat->pod_modes[index_row*4 + i*4 + 3][index_column + j] = p[index_row + i][index_column_p + j];
            }
        }

        for(int i = 0; i < n_internal_vertex_1stdd; i++){
            hlpod_mat->node_id[index] = index;
            hlpod_mat->subdomain_id_in_nodes_internal[index][m] = m + 1;
            index++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex_1stdd;

        index_row += n_internal_vertex_1stdd;
        index_column_p += num_modes_p[m];
        index_column += num_modes_p[m];
    }

    hlpod_vals->num_modes = index_column;

}

void ROM_std_hlpod_set_podmodes_local_para_diag_decoupled(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** snapmat_local;
    int* node_id_local;
    int n_internal_vertex_1stdd;

    hlpod_mat->pod_modes_decoupled_p = BB_std_calloc_2d_double(hlpod_mat->pod_modes_decoupled_p, (total_num_nodes)* 4, (num_modes_max_1 + num_modes_max_2));
    hlpod_mat->pod_modes_decoupled_v = BB_std_calloc_2d_double(hlpod_mat->pod_modes_decoupled_v, (total_num_nodes)* 4, (num_modes_max_1 + num_modes_max_2));
    //hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
	//for arbit dof ddecm
	//hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);
	/**/
    
    int index = 0;
    int index_row = 0;
    int index_column = 0;
    int index_column_v = 0;
    int index_column_p = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex_1stdd = hlpod_mat->n_internal_vertex_subd[m];

        /*
        for(int j = 0; j < num_modes_v[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                for(int k = 0; k < 3; k++){
                    hlpod_mat->pod_modes[index_row*4 + i*4 + k][index_column + j] = v[index_row*3 + i*3 + k][index_column_v + j];
                }
            }
        }
        */
    
        index_column_v += num_modes_v[m];
        index_column += num_modes_v[m];

        //hlpod_mat->num_modes_internal[m] = num_modes_v[m] + num_modes_p[m];

        for(int j = 0; j < num_modes_p[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                hlpod_mat->pod_modes_decoupled_p[index_row*4 + i*4 + 3][index_column + j] = p[index_row + i][index_column_p + j];
            }
        }
/*
        for(int i = 0; i < n_internal_vertex_1stdd; i++){
            hlpod_mat->node_id[index] = index;
            hlpod_mat->subdomain_id_in_nodes_internal[index][m] = m + 1;
            index++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex_1stdd;
*/
        index_row += n_internal_vertex_1stdd;
        index_column_p += num_modes_p[m];
        index_column += num_modes_p[m];
    }

    index = 0;
    index_row = 0;
    index_column = 0;
    index_column_v = 0;
    index_column_p = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex_1stdd = hlpod_mat->n_internal_vertex_subd[m];

        
        for(int j = 0; j < num_modes_v[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                for(int k = 0; k < 3; k++){
                    hlpod_mat->pod_modes_decoupled_v[index_row*4 + i*4 + k][index_column + j] = v[index_row*3 + i*3 + k][index_column_v + j];
                }
            }
        }
        
    
        index_column_v += num_modes_v[m];
        index_column += num_modes_v[m];

        //hlpod_mat->num_modes_internal[m] = num_modes_v[m] + num_modes_p[m];
        /*
        for(int j = 0; j < num_modes_p[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                hlpod_mat->pod_modes_decoupled_p[index_row*4 + i*4 + 3][index_column + j] = p[index_row + i][index_column_p + j];
            }
        }
        */
/*
        for(int i = 0; i < n_internal_vertex_1stdd; i++){
            hlpod_mat->node_id[index] = index;
            hlpod_mat->subdomain_id_in_nodes_internal[index][m] = m + 1;
            index++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex_1stdd;
*/
        index_row += n_internal_vertex_1stdd;
        index_column_p += num_modes_p[m];
        index_column += num_modes_p[m];
    }

    //hlpod_vals->num_modes = index_column;

}

void ROM_std_hlpod_set_podmodes_global_para_diag(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, (total_num_nodes)*4, num_modes_max_1+num_modes_max_2);

    for(int j = 0; j < num_modes_max_1; j++){
        for(int i = 0; i < total_num_nodes; i++){
            for(int k = 0; k < 3; k++){
                hlpod_mat->pod_modes[i*4 + k][j] = v[i*3 + k][j];
            }
        }
    }

    for(int j = 0; j < num_modes_max_2; j++){
        for(int i = 0; i < total_num_nodes; i++){
            hlpod_mat->pod_modes[i*4 + 3][j + num_modes_max_1] = p[i][j];
        }
    }

    hlpod_vals->num_modes = num_modes_max_1+num_modes_max_2;
}

void ROM_std_hlpod_set_podmodes_global_para_Aphi(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const char*     directory)
{
    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, total_num_nodes, num_modes_max_1+num_modes_max_2);

    for(int j = 0; j < num_modes_max_1; j++){
        for(int i = 0; i < total_num_nodes; i++){
            hlpod_mat->pod_modes[i][j] = v[i][j];
        }
    }

    for(int j = 0; j < num_modes_max_2; j++){
        for(int i = 0; i < total_num_nodes; i++){
            hlpod_mat->pod_modes[i][j + num_modes_max_1] = p[i][j];
        }
    }

    hlpod_vals->num_modes = num_modes_max_1+num_modes_max_2;
}

/*for hyper-reduction*/
void HROM_ddecm_set_podbasis_ovl(
    MONOLIS_COM*  	monolis_com,
	HLPOD_VALUES*		hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
	const int		total_num_nodes,
    const int       dof)
{
	hlpod_mat->pod_basis_hr = BB_std_calloc_2d_double(hlpod_mat->pod_basis_hr, total_num_nodes*dof, hlpod_vals->num_modes_max);

    const int NDOF  = total_num_nodes*dof;

    double* monolis_in;
    monolis_in = BB_std_calloc_1d_double(monolis_in, NDOF);

    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
		for(int k = 0; k < NDOF; k++){
			monolis_in[k] = 0;
		}
		for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
			monolis_in[j] = hlpod_mat->pod_modes[j][l];
		}
		monolis_mpi_update_R(monolis_com, NDOF, dof, monolis_in);
		for(int k = 0; k < NDOF; k++){
            hlpod_mat->pod_basis_hr[k][l] = monolis_in[k];
        }
	}

	BB_std_free_1d_double(monolis_in, NDOF);

}

void HROM_ddecm_set_podbasis_ovl_decoupled(
    MONOLIS_COM*  	monolis_com,
	HLPOD_VALUES*		hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
	const int		total_num_nodes,
    const int       dof)
{
	hlpod_mat->pod_basis_hr_decoupled_p = BB_std_calloc_2d_double(hlpod_mat->pod_basis_hr_decoupled_p, total_num_nodes*dof, hlpod_vals->num_modes_max);

    const int NDOF  = total_num_nodes*dof;

    double* monolis_in;
    monolis_in = BB_std_calloc_1d_double(monolis_in, NDOF);

    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
		for(int k = 0; k < NDOF; k++){
			monolis_in[k] = 0;
		}
		for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
			monolis_in[j] = hlpod_mat->pod_modes_decoupled_p[j][l];
		}
		monolis_mpi_update_R(monolis_com, NDOF, dof, monolis_in);
		for(int k = 0; k < NDOF; k++){
            hlpod_mat->pod_basis_hr_decoupled_p[k][l] = monolis_in[k];
        }
	}

    hlpod_mat->pod_basis_hr_decoupled_v = BB_std_calloc_2d_double(hlpod_mat->pod_basis_hr_decoupled_v, total_num_nodes*dof, hlpod_vals->num_modes_max);


    for(int l = 0; l < hlpod_vals->num_modes_max; l++){
		for(int k = 0; k < NDOF; k++){
			monolis_in[k] = 0;
		}
		for(int j = 0; j < monolis_com->n_internal_vertex * dof; j++){
			monolis_in[j] = hlpod_mat->pod_modes_decoupled_v[j][l];
		}
		monolis_mpi_update_R(monolis_com, NDOF, dof, monolis_in);
		for(int k = 0; k < NDOF; k++){
            hlpod_mat->pod_basis_hr_decoupled_v[k][l] = monolis_in[k];
        }
	}

	BB_std_free_1d_double(monolis_in, NDOF);

}

void ROM_std_hlpod_set_podmodes_local_para_Aphi(
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*    	hlpod_meta,
    const int		total_num_nodes,
    const int 		n_internal_vertex,
    double** 		v,
    double** 		p,
    int* 			num_modes_v,
    int* 			num_modes_p,
    const int       num_modes_max_1,
    const int       num_modes_max_2,
    const int 		dof_1,
    const int 		dof_2,
    const int		num_2nd_subdomains,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    double** snapmat_local;
    int* node_id_local;
    int n_internal_vertex_1stdd;

    hlpod_mat->pod_modes = BB_std_calloc_2d_double(hlpod_mat->pod_modes, (total_num_nodes)* 4, (num_modes_max_1 + num_modes_max_2));
    hlpod_mat->num_modes_internal = BB_std_calloc_1d_int(hlpod_mat->num_modes_internal, num_2nd_subdomains);
	//for arbit dof ddecm
	hlpod_mat->subdomain_id_in_nodes_internal = BB_std_calloc_2d_int(hlpod_mat->subdomain_id_in_nodes_internal, total_num_nodes, num_2nd_subdomains);
	/**/
    
    int index = 0;
    int index_row = 0;
    int index_column = 0;
    int index_column_v = 0;
    int index_column_p = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        int subdomain_id = hlpod_meta->subdomain_id[m];
        n_internal_vertex_1stdd = hlpod_mat->n_internal_vertex_subd[m];

        for(int j = 0; j < num_modes_v[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                //for(int k = 0; k < 3; k++){
                hlpod_mat->pod_modes[index_row+ i][index_column + j] = v[index_row + i][index_column_v + j];
                //}
            }
        }
    
        index_column_v += num_modes_v[m];
        index_column += num_modes_v[m];

        hlpod_mat->num_modes_internal[m] = num_modes_v[m] + num_modes_p[m];

        for(int j = 0; j < num_modes_p[m]; j++){
            for(int i = 0; i < n_internal_vertex_1stdd; i++){
                hlpod_mat->pod_modes[index_row + i][index_column + j] = p[index_row + i][index_column_p + j];
            }
        }

        for(int i = 0; i < n_internal_vertex_1stdd; i++){
            hlpod_mat->node_id[index] = index;
            hlpod_mat->subdomain_id_in_nodes_internal[index][m] = m + 1;
            index++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex_1stdd;

        index_row += n_internal_vertex_1stdd;
        index_column_p += num_modes_p[m];
        index_column += num_modes_p[m];
    }

    hlpod_vals->num_modes = index_column;

}
