
#include "hlpod_comm.h"

void ROM_std_hlpod_get_neib_vec(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*	hlpod_vals,
    HLPOD_MAT* 	    hlpod_mat,
    const int 		num_modes,
    const int       ndof)
{
    int n_neib_vec;
    int n_vec = num_modes;    //自領域のベクトル数
    
    monolis_mpi_get_n_neib_vector(
        monolis_com,
        n_vec,
        &n_neib_vec);       //出力：自領域と隣接領域の合計ベクトル数

    hlpod_vals->n_neib_vec = n_neib_vec;

    const int np = monolis_com->n_internal_vertex + monolis_com->recv_index[monolis_com->recv_n_neib];
    const int n_internal_vertex = monolis_com->n_internal_vertex;

    double** my_vec;
    my_vec = BB_std_calloc_2d_double(my_vec, np*ndof, n_vec);

    for(int i = 0; i < n_vec; i++){	
        for(int j = 0; j < n_internal_vertex; j++){
            for(int k = 0; k < ndof; k++){
                my_vec[j*ndof + k][i] = hlpod_mat->pod_modes[j*ndof + k][i];
            }
        }
    }

    hlpod_mat->neib_vec = BB_std_calloc_2d_double(hlpod_mat->neib_vec, np*ndof, n_neib_vec);

    monolis_mpi_get_neib_vector_R(
        monolis_com,
        np,						//配列サイズ
        ndof,					//計算点が持つ自由度
        n_vec,					//自領域のベクトル数
        n_neib_vec,				//自領域と隣接領域の合計ベクトル数
        my_vec,					//自領域のベクトル
        hlpod_mat->neib_vec);	//自領域と隣接領域が並んだベクトル

}


void ROM_std_hlpod_get_neib_vec_save_memory(
    MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES*	hlpod_vals,
    const int 		num_modes)
{
    int n_neib_vec;
    int n_vec = num_modes;    //自領域のベクトル数
    
    monolis_mpi_get_n_neib_vector(
        monolis_com,
        n_vec,
        &n_neib_vec);       //出力：自領域と隣接領域の合計ベクトル数

    hlpod_vals->n_neib_vec = n_neib_vec;

}


void ROM_std_hlpod_set_comm_table_para_subd(
        MONOLIS_COM*	monolis_com,
        MONOLIS_COM*	mono_com_rom,
        const int 		n_neib)
{
    mono_com_rom->comm 		= monolis_mpi_get_global_comm();
    mono_com_rom->comm_size 	= monolis_mpi_get_global_comm_size();
    mono_com_rom->my_rank 	= monolis_mpi_get_global_my_rank();

    mono_com_rom->recv_n_neib = monolis_com->recv_n_neib;
    mono_com_rom->send_n_neib = monolis_com->send_n_neib;

    mono_com_rom->n_internal_vertex = 1;

    mono_com_rom->recv_index		= BB_std_calloc_1d_int(mono_com_rom->recv_index,	n_neib + 1);
    mono_com_rom->recv_item		= BB_std_calloc_1d_int(mono_com_rom->recv_item,  n_neib);	
    mono_com_rom->send_index		= BB_std_calloc_1d_int(mono_com_rom->send_index,	n_neib + 1);
    mono_com_rom->send_item		= BB_std_calloc_1d_int(mono_com_rom->send_item, 	n_neib);	

    mono_com_rom->recv_neib_pe	= BB_std_calloc_1d_int(mono_com_rom->recv_neib_pe,	n_neib);
    mono_com_rom->send_neib_pe	= BB_std_calloc_1d_int(mono_com_rom->send_neib_pe, 	n_neib);	

    for(int i = 0; i < n_neib; i++){
        mono_com_rom->recv_neib_pe[i] = monolis_com->recv_neib_pe[i];
        mono_com_rom->send_neib_pe[i] = monolis_com->send_neib_pe[i];
    }
    
    for(int i = 0; i < n_neib +1; i++){
        mono_com_rom->recv_index[i] = i;
    }

    for(int i = 0; i < n_neib; i++){
        mono_com_rom->recv_item[i] = i + 1;
    }

    for(int i = 0; i < n_neib +1; i++){
        mono_com_rom->send_index[i] = i;
    }

    for(int i = 0; i < n_neib; i++){
        mono_com_rom->send_item[i] = 0;
    }
}

//level1領域の選択された基底(p-adaptive)本数の共有
/*
void ROM_std_hlpod_get_neib_num_modes_para_subd(
        MONOLIS_COM*  	monolis_com,
        HLPOD_VALUES* 	hlpod_vals,
        HLPOD_MAT* 	    hlpod_mat,
        const int       np,
        const int       num_my_modes)
{
    hlpod_mat->num_modes_2nddd = BB_std_calloc_1d_int(hlpod_mat->num_modes_2nddd, np);
    hlpod_mat->num_modes_2nddd[0] = num_my_modes;

    monolis_mpi_update_I(monolis_com, np, 1, hlpod_mat->num_modes_2nddd);

    hlpod_vals->num_modes_max = ROM_BB_findMax (hlpod_mat->num_modes_2nddd, 1 + monolis_com->recv_n_neib);
}
*/

void ROM_std_hlpod_get_neib_num_modes_para_subd(
        MONOLIS_COM*    monolis_com,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_MAT*          hlpod_mat,
        const int       np,
        const int       num_my_modes)
{
    hlpod_mat->num_modes_2nddd = BB_std_calloc_1d_int(hlpod_mat->num_modes_2nddd, np);
    hlpod_mat->num_modes_2nddd[0] = num_my_modes;

    monolis_mpi_update_I(monolis_com, np, 1, hlpod_mat->num_modes_2nddd);

    int num_my_modes_in = num_my_modes;

    monolis_allreduce_I(
        1,
        &num_my_modes_in,
        MONOLIS_MPI_MAX,
        monolis_com->comm);

    hlpod_vals->num_modes_max = num_my_modes_in;
}

void ROM_std_hlpod_get_neib_num_modes_mode_subd(
    MONOLIS_COM*  	mono_com,
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		np,
    const char*     directory)
{
    int max = monolis_com->n_internal_vertex;
    const int n_internal_vertex = monolis_com->n_internal_vertex;
    const int n_internal_sum = hlpod_meta->n_internal_sum;
    
    monolis_allreduce_I(
        1,
        &max,
        MONOLIS_MPI_MAX,
        mono_com->comm);

    int* array = BB_std_calloc_1d_int(array, max * np);
    for(int j = 0; j < n_internal_vertex; j++){
        array[j] = hlpod_mat->num_modes_internal[j];
    }
    monolis_mpi_update_I(mono_com, np, max, array);

    /* 隣接を含めた全頂点(自分＋隣接)分のモード数配列を確保 */
    int index=0;

    int* num_neib_modes_1stdd;
    num_neib_modes_1stdd = BB_std_calloc_1d_int(num_neib_modes_1stdd, n_internal_sum + n_internal_vertex);

    for(int i = 0; i < max * np; i++) {
        if(array[i] != 0){
            num_neib_modes_1stdd[index] = array[i];
            index++;
        }
    }

    /* hlpod_mat->num_modes_1stdd_neib を確保して計算 */
    //hlpod_mat->num_modes_1stdd_neib = BB_std_calloc_1d_int(hlpod_mat->num_modes_1stdd_neib, n_internal_vertex);
    int* num_modes_1stdd_neib = BB_std_calloc_1d_int(num_modes_1stdd_neib, n_internal_vertex);

    /* 1stDDの隣接領域を含めた総基底本数を計算 */
    for (int k = 0; k < n_internal_vertex; k++) {
        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for (int i = iS; i < iE; i++) {
            int item_index = hlpod_meta->item[i];
            int global_id_value = hlpod_meta->subdomain_id[item_index];

            /* グローバルIDに対応する場所を探して基底数を加算 */
            for (int j = 0; j < n_internal_sum + n_internal_vertex; j++){
                if (global_id_value == hlpod_meta->subdomain_id_neib[j]) {
                    //hlpod_mat->num_modes_1stdd_neib[k] += num_neib_modes_1stdd[j];
                    num_modes_1stdd_neib[k] += num_neib_modes_1stdd[j];
                }
            }
        }
    }

    /* 自領域の基底数を加算 (num_neib_modes_1stdd[k] を自分に足す) */
    for(int k = 0; k < n_internal_vertex; k++){
        //hlpod_mat->num_modes_1stdd_neib[k] += num_neib_modes_1stdd[k];
        num_modes_1stdd_neib[k] += num_neib_modes_1stdd[k];
    }

    /* 隣接含めた配列長+1 の配列を確保し、累積和を作成 */
    int total_size = n_internal_sum + n_internal_vertex;
    hlpod_mat->num_modes_1stdd = BB_std_calloc_1d_int(hlpod_mat->num_modes_1stdd, total_size + 1);
    hlpod_mat->num_modes_1stdd[0] = 0;

    for(int i = 0; i < total_size; i++){
        hlpod_mat->num_modes_1stdd[i + 1] =
            hlpod_mat->num_modes_1stdd[i] + num_neib_modes_1stdd[i];
    }

    BB_std_free_1d_int(array, max * np);
    BB_std_free_1d_int(num_neib_modes_1stdd, n_internal_sum + n_internal_vertex);
}
