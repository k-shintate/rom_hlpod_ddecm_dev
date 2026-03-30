#include "core_mpi.h"

//テスト用関数
void lb_neib_test(
    const int myrank,
    const int Ddof,
    LPOD_LB* lpod_lb)
{
    if(myrank == 0){
        printf("meta_list_n_neib\n");       //探す対象1
        for(int m = 0; m < Ddof; m++){
            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                printf("%d  ",lpod_lb->meta_list_n_neib[j][m]);
            }
            printf("\n");
        }
        printf("\n\n\n");

        printf("domain_list_neib\n");   //探す対象2
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int i = 0; i < lpod_lb->max_Ddof; i++){
                printf("%d  ",lpod_lb->domain_list_neib[i][k]);
            }
            printf("\n");
        }
        printf("\n\n\n");
    }

    if(myrank == 0){    //1,2が一致したときtrue
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int m = 0; m < Ddof; m++){
                for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                    if(lpod_lb->overlap_exists[j][m][k] == true){
                        printf("1 ");
                    }
                    else{
                        printf("0 ");
                    }
                }
                printf("\n");
            }
            printf("\n\n\n");
        }
    }
}

void lb_set_node_internal(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];

    int total_num_nodes_local;
    int ndof;

    const int Ddof = lpod_lb->Ddof;
    lpod_lb->n_internal_vertex = BB_std_calloc_1d_int(lpod_lb->n_internal_vertex, Ddof);  
    lpod_lb->total_num_nodes_local = BB_std_calloc_1d_int(lpod_lb->total_num_nodes_local, Ddof);

/*合計内部節点数の計算*/
    for(int d = 0; d < lpod_lb->Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_NODE, lpod_lb->domain_list[d]);
        
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s %d", char_n_internal, &(ndof));
        fscanf(fp, "%d", &(lpod_lb->n_internal_vertex[d]));
        fclose(fp);
    }

    lpod_lb->total_n_internal = 0;
    for(int i = 0; i < lpod_lb->Ddof; i++){
        lpod_lb->total_n_internal += lpod_lb->n_internal_vertex[i];
    }
    const int total_n_internal = lpod_lb->total_n_internal;

   lpod_lb->node_id_internal = BB_std_calloc_1d_int(lpod_lb->node_id_internal, total_n_internal);

    int* node_id_local;
    int iS = 0;
    int iE = lpod_lb->n_internal_vertex[0];
/*global idとの対応*/
    for(int d = 0; d < lpod_lb->Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(lpod_lb->total_num_nodes_local[d]),&(ndof));

        node_id_local = BB_std_calloc_1d_int(node_id_local, lpod_lb->total_num_nodes_local[d]);

        for(int i = 0; i < lpod_lb->n_internal_vertex[d]; i++) {
            fscanf(fp, "%d", &(node_id_local[i]));
        }
        fclose(fp);

        for(int i = 0; i < lpod_lb->n_internal_vertex[d]; i++){
            lpod_lb->node_id_internal[iS + i] = node_id_local[i];
        }
        iS += lpod_lb->n_internal_vertex[d];
        iE += lpod_lb->n_internal_vertex[d + 1];

        BB_std_free_1d_int(node_id_local, lpod_lb->total_num_nodes_local[d]);
    }

    lpod_lb->node_id_internal_for_sort = BB_std_calloc_1d_int(lpod_lb->node_id_internal_for_sort, total_n_internal);    //各MPI計算領域に属するメタメッシュの領域数

    for(int i = 0; i < total_n_internal; i++){
        lpod_lb->node_id_internal_for_sort[i] = i;
    }

    //ROM_BB_bubble_sort_with_id(lpod_lb->node_id_internal, lpod_lb->node_id_internal_for_sort, total_n_internal);
    double t = monolis_get_time_global_sync();
}


//隣接領域リストの作成
void lb_set_neib_list(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int myrank = monolis_mpi_get_global_my_rank();
    lpod_lb->myrank = myrank;


    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    char char_n_internal[BUFFER_SIZE];

    int graph_ndof;
    int ndof;
    int total_num_nodes;
    int recv_num_nodes;
    int n_neib;
    int total_n_neib;

/*meta_domainの隣接領域数を読み込む***/
    char fname_meta_domain_neib[BUFFER_SIZE];
    snprintf(fname_meta_domain_neib, BUFFER_SIZE, "metagraph_parted.0/%s.recv.%d", OUTPUT_FILENAME_DOMAIN_GRAPH, myrank);
    fp = BBFE_sys_read_fopen(fp, fname_meta_domain_neib, directory);
    fscanf(fp, "%d %d", &(lpod_lb->meta_domain_neib), &(n_neib));


    lpod_lb->meta_domain_list_neib = BB_std_calloc_1d_int(lpod_lb->meta_domain_list_neib, lpod_lb->meta_domain_neib);

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        fscanf(fp, "%d", &(lpod_lb->meta_domain_list_neib[i]));
    }

    fclose(fp);
/************************************/

    lpod_lb->domain_list_neib = BB_std_calloc_2d_int(lpod_lb->domain_list_neib, lpod_lb->max_Ddof, lpod_lb->meta_domain_neib);

    int Ddof_neib;
    int tmp;

    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){
        int list = lpod_lb->meta_domain_list_neib[j];

        snprintf(fname, BUFFER_SIZE, 
            "metagraph_parted.0/%s.n_internal.%d", 
            OUTPUT_FILENAME_DOMAIN_GRAPH, list);
            
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s %d", char_n_internal, &(graph_ndof));
        fscanf(fp, "%d", &(Ddof_neib));
        fclose(fp);

        snprintf(fname, BUFFER_SIZE, 
            "metagraph_parted.0/%s.id.%d", 
            OUTPUT_FILENAME_DOMAIN_GRAPH, list);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_n_neib), &(graph_ndof));

        for(int i = 0; i < Ddof_neib; i++){
            fscanf(fp, "%d", &(lpod_lb->domain_list_neib[i][j]));
        }
        for(int i = Ddof_neib; i<lpod_lb->max_Ddof; i++){
            lpod_lb->domain_list_neib[i][j] = -1;   //探索で重複しないよう-1を代入する
        }
        fclose(fp);
    }
    
    lpod_lb->meta_n_neib = BB_std_calloc_1d_int(lpod_lb->meta_n_neib, Ddof);
    int max_meta_n_neib = 0;
        
    for(int j = 0; j < Ddof; j++){
        int list = lpod_lb->domain_list[j];
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, list);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(lpod_lb->meta_n_neib[j]), &(ndof));

        lpod_lb->meta_list_neib = BB_std_calloc_1d_int(lpod_lb->meta_list_neib, lpod_lb->meta_n_neib[j]);

        for(int i = 0; i < lpod_lb->meta_n_neib[j]; i++){
            fscanf(fp, "%d", &(lpod_lb->meta_list_neib[i]));
        }
        fclose(fp);
        if(max_meta_n_neib < lpod_lb->meta_n_neib[j]){
            max_meta_n_neib = lpod_lb->meta_n_neib[j];
        }
    }

//事前にcallocをしておく
    lpod_lb->meta_list_n_neib = BB_std_calloc_2d_int(lpod_lb->meta_list_n_neib, max_meta_n_neib,Ddof);
    lpod_lb->meta_list_num_n_neib = BB_std_calloc_2d_int(lpod_lb->meta_list_num_n_neib, max_meta_n_neib + 1, Ddof);    

//send用追加
    lpod_lb->meta_list_num_n_neib_send = BB_std_calloc_2d_int(lpod_lb->meta_list_num_n_neib_send, max_meta_n_neib + 1, Ddof);

//recv用の配列
    lpod_lb->overlap_exists = std_calloc_3d_bool(lpod_lb->overlap_exists, max_meta_n_neib, lpod_lb->max_Ddof, lpod_lb->meta_domain_neib);

    for(int m = 0; m < Ddof; m++){
        int list = lpod_lb->domain_list[m];
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, list);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(lpod_lb->meta_n_neib[m]), &(recv_num_nodes));

        for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
            fscanf(fp, "%d", &(lpod_lb->meta_list_n_neib[j][m]));
        }

        for(int j = 0; j < lpod_lb->meta_n_neib[m] + 1; j++){
            fscanf(fp, "%d", &(lpod_lb->meta_list_num_n_neib[j][m]));
        }
        fclose(fp);

        //領域グラフからオーバーラップ領域の隣接関係を取得
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                for(int i = 0; i < lpod_lb->max_Ddof; i++){
                    if(lpod_lb->meta_list_n_neib[j][m] == lpod_lb->domain_list_neib[i][k]){
                        lpod_lb->overlap_exists[j][m][k] = true;
                    }
                }
            }
        }
    }
    double t = monolis_get_time_global_sync();
}

/**************************/
void lb_set_node_overlap(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int myrank = monolis_mpi_get_global_my_rank();

    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int total_num_nodes_local;
    int* node_id_local_overlap;

    int tmp = 0;
    int recv_num_nodes_neib = 0;

    int* recv_index;
    recv_index = BB_std_calloc_1d_int(recv_index, lpod_lb->meta_domain_neib);

    int* index_i;
    index_i = BB_std_calloc_1d_int(index_i, lpod_lb->meta_domain_neib);

    int sum_overlap = 0;

    int* meta_recv_local_id;

/*オーバーラップの合計を計算*/
    for(int m = 0; m < Ddof; m++){
        //int sum = 0;
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){

            int iS = 0;
            int iE = 0;

            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                    if(lpod_lb->overlap_exists[j][m][k] == true){
                        iS = lpod_lb->meta_list_num_n_neib[j][m];
                        iE = lpod_lb->meta_list_num_n_neib[j + 1][m];
                        sum_overlap += iE - iS;
                    }

            }

        }
    }
/*******************/
    lpod_lb->sum_overlap = sum_overlap;

    //重複を含む行列サイズでcalloc
    lpod_lb->recv_id = BB_std_calloc_2d_int(lpod_lb->recv_id, lpod_lb->meta_domain_neib, sum_overlap);

    double*** global_x;
    global_x = BB_std_calloc_3d_double(global_x, sum_overlap, lpod_lb->meta_domain_neib, 3);

    double** local_x;

/*****************/

    //lb_neib_test(myrank,Ddof,lpod_lb);

/**************************/

    double** global_x_internal;
    global_x_internal = BB_std_calloc_2d_double(global_x_internal, lpod_lb->total_n_internal, 3);

    int index_x = 0;

    //領域内自由度のうち、対応する領域のみのオーバーラップ領域から読み込む
    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(tmp), &(recv_num_nodes_neib));
        for(int i = 0; i < lpod_lb->meta_n_neib[m]; i++){
            fscanf(fp, "%d", &(tmp));
        }
        for(int i = 0; i < lpod_lb->meta_n_neib[m] + 1; i++){
            fscanf(fp, "%d", &(tmp));
        }
        int val = recv_num_nodes_neib;
        meta_recv_local_id = BB_std_calloc_1d_int(meta_recv_local_id, val);

        for(int i = 0; i < val; i++){
            fscanf(fp, "%d", &(meta_recv_local_id[i]));
        }

        fclose(fp);
        /*読み込み終了*************/

        int iS = 0;
        int iE = 0;
        int ndof;

        //オーバーラップした場合と場合分け
        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        iE = total_num_nodes_local;
        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, iE);
        for(int i = 0; i < iE; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/

        /*node.datの読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp,"%d", &(total_num_nodes_local));

        local_x = BB_std_calloc_2d_double(local_x, total_num_nodes_local, 3);

	    for(int i = 0; i < total_num_nodes_local; i++) {
		    fscanf(fp,"%lf %lf %lf", 
                    &(local_x[i][0]), &(local_x[i][1]), &(local_x[i][2]));

	    }
	    fclose(fp);

        for(int i = 0; i < lpod_lb->n_internal_vertex[m]; i++) {
            for(int n = 0; n < 3 ; n++){
                global_x_internal[index_x][n] = local_x[i][n];
            }
            index_x++;
        }

        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                if(lpod_lb->overlap_exists[j][m][k] == true){
                    int list = lpod_lb->meta_list_n_neib[j][m];

                    //local_idの読み込み
                    index_i[k] = 0;

                    iS = lpod_lb->meta_list_num_n_neib[j][m];
                    iE = lpod_lb->meta_list_num_n_neib[j + 1][m];

                    for(int i = iS; i < iE; i++){
                        int index = meta_recv_local_id[i];
                        lpod_lb->recv_id[k][index_i[k] + recv_index[k]] = node_id_local_overlap[index];
                        for(int n = 0; n < 3 ; n++){
                            global_x[index_i[k] + recv_index[k]][k][n] = local_x[index][n];
                        }
                        index_i[k]++;
                    }
                    recv_index[k] += iE-iS;
                }
            }
        }

        BB_std_free_1d_int(meta_recv_local_id, val);
        BB_std_free_1d_int(node_id_local_overlap, iE);
        BB_std_free_2d_double(local_x, total_num_nodes_local, 3);
    }


/***********************/
    int** id_for_sort;
    id_for_sort = BB_std_calloc_2d_int(id_for_sort, lpod_lb->meta_domain_neib, sum_overlap);
    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        for(int j = 0; j < sum_overlap; j++){
            id_for_sort[i][j] = j;
        }
    }

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        ROM_BB_bubble_sort_with_id(lpod_lb->recv_id[i], id_for_sort[i], recv_index[i]);
    }

    int return_val;
    int* num_search;
    num_search = BB_std_calloc_1d_int(num_search, lpod_lb->meta_domain_neib + 1);

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        std_bubble_sort_remove_duplicate_node_return(lpod_lb->recv_id[i], recv_index[i], &return_val);
        num_search[i + 1] = num_search[i] + recv_index[i] - return_val;
    }

    lpod_lb->local_id = BB_std_calloc_1d_int(lpod_lb->local_id, lpod_lb->total_n_internal + sum_overlap);
    lpod_lb->local_id_for_sort = BB_std_calloc_1d_int(lpod_lb->local_id_for_sort, lpod_lb->total_n_internal + sum_overlap);

    for(int i = 0; i < lpod_lb->total_n_internal; i++){   
        int val = lpod_lb->node_id_internal[i];
        lpod_lb->local_id[i] = val;
        lpod_lb->local_id_for_sort[i] = i;
    }

/*逐次版との同じidの並びにするための操作 過剰なため削除可能*/
    int* overlap_id_for_sort;
    int* overlap_id;
    overlap_id_for_sort = BB_std_calloc_1d_int(overlap_id_for_sort, sum_overlap);
    overlap_id = BB_std_calloc_1d_int(overlap_id, sum_overlap);
    
    int index = 0;
    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(lpod_lb->recv_id[j][i] != -1){
                overlap_id[index] = lpod_lb->recv_id[j][i];
                overlap_id_for_sort[index] = index + lpod_lb->total_n_internal;
                index++;
            }
        }
    }

    ROM_BB_bubble_sort_with_id(overlap_id, overlap_id_for_sort, index);

    int ie = index;
    index = lpod_lb->total_n_internal;
    for(int i = 0; i < ie; i++){
        lpod_lb->local_id[index] = overlap_id[i];
        lpod_lb->local_id_for_sort[index] = index;
        index++;
    }

    lpod_lb->total_num_nodes_lb = lpod_lb->total_n_internal + sum_overlap;

    ROM_BB_bubble_sort_with_id(lpod_lb->local_id, lpod_lb->local_id_for_sort, lpod_lb->total_num_nodes_lb);
    return_val = 0;
    std_bubble_sort_remove_duplicate_node_return(lpod_lb->local_id, lpod_lb->total_num_nodes_lb, &return_val);

/*recv配列の作成***********/

    LB_write_node_recv(recv_index, num_search,lpod_lb, directory);
/**************************/

    int ndof = 1;

/*id配列の作成***********/
    LB_write_node_id(ndof,ie,overlap_id,num_search,lpod_lb,directory);

/**************************/
    double** global_x_overlap;
    global_x_overlap = BB_std_calloc_2d_double(global_x_overlap, sum_overlap, 3);
    
    index = 0;
    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(lpod_lb->recv_id[j][i] != -1){
                int index1 = id_for_sort[j][i];
                for(int k = 0; k < 3 ; k++){
                    global_x_overlap[index][k] = global_x[index1][j][k];
                }
                index++;
            }
        }
    }

    //nodeの書き込み
    LB_write_node(
        ie,lpod_lb, num_search, global_x_internal, global_x_overlap, overlap_id_for_sort,directory);

    //internal_vertexの書き出し
    LB_write_node_internal(
        ndof, lpod_lb->total_n_internal, lpod_lb, directory);

    BB_std_calloc_3d_double(global_x, sum_overlap, lpod_lb->meta_domain_neib, 3);
    BB_std_calloc_2d_double(global_x_internal, lpod_lb->total_n_internal, 3);
    BB_std_calloc_2d_double(global_x_overlap, sum_overlap, 3);
    BB_std_calloc_1d_int(num_search, lpod_lb->meta_domain_neib + 1);
    BB_std_calloc_1d_int(recv_index, lpod_lb->meta_domain_neib);
    BB_std_calloc_1d_int(index_i, lpod_lb->meta_domain_neib);
    BB_std_free_2d_int(id_for_sort, lpod_lb->meta_domain_neib, sum_overlap);
    BB_std_free_1d_int(overlap_id, sum_overlap);
    BB_std_free_1d_int(overlap_id_for_sort, sum_overlap);

    double t = monolis_get_time_global_sync();
}


void lb_set_node_overlap_without_write(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int myrank = monolis_mpi_get_global_my_rank();

    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int total_num_nodes_local;
    int* node_id_local_overlap;

    int tmp = 0;
    int recv_num_nodes_neib = 0;

    int* recv_index;
    recv_index = BB_std_calloc_1d_int(recv_index, lpod_lb->meta_domain_neib);

    int* index_i;
    index_i = BB_std_calloc_1d_int(index_i, lpod_lb->meta_domain_neib);

    int sum_overlap = 0;

    int* meta_recv_local_id;

/*オーバーラップの合計を計算*/
    for(int m = 0; m < Ddof; m++){
        //int sum = 0;
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){

            int iS = 0;
            int iE = 0;

            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                    if(lpod_lb->overlap_exists[j][m][k] == true){
                        iS = lpod_lb->meta_list_num_n_neib[j][m];
                        iE = lpod_lb->meta_list_num_n_neib[j + 1][m];
                        sum_overlap += iE - iS;
                    }

            }

        }
    }
/*******************/
    lpod_lb->sum_overlap = sum_overlap;

    //重複を含む行列サイズでcalloc
    lpod_lb->recv_id = BB_std_calloc_2d_int(lpod_lb->recv_id, lpod_lb->meta_domain_neib, sum_overlap);

    double*** global_x;
    global_x = BB_std_calloc_3d_double(global_x, sum_overlap, lpod_lb->meta_domain_neib, 3);

    double** local_x;

/*****************/

    //lb_neib_test(myrank,Ddof,lpod_lb);

/**************************/

    double** global_x_internal;
    global_x_internal = BB_std_calloc_2d_double(global_x_internal, lpod_lb->total_n_internal, 3);

    int index_x = 0;

    //領域内自由度のうち、対応する領域のみのオーバーラップ領域から読み込む
    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(tmp), &(recv_num_nodes_neib));
        for(int i = 0; i < lpod_lb->meta_n_neib[m]; i++){
            fscanf(fp, "%d", &(tmp));
        }
        for(int i = 0; i < lpod_lb->meta_n_neib[m] + 1; i++){
            fscanf(fp, "%d", &(tmp));
        }
        int val = recv_num_nodes_neib;
        meta_recv_local_id = BB_std_calloc_1d_int(meta_recv_local_id, val);

        for(int i = 0; i < val; i++){
            fscanf(fp, "%d", &(meta_recv_local_id[i]));
        }

        fclose(fp);
        /*読み込み終了*************/

        int iS = 0;
        int iE = 0;
        int ndof;

        //オーバーラップした場合と場合分け
        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        iE = total_num_nodes_local;
        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, iE);
        for(int i = 0; i < iE; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/

        /*node.datの読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp,"%d", &(total_num_nodes_local));

        local_x = BB_std_calloc_2d_double(local_x, total_num_nodes_local, 3);

	    for(int i = 0; i < total_num_nodes_local; i++) {
		    fscanf(fp,"%lf %lf %lf", 
                    &(local_x[i][0]), &(local_x[i][1]), &(local_x[i][2]));

	    }
	    fclose(fp);

        for(int i = 0; i < lpod_lb->n_internal_vertex[m]; i++) {
            for(int n = 0; n < 3 ; n++){
                global_x_internal[index_x][n] = local_x[i][n];
            }
            index_x++;
        }

        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                if(lpod_lb->overlap_exists[j][m][k] == true){
                    int list = lpod_lb->meta_list_n_neib[j][m];

                    //local_idの読み込み
                    index_i[k] = 0;

                    iS = lpod_lb->meta_list_num_n_neib[j][m];
                    iE = lpod_lb->meta_list_num_n_neib[j + 1][m];

                    for(int i = iS; i < iE; i++){
                        int index = meta_recv_local_id[i];
                        lpod_lb->recv_id[k][index_i[k] + recv_index[k]] = node_id_local_overlap[index];
                        for(int n = 0; n < 3 ; n++){
                            global_x[index_i[k] + recv_index[k]][k][n] = local_x[index][n];
                        }
                        index_i[k]++;
                    }
                    recv_index[k] += iE-iS;
                }
            }
        }

        BB_std_free_1d_int(meta_recv_local_id, val);
        BB_std_free_1d_int(node_id_local_overlap, iE);
        BB_std_free_2d_double(local_x, total_num_nodes_local, 3);
    }


/***********************/
    int** id_for_sort;
    id_for_sort = BB_std_calloc_2d_int(id_for_sort, lpod_lb->meta_domain_neib, sum_overlap);
    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        for(int j = 0; j < sum_overlap; j++){
            id_for_sort[i][j] = j;
        }
    }

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        ROM_BB_bubble_sort_with_id(lpod_lb->recv_id[i], id_for_sort[i], recv_index[i]);
    }

    int return_val;
    int* num_search;
    num_search = BB_std_calloc_1d_int(num_search, lpod_lb->meta_domain_neib + 1);

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        std_bubble_sort_remove_duplicate_node_return(lpod_lb->recv_id[i], recv_index[i], &return_val);
        num_search[i + 1] = num_search[i] + recv_index[i] - return_val;
    }

    lpod_lb->local_id = BB_std_calloc_1d_int(lpod_lb->local_id, lpod_lb->total_n_internal + sum_overlap);
    lpod_lb->local_id_for_sort = BB_std_calloc_1d_int(lpod_lb->local_id_for_sort, lpod_lb->total_n_internal + sum_overlap);

    for(int i = 0; i < lpod_lb->total_n_internal; i++){   
        int val = lpod_lb->node_id_internal[i];
        lpod_lb->local_id[i] = val;
        lpod_lb->local_id_for_sort[i] = i;
    }

/*逐次版との同じidの並びにするための操作 過剰なため削除可能*/
    int* overlap_id_for_sort;
    int* overlap_id;
    overlap_id_for_sort = BB_std_calloc_1d_int(overlap_id_for_sort, sum_overlap);
    overlap_id = BB_std_calloc_1d_int(overlap_id, sum_overlap);
    
    int index = 0;
    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(lpod_lb->recv_id[j][i] != -1){
                overlap_id[index] = lpod_lb->recv_id[j][i];
                overlap_id_for_sort[index] = index + lpod_lb->total_n_internal;
                index++;
            }
        }
    }

    ROM_BB_bubble_sort_with_id(overlap_id, overlap_id_for_sort, index);

    int ie = index;
    index = lpod_lb->total_n_internal;
    for(int i = 0; i < ie; i++){
        lpod_lb->local_id[index] = overlap_id[i];
        lpod_lb->local_id_for_sort[index] = index;
        index++;
    }

    lpod_lb->total_num_nodes_lb = lpod_lb->total_n_internal + sum_overlap;

    ROM_BB_bubble_sort_with_id(lpod_lb->local_id, lpod_lb->local_id_for_sort, lpod_lb->total_num_nodes_lb);
    return_val = 0;
    std_bubble_sort_remove_duplicate_node_return(lpod_lb->local_id, lpod_lb->total_num_nodes_lb, &return_val);

    int ndof = 1;

    double** global_x_overlap;
    global_x_overlap = BB_std_calloc_2d_double(global_x_overlap, sum_overlap, 3);
    
    index = 0;
    for(int j = 0; j < lpod_lb->meta_domain_neib; j++){   
        for(int i = 0; i < recv_index[j]; i++){
            if(lpod_lb->recv_id[j][i] != -1){
                int index1 = id_for_sort[j][i];
                for(int k = 0; k < 3 ; k++){
                    global_x_overlap[index][k] = global_x[index1][j][k];
                }
                index++;
            }
        }
    }

    BB_std_calloc_3d_double(global_x, sum_overlap, lpod_lb->meta_domain_neib, 3);
    BB_std_calloc_2d_double(global_x_internal, lpod_lb->total_n_internal, 3);
    BB_std_calloc_2d_double(global_x_overlap, sum_overlap, 3);
    BB_std_calloc_1d_int(num_search, lpod_lb->meta_domain_neib + 1);
    BB_std_calloc_1d_int(recv_index, lpod_lb->meta_domain_neib);
    BB_std_calloc_1d_int(index_i, lpod_lb->meta_domain_neib);
    BB_std_free_2d_int(id_for_sort, lpod_lb->meta_domain_neib, sum_overlap);
    BB_std_free_1d_int(overlap_id, sum_overlap);
    BB_std_free_1d_int(overlap_id_for_sort, sum_overlap);

    double t = monolis_get_time_global_sync();
}


/*send行列の作成***************/
void lb_set_node_send(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int myrank = monolis_mpi_get_global_my_rank();

    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int send_num_nodes;
    
    for(int m = 0; m < Ddof; m++){
        int list = lpod_lb->domain_list[m];
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.send.%d", INPUT_FILENAME_NODE, list);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(lpod_lb->meta_n_neib[m]), &(send_num_nodes));

        for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
            fscanf(fp, "%d", &(lpod_lb->meta_list_n_neib[j][m]));
        }
        for(int j = 0; j < lpod_lb->meta_n_neib[m] + 1; j++){
            fscanf(fp, "%d", &(lpod_lb->meta_list_num_n_neib_send[j][m]));
        }
        fclose(fp);
    }

/*calc_sum_ovelap*/
    int sum_overlap_send = 0;
    for(int m = 0; m < Ddof; m++){
        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){

            int iS = 0;
            int iE = 0;

            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                    if(lpod_lb->overlap_exists[j][m][k] == true){
                        iS = lpod_lb->meta_list_num_n_neib_send[j][m];
                        iE = lpod_lb->meta_list_num_n_neib_send[j + 1][m];
                        sum_overlap_send += iE - iS;
                    }

            }

        }
    }
/*******************/
    
    //重複を含む行列サイズでcalloc
    int** send_id;
    send_id = BB_std_calloc_2d_int(send_id, lpod_lb->meta_domain_neib, sum_overlap_send);

    int* node_id_local_overlap_send;
    int iS = 0;
    int iE = 0;
    int tmp;
    int ndof;
    int total_num_nodes_local;

    int* send_index;
    send_index = BB_std_calloc_1d_int(send_index, lpod_lb->meta_domain_neib);
    int* index_i;
    index_i = BB_std_calloc_1d_int(index_i, lpod_lb->meta_domain_neib);

    int* meta_send_local_id;
    int send_num_nodes_neib;

    for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
        send_index[k] = 0;
    }

    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.send.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d %d", &(tmp), &(send_num_nodes_neib));
        for(int i = 0; i < lpod_lb->meta_n_neib[m]; i++){
            fscanf(fp, "%d", &(tmp));
        }
        for(int i = 0; i < lpod_lb->meta_n_neib[m] + 1; i++){
            fscanf(fp, "%d", &(tmp));
        }
        int val = send_num_nodes_neib;
        meta_send_local_id = BB_std_calloc_1d_int(meta_send_local_id, val);

        for(int i = 0; i < val; i++){
            fscanf(fp, "%d", &(meta_send_local_id[i]));
        }

        fclose(fp);
        /*読み込み終了*************/

        //オーバーラップした場合と場合分け
        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        node_id_local_overlap_send = BB_std_calloc_1d_int(node_id_local_overlap_send, total_num_nodes_local);
        for(int i = 0; i < total_num_nodes_local; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap_send[i]));
        }
        fclose(fp);
        /**********************/

        for(int k = 0; k < lpod_lb->meta_domain_neib; k++){
            for(int j = 0; j < lpod_lb->meta_n_neib[m]; j++){
                if(lpod_lb->overlap_exists[j][m][k] == true){
                    int list = lpod_lb->meta_list_n_neib[j][m];

                    //local_idの読み込み
                    index_i[k] = 0;

                    iS = lpod_lb->meta_list_num_n_neib_send[j][m];
                    iE = lpod_lb->meta_list_num_n_neib_send[j + 1][m];
     
                    for(int i = iS; i < iE; i++){
                        int index = meta_send_local_id[i];
                        send_id[k][index_i[k] + send_index[k]] = node_id_local_overlap_send[index];
                        index_i[k]++;
                    }
                    send_index[k] += iE-iS;

                }

            }
        }

        BB_std_free_1d_int(meta_send_local_id, val);
        BB_std_free_1d_int(node_id_local_overlap_send, total_num_nodes_local);
    }

/*sendについてソートと書き出し**********************/
    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        std_bubble_sort(send_id[i], send_index[i]);
    }

    int return_val;
    int* num_send;
    num_send = BB_std_calloc_1d_int(num_send, lpod_lb->meta_domain_neib + 1);
    num_send[0] = 0;

    for(int i = 0; i < lpod_lb->meta_domain_neib; i++){
        std_bubble_sort_remove_duplicate_node_return(send_id[i], send_index[i], &return_val);
        num_send[i + 1] = num_send[i] + send_index[i] - return_val;
    }

    LB_write_node_send(send_index, num_send, send_id, lpod_lb, directory);

    BB_std_free_2d_int(send_id, lpod_lb->meta_domain_neib, sum_overlap_send);
    BB_std_free_1d_int(num_send, lpod_lb->meta_domain_neib + 1);
    BB_std_free_1d_int(send_index, lpod_lb->meta_domain_neib);
    BB_std_free_1d_int(index_i, lpod_lb->meta_domain_neib);
    double t = monolis_get_time_global_sync();
}



void lb_set_D_bc(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int sum_overlap = lpod_lb->sum_overlap;
    const int myrank = monolis_mpi_get_global_my_rank();

    FILE* fp;
    char id[BUFFER_SIZE];
    char fname[BUFFER_SIZE];

    bool* D_bc_exists;
    double* imposed_D_val;
    int* D_bc_index;

    int* node_id_local_overlap;

    lpod_lb->total_num_D_bc = read_total_num_D_bc(Ddof, lpod_lb->domain_list, directory);


    int* global_D_bc_index;
    double* global_imposed_D_val;
    global_D_bc_index = BB_std_calloc_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    global_imposed_D_val = BB_std_calloc_1d_double(global_imposed_D_val, lpod_lb->total_num_D_bc);

    int iS = 0;
    int iE = 0;

    int index_D = 0;
    int total_num_nodes_local;  int num_D_bcs;  int block_size;
    int node_id;  int block_id;  double val;
    int ndof;

    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*D_bcの読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_D_BC, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

        imposed_D_val = BB_std_calloc_1d_double(imposed_D_val, num_D_bcs);
        D_bc_index = BB_std_calloc_1d_int(D_bc_index, num_D_bcs);

        double value = 0;
    	for(int i = 0; i < num_D_bcs; i++) {
            fscanf(fp, "%d %d %lf", &(node_id), &(block_id), &(value));
            int index = (block_size)*node_id + block_id;
		    imposed_D_val[i] = value;
            D_bc_index[i] = index;    //local node id
	    }
        fclose(fp);
        /**********************/

        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, total_num_nodes_local);
        for(int i = 0; i < total_num_nodes_local; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/

        for(int i = 0; i < num_D_bcs; i++) {
            global_imposed_D_val[index_D] = imposed_D_val[i];
            global_D_bc_index[index_D] = node_id_local_overlap[D_bc_index[i]];
            index_D++;
        }
	BB_std_free_1d_double(imposed_D_val, num_D_bcs);
        BB_std_free_1d_int(node_id_local_overlap, total_num_nodes_local);
    }

    int* id_for_sort_D_bc;
    id_for_sort_D_bc = BB_std_calloc_1d_int(id_for_sort_D_bc, index_D);
    for(int i = 0; i < index_D; i++) {
        id_for_sort_D_bc[i] = i;
    }
    ROM_BB_bubble_sort_with_id(global_D_bc_index, id_for_sort_D_bc, index_D);
    int return_val = 0;
    std_bubble_sort_remove_duplicate_node_return(global_D_bc_index, index_D, &return_val);

    //D_bcの書き込み
    LB_write_D_bc(global_D_bc_index, return_val, global_imposed_D_val, lpod_lb, directory);

    BB_std_free_1d_int(id_for_sort_D_bc, index_D);
    BB_std_free_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    BB_std_free_1d_double(global_imposed_D_val, lpod_lb->total_num_D_bc);
    
    double t = monolis_get_time_global_sync();
}


void lb_set_D_bc_p(
    const int       nd1,
    const char*     directory,
    const char*     label,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int sum_overlap = lpod_lb->sum_overlap;
    const int myrank = monolis_mpi_get_global_my_rank();

    FILE* fp;
    char id[BUFFER_SIZE];
    char fname[BUFFER_SIZE];

    bool* D_bc_exists;
    double* imposed_D_val;
    int* D_bc_index;

    int* node_id_local_overlap;

    lpod_lb->total_num_D_bc = read_total_num_D_bc_p(Ddof, lpod_lb->domain_list, label, directory);


    int* global_D_bc_index;
    double* global_imposed_D_val;
    global_D_bc_index = BB_std_calloc_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    global_imposed_D_val = BB_std_calloc_1d_double(global_imposed_D_val, lpod_lb->total_num_D_bc);

    int iS = 0;
    int iE = 0;

    int index_D = 0;
    int total_num_nodes_local;  int num_D_bcs;  int block_size;
    int node_id;  int block_id;  double val;
    int ndof;

    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*D_bcの読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.dat.%d", label, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

        imposed_D_val = BB_std_calloc_1d_double(imposed_D_val, num_D_bcs);
        D_bc_index = BB_std_calloc_1d_int(D_bc_index, num_D_bcs);

        double value = 0;
    	for(int i = 0; i < num_D_bcs; i++) {
            fscanf(fp, "%d %d %lf", &(node_id), &(block_id), &(value));
            int index = (block_size)*node_id + block_id;
		    imposed_D_val[i] = value;
            D_bc_index[i] = index;    //local node id
	    }
        fclose(fp);
        /**********************/

        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, total_num_nodes_local);
        for(int i = 0; i < total_num_nodes_local; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/

        for(int i = 0; i < num_D_bcs; i++) {
            global_imposed_D_val[index_D] = imposed_D_val[i];
            global_D_bc_index[index_D] = node_id_local_overlap[D_bc_index[i]];
            index_D++;
        }
	    BB_std_free_1d_double(imposed_D_val, num_D_bcs);
        BB_std_free_1d_int(node_id_local_overlap, total_num_nodes_local);
    }

    int* id_for_sort_D_bc;
    id_for_sort_D_bc = BB_std_calloc_1d_int(id_for_sort_D_bc, index_D);
    for(int i = 0; i < index_D; i++) {
        id_for_sort_D_bc[i] = i;
    }
    ROM_BB_bubble_sort_with_id(global_D_bc_index, id_for_sort_D_bc, index_D);
    int return_val = 0;
    std_bubble_sort_remove_duplicate_node_return(global_D_bc_index, index_D, &return_val);

    //D_bcの書き込み
    LB_write_D_bc_p(global_D_bc_index, return_val, global_imposed_D_val, lpod_lb, label, directory);

    BB_std_free_1d_int(id_for_sort_D_bc, index_D);
    BB_std_free_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    BB_std_free_1d_double(global_imposed_D_val, lpod_lb->total_num_D_bc);
    
    double t = monolis_get_time_global_sync();
}


void lb_set_D_bc_v(
    const int       nd1,
    const char*     directory,
    const int       dof,
    const char*     label,
    LPOD_LB*        lpod_lb)
{
    const int Ddof = lpod_lb->Ddof;
    const int sum_overlap = lpod_lb->sum_overlap;
    const int myrank = monolis_mpi_get_global_my_rank();
    int tmp;

    FILE* fp;
    char id[BUFFER_SIZE];
    char fname[BUFFER_SIZE];

    bool* D_bc_exists;
    double** imposed_D_val;
    int* D_bc_index;

    int* node_id_local_overlap;

    tmp = read_total_num_D_bc_v(Ddof, lpod_lb->domain_list, directory);

    lpod_lb->total_num_D_bc = tmp / dof;


    int* global_D_bc_index;
    double** global_imposed_D_val;
    global_D_bc_index = BB_std_calloc_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    global_imposed_D_val = BB_std_calloc_2d_double(global_imposed_D_val, lpod_lb->total_num_D_bc, dof);

    int iS = 0;
    int iE = 0;

    int index_D = 0;
    int total_num_nodes_local;  int num_D_bcs;  int block_size;
    int node_id;  int block_id;  double val;
    int ndof;

    for(int m = 0; m < Ddof; m++){
        int mydomain = lpod_lb->domain_list[m];

        /*D_bcの読み込み*****************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.dat.%d", label ,mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp,"%d %d", &(tmp), &(block_size));

        num_D_bcs = tmp/ dof;

        imposed_D_val = BB_std_calloc_2d_double(imposed_D_val, num_D_bcs, dof);
        D_bc_index = BB_std_calloc_1d_int(D_bc_index, num_D_bcs);

        double value = 0;

    	for(int i = 0; i < num_D_bcs; i++) {
            for(int j = 0; j < dof; j++){
                fscanf(fp, "%d %d %lf", &(node_id), &(block_id), &(value));
                imposed_D_val[i][j] = value;
            }
            D_bc_index[i] = node_id;    //local node id
	    }
        /**********************/

        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, total_num_nodes_local);
        for(int i = 0; i < total_num_nodes_local; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/

        for(int i = 0; i < num_D_bcs; i++) {
            for(int j = 0; j < dof; j++){
                global_imposed_D_val[index_D][j] = imposed_D_val[i][j];
            }
            global_D_bc_index[index_D] = node_id_local_overlap[D_bc_index[i]];
            index_D++;
        }
	    BB_std_free_2d_double(imposed_D_val, num_D_bcs, dof);
        BB_std_free_1d_int(node_id_local_overlap, total_num_nodes_local);
    }

    int* id_for_sort_D_bc;
    id_for_sort_D_bc = BB_std_calloc_1d_int(id_for_sort_D_bc, lpod_lb->total_num_D_bc);
    for(int i = 0; i < lpod_lb->total_num_D_bc; i++) {
        id_for_sort_D_bc[i] = i;
    }
        
    int return_val = 0;

    char fname_out[BUFFER_SIZE];
    snprintf(fname_out, BUFFER_SIZE, "merged_graph/%s.dat.%d", label ,lpod_lb->myrank);    
    fp = BBFE_sys_write_fopen(fp, fname_out, directory);
    int count = 0;

    fprintf(fp,"%d %d\n", (lpod_lb->total_num_D_bc - return_val)*dof, block_size);

    for(int i = 0; i < lpod_lb->total_num_D_bc; i++){
        if(global_D_bc_index[i] != -1){
            int index = global_D_bc_index[i];
            int value = ROM_BB_binarySearch(lpod_lb->local_id, index, lpod_lb->total_num_nodes_lb);
            int local_id = lpod_lb->local_id_for_sort[value];
            
            for(int j = 0; j < dof; j++){
                fprintf(fp,"%d %d %e\n", local_id, j, global_imposed_D_val[id_for_sort_D_bc[i]][j]);
            }
        }
    }


    double t = monolis_get_time_global_sync();

    BB_std_free_1d_int(id_for_sort_D_bc, lpod_lb->total_num_D_bc);
    BB_std_free_1d_int(global_D_bc_index, lpod_lb->total_num_D_bc);
    BB_std_free_2d_double(global_imposed_D_val, lpod_lb->total_num_D_bc, dof);

    t = monolis_get_time_global_sync();
}



void lb_set_elem_internal(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    const int Ddof = lpod_lb->Ddof;
    const int myrank = monolis_mpi_get_global_my_rank();

    int ndof;

    int* total_num_elems_local;
    total_num_elems_local = BB_std_calloc_1d_int(total_num_elems_local, Ddof);  

    int* n_internal_vertex_elem;
    n_internal_vertex_elem = BB_std_calloc_1d_int(n_internal_vertex_elem, Ddof);  

//合計内部節点数の計算
    int total_n_internal = 0;
    int total_n_overlap = 0;
    for(int d = 0; d < Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);
        
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s %d", char_n_internal, &(ndof));
        fscanf(fp, "%d", &(n_internal_vertex_elem[d]));
        fclose(fp);
    }

    for(int d = 0; d < Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_elems_local[d]),&(ndof));
    }

    for(int i = 0; i < Ddof; i++){
        total_n_internal += n_internal_vertex_elem[i];
    }

    for(int i = 0; i < Ddof; i++){
        total_n_overlap += total_num_elems_local[i] - n_internal_vertex_elem[i];
    }


    int* elem_id_bef_sort;
    elem_id_bef_sort = BB_std_calloc_1d_int(elem_id_bef_sort, total_n_internal);

    int* elem_id_overlap;
    elem_id_overlap = BB_std_calloc_1d_int(elem_id_overlap, total_n_overlap);

    /*内部idの読み込みと足し込み*/

    int* elem_id_local;
    int* elem_id_local_overlap;

    int iS_i = 0;
    int iS_o = 0;
    int iE = n_internal_vertex_elem[0];

    for(int d = 0; d < Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_elems_local[d]),&(ndof));

        elem_id_local = BB_std_calloc_1d_int(elem_id_local, n_internal_vertex_elem[d]);
        elem_id_local_overlap = BB_std_calloc_1d_int(elem_id_local_overlap, total_num_elems_local[d] - n_internal_vertex_elem[d]);

        for(int i = 0; i < n_internal_vertex_elem[d]; i++) {
            fscanf(fp, "%d", &(elem_id_local[i]));
        }

        int index = 0;
        for(int i = n_internal_vertex_elem[d]; i < total_num_elems_local[d]; i++){
            fscanf(fp, "%d", &(elem_id_local_overlap[index]));
            index++;
        }
        fclose(fp);

        for(int i = 0; i < n_internal_vertex_elem[d]; i++){
            elem_id_bef_sort[iS_i + i] = elem_id_local[i];
        }
        iS_i += n_internal_vertex_elem[d];

        index = 0;
        for(int i = n_internal_vertex_elem[d]; i < total_num_elems_local[d]; i++){
            elem_id_overlap[iS_o + index] = elem_id_local_overlap[index];
            index++;
        }

        iS_o += total_num_elems_local[d] - n_internal_vertex_elem[d];
    }

    /*elem.datの書き出しのための読み込み*/
    int local_node_dof;

    int** global_conn;
    global_conn = BB_std_calloc_2d_int(global_conn, total_n_internal + total_n_overlap, 8);

    int** global_conn_internal;
    int** global_conn_overlap;

    int total_num_elems;

    iS_i = 0;
    iS_o = 0;
    int total_num_nodes_local;
    int* node_id_local_overlap;

    for(int d = 0; d < Ddof; d++){
        int mydomain = lpod_lb->domain_list[d];
        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        iE = total_num_nodes_local;
        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, iE);
        for(int i = 0; i < iE; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp,"%d %d", &(total_num_elems), &(local_node_dof));
        

        global_conn_internal = BB_std_calloc_2d_int(global_conn_internal, n_internal_vertex_elem[d], 8);
        global_conn_overlap = BB_std_calloc_2d_int(global_conn_overlap, total_num_elems_local[d] - n_internal_vertex_elem[d], 8);

        for(int e = 0; e < n_internal_vertex_elem[d]; e++) {
            for(int i = 0; i < local_node_dof ; i++){
                fscanf(fp, "%d", &(global_conn_internal[e][i]));
            }
        }

        int index = 0;
        for(int e = n_internal_vertex_elem[d]; e < total_num_elems_local[d]; e++){
            for(int i = 0; i < local_node_dof ; i++){
                fscanf(fp, "%d", &(global_conn_overlap[index][i]));
            }
            index++;
        }
        fclose(fp);

        for(int e = 0; e < n_internal_vertex_elem[d]; e++){
            for(int i = 0; i < local_node_dof ; i++) {
                index = global_conn_internal[e][i];
                global_conn[iS_i + e][i] = node_id_local_overlap[index];
            }
        }
        iS_i += n_internal_vertex_elem[d];

        int index1;
        index = 0;
        for(int e = n_internal_vertex_elem[d]; e < total_num_elems_local[d]; e++){
            int index_row = total_n_internal + iS_o + index;
            for(int i = 0; i < local_node_dof ; i++){
                index1 = global_conn_overlap[index][i];
                global_conn[index_row][i] = node_id_local_overlap[index1];
            }
            index++;
        }
        iS_o += total_num_elems_local[d] - n_internal_vertex_elem[d];
    }

//読み込み後にバブルソートを行う

//idの代入
    int* id_overlap;
    id_overlap = BB_std_calloc_1d_int(id_overlap, total_n_overlap);

    for(int i = 0; i < total_n_overlap; i++){
        id_overlap[i] = i + total_n_internal;
    }

    ROM_BB_bubble_sort_with_id(elem_id_overlap, id_overlap, total_n_overlap);

    int* elem_internal_id_savings;
    elem_internal_id_savings = BB_std_calloc_1d_int(elem_internal_id_savings, total_n_overlap);
    int num_internal_elem_add;

    int* id_internal_savings;
    id_internal_savings = BB_std_calloc_1d_int(id_internal_savings, total_n_overlap);

//overlap領域同士での重なりを判定する
    std_bubble_sort_remove_duplicate_elem_with_id(
        elem_id_overlap, elem_internal_id_savings, 
        &num_internal_elem_add, total_n_overlap,
        id_overlap, id_internal_savings);

//内部要素と判定された数
    const int total_num_internal_bef_sort = total_n_internal + num_internal_elem_add;

//internalへ足しこみをおこなう
    int* sum_internal_global_id;
    sum_internal_global_id = BB_std_calloc_1d_int(sum_internal_global_id, total_num_internal_bef_sort);
    int* sum_internal_conn_id;
    sum_internal_conn_id = BB_std_calloc_1d_int(sum_internal_conn_id, total_num_internal_bef_sort);

//内部要素の足しこみ
    for(int i = 0; i < total_n_internal; i++){
        sum_internal_global_id[i] = elem_id_bef_sort[i];    //内部要素のglobal id
        sum_internal_conn_id[i] = i;
    }

//overlap要素のうち内部要素と判定された要素の足しこみ
    for(int i = 0; i < num_internal_elem_add; i++){ //overlap要素の内、内部要素と判定されたglobal id
        sum_internal_global_id[i + total_n_internal] = elem_internal_id_savings[i];
        sum_internal_conn_id[i + total_n_internal] = id_internal_savings[i];
    }

//バブルソートを行う
//internalに関して
    ROM_BB_bubble_sort_with_id(
        sum_internal_global_id, sum_internal_conn_id, total_num_internal_bef_sort);
    int return_val_internal = 0;
    std_count_negative_val(
        sum_internal_global_id,  total_num_internal_bef_sort, &return_val_internal);

//overlapに関して
    ROM_BB_bubble_sort_with_id(
        elem_id_overlap, id_overlap, total_n_overlap);
    int return_val_overlap = 0;
    std_count_negative_val(
        elem_id_overlap, total_n_overlap, &return_val_overlap);

//idの書き出し
    LB_write_elem_id(
        ndof, 
        total_n_internal + total_n_overlap + num_internal_elem_add - return_val_internal - return_val_overlap,
        total_n_internal + num_internal_elem_add,
        total_n_overlap,
        sum_internal_global_id,
        elem_id_overlap,
        lpod_lb,
        directory);

//elem.datの書き出し
    LB_write_elem(
        local_node_dof, 
        total_n_internal + num_internal_elem_add,
        total_n_internal + total_n_overlap + num_internal_elem_add - return_val_internal - return_val_overlap,
        total_n_overlap,
        elem_id_overlap, 
        id_overlap,
        sum_internal_global_id, 
        global_conn,
        sum_internal_conn_id,
        lpod_lb,
        directory);

//internal_vertexの書き出し
    LB_write_elem_internal(
        ndof,
        total_n_internal + num_internal_elem_add, 
        lpod_lb, 
        directory);

}


//gedatsuの最新コミットに対応
//internal→内部要素 (従来定義) + マスター要素を含む連結後の内部要素 (従来定義)
//overlap→overlapする節点を含む要素 マスター・スレイブ判定していない (FEM解析内で判定するため)
void lb_set_elem_internal_2(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    const int Ddof = lpod_lb->Ddof;
    int ndof;
    const int myrank = monolis_mpi_get_global_my_rank();

    int* total_num_elems_local;
    total_num_elems_local = BB_std_calloc_1d_int(total_num_elems_local, Ddof);  

    int* n_internal_vertex_elem;
    n_internal_vertex_elem = BB_std_calloc_1d_int(n_internal_vertex_elem, Ddof);  

//合計内部要素数の計算
    int total_n_internal = 0;
    int total_n_overlap = 0;
    int tmp;
    for(int d = 0; d < Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp, "%d %d", &(n_internal_vertex_elem[d]), &(tmp));
        //printf("%d\n",n_internal_vertex_elem[d]);
        fclose(fp);
    }

    for(int i = 0; i < Ddof; i++){
        total_n_internal += n_internal_vertex_elem[i];
    }

    int* elem_id_bef_sort;
    elem_id_bef_sort = BB_std_calloc_1d_int(elem_id_bef_sort, total_n_internal);

    int* elem_id_overlap;
    elem_id_overlap = BB_std_calloc_1d_int(elem_id_overlap, total_n_overlap);

    /*内部idの読み込みと足し込み*/
    int* elem_id_local;
    int* elem_id_local_overlap;

    int iS_i = 0;
    int iS_o = 0;
    int iE = n_internal_vertex_elem[0];

    for(int d = 0; d < Ddof; d++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_elems_local[d]),&(ndof));

        elem_id_local = BB_std_calloc_1d_int(elem_id_local, n_internal_vertex_elem[d]);
        elem_id_local_overlap = BB_std_calloc_1d_int(elem_id_local_overlap, total_num_elems_local[d] - n_internal_vertex_elem[d]);

        for(int i = 0; i < n_internal_vertex_elem[d]; i++) {
            fscanf(fp, "%d", &(elem_id_local[i]));
        }

        fclose(fp);

        for(int i = 0; i < n_internal_vertex_elem[d]; i++){
            elem_id_bef_sort[iS_i + i] = elem_id_local[i];
        }
        iS_i += n_internal_vertex_elem[d];

    }

    /*elem.datの書き出しのための読み込み*/
    int local_node_dof;

    int** global_conn;
    global_conn = BB_std_calloc_2d_int(global_conn, total_n_internal, 8);

    int** global_conn_internal;
    int** global_conn_overlap;

    //int total_num_elems;

    iS_i = 0;
    iS_o = 0;
    int total_num_nodes_local;
    int* node_id_local_overlap;

    for(int d = 0; d < Ddof; d++){
        int mydomain = lpod_lb->domain_list[d];
        /*idの読み込み**********/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, mydomain);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

        iE = total_num_nodes_local;
        node_id_local_overlap = BB_std_calloc_1d_int(node_id_local_overlap, iE);
        for(int i = 0; i < iE; i++) {
            fscanf(fp, "%d", &(node_id_local_overlap[i]));
        }
        fclose(fp);
        /**********************/
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM, lpod_lb->domain_list[d]);

        fp = BBFE_sys_read_fopen(fp, fname, directory);

        fscanf(fp,"%d %d", &(tmp), &(local_node_dof));

        global_conn_internal = BB_std_calloc_2d_int(global_conn_internal, n_internal_vertex_elem[d], 8);

        for(int e = 0; e < n_internal_vertex_elem[d]; e++) {
            for(int i = 0; i < local_node_dof ; i++){
                fscanf(fp, "%d", &(global_conn_internal[e][i]));
            }
        }

        fclose(fp);

        int index = 0;
        for(int e = 0; e < n_internal_vertex_elem[d]; e++){
            for(int i = 0; i < local_node_dof ; i++) {
                index = global_conn_internal[e][i];
                global_conn[iS_i + e][i] = node_id_local_overlap[index];
            }
        }
        iS_i += n_internal_vertex_elem[d];
    }


    //gedatsu変更後追加した関数
    int* id_internal;
    id_internal = BB_std_calloc_1d_int(id_internal, total_n_internal);

    for(int i = 0; i < total_n_internal; i++){
        id_internal[i] = i;
    }

    ROM_BB_bubble_sort_with_id(elem_id_bef_sort, id_internal, total_n_internal);
    int return_val = 0;
    std_bubble_sort_remove_duplicate_node_return(elem_id_bef_sort, total_n_internal, &return_val);


    const int total_num_internal_bef_sort = total_n_internal;    

    int* sum_internal_global_id;
    sum_internal_global_id = BB_std_calloc_1d_int(sum_internal_global_id, total_num_internal_bef_sort);
    int* sum_internal_conn_id;
    sum_internal_conn_id = BB_std_calloc_1d_int(sum_internal_conn_id, total_num_internal_bef_sort);

    for(int i = 0; i < total_n_internal; i++){
        sum_internal_global_id[i] = elem_id_bef_sort[i];
        sum_internal_conn_id[i] = id_internal[i];
    }

//ファイル書き出し
    LB_write_elem_2(
        local_node_dof,
        total_n_internal - return_val,
        total_n_internal,
        sum_internal_global_id, 
        global_conn,
        sum_internal_conn_id,
        lpod_lb,
        directory);

    LB_write_elem_id_2(
        ndof,
        total_n_internal - return_val,
        total_n_internal,
        sum_internal_global_id,
        lpod_lb,
        directory);

    LB_write_elem_internal(
        ndof,
        total_n_internal, 
        lpod_lb, 
        directory);

    double t = monolis_get_time_global_sync();
}

