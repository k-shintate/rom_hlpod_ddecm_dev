//全ノードで保持するプログラム
//デバッグ用

#include "core.h"

void lb_read_meta_graph(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    char fname_id_graph[BUFFER_SIZE];
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    FILE* fp;

    int total_n_neib;
    int graph_ndof;     //graph_weightの節点自由度
    int** list_n_neib;  //メタメッシュの隣接領域番号リスト
    
    pod_lb->Ddof = BB_std_calloc_1d_int(pod_lb->Ddof, nd1);    //各MPI計算領域に属するメタメッシュの領域数

    //内部領域リストの取得
    pod_lb->max_n_neib = 0;
    for(int j = 0; j < nd1; j++){
        snprintf(fname_n_internal_graph, BUFFER_SIZE, "meta_graph_parted.0/%s.n_internal.%d", OUTPUT_FILENAME_DOMAIN_GRAPH, j);
        
        fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
        fscanf(fp, "%s %d", char_n_internal, &(graph_ndof));
        fscanf(fp, "%d", &(pod_lb->Ddof[j]));
        fclose(fp);

        printf("n_neib = %d\n", pod_lb->Ddof[j]);

        if(pod_lb->Ddof[j] > pod_lb->max_n_neib){
            pod_lb->max_n_neib = pod_lb->Ddof[j];
        }
    }
    printf("pod_lb->max_n_neib = %d\n", pod_lb->max_n_neib);

    pod_lb->list_n_neib = BB_std_calloc_2d_int(list_n_neib, pod_lb->max_n_neib, nd1);

    for(int j = 0; j < nd1; j++){
        snprintf(fname_id_graph, BUFFER_SIZE, "meta_graph_parted.0/%s.id.%d", OUTPUT_FILENAME_DOMAIN_GRAPH, j);
        fp = BBFE_sys_read_fopen(fp, fname_id_graph, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(total_n_neib), &(graph_ndof));

        for(int i = 0; i < pod_lb->Ddof[j]; i++){
            fscanf(fp, "%d", &(pod_lb->list_n_neib[i][j]));
            printf("list_n_neib = %d\n", pod_lb->list_n_neib[i][j]);
        }
        fclose(fp);
    }
}

void lb_set_node_global_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    char fname_n_internal[BUFFER_SIZE];
    char fname_id_out[BUFFER_SIZE];
    char fname_id_in[BUFFER_SIZE];
    char fname_n_internal_out[BUFFER_SIZE];

	FILE* fp;

    char char_n_internal[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    int* sum_n_overlap;

    int*** node_id_local;
    int** n_internal_vertex;

    int ndof;
    int total_num_nodes_local;

    pod_lb->sum_n_internal = BB_std_calloc_1d_int(pod_lb->sum_n_internal, nd1);
    sum_n_overlap = BB_std_calloc_1d_int(sum_n_overlap, nd1);


    fp = BBFE_sys_read_fopen(fp, INPUT_FILENAME_NODE, directory);
    fscanf(fp, "%d", &(pod_lb->total_num_nodes));
    printf("total_num_nodes = %d\n",pod_lb->total_num_nodes);
    fclose(fp);

    bool** bool_global_id_internal;
    bool** bool_global_id_overlap;

    pod_lb->bool_global_id_internal = BB_std_calloc_2d_bool(pod_lb->bool_global_id_internal, pod_lb->total_num_nodes, nd1);            
    pod_lb->bool_global_id_overlap = BB_std_calloc_2d_bool(pod_lb->bool_global_id_overlap, pod_lb->total_num_nodes, nd1);

    //globalとinternalを対応させるための配列
    pod_lb->global_id = BB_std_calloc_2d_int(pod_lb->global_id, pod_lb->total_num_nodes, nd1);

    int  max_total_num_nodes_local = 0;
    for(int n = 0; n < nd1; n++){
        for(int m = 0; m < pod_lb->Ddof[n]; m++){
            snprintf(fname_id_in, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, pod_lb->list_n_neib[m][n]);

            fp = BBFE_sys_read_fopen(fp, fname_id_in, directory);
            fscanf(fp, "%s", id);
            fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

            if(max_total_num_nodes_local <= total_num_nodes_local)
            {
                max_total_num_nodes_local = total_num_nodes_local;
            }
        }
    }

    node_id_local = ROM_BB_calloc_3d_int(node_id_local, max_total_num_nodes_local, nd1, pod_lb->max_n_neib);
    n_internal_vertex = BB_std_calloc_2d_int(n_internal_vertex, pod_lb->max_n_neib, nd1);

    int local_num_nodes = 0;
    int total_n_internal = 0;
    for(int n = 0; n < nd1; n++){
        int sum_internal = 0;
        int sum_overlap = 0;
        int doublecount = 0;
        
        for(int m = 0; m < pod_lb->Ddof[n]; m++){
            snprintf(fname_n_internal, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_NODE, pod_lb->list_n_neib[m][n]);
        
            fp = BBFE_sys_read_fopen(fp, fname_n_internal, directory);
            fscanf(fp, "%s %d", char_n_internal, &(ndof));
            fscanf(fp, "%d", &(n_internal_vertex[m][n]));
            fclose(fp);

            printf("n_internal_vertex[n] = %d\n", n_internal_vertex[m][n]);

            snprintf(fname_id_in, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, pod_lb->list_n_neib[m][n]);

            fp = BBFE_sys_read_fopen(fp, fname_id_in, directory);
            fscanf(fp, "%s", id);
            fscanf(fp, "%d %d", &(total_num_nodes_local),&(ndof));

            printf("total_num_nodes_local = %d\n", total_num_nodes_local);

            for(int i = 0; i < total_num_nodes_local; i++) {
                fscanf(fp, "%d", &(node_id_local[i][n][m]));
            }
            fclose(fp);

            for(int i = 0; i < n_internal_vertex[m][n]; i++){
                pod_lb->bool_global_id_internal[node_id_local[i][n][m]][n] = true;
            }
            sum_internal += n_internal_vertex[m][n];

            int count = 0;
            for(int i = n_internal_vertex[m][n]; i < total_num_nodes_local; i++){
                if(pod_lb->bool_global_id_overlap[node_id_local[i][n][m]][n] == true){
                    count++;
                }
                pod_lb->bool_global_id_overlap[node_id_local[i][n][m]][n] = true;
            }
            sum_overlap += total_num_nodes_local - n_internal_vertex[m][n] - count;
        }

        for(int i=0; i < pod_lb->total_num_nodes; i++){
            if(pod_lb->bool_global_id_internal[i][n] == true && pod_lb->bool_global_id_overlap[i][n] == true){
                pod_lb->bool_global_id_overlap[i][n] = false;
                doublecount ++;
            }
        }

        pod_lb->sum_n_internal[n] = sum_internal;
        sum_n_overlap[n] = sum_overlap - doublecount;

        /*********ファイル書き込み********/
        //node.dat.idファイル
        snprintf(fname_id_out, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, n);
        
        fp = BBFE_sys_write_fopen(fp, fname_id_out, directory);
        fprintf(fp,"#id\n");
        int local_num_nodes = pod_lb->sum_n_internal[n] + sum_n_overlap[n];
        fprintf(fp,"%d %d\n", local_num_nodes ,ndof);

        int global_id_count = 0;

        for(int i=0; i < pod_lb->total_num_nodes; i++){
            if(pod_lb->bool_global_id_internal[i][n] == true ){
		        fprintf(fp,"%d\n", i);
                pod_lb->global_id[i][n] = global_id_count;
                global_id_count++;    //globalとlocalを対応させる配列
            }
        }

        for(int i=0; i < pod_lb->total_num_nodes; i++){
            if(pod_lb->bool_global_id_overlap[i][n] == true ){
		        fprintf(fp,"%d\n", i);
                pod_lb->global_id[i][n] = global_id_count;
                global_id_count++;
            }
        }
        fclose(fp);

        //node.dat.n_internal.datファイル
        snprintf(fname_n_internal_out, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_NODE, n);
            
        fp = BBFE_sys_write_fopen(fp, fname_n_internal_out, directory);
        fprintf(fp,"#n_internal %d\n", ndof);
        fprintf(fp,"%d", pod_lb->sum_n_internal[n]);

        fclose(fp);
    }

    //メタメッシュの分割メタメッシュに対するid

    for(int n = 0; n < nd1; n++){
        for(int m = 0; m < pod_lb->Ddof[n]; m++){
            int domain_id = pod_lb->list_n_neib[m][n];
            snprintf(fname_n_internal_out, BUFFER_SIZE, "meta_mesh/parted.0/mm_to_Dmm_internal_id.%d", domain_id);    
            fp = BBFE_sys_write_fopen(fp, fname_n_internal_out, directory);
            fprintf(fp,"#internal_id\n");   //内部節点のみのid　gedatsuによるparted.0ファイルとは形式が違う

            for(int i = 0; i < n_internal_vertex[m][n]; i++){
                int index = node_id_local[i][n][m];
                int row_index = pod_lb->global_id[index][n];
                fprintf(fp,"%d\n",row_index);
            }

            fclose(fp);
        }
    }
}


void lb_write_node_local_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    char fname_meta_node_in[BUFFER_SIZE];
    char fname_D_bc_out[BUFFER_SIZE];
    char fname_node_out[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    FILE* fp;

    int num_D_bcs;
    int block_size;
    int ndof;
    bool* D_bc_exists;
    double* imposed_D_val;
    double** global_x;
    int local_num_nodes;

    printf("D_bc, nodeの書き込み\n");

/***global D_bcの読み込み***/
    fp = BBFE_sys_read_fopen(fp, INPUT_FILENAME_D_BC, directory);

    fscanf(fp,"%d %d", &(num_D_bcs), &(block_size));

    D_bc_exists = BB_std_calloc_1d_bool(D_bc_exists, pod_lb->total_num_nodes);            
    imposed_D_val = BB_std_calloc_1d_double(imposed_D_val, pod_lb->total_num_nodes);

    int node_id;  int block_id;  double val;
	for(int i = 0; i < num_D_bcs; i++) {
        fscanf(fp, "%d %d %lf", &(node_id), &(block_id), &(val));

        int index = (block_size)*node_id + block_id;
        D_bc_exists[ index ]   = true;
		imposed_D_val[ index ] = val;
	}
/***********************/


/***global node.datの読み込み***/
    fp = BBFE_sys_read_fopen(fp, INPUT_FILENAME_NODE, directory);

    fscanf(fp,"%d", &(pod_lb->total_num_nodes));

    global_x = BB_std_calloc_2d_double(global_x, pod_lb->total_num_nodes, 3);

	for(int i = 0; i < pod_lb->total_num_nodes; i++) {
		fscanf(fp,"%lf %lf %lf", 
                    &(global_x[i][0]), &(global_x[i][1]), &(global_x[i][2]));
	}

	fclose(fp);
/***********************/
    int* local_node_id;
    int* D_bc_id;

    for(int j = 0; j < nd1; j++) {
        //node_idの読み込み
        printf("node_idの読み込み\n");
        snprintf(fname_meta_node_in, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_NODE, j);
        fp = BBFE_sys_read_fopen(fp, fname_meta_node_in, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(local_num_nodes), &ndof);

        printf("local_num_nodes = %d", local_num_nodes);

        local_node_id = BB_std_calloc_1d_int(local_node_id, local_num_nodes);
        D_bc_id = BB_std_calloc_1d_int(D_bc_id, local_num_nodes);

        for(int i = 0; i < local_num_nodes; i++){
            fscanf(fp, "%d", &(local_node_id[i]));
        }
        fclose(fp);

        //D_bcの書き込み
        printf("D_bcの書き込み\n");
        snprintf(fname_D_bc_out, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_D_BC, j);    
        fp = BBFE_sys_write_fopen(fp, fname_D_bc_out, directory);

        int count = 0;
        for(int i = 0; i< local_num_nodes; i++) {
            int index = local_node_id[i];
            if(D_bc_exists[index] == true){
                D_bc_id[count] = i;
                count++;
            }
        }

        fprintf(fp,"%d %d\n",count, block_size);

        for(int i = 0; i < count; i++) {
            int index = D_bc_id[i];
            int global_index = local_node_id[index];
            fprintf(fp,"%d %d   %e\n",index, block_id, imposed_D_val[global_index]);
        }

        fclose(fp);

        //nodeの書き込み
        printf("nodeの書き込み\n");
        snprintf(fname_node_out, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_NODE, j);    
        fp = BBFE_sys_write_fopen(fp, fname_node_out, directory);

        fprintf(fp,"%d\n",local_num_nodes);

        for(int i = 0; i< local_num_nodes; i++) {
            int index = local_node_id[i];
            fprintf(fp,"%lf %lf %lf\n", 
                    global_x[index][0], global_x[index][1], global_x[index][2]);
        }

        fclose(fp);
    }

    BB_std_free_1d_bool(D_bc_exists, pod_lb->total_num_nodes);
	BB_std_free_1d_double(imposed_D_val, pod_lb->total_num_nodes);
	BB_std_free_2d_double(global_x, pod_lb->total_num_nodes, 3);
}

void lb_write_sendrecv(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    //隣接リストの取得

    printf("\n 隣接リストの取得\n");
    char fname_recv_in[BUFFER_SIZE];
    char fname_recv_out[BUFFER_SIZE];
    char fname_send_out[BUFFER_SIZE];

    FILE* fp;

    int** list_neib_graph;
    int* n_neib_graph;
    int* num_send_nodes;
    int* tmp_send;
    int* num_recv_nodes;
    int* tmp_recv;

    int* sum_n_internal;
    int ndof;

    list_neib_graph = BB_std_calloc_2d_int(list_neib_graph, pod_lb->max_n_neib, nd1);
    n_neib_graph = BB_std_calloc_1d_int(n_neib_graph, nd1);     //グラフとしての隣接領域番号

    for(int j = 0; j < nd1; j++){
        snprintf(fname_recv_in, BUFFER_SIZE, "meta_graph_parted.0/%s.recv.%d", OUTPUT_FILENAME_DOMAIN_GRAPH, j);
        fp = BBFE_sys_read_fopen(fp, fname_recv_in, directory);
        fscanf(fp, "%d %d", &(n_neib_graph[j]), &(ndof));

        printf("n_neib_graph = %d\n", n_neib_graph[j]);

        for(int i = 0; i < n_neib_graph[j]; i++){
            fscanf(fp, "%d", &(list_neib_graph[i][j]));
            printf("list_neib_graph = %d\n", list_neib_graph[i][j]);
        }
        fclose(fp);
    }

//recv,sendファイル
    printf("\n sendファイルの作成\n");
    for(int n = 0; n < nd1; n++){
        snprintf(fname_send_out, BUFFER_SIZE, "parted.0/%s.send.%d", INPUT_FILENAME_NODE, n);

        tmp_send = BB_std_calloc_1d_int(tmp_send, pod_lb->sum_n_internal[n]);
        num_send_nodes = BB_std_calloc_1d_int(num_send_nodes, nd1);

        fp = BBFE_sys_write_fopen(fp, fname_send_out, directory);

        //隣接する領域リストをraw_index,出力するノードをcolumn_index,tmpとして一時保管
        int column_index = 0;
        for(int j = 0; j < n_neib_graph[n]; j++){
            int count = 0;
            int row_index = list_neib_graph[j][n];
            for(int i = 0; i < pod_lb->total_num_nodes; i++){        //全通信を減らす
                if(pod_lb->bool_global_id_internal[i][n] == true 
                && pod_lb->bool_global_id_overlap[i][row_index] == true){
                    tmp_send[column_index] = i;
                    count++;
                    column_index++;
                }
            }
            num_send_nodes[j] = count;
            printf("count = %d\n", count);
            printf("index = %d\n", row_index);
        }


        fprintf(fp,"%d %d\n",n_neib_graph[n], column_index);

        for(int i = 0; i < n_neib_graph[n]; i++){
	        fprintf(fp,"%d\n", list_neib_graph[i][n]);
        }
        fprintf(fp,"%d\n", 0);

        int sum = 0;
        for(int i = 0; i < n_neib_graph[n]; i++){
            sum += num_send_nodes[i];
	        fprintf(fp,"%d\n", sum);
        }

        for(int i = 0; i < sum; i++){
            int index_conn = tmp_send[i];
            fprintf(fp,"%d\n",pod_lb->global_id[index_conn][n]+1);
        }

        fclose(fp);


        printf("\n recvファイルの作成\n");

        snprintf(fname_recv_out, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, n);

        tmp_recv = BB_std_calloc_1d_int(tmp_recv, pod_lb->sum_n_internal[n]);
        num_recv_nodes = BB_std_calloc_1d_int(num_recv_nodes, nd1);

        fp = BBFE_sys_write_fopen(fp, fname_recv_out, directory);

        column_index = 0;
        for(int j = 0; j < n_neib_graph[n]; j++){
            int count = 0;
            int row_index = list_neib_graph[j][n];
            for(int i = 0; i < pod_lb->total_num_nodes; i++){        //全通信を減らす
                if(pod_lb->bool_global_id_overlap[i][n] == true 
                && pod_lb->bool_global_id_internal[i][row_index] == true){
                    tmp_recv[column_index] = i;
                    count++;
                    column_index++;
                }
            }
            num_recv_nodes[j] = count;
            printf("count = %d\n", count);
            printf("index = %d\n", row_index);
        }

        fprintf(fp,"%d %d\n",n_neib_graph[n], column_index);

        for(int i =0; i < n_neib_graph[n]; i++){
	        fprintf(fp,"%d\n", list_neib_graph[i][n]);
        }
        fprintf(fp,"%d\n", 0);

        sum = 0;
        for(int i =0; i < n_neib_graph[n]; i++){
            sum += num_recv_nodes[i];
	        fprintf(fp,"%d\n", sum);
        }

        for(int i = 0; i < sum; i++){
            int index_conn = tmp_recv[i];
            fprintf(fp,"%d\n",pod_lb->global_id[index_conn][n]+1);
        }

        fclose(fp);

    }

	BB_std_free_2d_bool(pod_lb->bool_global_id_internal, pod_lb->total_num_nodes, nd1);
	BB_std_free_2d_bool(pod_lb->bool_global_id_overlap, pod_lb->total_num_nodes, nd1);
}


void lb_set_elem_global_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    //elemファイルの作成
    printf("\n elemファイルの作成\n");
    char fname_elem_n_internal[BUFFER_SIZE];
    char fname_n_internal_out[BUFFER_SIZE];
    char fname_id_in[BUFFER_SIZE];
    char fname_id_out[BUFFER_SIZE];
    char fname_id[BUFFER_SIZE];
    FILE* fp;

    char id[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];

	int recv_num_elems;
    int total_num_elems_local;
    int total_num_elems;
    int** elem_id_local;

    fp = BBFE_sys_read_fopen(fp, INPUT_FILENAME_ELEM, directory);
    fscanf(fp, "%d", &(pod_lb->total_num_elems));
    fclose(fp);

    int* sum_n_internal;
    int* sum_n_overlap;
    int* n_internal_vertex;

    sum_n_internal = BB_std_calloc_1d_int(sum_n_internal, nd1);
    sum_n_overlap = BB_std_calloc_1d_int(sum_n_overlap, nd1);
    n_internal_vertex = BB_std_calloc_1d_int(n_internal_vertex, nd1);

    bool** bool_global_id_internal;
    bool** bool_global_id_overlap;

    pod_lb->bool_global_id_internal = BB_std_calloc_2d_bool(bool_global_id_internal, pod_lb->total_num_elems, nd1);            
    pod_lb->bool_global_id_overlap = BB_std_calloc_2d_bool(bool_global_id_overlap, pod_lb->total_num_elems, nd1);

    int local_num_elems = 0;

    int count;
    int count_overlap;
    int doublecount_internal;
    int ndof;

    for(int n = 0; n < nd1; n++){
        int sum_internal = 0;
        int sum_overlap = 0;
        int doublecount = 0;
        for(int m = 0; m < pod_lb->Ddof[n]; m++){
            snprintf(fname_elem_n_internal, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_ELEM, pod_lb->list_n_neib[m][n]);
        
            fp = BBFE_sys_read_fopen(fp, fname_elem_n_internal, directory);
            fscanf(fp, "%s %d", char_n_internal, &(ndof));
            fscanf(fp, "%d", &(n_internal_vertex[n]));
            fclose(fp);

            printf("n_internal_vertex = %d\n", n_internal_vertex[n]);

            snprintf(fname_id_in, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, pod_lb->list_n_neib[m][n]);

            fp = BBFE_sys_read_fopen(fp, fname_id_in, directory);
            fscanf(fp, "%s", id);
            fscanf(fp, "%d %d", &(total_num_elems_local),&(ndof));

            printf("total_num_elems_local = %d\n", total_num_elems_local);

            elem_id_local = BB_std_calloc_2d_int(elem_id_local, total_num_elems_local, nd1);

            for(int i = 0; i < total_num_elems_local; i++) {
                fscanf(fp, "%d", &(elem_id_local[i][n]));
            }
            fclose(fp);

            for(int i = 0; i < n_internal_vertex[n]; i++){
                pod_lb->bool_global_id_internal[elem_id_local[i][n]][n] = true;
            }
            sum_internal += n_internal_vertex[n];

            count = 0;
            count_overlap = 0;
            doublecount_internal = 0;
            //overlap同士で重なった場合、内部要素、重ならない場合、overlap要素と場合分け
            for(int i = n_internal_vertex[n]; i < total_num_elems_local; i++){
                if(pod_lb->bool_global_id_overlap[elem_id_local[i][n]][n] == true){
                    if(pod_lb->bool_global_id_internal[elem_id_local[i][n]][n] == true){
                        doublecount_internal++;
                    }
                    pod_lb->bool_global_id_overlap[elem_id_local[i][n]][n] = false;
                    pod_lb->bool_global_id_internal[elem_id_local[i][n]][n] = true;
                    count++;
                }
                else{
                    pod_lb->bool_global_id_overlap[elem_id_local[i][n]][n] = true;
                    count_overlap++;
                }
            }

            sum_overlap += count_overlap - count;
            sum_internal += count - doublecount_internal;
            //total_num_elems_localファイル読み込みの値
            //n_internal_vertexファイル読み込み
            //count_overlap : overlap要素同士で重なった要素 すなわち内部要素
            //count : overlap要素同士で重ならなかった要素と重なった要素の総和
            //doublecount_internal : 内部要素同士で重なった要素の総和
        }
        
        for(int i = 0; i < pod_lb->total_num_elems; i++){
            if(pod_lb->bool_global_id_internal[i][n] == true && pod_lb->bool_global_id_overlap[i][n] == true){
                pod_lb->bool_global_id_overlap[i][n] = false;
                doublecount++;
            }
        }

        sum_n_internal[n] = sum_internal;
        sum_n_overlap[n] = sum_overlap - doublecount;

        printf("sum internal = %d\n", sum_n_internal[n]);
        printf("sum overlap = %d\n", sum_n_overlap[n]);

        //ファイル書き込み

        //elem.datファイル
        snprintf(fname_id_out, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, n);
        
        fp = BBFE_sys_write_fopen(fp, fname_id_out, directory);
        fprintf(fp,"#id\n");
        int local_num_elems = sum_n_internal[n] + sum_n_overlap[n];
        fprintf(fp,"%d %d\n", local_num_elems ,ndof);

        for(int i = 0; i < pod_lb->total_num_elems; i++){
            if(pod_lb->bool_global_id_internal[i][n] == true ){
		        fprintf(fp,"%d\n", i);
            }
        }

        for(int i = 0; i < pod_lb->total_num_elems; i++){
            if(pod_lb->bool_global_id_overlap[i][n] == true ){
		        fprintf(fp,"%d\n", i);
            }
        }
        fclose(fp);

        snprintf(fname_n_internal_out, BUFFER_SIZE, "parted.0/%s.n_internal.%d", INPUT_FILENAME_ELEM, n);
        
        fp = BBFE_sys_write_fopen(fp, fname_n_internal_out, directory);
        fprintf(fp,"#n_internal %d\n", ndof);
        fprintf(fp,"%d", sum_n_internal[n]);

        fclose(fp);

        BB_std_free_2d_int(elem_id_local, total_num_elems_local, nd1);
    }

	BB_std_free_2d_bool(pod_lb->bool_global_id_internal, pod_lb->total_num_elems, nd1);
	BB_std_free_2d_bool(pod_lb->bool_global_id_overlap, pod_lb->total_num_elems, nd1);
}

void lb_write_elem_local_id(
    const int       nd1,
    const char*     directory,
    POD_LB*        pod_lb)
{
    char fname_elem_out[BUFFER_SIZE];
    char fname_meta_elem_in[BUFFER_SIZE];
    char id[BUFFER_SIZE];

    FILE* fp;

    int** global_conn;
    int local_node_dof;
    int local_num_elems;
    int ndof;

/***global elem.datの読み込み***/
	fp = BBFE_sys_read_fopen(fp, INPUT_FILENAME_ELEM, directory);

	// read the number of elements
    fscanf(fp,"%d %d", &(pod_lb->total_num_elems), &(local_node_dof));
    global_conn = BB_std_calloc_2d_int(global_conn, pod_lb->total_num_elems, 8);

    printf("local_node_dof = %d\n", local_node_dof);

	// read the connectivities of elements
	for(int e = 0; e < pod_lb->total_num_elems ; e++) {
		for(int i = 0; i < local_node_dof ; i++) {
			fscanf(fp, "%d", &(global_conn[e][i]));
		}
	}

	fclose(fp);
/***********************/
    int* local_elem_id;

    for(int n = 0; n < nd1; n++) {
        //elem_idの読み込み
        snprintf(fname_meta_elem_in, BUFFER_SIZE, "parted.0/%s.id.%d", INPUT_FILENAME_ELEM, n);
        fp = BBFE_sys_read_fopen(fp, fname_meta_elem_in, directory);
        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &(local_num_elems), &ndof);

        printf("local_num_elems = %d", local_num_elems);

        local_elem_id = BB_std_calloc_1d_int(local_elem_id, local_num_elems);

        for(int i = 0; i < local_num_elems; i++){
            fscanf(fp, "%d", &(local_elem_id[i]));
        }
        fclose(fp);

        //elemの書き込み
        snprintf(fname_elem_out, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM, n);
        fp = BBFE_sys_write_fopen(fp, fname_elem_out, directory);

        fprintf(fp,"%d %d\n",local_num_elems, local_node_dof);

        for(int j = 0; j < local_num_elems; j++) {
            for(int i = 0; i < local_node_dof; i++){
                int index = local_elem_id[j];
                int index_conn = global_conn[index][i];
                fprintf(fp,"%d ",pod_lb->global_id[index_conn][n]);
            }
            fprintf(fp,"\n");
        }

        fclose(fp);
    }

	BB_std_free_2d_int(pod_lb->global_id, pod_lb->total_num_nodes, nd1);
    BB_std_free_2d_int(global_conn, pod_lb->total_num_elems, 8);
}