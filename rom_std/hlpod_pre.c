
#include "hlpod_pre.h"
    
static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_NODE        = "node.dat";

void ROM_std_hlpod_get_meta_neib(
    MONOLIS_COM*  	monolis_com,
    HLPOD_META*		hlpod_meta,
    const char*     metagraph,
    const char*     directory)
{
    int num_meta_nodes;
    int tmp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    FILE* fp = NULL;
    
    const int n_internal_vertex = monolis_com->n_internal_vertex;
    const int myrank = monolis_mpi_get_global_my_rank();

    int num_neib;

    /* neib_id 配列を確保して隣接IDを読む */
    snprintf(fname, BUFFER_SIZE, "%s.recv.%d", metagraph, myrank);
    fp = ROM_BB_read_fopen(fp, fname, directory);

    fscanf(fp, "%d %d", &(num_neib), &tmp);

    int* neib_id = BB_std_calloc_1d_int(neib_id, num_neib);
    for (int i = 0; i < num_neib; i++) {
        fscanf(fp, "%d", &(neib_id[i]));
    }
    fclose(fp);

    /* 各隣接の内部ノード数を読み込み */
    int* n_internal = BB_std_calloc_1d_int(n_internal, num_neib);
    hlpod_meta->n_internal_sum = 0;
    for (int m = 0; m < num_neib; m++) {
        snprintf(fname, BUFFER_SIZE, "%s.n_internal.%d", metagraph, neib_id[m]);
        n_internal[m] = ROM_std_hlpod_read_n_internal(fp, fname, directory);
        hlpod_meta->n_internal_sum += n_internal[m];
    }

    /* 全隣接分の内部ノード数 + 自分の内部頂点数 */
    hlpod_meta->subdomain_id_neib = BB_std_calloc_1d_int(hlpod_meta->subdomain_id_neib, hlpod_meta->n_internal_sum + n_internal_vertex);

    /* 自身の .id ファイルからノードIDを読み込み */
    snprintf(fname, BUFFER_SIZE, "%s.id.%d", metagraph, myrank);
    ROM_std_hlpod_read_node_id(fp, hlpod_meta->subdomain_id_neib, n_internal_vertex, fname, directory);

    int index_internal = 0;
    index_internal += n_internal_vertex;

    /* 隣接ノードIDを読み込み */
    for (int m = 0; m < num_neib; m++) {
        snprintf(fname, BUFFER_SIZE, "%s.id.%d", metagraph, neib_id[m]);
        fp = ROM_BB_read_fopen(fp, fname, directory);

        fscanf(fp, "%s", id);
        fscanf(fp, "%d %d", &num_meta_nodes, &tmp);

        for (int i = 0; i < n_internal[m]; i++) {
            fscanf(fp, "%d", &(hlpod_meta->subdomain_id_neib[index_internal]));
            index_internal++;
        }
        fclose(fp);
    }

}


void ROM_std_hlpod_get_subdomain_id(
    HLPOD_VALUES* 	hlpod_vals,
    HLPOD_META*		hlpod_meta,
    const int       num,
    const char*     metagraph,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    int num_meta_nodes;
    int tmp;
    const int myrank = monolis_mpi_get_global_my_rank();

    snprintf(fname, BUFFER_SIZE, "%s.n_internal.%d", metagraph, num);
    hlpod_vals->num_2nd_subdomains = ROM_std_hlpod_read_n_internal(fp, fname, directory);

    snprintf(fname, BUFFER_SIZE, "%s.id.%d", metagraph, myrank);
    fp = ROM_BB_read_fopen(fp, fname, directory);

    fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &num_meta_nodes, &tmp);

    hlpod_meta->subdomain_id =
        BB_std_calloc_1d_int(hlpod_meta->subdomain_id, num_meta_nodes);

    for (int i = 0; i < num_meta_nodes; i++) {
        fscanf(fp, "%d", &(hlpod_meta->subdomain_id[i]));
    }
    fclose(fp);

}


void ROM_std_hlpod_read_node_id_pod_subd(
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int		num_2nd_subdomains,
    const char*     label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    int* node_id_local;
    int n_internal_vertex;

    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd, num_2nd_subdomains);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, total_num_nodes);

    int index = 0;
    int index_row = 0;

    for(int m = 0; m < num_2nd_subdomains; m++){
        snprintf(fname, BUFFER_SIZE, "%s/%s.n_internal.%d", label, INPUT_FILENAME_NODE, m);
        n_internal_vertex = ROM_std_hlpod_read_n_internal(fp, fname, directory);

        node_id_local = BB_std_calloc_1d_int(node_id_local, n_internal_vertex);

        snprintf(fname, BUFFER_SIZE, "%s/%s.id.%d", label, INPUT_FILENAME_NODE, m);
        ROM_std_hlpod_read_node_id(fp, node_id_local, n_internal_vertex, fname, directory);
        int sum = 0;
        for(int i = 0; i < n_internal_vertex; i++){
            hlpod_mat->node_id[index] = node_id_local[i];
            index++;
            sum++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = sum;

        BB_std_free_1d_int(node_id_local, n_internal_vertex);
    }

}


void ROM_std_hlpod_read_node_id_pod_subd_para(
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		total_num_nodes,
    const int 		N_internal_vertex,
    const int		num_2nd_subdomains,
    const char*     label,
    const char*     directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    int n_internal_vertex;

    hlpod_mat->n_internal_vertex_subd = BB_std_calloc_1d_int(hlpod_mat->n_internal_vertex_subd , num_2nd_subdomains);
    hlpod_mat->node_id = BB_std_calloc_1d_int(hlpod_mat->node_id, N_internal_vertex);

    int index = 0;
    for(int m = 0; m < num_2nd_subdomains; m++){
        snprintf(fname, BUFFER_SIZE, "%s/%s.n_internal.%d", label, INPUT_FILENAME_NODE, hlpod_meta->subdomain_id[m]);
        n_internal_vertex = ROM_std_hlpod_read_n_internal(fp, fname, directory);

        for(int i = 0; i < n_internal_vertex; i++){
            hlpod_mat->node_id[index] = index;
            index++;
        }
        hlpod_mat->n_internal_vertex_subd[m] = n_internal_vertex;

    }

}

/*solve as dense matrix*/
void ROM_std_hlpod_set_nonzero_pattern(
    MONOLIS*     	monolis,
    HLPOD_MAT*      hlpod_mat,
    const int 		num_modes)
{
    const int k = num_modes;
    int* index;
    int* item;
    int* connectivity;

    index = BB_std_calloc_1d_int(index, 1);
    item = BB_std_calloc_1d_int(item, 1);
    
    index[0] = 0;
    item[0] = 1;

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        monolis,
        1,					//nnode:節点数
        k,					//ndof:節点あたりの自由度
        index,				//節点グラフを定義するindex配列
        item);				//設定グラフを定義するitem配列

    connectivity = BB_std_calloc_1d_int(connectivity, 1);
    connectivity[0] = 0;

    BB_std_free_1d_int(index, 1);
    BB_std_free_1d_int(item, 1);
    BB_std_free_1d_int(connectivity, 1);
}

void ROM_std_hlpod_set_nonzero_pattern_bcsr(
    MONOLIS*     	monolis,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int		num_modes,
    const char*     label,
    const char*		directory)
{
    const char* fname;
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);

    fp = ROM_BB_read_fopen(fp, fname, directory);

    fscanf(fp, "%d", &(num_nodes));

    for(int e = 0; e < num_nodes; e++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes) );

        for(int i = 0; i < num_adj_nodes; i++) {
            fscanf(fp, "%d", &(tmp));
        }
        sum += num_adj_nodes;
    }
    fclose(fp);
    hlpod_meta->index = BB_std_calloc_1d_int(hlpod_meta->index, num_nodes + 1);
    hlpod_meta->item = BB_std_calloc_1d_int(hlpod_meta->item, sum);

    sum = 0;
    int index_sum = 0;
    hlpod_meta->index[0] = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);
    fp = ROM_BB_read_fopen(fp, fname, directory);
    fscanf(fp, "%d", &(num_nodes));

    for(int i = 0; i < num_nodes; i++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes));

        for(int j = 0; j < num_adj_nodes; j++) {
            fscanf(fp, "%d", &(tmp));
            hlpod_meta->item[index_sum] = tmp;
            index_sum++;
        }
        sum += num_adj_nodes;
        hlpod_meta->index[i + 1] = sum;
    }	
    fclose(fp);

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
        monolis,
        num_nodes,
        hlpod_mat->num_modes_internal,
        hlpod_meta->index,
        hlpod_meta->item);

}


void ROM_std_hlpod_set_nonzero_pattern_bcsr_para(
    MONOLIS*     	monolis,
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const char*     metagraph,
    const char*		directory)
{
    char fname[BUFFER_SIZE];
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;

    const int myrank = monolis_mpi_get_global_my_rank();

    snprintf(fname, BUFFER_SIZE, "%s.%d", metagraph, myrank);
    fp = ROM_BB_read_fopen(fp, fname, directory);

    fscanf(fp, "%d", &(num_nodes));

    for(int e = 0; e < num_nodes; e++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes) );

        for(int i = 0; i < num_adj_nodes; i++) {
            fscanf(fp, "%d", &(tmp));
        }
        sum += num_adj_nodes;
    }

    hlpod_meta->index = BB_std_calloc_1d_int(hlpod_meta->index, num_nodes + 1);
    hlpod_meta->item = BB_std_calloc_1d_int(hlpod_meta->item, sum);

    sum = 0;
    int index = 0;
    hlpod_meta->index[0] = 0;

    snprintf(fname, BUFFER_SIZE, "%s.%d", metagraph, myrank);
    fp = ROM_BB_read_fopen(fp, fname, directory);
    fscanf(fp, "%d", &(num_nodes));

    for(int i = 0; i < num_nodes; i++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes));

        for(int j = 0; j < num_adj_nodes; j++) {
            fscanf(fp, "%d", &(hlpod_meta->item[index]));
            index++;
        }
        sum += num_adj_nodes;
        hlpod_meta->index[i + 1] = sum;
    }	
    fclose(fp);

    hlpod_meta->num_meta_nodes = num_nodes;
}

//for arbit_dof_monolis_solver
void ROM_std_hlpod_get_n_dof_list(
    MONOLIS_COM*  	monolis_com,
    HLPOD_MAT* 	    hlpod_mat,
    HLPOD_META*		hlpod_meta,
    const int 		max_num_bases)
{
    const int M = max_num_bases;
    const int M_all = M * hlpod_meta->num_meta_nodes;

    int length1;
    int length2;

    int* node_weight;

    hlpod_meta->n_dof_list = BB_std_calloc_1d_int(hlpod_meta->n_dof_list, hlpod_meta->num_meta_nodes);

    node_weight = BB_std_calloc_1d_int(node_weight, monolis_com->n_internal_vertex);

    int index1 = 0;
    int index2 = 0;

    for(int k = 0; k < monolis_com->n_internal_vertex; k++){
        node_weight[k] = 0;

        int iS = hlpod_meta->index[k];
        int iE = hlpod_meta->index[k + 1];

        for(int i = iS; i < iE; i++){
            for(int j = 0; j < hlpod_meta->n_internal_sum + monolis_com->n_internal_vertex; j++){
                if (hlpod_meta->subdomain_id[hlpod_meta->item[i]] == hlpod_meta->subdomain_id_neib[j]){
                    int IS = hlpod_mat->num_modes_1stdd[k];
                    int IE = hlpod_mat->num_modes_1stdd[k+1];

                    int IIS = hlpod_mat->num_modes_1stdd[j];
                    int IIE = hlpod_mat->num_modes_1stdd[j + 1];

                    length1 = IE - IS;
                    index1 = 0;

                    for(int m = IS; m < IE; m++){					
                        length2 = IIE - IIS;
                        index2 = 0;
                        for(int n = IIS; n < IIE; n++){		
                            index2++;
                        }
                        index1++;

                        hlpod_meta->n_dof_list[hlpod_meta->item[i]] = IIE - IIS;
                    }
                    hlpod_meta->n_dof_list[k] = IE - IS;
                    node_weight[k] += length1*length2;
                }
            }
        }
        node_weight[k] += length1*length1;
    }
}
