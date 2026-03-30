
#include "graph.h"

const char* get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename)
{
	const char* dir_name;

	if(argc < 2) { dir_name = "."; }
	else         { dir_name = argv[1]; }

	printf("%s Main directory: %s\n", codename, dir_name);

	return dir_name;
}


//ノードの重みを1ずつ増やす関数(テスト用)
void lb_set_node_weight(
        const int       num_subdomains,
	    const char*     directory)
{
	char fname[BUFFER_SIZE];

	snprintf(fname, BUFFER_SIZE, "%s", OUTPUT_FILENAME_NODAL_WEIGHT);

    FILE* fp;
    int ndof = 1;   //nodeの自由度

    fp = BBFE_sys_write_fopen(fp, fname, directory);
    fprintf(fp,"#node_weight\n");
    fprintf(fp,"%d %d\n",num_subdomains, ndof);

    for(int j = 0; j < num_subdomains; j++){
        for(int i = 0; i < ndof; i++){
            fprintf(fp,"%d", 1);
        }
        fprintf(fp,"\n");
    }

    fclose(fp);
}

void lb_set_node_weight_for_hyperreduction(
        const int       num_subdomains,
	    const char*     directory)
{
	char fname[BUFFER_SIZE];
    FILE* fp;
    int ndof = 1;   //nodeの自由度

    int* num_selected_elems;
    int* num_modes;

    num_selected_elems = BB_std_calloc_1d_int(num_selected_elems, num_subdomains);
    num_modes = BB_std_calloc_1d_int(num_modes, num_subdomains);

    for(int i = 0; i < num_subdomains; i++){
    	snprintf(fname, BUFFER_SIZE,"pod_modes/num_modes.dat.%d", i);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d", &(num_modes[i]));
        printf("%d", num_modes[i]);
        fclose(fp);
    }

    for(int i = 0; i < num_subdomains; i++){
    	snprintf(fname, BUFFER_SIZE,"DDECM/num_selected_node.%d.txt", i);
        fp = BBFE_sys_read_fopen(fp, fname, directory);
        fscanf(fp, "%d", &(num_selected_elems[i]));
        printf("%d", num_selected_elems[i]);
        fclose(fp);
    }

   	snprintf(fname, BUFFER_SIZE, "merged_graph/%s", OUTPUT_FILENAME_NODAL_WEIGHT);
    fp = BBFE_sys_write_fopen(fp, fname, directory);
    fprintf(fp,"#node_weight\n");
    fprintf(fp,"%d %d\n",num_subdomains, ndof);

    for(int j = 0; j < num_subdomains; j++){
        for(int i = 0; i < ndof; i++){
           fprintf(fp,"%d", num_selected_elems[j]);
	}
        fprintf(fp,"\n");
    }

    fclose(fp);

    BB_std_free_1d_int(num_selected_elems, num_subdomains);
    BB_std_free_1d_int(num_modes, num_subdomains);
}


const char* hlpod_get_directory_name(
    const char* directory)
{
    static const char* dir_name;
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE, "%s/load_balancing", directory);
    dir_name = fname;
    printf("fname = %s", dir_name);

    return dir_name;
}

//グラフファイルの作成
void hlpod_set_subdomain_graph(
        const int       num_subdomains,
	    const char*     directory)
{
    char fname_n_internal[BUFFER_SIZE];
	char fname_out[BUFFER_SIZE];

	FILE* fp_in;
    FILE* fp_out;

	int* num_n_neib;
	int recv_num_nodes;
	int n_neib;

	snprintf(fname_out, BUFFER_SIZE, "%s", OUTPUT_FILENAME_DOMAIN_GRAPH);

    fp_out = BBFE_sys_write_fopen(fp_out, fname_out, directory);
    fprintf(fp_out,"%d\n",num_subdomains);
	fclose(fp_out);

    for(int m = 0; m < num_subdomains; m++){
        snprintf(fname_n_internal, BUFFER_SIZE, "parted.0/%s.recv.%d", INPUT_FILENAME_NODE, m);
    
        fp_in = BBFE_sys_read_fopen(fp_in, fname_n_internal, directory);
        fscanf(fp_in, "%d %d", &(n_neib), &(recv_num_nodes));

        num_n_neib = BB_std_calloc_1d_int(num_n_neib, n_neib);

        for(int i = 0; i < n_neib; i++){
            fscanf(fp_in, "%d", &(num_n_neib[i]));
        }
        fclose(fp_in);

        fp_out = ROM_BB_write_add_fopen(fp_out, fname_out, directory);
        
        fprintf(fp_out,"%d %d", m, n_neib);
	    for(int i = 0; i < n_neib; i++){
		    fprintf(fp_out," %d",num_n_neib[i]);
	    }
        fprintf(fp_out,"\n");
        fclose(fp_out);

        BB_std_free_1d_int(num_n_neib, n_neib);
    }
}

void open_gedatsu_external_file(
    const int       nd1,
    const char*  	directory)
{
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE, 
        "cd %s/load_balancing && cp -r ./../parted.0/metagraph.dat ./ && gedatsu_nodal_graph_partitioner -n %d -i %s -d %s -inw %s",
        directory, nd1, OUTPUT_FILENAME_DOMAIN_GRAPH, OUTPUT_FILENAME_PARTED_METAGRAPH, OUTPUT_FILENAME_NODAL_WEIGHT);

    system(fname);
}


//分割されたメタグラフを対応するMPIプロセスで読み込む
void lb_read_meta_graph_mpi(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb)
{
    char fname_id_graph[BUFFER_SIZE];
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    char id[BUFFER_SIZE];
    
    FILE* fp;

    int total_n_neib;
    int graph_ndof;     //graph_weightの節点自由度
    int* domain_list;  //メタメッシュの隣接領域番号リスト
    
/*max_n_neibの共有*/
    //内部領域リストの取得
    lpod_lb->max_Ddof = 0;
    for(int j = 0; j < nd1; j++){
        snprintf(fname_n_internal_graph, BUFFER_SIZE, "metagraph_parted.0/%s.n_internal.%d", OUTPUT_FILENAME_DOMAIN_GRAPH, j);
        
        fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
        fscanf(fp, "%s %d", char_n_internal, &(graph_ndof));
        fscanf(fp, "%d", &(lpod_lb->Ddof));
        fclose(fp);

        if(lpod_lb->Ddof > lpod_lb->max_Ddof){
            lpod_lb->max_Ddof = lpod_lb->Ddof;
        }
    }

    const int myrank = monolis_mpi_get_global_my_rank();
    //分割メタグラフの内部節点数, 領域自由度の取得
    snprintf(fname_n_internal_graph, BUFFER_SIZE, 
        "metagraph_parted.0/%s.n_internal.%d", 
        OUTPUT_FILENAME_DOMAIN_GRAPH, myrank);
        
    fp = BBFE_sys_read_fopen(fp, fname_n_internal_graph, directory);
    fscanf(fp, "%s %d", char_n_internal, &(graph_ndof));
    fscanf(fp, "%d", &(lpod_lb->Ddof));
    fclose(fp);


/*分割メタグラフidから連結関係のリストの作成*/
    //Ddofのループでリストを取得
    lpod_lb->domain_list = BB_std_calloc_1d_int(domain_list, lpod_lb->max_Ddof);
    snprintf(fname_id_graph, BUFFER_SIZE, 
        "metagraph_parted.0/%s.id.%d", 
        OUTPUT_FILENAME_DOMAIN_GRAPH, myrank);
    fp = BBFE_sys_read_fopen(fp, fname_id_graph, directory);
    fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_n_neib), &(graph_ndof));

    for(int i = 0; i < lpod_lb->Ddof; i++){
        fscanf(fp, "%d", &(lpod_lb->domain_list[i]));
    }
    fclose(fp);
}
