
#include "core.h"
#include "core_mpi.h"
#include "lb_dataset.h"
#include "graph.h"

static const char* OPTION_NUM_1STDD = "-np1";
static const char* OPTION_NUM_2NDDD = "-np2";
static const char* OPTION_DOF_BC    = "-nbc";
static const char* OPTION_FILENAME_BC = "-fbc";


typedef struct
{
	const char* directory;
    int num_1st_subdomains;
    int num_2nd_subdomains;
    char* bc_label;
    int bc_dof;

} CONDITIONS;


typedef struct
{
	CONDITIONS   cond;
    POD_LB			pod_lb;
	LPOD_LB			lpod_lb;

} FE_SYSTEM;


void set_condition(
    int argc,
    char* argv[],
    CONDITIONS* cond)
{
	int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_1STDD);
    if(num == -1) {
		printf("\nargs error num_1st_subdomains");
		exit(1);
    }
    else {
        cond->num_1st_subdomains = atoi(argv[num + 1]);	//hddにおける1層目のメッシュ分割数 (POD計算領域数)
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_2NDDD);
    if(num == -1) {
		printf("\nargs error num_2nd_subdomains");
		exit(1);
    }
    else {
		cond->num_2nd_subdomains = atoi(argv[num + 1]);	//hddにおける2層目のメッシュ分割数 (並列計算領域数)
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_DOF_BC);
    if(num == -1) {
		printf("\nargs error num_2nd_subdomains");
		exit(1);
    }
    else {
		cond->bc_dof = atoi(argv[num + 1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_FILENAME_BC);
    if(num == -1) {
		printf("\nargs error num_2nd_subdomains");
		exit(1);
    }
    else {
		cond->bc_label = argv[num + 1];
    }

	//printf("\nnum_1st_dd = %d\n", cond->num_1st_subdomains);
	//printf("\nnum_2nd_dd = %d\n", cond->num_2nd_subdomains);
}

void set_graph_to_communicator(
        POD_LB*        pod_lb,
        const int       num_2nd_subdomains,
	    const char*     directory)
{
   //hddにおける1層目のメッシュ分割数 (POD計算領域数)
    const int  nd1 = num_2nd_subdomains;
    //分割されたメタグラフからの読み取り
    lb_read_meta_graph(
        nd1,
        directory,
        pod_lb);

    lb_set_node_global_id(
        nd1,
        directory,
        pod_lb);

    lb_write_node_local_id(
        nd1,
        directory,
        pod_lb);

    lb_write_sendrecv(
        nd1,
        directory,
        pod_lb);   

    lb_set_elem_global_id(
        nd1,
        directory,
        pod_lb);

    lb_write_elem_local_id(
        nd1,
        directory,
        pod_lb);
}


void set_graph_to_communicator_mpi(
        LPOD_LB*        lpod_lb,
        const int       num_2nd_subdomains,
        const int       dof,
        const char*     label,
	    const char*     directory)
{
    const int  nd1 = num_2nd_subdomains;

    lb_read_meta_graph_mpi(
        nd1,
        directory,
        lpod_lb);

    lb_set_node_internal(
        nd1,
        directory,
        lpod_lb);

    lb_set_neib_list(
        nd1,
        directory,
        lpod_lb);
    
    lb_set_node_overlap_without_write(
        nd1,
        directory,
        lpod_lb);

    if(dof==3){
        lb_set_D_bc_v(
            nd1,
            directory,
            3,
            label,
            lpod_lb);    
    }
    else if(dof==1){
        lb_set_D_bc_p(
            nd1,
            directory,
            label,
            lpod_lb);
    }
    else{
        printf("Error : dof is not set\n");
    }
}

int main (
		int argc,
		char* argv[])
{
	printf("\n");

	BB_vtk_void();

	FE_SYSTEM sys;

	monolis_global_initialize();

	sys.cond.directory = get_directory_name(argc, argv, CODENAME);

    set_condition(argc, argv, &(sys.cond));
    
    set_graph_to_communicator_mpi(&(sys.lpod_lb), sys.cond.num_2nd_subdomains, sys.cond.bc_dof, sys.cond.bc_label, sys.cond.directory);

	monolis_global_finalize();

	printf("\n");

	return 0;
}
