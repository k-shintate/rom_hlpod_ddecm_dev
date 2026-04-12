
#include "./../mag_core/convdiff_core.h"
#include "./../mag_core/nedelec_core.h"

#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include <complex.h>

#include "mag_dataset.h"
#include "core_ROM.h"

#include "3ph_tr_NR_JA.h"
#include "team21c.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 3;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-6;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 100000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.001;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 10;

const double DELTA    = 1.0E-06;
const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

static const char* OUTPUT_FILENAME_NEDELEC_VTK          = "result_nedelec_%06d.vtk";

const double mu0 = 4.0*M_PI*1e-7; // H/m
const double Nu  = 1.0 / mu0;     // ← ここを μ0 に合わせる
const double Sigma = 5.0*1.0e7;                  // 無次元化（典型導電率）
const double Sigma_cop =  5.0*1.0e8;                  // 無次元化（典型導電率）
const double Sigma_air = 0.00000;                  // 無次元化（典型導電率）

const double freq  = 10.0;                 // [Hz]
//const double Omega = 2.0 * M_PI * freq;   // 角周波数
const double Omega = 0;   // 角周波数

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";



void ROM_read_args(
    int 		argc,
    char* 		argv[],
    ROM_PRM*    rom_prm)
{
	int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_1STDD);
    if(num == -1) {
		printf("\nargs error num_subdomains");
		exit(1);
    }
    else {
        rom_prm->num_subdomains = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_MODES);
    if(num == -1) {
		printf("\nargs error num_modes");
		exit(1);
    }
    else {
		rom_prm->num_modes = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_PADAPTIVE);
    if(num == -1) {
		printf("\nargs error rom_epsilon");
		exit(1);
    }
    else {
        rom_prm->rom_epsilon = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_SOLVER_TYPE);
    if(num == -1) {
		printf("\nargs error solver_type");
		exit(1);
    }
    else {
        rom_prm->solver_type = atof(argv[num+1]);
    }

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);

}



static int cmp_edgecand(const void *pa, const void *pb) {
    const EdgeCand *A = (const EdgeCand*)pa;
    const EdgeCand *B = (const EdgeCand*)pb;
    if (A->a != B->a) return (A->a < B->a) ? -1 : 1;
    if (A->b != B->b) return (A->b < B->b) ? -1 : 1;
    // 同じキー内では順序はどうでもいい
    if (A->elem != B->elem) return (A->elem < B->elem) ? -1 : 1;
    return (A->ledge < B->ledge) ? -1 : (A->ledge > B->ledge);
}

void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->phi        = BB_std_calloc_1d_double(vals->phi,     total_num_nodes);
	vals->error    = BB_std_calloc_1d_double(vals->error, total_num_nodes);
	vals->theo_sol = BB_std_calloc_1d_double(vals->theo_sol, total_num_nodes);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}

void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->phi, fe->total_num_nodes, "phi");

	fclose(fp);

}

void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			filename,
			sys->cond.directory,
			t);

}

void ref_hex8_shape_grads(const double xi[3],
                                        double dN_dxi[8],
                                        double dN_deta[8],
                                        double dN_dzeta[8])
{
    const double r = xi[0], s = xi[1], t = xi[2];
    const double c = 1.0/8.0;

    // 各ノードの符号（参照節点座標の ±1）
    const int sx[8] = {-1,+1,+1,-1,-1,+1,+1,-1};
    const int sy[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
    const int sz[8] = {-1,-1,-1,-1,+1,+1,+1,+1};

    for (int a=0; a<8; ++a) {
        dN_dxi[a]   = c * sx[a] * (1 + sy[a]*s) * (1 + sz[a]*t);
        dN_deta[a]  = c * sy[a] * (1 + sx[a]*r) * (1 + sz[a]*t);
        dN_dzeta[a] = c * sz[a] * (1 + sx[a]*r) * (1 + sy[a]*s);
    }
}

static inline int invert3x3(const double A[3][3], double Ainv[3][3], double *det_out)
{
    const double a=A[0][0], b=A[0][1], c=A[0][2];
    const double d=A[1][0], e=A[1][1], f=A[1][2];
    const double g=A[2][0], h=A[2][1], i=A[2][2];

    const double co00 =  e*i - f*h;
    const double co01 =  c*h - b*i;
    const double co02 =  b*f - c*e;
    const double co10 =  f*g - d*i;
    const double co11 =  a*i - c*g;
    const double co12 =  c*d - a*f;
    const double co20 =  d*h - e*g;
    const double co21 =  b*g - a*h;
    const double co22 =  a*e - b*d;

    const double det = a*co00 + b*co10 + c*co20;
    if (det_out) *det_out = det;
    if (fabs(det) < 1e-30) return 0;

    const double invdet = 1.0/det;
    
    Ainv[0][0] = co00*invdet;  Ainv[0][1] = co10*invdet;  Ainv[0][2] = co20*invdet;
    Ainv[1][0] = co01*invdet;  Ainv[1][1] = co11*invdet;  Ainv[1][2] = co21*invdet;
    Ainv[2][0] = co02*invdet;  Ainv[2][1] = co12*invdet;  Ainv[2][2] = co22*invdet;
    return 1;
}

static inline int compute_jacobian_and_inverse_at(
    const double xnod[8][3],
    const double xi[3],
    double J[3][3],
    double J_inv[3][3],
    double *detJ)
{
    double dN_dxi[8], dN_deta[8], dN_dzeta[8];
    ref_hex8_shape_grads(xi, dN_dxi, dN_deta, dN_dzeta);

    // J(m,n) = Σ_a x_a(m) * ∂N_a/∂(xi,eta,zeta)_n
    // m = 0:x,1:y,2:z   n = 0:xi,1:eta,2:zeta
    for (int m=0; m<3; ++m) J[m][0]=J[m][1]=J[m][2]=0.0;

    for (int a=0; a<8; ++a) {
        const double xa = xnod[a][0];
        const double ya = xnod[a][1];
        const double za = xnod[a][2];

        J[0][0] += xa * dN_dxi[a];   J[0][1] += xa * dN_deta[a];   J[0][2] += xa * dN_dzeta[a];
        J[1][0] += ya * dN_dxi[a];   J[1][1] += ya * dN_deta[a];   J[1][2] += ya * dN_dzeta[a];
        J[2][0] += za * dN_dxi[a];   J[2][1] += za * dN_deta[a];   J[2][2] += za * dN_dzeta[a];
    }

    return invert3x3(J, J_inv, detJ);
}

void edge_midpoint_ref_coords(int p, double xi[3]) {
    static const double P[12][3] = {
        // x-edges 0..3 : (eta,zeta) = (--),(+-),(-+),(++)
        { 0,-1,-1}, { 0,+1,-1}, { 0,-1,+1}, { 0,+1,+1},

        // y-edges 4..7 : (xi,zeta) = (--),(+-),(-+),(++)
        {-1, 0,-1}, {+1, 0,-1}, {-1, 0,+1}, {+1, 0,+1},

        // z-edges 8..11: (xi,eta)  = (--),(+-),(-+),(++)
        {-1,-1, 0}, {+1,-1, 0}, {-1,+1, 0}, {+1,+1, 0}
    };
    xi[0]=P[p][0]; xi[1]=P[p][1]; xi[2]=P[p][2];
}

// 各要素 e について：12 点の N[j] を保存（あなたの ned->N_edge[e][p][j][3]）
void precompute_N_edge_for_element(
    int e,
    const double xnod[8][3],            // 要素 e の 8 節点座標
    double N_edge_tab[12][12][3])       // [p][j][comp] = ned->N_edge[e][p][j][comp]
{
    for (int p = 0; p < 12; ++p) {
        double xi[3]; edge_midpoint_ref_coords(p, xi);
        double J[3][3], J_inv[3][3], detJ;
        compute_jacobian_and_inverse_at(xnod, xi, J, J_inv, &detJ); // ←既存の J 計算
        double *N_edge[12]; double Nbuf[12][3];
        for (int j=0;j<12;++j) N_edge[j]=Nbuf[j];

        BBFE_std_shapefunc_hex1st_nedelec_get_val(xi, N_edge, J_inv);

        for (int j=0;j<12;++j){
            for (int d=0; d<3; ++d)
                N_edge_tab[p][j][d] = N_edge[j][d];
        }
    }

}


void reconstruct_field_at_edge_points(
    BBFE_DATA* fe,
    NEDELEC* ned,
    const double* u_coeff,      // [num_global_edges]
    double** field_at_edges,    // [num_global_edges][3]
    int total_global_elems)
{
    // 3方向の和で総エッジ数を取る（構造Hexの例）
    
    int num_global_edges = fe->total_num_nodes;

    double N_edge_tab[12][12][3];
    double xnod[8][3];

    for (int ge=0; ge<num_global_edges; ++ge)
        field_at_edges[ge][0]=field_at_edges[ge][1]=field_at_edges[ge][2]=0.0;

    int *hits = (int*)calloc(num_global_edges,sizeof(int));

    for (int e=0; e<total_global_elems; ++e) {
        
        for(int i = 0; i < 8; i++){
            xnod[i][0] = fe->x[fe->conn[e][i]][0];
            xnod[i][1] = fe->x[fe->conn[e][i]][1];
            xnod[i][2] = fe->x[fe->conn[e][i]][2];
        }

        precompute_N_edge_for_element(
            e,
            xnod,            // 要素 e の 8 節点座標
            N_edge_tab);       // [p][j][comp] = ned->N_edge[e][p][j][comp]

        double c[12];
        for (int j=0; j<12; ++j) {
            int ge_j = ned->nedelec_conn[e][j];
            c[j] = u_coeff[ge_j] * ned->edge_sign[e][j];
        }

        for (int k=0; k<12; ++k) {
            int ge_k = ned->nedelec_conn[e][k];
            int p = k;
            double v0=0, v1=0, v2=0;
            for (int j=0; j<12; ++j) {
                v0 += c[j] * N_edge_tab[p][j][0];
                v1 += c[j] * N_edge_tab[p][j][1];
                v2 += c[j] * N_edge_tab[p][j][2];
            }
            field_at_edges[ge_k][0] += v0;
            field_at_edges[ge_k][1] += v1;
            field_at_edges[ge_k][2] += v2;
            ++hits[ge_k];
        }
    }

    for (int ge=0; ge<num_global_edges; ++ge) if (hits[ge]) {
        field_at_edges[ge][0] /= hits[ge];
        field_at_edges[ge][1] /= hits[ge];
        field_at_edges[ge][2] /= hits[ge];
    }
    free(hits);
}


void output_result_file_vtk_nedelec_edge(
		BBFE_DATA*      fe,
		VALUES*         vals,
        NEDELEC*        ned,
		const char*     filename,
		const char*     directory,
        int total_edges,
		double          t)
{

    double** node_result = BB_std_calloc_2d_double(node_result, fe->total_num_nodes, 3);
    double** edge_result = BB_std_calloc_2d_double(edge_result, fe->total_num_nodes, 3);

    reconstruct_field_at_edge_points(fe, ned, vals->V, edge_result, fe->total_num_elems);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    BB_vtk_write_points_3d(fp, total_edges, ned->nedelec_coords);
	BB_vtk_write_cells(fp, fe->total_num_elems, 12, ned->nedelec_conn);
    BB_vtk_write_cell_types(fp, fe->total_num_elems, 4);

    fprintf(fp, "POINT_DATA %d\n", total_edges);
    BB_vtk_write_point_vals_scalar(fp, vals->V, total_edges, "coordinates");
    BB_vtk_write_point_vals_vector(fp, edge_result, total_edges, "V_edge");

	fclose(fp);
}

void output_B_cell_vtk(
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    NEDELEC*    ned,
    const double* Aphi,
    const char*  filename,
    const char*  directory)
{
    // 要素ごとの B を作る
    double** B_cell = BB_std_calloc_2d_double(B_cell, fe->total_num_elems, 3);
    compute_B_cell_average(fe, basis, ned, Aphi, B_cell);

    FILE* fp = BBFE_sys_write_fopen(fp, filename, directory);

    // 形状・接続は既存関数に任せる
    switch(fe->local_num_nodes){
        case 4: BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);       break;
        case 8: BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);  break;
        default: fprintf(stderr,"unsupported local_num_nodes\n"); break;
    }

    // 既存の点データを書きたければここで POINT_DATA ... を挿む
    // （今回は B を CELL_DATA として書く）
    fprintf(fp, "CELL_DATA %d\n", fe->total_num_elems);
    fprintf(fp, "VECTORS B double\n");
    for(int e=0; e<fe->total_num_elems; ++e){
        fprintf(fp, "%.16e %.16e %.16e\n",
            B_cell[e][0], B_cell[e][1], B_cell[e][2]);
    }

    fclose(fp);
    BB_std_free_2d_double(B_cell, fe->total_num_elems, 3);
}

void output_files_nedelec(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
    char fname_B[BUFFER_SIZE];

	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_NEDELEC_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);
    snprintf(fname_B  , BUFFER_SIZE, "B_cell_%06d.vtk",          file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
	output_result_file_vtk_nedelec_edge(
			&(sys->fe),
			&(sys->vals),
            &(sys->ned),
			filename,
			sys->cond.directory,
            //sys->ned.num_edges,
            sys->fe.total_num_nodes,
			t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.V,
			filename,
			sys->cond.directory);
    
    filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_B);
    output_B_cell_vtk(
            &(sys->fe), &(sys->basis), &(sys->ned),
            sys->vals.Aphi_time,
            filename, sys->cond.directory);

}

double hot_start_read_initialize_val(
    double*     int_val,
    const char* input_fname,
    const char* directory)
{
    int BUFFER_SIZE = 1024;
    int total_num_nodes;
    int ndof;
    double t = 0.0;
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];

	fp = BBFE_sys_read_fopen(fp, "hot_start/start_time.dat", directory);
	fscanf(fp, "%s", id);
    fscanf(fp, "%lf", &(t));
    fclose(fp);

	fp = BBFE_sys_read_fopen(fp, input_fname, directory);
	fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes), &(ndof));
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fscanf(fp, "%lf", &(int_val[i * ndof + j]));
        }
    }
	fclose(fp);

    return t;
}

void hot_start_write_initialize_val(
    double*         int_val,
    const int       total_num_nodes,
    const int       ndof,
    const double    time,
    const char*     output_fname,
    const char*     directory)
{
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];

	fp = BBFE_sys_write_fopen(fp, output_fname, directory);
	fprintf(fp, "initialization\n");
    fprintf(fp, "%d %d\n", total_num_nodes, ndof);
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fprintf(fp, "%e ", int_val[i * ndof + j]);
        }
        fprintf(fp, "\n");        
    }
	fclose(fp);

    if(monolis_mpi_get_global_my_rank()==0){
        fp = BBFE_sys_write_fopen(fp, "hot_start/start_time.dat", directory);
        fprintf(fp, "start_time\n");
        fprintf(fp, "%lf", time);
        fclose(fp);
    }
}


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

    monolis_initialize(&(sys.monolis));
    
	BBFE_mag_pre_C(
			&(sys.fe), &(sys.basis), &(sys.ned), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

    const char* filename;
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem_bool.dat");
    set_elem_prop(
			&(sys.fe),
            &(sys.ned), 
            sys.cond.directory,
            filename);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

    FILE* fp;
	if(monolis_mpi_get_global_my_rank() == 0){
	    fp = BBFE_sys_write_fopen(fp, "NR_prop.txt", sys.cond.directory);
		fclose(fp);
	}

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));

	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));
    
    BBFE_mag_set_basis(
        &(sys.fe),
        &(sys.basis),
        &(sys.ned),
        sys.fe.local_num_nodes,
        sys.vals.num_ip_each_axis);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			BLOCK_SIZE);
   
    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
    BBFE_fluid_sups_read_Dirichlet_bc_NR(
            &(sys.bc_NR),
            filename,
            sys.cond.directory,
            sys.fe.total_num_nodes,
            BLOCK_SIZE);

    sys.vals.Aphi_time = (double *)calloc(sys.fe.total_num_nodes, sizeof(double));
    sys.vals.Aphi_time_curr = (double *)calloc(sys.fe.total_num_nodes, sizeof(double));

    /*for ROM *****************************************/
	
    /*for ROM input data*/

    ROM_read_args(argc, argv, &(sys.rom_prm_p));
	ROM_read_args(argc, argv, &(sys.rom_prm_v));

    ROM_set_param(
            &(sys.rom_p),
            sys.rom_prm_p.num_subdomains,
            sys.rom_prm_p.num_modes,
            sys.rom_prm_p.rom_epsilon,
            sys.rom_prm_p.solver_type);
    
    ROM_set_param(
            &(sys.rom_v),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);
    
    ROM_set_param(
            &(sys.rom_sups),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_p.num_modes + sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm_v.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm_v.solver_type);
/*
	ROM_std_hlpod_pre(
            &(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

	ROM_std_hlpod_pre(
            &(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

	ROM_std_hlpod_pre(
            &(sys.rom_sups),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
*/
            /******************/

    /*for offline******/
    /*
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);

    ROM_std_hlpod_offline_set_num_snapmat(
        &(sys.rom_p),
        //sys.vals.finish_time,
        1.0e-6,
        sys.vals.dt,
        sys.vals.snapshot_interval,
        1);
    
    ROM_std_hlpod_offline_set_num_snapmat(
        &(sys.rom_v),
        //sys.vals.finish_time,
        1.0e-6,
        sys.vals.dt,
        sys.vals.snapshot_interval,
        1);
*/
        /*
    ROM_std_hlpod_offline_set_num_snapmat(
        &(sys.rom_sups),
        //sys.vals.finish_time,
        1.0e-6,
        sys.vals.dt,
        sys.vals.snapshot_interval,
        1);
*/
	/****************** solver ********************/
    int step = 0;
    //double t = 1.5e-3;
    double t = 0.0;
	int file_num = 0;
    int count = 0;  //for ROM
    double t_hs = 0.0; //for hot start
    int step_hs = 0; //for hot start

    int nsteps = (int)ceil(sys.vals.finish_time / sys.vals.dt);

    /*
    char fnode[BUFFER_SIZE];
    snprintf(fnode, BUFFER_SIZE, "B_node_%06d.vtk", step);
    const char* fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
    output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, fn1, sys.cond.directory);
*/
    
/*
    char fname[BUFFER_SIZE];         
    snprintf(fname, BUFFER_SIZE, "hot_start/%s.%d.dat", "velosity_pressure", monolis_mpi_get_global_my_rank());
    t_hs = hot_start_read_initialize_val(sys.vals.Aphi_time, fname, sys.cond.directory);
    step_hs = t_hs / sys.vals.dt + 2;
    printf("Hot start time: %lf\n", t_hs);
    printf("Hot start step: %d\n", step_hs);
    t = t_hs;
    printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hs),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            1,
			1);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hs),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            1,
			1);
    */
    //monolis_copy_mat_nonzero_pattern_R(&(sys.monolis), &(sys.monolis_mass));

    //ja_init_states_zero();
    //ja_allocate_states(sys.fe.total_num_elems, sys.basis.num_integ_points);

    for (step = step_hs; step <= nsteps; ++step) {
        t += sys.vals.dt;

        printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
  
        solver_fom_NR_Aphi_team7(
            sys, t, count, 
            sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr,
            sys.fe.total_num_nodes);
/*
        log_accuracy_metrics_JA(
            &sys, sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr, step, t, sys.vals.dt, CURRENT_AMP);
*/    
        for(int i=0; i<sys.fe.total_num_nodes; ++i){
            sys.vals.Aphi_time[i] = sys.vals.Aphi_time_curr[i];
        }

        count ++;

        if(step%sys.vals.output_interval == 0) {
            //output_files(&sys, step, t);
            file_num += 1;
        }

        // 可視化・ログ
        
        if (step % sys.vals.output_interval == 0) {
            sys.vals.V   = BB_std_calloc_1d_double(sys.vals.V,   sys.fe.total_num_nodes);
            sys.vals.phi = BB_std_calloc_1d_double(sys.vals.phi, sys.fe.total_num_nodes);
            copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), sys.vals.Aphi_time, sys.vals.V, sys.vals.phi, sys.fe.total_num_elems);

            BB_std_free_1d_double(sys.vals.V, sys.fe.total_num_nodes);
            BB_std_free_1d_double(sys.vals.phi, sys.fe.total_num_nodes);

            char fnode[BUFFER_SIZE];
            snprintf(fnode, BUFFER_SIZE, "B_node_%06d.vtk", step);
            const char* fn1 = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fnode);
            output_B_node_vtk(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, fn1, sys.cond.directory);

            //output_files_nedelec(&sys, step, t);
            //output_files(&sys, step, t);
        }
        
/*
        if(step%1 == 0){
                char fname[BUFFER_SIZE];
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%d.dat", "velosity_pressure", monolis_mpi_get_global_my_rank());
                hot_start_write_initialize_val(sys.vals.Aphi_time, sys.fe.total_num_nodes, 1, t, fname, sys.cond.directory);
            }
        else{
        }
*/
    }

    
    ROM_std_hlpod_set_pod_modes_diag(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.monolis_com.n_internal_vertex,
		1,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);
    
    /*
    monolis_com_initialize_by_self(&(sys.mono_com0));
    ROM_std_hlpod_set_pod_modes_diag_mass(
        &(sys.monolis_mass),
        &(sys.mono_com0),
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.monolis_com.n_internal_vertex,
		1,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);
*/
    /*for writing vtk*/
    ROM_std_hlpod_read_pod_modes_diag(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.monolis_com.n_internal_vertex,
		1,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);
/*
	ROM_sys_hlpod_fe_write_pod_modes_vtk_diag(
		&(sys.fe),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		10,
		10,
		1,
		1,
		"pod_modes_vtk/pod_modes_v.vtk",
		"pod_modes_vtk/pod_modes_p.vtk",
		sys.cond.directory);
*/

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
