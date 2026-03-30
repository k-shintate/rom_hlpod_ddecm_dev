
#include "./../src/convdiff_core.h"
#include "./../src/nedelec_core.h"

#include "./../src/elemmat.h"
#include "./../src/shapefunc.h"
#include "./../src/std.h"

#include <complex.h>

#include "./../src/team11_core.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 4;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 200000;
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

/*
const double mu0 = 4.0*M_PI*1e-7; // H/m
const double Nu  = 1.0 / mu0;     // ← ここを μ0 に合わせる
const double Sigma_cop =  5.0*1.0e8;                  // 無次元化（典型導電率）
const double Sigma_air = 0.00000;                  // 無次元化（典型導電率）

const double freq  = 10.0;                 // [Hz]
const double Omega = 0;   // 角周波数
*/

const double mu0   = 4.0*M_PI*1e-7;      // H/m
const double Nu    = 1.0 / mu0;          // 1/H/m = A/(T·m)
const double Sigma_cop = 5.0e8;          // S/m
const double Sigma_air = 0.0;            // S/m

const double freq  = 50.0;               // Hz
const double Omega = 2.0*M_PI*freq;      // rad/s  ★これが正しい

typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double dt;
	double finish_time;
	int    output_interval;

	//double* T;
	double* error;
	double* theo_sol;

    double * V;
    double * phi;
    double _Complex * Aphi;
    double* Aphi_time;
    double* Aphi_time_curr;

} VALUES;

typedef struct
{
	const char* directory;

} CONDITIONS;

typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;
	BBFE_BC      bc;
    BBFE_BC      bc_internal;
    BBFE_BC      bc_NR;
	MONOLIS      monolis;
	MONOLIS_COM  monolis_com;

	CONDITIONS   cond;
	VALUES       vals;

    NEDELEC     ned;

} FE_SYSTEM;


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

static inline void ref_hex8_shape_grads(const double xi[3],
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

static inline void edge_midpoint_ref_coords(int p, double xi[3]) {
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

/* ============================================================
 * Constants & Material Properties (SI Units)
 * ============================================================ */
static const double SIGMA_COPPER = 5.77e7; /* [S/m] */
static const double SIGMA_CORE   = 3.72e3; /* [S/m] */
static const double FREQ_HZ      = 50.0;   /* [Hz] */
static const double CURRENT_AMP  = 40;    /* [A] */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double MU0    = 4.0 * M_PI * 1.0e-7;
static const double NU_LIN = 1.0 / (4.0 * M_PI * 1.0e-7); /* 1/mu0 */


/* ============================================================
 * Helper Functions & Structures
 * ============================================================ */

/* dot */
static inline double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* norm */
static inline double norm3(const double a[3]) {
    return sqrt(dot3(a,a));
}

/* --- Nonlinear Material Property (Brauer Law) --- */
/* nu(B) = 100 + 10 * exp(2 * |B|^2) */
static inline void eval_nu_and_dnudB(double Bmag, double* nu, double* dnudB){
    /* 真空の磁気抵抗率 (約 8.0e5) */
    const double NU0 = 1.0 / (4.0 * M_PI * 1.0e-7); 
    
    double B2 = Bmag * Bmag;
    
    /* 指数部が爆発しないようにクランプ (例: B=2.5T程度で制限) */
    /* exp(2 * 6.25) = exp(12.5) approx 2.6e5 */
    if (B2 > 6.25) B2 = 6.25; 

    double exp_val = exp(2.0 * B2);
    *nu = 100.0 + 10.0 * exp_val;
    *dnudB = 40.0 * Bmag * exp_val;

    /* オプション: nu が真空の透磁率を超えないように物理的なリミッタをかける
       (飽和後は真空の傾きに近づくのが一般的であるため) */
    if(*nu > NU0) {
        *nu = NU0;
        *dnudB = 0.0; // あるいは飽和域の傾き
    }
}
static inline double get_reluctivity_nu(double Bmag) {
    double B2 = Bmag * Bmag;
    return 100.0 + 10.0 * exp(2.0 * B2);
}

/* --- Coil Configuration Structure --- */
typedef struct {
    double axis[3];
    double center[3];
    double turns;
    double area;
    double phase_shift;
} COIL_INFO;

static const double MM_TO_M = 0.001;

/* ============================================================
 * Coil info (prop 1..3 are coils)
 * ============================================================ */
static inline int get_coil_info(int elem_prop, COIL_INFO* info) {
    info->axis[0] = 0.0; info->axis[1] = 1.0; info->axis[2] = 0.0; /* Y-axis */
    info->turns = 100.0;

    /* PI*(9^2 - 4.5^2) = 190.85 mm^2 -> m^2 */
    info->area  = 190.85 * (MM_TO_M * MM_TO_M);

    double cy = 30.0 * MM_TO_M;
    double cz = 20.0 * MM_TO_M;

    if (elem_prop == 1) {
        info->center[0] = 2.5 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = 0.0;
        return 1;
    } else if (elem_prop == 2) {
        info->center[0] = 35.0 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = -2.0 * M_PI / 3.0;
        return 1;
    } else if (elem_prop == 3) {
        info->center[0] = 67.5 * MM_TO_M;
        info->center[1] = cy;
        info->center[2] = cz;
        info->phase_shift = +2.0 * M_PI / 3.0;
        return 1;
    }
    return 0;
}

/* Air stabilization factors (tune as needed) */
static const double AIR_PHI_SIGMA_FACTOR  = 1.0e-5; /* phi-phi stabilization relative to copper sigma */
static const double AIR_A_MASS_SIGMA_FACTOR = 1.0e-5; /* A-A mass stabilization relative to copper sigma */

/* ============================================================
 * Material coefficient policy (prop: 1..3 coil, 4 core, 5 air)
 * ============================================================ */
static inline void get_sigmas_for_prop(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2 || prop == 3){
        /* coil (conductor) */
        *sigma_mass_A = SIGMA_COPPER * AIR_A_MASS_SIGMA_FACTOR;
        *sigma_cpl    = 0.0;
        *sigma_phi    = SIGMA_COPPER * AIR_PHI_SIGMA_FACTOR;
    } else if(prop == 4){
        /* core (A-only by default; phi is not solved here) */
        *sigma_mass_A = SIGMA_CORE;
        *sigma_cpl    = SIGMA_CORE;
        *sigma_phi    = SIGMA_CORE;
    } else if(prop == 5){
        /* air: coupling cut; add small penalties for solvability */
        *sigma_mass_A = SIGMA_COPPER * AIR_A_MASS_SIGMA_FACTOR;
        *sigma_cpl    = 0.0;
        *sigma_phi    = SIGMA_COPPER * AIR_PHI_SIGMA_FACTOR;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


/* ============================================================
 * B = curl(A) at integration point
 * ============================================================ */
static inline void compute_B_ip(
    const NEDELEC* ned, int e, int p,
    const double* x_curr, int nEdge,
    double B[3]
){
    B[0]=B[1]=B[2]=0.0;
    for(int m=0; m<ned->local_num_edges; ++m){
        int gm = ned->nedelec_conn[e][m];
        int sm = ned->edge_sign[e][m];
        const double* curlNm = ned->curl_N_edge[e][p][m];
        double am = x_curr[gm] * (double)sm;
        B[0] += am * curlNm[0];
        B[1] += am * curlNm[1];
        B[2] += am * curlNm[2];
    }
}

/* integration point coordinates */
static inline void get_interp_coords(
    int e, int p,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    double x_ip[3]
){
    x_ip[0] = 0.0; x_ip[1] = 0.0; x_ip[2] = 0.0;
    for(int n=0; n < fe->local_num_nodes; ++n){
        int global_node_id = fe->conn[e][n];
        double N_val = basis->N[p][n];
        x_ip[0] += N_val * fe->x[global_node_id][0];
        x_ip[1] += N_val * fe->x[global_node_id][1];
        x_ip[2] += N_val * fe->x[global_node_id][2];
    }
}

static inline void update_Aphi_NR(double* x_curr, const double* delta, int n_dof_total, double relaxation){
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] += relaxation * delta[i];
    }
}

/* ============================================================
 * Local utility functions
 * ============================================================ */
static double local_sum_sq(const double* v, int n) {
    double s = 0.0;
    for(int i=0; i<n; ++i) s += v[i]*v[i];
    return s;
}
static double local_max_abs(const double* v, int n) {
    double m = 0.0;
    for(int i=0; i<n; ++i){
        double a = fabs(v[i]);
        if(a > m) m = a;
    }
    return m;
}

/* ============================================================
 * Jacobian Assembly (Newton)
 * ============================================================ */
void set_element_mat_NR_Aphi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0/dt;
    const double epsB   = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];
        double sigma_mass_A, sigma_cpl, sigma_phi;
        //get_sigmas_for_prop(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);
        if(ned->elem_prop[e] == 1){ // 導体
            sigma_cpl = Sigma_cop;
            sigma_phi  = Sigma_cop;
            sigma_mass_A   = Sigma_cop;
        }
        else if(ned->elem_prop[e] == 2){ // 空気
            // 【重要1】連成項は完全カット (物理的に正しい)
            sigma_cpl = 0.0; 
            // 【重要2】φの行列が特異にならないようダミー値を入れる
            //  小さすぎると悪条件になるため、導体の 1/1000 程度を推奨
            //sigma_phi = Sigma_cop * 1.0e-12; 
            sigma_phi = 0,0; 
            // 【重要3】Aの行列安定化のために微小な質量項を入れる
            //  これによりBiCGStabの挙動が改善する
            //sigma_mass_A   = Sigma_cop * 1.0e-12;
            sigma_mass_A   = 0.0;
        }
        else {
            exit(1);
            continue; // その他の要素
        }

        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        int nonlinear_mu = 0;
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ========= [1] A-A : Curl-Curl (tangent stiffness) ========= */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    if(nonlinear_mu){
                        double B[3]; compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        double Bmag = fmax(norm3(B), epsB);
                        double nu, dnudB;
                        eval_nu_and_dnudB(Bmag, &nu, &dnudB);
                        double alpha = dnudB / Bmag;
                        val_ip_C[p] = nu * dot3(ci,cj) + alpha * dot3(ci,B) * dot3(cj,B);
                    }else{
                        val_ip_C[p] = NU_LIN * dot3(ci,cj);
                    }
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, (double)(si*sj) * v);
            }
        }

        /* ========= [2] A-A : Mass (sigma_mass_A/dt) ========= */
        //if(sigma_mass_A > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0,
                        (double)(si*sj) * v * inv_dt);
                }
            }
        //}

        /* ========= [3] Phi-Phi : Laplace (sigma_phi * gradN·gradN) ========= */
        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;

        //if(sigma_laplace > 0.0){
            for(int i=0; i<fe->local_num_nodes; ++i){
                int gi = fe->conn[e][i];
                for(int j=0; j<fe->local_num_nodes; ++j){
                    int gj = fe->conn[e][j];
                    for(int p=0; p<np; ++p){
                        /* sigma_phi -> sigma_laplace に変更 */
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], sigma_laplace);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v);
                }
            }
        //}

        /* ========= [4] A-Phi : Coupling C ========= */
        //if(sigma_cpl > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){      /* col: node (phi) */
                int gn = fe->conn[e][n];
                for(int j=0; j<ned->local_num_edges; ++j){ /* row: edge (A) */
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], ned->N_edge[e][p][j], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v);
                }
            }
        //}

        /* ========= [5] Phi-A : Coupling C^T / dt ========= */
        //if(sigma_cpl > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){      /* col: edge (A) */
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int n=0; n<fe->local_num_nodes; ++n){  /* row: node (phi) */
                    int gn = fe->conn[e][n];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], fe->geo[e][p].grad_N[n], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gn, gi, 0, 0, (double)si * v * inv_dt);
                }
            }
        //}
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * ============================================================ */
void set_element_vec_NR_Aphi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0 / dt;
    const double epsB   = 1.0e-14;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;

        //get_sigmas_for_prop(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);
        //int nonlinear_mu = (prop == 4) ? 1 : 0;
        //int nonlinear_mu = (prop == 4) ? 1 : 0;

        if(ned->elem_prop[e] == 1){ // 導体
            sigma_cpl = Sigma_cop;
            sigma_phi  = Sigma_cop;
            sigma_mass_A   = Sigma_cop;
        }
        else if(ned->elem_prop[e] == 2){ // 空気
            // 【重要1】連成項は完全カット (物理的に正しい)
            sigma_cpl = 0.0; 
            // 【重要2】φの行列が特異にならないようダミー値を入れる
            //  小さすぎると悪条件になるため、導体の 1/1000 程度を推奨
            sigma_phi = 0.0;
            //sigma_phi = Sigma_cop * 1.0e-11; 
            // 【重要3】Aの行列安定化のために微小な質量項を入れる
            //  これによりBiCGStabの挙動が改善する
            sigma_mass_A   = 0.0;
            //sigma_mass_A   = Sigma_cop * 1.0e-11;
        }
        else {
            exit(1);
            continue; // その他の要素
        }

        int nonlinear_mu = 0;

        /* coil source uses prop 1..3 */
        COIL_INFO coil;
        int is_coil = 0;
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation) -> contributes to A residual (+)
         * ========================================================== */
        /*
        if(is_coil){
            double omega = 2.0 * M_PI * FREQ_HZ;
            double I_t   = CURRENT_AMP * sin(omega * current_time + coil.phase_shift);
            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p=0; p<np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double d[3] = { x_ip[0]-coil.center[0], x_ip[1]-coil.center[1], x_ip[2]-coil.center[2] };
                    double da = dot3(d, coil.axis);
                    double r[3] = { d[0]-da*coil.axis[0], d[1]-da*coil.axis[1], d[2]-da*coil.axis[2] };

                    double tdir[3] = {
                        coil.axis[1]*r[2] - coil.axis[2]*r[1],
                        coil.axis[2]*r[0] - coil.axis[0]*r[2],
                        coil.axis[0]*r[1] - coil.axis[1]*r[0]
                    };
                    double n_tdir = norm3(tdir);

                    double Js[3] = {0.0, 0.0, 0.0};
                    if(n_tdir > eps_r){
                        double inv_n = 1.0/n_tdir;
                        Js[0] = J_mag * tdir[0] * inv_n;
                        Js[1] = J_mag * tdir[1] * inv_n;
                        Js[2] = J_mag * tdir[2] * inv_n;
                    }
                    val_ip_C[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis->mat.R.B[gi] += (double)si * integ;
            }
        }
        */
        

        /* ==========================================================
         * [2] A-mass term:
         *   - conductor/core:  -(sigma/dt) M (A_curr - A_prev)
         *   - air stabilization: -(sigma_stab/dt) M (A_curr)  (NO history)
         * ========================================================== */
        //if(sigma_mass_A > 0.0){
            int use_history = (prop != 1); /* air uses NO history (penalty only) */

            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                double acc = 0.0;
                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj_e = ned->nedelec_conn[e][j];
                    int gj = gj_e;
                    int sj = ned->edge_sign[e][j];

                    double coeffA;
                    if(use_history){
                        coeffA = (x_curr[gj] - x_prev[gj]); /* (A^{n+1}-A^n) */
                    }else{
                        coeffA = (x_curr[gj] - x_prev[gj]);
                        //coeffA = x_curr[gj];
                    }

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double mij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)(si*sj) * mij * coeffA * inv_dt;
                }
                monolis->mat.R.B[gi] -= acc; /* move to residual */
            }
        //}

        /* ==========================================================
         * [3] Curl-curl stiffness term: -(K(A_curr) * A_curr)
         * ========================================================== */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            double acc = 0.0;
            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];
                double a_val = x_curr[gj];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    double nu_val;
                    if(nonlinear_mu){
                        double B[3]; compute_B_ip(ned, e, p, x_curr, nEdge, B);
                        nu_val = get_reluctivity_nu(fmax(norm3(B), epsB));
                    }else{
                        nu_val = NU_LIN;
                    }
                    val_ip_C[p] = nu_val * dot3(ci, cj);
                }
                double kij = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                acc += (double)(si*sj) * kij * a_val;
            }
            monolis->mat.R.B[gi] -= acc;
        }

        /* ==========================================================
         * [4] A-Phi coupling in A-equation: -( sigma_cpl * C * phi )
         * ========================================================== */
        //if(sigma_cpl > 0.0){
            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                double acc = 0.0;
                for(int n=0; n<fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_val = x_curr[gn];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], ned->N_edge[e][p][j], sigma_cpl);
                    }
                    double cjn = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)sj * cjn * phi_val;
                }
                monolis->mat.R.B[gj] -= acc;
            }
        //}

        /* ==========================================================
         * [5] Phi-equation:
         *    (a) - sigma_cpl/dt * C^T * (A_curr - A_prev)  (air: sigma_cpl=0)
         *    (b) - sigma_phi * Kphi * phi_curr             (air: sigma_phi>0 stabilization)
         * ========================================================== */

        /* (a) dA/dt coupling term */
        //if(sigma_cpl > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int i=0; i<ned->local_num_edges; ++i){
                    int gi = ned->nedelec_conn[e][i];
                    int si = ned->edge_sign[e][i];
                    double da = x_curr[gi] - x_prev[gi];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], fe->geo[e][p].grad_N[n], sigma_cpl);
                    }
                    double gin = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += (double)si * gin * da * inv_dt;
                }
                monolis->mat.R.B[gn] -= acc;
            }
        //}

        /* (b) phi-phi Laplace/penalty term */
        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;
        //if(sigma_laplace > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int m=0; m<fe->local_num_nodes; ++m){
                    int gm = fe->conn[e][m];
                    double phi_m = x_curr[gm];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], fe->geo[e][p].grad_N[m], sigma_laplace);
                    }
                    double knm = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    acc += knm * phi_m;
                }
                monolis->mat.R.B[gn] -= acc;
            }
        //}
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

void apply_dirichlet_bc_for_A_and_phi(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned)
{
    int num_nodes = fe->total_num_nodes;

    int* node_is_conductor = (int*)calloc(num_nodes, sizeof(int));

    /* 要素をループして、導体要素(prop=1,2,3,4)に含まれるノードにフラグを立てる */
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop == 2){
            /* 空気：まだ導体(1)と判定されていない場合のみ 空気(2) をセット */
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                if(node_is_conductor[gn] != 1){ 
                    node_is_conductor[gn] = 2;
                }
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 1){
            /* 導体：優先度高 (1) をセット */
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 1; 
            }
        }
    }

    //bool* is_dir_edge = NULL;
    static int  is_dir_edge_n = 0;
    //if(!is_dir_edge){
        is_dir_edge_n  = fe->total_num_nodes;
        //is_dir_edge   = calloc((size_t)is_dir_edge_n, sizeof(bool));
        bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
        build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);
    //}

    /* --------------------------------------------------------
       2. A (Edge) に対する境界条件: A_tan = 0
       外部境界(bc->D_bc_existsが両端でtrue)のエッジを0固定
       -------------------------------------------------------- */
    
    /* 要素タイプに応じたエッジテーブル (Tet/Hex) */
    const int nen = fe->local_num_nodes;
    int n_local_edges = 0;
    const int (*edge_tbl)[2] = NULL;
    
    if (nen == 4) {
        n_local_edges = 6;
        edge_tbl = tet_edge_conn; /* bbfe_std内で定義されていると仮定 */
    } else if (nen == 8) {
        n_local_edges = 12;
        edge_tbl = hex_edge_conn;
    }

    
    for (int e = 0; e < fe->total_num_elems; ++e){
        for (int i = 0; i < n_local_edges; ++i){
            int n1_local = edge_tbl[i][0];
            int n2_local = edge_tbl[i][1];

            int gn1 = fe->conn[e][n1_local];
            int gn2 = fe->conn[e][n2_local];
            if (bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2]) {
                int ged = ned->nedelec_conn[e][i];                 

                if(!is_dir_edge[ged]) continue; 

                monolis_set_Dirichlet_bc_R(
                    monolis, 
                    monolis->mat.R.B, 
                    ged, 
                    0,    
                    0.0  
                );
            }
        }
    }

    /* --------------------------------------------------------
       3. phi (Node) に対する境界条件
       
       Case A: 空気領域 (絶縁体) -> phi = 0 で固定 (特異性回避)
       Case B: 導体領域の境界 -> 接地などが必要なら設定
       -------------------------------------------------------- */
/*
       for (int i = 0; i < num_nodes; ++i){
        if (node_is_conductor[i] == 2) {
            monolis_set_Dirichlet_bc_R(
                monolis, 
                monolis->mat.R.B, 
                i, 
                0, 
                0.0
            );
        }
        
    }
*/
    free(node_is_conductor);
}

static void debug_max_B_and_nu_core(
    FE_SYSTEM* sys,
    const double* x_curr,
    int n_dof_total,
    int it, int step, double t,
    const char* directory
){
    BBFE_DATA*  fe   = &(sys->fe);
    BBFE_BASIS* basis= &(sys->basis);
    NEDELEC*    ned  = &(sys->ned);

    const int np   = basis->num_integ_points;
    const double epsB = 1.0e-14;

    /* local maxima */
    double local_max_B  = 0.0;
    double local_max_nu = 0.0;

    /* NaN/Inf counter (doubleでreduce) */
    double local_bad_nu = 0.0;

    double local_max_x = 0.0;
    for(int i=0;i<n_dof_total;++i){
        double a = fabs(x_curr[i]);
        if(a > local_max_x) local_max_x = a;
    }

    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop != 4) continue; /* coreだけ */

        for(int p=0; p<np; ++p){
            double B[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B);

            double Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
            if(Bmag > local_max_B) local_max_B = Bmag;

            double nu = get_reluctivity_nu(fmax(Bmag, epsB));
            if(!isfinite(nu)){
                local_bad_nu += 1.0;
            } else {
                if(nu > local_max_nu) local_max_nu = nu;
            }
        }
    }

    /* MPI reduce (MAX/SUM) */
    double g_max_B  = local_max_B;
    double g_max_nu = local_max_nu;
    double g_bad_nu = local_bad_nu;
    double g_max_x  = local_max_x;

    monolis_allreduce_R(1, &g_max_B,  MONOLIS_MPI_MAX, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_max_nu, MONOLIS_MPI_MAX, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_bad_nu, MONOLIS_MPI_SUM, sys->monolis_com.comm);
    monolis_allreduce_R(1, &g_max_x,  MONOLIS_MPI_MAX, sys->monolis_com.comm);

    if(sys->monolis_com.my_rank == 0){
        printf("[DBG step=%d it=%d t=%.6e] core: max|x|=%.6e  max|B|=%.6e  max(nu)=%.6e  bad_nu(count)=%.0f\n",
               step, it, t, g_max_x, g_max_B, g_max_nu, g_bad_nu);
    }

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = BBFE_sys_write_add_fopen(fp, "NR_prop.txt", sys->cond.directory);
		fprintf(fp, "[DBG step=%d it=%d t=%.6e] core: max|x|=%.6e  max|B|=%.6e  max(nu)=%.6e  bad_nu(count)=%.0f\n",
               step, it, t, g_max_x, g_max_B, g_max_nu, g_bad_nu);
		fclose(fp);
	}
}

void compute_lumped_mass_air(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    double* M_lump_air   // [fe->total_num_nodes] 事前にcallocで0初期化
){
    const int np = basis->num_integ_points;
    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    for(int e=0; e<fe->total_num_elems; ++e){
        //if(ned->elem_prop[e] == 1){

            BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

            for(int a=0; a<fe->local_num_nodes; ++a){
                const int g = fe->conn[e][a];
                double acc = 0.0;
                for(int p=0; p<np; ++p){
                    acc += basis->integ_weight[p] * J_ip[p] * basis->N[p][a]; // ∫ N_a dV
                }
                M_lump_air[g] += acc;
            }
        //}
    }

    BB_std_free_1d_double(J_ip, np);
}

void add_coulomb_gauge_penalty_air(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* M_lump_air, // [total_num_nodes]
    double gamma              // 係数
){
    const int np = basis->num_integ_points;
    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    // local buffers
    double D[fe->local_num_nodes][ned->local_num_edges]; // D[a][i] = -∫ gradN_a · w_i dV  (air element)

    for(int e=0; e<fe->total_num_elems; ++e){
        if(ned->elem_prop[e] == 1||ned->elem_prop[e] == 2||ned->elem_prop[e] == 3||ned->elem_prop[e] == 5){

            BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

            // --- D[a][i] を作る ---
            for(int a=0; a<fe->local_num_nodes; ++a){
                for(int i=0; i<ned->local_num_edges; ++i){
                    double acc = 0.0;
                    for(int p=0; p<np; ++p){
                        const double* gNa = fe->geo[e][p].grad_N[a];      // [3]
                        const double* wi  = ned->N_edge[e][p][i];         // [3]
                        const double dot  = gNa[0]*wi[0] + gNa[1]*wi[1] + gNa[2]*wi[2];
                        acc += basis->integ_weight[p] * J_ip[p] * dot;
                    }
                    D[a][i] = -acc;
                }
            }

            // --- K_ij += gamma Σ_a D[a][i]D[a][j]/M_lump_air[g(a)] ---
            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = ned->edge_sign[e][i];

                for(int j=0; j<ned->local_num_edges; ++j){
                    const int gj = ned->nedelec_conn[e][j];
                    const int sj = ned->edge_sign[e][j];

                    double kij = 0.0;
                    for(int a=0; a<fe->local_num_nodes; ++a){
                        const int gnode = fe->conn[e][a];
                        const double m  = M_lump_air[gnode];
                        if(m > 0.0){
                            kij += (D[a][i] * D[a][j]) / m;
                        }
                    }

                    kij *= gamma;

                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis, gi, gj, 0, 0, (double)(si * sj) * kij
                    );
                }
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
}

void add_coulomb_gauge_penalty_air_matrix_and_residual(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* M_lump_air,
    double gamma,
    const double* x_curr, // [total_num_edges] 現在の解ベクトル (A)
    double* R_vec         // [total_num_edges] 残差ベクトル (b - Ax)
){
    const int np = basis->num_integ_points;
    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    // D[a][i] バッファ
    double D[fe->local_num_nodes][ned->local_num_edges]; 

    for(int e=0; e<fe->total_num_elems; ++e){
        if(ned->elem_prop[e] == 1||ned->elem_prop[e] == 2||ned->elem_prop[e] == 3||ned->elem_prop[e] == 5){ // Air only

            BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

            // --- D[a][i] の計算 (変更なし) ---
            for(int a=0; a<fe->local_num_nodes; ++a){
                for(int i=0; i<ned->local_num_edges; ++i){
                    double acc = 0.0;
                    for(int p=0; p<np; ++p){
                        const double* gNa = fe->geo[e][p].grad_N[a];
                        const double* wi  = ned->N_edge[e][p][i];
                        const double dot  = gNa[0]*wi[0] + gNa[1]*wi[1] + gNa[2]*wi[2];
                        acc += basis->integ_weight[p] * J_ip[p] * dot;
                    }
                    D[a][i] = -acc;
                }
            }

            // --- 行列と残差への加算 ---
            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = ned->edge_sign[e][i];
                
                double row_sum_val = 0.0; // 残差計算用： (K_local * x_local)_i

                for(int j=0; j<ned->local_num_edges; ++j){
                    const int gj = ned->nedelec_conn[e][j];
                    const int sj = ned->edge_sign[e][j];

                    // K_ij の計算
                    double kij = 0.0;
                    for(int a=0; a<fe->local_num_nodes; ++a){
                        const int gnode = fe->conn[e][a];
                        const double m  = M_lump_air[gnode];
                        if(m > 0.0){
                            kij += (D[a][i] * D[a][j]) / m;
                        }
                    }
                    kij *= gamma;

                    // 1. 行列への加算 (Jacobian)
                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis, gi, gj, 0, 0, (double)(si * sj) * kij
                    );

                    // 2. 残差計算のための準備
                    // 現在の解 x_j を取得 (MPIの場合は通信済みの値が必要)
                    double x_val = x_curr[gj]; 
                    
                    // 行列ベクトル積 K_ij * x_j を加算
                    // 注意: edge_sign (si, sj) を考慮する
                    row_sum_val += ((double)(si * sj) * kij) * x_val;
                }

                // 3. 残差ベクトルへの加算 (Residual)
                // 通常、Newton法の式は J*dx = -R なので、
                // R = (Internal Force) - (External Force)
                // ペナルティ項は Internal Force (剛性) の一部として「加算」します。
                
                //#pragma omp atomic
                R_vec[gi] += row_sum_val; 
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
}

// エア領域に「発熱しないバネ（微小質量）」を付加する
void add_stiffness_penalty_air(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    double alpha,         // ペナルティ係数（後述）
    const double* x_curr, // 現在の解
    double* R_vec         // 残差ベクトル
){
    const int np = basis->num_integ_points;
    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    for(int e=0; e<fe->total_num_elems; ++e){
        // エア領域のみ対象
        if(ned->elem_prop[e] == 1 || ned->elem_prop[e] == 2 || ned->elem_prop[e] == 3 || ned->elem_prop[e] == 5){

            BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = ned->edge_sign[e][i];
                
                double row_sum = 0.0;

                for(int j=0; j<ned->local_num_edges; ++j){
                    const int gj = ned->nedelec_conn[e][j];
                    const int sj = ned->edge_sign[e][j];

                    // --- 通常の質量行列項の計算 ∫(wi・wj)dV ---
                    double mij = 0.0;
                    for(int p=0; p<np; ++p){
                        const double* wi = ned->N_edge[e][p][i];
                        const double* wj = ned->N_edge[e][p][j];
                        double dot = wi[0]*wj[0] + wi[1]*wj[1] + wi[2]*wj[2];
                        
                        mij += basis->integ_weight[p] * J_ip[p] * dot;
                    }
                    
                    // 係数をかける
                    double kij_penalty = alpha * mij;

                    // 1. 剛性行列(K)に加算
                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis, gi, gj, 0, 0, (double)(si * sj) * kij_penalty
                    );

                    // 2. 残差計算 (K * x)
                    // 注: あなたのソルバーの流儀に合わせて += か -= を選ぶ
                    // 一般的には、K項と同じ符号で足す
                    row_sum += ((double)(si * sj) * kij_penalty) * x_curr[gj];
                }

                // 残差ベクトルへの加算
                // R = F - Kx なら -=
                // R = Kx - F なら +=
                // 前回のテストで「-=」の方が正しそうならそちらに合わせてください
                // #pragma omp atomic
                R_vec[gi] += row_sum; 
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
}

double Bz_of_time(
    double t)
{
    return 1;
}

void solver_fom_NR_Aphi_team11(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            dx[i] = 0.0;
        }


        // ディリクレ境界（毎ステップ適用）

        double val = Bz_of_time(t);
        apply_nedelec_boundary_conditions_TEAM11_NR(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned), val, t, x_curr);

        /* 組み立て */
        set_element_mat_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        /* 線形解は dx に入れる（x_curr と分ける！） */
        monowrap_solve_R(
            &(sys.monolis),
            &(sys.monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        /* Newton update */
        update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), dx, A_delta, phi_delta, sys.fe.total_num_elems);

        /*
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));
        
        apply_nedelec_boundary_conditions_TEAM11_time(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned), val, t);
        */

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        //int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v);
        //int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        if (conv_v && conv_p) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }

    free(dx);
    free(rvec);


    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
}


void solver_fom_NR_Aphi(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total
){
    const int max_iter = 20;
    const double relaxation = 1.0;

    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double* rvec_old = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double* A_delta   = BB_std_calloc_1d_double(A_delta,   sys.fe.total_num_nodes);
    double* phi_delta = BB_std_calloc_1d_double(phi_delta, sys.fe.total_num_nodes);
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);
    double* A_old   = BB_std_calloc_1d_double(A_old,   sys.fe.total_num_nodes);
    double* phi_old = BB_std_calloc_1d_double(phi_old, sys.fe.total_num_nodes);
    copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_prev, A_old, phi_old, sys.fe.total_num_elems);

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }

    double r0 = -1.0;
    const double rel_tol_v = 1.0e-6;
    const double abs_tol_v = 1.0e-12;
    const double rel_tol_p = 1.0e-6;
    const double abs_tol_p = 1.0e-12;
    const double rel_tol_r = 1.0e-6;
    const double tiny      = 1.0e-30;
    int max_iter_NR = 20;

    for(int it=0; it<max_iter; ++it){

        debug_max_B_and_nu_core(&sys, x_curr, n_dof_total, it, step, t, sys.cond.directory);

        monolis_clear_mat_value_R(&(sys.monolis));
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
            dx[i] = 0.0;
        }

        /* 組み立て */
        set_element_mat_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_curr, sys.vals.dt);
        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        double* M_lump_air = BB_std_calloc_1d_double(M_lump_air, sys.fe.total_num_nodes);

        compute_lumped_mass_air(
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            M_lump_air);

        add_coulomb_gauge_penalty_air(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            M_lump_air,
            0); /* gamma */
        
        add_coulomb_gauge_penalty_air_matrix_and_residual(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            M_lump_air,
            0,
            x_curr,
            sys.monolis.mat.R.B);
        
        add_stiffness_penalty_air(
            &(sys.monolis),
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            0,
            x_curr,
            sys.monolis.mat.R.B);

        /* 残差ベクトルを保存（B = -F） */
        if(it==0){
            for(int i=0; i<n_dof_total; ++i) rvec_old[i] = sys.monolis.mat.R.B[i];
        }

        monowrap_solve_R(
            &(sys.monolis),
            &(sys.monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

        /* Newton update */
        update_Aphi_NR(x_curr, dx, n_dof_total, relaxation);

        /* 残差ベクトルを保存（B = -F） */
        for(int i=0; i<n_dof_total; ++i){
            sys.monolis.mat.R.B[i] = 0.0;
            sys.monolis.mat.R.X[i] = 0.0;
        }

        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), x_curr, A, phi, sys.fe.total_num_elems);
        copy_Aphi_to_V_phi_time2(&(sys.fe), &(sys.ned), dx, A_delta, phi_delta, sys.fe.total_num_elems);

        set_element_vec_NR_Aphi(&(sys.monolis), &(sys.fe), &(sys.basis), &(sys.ned),
                                x_prev, x_curr, sys.vals.dt, t);
        
        apply_dirichlet_bc_for_A_and_phi(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned));

        for(int i=0; i<n_dof_total; ++i) rvec[i] = sys.monolis.mat.R.B[i];

        /*収束判定 別の関数にまとめたい*/
        double norm_v = calc_internal_norm_1d(
            A,
            sys.monolis_com.n_internal_vertex,
            1);

        double norm_delta_v = calc_internal_norm_1d(
            A_delta,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r_old = calc_internal_norm_1d(
            rvec_old,
            sys.monolis_com.n_internal_vertex,
            1);
        
        double norm_r = calc_internal_norm_1d(
            rvec,
            sys.monolis_com.n_internal_vertex,
            1);

        /* 圧力の L2 ノルム（内部自由度のみ） */
        double norm_p = 0.0, norm_delta_p = 0.0;
        for (int ii = 0; ii < sys.monolis_com.n_internal_vertex; ++ii) {
            double pv  = phi[ii];
            double dpv = phi_delta[ii];
            norm_p       += pv  * pv;
            norm_delta_p += dpv * dpv;
        }

        /* L∞（最大変化量）：速度は3成分、圧力は1成分 */
        double linf_delta_v_local = 0.0;
        double linf_delta_p_local = 0.0;
        for (int i_node = 0; i_node < sys.monolis_com.n_internal_vertex; ++i_node) {
            double av = fabs(A_delta[i_node]);
            if (av > linf_delta_v_local) linf_delta_v_local = av;
        
            double ap = fabs(phi_delta[i_node]);
            if (ap > linf_delta_p_local) linf_delta_p_local = ap;
        }

        monolis_allreduce_R(1, &norm_v,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_v, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_p,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_delta_p, MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r,       MONOLIS_MPI_SUM, sys.monolis_com.comm);
        monolis_allreduce_R(1, &norm_r_old, MONOLIS_MPI_SUM, sys.monolis_com.comm);

        double linf_delta_v = linf_delta_v_local;
        double linf_delta_p = linf_delta_p_local;
        monolis_allreduce_R(1, &linf_delta_v, MONOLIS_MPI_MAX, sys.monolis_com.comm);
        monolis_allreduce_R(1, &linf_delta_p, MONOLIS_MPI_MAX, sys.monolis_com.comm);

        double nrm_v        = sqrt(norm_v);
        double nrm_dv       = sqrt(norm_delta_v);
        double nrm_p        = sqrt(norm_p);
        double nrm_dp       = sqrt(norm_delta_p);
        double nrm_r        = sqrt(norm_r);
        double nrm_r_old       = sqrt(norm_r_old);

        double denom_v = fmax(nrm_v,  tiny);
        double denom_p = fmax(nrm_p,  tiny);

        int conv_v = (nrm_dv <= abs_tol_v) || (nrm_dv/denom_v <= rel_tol_v) || (linf_delta_v <= abs_tol_v) || (nrm_r/nrm_r_old <= rel_tol_r);
        int conv_p = (nrm_dp <= abs_tol_p) || (nrm_dp/denom_p <= rel_tol_p) || (linf_delta_p <= abs_tol_p) || (nrm_r/nrm_r_old <= rel_tol_r);

        if(monolis_mpi_get_global_my_rank() == 0){
            printf("[NR %2d] ||dv||2=%.3e  ||v||2=%.3e  rel=%.3e  Linf(dv)=%.3e\n",
                it, nrm_dv, nrm_v, nrm_dv/denom_v, linf_delta_v);
            printf("[NR %2d] ||dp||2=%.3e  ||p||2=%.3e  rel=%.3e  Linf(dp)=%.3e  |r_old| = %.3e |r| = %.3e  |r|/|r_old| = %.3e\n", 
                it, nrm_dp, nrm_p, nrm_dp/denom_p, linf_delta_p, nrm_r_old, nrm_r, nrm_r / nrm_r_old);
        }

        if (conv_v && conv_p) {
            double max_du = 0.0;
            for (int ii = 0; ii < sys.fe.total_num_nodes; ++ii) {    
                double du = fabs(A[ii] - A_old[ii]);
                    if (du > max_du) max_du = du;
            }
            if(monolis_mpi_get_global_my_rank() == 0){
                printf("[step %d] max|v^{n+1}-v^{n}| = %.6e\n", step, max_du);
            }

            ROM_BB_vec_copy_1d(
                A,
                A_old,
                sys.fe.total_num_nodes,
                1);

            break;
        }

    }

/*
    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_std_hlpod_set_snapmat_nobc(
                    phi,
                    &(sys.rom_p.hlpod_mat),
                    sys.fe.total_num_nodes,
                    1,
                    (int)(step/sys.vals.snapshot_interval));
        }
        else{
            ROM_std_hlpod_set_snapmat_nobc(
                    phi,
                    &(sys.rom_p.hlpod_mat),
                    sys.mono_com.n_internal_vertex,
                    1,
                    (int)(step/sys.vals.snapshot_interval));
        }

    }

    if(step%sys.vals.snapshot_interval == 0) {
        printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));

        if(monolis_mpi_get_global_comm_size() == 1){
            ROM_sys_hlpod_fe_set_snap_mat_para_ned(
                    A,
                    &(sys.fe),
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc),
                    &(sys.ned),
                    sys.fe.total_num_nodes,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
        }
        else{
            ROM_sys_hlpod_fe_set_snap_mat_para_ned(
                    A,
                    &(sys.fe),
                    &(sys.rom_v.hlpod_mat),
                    &(sys.bc),
                    &(sys.ned),
                    sys.fe.total_num_nodes,
                    1,
                    ((int)step/sys.vals.snapshot_interval));
            
        }

    }
    */

    free(dx);
    free(rvec);

    BB_std_free_1d_double(A_delta,   sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_delta, sys.fe.total_num_nodes);
    BB_std_free_1d_double(A,         sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi,       sys.fe.total_num_nodes);
    BB_std_free_1d_double(A_old,     sys.fe.total_num_nodes);
    BB_std_free_1d_double(phi_old,   sys.fe.total_num_nodes);
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

static void log_accuracy_metrics(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    int step,
    double t,
    double dt,
    double P_input
){
    BBFE_DATA* fe     = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC* ned      = &(sys->ned);

    const int np    = basis->num_integ_points;
    const int nEdge = fe->total_num_nodes;
    const double inv_dt = 1.0 / dt;

    /* --- Probe Settings (Monitoring Points) --- */
    /* 鉄心内の飽和を確認したい座標を指定 (単位: m) */
    const double search_radius = 2.0e-3; // 探索半径 2mm
    const double target_pos1[3] = { 0.0025, 0.0350, 0.0200 }; // 測定点1
    const double target_pos2[3] = { 0.0350, 0.0350, 0.0200 }; // 測定点2
    const double target_pos3[3] = { 0.0675, 0.0350, 0.0200 }; // 測定点3

    /* --- Data Accumulators Definition ---
     * Size: 29
     *
     * [MPI_SUM target] (Integrals)
     * 0-2:  (Reserved/Unused)
     * 3-5:  V_fem (EMF) for U, V, W
     * 6:    Total Loss
     * 7:    Total Mag Energy
     * 8-10: Mag Energy (Coil, Core, Air)
     * 11-13: Volume (Coil, Core, Air)
     * 16-18: Loss Breakdown (Coil, Core, Air)
     *
     * [MPI_MAX target] (Point Values)
     * 14:   Max B_mag (Global)
     * 15:   Max Nu (Air Region check)
     * 20:   Probe1 Nu
     * 21:   Probe1 Distance
     * 22:   Probe1 B_mag
     * 23-25: Probe2 (Nu, Dist, B)
     * 26-28: Probe3 (Nu, Dist, B)
     */
    const int N_LOG = 29;
    double local_vals[29];
    for(int i=0; i<N_LOG; i++) local_vals[i] = 0.0;

    /* 真空の磁気抵抗率 */
    const double NU0 = 1.0 / (4.0 * M_PI * 1.0e-7);

    /* 電流源計算 (Source Current) */
    double I_src[3] = {0.0};
    double omega = 2.0 * M_PI * FREQ_HZ;
    for(int k=0; k<3; k++){
        COIL_INFO info;
        get_coil_info(k+1, &info);
        I_src[k] = CURRENT_AMP * sin(omega * t + info.phase_shift);
    }

    /* 内部要素数読み込み (Auxiliary Data) */
    const char* fname;
    FILE* fp_in;
    char id[128];
    int tmp;
    int num_internal_elems = 0;
    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "graph_elem.dat.n_internal");
    fp_in = BBFE_sys_read_fopen(fp_in, fname, sys->cond.directory);
    fscanf(fp_in, "%s %d", id, &(tmp));
    fscanf(fp_in, "%d", &(num_internal_elems));
    fclose(fp_in);

    /* --- Element Loop --- */
    for(int e=0; e<num_internal_elems; ++e){
        int prop = ned->elem_prop[e];

        /* Region Identification */
        int region_type = 0; // 0:Other, 1:Coil, 2:Core, 3:Air
        if(prop == 1 || prop == 2 || prop == 3) region_type = 1;
        else if(prop == 4) region_type = 2;
        else if(prop == 5) region_type = 3;

        double sigma_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop(prop, &sigma_A, &sigma_cpl, &sigma_phi);

        /* Coil Info */
        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);
        int phase_idx = -1;
        if(prop == 1) phase_idx = 0;
        if(prop == 2) phase_idx = 1;
        if(prop == 3) phase_idx = 2;

        /* Calculate Element Centroid for Probing */
        double center[3] = {0.0};
        for(int n=0; n<fe->local_num_nodes; ++n){
            int node_id = fe->conn[e][n];
            center[0] += fe->x[node_id][0];
            center[1] += fe->x[node_id][1];
            center[2] += fe->x[node_id][2];
        }
        center[0] /= fe->local_num_nodes;
        center[1] /= fe->local_num_nodes;
        center[2] /= fe->local_num_nodes;

        /* Distance to Probes */
        double dist1 = sqrt(pow(center[0]-target_pos1[0],2) + pow(center[1]-target_pos1[1],2) + pow(center[2]-target_pos1[2],2));
        double dist2 = sqrt(pow(center[0]-target_pos2[0],2) + pow(center[1]-target_pos2[1],2) + pow(center[2]-target_pos2[2],2));
        double dist3 = sqrt(pow(center[0]-target_pos3[0],2) + pow(center[1]-target_pos3[1],2) + pow(center[2]-target_pos3[2],2));

        /* Integration Loop */
        for(int p=0; p<np; ++p){
            double w_detJ = basis->integ_weight[p] * fe->geo[e][p].Jacobian;
            
            /* Volume Accumulation */
            if(region_type == 1) local_vals[11] += w_detJ;
            if(region_type == 2) local_vals[12] += w_detJ;
            if(region_type == 3) local_vals[13] += w_detJ;

            /* B-Field Calculation */
            double B_curr[3], B_prev[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B_curr);
            compute_B_ip(ned, e, p, x_prev, nEdge, B_prev);

            double Bmag = norm3(B_curr);
            
            /* Max B Check (Global) */
            if(Bmag > local_vals[14]) local_vals[14] = Bmag;

            /* Reluctivity (Nu) Calculation */
            double current_nu;
            if (region_type == 2) { 
                /* Core: Non-linear */
                current_nu = get_reluctivity_nu(fmax(Bmag, 1.0e-14));
            } else { 
                /* Air/Coil: Linear (Vacuum) */
                current_nu = NU0;
            }

            /* Max Nu Check (Air only, for stability check) */
            if(region_type == 3) {
                if(current_nu > local_vals[15]) local_vals[15] = current_nu;
            }

            /* --- Probe Data Capture (Only in Core) --- */
            if(region_type == 2){
                // Probe 1
                if(dist1 < search_radius){
                    local_vals[20] = current_nu;
                    local_vals[21] = dist1;
                    local_vals[22] = Bmag;
                }
                // Probe 2
                if(dist2 < search_radius){
                    local_vals[23] = current_nu;
                    local_vals[24] = dist2;
                    local_vals[25] = Bmag;
                }
                // Probe 3
                if(dist3 < search_radius){
                    local_vals[26] = current_nu;
                    local_vals[27] = dist3;
                    local_vals[28] = Bmag;
                }
            }

            /* A-field Reconstruction & E-field Calculation */
            double A_curr[3]={0}, A_prev[3]={0};
            double grad_phi[3]={0};
            
            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                const double* N = ned->N_edge[e][p][i];
                double ac = x_curr[gi]*si;
                double ap = x_prev[gi]*si;
                for(int d=0; d<3; d++){ A_curr[d]+=ac*N[d]; A_prev[d]+=ap*N[d]; }
            }
            if(sigma_cpl > 0.0 || sigma_phi > 0.0){
                for(int n=0; n<fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double ph = x_curr[gn];
                    const double* dN = fe->geo[e][p].grad_N[n];
                    for(int d=0; d<3; d++) grad_phi[d] += ph*dN[d];
                }
            }

            double E_vec[3];
            for(int d=0; d<3; d++) E_vec[d] = -(A_curr[d]-A_prev[d])*inv_dt - grad_phi[d];

            /* Coil Source Term (Joule Heating Source Term & EMF Integration) */
            if(is_coil && phase_idx >= 0){
                double x_ip[3]; 
                get_interp_coords(e, p, fe, basis, x_ip);
                
                double d_vec[3] = { x_ip[0]-coil.center[0], x_ip[1]-coil.center[1], x_ip[2]-coil.center[2] };
                double da = dot3(d_vec, coil.axis);
                double r_vec[3] = { d_vec[0]-da*coil.axis[0], d_vec[1]-da*coil.axis[1], d_vec[2]-da*coil.axis[2] };
                double t_vec[3] = { coil.axis[1]*r_vec[2]-coil.axis[2]*r_vec[1], coil.axis[2]*r_vec[0]-coil.axis[0]*r_vec[2], coil.axis[0]*r_vec[1]-coil.axis[1]*r_vec[0] };
                double t_norm = norm3(t_vec);
                
                if(t_norm > 1.0e-12){
                    double J0_mag = coil.turns / coil.area;
                    double J0[3] = { J0_mag*t_vec[0]/t_norm, J0_mag*t_vec[1]/t_norm, J0_mag*t_vec[2]/t_norm };
                    /* EMF Integration (V_fem) */
                    local_vals[3 + phase_idx] += dot3(E_vec, J0) * w_detJ;
                }
            }

            /* --- Joule Loss Calculation --- */
            double dLoss = 0.0;
            if (sigma_A > 0.0) {
                dLoss = sigma_A * dot3(E_vec, E_vec) * w_detJ;
            }

            local_vals[6] += dLoss; /* Total Loss */

            /* Loss Breakdown */
            if(region_type == 1) local_vals[16] += dLoss; // Coil
            if(region_type == 2) local_vals[17] += dLoss; // Core
            if(region_type == 3) local_vals[18] += dLoss; // Air

            /* --- Magnetic Energy Calculation --- */
            double H_vec[3] = { current_nu*B_curr[0], current_nu*B_curr[1], current_nu*B_curr[2] };
            double dB_dt[3] = { (B_curr[0]-B_prev[0])*inv_dt, (B_curr[1]-B_prev[1])*inv_dt, (B_curr[2]-B_prev[2])*inv_dt };
            double dW_val = dot3(H_vec, dB_dt) * w_detJ;

            local_vals[7] += dW_val; /* Total Mag Energy Change */
            
            if(region_type == 1) local_vals[8] += dW_val; // Coil
            if(region_type == 2) local_vals[9] += dW_val; // Core
            if(region_type == 3) local_vals[10] += dW_val; // Air

        } // End Integ Loop
    } // End Element Loop

    /* --- MPI Reduce (Split logic for SUM and MAX) --- */
    double global_vals[29];

    /* 1. 集計 (SUM) 用の一時配列作成 */
    /* 対象: 0-13, 16-19 (積分値) */
    double val_sum[29];
    for(int i=0; i<29; i++) val_sum[i] = local_vals[i];
    
    // MAX対象の変数は 0 にして SUM に混ざらないようにする
    val_sum[14] = 0.0; val_sum[15] = 0.0; 
    for(int i=20; i<29; i++) val_sum[i] = 0.0; 
    
    monolis_allreduce_R(29, val_sum, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    double val_max[29];
    for(int i=0; i<29; i++) val_max[i] = 0.0;

    val_max[14] = local_vals[14];
    val_max[15] = local_vals[15];
    for(int i=20; i<29; i++) val_max[i] = local_vals[i];

    monolis_allreduce_R(29, val_max, MONOLIS_MPI_MAX, sys->monolis_com.comm);

    /* 3. 結果の統合 */
    for(int i=0; i<29; i++) global_vals[i] = val_sum[i];
    
    // MAXで計算した部分を上書き
    global_vals[14] = val_max[14];
    global_vals[15] = val_max[15];
    for(int i=20; i<29; i++) global_vals[i] = val_max[i];


    /* --- Output Processing (Rank 0 Only) --- */
    if(sys->monolis_com.my_rank == 0){
        
        /* 誘導起電力 (Back EMF) */
        double V_fem[3]= { global_vals[3], global_vals[4], global_vals[5] };
        
        /* 端子電圧計算 */
        const double R_DC_COIL = 0.1; // [Ohm]
        double V_term[3];
        for(int k=0; k<3; k++){
            V_term[k] = R_DC_COIL * I_src[k] - V_fem[k]; 
        }

        /* プローブデータの取り出し */
        double p1_nu = global_vals[20]; double p1_B = global_vals[22];
        double p2_nu = global_vals[23]; double p2_B = global_vals[25];
        double p3_nu = global_vals[26]; double p3_B = global_vals[28];

        /* コンソール表示 (デバッグ用) */
        printf("  [Monitor] Step:%d t=%.4e | MaxB: %.2f [T] | P1(Core): B=%.2f Nu=%.2e| P2(Core): B=%.2f Nu=%.2e| P3(Core): B=%.2f Nu=%.2e\n", 
               step, t, global_vals[14], p1_B, p1_nu, p2_B, p2_nu, p3_B, p3_nu);

        /* CSV Output */
        FILE* fp_out;
        fp_out = BBFE_sys_write_add_fopen(fp_out, "terminal_voltage_log.csv", sys->cond.directory);
        
        /* Header Output (初回のみ) */
        if(step == 0){
             fprintf(fp_out, "Step,Time,I_u,I_v,I_w,V_t_u,V_t_v,V_t_w,EMF_u,EMF_v,EMF_w,Power_In,Max_B,Max_Nu,Loss_Coil,Loss_Core,Loss_Air,P1_Nu,P1_B,P2_Nu,P2_B,P3_Nu,P3_B\n");
        }

        /* Data Output */
        fprintf(fp_out, "%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                step, t,
                I_src[0], I_src[1], I_src[2],    // Current
                V_term[0], V_term[1], V_term[2], // Terminal Voltage
                -V_fem[0], -V_fem[1], -V_fem[2], // EMF (Induced Voltage)
                P_input,                         // Input Power
                global_vals[14],                 // Max B (Global)
                global_vals[15],                 // Max Nu (Air)
                global_vals[16], global_vals[17], global_vals[18], // Loss Breakdown
                p1_nu, p1_B,                     // Probe 1 Data
                p2_nu, p2_B,                     // Probe 2 Data
                p3_nu, p3_B                      // Probe 3 Data
        );
        fclose(fp_out);
    }
}

void output_result_B_node_sim(
    BBFE_DATA* fe, VALUES* vals, BBFE_BASIS* basis, NEDELEC* ned, const double* Aphi, double t,
    const char* directory)
{
    double** B_cell = BB_std_calloc_2d_double(B_cell, fe->total_num_elems, 3);
    double** B_node = BB_std_calloc_2d_double(B_node, fe->total_num_nodes, 3);

    compute_B_cell_average(fe, basis, ned, Aphi, B_cell);
    accumulate_B_cell_to_nodes(fe, B_cell, B_node);

    output_result_file_B_node_sim(fe, t, B_node, directory);

    BB_std_free_2d_double(B_cell, fe->total_num_elems, 3);
    BB_std_free_2d_double(B_node, fe->total_num_nodes, 3);
}

void output_result_B_node_sim_I(
    BBFE_DATA* fe, VALUES* vals, BBFE_BASIS* basis, NEDELEC* ned, const double* Aphi, double t,
    const char* directory)
{
    double** B_cell = BB_std_calloc_2d_double(B_cell, fe->total_num_elems, 3);
    double** B_node = BB_std_calloc_2d_double(B_node, fe->total_num_nodes, 3);

    compute_B_cell_average(fe, basis, ned, Aphi, B_cell);
    accumulate_B_cell_to_nodes(fe, B_cell, B_node);

    output_result_file_B_node_sim_I(fe, t, B_node, directory);

    BB_std_free_2d_double(B_cell, fe->total_num_elems, 3);
    BB_std_free_2d_double(B_node, fe->total_num_nodes, 3);
}


static inline void A_uniform(const double B[3], const double x[3], double A[3]){
    // 0.5 * (B x x)
    A[0] = 0.5*( B[1]*x[2] - B[2]*x[1] );
    A[1] = 0.5*( B[2]*x[0] - B[0]*x[2] );
    A[2] = 0.5*( B[0]*x[1] - B[1]*x[0] );
}

static inline void A_dipole_from_C(double C, const double x[3], double A[3]){
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r  = sqrt(r2);
    double r3 = r2*r;
    if(r3 < 1e-30){ A[0]=A[1]=A[2]=0.0; return; }

    // A_dip = -C * (ez x x)/r^3 = -C*(-y, x, 0)/r^3 = (C*y/r^3, -C*x/r^3, 0)
    A[0] =  C * x[1] / r3;
    A[1] = -C * x[0] / r3;
    A[2] =  0.0;
}

static inline double edge_integral_A_dot_dl(
    const double B[3], double C,
    const double x1[3], const double x2[3])
{
    double dx[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };

    // 2点Gauss on s in [0,1]
    const double gp[2] = { 0.5 - 0.2886751345948129, 0.5 + 0.2886751345948129 };
    const double gw[2] = { 0.5, 0.5 };

    double In = 0.0;
    for(int q=0;q<2;++q){
        double s = gp[q];
        double x[3] = { x1[0] + s*dx[0], x1[1] + s*dx[1], x1[2] + s*dx[2] };

        double Au[3], Ad[3], A[3];
        A_uniform(B, x, Au);
        A_dipole_from_C(C, x, Ad);
        A[0]=Au[0]+Ad[0]; A[1]=Au[1]+Ad[1]; A[2]=Au[2]+Ad[2];

        In += gw[q] * (A[0]*dx[0] + A[1]*dx[1] + A[2]*dx[2]);
    }
    return In;
}

void apply_nedelec_boundary_conditions_team6(
    MONOLIS*    monolis,
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    double      Bz_t,
    double      t_sec,
    double*      val,
    const char*    directory)
{
    const double b = 0.055;
    //const double C = team11_C_unitstep(t_sec, b) * (Bz_t / 1.0);
    const double B_vec[3] = {0.0, 0.0, Bz_t};

    const int nen = fe->local_num_nodes;
    int num_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if (nen == 4) {
        num_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if (nen == 8) {
        num_edges = 12;
        edge_tbl = hex_edge_conn;
    } else {
        fprintf(stderr, "Unsupported element type: local_num_nodes=%d\n", nen);
        exit(EXIT_FAILURE);
    }

    static int  is_dir_edge_n = 0;
    is_dir_edge_n  = fe->total_num_nodes;
    bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
    build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);


    for (int e = 0; e < fe->total_num_elems; ++e){

        for (int ed = 0; ed < num_edges; ++ed){

            const int ln1 = edge_tbl[ed][0];
            const int ln2 = edge_tbl[ed][1];
            const int gn1 = fe->conn[e][ln1];
            const int gn2 = fe->conn[e][ln2];

            const int gid = ned->nedelec_conn[e][ed];

            if (!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;

            if(!is_dir_edge[gid]) continue; 

            const double x1[3] = { fe->x[gn1][0], fe->x[gn1][1], fe->x[gn1][2] };
            const double x2[3] = { fe->x[gn2][0], fe->x[gn2][1], fe->x[gn2][2] };

            double Aint = edge_integral_A_dot_dl(B_vec, 0, x1, x2);

            //if (ned->edge_sign) {
                Aint *= ned->edge_sign[e][ed]; /* ±1 */
            //}

            const int ge = gid;

            //val[ge] = Aint;
            double _Complex val = Aint;

            //printf("%lf ", Aint);

            monolis_set_Dirichlet_bc_C(
                            monolis,
                            monolis->mat.C.B,
                            ge,
                            0,
                            val);
        }
    }


    int num_nodes = fe->total_num_nodes;
    double* node_is_conductor = (double*)calloc(num_nodes, sizeof(double));

    /* 要素をループして、導体要素(prop=1,2,3,4)に含まれるノードにフラグを立てる */
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop == 2){
            /* 空気：まだ導体(1)と判定されていない場合のみ 空気(2) をセット */
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                //if(node_is_conductor[gn] != 1){ 
                    node_is_conductor[gn] = 2;
                //}
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 1){
            /* 導体：優先度高 (1) をセット */
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 1; 
            }
        }
    }

    FILE* fp;
    const char* filename;
    filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", "result_elem_0.00001.vtk");
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
	BB_vtk_write_point_vals_scalar(fp, node_is_conductor, fe->total_num_nodes, "phi");

	fclose(fp);


    /* --------------------------------------------------------
       3. phi (Node) に対する境界条件
       
       Case A: 空気領域 (絶縁体) -> phi = 0 で固定 (特異性回避)
       Case B: 導体領域の境界 -> 接地などが必要なら設定
       -------------------------------------------------------- */

    int count = 0;
    if(monolis_mpi_get_global_my_rank()==2){
        for (int i = 0; i < num_nodes; ++i){
            if (node_is_conductor[i] == 1) {
                count++;
                if(count == 1){
                    monolis_set_Dirichlet_bc_C(
                    monolis, 
                    monolis->mat.C.B, 
                    i, 
                    0, 
                    0.0 + 0.0*I
                );

                printf("\n\n\n\n\nadd D_bc\n\n\n\n\n");

                }
                else{
                }
            }
        }
    }

    for (int i = 0; i < num_nodes; ++i){
        if (node_is_conductor[i] == 2) {
            monolis_set_Dirichlet_bc_C(
                monolis, 
                monolis->mat.C.B, 
                i, 
                0, 
                0.0 + 0.0*I
            );
        }
      
    }


}

void set_element_mat_nedelec_Aphi_team6(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    (void)bc;

    const int np = basis->num_integ_points;

    double* J_ip = BB_std_calloc_1d_double(J_ip, np);
    double _Complex* val_ip_C = BB_std_calloc_1d_double_C(val_ip_C, np);

    for(int e=0; e<fe->total_num_elems; ++e){
        const int prop = ned->elem_prop[e];
        double sigma;
        
        if(prop == 1){ // 導体
            sigma = Sigma_cop;
        }
        else if(prop == 2){ // 空気
            //sigma = 0.00001;
            sigma = 0.0;
        }
        else {
            exit(1);
            continue; // その他の要素
        }

        if(prop != 1 && prop != 2){
            fprintf(stderr, "TEAM6: unsupported elem_prop=%d at elem %d\n", prop, e);
            exit(EXIT_FAILURE);
        }

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        /* --------------------
         * (1) A-A : K = ∫ nu curlNi·curlNj
         * sign: si*sj
         * -------------------- */
        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j=0; j<ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += Nu * BBFE_elemmat_mag_mat_curl(
                        ned->curl_N_edge[e][p][i],
                        ned->curl_N_edge[e][p][j],
                        1.0
                    );
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)(si*sj);  // ★ edge_sign

                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        /* 空気(prop=2, sigma=0)なら、以下のσ項は不要 */
        if(sigma <= 0.0){
            continue;
        }

        /* --------------------
         * (2) A-A : i w M = i w ∫ sigma Ni·Nj
         * sign: si*sj
         * -------------------- */
        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j=0; j<ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        ned->N_edge[e][p][j],
                        sigma
                    ) * Omega * I;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)(si*sj);  // ★ edge_sign

                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

        /* --------------------
         * (3) phi-phi : S = ∫ sigma gradψm·gradψn   (★ iω無し)
         * -------------------- */
        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm =  fe->conn[e][m];

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn =  fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[m],
                        fe->geo[e][p].grad_N[n],
                        sigma
                    );
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
            }
        }

        /* --------------------
         * (4) A-phi : C = ∫ sigma N_edge · gradψn   (★ iω無し)
         * row=edge, col=node
         * sign: (row edge sign) = sj
         * -------------------- */
        for(int j=0; j<ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int n=0; n<fe->local_num_nodes; ++n){
                const int gn =  fe->conn[e][n];

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] -= BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[n],
                        ned->N_edge[e][p][j],
                        sigma
                    )* I*Omega;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)sj;  // ★ edge_sign（row edge）

                monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
            }
        }

        /* --------------------
         * (5) phi-A : i w C^T = i w ∫ sigma gradψm · N_edge
         * row=node, col=edge
         * sign: (col edge sign) = si
         * -------------------- */
        for(int m=0; m<fe->local_num_nodes; ++m){
            const int gm =  fe->conn[e][m];

            for(int i=0; i<ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int p=0; p<np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] -= BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        fe->geo[e][p].grad_N[m],
                        sigma
                    ) * Omega * I;
                }

                double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                v *= (double)si;  // ★ edge_sign（col edge）

                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gi, 0, 0, v);
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
    BB_std_free_1d_double_C(val_ip_C, np);
}

/* ============================================================
 * TEAM6 : Complex RHS assembly
 *
 * TEAM6 は通常「外部一様磁場を境界 Dirichlet で与える」ので、
 * 体積 RHS は 0（Js無し）でOKです。
 *
 * もし体積電流 Js を入れるなら、edge_sign は “行側”に掛けます：
 *   rhs(edge_i) += si * ∫ N_i · Js dV
 * ============================================================ */
void set_element_vec_nedelec_Aphi_team6(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    (void)basis; (void)bc; (void)ned;

    /* TEAM6の標準：体積RHSなし（境界Dirichletで励磁） */
    /* ここでは何もしない（ゼロのまま） */
    (void)monolis;
    (void)fe;
}

void copy_Aphi_to_V_phi_C(
    BBFE_DATA* fe,
    NEDELEC* ned,
    double _Complex * Aphi,
    double * V,
    double * phi,
    const int total_num_elems)
{

    for(int i = 0; i < total_num_elems; i++){
        for(int j = 0; j < ned->local_num_edges; j++){
            V[i] = creal(Aphi[ned->nedelec_conn[i][j]]);
        }
    }

    for(int i = 0; i < total_num_elems; i++){
        for(int j = 0; j < fe->local_num_nodes; j++){
            phi[i] = creal(Aphi[fe->conn[i][j]]);
        }
    }

}

void solver_fom_NR_Aphi_team6(
    FE_SYSTEM sys,
    double t,
    int step,
    double* x_prev,
    double* x_curr,
    int n_dof_total)
{
    //double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));
    double _Complex * Aphi = (double _Complex *)calloc(sys.fe.total_num_nodes, sizeof(double _Complex));
    double* A   = BB_std_calloc_1d_double(A,   sys.fe.total_num_nodes);
    double* phi = BB_std_calloc_1d_double(phi, sys.fe.total_num_nodes);

    set_element_mat_nedelec_Aphi_team6(&(sys.monolis), &(sys.fe), &(sys.basis),
                            &(sys.bc), &(sys.ned));
    set_element_vec_nedelec_Aphi_team6(&(sys.monolis), &(sys.fe), &(sys.basis),
                            &(sys.bc), &(sys.ned));

    double val = Bz_of_time(t);
    apply_nedelec_boundary_conditions_team6(&(sys.monolis), &(sys.fe), &(sys.bc), &(sys.ned), val, t, x_curr, sys.cond.directory);

    monowrap_solve_C(
        &(sys.monolis),
        &(sys.monolis_com),
        Aphi,
        MONOLIS_ITER_COCG,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

    for(int i = 0; i < sys.fe.total_num_elems; i++){
        for(int j = 0; j < sys.fe.local_num_nodes; j++){
            //Aphi[sys.fe.conn[i][j]] = I * Omega * Aphi[sys.fe.conn[i][j]];
        }
    }

    for(int i = 0; i < sys.fe.total_num_nodes; i++){
        x_prev[i] = creal(Aphi[i]);
        x_curr[i] = cimag(Aphi[i]);
    }

    //copy_Aphi_to_V_phi_C(&(sys.fe), &(sys.ned), Aphi, A, phi, sys.fe.total_num_elems);
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
    /*
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
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            1,
            &(sys.rom_sups.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm_v.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm_v.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            3,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
	ROM_std_hlpod_pre(
            &(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

	ROM_std_hlpod_pre(
            &(sys.rom_sups),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
*/
    /******************/

    /*for offline******/
    //ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);

	/****************** solver ********************/
    int step = 0;
    double t = 0;
	int file_num = 0;
    int count = 0;  //for ROM
    double t_hs = 0.0; //for hot start
    int step_hs = 0; //for hot start

    int nsteps = (int)ceil(sys.vals.finish_time / sys.vals.dt);

    
    //output_result_B_node_analysis(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, t, sys.cond.directory);       
    //output_team6_analytic_series(&(sys.fe), sys.cond.directory); 
    output_team6_benchmark(sys.cond.directory); 
/*
    char fname[BUFFER_SIZE];         
    snprintf(fname, BUFFER_SIZE, "hot_start/%s.%d.dat", "velosity_pressure", monolis_mpi_get_global_my_rank());
    t_hs = hot_start_read_initialize_val(sys.vals.Aphi_time, fname, sys.cond.directory);
    step_hs = t_hs / sys.vals.dt;
    printf("Hot start time: %lf\n", t_hs);
    printf("Hot start step: %d\n", step_hs);
    t = t_hs;
    printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

    output_result_B_node_analysis(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, t, sys.cond.directory);
*/

/*
	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hs),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			1);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            ((double)sys.vals.finish_time - t_hs),
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			1);
*/

    for (step = step_hs; step <= nsteps; ++step) {
        t += sys.vals.dt;

        printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
 /*
        solver_fom_NR_Aphi_team11(
            sys, t, count, 
            sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr,
            sys.fe.total_num_nodes);
*/     
        solver_fom_NR_Aphi_team6(
            sys, t, count, 
            sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr,
            sys.fe.total_num_nodes);

/*
        solver_fom_NR_Aphi_collect_snapmat(
            sys, t, count, 
            sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr,
            sys.fe.total_num_nodes);
*/
        log_accuracy_metrics(
            &sys, sys.vals.Aphi_time,
            sys.vals.Aphi_time_curr, step, t, sys.vals.dt, CURRENT_AMP);
        
        for(int i=0; i<sys.fe.total_num_nodes; ++i){
            //sys.vals.Aphi_time[i] = sys.vals.Aphi_time_curr[i];
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
             
            output_result_B_node_sim(&(sys.fe), &(sys.vals), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, t, sys.cond.directory);  
            output_result_B_node_sim_I(&(sys.fe), &(sys.vals), &(sys.basis), &(sys.ned), sys.vals.Aphi_time_curr, t, sys.cond.directory);            
            //output_files_nedelec(&sys, step, t);
            //output_files(&sys, step, t);

            //output_result_file_B_integration_point_sim_axis(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time, t, sys.cond.directory);  
            //output_result_file_B_integration_point_sim_axis_I(&(sys.fe), &(sys.basis), &(sys.ned), sys.vals.Aphi_time_curr, t, sys.cond.directory);  
        }

        if(step%1 == 0){
                char fname[BUFFER_SIZE];
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%d.dat", "velosity_pressure", monolis_mpi_get_global_my_rank());
                hot_start_write_initialize_val(sys.vals.Aphi_time, sys.fe.total_num_nodes, 1, t, fname, sys.cond.directory);
            }
        else{
        }
    }

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
