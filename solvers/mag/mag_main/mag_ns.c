
#include "./../mag_core/convdiff_core.h"
#include "./../mag_core/nedelec_core.h"

#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include <complex.h>

#include "mag_dataset.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 4;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-10;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;

const double DELTA    = 1.0E-06;
const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";


const double Nu    = 1.0;                  // 無次元化 or μ=1
const double Sigma = 1.0;                  // 無次元化（典型導電率）
const double freq  = 10.0;                 // [Hz]
const double Omega = 1.0;   // 角周波数

double manusol_get_mass_coef(
		double x[3])
{
	double val = 1.0;

	return val;
}


double manusol_get_diff_coef(
		double x[3])
{
    double val = 1.0;

	return val;
}

void set_ansvec(
    BBFE_DATA* fe,
    double** source)
{
    
    for(int i = 0; i < fe->total_num_nodes; i++) {
        double X = fe->x[i][0];
        double Y = fe->x[i][1];
        double Z = fe->x[i][2];  // ノードのZ座標

        source[i][0] = 0.0;
        source[i][1] = 0.0;
        source[i][2] = sin(M_PI * X) * sin(M_PI * Y);

    }
}

void manusol_set_theo_sol(
		BBFE_DATA* fe,
		double*  theo_sol,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		//theo_sol[i] = cos(M_PI*fe->x[i][0]) * cos(M_PI*fe->x[i][1]) * cos(M_PI*fe->x[i][2]);
        theo_sol[i] = fe->x[i][0]*fe->x[i][0] - fe->x[i][1]*fe->x[i][1];
	}
}


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
    double val = x*x - y*y;

	return val;
}

void set_ansvec_edge(
    NEDELEC*    ned,
    double** source,
    const int   total_num_edge)
{
    
    for(int i = 0; i < total_num_edge; i++) {
        double X = ned->nedelec_coords[i][0];  // ノードのX座標
        double Y = ned->nedelec_coords[i][1];  // ノードのY座標
        double Z = ned->nedelec_coords[i][2];  // ノードのZ座標

        source[i][0] = 0.0;
        source[i][1] = 0.0;
        source[i][2] = sin(M_PI * X) * sin(M_PI * Y);

    }
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
            c[j] = u_coeff[ge_j];
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

	// for manufactured solution
    manusol_set_theo_sol(fe, vals->theo_sol, t);
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->theo_sol, vals->phi);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "abs_error");
	BB_vtk_write_point_vals_scalar(fp, vals->theo_sol, fe->total_num_nodes, "theoretical");

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

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.phi,
			filename,
			sys->cond.directory);

	/**** for manufactured solution ****/
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_theo_sol(&(sys->fe), source, t);
    
	double L2_error = BBFE_elemmat_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis_com),
			t,
			sys->vals.phi,
			manusol_get_sol);

	printf("%s L2 error: %e\n", CODENAME, L2_error);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = BBFE_sys_write_add_fopen(fp, "l2_error.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error);
		fclose(fp);
	}

	//BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}

void output_result_file_vtk_nedelec_edge(
		BBFE_DATA*      fe,
		VALUES*         vals,
        NEDELEC*        ned,
		const char*     filename,
		const char*     directory,
		double          t)
{

    double** node_result = BB_std_calloc_2d_double(node_result, fe->total_num_nodes, 3);
    double** edge_result = BB_std_calloc_2d_double(edge_result, fe->total_num_nodes, 3);

    reconstruct_field_at_edge_points(fe, ned, vals->V, edge_result, fe->total_num_elems);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    BB_vtk_write_points_3d(fp, fe->total_num_nodes, ned->nedelec_coords);
	BB_vtk_write_cells(fp, fe->total_num_elems, 12, ned->nedelec_conn);
    BB_vtk_write_cell_types(fp, fe->total_num_elems, 4);

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    BB_vtk_write_point_vals_scalar(fp, vals->V, fe->total_num_nodes, "coordinates");
    BB_vtk_write_point_vals_vector(fp, edge_result, fe->total_num_nodes, "V_edge");

	double** source;
	source = BB_std_calloc_2d_double(source, fe->total_num_nodes, 3);
    set_ansvec_edge(ned, source, fe->total_num_nodes);
	BB_vtk_write_point_vals_vector(fp, source, fe->total_num_nodes, "theo_sol_edge");

    double** edge_error;
	edge_error = BB_std_calloc_2d_double(edge_error, fe->total_num_nodes, 3);
    double error = calc_nodal_error_vector(edge_result, source, edge_error, fe->total_num_nodes);
    printf("\n\nL2 error: %e\n\n", error);
    BB_vtk_write_point_vals_vector(fp, edge_error, fe->total_num_nodes, "edge_error");
    BB_std_free_2d_double(source, fe->total_num_nodes, 3);
    BB_std_free_2d_double(edge_error, fe->total_num_nodes, 3);
    
	fclose(fp);
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
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", "result_nedelec_edge.vtk");
	output_result_file_vtk_nedelec_edge(
			&(sys->fe),
			&(sys->vals),
            &(sys->ned),
			filename,
			sys->cond.directory,
			t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.V,
			filename,
			sys->cond.directory);

}


void apply_nedelec_boundary_conditions(
        MONOLIS* monolis,
        BBFE_DATA* fe,
        BBFE_BC*      bc,
        NEDELEC*     ned)
{
    double val_points[2];
    for (int e = 0; e < fe->total_num_elems; e++) {

        // 各エッジごとに、両端の節点が境界上にあるかをチェック
        for (int ed = 0; ed < ned->local_num_edges; ed++) {
            int local_node1 = hex_edge_conn[ed][0];
            int local_node2 = hex_edge_conn[ed][1];

            int global_node1 = fe->conn[e][local_node1];
            int global_node2 = fe->conn[e][local_node2];

            if (bc->D_bc_exists[global_node1] && bc->D_bc_exists[global_node2]) {
                int global_edge = ned->nedelec_conn[e][ed];

                monolis_set_Dirichlet_bc_C(
                        monolis,
                        monolis->mat.C.B,
                        global_edge,
                        0,
                        0.0);
                    
            }
        }
    }


    for (int e = 0; e < fe->total_num_elems; e++) {

        for (int ed = 0; ed < fe->local_num_nodes; ed++) {

            int global_node1 = fe->conn[e][ed];

            //case 1
            double _Complex val = fe->x[global_node1][0]*fe->x[global_node1][0] - fe->x[global_node1][1]* fe->x[global_node1][1];
            //case 2
            //double _Complex val = 0.0;
            
            if (bc->D_bc_exists[global_node1]) { 
                    monolis_set_Dirichlet_bc_C(
                            monolis,
                            monolis->mat.C.B,
                            global_node1,
                            0,
                            val);

            }
        }
    }

}


void set_element_mat_nedelec_Aphi(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
        BBFE_BC*      bc,
        NEDELEC*     ned,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	double _Complex* val_ip_C;
	val_ip_C      = BB_std_calloc_1d_double_C(val_ip_C, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

    double*** N_edge = BB_std_calloc_3d_double(N_edge, np, ned->local_num_edges, 3);
	double*** curl_N_edge = BB_std_calloc_3d_double(curl_N_edge, np, ned->local_num_edges, 3);

    for (int e = 0; e < fe->total_num_elems; e++) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for (int p = 0; p < np; p++) {
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]); 
            a_ip[p] = manusol_get_mass_coef(x_ip[p]);
            k_ip[p] = manusol_get_diff_coef(x_ip[p]);
        }


        for (int i = 0; i < ned->local_num_edges; i++) {
            for (int j = 0; j < ned->local_num_edges; j++) {
                for (int p = 0; p < np; p++) {
                    val_ip_C[p] = 0.0 + 0.0 * I;

                    val_ip_C[p] += BBFE_elemmat_mag_mat_curl(
                                    ned->curl_N_edge[e][p][i], 
                                    ned->curl_N_edge[e][p][j], 1) * Nu;
                }

                double _Complex integ_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_C(
                    monolis,
                    ned->nedelec_conn[e][i],
                    ned->nedelec_conn[e][j],
                    0, 0,
                    integ_val);

            }
        }

        for (int i = 0; i < ned->local_num_edges; i++) {
            for (int j = 0; j < ned->local_num_edges; j++) {
                for (int p = 0; p < np; p++) {
                    val_ip_C[p] = 0.0 + 0.0 * I;
                    
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                                    ned->N_edge[e][p][i],
                                    ned->N_edge[e][p][j], Sigma) * Omega * I;
                }

                double _Complex integ_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_C(
                    monolis,
                    ned->nedelec_conn[e][i],
                    ned->nedelec_conn[e][j],
                    0, 0,
                    integ_val);

            }
        }

        for (int i = 0; i < fe->local_num_nodes; i++) {
            for (int j = 0; j < fe->local_num_nodes; j++) {
                for (int p = 0; p < np; p++) {
                    val_ip_C[p] = 0.0 + 0.0 * I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i],
                            fe->geo[e][p].grad_N[j], Sigma) * Omega * I;
                }

                double _Complex integ_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_C(
                    monolis,
                    fe->conn[e][i],
                    fe->conn[e][j],
                    0, 0,
                    integ_val);

            }
        }

        for (int i = 0; i < fe->local_num_nodes; i++) {
            for (int j = 0; j < ned->local_num_edges; j++) {
                for (int p = 0; p < np; p++) {
                    val_ip_C[p] = 0.0 + 0.0 * I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i],
                            ned->N_edge[e][p][j], Sigma) * Omega * I;
                }

                double _Complex integ_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_C(
                    monolis,
                    fe->conn[e][i],
                    ned->nedelec_conn[e][j],
                    0, 0,
                    integ_val);
            }
        }


        for (int i = 0; i < ned->local_num_edges; i++) {
            for (int j = 0; j < fe->local_num_nodes; j++) {
                for (int p = 0; p < np; p++) {
                    val_ip_C[p] = 0.0 + 0.0 * I;
                    val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            fe->geo[e][p].grad_N[j], Sigma) * Omega * I;
                }

                double _Complex integ_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_C(
                    monolis,
                    ned->nedelec_conn[e][i],
                    fe->conn[e][j], 
                    0, 0,
                    integ_val);

            }
        }

    }


	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}


void set_element_vec_nedelec_Aphi(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
        BBFE_BC*      bc,
        NEDELEC*     ned,
		VALUES*      vals,
		double       t)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
	double _Complex* val_ip_C;
	val_ip_C      = BB_std_calloc_1d_double_C(val_ip_C, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, 3);

	double _Complex* f_ip_C;
	f_ip_C      = BB_std_calloc_1d_double_C(f_ip_C, np);

    for (int e = 0; e < fe->total_num_elems; e++) {
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
        BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

        for (int p = 0; p < np; p++) {
            BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]); 
        }

        for (int i = 0; i < ned->local_num_edges; i++) {
            for (int p = 0; p < np; p++) {
                compute_Js(x_ip[p][0], x_ip[p][1], x_ip[p][2], Nu, Sigma, Omega, f_ip_C);            
                
                val_ip_C[p] = ned->N_edge[e][p][i][0] * f_ip_C[0]
                        + ned->N_edge[e][p][i][1] * f_ip_C[1]
                        + ned->N_edge[e][p][i][2] * f_ip_C[2];
    
            }
            double _Complex rhs_val = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, Jacobian_ip);
            
            monolis->mat.C.B[ned->nedelec_conn[e][i]] += rhs_val;
        }
    }


	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
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
    
	BBFE_mag_pre(
			&(sys.fe), &(sys.basis), &(sys.ned), (&sys.bc), (&sys.monolis0), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

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

    compute_nedelec_edge_coords(
            &(sys.fe),
            &(sys.ned),
           sys.fe.total_num_elems);

    const char* filename;
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			BLOCK_SIZE);

    sys.vals.Aphi = (double _Complex *)calloc(sys.fe.total_num_nodes, sizeof(double _Complex));

    printf("set_mat\n\n");
	set_element_mat_nedelec_Aphi(
			&(sys.monolis0),
			&(sys.fe),
			&(sys.basis),
            &(sys.bc),
            &(sys.ned),
			&(sys.vals));

	/****************** solver ********************/
	int file_num = 0;

    set_element_vec_nedelec_Aphi(
            &(sys.monolis0),
            &(sys.fe),
            &(sys.basis),
            &(sys.bc),
            &(sys.ned),
            &(sys.vals),
            0.0);

    apply_nedelec_boundary_conditions(
            &(sys.monolis0),
            &(sys.fe),
            &(sys.bc),
            &(sys.ned));

    monolis_show_iterlog (&(sys.monolis0), true);
    
    monowrap_solve_C(
            &(sys.monolis0),
            &(sys.monolis_com),
            sys.vals.Aphi,
            MONOLIS_ITER_COCG,
            MONOLIS_PREC_DIAG,
            sys.vals.mat_max_iter,
            sys.vals.mat_epsilon);

    /**********************************************/

    sys.vals.V = BB_std_calloc_1d_double(sys.vals.V, sys.fe.total_num_nodes);
    sys.vals.phi = BB_std_calloc_1d_double(sys.vals.phi, sys.fe.total_num_nodes);

    copy_Aphi_to_V_phi(
            sys.vals.Aphi,
            sys.vals.V,
            sys.vals.phi,
            sys.fe.total_num_nodes,
            sys.fe.total_num_nodes);

    //output_files_nedelec(&sys, file_num, 0.0);
    output_files(&sys, file_num, 0.0);

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis0));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
