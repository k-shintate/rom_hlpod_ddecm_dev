#pragma once

#include "convdiff_core.h"

const int hex_edge_conn[12][2] = {
    {0, 1},
    {3, 2},
    {4, 5},
    {7, 6},
    {0, 3},
    {1, 2},
    {4, 7},
    {5, 6},
    {0, 4},
    {1, 5},
    {3, 7},
    {2, 6}
};

static const int tet_edge_conn[6][2] = {
    {0, 1}, {1, 2}, {2, 0},
    {0, 3}, {1, 3}, {2, 3}
};

/* ====== 候補エッジ構造体 ====== */
typedef struct {
    int a, b;      /* 無向キー: a=min(u,v), b=max(u,v) */
    int u, v;      /* 要素ローカルの向き付き端点 (u->v) */
    int elem;      /* 要素ID */
    int ledge;     /* ローカルエッジID */
} EdgeCand;


typedef struct
{
	int**    nedelec_conn;
    int**    edge_sign; 
    double**  nedelec_coords;
    double**** N_edge;
    double**** curl_N_edge;
    
    int num_edges;

    int*    elem_prop;
    int local_num_edges;

} NEDELEC;



// Aphi -> 要素ごとの B を作る（ガウス点のヤコビアン重み付き平均）
void compute_B_cell_average(
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    NEDELEC*     ned,
    const double* Aphi,      // [total_num_nodes + ned->num_edges]
    double**      B_cell     // [fe->total_num_elems][3] 事前に確保 or NULLなら内部で確保
);

void copy_Aphi_to_V_phi_time2(
    BBFE_DATA* fe,
    NEDELEC* ned,
    double * Aphi,
    double * V,
    double * phi,
    const int total_num_elems);

void build_dirichlet_edge_mask_from_boundary_faces_tet(
    const BBFE_DATA* fe,
    const BBFE_BC* bc,
    const NEDELEC* ned,
    bool* is_dir_edge,          /* size = num_global_edges, zero-initialized */
    int  num_global_edges);


void accumulate_B_cell_to_nodes(
    BBFE_DATA* fe,
    double** B_cell,
    double** B_node_out);

// VTK 出力（節点ベクトル）
void output_B_node_vtk(
    BBFE_DATA* fe, BBFE_BASIS* basis, NEDELEC* ned, const double* Aphi,
    const char* filename, const char* directory);

void compute_nedelec_edge_coords(
        BBFE_DATA* fe,
        NEDELEC* ned,
        int total_elems);


void BBFE_mag_set_basis(
        BBFE_DATA*     fe,
        BBFE_BASIS*    basis,
        NEDELEC*       ned,
        int            local_num_nodes,
        int            num_integ_points_each_axis);


void read_edge_sign(
        BBFE_DATA *fe,
        NEDELEC *ned,
        const char *filename,
        const char *directory,
        int num_integ_points);


void ROM_std_hlpod_set_nonzero_pattern_bcsr(
        MONOLIS*     	monolis,
        const char*     label,
        const char*		directory);


void set_elem_prop(
		BBFE_DATA*      fe,
        NEDELEC*        ned,
		const char*     directory,
        const char*     filename);


void read_elem_types(
        BBFE_DATA*    fe,
        NEDELEC *ned,
        const char*     label,
        const char*		directory);


void read_num_nodes(
        BBFE_DATA*    fe,
        const char*     label,
        const char*		directory);


void read_connectivity_graph_lag_nedelec(
        BBFE_DATA *fe,
        NEDELEC *ned,
        const char *filename,
        const char *directory,
        int num_integ_points);


void monowrap_solve_C(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double _Complex*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon);


void monowrap_solve_R(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon);

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
        BBFE_BC*     bc,
        const char*  filename,
        const char*  directory,
        const int    total_num_nodes,
        const int    block_size);

void BBFE_mag_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
        NEDELEC*      ned,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution);

void set_element_mat_nedelec_Aphi_team21a2(
        MONOLIS*     monolis,
        BBFE_DATA*   fe,
        BBFE_BASIS*  basis,
        BBFE_BC*     bc,
        NEDELEC*     ned);

void set_element_vec_nedelec_Aphi_team21a2(
        MONOLIS*     monolis,
        BBFE_DATA*   fe,
        BBFE_BASIS*  basis,
        BBFE_BC*     bc,
        NEDELEC*     ned);

void ROM_std_hlpod_set_nonzero_pattern_bcsr_C(
        MONOLIS*     	monolis,
        const char*     label,
        const char*		directory);

void BBFE_mag_pre_C(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
        NEDELEC*      ned,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution);