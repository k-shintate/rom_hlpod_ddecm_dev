#pragma once

#include "fluid_core.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

#define T4_WALL_DIAG_NUM_REGIONS 8

enum {
    T4_WALL_DIAG_FRONT_POS_DRAG = 0,
    T4_WALL_DIAG_REAR_NEG_DRAG  = 1,
    T4_WALL_DIAG_TOP_POS_LIFT   = 2,
    T4_WALL_DIAG_BOTTOM_NEG_LIFT= 3,
    T4_WALL_DIAG_SIDE_POS       = 4,
    T4_WALL_DIAG_SIDE_NEG       = 5,
    T4_WALL_DIAG_SLANT_MIXED    = 6,
    T4_WALL_DIAG_OTHER          = 7
};

static const char* t4_wall_diag_region_name[T4_WALL_DIAG_NUM_REGIONS] = {
    "front_pos_drag",
    "rear_neg_drag",
    "top_pos_lift",
    "bottom_neg_lift",
    "side_pos",
    "side_neg",
    "slant_mixed",
    "other"
};


typedef struct {
    double area;
    double h_sum;
    double h_mean;
    double h_min;
    double h_max;
} T4WallPenaltyScaleDiagnostics;

typedef struct {
    double area;
    double int_n[3];

    double max_abs_un_over_Uref;
    double int_abs_un;
    double int_un;
    double int_un2;
    double mean_abs_un_over_Uref;
    double mean_un_over_Uref;
    double rms_un_over_Uref;
    double Uref_scale;

    double p_min;
    double p_max;
    double p_mean;
    double p_rms;
    double p_area_sum;
    double p2_area_sum;

    double cp_min;
    double cp_max;
    double cp_mean;
    double cp_rms;
    double cp_area_sum;
    double cp2_area_sum;

    /* Pressure force in the drag direction: int_Gamma p (n.eU) dGamma.
       If eU=(1,0,0), this is int_Gamma p n_x dGamma. */
    double int_pnU;
    double area_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_pnU_region[T4_WALL_DIAG_NUM_REGIONS];
    double int_abs_un_region[T4_WALL_DIAG_NUM_REGIONS];

    double max_norm_r_over_Uref;
    double int_norm_r;
    double mean_norm_r_over_Uref;

} T4WallDiagnostics;

/* ========================================================================== */
/* TET4/TRI3 wall mesh and y+ diagnostics                                     */
/* ========================================================================== */

#define T4_WALL_MESH_DIAG_NUM_STATS 18
#define T4_WALL_MESH_DIAG_NUM_YPLUS_BINS 5

enum {
    T4_WMD_YPLUS_GRAD_CENTROID = 0,
    T4_WMD_YPLUS_SPALDING_CENTROID,
    T4_WMD_YPLUS_SPALDING_OVER_GRAD,
    T4_WMD_UTAU,
    T4_WMD_TAU_W,
    T4_WMD_CF,
    T4_WMD_H_N_PENALTY,
    T4_WMD_OWNER_FACE_HEIGHT,
    T4_WMD_H_T,
    T4_WMD_ASPECT_RATIO,
    T4_WMD_NONORTHOGONALITY_DEG,
    T4_WMD_OWNER_TET_QUALITY,
    T4_WMD_MU_EFF_OVER_MU,
    T4_WMD_BETA_DT_OVER_BETA_MU,
    T4_WMD_PENALTY_TO_SHEAR,
    T4_WMD_OWNER_CFL,
    T4_WMD_OWNER_RE_H,
    T4_WMD_OWNER_PE_EFF
};

typedef struct {
    double force_row_sum[3];
    double force_reference[3];

    double force_pressure[3];
    double force_weak_viscous_eff[3];
    double force_penalty_file[3];

    double max_qp_row_reference_abs_error;
    double max_partition_unity_abs_error;
    double max_gradient_sum_abs_error;

    /* doubles so monolis_allreduce_R can reduce them directly */
    double mapped_surface_faces;
    double skipped_surface_faces;
    double quadrature_points;
    double invalid_faces;
} T4NitscheRowSumDiagnostics;



typedef struct {
    /*
     * Raw area-weighted accumulators.
     * MPI reduction後も保持する。
     */
    double area;

    double sum[T4_WALL_MESH_DIAG_NUM_STATS];
    double sum2[T4_WALL_MESH_DIAG_NUM_STATS];
    double min[T4_WALL_MESH_DIAG_NUM_STATS];
    double max[T4_WALL_MESH_DIAG_NUM_STATS];

    /*
     * Finalized statistics.
     */
    double mean[T4_WALL_MESH_DIAG_NUM_STATS];
    double rms[T4_WALL_MESH_DIAG_NUM_STATS];

    /*
     * y+ area bins:
     *   0: y+ < 1
     *   1: 1 <= y+ < 5
     *   2: 5 <= y+ < 30
     *   3: 30 <= y+ < 100
     *   4: y+ >= 100
     */
    double yplus_bin_area[T4_WALL_MESH_DIAG_NUM_YPLUS_BINS];
    double yplus_bin_fraction[T4_WALL_MESH_DIAG_NUM_YPLUS_BINS];

    /*
     * Separation / reverse shear.
     */
    double streamwise_shear_valid_area;
    double reverse_shear_area;
    double reverse_shear_area_fraction;

    /*
     * Mesh-quality area fractions.
     */
    double nonorth_gt_30_area;
    double nonorth_gt_60_area;
    double nonorth_gt_30_fraction;
    double nonorth_gt_60_fraction;

    double quality_lt_0p1_area;
    double quality_lt_0p2_area;
    double quality_lt_0p1_fraction;
    double quality_lt_0p2_fraction;

    /*
     * VMS effective-viscosity area fractions.
     */
    double mu_ratio_gt_2_area;
    double mu_ratio_gt_10_area;
    double mu_ratio_gt_2_fraction;
    double mu_ratio_gt_10_fraction;

    /*
     * Nitsche time-penalty dominance.
     */
    double beta_dt_dominant_area;
    double beta_dt_dominant_fraction;

    /*
     * Invalid or degenerate faces.
     */
    double invalid_face_count;
} T4WallMeshDiagnostics;


/* ========================================================================== */
/* Domain continuity diagnostics                                              */
/* ========================================================================== */

typedef struct {
    double volume;

    double int_abs_div_u;
    double int_div_u2;
    double max_abs_div_u;

    double mean_abs_div_u;
    double rms_div_u;

    double invalid_element_count;
} T4DomainContinuityDiagnostics;


/* ========================================================================== */
/* Public functions                                                           */
/* ========================================================================== */
/*
void BBFE_fluid_sups_impose_surface_velocity_on_values(
    const BBFE_DATA* surf,
    double** v,
    const double imposed_velocity[3]);

int BBFE_fluid_sups_add_surface_velocity_Dirichlet(
    BBFE_BC* bc,
    const BBFE_DATA* surf,
    const double imposed_velocity[3]);
*/
int calc_wall_mesh_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    VALUES* vals,
    double rho,
    double mu_molecular,
    double Uref,
    const double Uw[3],
    const double eU[3],
    double gamma_n,
    int write_nodal_yplus,
    T4WallMeshDiagnostics* diag);

void wall_mesh_diagnostics_allreduce(
    T4WallMeshDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void output_wall_mesh_diagnostics_tet4_tri3(
    const T4WallMeshDiagnostics* diag,
    double t,
    const char* directory);


int calc_domain_continuity_diagnostics_tet4(
    const BBFE_DATA* fe,
    const VALUES* vals,
    T4DomainContinuityDiagnostics* diag);

void domain_continuity_diagnostics_allreduce(
    T4DomainContinuityDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void output_domain_continuity_diagnostics_tet4(
    const T4DomainContinuityDiagnostics* diag,
    double t,
    const char* directory);

void output_cavity_center_vx(
    VALUES*        vals,
    const char*    method,
    const char*    directory);

void initialize_velocity_pressure_ahmedbody(
        BBFE_DATA*     fe,
        double** v,
        double* p,
        const int total_num_nodes);

int BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
    BBFE_BC* bc,
    const BBFE_DATA* fe,
    const double z0,
    const double tolerance);

void initialize_velocity_pressure_karman_vortex(
        double** v,
        double* p,
        const int total_num_nodes);

void initialize_velocity_pressure_backstep(
    BBFE_DATA*     fe,
        double**    v,
        double*     p,
        const int   total_num_nodes);

void initialize_velocity_pressure_cavity(
        double** v,
    double** v_old,
        double* p,
        const int total_num_nodes);

void BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void output_result_file_karman_vortex(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure_inf(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void memory_allocation_nodal_values(
    VALUES*         vals,
    const int       total_num_nodes);

void assign_default_values(
    VALUES*         vals);

void print_all_values(
    VALUES*         vals);

void read_calc_conditions(
    VALUES*         vals,
    const char*     directory);

void output_result_file_vtk(
    BBFE_DATA*     fe,
    VALUES*        vals,
    const char*    filename,
    const char*    directory,
    double         t);

void output_files(
    FE_SYSTEM*  sys,
    int         file_num,
    double      t);

void set_element_mat_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_corr(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_ppe(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void BBFE_fluid_sups_read_Dirichlet_bc(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
    BBFE_BC*     bc,
    VALUES*      vals,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void set_element_mat(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void solver_fom_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void memory_allocation_nodal_values_AB2(
    VALUES*         vals,
    const int       total_num_nodes);

void ROM_BB_vec_copy_2d(
    double**  in,   //input
    double**  out,  //output
    const int num1,
    const int num2);

void solver_fom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void pre_surface(
    BBFE_DATA*    surf,
    BBFE_BASIS*   basis,
    const char*   directory,
    const char*   fname,
    int           num_integ_points_each_axis);

void pre_surface_internal(
                //MONOLIS_COM*  monolis_com_surf,
                BBFE_DATA*    surf,
                BBFE_BASIS*   basis,
                const char*   directory,
                const char*   fname,
                int           num_integ_points_each_axis);

void calc_Cd_p(
    BBFE_DATA*  surf,
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    VALUES*     vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out);

int calc_Cd_v(const BBFE_DATA* surf, const BBFE_DATA* fe, const BBFE_BASIS* basis, MONOLIS_COM* monolis_com,
                       const VALUES* vals,
                       double rho, double mu, double Uinf, double Aref,
                       const double eU_in[3], const double eP_in[3],
                       double* D_out, double* L_out, double* Cd_out, double* Cl_out);

void update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes);

double calc_internal_norm_2d(
    double**  in,       //input
        const int num1,
    const int num2);

    double calc_internal_norm_1d(
        double*  in,       //input
        const int num1,
        const int num2);

void set_element_mat_NR_Tezuer(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

double BBFE_elemmat_fluid_sups_coef_metric_tensor_LSIC(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

void solver_fom_NR(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
                BBFE_BC*     bc,
                const char*  filename,
                const char*  directory,
                const int    total_num_nodes,
                const int    block_size);

void solver_fom_NR_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

double BBFE_elemmat_fluid_sups_coef_metric_tensor(
    const double J_inv[3][3],
    const double Jacobian,
    const double density,
    const double viscosity,
    const double v[3],
    const double dt);

void ROM_ddecm_set_reduced_mat_para_save_memory(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);


void ROM_ddecm_set_reduced_vec_para(
    MONOLIS*        monolis,
    BBFE_DATA*      fe,
    VALUES*         vals,
    BBFE_BASIS*     basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int       num_modes,
    const int       num_subdomains,
    const double    dt,
    double          t);

void ROM_set_reduced_vec_para_debug(
        MONOLIS*        monolis,
        BBFE_DATA*      fe,
        VALUES*         vals,
        BBFE_BASIS*     basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
        HLPOD_MAT*      hlpod_mat,
        const int       num_modes,
        const int       num_subdomains,
        const double    dt,
        double          t);

void ROM_set_reduced_mat_para_save_memory_debug(
        MONOLIS*        monolis,
        BBFE_DATA*      fe,
        VALUES*         vals,
        BBFE_BASIS*     basis,
        BBFE_BC*        bc,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int       num_modes,
        const int       num_subdomains,
        const double    dt);

void HROM_ddecm_set_residuals_NR_Tezuer(
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals,
    BBFE_BC*        bc,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int               num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_reduced_mat_para_save_memory_NR(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);

void HROM_ddecm_set_reduced_vec_para_NR(
        MONOLIS*         monolis,
        BBFE_DATA*       fe,
        VALUES*          vals,
        BBFE_BASIS*      basis,
        HR_VALUES*       hr_vals,
        HLPOD_VALUES*    hlpod_vals,
        HLPOD_DDHR*      hlpod_ddhr,
        HLPOD_MAT*       hlpod_mat,
        const int        num_modes,
        const int        num_subdomains,
        const double     dt,
        double           t);


void HROM_ddecm_set_reduced_mat_para_save_memory_NR2(
    MONOLIS*      monolis,
    BBFE_DATA*    fe,
    VALUES*       vals,
    BBFE_BASIS*   basis,
    BBFE_BC*      bc,
    HLPOD_VALUES* hlpod_vals,
    HLPOD_MAT*    hlpod_mat,
    HLPOD_DDHR*   hlpod_ddhr,
    const int     num_modes,
    const int     num_subdomains,
    const double  dt);

void HROM_ddecm_set_reduced_vec_para_NR2(
        MONOLIS*         monolis,
        BBFE_DATA*       fe,
        VALUES*          vals,
        BBFE_BASIS*      basis,
        HR_VALUES*       hr_vals,
        HLPOD_VALUES*    hlpod_vals,
        HLPOD_DDHR*      hlpod_ddhr,
        HLPOD_MAT*       hlpod_mat,
        const int        num_modes,
        const int        num_subdomains,
        const double     dt,
        double           t);



void BBFE_fluid_initialize_velocity_pressure(
        double**    v,
        double*     p,
        const int   total_num_nodes);

void BBFE_fluid_initialize_velocity_pressure_karman_vortex(
	double**    v,
    double**    v_old,
	double*     p,
	const int   total_num_nodes);

void BBFE_fluid_sups_add_velocity_pressure(
        double**  v,
        double*   p,
        double*   sol_vec,
        const int total_num_nodes);

void BBFE_fluid_add_Dbc(
        double*         ansvec,
        BBFE_BC*        bc,
        const int       total_num_nodes,
        const int       dof);

void BBFE_fluid_update_velocity_pressure_NR(
        double**    v,
        double**    delta_v,
        double*     p,
        double*     delta_p,
        const int   total_num_nodes);

double BBFE_fluid_calc_internal_norm_2d(
        double**  in,
        const int num1,
        const int num2);

double BBFE_fluid_calc_internal_norm_1d(
        double*  in,
        const int num1,
        const int num2);

void BBFE_fluid_vec_copy_2d(
        double**  in,
        double**  out,
        const int num1,
        const int num2);

void BBFE_fluid_vec_copy_1d(
        double* dst,
        const double* src,
        int size);

double BBFE_fluid_calc_internal_norm_total(
        const double* rvec,
        int n_internal_vertex);

double BBFE_fluid_calc_internal_norm_momentum(
        const double* rvec,
        int n_internal_vertex);

double BBFE_fluid_calc_internal_norm_mass(
        const double* rvec,
        int n_internal_vertex);


void set_wall_face_vecmat_symmetric_nitsche(
    MONOLIS* monolis,
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_eff,
    const double Uw[3],
    double gamma_n);

int calc_Cd_v_tet4_tri3_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out);

void calc_Cd_p_hex_or_tet(
    BBFE_DATA*  surf,
    BBFE_DATA*  fe,
    BBFE_BASIS* basis,
    VALUES*     vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out);

int calc_wall_nitsche_gamma_sweep_split_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma);


void wall_gamma_sweep_split_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    MONOLIS_COM* monolis_com);


void output_wall_gamma_sweep_split_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_by_gamma,
    const double* L_penalty_by_gamma,
    const double* D_tangent_by_gamma,
    const double* L_tangent_by_gamma,
    const double* D_total_by_gamma,
    const double* L_total_by_gamma,
    double t,
    const char* directory);

int calc_wall_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    T4WallDiagnostics* diag);

void wall_diagnostics_allreduce(
    T4WallDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void output_wall_diagnostics_tet4_tri3(
    const T4WallDiagnostics* diag,
    double t,
    const char* directory);


static void t4_wall_penalty_scale_finalize(
    T4WallPenaltyScaleDiagnostics* diag);

int calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    T4WallPenaltyScaleDiagnostics* scale_diag);

int calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    T4WallPenaltyScaleDiagnostics* scale_diag);

int calc_wall_nitsche_gamma_sweep_decomposed_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu,
    const double Uw[3],
    const double eU_in[3],
    const double eP_in[3],
    const double* gamma_list,
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    T4WallPenaltyScaleDiagnostics* scale_diag);

void wall_gamma_sweep_decomposed_allreduce(
    int num_gamma,
    double* D_penalty_by_gamma,
    double* L_penalty_by_gamma,
    double* D_tangent_by_gamma,
    double* L_tangent_by_gamma,
    double* D_total_by_gamma,
    double* L_total_by_gamma,
    double* D_penalty_mu_by_gamma,
    double* L_penalty_mu_by_gamma,
    double* D_penalty_dt_by_gamma,
    double* L_penalty_dt_by_gamma,
    MONOLIS_COM* monolis_com);

void wall_penalty_scale_diagnostics_allreduce(
    T4WallPenaltyScaleDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void output_wall_gamma_sweep_penalty_decomposition_tet4_tri3(
    int num_gamma,
    const double* gamma_list,
    const double* D_penalty_mu_by_gamma,
    const double* L_penalty_mu_by_gamma,
    const double* D_penalty_dt_by_gamma,
    const double* L_penalty_dt_by_gamma,
    double t,
    const char* directory);

void output_wall_penalty_scale_diagnostics_tet4_tri3(
    const T4WallPenaltyScaleDiagnostics* diag,
    double t,
    const char* directory);

int calc_wall_diagnostics_tet4_tri3_with_Uw(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    T4WallDiagnostics* diag);

int calc_wall_diagnostics_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double Uref,
    double p_ref,
    const double eU_in[3],
    const double eP_in[3],
    T4WallDiagnostics* diag);


int output_wall_max_error_locations_tet4_tri3(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double Uref,
    const double Uw_in[3],
    MONOLIS_COM* monolis_com,
    double t,
    const char* directory);

int output_domain_max_div_u_location_tet4(
    const BBFE_DATA* fe,
    const VALUES* vals,
    MONOLIS_COM* monolis_com,
    double t,
    const char* directory);

int calc_nitsche_row_sum_reaction_tet4_tri3_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw_in[3],
    double gamma_n,
    T4NitscheRowSumDiagnostics* diag);

static int t4_all_ranks_status_ok(
    int local_status,
    MONOLIS_COM* monolis_com,
    const char* label);

void nitsche_row_sum_diagnostics_allreduce(
    T4NitscheRowSumDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void solver_fom_VMS(
    FE_SYSTEM        sys,
    const BBFE_DATA* surf_wall,
    const BBFE_BASIS* basis_wall,
    double           t,
    const int        step);
