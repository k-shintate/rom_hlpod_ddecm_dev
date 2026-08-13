#pragma once

#include "core_FOM.h"


#define ST_SCHEME_CG1      1

#define ST_TIME_FAMILY_DG  0
#define ST_TIME_FAMILY_CG  1

#define ST_TIME_MAX_DOF    8
#define ST_TIME_MAX_IP     8


typedef struct {
	int family;
	int degree;
	int n_dof;

	/*
	 * is_nodal = 1:
	 *   time coefficient alpha is the value at dof_tau[alpha].
	 */
	int is_nodal;
	double dof_tau[ST_TIME_MAX_DOF];

	/*
	 * Temporal matrices:
	 *
	 * Mt[alpha][beta]
	 *   = integral psi_alpha * psi_beta d tau
	 *
	 * Dt[alpha][beta]
	 *   = integral psi_alpha * d(psi_beta)/d tau d tau
	 */
	double Mt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double Dt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];

	/*
	 * Endpoint values:
	 *
	 * e_minus[alpha] = psi_alpha(0)
	 * e_plus [alpha] = psi_alpha(1)
	 */
	double e_minus[ST_TIME_MAX_DOF];
	double e_plus [ST_TIME_MAX_DOF];

	/*
	 * dG flux matrices:
	 *
	 * E_self = e_minus * e_minus^T
	 * E_prev = e_minus * e_plus^T
	 */
	double E_self[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double E_prev [ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];

	/*
	 * Temporal quadrature for the source term.
	 */
	int num_time_ip;
	double tau[ST_TIME_MAX_IP];
	double weight[ST_TIME_MAX_IP];
	double psi[ST_TIME_MAX_DOF][ST_TIME_MAX_IP];

} ST_TIME_ELEMENT;


typedef struct {
	int type;

	double lhs_mass;
	double lhs_operator;

	double rhs_old_mass;
	double rhs_old_operator;

	double rhs_source_old;
	double rhs_source_new;

} ST_TIME_SCHEME;



/* ------------------------------------------------------------------------- */
/* Time-element initialization                                                */
/* ------------------------------------------------------------------------- */

void ST_TimeElement_init_dG(
		ST_TIME_ELEMENT* te,
		const int        degree);


/* ------------------------------------------------------------------------- */
/* Single-slab dG block                                                       */
/* ------------------------------------------------------------------------- */

typedef double (*ST_DIRICHLET_VALUE_FUNC)(
		const double* x,
		int           component,
		double        time,
		void*         user_data);

void BBFE_sys_monowrap_set_Dirichlet_bc_space_time(
		MONOLIS*                monolis,
		BBFE_DATA*              fe,
		const ST_TIME_ELEMENT*  te,
		const int               num_window_slabs,
		const int               num_components,
		const BBFE_BC*          bc,
		const double            t_window_start,
		const double            dt,
		ST_DIRICHLET_VALUE_FUNC value_func,
		void*                   user_data,
		double*                 g_rhs);

void set_element_mat_ST_dG_block(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		BBFE_BASIS*            basis,
		VALUES*                vals,
		const ST_TIME_ELEMENT* te);

void set_element_vec_ST_dG_block(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		BBFE_BASIS*            basis,
		VALUES*                vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st,
		double                 t_old,
		double                 t_new);

void set_element_mat_ST_dG0_block(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals);

void set_element_vec_ST_dG0_block(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_old,
		double        t_new);

void set_element_mat_ST_dG1_block(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals);

void set_element_vec_ST_dG1_block(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_old,
		double        t_new);

void set_element_mat_ST_dG2_block(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals);

void set_element_vec_ST_dG2_block(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_old,
		double        t_new);

static void check_ST_Dirichlet_values(
        const FE_SYSTEM*       sys,
        const ST_TIME_ELEMENT* te,
        const int              num_window_slabs,
        const double           t_window_start,
        const double*          T_window);

/* ------------------------------------------------------------------------- */
/* Temporal trace conversion                                                  */
/* ------------------------------------------------------------------------- */

void ST_TimeElement_set_prev_trace_from_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const double*          T_nodal,
		double*                T_prev_st);

void ST_TimeElement_extract_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const double*          T_st,
		double*                T_nodal);


/* ------------------------------------------------------------------------- */
/* cG(1) / Crank-Nicolson                                                     */
/* ------------------------------------------------------------------------- */

void set_element_mat_ST_1step(
		MONOLIS*        monolis,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
		ST_TIME_SCHEME* scheme);

void set_element_vec_ST_1step(
		MONOLIS*        monolis,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
		ST_TIME_SCHEME* scheme,
		double          t_old,
		double          t_new);

void set_element_mat_ST_cG1(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals);

void set_element_vec_ST_cG1(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		double      t_new);

void solver_fom_ST_1step(
		FE_SYSTEM       sys,
		ST_TIME_SCHEME* scheme,
		double          t_old,
		double          t_new,
		const int       step);

void solver_fom_ST_cG1(
		FE_SYSTEM sys,
		double    t_new,
		const int step);


/* ------------------------------------------------------------------------- */
/* Windowed space-time dG                                                     */
/* ------------------------------------------------------------------------- */

int ST_window_gid(
		const int slab_id,
		const int time_dof,
		const int space_id,
		const int num_time_dofs,
		const int num_space_nodes);

void set_element_mat_ST_dG_window(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		BBFE_BASIS*            basis,
		VALUES*                vals,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs);

void set_element_vec_ST_dG_window(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		BBFE_BASIS*            basis,
		VALUES*                vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st,
		double                 t_window_start,
		const int              num_window_slabs);

void set_element_mat_ST_dG_window_degree(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		const int   degree,
		const int   num_window_slabs);

void set_element_vec_ST_dG_window_degree(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const int     degree,
		const double* T_prev_st,
		double        t_window_start,
		const int     num_window_slabs);

void set_element_mat_ST_dG0_window(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		const int   num_window_slabs);

void set_element_vec_ST_dG0_window(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_window_start,
		const int     num_window_slabs);

void set_element_mat_ST_dG1_window(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		const int   num_window_slabs);

void set_element_vec_ST_dG1_window(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_window_start,
		const int     num_window_slabs);

void set_element_mat_ST_dG2_window(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		VALUES*     vals,
		const int   num_window_slabs);

void set_element_vec_ST_dG2_window(
		MONOLIS*      monolis,
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		VALUES*       vals,
		const double* T_prev_st,
		double        t_window_start,
		const int     num_window_slabs);


/* ------------------------------------------------------------------------- */
/* Window extraction                                                         */
/* ------------------------------------------------------------------------- */

void ST_Window_extract_slab_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const int              slab,
		const double*          T_window,
		double*                T_nodal);

void ST_Window_extract_last_slab_state(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const double*          T_window,
		double*                T_last_st);

void ST_Window_extract_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const double*          T_window,
		double*                T_nodal);


/* ------------------------------------------------------------------------- */
/* ST graph, boundary condition and solver                                    */
/* ------------------------------------------------------------------------- */

void set_nonzero_pattern_space_time(
		MONOLIS*               monolis,
		const char*            label,
		const char*            directory,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs);

void BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		VALUES*                vals,
		BBFE_BC*               bc,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs,
		const double           t_window_start,
		double*                g_rhs);

void solver_fom_space_time(
		FE_SYSTEM*             sys,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs,
		const double           t_window_start,
		const int              window_step,
		const double*          T_prev_st,
		double*                T_window,
		double*                T_last_st);
