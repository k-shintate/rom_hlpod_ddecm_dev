#include "core_FOM.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Unified time schemes for heat equation
 *
 * Semi-discrete form:
 *     M dT/dt + K T = F
 *
 * One-step wrapper:
 *     cG(1) / Crank-Nicolson only.
 *
 * Space-time dG schemes, including dG(0), are assembled through the
 * time-block and windowed space-time assemblers below.
 *----------------------------------------------------------------------------*/

#ifndef ST_TIME_SCHEME_H
#define ST_TIME_SCHEME_H

#define ST_SCHEME_CG1 1

#define ST_TIME_FAMILY_DG 0
#define ST_TIME_FAMILY_CG 1
#define ST_TIME_MAX_DOF  8
#define ST_TIME_MAX_IP   8

typedef struct {
	int family;
	int degree;
	int n_dof;

	/*
	 * Reference-time matrices on tau in [0,1].
	 *
	 * Mt[alpha][beta] = int_0^1 psi_alpha psi_beta d tau
	 * Dt[alpha][beta] = int_0^1 psi_alpha d(psi_beta)/d tau d tau
	 */
	double Mt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double Dt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];

	/* Endpoint vectors. */
	double e_minus[ST_TIME_MAX_DOF];  /* psi_alpha(0) */
	double e_plus [ST_TIME_MAX_DOF];  /* psi_alpha(1) */

	/*
	 * dG upwind time-flux matrices.
	 *
	 * E_self = e_minus e_minus^T
	 * E_prev = e_minus e_plus^T
	 */
	double E_self[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double E_prev [ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];


	/*
	 * Time quadrature for the source term on tau in [0,1].
	 * psi[alpha][q] = psi_alpha(tau_q).
	 */
	int num_time_ip;
	double tau[ST_TIME_MAX_IP];
	double weight[ST_TIME_MAX_IP];
	double psi[ST_TIME_MAX_DOF][ST_TIME_MAX_IP];

} ST_TIME_ELEMENT;


static void ST_TimeElement_clear(
		ST_TIME_ELEMENT* te)
{
	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof = 0;
	te->num_time_ip = 0;

	for(int q=0; q<ST_TIME_MAX_IP; q++) {
		te->tau[q] = 0.0;
		te->weight[q] = 0.0;
	}

	for(int i=0; i<ST_TIME_MAX_DOF; i++) {
		te->e_minus[i] = 0.0;
		te->e_plus [i] = 0.0;
		for(int q=0; q<ST_TIME_MAX_IP; q++) {
			te->psi[i][q] = 0.0;
		}

		for(int j=0; j<ST_TIME_MAX_DOF; j++) {
			te->Mt[i][j] = 0.0;
			te->Dt[i][j] = 0.0;
			te->E_self[i][j] = 0.0;
			te->E_prev [i][j] = 0.0;
		}
	}
}


/*
 * dG(0) TimeElement:
 *
 *   psi_0(tau) = 1
 *
 * With the validated code scaling, the resulting single-slab block equation is
 *
 *   ((D_t + E_self)/dt \otimes M + M_t \otimes K) U
 *   =
 *   (E_prev/dt \otimes M) U_prev + source
 *
 * Since dG(0) has only one time DOF:
 *
 *   Mt       = [1]
 *   Dt       = [0]
 *   E_self   = [1]
 *   E_prev   = [1]
 */
static void ST_TimeElement_init_dG0(
		ST_TIME_ELEMENT* te)
{
	ST_TimeElement_clear(te);

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof  = 1;

	te->Mt[0][0] = 1.0;
	te->Dt[0][0] = 0.0;

	te->e_minus[0] = 1.0;
	te->e_plus [0] = 1.0;

	te->E_self[0][0] = 1.0;
	te->E_prev [0][0] = 1.0;

	/* Right-endpoint quadrature gives BE-compatible F(t_{n+1}). */
	te->num_time_ip = 1;
	te->tau[0] = 1.0;
	te->weight[0] = 1.0;
	te->psi[0][0] = 1.0;
}


/*
 * dG(1) TimeElement with nodal basis on [0,1]:
 *
 *   psi_0(tau) = 1 - tau
 *   psi_1(tau) = tau
 *
 * Matrix convention:
 *   Mt[alpha][beta] = int psi_alpha psi_beta d tau
 *   Dt[alpha][beta] = int psi_alpha d(psi_beta)/d tau d tau
 *
 * dG upwind flux:
 *   E_self = e_minus e_minus^T
 *   E_prev = e_minus e_plus^T
 */
static void ST_TimeElement_init_dG1(
		ST_TIME_ELEMENT* te)
{
	ST_TimeElement_clear(te);

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 1;
	te->n_dof  = 2;

	/* Mass matrix on reference time interval. */
	te->Mt[0][0] = 1.0 / 3.0;
	te->Mt[0][1] = 1.0 / 6.0;
	te->Mt[1][0] = 1.0 / 6.0;
	te->Mt[1][1] = 1.0 / 3.0;

	/* Derivative matrix: int psi_alpha * d psi_beta / d tau. */
	te->Dt[0][0] = -0.5;
	te->Dt[0][1] =  0.5;
	te->Dt[1][0] = -0.5;
	te->Dt[1][1] =  0.5;

	te->e_minus[0] = 1.0;
	te->e_minus[1] = 0.0;
	te->e_plus [0] = 0.0;
	te->e_plus [1] = 1.0;

	for(int a=0; a<te->n_dof; a++) {
		for(int b=0; b<te->n_dof; b++) {
			te->E_self[a][b] = te->e_minus[a] * te->e_minus[b];
			te->E_prev [a][b] = te->e_minus[a] * te->e_plus [b];
		}
	}

	/* Two-point Gauss quadrature on [0,1] for source integration. */
	const double s = 0.57735026918962576451;
	te->num_time_ip = 2;
	te->tau[0] = 0.5 * (1.0 - s);
	te->tau[1] = 0.5 * (1.0 + s);
	te->weight[0] = 0.5;
	te->weight[1] = 0.5;

	for(int q=0; q<te->num_time_ip; q++) {
		double tau = te->tau[q];
		te->psi[0][q] = 1.0 - tau;
		te->psi[1][q] = tau;
	}
}


/*
 * dG(2) TimeElement with nodal quadratic basis on [0,1]
 * at tau = 0, 1/2, 1:
 *
 *   psi_0(tau) = 2 tau^2 - 3 tau + 1
 *   psi_1(tau) = -4 tau^2 + 4 tau
 *   psi_2(tau) = 2 tau^2 - tau
 *
 * Matrix convention:
 *   Mt[alpha][beta] = int psi_alpha psi_beta d tau
 *   Dt[alpha][beta] = int psi_alpha d(psi_beta)/d tau d tau
 *
 * dG upwind flux:
 *   E_self = e_minus e_minus^T
 *   E_prev = e_minus e_plus^T
 */
static void ST_TimeElement_init_dG2(
		ST_TIME_ELEMENT* te)
{
	ST_TimeElement_clear(te);

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 2;
	te->n_dof  = 3;

	/* Mass matrix on reference time interval. */
	te->Mt[0][0] =  2.0 / 15.0;
	te->Mt[0][1] =  1.0 / 15.0;
	te->Mt[0][2] = -1.0 / 30.0;

	te->Mt[1][0] =  1.0 / 15.0;
	te->Mt[1][1] =  8.0 / 15.0;
	te->Mt[1][2] =  1.0 / 15.0;

	te->Mt[2][0] = -1.0 / 30.0;
	te->Mt[2][1] =  1.0 / 15.0;
	te->Mt[2][2] =  2.0 / 15.0;

	/* Derivative matrix: int psi_alpha * d psi_beta / d tau. */
	te->Dt[0][0] = -1.0 / 2.0;
	te->Dt[0][1] =  2.0 / 3.0;
	te->Dt[0][2] = -1.0 / 6.0;

	te->Dt[1][0] = -2.0 / 3.0;
	te->Dt[1][1] =  0.0;
	te->Dt[1][2] =  2.0 / 3.0;

	te->Dt[2][0] =  1.0 / 6.0;
	te->Dt[2][1] = -2.0 / 3.0;
	te->Dt[2][2] =  1.0 / 2.0;

	te->e_minus[0] = 1.0;
	te->e_minus[1] = 0.0;
	te->e_minus[2] = 0.0;

	te->e_plus [0] = 0.0;
	te->e_plus [1] = 0.0;
	te->e_plus [2] = 1.0;

	for(int a=0; a<te->n_dof; a++) {
		for(int b=0; b<te->n_dof; b++) {
			te->E_self[a][b] = te->e_minus[a] * te->e_minus[b];
			te->E_prev [a][b] = te->e_minus[a] * te->e_plus [b];
		}
	}

	/* Three-point Gauss quadrature on [0,1] for source integration. */
	const double s = 0.77459666924148337704; /* sqrt(3/5) */
	te->num_time_ip = 3;
	te->tau[0] = 0.5 * (1.0 - s);
	te->tau[1] = 0.5;
	te->tau[2] = 0.5 * (1.0 + s);
	te->weight[0] = 5.0 / 18.0;
	te->weight[1] = 4.0 / 9.0;
	te->weight[2] = 5.0 / 18.0;

	for(int q=0; q<te->num_time_ip; q++) {
		double tau = te->tau[q];
		te->psi[0][q] =  2.0 * tau * tau - 3.0 * tau + 1.0;
		te->psi[1][q] = -4.0 * tau * tau + 4.0 * tau;
		te->psi[2][q] =  2.0 * tau * tau - tau;
	}
}


static void ST_TimeElement_init_dG(
		ST_TIME_ELEMENT* te,
		const int degree)
{
	if(degree == 0) {
		ST_TimeElement_init_dG0(te);
	}
	else if(degree == 1) {
		ST_TimeElement_init_dG1(te);
	}
	else if(degree == 2) {
		ST_TimeElement_init_dG2(te);
	}
	else {
		printf("ERROR: unsupported dG degree %d in ST_TimeElement_init_dG.\n", degree);
		exit(1);
	}
}


typedef struct {
	int type;

	/* Left hand side:
	 * A = lhs_mass * M + lhs_operator * K
	 */
	double lhs_mass;
	double lhs_operator;

	/* Right hand side old solution terms:
	 * rhs = rhs_old_mass * M*T_old
	 *     + rhs_old_operator * K*T_old
	 */
	double rhs_old_mass;
	double rhs_old_operator;

	/* Source terms:
	 * rhs += rhs_source_old * F(t_old)
	 *      + rhs_source_new * F(t_new)
	 */
	double rhs_source_old;
	double rhs_source_new;

} ST_TIME_SCHEME;


/*
 * One-step coefficients are intentionally kept cG-only.
 * dG(0) is handled by the space-time dG block/window assemblers.
 */
static void ST_TimeScheme_set_cG1(
		ST_TIME_SCHEME* scheme,
		const double dt)
{
	scheme->type = ST_SCHEME_CG1;

	/*
	 * cG(1) / Crank-Nicolson:
	 *
	 * (1/dt M + 1/2 K) T^{n+1}
	 * =
	 * (1/dt M - 1/2 K) T^n
	 * + 1/2 F^n + 1/2 F^{n+1}
	 */
	scheme->lhs_mass     = 1.0 / dt;
	scheme->lhs_operator = 0.5;

	scheme->rhs_old_mass     = 1.0 / dt;
	scheme->rhs_old_operator = -0.5;

	scheme->rhs_source_old = 0.5;
	scheme->rhs_source_new = 0.5;
}


static const char* ST_TimeScheme_name(
		const ST_TIME_SCHEME* scheme)
{
	if(scheme->type == ST_SCHEME_CG1) return "cG(1) / Crank-Nicolson";

	return "unsupported one-step scheme";
}

#endif


/*----------------------------------------------------------------------------
 * Matrix terms
 *----------------------------------------------------------------------------*/

static void add_heat_mass_matrix(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double*  a_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	(void)vals;

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {
				for(int p=0; p<np; p++) {
					val_ip[p] =
						coef *
						BBFE_elemmat_convdiff_mat_mass(
							basis->N[p][i],
							basis->N[p][j],
							a_ip[p]);
				}

				double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
					monolis,
					fe->conn[e][i], fe->conn[e][j],
					0, 0,
					integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
}


static void add_heat_diffusion_matrix(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	(void)vals;

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {
				for(int p=0; p<np; p++) {
					/*
					 * Keep the same sign convention as the validated code.
					 */
					val_ip[p] =
						coef *
						(
							- BBFE_elemmat_convdiff_mat_diff(
								fe->geo[e][p].grad_N[i],
								fe->geo[e][p].grad_N[j],
								k_ip[p])
						);
				}

				double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
					monolis,
					fe->conn[e][i], fe->conn[e][j],
					0, 0,
					integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
}


/*----------------------------------------------------------------------------
 * RHS terms
 *----------------------------------------------------------------------------*/

static void add_heat_old_mass_vector(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	double*  local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;
	double*  a_ip;
	double*  T_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, vals->T, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] =
					coef *
					BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i],
						T_ip[p],
						a_ip[p]);
			}

			double integ_val = BBFE_std_integ_calc(
				np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[fe->conn[e][i]] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_1d_double(local_T, nl);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
}


static void add_heat_old_diffusion_vector(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;

	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	double*  local_T;

	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;
	double*  k_ip;

	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(
				local_x, fe, fe->x, e, 3);

		BBFE_elemmat_set_local_array_scalar(
				local_T, fe, vals->T, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(
					x_ip[p], nl, local_x, basis->N[p]);

			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;

				for(int j=0; j<nl; j++) {
					/*
					 * Same sign convention as add_heat_diffusion_matrix().
					 */
					double Kij =
						- BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i],
							fe->geo[e][p].grad_N[j],
							k_ip[p]);

					val_ip[p] += coef * Kij * local_T[j];
				}
			}

			double integ_val = BBFE_std_integ_calc(
					np,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);

			monolis->mat.R.B[fe->conn[e][i]] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_1d_double(local_T, nl);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
}


static void add_heat_source_vector(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       t,
		double       coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double** v_ip;
	double*  a_ip;
	double*  k_ip;
	double*  f_ip;

	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	(void)vals;

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);

			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

			f_ip[p] = manusol_get_source(
				x_ip[p],
				t,
				a_ip[p],
				v_ip[p],
				k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] =
					coef *
					BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i],
						f_ip[p]);
			}

			double integ_val = BBFE_std_integ_calc(
				np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[fe->conn[e][i]] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(f_ip, np);
}


/*----------------------------------------------------------------------------
 * Time-block dG slab assemblers
 *
 * These routines assemble one time slab with n_t = degree + 1 time DOFs.
 * They include dG(0) as n_t = 1, dG(1) as n_t = 2, and dG(2) as n_t = 3.
 *
 * Unknown ordering for one slab:
 *
 *   T_st[alpha * Nx + i]
 *
 * where alpha is the local time DOF and i is the spatial node ID.
 * The MONOLIS matrix must therefore be allocated for n_t * Nx scalar DOFs.
 *----------------------------------------------------------------------------*/

static int ST_block_gid(
		const int time_dof,
		const int space_id,
		const int num_space_nodes)
{
	return time_dof * num_space_nodes + space_id;
}


void set_element_mat_ST_dG_block(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double*  a_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int alpha=0; alpha<nt; alpha++) {
			for(int beta=0; beta<nt; beta++) {
				double coef_mass =
					(te->Dt[alpha][beta] + te->E_self[alpha][beta]) / vals->dt;

				double coef_diff = te->Mt[alpha][beta];

				for(int i=0; i<nl; i++) {
					for(int j=0; j<nl; j++) {
						for(int p=0; p<np; p++) {
							double Mij =
								BBFE_elemmat_convdiff_mat_mass(
									basis->N[p][i],
									basis->N[p][j],
									a_ip[p]);

							double Kij =
								- BBFE_elemmat_convdiff_mat_diff(
									fe->geo[e][p].grad_N[i],
									fe->geo[e][p].grad_N[j],
									k_ip[p]);

							val_ip[p] = coef_mass * Mij + coef_diff * Kij;
						}

						double integ_val = BBFE_std_integ_calc(
							np,
							val_ip,
							basis->integ_weight,
							Jacobian_ip);

						int row = ST_block_gid(alpha, fe->conn[e][i], nx);
						int col = ST_block_gid(beta,  fe->conn[e][j], nx);

						monolis_add_scalar_to_sparse_matrix_R(
							monolis,
							row,
							col,
							0,
							0,
							integ_val);
					}
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}


static void add_heat_prev_mass_vector_ST_dG_block(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** local_T_prev;
	local_T_prev = BB_std_calloc_2d_double(local_T_prev, nt, nl);

	double** x_ip;
	double*  a_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	(void)vals;

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int beta=0; beta<nt; beta++) {
			BBFE_elemmat_set_local_array_scalar(
					local_T_prev[beta],
					fe,
					&(T_prev_st[beta * nx]),
					e);
		}

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
		}

		for(int alpha=0; alpha<nt; alpha++) {
			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					for(int beta=0; beta<nt; beta++) {
						double coef = te->E_prev[alpha][beta] / vals->dt;

						if(fabs(coef) < 1.0e-30) continue;

						double T_beta_ip = BBFE_std_mapping_scalar(
								nl,
								local_T_prev[beta],
								basis->N[p]);

						val_ip[p] +=
							coef *
							BBFE_elemmat_convdiff_vec_mass(
								basis->N[p][i],
								T_beta_ip,
								a_ip[p]);
					}
				}

				double integ_val = BBFE_std_integ_calc(
						np,
						val_ip,
						basis->integ_weight,
						Jacobian_ip);

				int row = ST_block_gid(alpha, fe->conn[e][i], nx);
				monolis->mat.R.B[row] += integ_val;
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(local_T_prev, nt, nl);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
}


static void add_heat_source_vector_ST_dG_block(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		double                 t_old,
		double                 t_new)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double** v_ip;
	double*  a_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	(void)vals;

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int alpha=0; alpha<nt; alpha++) {
			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					for(int q=0; q<te->num_time_ip; q++) {
						double tau = te->tau[q];
						double tq  = t_old + tau * (t_new - t_old);
						double coef = te->weight[q] * te->psi[alpha][q];

						if(fabs(coef) < 1.0e-30) continue;

						double f_ip = manusol_get_source(
								x_ip[p],
								tq,
								a_ip[p],
								v_ip[p],
								k_ip[p]);

						val_ip[p] +=
							coef *
							BBFE_elemmat_convdiff_vec_source(
								basis->N[p][i],
								f_ip);
					}
				}

				double integ_val = BBFE_std_integ_calc(
						np,
						val_ip,
						basis->integ_weight,
						Jacobian_ip);

				int row = ST_block_gid(alpha, fe->conn[e][i], nx);
				monolis->mat.R.B[row] += integ_val;
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}


void set_element_vec_ST_dG_block(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st,
		double                 t_old,
		double                 t_new)
{
	add_heat_prev_mass_vector_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			te,
			T_prev_st);

	add_heat_source_vector_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			te,
			t_old,
			t_new);
}


void ST_TimeElement_set_prev_trace_from_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const double*          T_nodal,
		double*                T_prev_st)
{
	for(int alpha=0; alpha<te->n_dof; alpha++) {
		for(int i=0; i<num_space_nodes; i++) {
			T_prev_st[alpha * num_space_nodes + i] = 0.0;
		}
	}

	/*
	 * For the current nodal dG(0)/dG(1) bases, the right endpoint is the last
	 * local time DOF. This gives sum_beta e_plus[beta] T_prev_beta = T_nodal.
	 */
	int beta_right = te->n_dof - 1;
	for(int i=0; i<num_space_nodes; i++) {
		T_prev_st[beta_right * num_space_nodes + i] = T_nodal[i];
	}
}


void ST_TimeElement_extract_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const double*          T_st,
		double*                T_nodal)
{
	for(int i=0; i<num_space_nodes; i++) {
		T_nodal[i] = 0.0;
		for(int alpha=0; alpha<te->n_dof; alpha++) {
			T_nodal[i] += te->e_plus[alpha] * T_st[alpha * num_space_nodes + i];
		}
	}
}


/*----------------------------------------------------------------------------
 * cG(1) one-step assemblers
 *----------------------------------------------------------------------------*/

void set_element_mat_ST_1step(
		MONOLIS*        monolis,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
		ST_TIME_SCHEME* scheme)
{
	/*
	 * A = lhs_mass * M + lhs_operator * K
	 */
	add_heat_mass_matrix(
			monolis,
			fe,
			basis,
			vals,
			scheme->lhs_mass);

	add_heat_diffusion_matrix(
			monolis,
			fe,
			basis,
			vals,
			scheme->lhs_operator);
}


void set_element_vec_ST_1step(
		MONOLIS*        monolis,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
		ST_TIME_SCHEME* scheme,
		double          t_old,
		double          t_new)
{
	const double eps = 1.0e-30;

	/*
	 * rhs += rhs_old_mass * M * T_old
	 */
	if(fabs(scheme->rhs_old_mass) > eps) {
		add_heat_old_mass_vector(
				monolis,
				fe,
				basis,
				vals,
				scheme->rhs_old_mass);
	}

	/*
	 * rhs += rhs_old_operator * K * T_old
	 */
	if(fabs(scheme->rhs_old_operator) > eps) {
		add_heat_old_diffusion_vector(
				monolis,
				fe,
				basis,
				vals,
				scheme->rhs_old_operator);
	}

	/*
	 * rhs += rhs_source_old * F(t_old)
	 */
	if(fabs(scheme->rhs_source_old) > eps) {
		add_heat_source_vector(
				monolis,
				fe,
				basis,
				vals,
				t_old,
				scheme->rhs_source_old);
	}

	/*
	 * rhs += rhs_source_new * F(t_new)
	 */
	if(fabs(scheme->rhs_source_new) > eps) {
		add_heat_source_vector(
				monolis,
				fe,
				basis,
				vals,
				t_new,
				scheme->rhs_source_new);
	}
}


void set_element_mat_ST_dG0_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG0(&te);

	set_element_mat_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te);
}


void set_element_vec_ST_dG0_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_old,
		double       t_new)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG0(&te);

	set_element_vec_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te,
			T_prev_st,
			t_old,
			t_new);
}


void set_element_mat_ST_dG1_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG1(&te);

	set_element_mat_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te);
}


void set_element_vec_ST_dG1_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_old,
		double       t_new)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG1(&te);

	set_element_vec_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te,
			T_prev_st,
			t_old,
			t_new);
}


void set_element_mat_ST_dG2_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG2(&te);

	set_element_mat_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te);
}


void set_element_vec_ST_dG2_block(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_old,
		double       t_new)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG2(&te);

	set_element_vec_ST_dG_block(
			monolis,
			fe,
			basis,
			vals,
			&te,
			T_prev_st,
			t_old,
			t_new);
}


/* cG(1) one-step wrappers. */

void set_element_mat_ST_cG1(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	ST_TIME_SCHEME scheme;
	ST_TimeScheme_set_cG1(&scheme, vals->dt);

	set_element_mat_ST_1step(
			monolis,
			fe,
			basis,
			vals,
			&scheme);
}


void set_element_vec_ST_cG1(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       t_new)
{
	ST_TIME_SCHEME scheme;
	ST_TimeScheme_set_cG1(&scheme, vals->dt);

	set_element_vec_ST_1step(
			monolis,
			fe,
			basis,
			vals,
			&scheme,
			t_new - vals->dt,
			t_new);
}


/*----------------------------------------------------------------------------
 * cG(1) one-step solver
 *----------------------------------------------------------------------------*/

void solver_fom_ST_1step(
		FE_SYSTEM       sys,
		ST_TIME_SCHEME* scheme,
		double          t_old,
		double          t_new,
		const int       step)
{
	printf("\n%s ----------------- step %d : %s ----------------\n",
			CODENAME,
			step,
			ST_TimeScheme_name(scheme));

	if(scheme->type != ST_SCHEME_CG1) {
		printf("ERROR: solver_fom_ST_1step supports only cG(1).\n");
		exit(1);
	}

	/*
	 * monolis0 contains the constant left-hand-side matrix.
	 */
	monolis_copy_mat_value_R(
			&(sys.monolis0),
			&(sys.monolis));

	monolis_clear_mat_value_rhs_R(
			&(sys.monolis));

	/*
	 * Build cG(1) / Crank-Nicolson RHS:
	 *   b = (1/dt) M T_old
	 *     - (1/2) K T_old
	 *     + (1/2) F_old
	 *     + (1/2) F_new
	 */
	set_element_vec_ST_1step(
			&(sys.monolis),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals),
			scheme,
			t_old,
			t_new);

	/*
	 * Dirichlet BC at new time.
	 */
	manusol_set_theo_sol(
			&(sys.fe),
			sys.vals.theo_sol,
			t_new);

	BBFE_manusol_set_bc_scalar(
			&(sys.fe),
			&(sys.bc),
			sys.vals.theo_sol,
			t_new);

	BBFE_sys_monowrap_set_Dirichlet_bc(
			&(sys.monolis),
			sys.fe.total_num_nodes,
			BLOCK_SIZE,
			&(sys.bc),
			sys.monolis.mat.R.B);

	BBFE_sys_monowrap_solve(
			&(sys.monolis),
			&(sys.monolis_com),
			sys.vals.T,
			MONOLIS_ITER_CG,
			MONOLIS_PREC_DIAG,
			sys.vals.mat_max_iter,
			sys.vals.mat_epsilon);
}


void solver_fom_ST_cG1(
		FE_SYSTEM sys,
		double    t_new,
		const int step)
{
	ST_TIME_SCHEME scheme;
	ST_TimeScheme_set_cG1(&scheme, sys.vals.dt);

	solver_fom_ST_1step(
			sys,
			&scheme,
			t_new - sys.vals.dt,
			t_new,
			step);
}


/*----------------------------------------------------------------------------
 * Windowed space-time dG assemblers
 *
 * These routines assemble a window containing multiple consecutive dG slabs.
 * They reuse the same TimeElement used by the single-slab dG block assembler.
 *
 * Unknown ordering in one window:
 *
 *   T_window[((slab * nt + alpha) * Nx) + i]
 *
 * where
 *   slab  = 0, ..., num_window_slabs-1
 *   alpha = 0, ..., nt-1
 *   i     = 0, ..., Nx-1
 *
 * The MONOLIS matrix must be allocated for
 *
 *   num_window_slabs * te->n_dof * fe->total_num_nodes
 *
 * scalar unknowns.
 *
 * Window matrix form:
 *
 *   [ A_0                         ] [U_0]   [F_0 + B_0 U_prev]
 *   [-B_1  A_1                    ] [U_1] = [F_1             ]
 *   [      -B_2 A_2               ] [U_2]   [F_2             ]
 *   [             ...             ] [...]   [...             ]
 *
 * with
 *
 *   A_s = ((Dt + E_self)/dt) \otimes M + Mt \otimes K
 *   B_s = (E_prev/dt) \otimes M
 *
 * For constant dt and time-independent coefficients, A_s and B_s are the same
 * for all slabs. Source terms are still evaluated at each slab time.
 *----------------------------------------------------------------------------*/

static int ST_window_gid(
		const int slab_id,
		const int time_dof,
		const int space_id,
		const int num_time_dofs,
		const int num_space_nodes)
{
	return (slab_id * num_time_dofs + time_dof) * num_space_nodes + space_id;
}


void set_element_mat_ST_dG_window(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		const int             num_window_slabs)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	if(num_window_slabs <= 0) {
		printf("ERROR: num_window_slabs must be positive. got %d\n", num_window_slabs);
		exit(1);
	}

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double*  a_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		/* Diagonal slab blocks A_s. */
		for(int slab=0; slab<num_window_slabs; slab++) {
			for(int alpha=0; alpha<nt; alpha++) {
				for(int beta=0; beta<nt; beta++) {
					double coef_mass =
						(te->Dt[alpha][beta] + te->E_self[alpha][beta]) / vals->dt;

					double coef_diff = te->Mt[alpha][beta];

					for(int i=0; i<nl; i++) {
						for(int j=0; j<nl; j++) {
							for(int p=0; p<np; p++) {
								double Mij =
									BBFE_elemmat_convdiff_mat_mass(
										basis->N[p][i],
										basis->N[p][j],
										a_ip[p]);

								double Kij =
									- BBFE_elemmat_convdiff_mat_diff(
										fe->geo[e][p].grad_N[i],
										fe->geo[e][p].grad_N[j],
										k_ip[p]);

								val_ip[p] = coef_mass * Mij + coef_diff * Kij;
							}

							double integ_val = BBFE_std_integ_calc(
								np,
								val_ip,
								basis->integ_weight,
								Jacobian_ip);

							int row = ST_window_gid(slab, alpha, fe->conn[e][i], nt, nx);
							int col = ST_window_gid(slab, beta,  fe->conn[e][j], nt, nx);

							monolis_add_scalar_to_sparse_matrix_R(
								monolis,
								row,
								col,
								0,
								0,
								integ_val);
						}
					}
				}
			}
		}

		/* Subdiagonal inter-slab blocks: -B_s. */
		for(int slab=1; slab<num_window_slabs; slab++) {
			for(int alpha=0; alpha<nt; alpha++) {
				for(int beta=0; beta<nt; beta++) {
					double coef_mass = - te->E_prev[alpha][beta] / vals->dt;

					if(fabs(coef_mass) < 1.0e-30) continue;

					for(int i=0; i<nl; i++) {
						for(int j=0; j<nl; j++) {
							for(int p=0; p<np; p++) {
								double Mij =
									BBFE_elemmat_convdiff_mat_mass(
										basis->N[p][i],
										basis->N[p][j],
										a_ip[p]);

								val_ip[p] = coef_mass * Mij;
							}

							double integ_val = BBFE_std_integ_calc(
								np,
								val_ip,
								basis->integ_weight,
								Jacobian_ip);

							int row = ST_window_gid(slab,     alpha, fe->conn[e][i], nt, nx);
							int col = ST_window_gid(slab - 1, beta,  fe->conn[e][j], nt, nx);

							monolis_add_scalar_to_sparse_matrix_R(
								monolis,
								row,
								col,
								0,
								0,
								integ_val);
						}
					}
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}


static void add_heat_prev_mass_vector_ST_dG_window_first_slab(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** local_T_prev;
	local_T_prev = BB_std_calloc_2d_double(local_T_prev, nt, nl);

	double** x_ip;
	double*  a_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int beta=0; beta<nt; beta++) {
			BBFE_elemmat_set_local_array_scalar(
					local_T_prev[beta],
					fe,
					&(T_prev_st[beta * nx]),
					e);
		}

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
		}

		for(int alpha=0; alpha<nt; alpha++) {
			for(int i=0; i<nl; i++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					for(int beta=0; beta<nt; beta++) {
						double coef = te->E_prev[alpha][beta] / vals->dt;

						if(fabs(coef) < 1.0e-30) continue;

						double T_beta_ip = BBFE_std_mapping_scalar(
								nl,
								local_T_prev[beta],
								basis->N[p]);

						val_ip[p] +=
							coef *
							BBFE_elemmat_convdiff_vec_mass(
								basis->N[p][i],
								T_beta_ip,
								a_ip[p]);
					}
				}

				double integ_val = BBFE_std_integ_calc(
						np,
						val_ip,
						basis->integ_weight,
						Jacobian_ip);

				int row = ST_window_gid(0, alpha, fe->conn[e][i], nt, nx);
				monolis->mat.R.B[row] += integ_val;
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(local_T_prev, nt, nl);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
}


static void add_heat_source_vector_ST_dG_window(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		double                 t_window_start,
		const int             num_window_slabs)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;
	int nt = te->n_dof;
	int nx = fe->total_num_nodes;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip,      np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;
	double** v_ip;
	double*  a_ip;
	double*  k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	for(int e=0; e<fe->total_num_elems; e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
		}

		for(int slab=0; slab<num_window_slabs; slab++) {
			double t_old = t_window_start + ((double)slab) * vals->dt;
			double t_new = t_old + vals->dt;

			for(int alpha=0; alpha<nt; alpha++) {
				for(int i=0; i<nl; i++) {
					for(int p=0; p<np; p++) {
						val_ip[p] = 0.0;

						for(int q=0; q<te->num_time_ip; q++) {
							double tau = te->tau[q];
							double tq  = t_old + tau * (t_new - t_old);
							double coef = te->weight[q] * te->psi[alpha][q];

							if(fabs(coef) < 1.0e-30) continue;

							double f_ip = manusol_get_source(
									x_ip[p],
									tq,
									a_ip[p],
									v_ip[p],
									k_ip[p]);

							val_ip[p] +=
								coef *
								BBFE_elemmat_convdiff_vec_source(
									basis->N[p][i],
									f_ip);
						}
					}

					double integ_val = BBFE_std_integ_calc(
							np,
							val_ip,
							basis->integ_weight,
							Jacobian_ip);

					int row = ST_window_gid(slab, alpha, fe->conn[e][i], nt, nx);
					monolis->mat.R.B[row] += integ_val;
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}


void set_element_vec_ST_dG_window(
		MONOLIS*              monolis,
		BBFE_DATA*            fe,
		BBFE_BASIS*           basis,
		VALUES*               vals,
		const ST_TIME_ELEMENT* te,
		const double*          T_prev_st,
		double                 t_window_start,
		const int             num_window_slabs)
{
	if(num_window_slabs <= 0) {
		printf("ERROR: num_window_slabs must be positive. got %d\n", num_window_slabs);
		exit(1);
	}

	/* Initial inflow from the previous window/slab enters only the first slab. */
	add_heat_prev_mass_vector_ST_dG_window_first_slab(
			monolis,
			fe,
			basis,
			vals,
			te,
			T_prev_st);

	/* Source terms over all slabs in the window. */
	add_heat_source_vector_ST_dG_window(
			monolis,
			fe,
			basis,
			vals,
			te,
			t_window_start,
			num_window_slabs);
}


void ST_Window_extract_last_slab_state(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const double*          T_window,
		double*                T_last_st)
{
	if(num_window_slabs <= 0) {
		printf("ERROR: num_window_slabs must be positive. got %d\n", num_window_slabs);
		exit(1);
	}

	int nt = te->n_dof;
	int last = num_window_slabs - 1;

	for(int alpha=0; alpha<nt; alpha++) {
		for(int i=0; i<num_space_nodes; i++) {
			T_last_st[alpha * num_space_nodes + i] =
				T_window[ST_window_gid(last, alpha, i, nt, num_space_nodes)];
		}
	}
}


void ST_Window_extract_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const double*          T_window,
		double*                T_nodal)
{
	if(num_window_slabs <= 0) {
		printf("ERROR: num_window_slabs must be positive. got %d\n", num_window_slabs);
		exit(1);
	}

	int nt = te->n_dof;
	int last = num_window_slabs - 1;

	for(int i=0; i<num_space_nodes; i++) {
		T_nodal[i] = 0.0;
		for(int alpha=0; alpha<nt; alpha++) {
			T_nodal[i] +=
				te->e_plus[alpha] *
				T_window[ST_window_gid(last, alpha, i, nt, num_space_nodes)];
		}
	}
}


void set_element_mat_ST_dG_window_degree(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const int    degree,
		const int    num_window_slabs)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG(&te, degree);

	set_element_mat_ST_dG_window(
			monolis,
			fe,
			basis,
			vals,
			&te,
			num_window_slabs);
}


void set_element_vec_ST_dG_window_degree(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const int    degree,
		const double* T_prev_st,
		double       t_window_start,
		const int    num_window_slabs)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG(&te, degree);

	set_element_vec_ST_dG_window(
			monolis,
			fe,
			basis,
			vals,
			&te,
			T_prev_st,
			t_window_start,
			num_window_slabs);
}


void set_element_mat_ST_dG0_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const int    num_window_slabs)
{
	set_element_mat_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			0,
			num_window_slabs);
}


void set_element_vec_ST_dG0_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_window_start,
		const int    num_window_slabs)
{
	set_element_vec_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			0,
			T_prev_st,
			t_window_start,
			num_window_slabs);
}


void set_element_mat_ST_dG1_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const int    num_window_slabs)
{
	set_element_mat_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			1,
			num_window_slabs);
}


void set_element_vec_ST_dG1_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_window_start,
		const int    num_window_slabs)
{
	set_element_vec_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			1,
			T_prev_st,
			t_window_start,
			num_window_slabs);
}


void set_element_mat_ST_dG2_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const int    num_window_slabs)
{
	set_element_mat_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			2,
			num_window_slabs);
}


void set_element_vec_ST_dG2_window(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		const double* T_prev_st,
		double       t_window_start,
		const int    num_window_slabs)
{
	set_element_vec_ST_dG_window_degree(
			monolis,
			fe,
			basis,
			vals,
			2,
			T_prev_st,
			t_window_start,
			num_window_slabs);
}
