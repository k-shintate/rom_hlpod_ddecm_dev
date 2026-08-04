#include "core_FOM_ST.h"

static void ST_TimeElement_clear(
		ST_TIME_ELEMENT* te)
{
	if(te == NULL) {
		fprintf(
				stderr,
				"ERROR: te is NULL in ST_TimeElement_clear.\n");
		exit(EXIT_FAILURE);
	}

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof = 0;

	te->is_nodal = 0;
	te->num_time_ip = 0;

	for(int alpha=0; alpha<ST_TIME_MAX_DOF; alpha++) {
		te->dof_tau[alpha] = 0.0;

		te->e_minus[alpha] = 0.0;
		te->e_plus [alpha] = 0.0;

		for(int q=0; q<ST_TIME_MAX_IP; q++) {
			te->psi[alpha][q] = 0.0;
		}

		for(int beta=0; beta<ST_TIME_MAX_DOF; beta++) {
			te->Mt[alpha][beta] = 0.0;
			te->Dt[alpha][beta] = 0.0;

			te->E_self[alpha][beta] = 0.0;
			te->E_prev [alpha][beta] = 0.0;
		}
	}

	for(int q=0; q<ST_TIME_MAX_IP; q++) {
		te->tau[q] = 0.0;
		te->weight[q] = 0.0;
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

static void ST_TimeElement_set_flux_matrices(
		ST_TIME_ELEMENT* te)
{
	for(int alpha=0; alpha<te->n_dof; alpha++) {
		for(int beta=0; beta<te->n_dof; beta++) {
			te->E_self[alpha][beta] =
				te->e_minus[alpha] *
				te->e_minus[beta];

			te->E_prev[alpha][beta] =
				te->e_minus[alpha] *
				te->e_plus[beta];
		}
	}
}


static void ST_TimeElement_init_dG0(
		ST_TIME_ELEMENT* te)
{
	ST_TimeElement_clear(te);

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof = 1;

	/*
	 * 現在のdG(0)は後退Euler互換とし、
	 * 唯一の時間係数をスラブ右端値として扱う。
	 */
	te->is_nodal = 1;
	te->dof_tau[0] = 1.0;

	te->Mt[0][0] = 1.0;
	te->Dt[0][0] = 0.0;

	te->e_minus[0] = 1.0;
	te->e_plus [0] = 1.0;

	ST_TimeElement_set_flux_matrices(te);

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

	/*
	 * dG(1)の節点型時間基底
	 *
	 * alpha = 0 : スラブ左端 tau = 0
	 * alpha = 1 : スラブ右端 tau = 1
	 */
	te->is_nodal = 1;
	te->dof_tau[0] = 0.0;
	te->dof_tau[1] = 1.0;

	te->Mt[0][0] = 1.0 / 3.0;
	te->Mt[0][1] = 1.0 / 6.0;
	te->Mt[1][0] = 1.0 / 6.0;
	te->Mt[1][1] = 1.0 / 3.0;

	te->Dt[0][0] = -0.5;
	te->Dt[0][1] =  0.5;
	te->Dt[1][0] = -0.5;
	te->Dt[1][1] =  0.5;

	te->e_minus[0] = 1.0;
	te->e_minus[1] = 0.0;

	te->e_plus[0] = 0.0;
	te->e_plus[1] = 1.0;

	ST_TimeElement_set_flux_matrices(te);

	const double s =
		0.57735026918962576451;

	te->num_time_ip = 2;

	te->tau[0] =
		0.5 * (1.0 - s);

	te->tau[1] =
		0.5 * (1.0 + s);

	te->weight[0] = 0.5;
	te->weight[1] = 0.5;

	for(int q=0; q<te->num_time_ip; q++) {
		const double tau =
			te->tau[q];

		te->psi[0][q] =
			1.0 - tau;

		te->psi[1][q] =
			tau;
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
	te->n_dof = 3;

	te->is_nodal = 1;

	te->dof_tau[0] = 0.0;
	te->dof_tau[1] = 0.5;
	te->dof_tau[2] = 1.0;

	te->Mt[0][0] =  2.0 / 15.0;
	te->Mt[0][1] =  1.0 / 15.0;
	te->Mt[0][2] = -1.0 / 30.0;

	te->Mt[1][0] =  1.0 / 15.0;
	te->Mt[1][1] =  8.0 / 15.0;
	te->Mt[1][2] =  1.0 / 15.0;

	te->Mt[2][0] = -1.0 / 30.0;
	te->Mt[2][1] =  1.0 / 15.0;
	te->Mt[2][2] =  2.0 / 15.0;

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

	te->e_plus[0] = 0.0;
	te->e_plus[1] = 0.0;
	te->e_plus[2] = 1.0;

	ST_TimeElement_set_flux_matrices(te);

	const double s =
		0.77459666924148337704;

	te->num_time_ip = 3;

	te->tau[0] =
		0.5 * (1.0 - s);

	te->tau[1] =
		0.5;

	te->tau[2] =
		0.5 * (1.0 + s);

	te->weight[0] = 5.0 / 18.0;
	te->weight[1] = 4.0 / 9.0;
	te->weight[2] = 5.0 / 18.0;

	for(int q=0; q<te->num_time_ip; q++) {
		const double tau =
			te->tau[q];

		te->psi[0][q] =
			2.0 * tau * tau
			- 3.0 * tau
			+ 1.0;

		te->psi[1][q] =
			-4.0 * tau * tau
			+ 4.0 * tau;

		te->psi[2][q] =
			2.0 * tau * tau
			- tau;
	}
}

static void check_ST_Dirichlet_values(
        const FE_SYSTEM*       sys,
        const ST_TIME_ELEMENT* te,
        const int              num_window_slabs,
        const double           t_window_start,
        const double*          T_window)
{
    const int nx =
        sys->fe.total_num_nodes;

    const int nt =
        te->n_dof;

    for(int slab = 0;
            slab < num_window_slabs;
            slab++) {

        const double t_s =
            t_window_start
            + (double)slab
            * sys->vals.dt;

        for(int alpha = 0;
                alpha < nt;
                alpha++) {

            const double t_alpha =
                t_s
                + te->dof_tau[alpha]
                * sys->vals.dt;

            double max_error = 0.0;
            int max_node = -1;
            int num_bc_nodes = 0;

            for(int i = 0;
                    i < nx;
                    i++) {

                if(!sys->bc.D_bc_exists[i]) {
                    continue;
                }

                const int st_node =
                    ST_window_gid(
                            slab,
                            alpha,
                            i,
                            nt,
                            nx);

                const double exact =
                    manusol_get_sol(
                            sys->fe.x[i][0],
                            sys->fe.x[i][1],
                            sys->fe.x[i][2],
                            t_alpha);

                const double error =
                    fabs(
                            T_window[st_node]
                            -
                            exact);

                num_bc_nodes++;

                if(error > max_error) {
                    max_error = error;
                    max_node = i;
                }
            }

            printf(
                    "BC check: slab=%d alpha=%d "
                    "time=%.16e nodes=%d "
                    "max_error=%.16e node=%d\n",
                    slab,
                    alpha,
                    t_alpha,
                    num_bc_nodes,
                    max_error,
                    max_node);
        }
    }
}

void ST_TimeElement_init_dG(
		ST_TIME_ELEMENT* te,
		const int        degree)
{
	if(te == NULL) {
		fprintf(
				stderr,
				"ERROR: te is NULL in ST_TimeElement_init_dG.\n");
		exit(EXIT_FAILURE);
	}

	switch(degree) {
	case 0:
		ST_TimeElement_init_dG0(te);
		break;

	case 1:
		ST_TimeElement_init_dG1(te);
		break;

	case 2:
		ST_TimeElement_init_dG2(te);
		break;

	default:
		fprintf(
				stderr,
				"ERROR: unsupported dG degree %d.\n",
				degree);
		exit(EXIT_FAILURE);
	}
}

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
	if(te == NULL ||
	   T_nodal == NULL ||
	   T_prev_st == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"ST_TimeElement_set_prev_trace_from_nodal_value.\n");
		exit(EXIT_FAILURE);
	}

	if(num_space_nodes <= 0 ||
	   te->n_dof <= 0) {
		fprintf(
				stderr,
				"ERROR: invalid dimensions in "
				"ST_TimeElement_set_prev_trace_from_nodal_value.\n");
		exit(EXIT_FAILURE);
	}

	/*
	 * 次式を満たす時間係数を構築する。
	 *
	 * sum_alpha e_plus[alpha] T_prev_st[alpha,i]
	 *     = T_nodal[i]
	 *
	 * 現在の節点型基底では、実質的に最後の時間自由度へ
	 * T_nodalを設定することになる。
	 */
	double endpoint_norm2 = 0.0;

	for(int alpha=0; alpha<te->n_dof; alpha++) {
		endpoint_norm2 +=
			te->e_plus[alpha] *
			te->e_plus[alpha];
	}

	if(endpoint_norm2 <= 1.0e-30) {
		fprintf(
				stderr,
				"ERROR: zero right-end temporal trace.\n");
		exit(EXIT_FAILURE);
	}

	for(int alpha=0; alpha<te->n_dof; alpha++) {
		const double temporal_coef =
			te->e_plus[alpha] /
			endpoint_norm2;

		for(int i=0; i<num_space_nodes; i++) {
			T_prev_st[
				alpha * num_space_nodes + i
			] =
				temporal_coef *
				T_nodal[i];
		}
	}
}


void ST_Window_extract_slab_right_trace_to_nodal_value(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const int              slab,
		const double*          T_window,
		double*                T_nodal)
{
	if(te == NULL ||
	   T_window == NULL ||
	   T_nodal == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"ST_Window_extract_slab_right_trace_to_nodal_value.\n");
		exit(EXIT_FAILURE);
	}

	if(num_space_nodes <= 0 ||
	   num_window_slabs <= 0 ||
	   slab < 0 ||
	   slab >= num_window_slabs) {
		fprintf(
				stderr,
				"ERROR: invalid dimensions or slab ID in "
				"ST_Window_extract_slab_right_trace_to_nodal_value.\n");
		exit(EXIT_FAILURE);
	}

	const int nt =
		te->n_dof;

	const int nx =
		num_space_nodes;

	for(int i=0; i<nx; i++) {
		T_nodal[i] = 0.0;

		for(int alpha=0; alpha<nt; alpha++) {
			const int gid =
				ST_window_gid(
						slab,
						alpha,
						i,
						nt,
						nx);

			T_nodal[i] +=
				te->e_plus[alpha] *
				T_window[gid];
		}
	}
}


void ST_Window_extract_last_slab_state(
		const ST_TIME_ELEMENT* te,
		const int              num_space_nodes,
		const int              num_window_slabs,
		const double*          T_window,
		double*                T_last_st)
{
	if(te == NULL ||
	   T_window == NULL ||
	   T_last_st == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"ST_Window_extract_last_slab_state.\n");
		exit(EXIT_FAILURE);
	}

	if(num_space_nodes <= 0 ||
	   num_window_slabs <= 0) {
		fprintf(
				stderr,
				"ERROR: invalid dimensions in "
				"ST_Window_extract_last_slab_state.\n");
		exit(EXIT_FAILURE);
	}

	const int nt =
		te->n_dof;

	const int nx =
		num_space_nodes;

	const int last_slab =
		num_window_slabs - 1;

	for(int alpha=0; alpha<nt; alpha++) {
		for(int i=0; i<nx; i++) {
			T_last_st[alpha * nx + i] =
				T_window[
					ST_window_gid(
							last_slab,
							alpha,
							i,
							nt,
							nx)
				];
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
	ST_Window_extract_slab_right_trace_to_nodal_value(
			te,
			num_space_nodes,
			num_window_slabs,
			num_window_slabs - 1,
			T_window,
			T_nodal);
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

int ST_window_gid(
		const int slab_id,
		const int time_dof,
		const int space_id,
		const int num_time_dofs,
		const int num_space_nodes)
{
	return
		(slab_id * num_time_dofs + time_dof)
		* num_space_nodes
		+ space_id;
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



#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

static int ST_graph_gid(
    const int slab,
    const int alpha,
    const int space_id,
    const int nt,
    const int nx)
{
    return (slab * nt + alpha) * nx + space_id;
}



static int ST_temporal_same_slab_is_nonzero(
		const ST_TIME_ELEMENT* te,
		const int              alpha,
		const int              beta,
		const double           eps)
{
	return
		fabs(
			te->Dt[alpha][beta]
			+ te->E_self[alpha][beta]
		) > eps
		||
		fabs(te->Mt[alpha][beta]) > eps;
}


static int ST_temporal_prev_slab_is_nonzero(
		const ST_TIME_ELEMENT* te,
		const int              alpha,
		const int              beta,
		const double           eps)
{
	return
		fabs(te->E_prev[alpha][beta]) > eps;
}


void set_nonzero_pattern_space_time(
		MONOLIS*               monolis,
		const char*            label,
		const char*            directory,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs)
{
	if(monolis == NULL ||
	   label == NULL ||
	   directory == NULL ||
	   te == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"set_nonzero_pattern_space_time.\n");
		exit(EXIT_FAILURE);
	}

	if(te->n_dof <= 0 ||
	   te->n_dof > ST_TIME_MAX_DOF) {
		fprintf(
				stderr,
				"ERROR: invalid temporal DOF count %d.\n",
				te->n_dof);
		exit(EXIT_FAILURE);
	}

	if(num_window_slabs <= 0) {
		fprintf(
				stderr,
				"ERROR: num_window_slabs must be positive.\n");
		exit(EXIT_FAILURE);
	}

	const double eps =
		1.0e-30;

	const int nt =
		te->n_dof;

	const char* fname =
		monolis_get_global_input_file_name(
				MONOLIS_DEFAULT_TOP_DIR,
				MONOLIS_DEFAULT_PART_DIR,
				label);

	FILE* fp;

	fp =
		BBFE_sys_read_fopen(
				fp,
				fname,
				directory);

	int num_nodes = 0;

	if(fscanf(fp, "%d", &num_nodes) != 1 ||
	   num_nodes <= 0) {
		fprintf(
				stderr,
				"ERROR: failed to read spatial node count.\n");
		fclose(fp);
		exit(EXIT_FAILURE);
	}

	int64_t num_space_items_64 = 0;

	for(int i=0; i<num_nodes; i++) {
		int node_id;
		int num_adj_nodes;

		if(fscanf(fp, "%d", &node_id) != 1 ||
		   fscanf(fp, "%d", &num_adj_nodes) != 1 ||
		   num_adj_nodes < 0) {
			fprintf(
					stderr,
					"ERROR: invalid graph row %d.\n",
					i);
			fclose(fp);
			exit(EXIT_FAILURE);
		}

		for(int j=0; j<num_adj_nodes; j++) {
			int adjacent_node;

			if(fscanf(fp, "%d", &adjacent_node) != 1) {
				fprintf(
						stderr,
						"ERROR: failed to read adjacency item.\n");
				fclose(fp);
				exit(EXIT_FAILURE);
			}
		}

		num_space_items_64 +=
			num_adj_nodes;
	}

	fclose(fp);

	if(num_space_items_64 > INT_MAX) {
		fprintf(
				stderr,
				"ERROR: spatial graph is too large.\n");
		exit(EXIT_FAILURE);
	}

	const int num_space_items =
		(int)num_space_items_64;

	int* space_index =
		BB_std_calloc_1d_int(
				space_index,
				num_nodes + 1);

	int* space_item =
		BB_std_calloc_1d_int(
				space_item,
				num_space_items);

	fp =
		BBFE_sys_read_fopen(
				fp,
				fname,
				directory);

	int num_nodes_second = 0;

	if(fscanf(fp, "%d", &num_nodes_second) != 1 ||
	   num_nodes_second != num_nodes) {
		fprintf(
				stderr,
				"ERROR: inconsistent spatial node count.\n");
		fclose(fp);
		exit(EXIT_FAILURE);
	}

	int item_position = 0;

	space_index[0] = 0;

	for(int i=0; i<num_nodes; i++) {
		int node_id;
		int num_adj_nodes;

		if(fscanf(fp, "%d", &node_id) != 1 ||
		   fscanf(fp, "%d", &num_adj_nodes) != 1 ||
		   num_adj_nodes < 0) {
			fprintf(
					stderr,
					"ERROR: invalid graph row %d.\n",
					i);
			fclose(fp);
			exit(EXIT_FAILURE);
		}

		for(int j=0; j<num_adj_nodes; j++) {
			int adjacent_node;

			if(fscanf(fp, "%d", &adjacent_node) != 1) {
				fprintf(
						stderr,
						"ERROR: failed to read adjacency item.\n");
				fclose(fp);
				exit(EXIT_FAILURE);
			}

			if(adjacent_node < 0 ||
			   adjacent_node >= num_nodes) {
				fprintf(
						stderr,
						"ERROR: invalid adjacent node %d.\n",
						adjacent_node);
				fclose(fp);
				exit(EXIT_FAILURE);
			}

			space_item[item_position++] =
				adjacent_node;
		}

		space_index[i + 1] =
			item_position;
	}

	fclose(fp);

	int* has_self =
		BB_std_calloc_1d_int(
				has_self,
				num_nodes);

	int64_t sum_spatial_row_sizes = 0;

	for(int i=0; i<num_nodes; i++) {
		has_self[i] = 0;

		for(int p=space_index[i];
				p<space_index[i + 1];
				p++) {
			if(space_item[p] == i) {
				has_self[i] = 1;
				break;
			}
		}

		int row_size =
			space_index[i + 1]
			- space_index[i];

		if(!has_self[i]) {
			row_size++;
		}

		sum_spatial_row_sizes +=
			row_size;
	}

	int num_same_pairs = 0;
	int num_prev_pairs = 0;

	for(int alpha=0; alpha<nt; alpha++) {
		for(int beta=0; beta<nt; beta++) {
			if(ST_temporal_same_slab_is_nonzero(
					te,
					alpha,
					beta,
					eps)) {
				num_same_pairs++;
			}

			if(ST_temporal_prev_slab_is_nonzero(
					te,
					alpha,
					beta,
					eps)) {
				num_prev_pairs++;
			}
		}
	}

	const int64_t num_st_nodes_64 =
		(int64_t)num_window_slabs
		* (int64_t)nt
		* (int64_t)num_nodes;

	const int64_t num_st_items_64 =
		(int64_t)num_window_slabs
		* (int64_t)num_same_pairs
		* sum_spatial_row_sizes
		+
		(int64_t)(num_window_slabs - 1)
		* (int64_t)num_prev_pairs
		* sum_spatial_row_sizes;

	if(num_st_nodes_64 > INT_MAX ||
	   num_st_items_64 > INT_MAX) {
		fprintf(
				stderr,
				"ERROR: ST graph exceeds int range.\n");
		exit(EXIT_FAILURE);
	}

	const int num_st_nodes =
		(int)num_st_nodes_64;

	const int num_st_items =
		(int)num_st_items_64;

	int* st_index =
		BB_std_calloc_1d_int(
				st_index,
				num_st_nodes + 1);

	int* st_item =
		BB_std_calloc_1d_int(
				st_item,
				num_st_items);

	int st_position = 0;

	st_index[0] = 0;

	for(int slab=0;
			slab<num_window_slabs;
			slab++) {

		for(int alpha=0;
				alpha<nt;
				alpha++) {

			for(int i=0;
					i<num_nodes;
					i++) {

				const int row =
					ST_window_gid(
							slab,
							alpha,
							i,
							nt,
							num_nodes);

				/*
				 * 同一スラブ:
				 * (slab,alpha,i) -> (slab,beta,j)
				 */
				for(int beta=0;
						beta<nt;
						beta++) {

					if(!ST_temporal_same_slab_is_nonzero(
							te,
							alpha,
							beta,
							eps)) {
						continue;
					}

					/*
					 * 同一空間節点を必ず追加する。
					 */
					st_item[st_position++] =
						ST_window_gid(
								slab,
								beta,
								i,
								nt,
								num_nodes);

					for(int p=space_index[i];
							p<space_index[i + 1];
							p++) {

						const int j =
							space_item[p];

						if(j == i) {
							continue;
						}

						st_item[st_position++] =
							ST_window_gid(
									slab,
									beta,
									j,
									nt,
									num_nodes);
					}
				}

				/*
				 * 前スラブ:
				 * (slab,alpha,i) -> (slab-1,beta,j)
				 */
				if(slab > 0) {
					for(int beta=0;
							beta<nt;
							beta++) {

						if(!ST_temporal_prev_slab_is_nonzero(
								te,
								alpha,
								beta,
								eps)) {
							continue;
						}

						st_item[st_position++] =
							ST_window_gid(
									slab - 1,
									beta,
									i,
									nt,
									num_nodes);

						for(int p=space_index[i];
								p<space_index[i + 1];
								p++) {

							const int j =
								space_item[p];

							if(j == i) {
								continue;
							}

							st_item[st_position++] =
								ST_window_gid(
										slab - 1,
										beta,
										j,
										nt,
										num_nodes);
						}
					}
				}

				st_index[row + 1] =
					st_position;
			}
		}
	}

	if(st_position != num_st_items) {
		fprintf(
				stderr,
				"ERROR: inconsistent ST graph size: "
				"expected %d, generated %d.\n",
				num_st_items,
				st_position);
		exit(EXIT_FAILURE);
	}

	monolis_get_nonzero_pattern_by_nodal_graph_R(
			monolis,
			num_st_nodes,
			1,
			st_index,
			st_item);

	BB_std_free_1d_int(
			space_index,
			num_nodes + 1);

	BB_std_free_1d_int(
			space_item,
			num_space_items);

	BB_std_free_1d_int(
			has_self,
			num_nodes);

	BB_std_free_1d_int(
			st_index,
			num_st_nodes + 1);

	BB_std_free_1d_int(
			st_item,
			num_st_items);
}




void BBFE_sys_monowrap_set_Dirichlet_bc_space_time(
    MONOLIS*                     monolis,
    BBFE_DATA*                   fe,
    const ST_TIME_ELEMENT*       te,
    const int                    num_window_slabs,
    const int                    num_components,
    const BBFE_BC*               bc,
    const double                 t_window_start,
    const double                 dt,
    ST_DIRICHLET_VALUE_FUNC      value_func,
    void*                        user_data,
    double*                      g_rhs)
{
    if(monolis == NULL ||
       fe      == NULL ||
       te      == NULL ||
       bc      == NULL ||
       g_rhs   == NULL) {
        fprintf(
            stderr,
            "ERROR: NULL pointer in "
            "BBFE_sys_monowrap_set_Dirichlet_bc_space_time.\n");
        exit(EXIT_FAILURE);
    }

    if(num_window_slabs <= 0) {
        fprintf(
            stderr,
            "ERROR: num_window_slabs must be positive. got %d\n",
            num_window_slabs);
        exit(EXIT_FAILURE);
    }

    if(num_components <= 0) {
        fprintf(
            stderr,
            "ERROR: num_components must be positive. got %d\n",
            num_components);
        exit(EXIT_FAILURE);
    }

    if(dt <= 0.0) {
        fprintf(
            stderr,
            "ERROR: dt must be positive. got %.16e\n",
            dt);
        exit(EXIT_FAILURE);
    }

    if(te->n_dof <= 0) {
        fprintf(
            stderr,
            "ERROR: invalid number of temporal DOFs: %d\n",
            te->n_dof);
        exit(EXIT_FAILURE);
    }

    /*
     * This implementation evaluates the boundary value at each
     * temporal interpolation node.
     */
    if(!te->is_nodal) {
        fprintf(
            stderr,
            "ERROR: this Dirichlet routine requires a nodal "
            "temporal basis.\n");
        exit(EXIT_FAILURE);
    }

    const int nx = fe->total_num_nodes;
    const int nt = te->n_dof;

    /*
     * s:
     *   Time-slab index.
     */
    for(int s = 0; s < num_window_slabs; s++) {

        /*
         * Start time of slab s:
         *
         *   t_s = t_window_start + s * dt
         */
        const double t_s =
            t_window_start + ((double)s) * dt;

        /*
         * alpha:
         *   Temporal DOF index in slab s.
         */
        for(int alpha = 0; alpha < nt; alpha++) {

            /*
             * Physical time associated with temporal DOF alpha:
             *
             *   t_{s,alpha}
             *     = t_s + dof_tau[alpha] * dt
             */
            const double t_alpha =
                t_s + te->dof_tau[alpha] * dt;

            /*
             * i:
             *   Spatial-node index.
             */
            for(int i = 0; i < nx; i++) {

                /*
                 * k:
                 *   Physical component index at spatial node i.
                 *
                 * For the current scalar heat problem:
                 *   num_components = 1
                 *   k = 0
                 */
                for(int k = 0; k < num_components; k++) {

                    const int spatial_bc_id =
                        num_components * i + k;

                    if(!bc->D_bc_exists[spatial_bc_id]) {
                        continue;
                    }

                    /*
                     * Flattened space-time node:
                     *
                     *   (s, alpha, i)
                     *       -> st_node
                     */
                    const int st_node =
                        ST_window_gid(
                            s,
                            alpha,
                            i,
                            nt,
                            nx);

                    double imposed_value;

                    /*
                     * Time-dependent boundary value.
                     */
                    if(value_func != NULL) {
                        imposed_value =
                            value_func(
                                fe->x[i],
                                k,
                                t_alpha,
                                user_data);
                    }
                    /*
                     * Constant or previously assigned boundary value.
                     */
                    else {
                        imposed_value =
                            bc->imposed_D_val[spatial_bc_id];
                    }

                    /*
                     * st_node is the MONOLIS node index.
                     * k is the component index.
                     */
                    monolis_set_Dirichlet_bc_R(
                        monolis,
                        g_rhs,
                        st_node,
                        k,
                        imposed_value);
                }
            }
        }
    }
}

static double ST_TimeElement_get_dof_tau(
    const ST_TIME_ELEMENT* te,
    const int              alpha)
{
    if(te == NULL) {
        fprintf(stderr, "ERROR: te is NULL.\n");
        exit(EXIT_FAILURE);
    }

    if(alpha < 0 || alpha >= te->n_dof) {
        fprintf(
            stderr,
            "ERROR: invalid temporal DOF alpha = %d.\n",
            alpha);
        exit(EXIT_FAILURE);
    }

    /*
     * dG(0):
     *   The current implementation uses the right endpoint,
     *   consistently with backward Euler.
     */
    if(te->degree == 0) {
        return 1.0;
    }

    /*
     * dG(1):
     *   alpha = 0 -> tau = 0
     *   alpha = 1 -> tau = 1
     */
    if(te->degree == 1) {
        return (double)alpha;
    }

    /*
     * dG(2):
     *   alpha = 0 -> tau = 0
     *   alpha = 1 -> tau = 1/2
     *   alpha = 2 -> tau = 1
     */
    if(te->degree == 2) {
        return 0.5 * (double)alpha;
    }

    fprintf(
        stderr,
        "ERROR: unsupported temporal degree %d.\n",
        te->degree);
    exit(EXIT_FAILURE);
}


void BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar(
		MONOLIS*               monolis,
		BBFE_DATA*             fe,
		VALUES*                vals,
		BBFE_BC*               bc,
		const ST_TIME_ELEMENT* te,
		const int              num_window_slabs,
		const double           t_window_start,
		double*                g_rhs)
{
	if(monolis == NULL ||
	   fe == NULL ||
	   vals == NULL ||
	   bc == NULL ||
	   te == NULL ||
	   g_rhs == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar.\n");
		exit(EXIT_FAILURE);
	}

	if(num_window_slabs <= 0) {
		fprintf(
				stderr,
				"ERROR: num_window_slabs must be positive.\n");
		exit(EXIT_FAILURE);
	}

	if(vals->dt <= 0.0) {
		fprintf(
				stderr,
				"ERROR: dt must be positive.\n");
		exit(EXIT_FAILURE);
	}

	if(te->n_dof <= 0 ||
	   !te->is_nodal) {
		fprintf(
				stderr,
				"ERROR: nodal temporal basis is required.\n");
		exit(EXIT_FAILURE);
	}

	const int nx =
		fe->total_num_nodes;

	const int nt =
		te->n_dof;

	for(int slab=0;
			slab<num_window_slabs;
			slab++) {

		const double t_s =
			t_window_start
			+ (double)slab
			* vals->dt;

		for(int alpha=0;
				alpha<nt;
				alpha++) {

			const double t_alpha =
				t_s
				+ te->dof_tau[alpha]
				* vals->dt;

			/*
			 * 時間自由度alphaに対応する物理時刻で
			 * 厳密解と境界値を設定する。
			 */
			manusol_set_theo_sol(
					fe,
					vals->theo_sol,
					t_alpha);

			BBFE_manusol_set_bc_scalar(
					fe,
					bc,
					vals->theo_sol,
					t_alpha);

			for(int i=0; i<nx; i++) {
				if(!bc->D_bc_exists[i]) {
					continue;
				}

				const int st_node =
					ST_window_gid(
							slab,
							alpha,
							i,
							nt,
							nx);

				monolis_set_Dirichlet_bc_R(
						monolis,
						g_rhs,
						st_node,
						0,
						bc->imposed_D_val[i]);
			}
		}
	}
}

static void ST_Window_set_Dirichlet_values_to_solution_scalar(
		BBFE_DATA*              fe,
		VALUES*                 vals,
		BBFE_BC*                bc,
		const ST_TIME_ELEMENT*  te,
		const int               num_window_slabs,
		const double            t_window_start,
		double*                 T_window)
{
	if(fe == NULL ||
	   vals == NULL ||
	   bc == NULL ||
	   te == NULL ||
	   T_window == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in "
				"ST_Window_set_Dirichlet_values_to_solution_scalar.\n");
		exit(EXIT_FAILURE);
	}

	if(num_window_slabs <= 0) {
		fprintf(
				stderr,
				"ERROR: num_window_slabs must be positive.\n");
		exit(EXIT_FAILURE);
	}

	if(vals->dt <= 0.0) {
		fprintf(
				stderr,
				"ERROR: dt must be positive.\n");
		exit(EXIT_FAILURE);
	}

	if(te->n_dof <= 0 ||
	   !te->is_nodal) {
		fprintf(
				stderr,
				"ERROR: nodal temporal basis is required.\n");
		exit(EXIT_FAILURE);
	}

	const int nx =
		fe->total_num_nodes;

	const int nt =
		te->n_dof;

	for(int slab=0;
			slab<num_window_slabs;
			slab++) {

		const double t_s =
			t_window_start
			+
			(double)slab
			*
			vals->dt;

		for(int alpha=0;
				alpha<nt;
				alpha++) {

			const double t_alpha =
				t_s
				+
				te->dof_tau[alpha]
				*
				vals->dt;

			/*
			 * 時間自由度に対応する物理時刻で
			 * Dirichlet境界値を再計算する。
			 */
			manusol_set_theo_sol(
					fe,
					vals->theo_sol,
					t_alpha);

			BBFE_manusol_set_bc_scalar(
					fe,
					bc,
					vals->theo_sol,
					t_alpha);

			for(int i=0;
					i<nx;
					i++) {

				if(!bc->D_bc_exists[i]) {
					continue;
				}

				const int st_node =
					ST_window_gid(
							slab,
							alpha,
							i,
							nt,
							nx);

				T_window[st_node] =
					bc->imposed_D_val[i];
			}
		}
	}
}


void solver_fom_space_time(
		FE_SYSTEM*              sys,
		const ST_TIME_ELEMENT*  te,
		const int               num_window_slabs,
		const double            t_window_start,
		const int               window_step,
		const double*           T_prev_st,
		double*                 T_window,
		double*                 T_last_st)
{
	if(sys == NULL ||
	   te == NULL ||
	   T_prev_st == NULL ||
	   T_window == NULL) {
		fprintf(
				stderr,
				"ERROR: NULL pointer in solver_fom_space_time.\n");
		exit(EXIT_FAILURE);
	}

	if(num_window_slabs <= 0 ||
	   te->n_dof <= 0 ||
	   sys->vals.dt <= 0.0) {
		fprintf(
				stderr,
				"ERROR: invalid ST dimensions or dt.\n");
		exit(EXIT_FAILURE);
	}

	const int nx =
		sys->fe.total_num_nodes;

	const double t_window_end =
		t_window_start
		+
		(double)num_window_slabs
		*
		sys->vals.dt;

	printf(
			"\n%s ----------------- "
			"ST window %d : [%e, %e] ----------------\n",
			CODENAME,
			window_step,
			t_window_start,
			t_window_end);

	/*
	 * 作業行列をクリアする。
	 */
	monolis_clear_mat_value_rhs_R(
			&(sys->monolis));

	/*
	 * Dirichlet処理前の定数ST行列を復元する。
	 */
	monolis_copy_mat_value_R(
			&(sys->monolis0),
			&(sys->monolis));

	/*
	 * 前窓からの流入項とsource項。
	 */
	set_element_vec_ST_dG_window(
			&(sys->monolis),
			&(sys->fe),
			&(sys->basis),
			&(sys->vals),
			te,
			T_prev_st,
			t_window_start,
			num_window_slabs);

	/*
	 * 行列・右辺へDirichlet条件を設定する。
	 */
	BBFE_sys_monowrap_set_Dirichlet_bc_ST_window_scalar(
			&(sys->monolis),
			&(sys->fe),
			&(sys->vals),
			&(sys->bc),
			te,
			num_window_slabs,
			t_window_start,
			sys->monolis.mat.R.B);

	/*
	 * 拘束自由度を初期解ベクトルへ明示的に設定する。
	 *
	 * MONOLISが拘束自由度を更新対象から除外する場合、
	 * ここで設定した値が保持される。
	 */
	ST_Window_set_Dirichlet_values_to_solution_scalar(
			&(sys->fe),
			&(sys->vals),
			&(sys->bc),
			te,
			num_window_slabs,
			t_window_start,
			T_window);

	/*
	 * dG ST行列は一般に非対称。
	 */
	BBFE_sys_monowrap_solve(
			&(sys->monolis),
			&(sys->monolis_com),
			T_window,
			MONOLIS_ITER_BICGSTAB,
			MONOLIS_PREC_DIAG,
			sys->vals.mat_max_iter,
			sys->vals.mat_epsilon);

	/*
	 * 求解後にもDirichlet自由度を明示的に復元する。
	 *
	 * この処理をT_last_stや右端値の抽出より先に行うこと。
	 */
	ST_Window_set_Dirichlet_values_to_solution_scalar(
			&(sys->fe),
			&(sys->vals),
			&(sys->bc),
			te,
			num_window_slabs,
			t_window_start,
			T_window);

	/*
	 * 境界値を含んだ状態を次窓へ渡す。
	 */
	if(T_last_st != NULL) {
		ST_Window_extract_last_slab_state(
				te,
				nx,
				num_window_slabs,
				T_window,
				T_last_st);
	}

	/*
	 * 境界値を含んだ右端物理解を抽出する。
	 */
	ST_Window_extract_right_trace_to_nodal_value(
			te,
			nx,
			num_window_slabs,
			T_window,
			sys->vals.T);
}