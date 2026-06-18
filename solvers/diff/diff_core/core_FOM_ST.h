
#ifndef ST_TIME_ELEMENT_H
#define ST_TIME_ELEMENT_H

#define ST_TIME_FAMILY_DG 0
#define ST_TIME_FAMILY_CG 1

#define ST_TIME_MAX_DOF 8

typedef struct {
	int family;
	int degree;
	int n_dof;

	/* reference time matrices on tau in [0,1] */
	double Mt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double Dt[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];

	/* endpoint values */
	double e_minus[ST_TIME_MAX_DOF];
	double e_plus [ST_TIME_MAX_DOF];

	/* dG time flux matrices */
	double F_self[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];
	double F_prev[ST_TIME_MAX_DOF][ST_TIME_MAX_DOF];

} ST_TIME_ELEMENT;


static void ST_TimeElement_clear(ST_TIME_ELEMENT* te)
{
	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof = 0;

	for(int i=0; i<ST_TIME_MAX_DOF; i++) {
		te->e_minus[i] = 0.0;
		te->e_plus [i] = 0.0;

		for(int j=0; j<ST_TIME_MAX_DOF; j++) {
			te->Mt[i][j] = 0.0;
			te->Dt[i][j] = 0.0;
			te->F_self[i][j] = 0.0;
			te->F_prev[i][j] = 0.0;
		}
	}
}


/* dG(0): psi_0(tau)=1 */
static void ST_TimeElement_init_dG0(ST_TIME_ELEMENT* te)
{
	ST_TimeElement_clear(te);

	te->family = ST_TIME_FAMILY_DG;
	te->degree = 0;
	te->n_dof = 1;

	te->Mt[0][0] = 1.0;   /* int_0^1 psi psi dtau */
	te->Dt[0][0] = 0.0;   /* int_0^1 psi_t psi dtau */

	te->e_minus[0] = 1.0;
	te->e_plus [0] = 1.0;

	/*
	 * dG time-upwind flux:
	 * self term:  + M T^{n+1}
	 * prev term:  - M T^n
	 */
	te->F_self[0][0] = 1.0;
	te->F_prev[0][0] = 1.0;
}

#endif

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

void set_element_mat_ST_dG0(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG0(&te);

	/*
	 * dG(0):
	 *
	 * A = (1/dt) * F_self * M + Mt * K
	 *
	 * F_self = 1
	 * Mt     = 1
	 */
	double coef_mass = te.F_self[0][0] / vals->dt;
	double coef_diff = te.Mt[0][0];

	add_heat_mass_matrix(
		monolis,
		fe,
		basis,
		vals,
		coef_mass);

	add_heat_diffusion_matrix(
		monolis,
		fe,
		basis,
		vals,
		coef_diff);
}


void set_element_vec_ST_dG0(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals,
		double       t_new)
{
	ST_TIME_ELEMENT te;
	ST_TimeElement_init_dG0(&te);

	/*
	 * dG(0):
	 *
	 * b = (1/dt) * F_prev * M * T_old + Mt * F_new
	 *
	 * F_prev = 1
	 * Mt     = 1
	 */
	double coef_old_mass = te.F_prev[0][0] / vals->dt;
	double coef_source   = te.Mt[0][0];

	add_heat_old_mass_vector(
		monolis,
		fe,
		basis,
		vals,
		coef_old_mass);

	add_heat_source_vector(
		monolis,
		fe,
		basis,
		vals,
		t_new,
		coef_source);
}