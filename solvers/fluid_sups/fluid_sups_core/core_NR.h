
#pragma once

#include "fluid_core.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

void set_element_mat_NR_linear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_NR_linear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_mat_NR_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_NR_linear_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_NR_linear_nonlinear2(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_NR_nonlinear(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void HROM_ddecm_set_residuals_NR(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_residuals_NR_blas(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t);

void HROM_set_element_mat_NR(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS* 	basis,
    BBFE_BC*     	bc,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_modes,
    const int 		num_subdomains,
    const double    dt);

void HROM_set_element_vec_NR(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    VALUES*         vals,
    BBFE_BASIS*	 	basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt,
    double       	t);

void HROM_set_element_mat_NR2(
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

void set_element_mat_NR_mass(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void ROM_std_hlpod_reduced_rhs_to_monollis_linear2(
    MONOLIS*		monolis,
    MONOLIS*		monolis_mass,
    MONOLIS_COM*    monolis_com,
    HLPOD_MAT*      hlpod_mat,
    double*         mode_coeff,
    double*         mode_coeff_old,
    const int       num_2nd_subdomains);

void HROM_ddecm_set_residuals_NR_PSPG(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_residuals_NR_blas2(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_residuals_NR_vec(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,                       //1 + monolis_com->recv_n_neib
    const double    dt,
    double          t);

void HROM_ddecm_set_residuals_NR_blas2_decoupled(
    BBFE_DATA*      fe,
    BBFE_BASIS*     basis,
    VALUES*         vals,
    BBFE_BC*        bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int       num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,
    const double    dt,
    double          t);

void HROM_set_element_mat_NR_decoupled_v(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*   hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void HROM_set_element_mat_NR_decoupled_p(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*   hlpod_vals,
    	HLPOD_MAT*      hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);


void HROM_set_element_vec_NR_p(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);

void HROM_set_element_vec_NR_v(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		VALUES*         vals,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*      hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);