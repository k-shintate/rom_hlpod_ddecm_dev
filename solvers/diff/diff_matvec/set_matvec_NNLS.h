#pragma once

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_RH_for_NNLS_para(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_residuals_for_NNLS_para(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    BBFE_BC*     	bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_neib,			//1 + monolis_com->recv_n_neib
    const double    dt,
    double       	t);

void HROM_ddecm_set_RH_for_NNLS(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    BBFE_BC*     	bc,
    HLPOD_MAT*      hlpod_mat,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int 		num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);

//残差ベクトルのみをNNLSに使う
void HROM_ddecm_set_residuals_for_NNLS(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);

void ddhr_set_matvec_RH_for_NNLS2(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_subdomains,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);

void hr_set_matvec_for_NNLS(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_HR*       hlpod_hr,
    const int       index_snap,
    const int       num_modes,
    const double    dt,
    double       	t);

void HROM_ecm_set_RH_for_NNLS(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_HR*       hlpod_hr,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);

void HROM_ecm_set_residuals_for_NNLS(
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    BBFE_BC*     	bc,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_VALUES*     hlpod_vals,
    HLPOD_HR*       hlpod_hr,
    const int       index_snap,
    const int       num_snapshot,
    const int       num_modes,
    const double    dt,
    double       	t);