#ifndef INC_SVD_H_
#define INC_SVD_H_

#include <stddef.h>
#include "monolis_wrapper_scalapack_c.h"

typedef struct {
  int rows;     /* 行数 */
  int cols;     /* 列数 */
  double *a;    /* a[i + j*rows] = A(i,j) */
} CMat;

typedef struct {
  int K;          /* 行数 */
  int N;          /* 現在の列数 */
  int r;          /* 現在のランク */
  double *U;      /* [K x r] column-major */
  double *V;      /* [N x r] column-major */
  double *sigma;  /* [r] 特異値 */
} IncSVD;

void incsvd_update(IncSVD* S, const CMat *B, int r_max, double tol);

void incsvd_update_rows(IncSVD* S, const CMat *Brow, int r_max, double tol);

void HROM_ddecm_write_selected_elems_inc_svd_init_with_first_block(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      dof,
    const char*    directory);

void HROM_ddecm_write_selected_elems_inc_svd_update(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      dof,
    const char*    directory);

void HROM_ddecm_write_selected_elems_inc_svd(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      max_iter,   /* NNLS */
    const double   tol,        /* NNLS */
    const int      dof,
    const char*    directory);


#endif /* INC_SVD_H_ */
