#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "DDHR_para_lb.h"
#include "monolis_wrapper_scalapack_c.h"
#include "inc_svd.h"

static const char* OUTPUT_FILENAME_RANK_INC_SVD = "rank_inc_svd.txt";

static inline CMat cm_alloc(int rows, int cols){
    CMat M = {rows, cols, (double*)calloc((size_t)rows*cols, sizeof(double))};
    return M;
}


void cm_free(CMat *M){ free(M->a); M->a=NULL; M->rows=M->cols=0; }
static inline double cm_get(const CMat *M, int i, int j){ return M->a[(size_t)i + (size_t)j*M->rows]; }
void   cm_set(CMat *M, int i, int j, double v){ M->a[(size_t)i + (size_t)j*M->rows] = v; }
void   cm_copy(const CMat *A, CMat *B){
    assert(A->rows==B->rows && A->cols==B->cols);
    memcpy(B->a, A->a, (size_t)A->rows*A->cols*sizeof(double));
}

/* C = A*B (col-major) */
static CMat matmul_nn(const CMat *A, const CMat *B){
    assert(A->cols == B->rows);
    int m=A->rows, k=A->cols, n=B->cols;
    CMat C = cm_alloc(m,n);
    for(int j=0;j<n;++j){
        for(int p=0;p<k;++p){
        double bpj = cm_get(B,p,j);
        const double *ai = &A->a[(size_t)p*A->rows];
        double *ci = &C.a[(size_t)j*C.rows];
        for(int i=0;i<m;++i) ci[i] += ai[i]*bpj;
        }
    }
    return C;
}

/* C = A^T * B */
static CMat matmul_tn(const CMat *A, const CMat *B){
    assert(A->rows == B->rows);
    int m=A->rows, k=A->cols, n=B->cols;
    CMat C = cm_alloc(k,n);
    for(int j=0;j<n;++j){
        for(int p=0;p<k;++p){
            double s=0.0;
            for(int i=0;i<m;++i) s += cm_get(A,i,p)*cm_get(B,i,j);
            cm_set(&C,p,j,s);
        }
    }
    return C;
}


static CMat mat_sub(const CMat *A, const CMat *B){
    assert(A->rows==B->rows && A->cols==B->cols);
    CMat C = cm_alloc(A->rows, A->cols);
    int n=A->rows*A->cols;
    for(int t=0;t<n;++t) C.a[t] = A->a[t] - B->a[t];
    return C;
}


static CMat concat_cols(const CMat *A, const CMat *B){
    assert(A->rows == B->rows);
    CMat C = cm_alloc(A->rows, A->cols + B->cols);
    for(int j=0;j<A->cols;++j)
        memcpy(&C.a[(size_t)j*C.rows], &A->a[(size_t)j*A->rows], (size_t)A->rows*sizeof(double));
    for(int j=0;j<B->cols;++j)
        memcpy(&C.a[(size_t)(j+A->cols)*C.rows], &B->a[(size_t)j*B->rows], (size_t)B->rows*sizeof(double));
    return C;
}


static CMat block_diag(const CMat *A, const CMat *B){
    CMat C = cm_alloc(A->rows + B->rows, A->cols + B->cols);
    for(int j=0;j<A->cols;++j) for(int i=0;i<A->rows;++i) cm_set(&C,i,j, cm_get(A,i,j));
    for(int j=0;j<B->cols;++j) for(int i=0;i<B->rows;++i) cm_set(&C, A->rows+i, A->cols+j, cm_get(B,i,j));
    return C;
}


static CMat eye(int n){
    CMat Imat = cm_alloc(n,n);
    for(int i=0;i<n;++i) cm_set(&Imat, i, i, 1.0);
    return Imat;
}
static CMat zeros(int m, int n){ return cm_alloc(m,n); }


static CMat diag_from_vector(const double *sigma, int r){
    CMat S = cm_alloc(r,r);
    for(int i=0;i<r;++i) cm_set(&S, i, i, sigma[i]);
    return S;
}


static CMat block_2x2(const CMat *A11, const CMat *A12, const CMat *A21, const CMat *A22){
    assert(A11->rows==A12->rows && A21->rows==A22->rows);
    assert(A11->cols==A21->cols && A12->cols==A22->cols);
    int r1=A11->rows, c1=A11->cols, r2=A22->rows, c2=A22->cols;
    CMat K = cm_alloc(r1+r2, c1+c2);
    for(int j=0;j<c1;++j) for(int i=0;i<r1;++i) cm_set(&K,i,j, cm_get(A11,i,j));
    for(int j=0;j<c2;++j) for(int i=0;i<r1;++i) cm_set(&K,i,c1+j, cm_get(A12,i,j));
    for(int j=0;j<c1;++j) for(int i=0;i<r2;++i) cm_set(&K,r1+i,j, cm_get(A21,i,j));
    for(int j=0;j<c2;++j) for(int i=0;i<r2;++i) cm_set(&K,r1+i,c1+j, cm_get(A22,i,j));
    return K;
}


void qr_mgs2(const CMat *E, CMat *Q, CMat *R, double tol){
    int K=E->rows, D=E->cols;
    *Q = cm_alloc(K, D);
    *R = cm_alloc(D, D);
    int q=0;
    for(int j=0;j<D;++j){
        for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(E,i,j));
        for(int k=0;k<q;++k){
        double r = 0.0;
        for(int i=0;i<K;++i) r += cm_get(Q,i,k)*cm_get(Q,i,q);
        cm_set(R,k,j, r);
        for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(Q,i,q) - r*cm_get(Q,i,k));
        }
        for(int k=0;k<q;++k){
        double r = 0.0;
        for(int i=0;i<K;++i) r += cm_get(Q,i,k)*cm_get(Q,i,q);
        cm_set(R,k,j, cm_get(R,k,j)+r);
        for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(Q,i,q) - r*cm_get(Q,i,k));
        }
        double nrm=0.0; for(int i=0;i<K;++i){ double v=cm_get(Q,i,q); nrm += v*v; }
        nrm = sqrt(nrm);
        if(nrm > tol){
        cm_set(R,q,j, nrm);
        for(int i=0;i<K;++i) cm_set(Q,i,q, cm_get(Q,i,q)/nrm);
        ++q;
        }
    }
    Q->cols = q; R->rows = q;
}


void jacobi_eigensymm(const CMat *Gin, double *lam, CMat *Q, int max_sweeps, double tol){
    int n = Gin->rows; 
    CMat G = cm_alloc(n,n); cm_copy(Gin,&G);
    *Q = eye(n);
    for(int sweep=0; sweep<max_sweeps; ++sweep){
        double off=0.0;
        for(int j=0;j<n;++j) for(int i=0;i<j;++i){ double g=cm_get(&G,i,j); off+=g*g; }
        if(sqrt(off) < tol) break;
        for(int p=0;p<n-1;++p){
        for(int q=p+1;q<n;++q){
            double gpq=cm_get(&G,p,q); if(fabs(gpq)<=0.0) continue;
            double gpp=cm_get(&G,p,p), gqq=cm_get(&G,q,q);
            double tau=(gqq-gpp)/(2.0*gpq);
            double t = (tau>=0)? 1.0/(tau+sqrt(1.0+tau*tau)) : -1.0/(-tau+sqrt(1.0+tau*tau));
            double c=1.0/sqrt(1.0+t*t), s=t*c;
            for(int k=0;k<n;++k){
            if(k==p||k==q) continue;
            double gpk=cm_get(&G,p,k), gqk=cm_get(&G,q,k);
            double npk=c*gpk - s*gqk, nqk=s*gpk + c*gqk;
            cm_set(&G,p,k,npk); cm_set(&G,k,p,npk);
            cm_set(&G,q,k,nqk); cm_set(&G,k,q,nqk);
            }
            double gpp_new=c*c*gpp - 2*c*s*gpq + s*s*gqq;
            double gqq_new=s*s*gpp + 2*c*s*gpq + c*c*gqq;
            cm_set(&G,p,p,gpp_new); cm_set(&G,q,q,gqq_new);
            cm_set(&G,p,q,0.0); cm_set(&G,q,p,0.0);
            for(int i=0;i<n;++i){
            double vip=cm_get(Q,i,p), viq=cm_get(Q,i,q);
            cm_set(Q,i,p, c*vip - s*viq);
            cm_set(Q,i,q, s*vip + c*viq);
            }
        }
        }
    }
    for(int i=0;i<n;++i) lam[i]=cm_get(&G,i,i);
    cm_free(&G);
}


void sort_eigs_desc(double *lam, CMat *Q){
    int n=Q->cols;
    for(int i=0;i<n-1;++i){
        int p=i; for(int j=i+1;j<n;++j) if(lam[j]>lam[p]) p=j;
        if(p!=i){
        double t=lam[i]; lam[i]=lam[p]; lam[p]=t;
        for(int r=0;r<Q->rows;++r){
            double tmp=cm_get(Q,r,i);
            cm_set(Q,r,i, cm_get(Q,r,p));
            cm_set(Q,r,p, tmp);
        }
        }
    }
}


void svd_dense(const CMat *Ksmall, CMat *Uhat, double **sigma_out, CMat *Vhat){
    int m=Ksmall->rows, n=Ksmall->cols, r=(m<n?m:n);

    CMat G = matmul_tn(Ksmall, Ksmall);      /* n×n */
    CMat Q; double *lam = (double*)calloc((size_t)n,sizeof(double));
    jacobi_eigensymm(&G, lam, &Q, 50, 1e-12); cm_free(&G);
    sort_eigs_desc(lam, &Q);
    double *S = (double*)calloc((size_t)r,sizeof(double));
    for(int i=0;i<r;++i) S[i] = lam[i]>0.0? sqrt(lam[i]) : 0.0;
    *Vhat = Q;
    CMat KV  = matmul_nn(Ksmall, Vhat);      /* m×n */
    *Uhat = cm_alloc(m,m);
    for(int j=0;j<m;++j){
        double inv = (j<r && S[j]>1e-15)? 1.0/S[j] : 0.0;
        for(int i=0;i<m;++i) cm_set(Uhat,i,j, cm_get(&KV,i,j)*inv);
    }
    cm_free(&KV); free(lam);
    *sigma_out = S;
  
}


static int select_rank(const double *sigma, int len, int r_max, double tol){
    int r = 0; double s0 = (len>0? sigma[0]:0.0);
    for(int i=0;i<len && i<r_max; ++i){ if(sigma[i] >= tol*s0) ++r; else break; }
    if(r==0 && len>0) r=1;
    return r;

    printf("select_k_by_cum_threshold: N0=%d, r_max=%d, tol=%.2e, r=%d\n", len, r_max, tol, r);

}


void truncate_uv(CMat *U_new, double **sigma_new, CMat *V_new, int r_keep){
    int m=U_new->rows, nV=V_new->rows;
    CMat U2 = cm_alloc(m, r_keep);
    for(int j=0;j<r_keep;++j)
        memcpy(&U2.a[(size_t)j*U2.rows], &U_new->a[(size_t)j*U_new->rows], (size_t)m*sizeof(double));
    cm_free(U_new); *U_new = U2;
    double *s2 = (double*)calloc((size_t)r_keep,sizeof(double));
    memcpy(s2, *sigma_new, (size_t)r_keep*sizeof(double));
    free(*sigma_new); *sigma_new = s2;
    CMat V2 = cm_alloc(nV, r_keep);
    for(int j=0;j<r_keep;++j)
        memcpy(&V2.a[(size_t)j*V2.rows], &V_new->a[(size_t)j*V_new->rows], (size_t)nV*sizeof(double));
    cm_free(V_new); *V_new = V2;
}


int select_k_by_cum_threshold(const double* Svals, int N0, int r_max, double tol) {
    const int upto   = N0;
    double total = 0.0;
    for (int i = 0; i < upto; ++i) total += Svals[i];

    printf("select_k_by_cum_threshold: total = %lf\n", total);

    int r = 0;
    if (total > 0.0) {
        double cum = 0.0;
        for (int i = 0; i < upto; ++i) {
            cum += Svals[i];
            if (cum / total >= 1 - tol) {  // しきい値到達
                r = i + 1;             // 要素数としての k
                break;
            }
        }
        if (r == 0) r = (tol <= 0.0 ? 1 : upto);
    } else {
        r = 0;
    }

    if (r == 0 && upto > 0) r = 1;

    printf("select_k_by_cum_threshold: N0=%d, r_max=%d, tol=%.2e, r=%d\n", N0, r_max, tol, r);

    return r;
}


void incsvd_update(IncSVD* S, const CMat *B, int r_max, double tol){
    const int K = S->K;
    const int Delta = B->cols;
    const int r = S->r;
    if (Delta <= 0) return;

    CMat U = (CMat){K, r, S->U};
    CMat C = (r>0) ? matmul_tn(&U, B) : cm_alloc(0, Delta); /* r×Δ */

    CMat UC = (r>0) ? matmul_nn(&U, &C) : cm_alloc(K, Delta);
    CMat E  = mat_sub(B, &UC); cm_free(&UC);
    CMat Q,R; qr_mgs2(&E, &Q, &R, 1e-14); cm_free(&E);
    int q = Q.cols;

    CMat Sdiag = diag_from_vector(S->sigma, r);
    CMat K11 = Sdiag, K12 = C, K21 = zeros(q,r), K22 = R;
    CMat Ksmall = block_2x2(&K11,&K12,&K21,&K22);

    CMat Uhat, Vhat; double *sigma_new=NULL;
    svd_dense(&Ksmall, &Uhat, &sigma_new, &Vhat);
    cm_free(&K11); cm_free(&K21); cm_free(&K22); cm_free(&Ksmall);

    CMat Ucm  = (CMat){K, r, S->U};
    CMat U_aug = (q>0) ? concat_cols(&Ucm, &Q) : Ucm;
    CMat U_new = matmul_nn(&U_aug, &Uhat);

    CMat Vcm  = (CMat){S->N, r, S->V};
    CMat I_delta = eye(Delta);
    CMat V_aug = block_diag(&Vcm, &I_delta);
    CMat V_new = matmul_nn(&V_aug, &Vhat);

    int len_sigma = (Uhat.rows < Vhat.cols) ? Uhat.rows : Vhat.cols;
    int r_keep = select_rank(sigma_new, len_sigma, r_max, tol);
    truncate_uv(&U_new, &sigma_new, &V_new, r_keep);

    free(S->U); free(S->V); free(S->sigma);
    S->U = U_new.a; S->V = V_new.a; S->sigma = sigma_new;
    S->r = r_keep; S->N += Delta;

    if(q>0) cm_free(&Q);
    cm_free(&C); cm_free(&Uhat); cm_free(&Vhat); cm_free(&I_delta); cm_free(&V_aug);
}


void init_with_first_block_by_monolis(
    IncSVD *S, double **matrixRM, int K, int N0,
    int r_max, double tol, int comm, int scalapack_comm)
{
    if (K<=0 || N0<=0){ S->K=S->N=S->r=0; S->U=S->V=S->sigma=NULL; return; }

    double **A0 = BB_std_calloc_2d_double(A0, K, N0);
    for (int j=0;j<K;++j){
        for (int i=0;i<N0;++i){
            A0[j][i] = matrixRM[j][i];
        }
    }

    double **U = BB_std_calloc_2d_double(U, K, N0);
    double **V = BB_std_calloc_2d_double(V, N0, N0);
    double *Svals = BB_std_calloc_1d_double(Svals, N0);

    monolis_scalapack_gesvd_R(K, N0, A0, U, Svals, V, comm, scalapack_comm);
    int r = select_k_by_cum_threshold(Svals, N0, r_max, tol);

    S->K = K; S->N = N0; S->r = r;
    S->U = (double*)malloc((size_t)K*r*sizeof(double));
    S->V = (double*)malloc((size_t)N0*r*sizeof(double));
    S->sigma = (double*)malloc((size_t)r*sizeof(double));

    for (int j=0;j<r;++j){
        S->sigma[j] = Svals[j];
        for (int i=0;i<K;++i){
            S->U[i + (size_t)j*K] = U[i][j];
        }
        for (int i=0;i<N0;++i){ 
            S->V[i + (size_t)j*N0] = V[j][i]; /* V(i,j) */
        }
    }

    BB_std_free_2d_double(A0, K, N0);
    BB_std_free_2d_double(U,  K, N0);
    BB_std_free_2d_double(V,  N0, N0);
    BB_std_free_1d_double(Svals, N0);
}


static CMat copy_block_colmajor_from_rowmajor(double **matrixRM, int K, int N, int j0, int Delta){
    CMat B = cm_alloc(K, Delta);
    for(int j=0;j<Delta;++j){
        int jj = j0 + j;
        for(int i=0;i<K;++i) B.a[(size_t)i + (size_t)j*K] = matrixRM[i][jj];
    }
    return B;
}


static double** alloc2d_rowmajor(int rows, int cols){
    double **M = (double**)malloc((size_t)rows*sizeof(double*));
    for(int i=0;i<rows;++i) M[i] = (double*)calloc((size_t)cols,sizeof(double));
    return M;
}


void free2d_rowmajor(double **M, int rows){
    for(int i=0;i<rows;++i) free(M[i]);
    free(M);
}


void build_compressed_system_from_incsvd(
    const IncSVD* S,
    double*** G_k_out, double** b_k_out)
{
    int r=S->r, K=S->K, N=S->N;
    double **G_k = BB_std_calloc_2d_double(G_k, r, N);
    double *b_k = (double*)calloc((size_t)r,sizeof(double));
    for(int i=0;i<r;++i){
        double sii = S->sigma[i];
        for(int j=0;j<N;++j){
        double Vij = S->V[(size_t)j + (size_t)i*N];  /* V: N×r */
        G_k[i][j] = sii * Vij;
        b_k[i] += G_k[i][j];
        }
    }

    *G_k_out = G_k;
    *b_k_out = b_k;
}


IncSVD* g_svd;
int     g_svd_count = 0;

void ensure_svd_states(int num_subdomains){
    if (g_svd && g_svd_count == num_subdomains) return;
    if (g_svd){
        for (int i=0;i<g_svd_count;++i){
            free(g_svd[i].U); free(g_svd[i].V); free(g_svd[i].sigma);
        }
        free(g_svd);
    }
    g_svd = (IncSVD*)calloc((size_t)num_subdomains, sizeof(IncSVD));
    g_svd_count = num_subdomains;
}


static CMat cm_transpose_new(const CMat *A){
    CMat T = cm_alloc(A->cols, A->rows);
    for(int j=0;j<A->cols;++j)
        for(int i=0;i<A->rows;++i)
        cm_set(&T, j, i, cm_get(A, i, j));
    return T;
}


void incsvd_update_rows(IncSVD* S, const CMat *Brow, int r_max, double tol){
    /* 事前条件チェック */
    if (!S || !S->U || !S->V || !S->sigma) return;
    if (Brow->cols != S->N || Brow->rows <= 0) return;

    const int DeltaK = Brow->rows;
    const int N      = S->N;
    const int K      = S->K;
    const int r      = S->r;

    CMat B_T = cm_transpose_new(Brow);


    IncSVD St = {0};
    St.K = N;          /* 行数 ← 元の列数 */
    St.N = K;          /* 列数 ← 元の行数 */
    St.r = r;

    St.U = (double*)malloc((size_t)St.K * r * sizeof(double));   /* ← V(N×r) をコピー */
    St.V = (double*)malloc((size_t)St.N * r * sizeof(double));   /* ← U(K×r) をコピー */
    St.sigma = (double*)malloc((size_t)r * sizeof(double));

    /* V → St.U */
    for(int j=0;j<r;++j)
        memcpy(&St.U[(size_t)j*St.K], &S->V[(size_t)j*N], (size_t)N*sizeof(double));
    /* U → St.V */
    for(int j=0;j<r;++j)
        memcpy(&St.V[(size_t)j*St.N], &S->U[(size_t)j*K], (size_t)K*sizeof(double));
    /* sigma */
    memcpy(St.sigma, S->sigma, (size_t)r*sizeof(double));

    /* 3) 既存の列追加を利用（A^T に ΔK 列を追加） */
    incsvd_update(&St, &B_T, r_max, tol);

    /* 4) St を元に S を更新（再転置：U↔V） */
    double *oldU = S->U, *oldV = S->V, *oldS = S->sigma;

    S->K = K + DeltaK;  /* 行数が増える */
    S->N = N;           /* 列数は不変 */
    S->r = St.r;

    S->U = (double*)malloc((size_t)S->K * S->r * sizeof(double));
    S->V = (double*)malloc((size_t)S->N * S->r * sizeof(double));
    S->sigma = (double*)malloc((size_t)S->r * sizeof(double));

    /* St.V（= 新しい U^T の“右特異ベクトル”） → S.U（(K+ΔK)×r） */
    for(int j=0;j<S->r;++j)
        memcpy(&S->U[(size_t)j*S->K], &St.V[(size_t)j*(K+DeltaK)], (size_t)(K+DeltaK)*sizeof(double));
    /* St.U（= 新しい V^T の“左特異ベクトル”） → S.V（N×r） */
    for(int j=0;j<S->r;++j)
        memcpy(&S->V[(size_t)j*S->N], &St.U[(size_t)j*N], (size_t)N*sizeof(double));
    memcpy(S->sigma, St.sigma, (size_t)S->r*sizeof(double));

    free(oldU); free(oldV); free(oldS);
    free(St.U); free(St.V); free(St.sigma);
    cm_free(&B_T);
}


static CMat copy_rowblock_colmajor_from_rowmajor(double **matrixRM, int K, int N, int i0, int DeltaK){
    if (DeltaK <= 0) return cm_alloc(0,0);
    CMat B = cm_alloc(DeltaK, N);
    for(int j=0;j<N;++j){
        for(int i=0;i<DeltaK;++i){
        int src_i = i0 + i;
        B.a[(size_t)i + (size_t)j*DeltaK] = matrixRM[src_i][j];
        }
    }
    return B;
}


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
    const char*    directory)
{
    ensure_svd_states(num_subdomains);

    const int r_max = 1000;
    const double svd_tol = 1.0e-12;

    const int comm = monolis_mpi_get_self_comm();
    int scalapack_comm; monolis_scalapack_comm_initialize(comm, &scalapack_comm);

        const int K = hlpod_ddhr->num_modes_1stdd[0] * total_num_snapshot + 1;
        const int N = hlpod_ddhr->num_elems[0];

        double **matrix = BB_std_calloc_2d_double(matrix, K, N);
        for (int j=0; j<K; ++j){
            for (int i=0; i<N; ++i){
                matrix[j][i] = hlpod_ddhr->matrix[j][i][0];
            }
        }

        init_with_first_block_by_monolis(&g_svd[0], matrix, K, N, r_max, svd_tol, comm, scalapack_comm);

        BB_std_free_2d_double(matrix, K, N);

    int K0 = g_svd[0].K;

}


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
    const char*    directory)
{
    const int r_max = 300;
    const double svd_tol = 1.0e-6;

    const int K_total = hlpod_ddhr->num_modes_1stdd[0] * total_num_snapshot  + 1;
    const int N       = hlpod_ddhr->num_elems[0];

    int K0 = g_svd[0].K;

    printf("Incremental SVD: subdomain %d, K_total=%d, N=%d, K0=%d, r=%d\n",
            0, K_total, N, K0, g_svd[0].r);

    double **matrix = BB_std_calloc_2d_double(matrix, K_total, N);
    for (int j=0; j<K_total; ++j){
        for (int i=0; i<N; ++i){
            matrix[j][i] = hlpod_ddhr->matrix[j][i][0];
        }
    }
    
    int g_block_cols =  100;

    /* 列ブロックを一定幅ずつ追加（例：g_block_cols を行方向にも使う）*/
    for (int i0 = 0; i0 < K_total; i0 += g_block_cols){
        int DeltaK = (i0 + g_block_cols <= K_total ? g_block_cols : (K_total - i0));
        CMat Brow = copy_rowblock_colmajor_from_rowmajor(matrix, K_total, N, i0, DeltaK);
        incsvd_update_rows(&g_svd[0], &Brow, r_max, svd_tol);
        cm_free(&Brow);
    }
    g_svd[0].K = K0;

    BB_std_free_2d_double(matrix, K_total, N);

}


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
    const int      max_iter,
    const double   tol,
    const int      dof,
    const char*    directory)
{

    int BUFFER_SIZE = 1024;
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;
   	const int myrank = monolis_mpi_get_global_my_rank();

	int index_1 = 0;
	int index_2 = 0;

	int* subdomain_id;
	subdomain_id = BB_std_calloc_1d_int(subdomain_id, num_subdomains);

	snprintf(fname, BUFFER_SIZE, "metagraph_parted.0/metagraph.dat.id.%d", myrank);

	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s", id);
	fscanf(fp, "%d %d", &(ndof), &(ndof));
	for (int i = 0; i < num_subdomains; i++) {
		fscanf(fp, "%d", &(subdomain_id[i]));
	}
	fclose(fp);

	int nl = fe->local_num_nodes;
    double t1 = monolis_get_time_global_sync();

	double residual;

	bool** bool_elem;
	int* total_id_selected_elems;
	double* total_elem_weight;
	int* total_num_selected_elems;

    const int max_ITER = 1000;
    const double TOL   = 1.0e-10;

	/**/
	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);

	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
	hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);

	hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
	hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);

	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
	/**/

	//lbから追加
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);

	int Index1 = 0;
	int Index2 = 0;

	int index_NNLS1 = 0;
	int index_NNLS2 = 0;


    const int K = g_svd[0].K;
    const int N = g_svd[0].N;
    printf("Finalizing subdomain %d: K=%d, N=%d, r=%d\n", 0, K, N, g_svd[0].r);

    
	const char* filename;
	snprintf(fname, BUFFER_SIZE,"hr_prm/%s", OUTPUT_FILENAME_RANK_INC_SVD);
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname);
	fp = ROM_BB_write_fopen(fp, filename, directory);
    fprintf(fp, "RANK_INC_SVD\n", fe->total_num_nodes);
	fprintf(fp, "%d\n", g_svd[0].r);
	fclose(fp);


    double **G_k = NULL; double *b_k = NULL;
    build_compressed_system_from_incsvd(&g_svd[0], &G_k, &b_k); /* r×N, r */

    double *ans_vec = BB_std_calloc_1d_double(ans_vec, N);
    monolis_optimize_nnls_R_with_sparse_solution(
        G_k, b_k, ans_vec, g_svd[0].r, N, max_ITER, TOL, &residual);

    int index = 0;
    for (int i = 0; i < hlpod_ddhr->num_elems[0]; i++) {
        if (ans_vec[i] != 0.0) {
            index++;
        }
    }

    total_num_selected_elems[0] = index;

    hr_write_NNLS_residual(residual, myrank, 0, directory);
    hr_write_NNLS_num_elems(total_num_selected_elems[0], myrank, 0, directory);

    total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[0]);
    total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[0]);

    index = 0;
    for (int i = 0; i < hlpod_ddhr->num_elems[0]; i++) {
        if (ans_vec[i] != 0.0) {
            total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][0];
            total_elem_weight[index] = ans_vec[i];
            index++;
        }
    }

    for (int h = 0; h < (total_num_selected_elems[0]); h++) {
        int e = total_id_selected_elems[h];

        for (int i = 0; i < nl; i++) {       //六面体一次要素は8
            for (int j = 0; j < nl; j++) {
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                for(int k = 0; k < dof; k++) {
                    if (bc->D_bc_exists[index_j*dof + k]) {
                        bool_elem[h][0] = true;
                        hlpod_ddhr->D_bc_exists[index_j][0] = true;
                    }
                }
            }
        }
    }

    index = 0;
    for (int h = 0; h < (total_num_selected_elems[0]); h++) {
        if (bool_elem[h][0]) {
            index++;
        }
    }

    //index = D_bcが付与された要素数
    hlpod_ddhr->num_selected_elems[0] = total_num_selected_elems[0] - index;
    hlpod_ddhr->num_selected_elems_D_bc[0] = index;

    int index1 = 0;
    int index2 = 0;
    for (int h = 0; h < (total_num_selected_elems[0]); h++) {
        int e = total_id_selected_elems[h];

        if (bool_elem[h][0]) {
            hlpod_ddhr->id_selected_elems_D_bc[index1][0] = total_id_selected_elems[h];
            hlpod_ddhr->elem_weight_D_bc[index1][0] = total_elem_weight[h];

            index1++;
        }
        else {
            hlpod_ddhr->id_selected_elems[index2][0] = total_id_selected_elems[h];
            hlpod_ddhr->elem_weight[index2][0] = total_elem_weight[h];

            index2++;
        }
    }

    Index1 += index1;
    Index2 += index2;

    BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[0]);
    BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[0]);

    double t = monolis_get_time_global_sync();

    FILE* fp1;
    FILE* fp2;
    char fname1[BUFFER_SIZE];
    char fname2[BUFFER_SIZE];

    snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[0]);
    snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[0]);

    fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
    fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

    fprintf(fp1, "%d\n", index1);
    fprintf(fp2, "%d\n", index2);

    index_1 = 0;
    index_2 = 0;

    for (int h = 0; h < (total_num_selected_elems[0]); h++) {
        if (bool_elem[h][0]) {
            fprintf(fp1, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems_D_bc[index_1][0]], hlpod_ddhr->elem_weight_D_bc[index_1][0]);
            index_1++;
        }
        else {
            fprintf(fp2, "%d %.30e\n", hlpod_ddhr->parallel_elems_id[hlpod_ddhr->id_selected_elems[index_2][0]], hlpod_ddhr->elem_weight[index_2][0]);

            index_2++;
        }
    }

    fclose(fp1);
    fclose(fp2);

	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);

	int max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);

    BB_std_free_2d_double(G_k, g_svd[0].r, N);
    BB_std_free_1d_double(b_k, g_svd[0].r);
    BB_std_free_1d_double(ans_vec, N);

    /* SVD 状態も解放（以降使わない前提） */
    free(g_svd[0].U); free(g_svd[0].V); free(g_svd[0].sigma);
    g_svd[0].U=g_svd[0].V=g_svd[0].sigma=NULL;
    g_svd[0].K = g_svd[0].N = g_svd[0].r = 0;

    free(g_svd); g_svd=NULL; g_svd_count=0;
}
