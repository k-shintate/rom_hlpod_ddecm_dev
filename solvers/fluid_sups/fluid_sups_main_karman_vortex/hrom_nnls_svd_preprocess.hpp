#ifndef HROM_NNLS_SVD_PREPROCESS_HPP
#define HROM_NNLS_SVD_PREPROCESS_HPP

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

/*
 * LAPACK direct singular-value decomposition.
 * MONOLIS already depends on LAPACK for its NNLS linear solves, so the same
 * link line normally provides dsyev_.
 */
extern "C" {
void dgesvd_(char* jobu, char* jobvt, int* m, int* n,
             double* a, int* lda, double* s,
             double* u, int* ldu, double* vt, int* ldvt,
             double* work, int* lwork, int* info);
}

struct HROM_NNLS_SVD_REDUCTION {
    int original_rows = 0;
    int cols = 0;
    int reduced_rows = 0;

    double** A = NULL;
    double* Ablock = NULL;
    double* b = NULL;

    double requested_rtol = 0.0;
    double effective_cutoff = 0.0;
    double sigma_max = 0.0;
    double sigma_min_retained = 0.0;
    double condition_estimate = 0.0;

    double total_sigma_energy = 0.0;
    double discarded_sigma_energy = 0.0;
    double discarded_sigma_energy_ratio = 0.0;

    double b_norm = 0.0;
    double b_orthogonal_norm = 0.0;
    double b_truncated_norm = 0.0;
};

static inline double** hrom_svd_alloc_matrix(
    int m,
    int n,
    double** block)
{
    if (m <= 0 || n <= 0) return NULL;

    *block = (double*)calloc(
        (size_t)m * (size_t)n,
        sizeof(double));

    double** A = (double**)malloc(
        sizeof(double*) * (size_t)m);

    if (*block == NULL || A == NULL) {
        free(*block);
        free(A);
        *block = NULL;
        return NULL;
    }

    for (int i = 0; i < m; ++i) {
        A[i] = (*block) + (size_t)i * (size_t)n;
    }

    return A;
}

static inline void hrom_svd_free_reduction(
    HROM_NNLS_SVD_REDUCTION* s)
{
    if (s == NULL) return;
    free(s->Ablock);
    free(s->A);
    free(s->b);
    s->Ablock = NULL;
    s->A = NULL;
    s->b = NULL;
    s->reduced_rows = 0;
}

static inline double hrom_nnls_residual_l2(
    double** A,
    const double* b,
    const double* x,
    int m,
    int n)
{
    long double r2 = 0.0L;

    for (int i = 0; i < m; ++i) {
        long double ax = 0.0L;
        for (int j = 0; j < n; ++j) {
            ax += (long double)A[i][j] * (long double)x[j];
        }
        const long double r = ax - (long double)b[i];
        r2 += r * r;
    }

    return std::sqrt((double)r2);
}

/*
 * Row-space truncated SVD for NNLS preprocessing.
 *
 * A = U Sigma V^T,  A in R^{m x n}.
 * LAPACK computes a thin SVD directly; this avoids squaring the condition
 * number as would happen with an A^T A eigendecomposition.  Then build
 *
 *     A_red = Sigma_r V_r^T
 *     b_red = U_r^T b.
 *
 * The unknown vector x and therefore the NNLS non-negativity constraints are
 * unchanged.  Only row directions whose singular values are below the cutoff
 * are removed.
 */
static inline HROM_NNLS_SVD_REDUCTION hrom_nnls_svd_reduce_rows(
    double** A,
    const double* b,
    int m,
    int n,
    double svd_rtol)
{
    HROM_NNLS_SVD_REDUCTION out;
    out.original_rows = m;
    out.cols = n;
    out.requested_rtol = svd_rtol;

    if (A == NULL || b == NULL || m <= 0 || n <= 0) {
        fprintf(stderr, "ERROR: invalid matrix passed to SVD preprocessing\n");
        return out;
    }

    const int k = std::min(m, n);

    /* LAPACK uses column-major storage. */
    std::vector<double> a_col((size_t)m * (size_t)n, 0.0);
    for (int c = 0; c < n; ++c) {
        for (int r = 0; r < m; ++r) {
            a_col[(size_t)r + (size_t)c * (size_t)m] = A[r][c];
        }
    }

    std::vector<double> sigma((size_t)k, 0.0);
    std::vector<double> u((size_t)m * (size_t)k, 0.0);
    std::vector<double> vt((size_t)k * (size_t)n, 0.0);

    char jobu = 'S';
    char jobvt = 'S';
    int mm = m;
    int nn = n;
    int lda = m;
    int ldu = m;
    int ldvt = k;
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;

    dgesvd_(
        &jobu,
        &jobvt,
        &mm,
        &nn,
        a_col.data(),
        &lda,
        sigma.data(),
        u.data(),
        &ldu,
        vt.data(),
        &ldvt,
        &work_query,
        &lwork,
        &info);

    if (info != 0) {
        fprintf(
            stderr,
            "ERROR: LAPACK dgesvd workspace query failed: info=%d\n",
            info);
        return out;
    }

    lwork = std::max(1, (int)std::ceil(work_query));
    std::vector<double> work((size_t)lwork, 0.0);

    dgesvd_(
        &jobu,
        &jobvt,
        &mm,
        &nn,
        a_col.data(),
        &lda,
        sigma.data(),
        u.data(),
        &ldu,
        vt.data(),
        &ldvt,
        work.data(),
        &lwork,
        &info);

    if (info != 0) {
        fprintf(
            stderr,
            "ERROR: LAPACK dgesvd failed to converge: info=%d\n",
            info);
        return out;
    }

    const double sigma_max = sigma.empty() ? 0.0 : sigma[0];
    out.sigma_max = sigma_max;

    long double b2 = 0.0L;
    for (int r = 0; r < m; ++r) {
        b2 += (long double)b[r] * (long double)b[r];
    }
    out.b_norm = std::sqrt((double)b2);

    if (sigma_max <= DBL_MIN) {
        fprintf(stderr, "ERROR: SVD preprocessing found a zero matrix\n");
        return out;
    }

    const double user_rtol = (svd_rtol > 0.0) ? svd_rtol : 0.0;
    const double numerical_rtol =
        100.0 * DBL_EPSILON * (double)std::max(m, n);
    const double effective_rtol = std::max(user_rtol, numerical_rtol);
    const double cutoff = effective_rtol * sigma_max;
    out.effective_cutoff = cutoff;

    std::vector<int> retained;
    retained.reserve((size_t)k);

    long double sigma_energy_total = 0.0L;
    long double sigma_energy_discarded = 0.0L;
    long double bproj_full2 = 0.0L;
    long double bproj_retained2 = 0.0L;

    std::vector<double> ub((size_t)k, 0.0);

    for (int mode = 0; mode < k; ++mode) {
        long double dot = 0.0L;
        for (int r = 0; r < m; ++r) {
            dot +=
                (long double)u[(size_t)r + (size_t)mode * (size_t)m]
                * (long double)b[r];
        }
        ub[(size_t)mode] = (double)dot;

        const double sv = sigma[(size_t)mode];
        sigma_energy_total += (long double)sv * (long double)sv;

        if (sv > numerical_rtol * sigma_max) {
            bproj_full2 += dot * dot;
        }

        if (sv >= cutoff && sv > numerical_rtol * sigma_max) {
            retained.push_back(mode);
        }
        else {
            sigma_energy_discarded += (long double)sv * (long double)sv;
        }
    }

    /* Never return an empty system. */
    if (retained.empty()) {
        retained.push_back(0);
    }

    out.reduced_rows = (int)retained.size();
    out.A = hrom_svd_alloc_matrix(
        out.reduced_rows,
        n,
        &out.Ablock);
    out.b = (double*)calloc(
        (size_t)out.reduced_rows,
        sizeof(double));

    if (out.A == NULL || out.b == NULL) {
        fprintf(stderr, "ERROR: failed to allocate SVD-reduced NNLS system\n");
        hrom_svd_free_reduction(&out);
        return out;
    }

    for (int rr = 0; rr < out.reduced_rows; ++rr) {
        const int mode = retained[(size_t)rr];
        const double sv = sigma[(size_t)mode];

        out.b[rr] = ub[(size_t)mode];
        bproj_retained2 +=
            (long double)out.b[rr] * (long double)out.b[rr];

        if (rr == 0) out.sigma_max = sv;
        out.sigma_min_retained = sv;

        for (int c = 0; c < n; ++c) {
            /* VT has shape k x n in column-major storage. */
            const double vtc =
                vt[(size_t)mode + (size_t)c * (size_t)k];
            out.A[rr][c] = sv * vtc;
        }
    }

    out.condition_estimate =
        (out.sigma_min_retained > 0.0)
        ? out.sigma_max / out.sigma_min_retained
        : HUGE_VAL;

    out.total_sigma_energy = (double)sigma_energy_total;
    out.discarded_sigma_energy = (double)sigma_energy_discarded;
    out.discarded_sigma_energy_ratio =
        (sigma_energy_total > 0.0L)
        ? (double)(sigma_energy_discarded / sigma_energy_total)
        : 0.0;

    long double orth2 = b2 - bproj_full2;
    if (orth2 < 0.0L &&
        std::fabs((double)orth2) <
            1.0e-12 * std::max(1.0, (double)b2)) {
        orth2 = 0.0L;
    }

    long double trunc2 = bproj_full2 - bproj_retained2;
    if (trunc2 < 0.0L &&
        std::fabs((double)trunc2) <
            1.0e-12 * std::max(1.0, (double)bproj_full2)) {
        trunc2 = 0.0L;
    }

    out.b_orthogonal_norm =
        std::sqrt(std::max(0.0, (double)orth2));
    out.b_truncated_norm =
        std::sqrt(std::max(0.0, (double)trunc2));

    return out;
}

static inline void hrom_print_svd_summary(
    const char* label,
    const HROM_NNLS_SVD_REDUCTION& s)
{
    printf("------------------------------------------------------------\n");
    printf("[%s SVD preprocessing]\n", label);
    printf("original NNLS rows              = %d\n", s.original_rows);
    printf("candidate columns               = %d\n", s.cols);
    printf("SVD retained rows               = %d\n", s.reduced_rows);
    printf("SVD removed rows                = %d\n", s.original_rows - s.reduced_rows);
    printf("SVD requested relative cutoff   = %.3e\n", s.requested_rtol);
    printf("sigma max                       = %.15e\n", s.sigma_max);
    printf("sigma min retained              = %.15e\n", s.sigma_min_retained);
    printf("effective sigma cutoff          = %.15e\n", s.effective_cutoff);
    printf("retained condition estimate     = %.15e\n", s.condition_estimate);
    printf("discarded singular energy       = %.15e %%\n",
           100.0 * s.discarded_sigma_energy_ratio);
    printf("||b||                           = %.15e\n", s.b_norm);
    printf("||b_perp(col(A))||              = %.15e\n", s.b_orthogonal_norm);
    printf("||b_truncated||                  = %.15e\n", s.b_truncated_norm);
    printf("------------------------------------------------------------\n");
}

#endif
