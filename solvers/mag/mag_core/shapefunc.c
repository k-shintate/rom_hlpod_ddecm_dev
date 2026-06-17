
#include "shapefunc.h"


#define TET2_NEDGE        6
#define TET2_NFACE        4
#define TET2_EDGE_MODE    2
#define TET2_FACE_MODE    2
#define TET2_NDOF         20
#define TET2_FACE_OFFSET  12

static const int tet_face_conn[4][3] = {
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
};

static inline void swap_int_2nd(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

static void cross3_2nd(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void BBFE_std_shapefunc_hex1st_get_derivative_ned(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{	
	double coef = 1.0/4.0;
	dN_dxi[3] =  coef *               (1.0+xi[1]) * (1.0+xi[2]); //4
	dN_dxi[1] =  coef *               (1.0+xi[1]) * (1.0-xi[2]); //2
	dN_dxi[2] =  coef *               (1.0-xi[1]) * (1.0+xi[2]); //3
	dN_dxi[0] =  coef *               (1.0-xi[1]) * (1.0-xi[2]); //1

	dN_det[7] =  coef * (1.0+xi[0])               * (1.0+xi[2]); //8
	dN_det[6] =  coef * (1.0-xi[0])               * (1.0+xi[2]); //7
	dN_det[5] =  coef * (1.0+xi[0])               * (1.0-xi[2]); //6
	dN_det[4] =  coef * (1.0-xi[0])               * (1.0-xi[2]); //5

	dN_dze[11] =  coef * (1.0+xi[0]) * (1.0+xi[1])              ; //12
	dN_dze[10] =  coef * (1.0-xi[0]) * (1.0+xi[1])              ; //11
	dN_dze[9] =  coef * (1.0+xi[0]) * (1.0-xi[1])              ; //10
	dN_dze[8] =  coef * (1.0-xi[0]) * (1.0-xi[1])              ; //9

}

void BBFE_std_shapefunc_hex1st_get_derivative_ned_curl(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{	
	double coef = 1.0/4.0;
	dN_dxi[3] =  coef *               (1.0+xi[1]);
	dN_dxi[1] = -coef *               (1.0+xi[1]);
	dN_dxi[2] =  coef *               (1.0-xi[1]);
	dN_dxi[0] = -coef *               (1.0-xi[1]);

	dN_det[7] =  coef                         * (1.0+xi[2]);
	dN_det[6] = -coef                         * (1.0+xi[2]);
	dN_det[5] =  coef                         * (1.0-xi[2]);
	dN_det[4] = -coef                         * (1.0-xi[2]);

	dN_dze[11] =  coef * (1.0+xi[0]);
	dN_dze[10] =  coef * (1.0-xi[0]);
	dN_dze[9] = -coef * (1.0+xi[0]);
	dN_dze[8] = -coef * (1.0-xi[0]);
}


void BBFE_std_shapefunc_hex1st_get_derivative_ned_curl2(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{	
	double coef = 1.0/4.0;

	dN_dxi[3] =  coef *                           (1.0+xi[2]);
	dN_dxi[1] =  coef *                           (1.0-xi[2]);
	dN_dxi[2] = -coef *                           (1.0+xi[2]);
	dN_dxi[0] = -coef *                           (1.0-xi[2]);

	dN_det[7] =  coef * (1.0+xi[0]);
	dN_det[6] =  coef * (1.0-xi[0]);
	dN_det[5] = -coef * (1.0+xi[0]);
	dN_det[4] = -coef * (1.0-xi[0]);

	dN_dze[11] = coef               * (1.0+xi[1]);
	dN_dze[10] = -coef              * (1.0+xi[1]);
	dN_dze[9] = coef               * (1.0-xi[1]);
	dN_dze[8] = -coef              * (1.0-xi[1]);

}

void BBFE_std_shapefunc_hex1st_nedelec_get_val(
    const double xi[3],      // 参照座標 (例えば積分点) (xi, eta, zeta)
    double** N_edge,         // 出力：12エッジ分の基底関数値、各エッジは 3 成分のベクトル
    const double J_inv[3][3])
{
    double dphi_dxi[12], dphi_det[12], dphi_dze[12];

    BBFE_std_shapefunc_hex1st_get_derivative_ned(xi, dphi_dxi, dphi_det, dphi_dze);

    for (int e = 0; e < 12; e++) {
        if (e == 0 || e == 1 || e == 2 || e == 3) {
            for (int d = 0; d < 3; ++d) {
                N_edge[e][d] = J_inv[0][d] * dphi_dxi[e];  // ξ方向 → ref成分 0
            }
        }
        else if (e == 4 || e == 5 || e == 6 || e == 7) {
            for (int d = 0; d < 3; ++d) {
                N_edge[e][d] = J_inv[1][d] * dphi_det[e];  // η方向 → ref成分 1
            }
        }
        else { // 8..11
            for (int d = 0; d < 3; ++d) {
                N_edge[e][d] = J_inv[2][d] * dphi_dze[e];  // ζ方向 → ref成分 2
            }
        }

    }
}


void BBFE_std_shapefunc_hex1st_nedelec_get_curl(
    const double xi[3],
    double **curl_N_edge,          /* [12][3] */
    const double J[3][3],          /* fe->geo[e][p].J : dx/dxi */
    const double detJ)             /* fe->geo[e][p].Jacobian */
{
    const double ksi  = xi[0];
    const double eta  = xi[1];
    const double zeta = xi[2];

    const double c = 0.25;

    double curl_hat[12][3];
    for (int e = 0; e < 12; ++e) {
        curl_hat[e][0] = 0.0;
        curl_hat[e][1] = 0.0;
        curl_hat[e][2] = 0.0;
    }

    /* === 参照要素上の curl_{xi} \hat N_e を計算 === */

    /* --- エッジ 0..3 : ξ 方向 (x方向) --- */
    /* φ_e(η,ζ):
       e=0: φ0 = c*(1-η)*(1-ζ)
       e=1: φ1 = c*(1+η)*(1-ζ)
       e=2: φ2 = c*(1-η)*(1+ζ)
       e=3: φ3 = c*(1+η)*(1+ζ)
       curl \hat N = (0, ∂φ/∂ζ, -∂φ/∂η)
    */
    double dphi_deta[4], dphi_dzeta[4];

    dphi_deta[0]  = -c*(1.0 - zeta);
    dphi_deta[1]  =  c*(1.0 - zeta);
    dphi_deta[2]  = -c*(1.0 + zeta);
    dphi_deta[3]  =  c*(1.0 + zeta);

    dphi_dzeta[0] = -c*(1.0 - eta);
    dphi_dzeta[1] = -c*(1.0 + eta);
    dphi_dzeta[2] =  c*(1.0 - eta);
    dphi_dzeta[3] =  c*(1.0 + eta);

    for (int e = 0; e < 4; ++e) {
        curl_hat[e][0] = 0.0;
        curl_hat[e][1] = dphi_dzeta[e];      /* ∂φ/∂ζ */
        curl_hat[e][2] = - dphi_deta[e];      /* -∂φ/∂η */
    }

    /* --- エッジ 4..7 : η 方向 (y方向) --- */
    /* ψ_e(ξ,ζ):
       e=4: ψ4 = c*(1-ksi)*(1-zeta)
       e=5: ψ5 = c*(1+ksi)*(1-zeta)
       e=6: ψ6 = c*(1-ksi)*(1+zeta)
       e=7: ψ7 = c*(1+ksi)*(1+zeta)
       curl \hat N = (-∂ψ/∂ζ, 0, ∂ψ/∂ξ)
    */
    double dpsi_dksi[4], dpsi_dzeta[4];

    dpsi_dksi[0]  = -c*(1.0 - zeta);
    dpsi_dksi[1]  =  c*(1.0 - zeta);
    dpsi_dksi[2]  = -c*(1.0 + zeta);
    dpsi_dksi[3]  =  c*(1.0 + zeta);

    dpsi_dzeta[0] = -c*(1.0 - ksi);
    dpsi_dzeta[1] = -c*(1.0 + ksi);
    dpsi_dzeta[2] =  c*(1.0 - ksi);
    dpsi_dzeta[3] =  c*(1.0 + ksi);

    for (int k = 0; k < 4; ++k) {
        int e = 4 + k;
        curl_hat[e][0] = - dpsi_dzeta[k];     /* -∂ψ/∂ζ */
        curl_hat[e][1] =  0.0;
        curl_hat[e][2] = dpsi_dksi[k];      /*  ∂ψ/∂ξ */
    }

    /* --- エッジ 8..11 : ζ 方向 (z方向) --- */
    /* χ_e(ξ,η):
       e=8:  χ8  = c*(1-ksi)*(1-eta)
       e=9:  χ9  = c*(1+ksi)*(1-eta)
       e=10: χ10 = c*(1-ksi)*(1+eta)
       e=11: χ11 = c*(1+ksi)*(1+eta)
       curl \hat N = (-∂χ/∂η, -∂χ/∂ξ, 0)
    */
    double dchi_dksi[4], dchi_deta_[4];

    dchi_dksi[0]   = -c*(1.0 - eta);
    dchi_dksi[1]   =  c*(1.0 - eta);
    dchi_dksi[2]   = -c*(1.0 + eta);
    dchi_dksi[3]   =  c*(1.0 + eta);

    dchi_deta_[0]  = -c*(1.0 - ksi);
    dchi_deta_[1]  = -c*(1.0 + ksi);
    dchi_deta_[2]  =  c*(1.0 - ksi);
    dchi_deta_[3]  =  c*(1.0 + ksi);


    for (int k = 0; k < 4; ++k) {
        int e = 8 + k;
        curl_hat[e][0] = -dchi_deta_[k];   /*  ∂χ/∂η */
        curl_hat[e][1] = dchi_dksi[k];    /* -∂χ/∂ξ */
        curl_hat[e][2] =  0.0;
    }

    /* === 物理座標への変換： curl N = (1/detJ) * J * curl_hat === */

    const double inv_detJ = 1.0 / detJ;

    for (int e = 0; e < 12; ++e) {
        for (int i = 0; i < 3; ++i) {
            curl_N_edge[e][i] =
                ( J[0][i] * curl_hat[e][0]
                + J[1][i] * curl_hat[e][1]
                + J[2][i] * curl_hat[e][2] ) * inv_detJ;
        }
    }
}

static inline void cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/* --- Nedelec(1st kind, lowest order) の値 --- */
void BBFE_std_shapefunc_tet1st_nedelec_get_val(
    const double xi[3],
    double **N_edge,              /* [6][3] */
    const double J_inv[3][3])     /* [a(ref)][i(phys)] = dxi_a/dx_i */
{
    double lam[4];
    double dlam_dxi[4], dlam_det[4], dlam_dze[4];

    BBFE_std_shapefunc_tet1st_get_val(xi, lam);
    BBFE_std_shapefunc_tet1st_get_derivative(xi, dlam_dxi, dlam_det, dlam_dze);

    /* grad_lam[k][a] = ∂λ_k / ∂ξ_a  (a=0:xi, 1:eta, 2:zeta) */
    double grad_lam[4][3];
    for (int k = 0; k < 4; ++k) {
        grad_lam[k][0] = dlam_dxi[k];
        grad_lam[k][1] = dlam_det[k];
        grad_lam[k][2] = dlam_dze[k];
    }

    for (int e = 0; e < 6; ++e) {
        const int i = tet_edge_conn[e][0];
        const int j = tet_edge_conn[e][1];

        /* 参照要素上: hatN[a] */
        double hatN[3];
        for (int a = 0; a < 3; ++a) {
            hatN[a] = lam[i]*grad_lam[j][a] - lam[j]*grad_lam[i][a];
        }

        /* 物理要素へ: N = J^{-T} hatN
           N_i = sum_a (∂ξ_a/∂x_i) * hatN_a = sum_a J_inv[a][i] * hatN[a]
        */
	    for(int phys = 0; phys < 3; ++phys){
            N_edge[e][phys] =
                J_inv[phys][0]*hatN[0] +
                J_inv[phys][1]*hatN[1] +
                J_inv[phys][2]*hatN[2];
        }
    }
}


/* --- Nedelec(1st kind, lowest order) の curl --- */
void BBFE_std_shapefunc_tet1st_nedelec_get_curl(
    const double xi[3],           /* 形だけ合わせる（実際は一定） */
    double **curl_N_edge,         /* [6][3] */
    const double J[3][3],         /* [a(ref)][i(phys)] = ∂x_i/∂ξ_a */
    const double detJ)
{
    (void)xi; 

    double lam_dummy[4];
    double dlam_dxi[4], dlam_det[4], dlam_dze[4];
    const double xi_dummy[3] = {0.25, 0.25, 0.25};
    BBFE_std_shapefunc_tet1st_get_val(xi_dummy, lam_dummy);
    BBFE_std_shapefunc_tet1st_get_derivative(xi_dummy, dlam_dxi, dlam_det, dlam_dze);

    double grad_lam[4][3];
    for (int k = 0; k < 4; ++k) {
        grad_lam[k][0] = dlam_dxi[k];
        grad_lam[k][1] = dlam_det[k];
        grad_lam[k][2] = dlam_dze[k];
    }

    const double inv_detJ = 1.0 / detJ;

    for (int e = 0; e < 6; ++e) {
        const int i = tet_edge_conn[e][0];
        const int j = tet_edge_conn[e][1];

        /* curl_hat = 2 * (gradλ_i × gradλ_j) */
        double tmp[3];
        cross3(grad_lam[i], grad_lam[j], tmp);

        double curl_hat[3] = { 2.0*tmp[0], 2.0*tmp[1], 2.0*tmp[2] };

        /* 物理要素へ: curl N = (1/detJ) * J * curl_hat
           curl_i = (1/detJ) * sum_a (∂x_i/∂ξ_a) * curl_hat_a
                  = (1/detJ) * sum_a J[a][i] * curl_hat[a]
        */
        for (int phys = 0; phys < 3; ++phys) {
            curl_N_edge[e][phys] =
                ( J[0][phys]*curl_hat[0] +
                  J[1][phys]*curl_hat[1] +
                  J[2][phys]*curl_hat[2] ) * inv_detJ;
        }
    }
}

static void tet2_W_ij(
    int i,
    int j,
    const double lam[4],
    const double glam[4][3],
    double W[3],
    double curlW[3]
){
    for(int d = 0; d < 3; ++d){
        W[d] = lam[i]*glam[j][d] - lam[j]*glam[i][d];
    }

    double c[3];
    cross3_2nd(glam[i], glam[j], c);

    for(int d = 0; d < 3; ++d){
        curlW[d] = 2.0 * c[d];
    }
}

static void curl_scalar_times_vec_2nd(
    const double gradf[3],
    double f,
    const double V[3],
    const double curlV[3],
    double out[3]
){
    double tmp[3];
    cross3_2nd(gradf, V, tmp);

    for(int d = 0; d < 3; ++d){
        out[d] = tmp[d] + f * curlV[d];
    }
}

void BBFE_std_shapefunc_tet2nd_nedelec_get_ref(
    const double xi[3],
    const int elem_global_node_id[4],
    double N[20][3],
    double curlN[20][3]
){
    const double r = xi[0];
    const double s = xi[1];
    const double t = xi[2];

    /*
     * Reference Tet:
     *   lambda0 = 1 - r - s - t
     *   lambda1 = r
     *   lambda2 = s
     *   lambda3 = t
     */
    double lam[4] = {
        1.0 - r - s - t,
        r,
        s,
        t
    };

    /*
     * grad lambda on reference element.
     * These are indexed by local vertex id 0..3.
     */
    const double glam[4][3] = {
        {-1.0, -1.0, -1.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0}
    };

    for(int a = 0; a < TET2_NDOF; ++a){
        for(int d = 0; d < 3; ++d){
            N[a][d]     = 0.0;
            curlN[a][d] = 0.0;
        }
    }

    /*
     * Edge DOF:
     *   local layout:
     *     ed = 0..5
     *     mode 0: W_ij
     *     mode 1: (lambda_i - lambda_j) W_ij
     *
     * Important:
     *   elem_global_node_id is used ONLY for orientation sorting.
     *   i and j remain local vertex ids, so lam[i] and glam[i] are valid.
     */
    for(int ed = 0; ed < TET2_NEDGE; ++ed){
        int i = tet_edge_conn[ed][0];
        int j = tet_edge_conn[ed][1];

        /*
         * Make edge direction consistent with original global node ids:
         *   smaller global node id -> larger global node id
         *
         * Do not use fe->conn[e] here if it may contain partition-local ids.
         */
        if(elem_global_node_id[i] > elem_global_node_id[j]){
            swap_int_2nd(&i, &j);
        }

        double W[3], curlW[3];

        tet2_W_ij(
            i,
            j,
            lam,
            glam,
            W,
            curlW
        );

        const int id0 = TET2_EDGE_MODE * ed;
        const int id1 = TET2_EDGE_MODE * ed + 1;

        /*
         * mode 0: W_ij
         */
        for(int d = 0; d < 3; ++d){
            N[id0][d]     = W[d];
            curlN[id0][d] = curlW[d];
        }

        /*
         * mode 1: (lambda_i - lambda_j) W_ij
         */
        const double f = lam[i] - lam[j];
        double gradf[3];

        for(int d = 0; d < 3; ++d){
            gradf[d] = glam[i][d] - glam[j][d];
            N[id1][d] = f * W[d];
        }

        curl_scalar_times_vec_2nd(
            gradf,
            f,
            W,
            curlW,
            curlN[id1]
        );
    }

    /*
     * Face DOF:
     *   local layout:
     *     fc = 0..3
     *     mode 0: lambda_i W_jk
     *     mode 1: lambda_j W_ki
     *
     * Important:
     *   The face vertices are sorted by original global node ids.
     *   The sorted values v[0], v[1], v[2] are still local vertex ids.
     */
    for(int fc = 0; fc < TET2_NFACE; ++fc){
        int v[3] = {
            tet_face_conn[fc][0],
            tet_face_conn[fc][1],
            tet_face_conn[fc][2]
        };

        /*
         * Sort only by original global node ids.
         * Never use elem_global_node_id as an index into lam/glam.
         */
        for(int p = 0; p < 3; ++p){
            for(int q = p + 1; q < 3; ++q){
                if(elem_global_node_id[v[p]] > elem_global_node_id[v[q]]){
                    swap_int_2nd(&v[p], &v[q]);
                }
            }
        }

        const int i = v[0];
        const int j = v[1];
        const int k = v[2];

        const int id0 = TET2_FACE_OFFSET + TET2_FACE_MODE * fc;
        const int id1 = TET2_FACE_OFFSET + TET2_FACE_MODE * fc + 1;

        /*
         * mode 0: lambda_i W_jk
         */
        {
            double W[3], curlW[3];

            tet2_W_ij(
                j,
                k,
                lam,
                glam,
                W,
                curlW
            );

            for(int d = 0; d < 3; ++d){
                N[id0][d] = lam[i] * W[d];
            }

            curl_scalar_times_vec_2nd(
                glam[i],
                lam[i],
                W,
                curlW,
                curlN[id0]
            );
        }

        /*
         * mode 1: lambda_j W_ki
         */
        {
            double W[3], curlW[3];

            tet2_W_ij(
                k,
                i,
                lam,
                glam,
                W,
                curlW
            );

            for(int d = 0; d < 3; ++d){
                N[id1][d] = lam[j] * W[d];
            }

            curl_scalar_times_vec_2nd(
                glam[j],
                lam[j],
                W,
                curlW,
                curlN[id1]
            );
        }
    }
}


void BBFE_std_shapefunc_tet2nd_nedelec_get_val_curl(
    const double xi[3],
    const int elem_global_node_id[4],
    const double J[3][3],
    const double J_inv[3][3],
    const double detJ,
    double **N_phys,
    double **curlN_phys
){
    double N_ref[20][3];
    double curlN_ref[20][3];

    /*
     * elem_global_node_id is used inside get_ref()
     * only to sort edge/face orientation.
     */
    BBFE_std_shapefunc_tet2nd_nedelec_get_ref(
        xi,
        elem_global_node_id,
        N_ref,
        curlN_ref
    );

    for(int a = 0; a < TET2_NDOF; ++a){
        for(int d = 0; d < 3; ++d){

            /*
             * Covariant Piola transform for H(curl):
             *   N_phys = J^{-T} N_ref
             *
             * This follows the existing tet1st implementation convention,
             * where J_inv[d][ref] is used.
             */
            N_phys[a][d] =
                  J_inv[d][0] * N_ref[a][0]
                + J_inv[d][1] * N_ref[a][1]
                + J_inv[d][2] * N_ref[a][2];

            /*
             * Curl transform:
             *   curl N_phys = (1 / detJ) J curl N_ref
             */
            curlN_phys[a][d] =
                (
                  J[0][d] * curlN_ref[a][0]
                + J[1][d] * curlN_ref[a][1]
                + J[2][d] * curlN_ref[a][2]
                ) / detJ;
        }
    }
}