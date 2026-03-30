
#include "shapefunc.h"

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

void cross3(const double a[3], const double b[3], double c[3])
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
        for (int phys = 0; phys < 3; ++phys) {
            N_edge[e][phys] =
                J_inv[0][phys]*hatN[0] +
                J_inv[1][phys]*hatN[1] +
                J_inv[2][phys]*hatN[2];
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
