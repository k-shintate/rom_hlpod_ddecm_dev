#include "team7.h"
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>

/* =========================================================
   TEAM Problem 7 : T-Phi formulation inspired by
   NGSolve/TEAM-problems/TEAM-7/team7TPhiT0_withMS.py

   Unknowns:
     - edge DOF : T   in H(curl)
     - node DOF : Phi in H1

   Weak form (conceptual):
     ∫ rho curl(T)·curl(vT)
   + j*w*mu ∫ (T - grad Phi)·(vT - grad vPhi)
   = -j*w*mu*J0*t ∫ T0·(vT - grad vPhi) on coil

 prop meaning in THIS file:
     1 : coil_inner
     2 : coil_limb
     3 : coil_corner
     4 : aluminum
     5 : hole
     6 : air

   Notes:
   - The mesh is assumed to be already split geometrically into
     coil_inner / coil_limb / coil_corner.
   - Therefore, we do NOT classify subregions from coordinates anymore.
   - T0 is selected by prop and evaluated with simple local formulas
     only inside the corresponding already-meshed subregion.
   ========================================================= */

/* ---------------- material / frequency constants ---------------- */

const double Sigma_al  = 3.526e7;   /* aluminum plate [S/m] */
const double Sigma_air = 0.0;

const double mu0_team7   = 4.0 * M_PI * 1.0e-7;   /* [H/m] */
const double freq_team7  = 50.0;                  /* [Hz] */
const double omega_team7 = 2.0 * M_PI * 50.0;     /* [rad/s] */

/* TEAM benchmark rated current */
static const double NI_RMS_team7 = 2742.0;        /* [A.turn rms] */

/* Source scaling in the NGSolve code is:
     -j*w*mu*J0*0.025*T0*(...)
   where 0.025 is the radial thickness factor used there.
*/
static const double TEAM7_T0_THICKNESS = 0.025;   /* [m] */

static const double TEAM7_PHI_STAB = 1.0e-1;

/* =========================================================
   helpers
   ========================================================= */

static inline int is_team7_coil_prop(int prop)
{
    return (prop == 1||prop == 2 || prop == 3);
}

static inline double clamp01_team7(double x)
{
    if(x < 0.0) return 0.0;
    if(x > 1.0) return 1.0;
    return x;
}

static inline double dot3_team7(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* =========================================================
   material accessors for T-Phi
   ========================================================= */

void get_material_for_prop_team7_Tphi(
    int prop,
    double* rho,   /* resistivity-like coefficient in curl-curl term */
    double* mu     /* permeability in mass/coupling terms           */
){
    *mu = mu0_team7;

    if(prop == 4){
        /* aluminum plate */
        *rho = 1.0 / Sigma_al;
    }
    else if(prop == 5){
        *rho = 0.0;
    }
    else {
        /* coil + air: no conduction curl-curl coefficient */
        *rho = 0.0;
    }
}

/* Representative coil information.
   Center chosen from the meshed TEAM-7 geometry:
     outer x : 94 ... 294 mm
     outer y :  0 ... 200 mm
     z-center: zCoil0 + 50 mm = 19 + 30 + 50 mm
*/
static inline int get_coil_info_team7(int elem_prop, COIL_INFO* info)
{
    if(!is_team7_coil_prop(elem_prop)) return 0;

    info->axis[0] = 0.0;
    info->axis[1] = 0.0;
    info->axis[2] = 1.0;

    info->turns = 1.0;

    /* full coil cross-sectional area used for J0 = NI / area */
    info->area = 2500.0 * (MM_TO_M * MM_TO_M) / 40;
    //info->area = 1;

    info->center[0] = 0.5 * (94.0 + 294.0) * MM_TO_M;
    info->center[1] = 0.5 * (0.0  + 200.0) * MM_TO_M;
    info->center[2] = (19.0 + 30.0 + 50.0) * MM_TO_M;

    return 1;
}

static inline double get_coil_current_team7(int prop)
{
    if(is_team7_coil_prop(prop)){
        return NI_RMS_team7;
    }
    return 0.0;
}
/* =========================================================
   T0 construction based on coordinates

   T0 is a tangential vector field following the coil path
   in the xy-plane, with scalar amplitude varying linearly
   across the conductor thickness.

   We use a counterclockwise tangential direction:
     - top    : (-1, 0)
     - right  : ( 0, 1)
     - bottom : ( 1, 0)
     - left   : ( 0,-1)
     - corner : tangent of the local circular arc
                t = (-yr/r, xr/r)
   ========================================================= */


static double get_team7_T0_scalar_from_prop(
    const int prop,
    const double x_ip[3],
    const COIL_INFO* coil
){
    const double MM = MM_TO_M;

    /* coil center based local coordinates */
    const double x = x_ip[0] - coil->center[0];
    const double y = x_ip[1] - coil->center[1];

    /* geometry consistent with the provided .geo */
    const double rOuter = 50.0 * MM;
    const double rInner = 25.0 * MM;
    const double thickness = rOuter - rInner;   /* 25 mm */

    /* inner rectangle half sizes */
    const double inner_hx = 75.0 * MM;
    const double inner_hy = 75.0 * MM;

    /* straight limb half length */
    const double limb_half_len = 50.0 * MM;

    /* quarter-annulus centers = corners of the center rectangle */
    const double corner_hx = 50.0 * MM;
    const double corner_hy = 50.0 * MM;

    if(prop == 1){
        /* coil_inner */
        if(x >= -inner_hx && x <= inner_hx &&
           y >= -inner_hy && y <= inner_hy){
            return 1.0;
        }
        return 0.0;
    }

    if(prop == 2){
        /* top limb: x in [-50, 50], y in [75, 100] */
        if(x >= -limb_half_len && x <= limb_half_len &&
           y >=  inner_hy && y <= inner_hy + thickness){
            return clamp01_team7((inner_hy + thickness - y) / thickness);
        }

        /* bottom limb: x in [-50, 50], y in [-100, -75] */
        if(x >= -limb_half_len && x <= limb_half_len &&
           y >= -inner_hy - thickness && y <= -inner_hy){
            return clamp01_team7((y + inner_hy + thickness) / thickness);
        }

        /* right limb: x in [75, 100], y in [-50, 50] */
        if(x >= inner_hx && x <= inner_hx + thickness &&
           y >= -limb_half_len && y <= limb_half_len){
            return clamp01_team7((inner_hx + thickness - x) / thickness);
        }

        /* left limb: x in [-100, -75], y in [-50, 50] */
        if(x >= -inner_hx - thickness && x <= -inner_hx &&
           y >= -limb_half_len && y <= limb_half_len){
            return clamp01_team7((x + inner_hx + thickness) / thickness);
        }

        return 0.0;
    }

    if(prop == 3){
        /* top-right corner, center = ( 50,  50) mm */
        {
            const double cx =  corner_hx;
            const double cy =  corner_hy;
            const double xr = x - cx;
            const double yr = y - cy;
            const double r  = sqrt(xr*xr + yr*yr);

            if(x >= corner_hx && y >= corner_hy &&
               r >= rInner && r <= rOuter){
                return clamp01_team7((rOuter - r) / thickness);
            }
        }

        /* top-left corner, center = (-50,  50) mm */
        {
            const double cx = -corner_hx;
            const double cy =  corner_hy;
            const double xr = x - cx;
            const double yr = y - cy;
            const double r  = sqrt(xr*xr + yr*yr);

            if(x <= -corner_hx && y >= corner_hy &&
               r >= rInner && r <= rOuter){
                return clamp01_team7((rOuter - r) / thickness);
            }
        }

        /* bottom-left corner, center = (-50, -50) mm */
        {
            const double cx = -corner_hx;
            const double cy = -corner_hy;
            const double xr = x - cx;
            const double yr = y - cy;
            const double r  = sqrt(xr*xr + yr*yr);

            if(x <= -corner_hx && y <= -corner_hy &&
               r >= rInner && r <= rOuter){
                return clamp01_team7((rOuter - r) / thickness);
            }
        }

        /* bottom-right corner, center = ( 50, -50) mm */
        {
            const double cx =  corner_hx;
            const double cy = -corner_hy;
            const double xr = x - cx;
            const double yr = y - cy;
            const double r  = sqrt(xr*xr + yr*yr);

            if(x >= corner_hx && y <= -corner_hy &&
               r >= rInner && r <= rOuter){
                return clamp01_team7((rOuter - r) / thickness);
            }
        }

        return 0.0;
    }

    return 0.0;
}


static void get_team7_T0_vector(
    const int prop,
    const COIL_INFO* coil,
    const double x_ip[3],
    double T0[3]
){
    const double s = get_team7_T0_scalar_from_prop(prop, x_ip, coil);
    T0[0] = 0.0;
    T0[1] = 0.0;
    T0[2] = s;
}

/* =========================================================
   matrix assembly for T-Phi
   ========================================================= */

void set_element_mat_nedelec_Aphi_team7(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    (void)bc;

    const int np = basis->num_integ_points;

    double* J_ip = BB_std_calloc_1d_double(J_ip, np);
    double _Complex* val_ip_C = BB_std_calloc_1d_double_C(val_ip_C, np);

    for(int e = 0; e < fe->total_num_elems; ++e){

        const int prop = ned->elem_prop[e];

        double rho = 0.0, mu = mu0_team7;
        get_material_for_prop_team7_Tphi(prop, &rho, &mu);

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j = 0; j < ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p = 0; p < np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;

                    if(rho != 0.0){
                        val_ip_C[p] += rho * BBFE_elemmat_mag_mat_curl(
                            ned->curl_N_edge[e][p][i],
                            ned->curl_N_edge[e][p][j],
                            1.0
                        );
                    }

                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7 * mu
                        * BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            ned->N_edge[e][p][j],
                            1.0
                        );
                }

                double _Complex v = BBFE_std_integ_calc_C(
                    np, val_ip_C, basis->integ_weight, J_ip
                );
                v *= (double)(si * sj);

                monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
            }
        }

    if(prop==4)
    {
        for(int j = 0; j < ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p = 0; p < np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += -(0.0 + 1.0*I) * omega_team7 * mu
                        * BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][j],
                            fe->geo[e][p].grad_N[n],
                            1.0
                        );
                }

                double _Complex v = BBFE_std_integ_calc_C(
                    np, val_ip_C, basis->integ_weight, J_ip
                );
                v *= (double)sj;

                monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
            }
        }

        for(int m = 0; m < fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int i = 0; i < ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int p = 0; p < np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;
                    val_ip_C[p] += -(0.0 + 1.0*I) * omega_team7 * mu
                        * BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[m],
                            ned->N_edge[e][p][i],
                            1.0
                        );
                }

                double _Complex v = BBFE_std_integ_calc_C(
                    np, val_ip_C, basis->integ_weight, J_ip
                );
                v *= (double)si;

                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gi, 0, 0, v);
            }
        }
    }

        for(int m = 0; m < fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p = 0; p < np; ++p){
                    val_ip_C[p] = 0.0 + 0.0*I;

                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7 * mu
                        * BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[m],
                            fe->geo[e][p].grad_N[n],
                            1.0
                        );

                    if(prop == 4){
                        val_ip_C[p] += TEAM7_PHI_STAB
                            * BBFE_elemmat_convdiff_mat_mass(
                                basis->N[p][m],
                                basis->N[p][n],
                                1.0
                            );
                    }
                }

                double _Complex v = BBFE_std_integ_calc_C(
                    np, val_ip_C, basis->integ_weight, J_ip
                );

                monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
            }
        }

    }

    BB_std_free_1d_double(J_ip, np);
    BB_std_free_1d_double_C(val_ip_C, np);
}

/* =========================================================
   RHS assembly:
     -j*w*mu*J0*t ∫ T0·vT
     +j*w*mu*J0*t ∫ T0·grad(vPhi)
   on coil region only
   ========================================================= */

void set_element_vec_nedelec_Aphi_team7(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    BBFE_BC*     bc,
    NEDELEC*     ned)
{
    (void)bc;

    const int np = basis->num_integ_points;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double _Complex* val_ip_C = BB_std_calloc_1d_double_C(val_ip_C, np);

    for(int e = 0; e < fe->total_num_elems; ++e){

        const int prop = ned->elem_prop[e];

        COIL_INFO coil;
        const int is_coil = get_coil_info_team7(prop, &coil);
        if(!is_coil) continue;

        double rho = 0.0, mu = mu0_team7;
        get_material_for_prop_team7_Tphi(prop, &rho, &mu);
        (void)rho;

        const double NI_t = get_coil_current_team7(prop);
        const double J0   = NI_t / fmax(coil.area, 1.0e-30);
        const double _Complex rhs_scale =
            -(0.0 + 1.0*I) * omega_team7 * mu * J0 * TEAM7_T0_THICKNESS;

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* edge RHS: -j*w*mu*J0*t ∫ T0·Ni */
        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int p = 0; p < np; ++p){
                double x_ip[3];
                double T0[3];
                get_interp_coords(e, p, fe, basis, x_ip);
                get_team7_T0_vector(prop, &coil, x_ip, T0);

                val_ip_C[p] = rhs_scale * dot3_team7(T0, ned->N_edge[e][p][i]);
            }

            double _Complex integ = BBFE_std_integ_calc_C(
                np, val_ip_C, basis->integ_weight, Jacobian_ip
            );

            monolis->mat.C.B[gi] -= (double)si * integ;
        }

        /* node RHS: +j*w*mu*J0*t ∫ T0·grad(Nm)
           because rhs has -j*w*mu*J0*t * (vT - grad(vPhi))
         */
        for(int m = 0; m < fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int p = 0; p < np; ++p){
                double x_ip[3];
                double T0[3];
                get_interp_coords(e, p, fe, basis, x_ip);
                get_team7_T0_vector(prop, &coil, x_ip, T0);

                val_ip_C[p] = -rhs_scale * dot3_team7(T0, fe->geo[e][p].grad_N[m]);
            }

            double _Complex integ = BBFE_std_integ_calc_C(
                np, val_ip_C, basis->integ_weight, Jacobian_ip
            );

            monolis->mat.C.B[gm] -= integ;
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double_C(val_ip_C, np);
}

/* =========================================================
   Boundary conditions

   Edge BC:
     T_tan = 0 on outer boundary

   Node BC:
     Fix Phi only at ONE reference node inside the conducting plate
     to remove null space.
   ========================================================= */

void apply_dirichlet_bc_for_A_and_phi_team7(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BC* bc,
    NEDELEC* ned)
{
    const int nen = fe->local_num_nodes;

    int is_dir_edge_n = fe->total_num_nodes;
    bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
    build_dirichlet_edge_mask_from_boundary_faces_tet(
        fe, bc, ned, is_dir_edge, is_dir_edge_n
    );

    /* edge Dirichlet: T_tan = 0 on outer boundary */
    int n_local_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if(nen == 4){
        n_local_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if(nen == 8){
        n_local_edges = 12;
        edge_tbl = hex_edge_conn;
    }

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int i = 0; i < n_local_edges; ++i){
            const int gn1 = fe->conn[e][edge_tbl[i][0]];
            const int gn2 = fe->conn[e][edge_tbl[i][1]];
            const int ged = ned->nedelec_conn[e][i];

            if(!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;
            if(!is_dir_edge[ged]) continue;

            monolis_set_Dirichlet_bc_C(
                monolis,
                monolis->mat.C.B,
                ged,
                0,
                0.0 + 0.0*I
            );
        }

        for(int i = 0; i < fe->local_num_nodes; ++i){
            const int gn1 = fe->conn[e][i];
            if(!(bc->D_bc_exists[gn1])) continue;
            monolis_set_Dirichlet_bc_C(
                monolis,
                monolis->mat.C.B,
                gn1,
                0,
                0.0 + 0.0*I
            );
            //printf("add phi bc");
        }

        const int prop = ned->elem_prop[e];
        for(int i = 0; i < n_local_edges; ++i){
            const int ged = ned->nedelec_conn[e][i];
            if(prop==1||prop==2||prop==3 ||prop==5||prop==6){
                monolis_set_Dirichlet_bc_C(
                    monolis,
                    monolis->mat.C.B,
                    ged,
                    0,
                    0.0 + 0.0*I
                );
            }
        }


    }

    /* node Dirichlet for Phi:
       pin ONE node in the conducting plate (prop==4) */
    /*
    if(monolis_mpi_get_global_my_rank()==7){
        int ref_node = -1;

        for(int e = 0; e < fe->total_num_elems && ref_node < 0; ++e){
            const int prop = ned->elem_prop[e];
            if(prop != 4) continue;

            for(int k = 0; k < fe->local_num_nodes; ++k){
                ref_node = fe->conn[e][k];
                break;
            }
        }

        if(ref_node >= 0){
            printf("\n\nadd D_bc\n\n");
            monolis_set_Dirichlet_bc_C(
                monolis,
                monolis->mat.C.B,
                ref_node,
                0,
                0.0 + 0.0*I
            );
        }
    }
    */

    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);
}

double calc_aluminum_loss_team7_freq(
    FE_SYSTEM* sys,
    const double _Complex* x_c,
    double freq_hz,
    const char* directory
){
    (void)directory;   /* not needed if fe->total_num_elems is valid */
    (void)freq_hz;     /* loss formula below does not explicitly need omega */

    BBFE_DATA*  fe    = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC*    ned   = &(sys->ned);

    const int np = basis->num_integ_points;
    const int ne = fe->total_num_elems;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

    double loss_local = 0.0;

    for(int e = 0; e < ne; ++e){
        /* aluminum only */
        if(ned->elem_prop[e] != 4) continue;

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double _Complex Jp[3] = {0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};

            /* J = curl(T_h) */
            for(int i = 0; i < ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);
                const double _Complex ti = (double)si * x_c[gi];

                Jp[0] += ti * ned->curl_N_edge[e][p][i][0];
                Jp[1] += ti * ned->curl_N_edge[e][p][i][1];
                Jp[2] += ti * ned->curl_N_edge[e][p][i][2];
            }

            const double J2 = creal(
                Jp[0] * conj(Jp[0]) +
                Jp[1] * conj(Jp[1]) +
                Jp[2] * conj(Jp[2])
            );

            /* average Joule loss density
               If x_c is a peak phasor:   p = 0.5 * |J|^2 / sigma
               If x_c is an RMS phasor:   p =       |J|^2 / sigma
            */
            val_ip[p] = 0.5 * (1.0 / Sigma_al) * J2;
        }

        loss_local += BBFE_std_integ_calc(
            np, val_ip, basis->integ_weight, Jacobian_ip
        );
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);

    double loss_global = loss_local;
    monolis_allreduce_R(1, &loss_global, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    return loss_global;
}
