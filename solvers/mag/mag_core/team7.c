#include "team7.h"
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>

/* =========================================================
   TEAM Problem 7 : A-phi formulation

   Unknowns:
     - edge DOF : A   in H(curl)
     - node DOF : phi in H1

   Harmonic convention: exp(j*w*t)

   Governing form (typical A-phi form):

     curl( nu curl A )
     + j*w*sigma A
     + j*w*sigma grad(phi)
     = Js

     div[ j*w*sigma (A + grad(phi)) ] = 0

   Weak form used here:

     ∫ nu curl(A)·curl(vA)
   + j*w ∫ sigma A·vA
   + j*w ∫ sigma grad(phi)·vA
   + j*w ∫ sigma A·grad(vphi)
   + j*w ∫ sigma grad(phi)·grad(vphi)
   = ∫ Js·vA

   prop meaning in THIS file:
     1 : coil_inner
     2 : coil_limb
     3 : coil_corner
     4 : aluminum
     5 : hole
     6 : air
   ========================================================= */

/* ---------------- material / frequency constants ---------------- */

const double Sigma_al  = 3.526e7;                /* aluminum plate [S/m] */
const double Sigma_air = 0.0;

const double mu0_team7   = 4.0 * M_PI * 1.0e-7; /* [H/m] */
const double freq_team7  = 50.0;                /* [Hz] */
const double omega_team7 = 2.0 * M_PI * 50.0;   /* [rad/s] */

/* TEAM benchmark rated current */
static const double NI_RMS_team7 = 2742.0;      /* [A.turn rms] */

/* small gauge stabilization for phi in non-conducting region */
static const double TEAM7_PHI_STAB = 1.0e-12;

/* =========================================================
   helpers
   ========================================================= */

static inline int is_team7_coil_prop(int prop)
{
    return (prop == 2 || prop == 3);
}

static inline double dot3_team7(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* =========================================================
   material accessors for A-phi
   ========================================================= */

void get_material_for_prop_team7_Aphi(
    int prop,
    double* nu,
    double* sigma
){
    *nu = 1.0 / mu0_team7;

    if(prop == 4){
        *sigma = Sigma_al;
    }
    else if(is_team7_coil_prop(prop)){
        /* coil is impressed-current region; keep as nonconducting unless
           you explicitly want stranded-conductor conductivity included */
        *sigma = Sigma_al*1.0e-4;
    }
    else {
        *sigma = 0.0;
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
    info->area = 2500.0 * (MM_TO_M * MM_TO_M) / 40.0;

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
   Known source current density Js in coil region
   ========================================================= */


   static int get_team7_known_Js_by_coord(
    int prop,
    const double x_ip[3],
    double Jmag,
    double Js[3]
){
    const double MM = MM_TO_M;

    /* geometry from .geo */
    const double aluZ        = 19.0  * MM;
    const double coilGapZ    = 30.0  * MM;
    const double coilHeightZ = 100.0 * MM;

    const double zCoil0 = aluZ + coilGapZ;          /* 49 mm */
    const double zCoil1 = zCoil0 + coilHeightZ;     /* 149 mm */

    const double rOuter = 50.0 * MM;
    const double rInner = 25.0 * MM;

    const double ox0 =  94.0 * MM;
    const double ox1 = 294.0 * MM;
    const double oy0 =   0.0 * MM;
    const double oy1 = 200.0 * MM;

    const double ix0 = 119.0 * MM;
    const double ix1 = 269.0 * MM;
    const double iy0 =  25.0 * MM;
    const double iy1 = 175.0 * MM;

    /* corner centers from .geo */
    const double cxBL = ox0 + rOuter;   /* 144 mm */
    const double cyBL = oy0 + rOuter;   /*  50 mm */
    const double cxBR = ox1 - rOuter;   /* 244 mm */
    const double cyBR = oy0 + rOuter;   /*  50 mm */
    const double cxTR = ox1 - rOuter;   /* 244 mm */
    const double cyTR = oy1 - rOuter;   /* 150 mm */
    const double cxTL = ox0 + rOuter;   /* 144 mm */
    const double cyTL = oy1 - rOuter;   /* 150 mm */

    const double x = x_ip[0];
    const double y = x_ip[1];
    const double z = x_ip[2];

    Js[0] = 0.0;
    Js[1] = 0.0;
    Js[2] = 0.0;

    /* only coil regions */
    if(!(prop == 1 || prop == 2 || prop == 3)) return 0;

    /* z-range of the single extruded coil in the current .geo */
    if(z < zCoil0 || z > zCoil1) return 0;


    /* -----------------------------------------------------
       coil_limb : 4 straight outer limbs from .geo
       ----------------------------------------------------- */
    if(prop == 2){
        /* bottom limb */
        if(x >= ix0 + rInner && x <= ix1 - rInner &&
           y >= oy0 && y <= iy0){
            Js[0] =  Jmag; Js[1] = 0.0;  Js[2] = 0.0;
            return 1;
        }

        /* right limb */
        if(x >= ix1 && x <= ox1 &&
           y >= iy0 + rInner && y <= iy1 - rInner){
            Js[0] = 0.0;  Js[1] =  Jmag; Js[2] = 0.0;
            return 1;
        }

        /* top limb */
        if(x >= ix0 + rInner && x <= ix1 - rInner &&
           y >= iy1 && y <= oy1){
            Js[0] = -Jmag; Js[1] = 0.0;  Js[2] = 0.0;
            return 1;
        }

        /* left limb */
        if(x >= ox0 && x <= ix0 &&
           y >= iy0 + rInner && y <= iy1 - rInner){
            Js[0] = 0.0;  Js[1] = -Jmag; Js[2] = 0.0;
            return 1;
        }

        return 0;
    }

    /* -----------------------------------------------------
       coil_corner : 4 quarter-annuli from .geo
       ----------------------------------------------------- */
    if(prop == 3){
        struct Corner {
            double cx, cy;
            int sx, sy;
        } cs[4] = {
            { cxBL, cyBL, -1, -1 },
            { cxBR, cyBR, +1, -1 },
            { cxTR, cyTR, +1, +1 },
            { cxTL, cyTL, -1, +1 }
        };

        for(int k = 0; k < 4; ++k){
            const double dx = x - cs[k].cx;
            const double dy = y - cs[k].cy;
            const double r  = sqrt(dx*dx + dy*dy);

            if(r < rInner || r > rOuter) continue;
            if(cs[k].sx * dx < 0.0) continue;
            if(cs[k].sy * dy < 0.0) continue;

            /* CCW tangential direction */
            Js[0] = -Jmag * dy / r;
            Js[1] =  Jmag * dx / r;
            Js[2] =  0.0;
            return 1;
        }
    }

    return 0;
}

/* ===================

/* =========================================================
   Matrix assembly for A-phi
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
        double nu;
        double sigma;
        get_material_for_prop_team7_Aphi(prop, &nu, &sigma);

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        /* A-A block: ∫ nu curl A · curl vA + j*w*sigma ∫ A·vA */
        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int j = 0; j < ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p = 0; p < np; ++p){
                    val_ip_C[p]  = 0.0 + 0.0*I;
                    val_ip_C[p] += nu * BBFE_elemmat_mag_mat_curl(
                        ned->curl_N_edge[e][p][i],
                        ned->curl_N_edge[e][p][j],
                        1.0
                    );
                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7
                                 * BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        ned->N_edge[e][p][j],
                        sigma
                    );
                }

                {
                    double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                    v *= (double)(si * sj);
                    monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
                }
            }
        }

        /* phi-phi block: j*w*sigma ∫ grad(phi)·grad(vphi) */
        for(int m = 0; m < fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p = 0; p < np; ++p){
                    val_ip_C[p]  = 0.0 + 0.0*I;
                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7
                                 * BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[m],
                        fe->geo[e][p].grad_N[n],
                        sigma
                    );
                }

                {
                    double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                    monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
                }
            }
        }

        /* A-phi block: j*w*sigma ∫ grad(phi)·vA */
        for(int j = 0; j < ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p = 0; p < np; ++p){
                    val_ip_C[p]  = 0.0 + 0.0*I;
                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7
                                 * BBFE_elemmat_mag_mat_mass(
                        fe->geo[e][p].grad_N[n],
                        ned->N_edge[e][p][j],
                        sigma
                    );
                }

                {
                    double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                    v *= (double)sj;
                    monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
                }
            }
        }

        /* phi-A block: j*w*sigma ∫ A·grad(vphi) */
        for(int m = 0; m < fe->local_num_nodes; ++m){
            const int gm = fe->conn[e][m];

            for(int i = 0; i < ned->local_num_edges; ++i){
                const int gi = ned->nedelec_conn[e][i];
                const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

                for(int p = 0; p < np; ++p){
                    val_ip_C[p]  = 0.0 + 0.0*I;
                    val_ip_C[p] += (0.0 + 1.0*I) * omega_team7
                                 * BBFE_elemmat_mag_mat_mass(
                        ned->N_edge[e][p][i],
                        fe->geo[e][p].grad_N[m],
                        sigma
                    );
                }

                {
                    double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
                    v *= (double)si;
                    monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gi, 0, 0, v);
                }
            }
        }
    }

    BB_std_free_1d_double(J_ip, np);
    BB_std_free_1d_double_C(val_ip_C, np);
}

/* =========================================================
   RHS assembly for A-phi
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
    double* val_ip = BB_std_calloc_1d_double(val_ip, np);

    for(int e = 0; e < fe->total_num_elems; ++e){
        const int prop = ned->elem_prop[e];

        COIL_INFO coil;
        const int is_coil = get_coil_info_team7(prop, &coil);
        if(!is_coil) continue;

        const double Ii = get_coil_current_team7(prop);
        const double J_mag = (coil.turns * Ii) / coil.area;

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int p = 0; p < np; ++p){
                double x_ip[3];
                double Js[3] = {0.0, 0.0, 0.0};
                get_interp_coords(e, p, fe, basis, x_ip);

                if(get_team7_known_Js_by_coord(prop, x_ip, J_mag, Js)){
                    val_ip[p] = dot3_team7(Js, ned->N_edge[e][p][i]);
                }
                else {
                    val_ip[p] = 0.0;
                }
            }

            {
                const double integ = BBFE_std_integ_calc(np, val_ip, basis->integ_weight, Jacobian_ip);
                monolis->mat.C.B[gi] -= (double)si * integ;
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip, np);
}

/* =========================================================
   BC for A and phi
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
    build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);

    /* --------------------------------------------------------
       A (edge) boundary condition: A_tan = 0 on outer boundary
       -------------------------------------------------------- */
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
            int gn1 = fe->conn[e][edge_tbl[i][0]];
            int gn2 = fe->conn[e][edge_tbl[i][1]];
            int ged = ned->nedelec_conn[e][i];

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
    }

    
    int num_nodes = fe->total_num_nodes;
    int* node_is_conductor = (int*)calloc(num_nodes, sizeof(int));


    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop==6){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 6;
            }
        }
    }

    
    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop==5){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 5;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];
        if(prop==1){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 1;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 2){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 2; 
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 3){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 3; 
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        int prop = ned->elem_prop[e];

        if(prop == 4){
            for(int k=0; k<fe->local_num_nodes; ++k){
                int gn = fe->conn[e][k];
                node_is_conductor[gn] = 4; 
            }
        }
    }

    for (int i = 0; i < num_nodes; ++i){
        if (node_is_conductor[i] == 1||node_is_conductor[i] == 2||node_is_conductor[i] == 3||node_is_conductor[i] == 5||node_is_conductor[i] == 6) {
        //if (node_is_conductor[i] == 2||node_is_conductor[i] == 4 ||node_is_conductor[i] == 3) {
	    //if (node_is_conductor[i] == 4) {

	    monolis_set_Dirichlet_bc_C(
                monolis, 
                monolis->mat.C.B, 
                i, 
                0, 
                0.0 + 0.0*I
            );
        }  
    }
    

    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);
}