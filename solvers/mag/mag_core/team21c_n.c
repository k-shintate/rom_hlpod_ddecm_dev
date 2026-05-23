#include "team7.h"
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>

/* ---------------- material / frequency constants ---------------- */

const double Sigma_p21_plate = 1.3889e6;             /* TEAM Problem 21a-2 non-magnetic steel plate [S/m] */
const double Sigma_air       = 0.0;

const double mu0_team21a0   = 4.0 * M_PI * 1.0e-7; /* [H/m] */
const double freq_team21a0  = 50.0;                /* [Hz] */
const double omega_team21a0 = 2.0 * M_PI * 50.0;   /* [rad/s] */

/* TEAM benchmark rated current */
static const double NI_RMS_team21a0 = 3000.0;      /* TEAM Problem 21a-2: source current [A.turn rms] */

/* small gauge stabilization for phi in non-conducting region */
static const double TEAM7_PHI_STAB = 1.0e-12;

/* =========================================================material accessors for A-phi========================================================= */

static inline double dot3_team21a0(const double a[3], const double b[3]){return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];}


static inline int is_team21a0_coil_prop(int prop){
    /* Physical IDs: 1 = COIL_1, 2 = COIL_2 */
    return (prop == 1 || prop == 2);
}

void get_material_for_prop_team21a0_Aphi(int prop,double* nu,double* sigma){

    *nu = 1.0 / mu0_team21a0;

    /* Physical groups used by the current Problem21a-2 mesh/input:
       1: COIL_1, 2: COIL_2, 3: NONMAG_PLATE, 4: AIR */
    if(prop == 3){
        /* NONMAG_PLATE */
        *sigma = Sigma_p21_plate;
        //*sigma = 0.0;
    }
    else if(prop == 1 || prop == 2){
        /* COIL_1 / COIL_2: impressed current source regions. */
        //*sigma = 0.0;
        *sigma = Sigma_p21_plate*1.0e-4;
    }
    else {
        /* AIR and any other non-conducting region */
        *sigma = 0.0;
    }

}

/* Representative coil information.Center chosen from the meshed TEAM-7 geometry:outer x : 94 ... 294 mmouter y :  0 ... 200 mmz-center: zCoil0 + 50 mm = 19 + 30 + 50 mm*/

static inline int get_coil_info_team21a0(int elem_prop, COIL_INFO* info){
    if(!is_team21a0_coil_prop(elem_prop)) return 0;

    /* TEAM Problem 21a-2 racetrack coils.
       Geometry is taken from the attached .geo:
         outer 270 mm x 270 mm, inner 200 mm x 200 mm,
         outer corner radius 45 mm, inner corner radius 10 mm,
         coil height 217 mm, gap 24 mm.
       The current path is in the x-y plane; the two coils carry opposite currents. */
    info->axis[0] = 0.0;
    info->axis[1] = 0.0;
    info->axis[2] = 1.0;

    info->turns = 1.0;

    /* Cross-section normal to the impressed current:
       radial build (35 mm) x coil height (217 mm). */
    info->area = (35.0 * 217.0) * (MM_TO_M * MM_TO_M);

    info->center[0] = 0.0;
    info->center[1] = 0.0;
    if(elem_prop == 1){
        /* COIL_1: lower coil */
        info->center[2] = -0.5 * 24.0 * MM_TO_M - 0.5 * 217.0 * MM_TO_M;
    } else {
        /* COIL_2: upper coil */
        info->center[2] =  0.5 * 24.0 * MM_TO_M + 0.5 * 217.0 * MM_TO_M;
    }

    return 1;
}

static inline double get_coil_current_team21a0(int prop){
    if(is_team21a0_coil_prop(prop)){
        return NI_RMS_team21a0;
    }
    return 0.0;
}

/* =========================================================Known source current density Js in coil region========================================================= */

typedef enum {TEAM7_SRC_PP,TEAM7_SRC_DSC} TEAM7_SRC_TYPE;

typedef struct {TEAM7_SRC_TYPE type;

/* common: real-valued amplitude for current density */
double value;

/* Parallelepiped */
double p0[3];
double p1[3];
double p2[3];
double p3[3];
double vec_val[3];

/* DoubleSectorialCylinder / annular sector */
double base[3];
int axis;              /* 0:x, 1:y, 2:z */
double height;
double start_deg;
double sweep_deg;
double r_in;
double r_out;

} TEAM7_SOURCE_DEF;

static int point_in_parallelepiped_source(const TEAM7_SOURCE_DEF* s,const double x[3],double J[3]){double b[3], c[3], d[3], rhs[3];double M[3][3];double det, u, v, w;

for(int i = 0; i < 3; ++i){
    b[i]   = s->p1[i] - s->p0[i];
    c[i]   = s->p2[i] - s->p0[i];
    d[i]   = s->p3[i] - s->p0[i];
    rhs[i] = x[i]     - s->p0[i];
}

M[0][0]=b[0]; M[0][1]=c[0]; M[0][2]=d[0];
M[1][0]=b[1]; M[1][1]=c[1]; M[1][2]=d[1];
M[2][0]=b[2]; M[2][1]=c[2]; M[2][2]=d[2];

det =
    M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) -
    M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) +
    M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);

if(fabs(det) < 1.0e-20) return 0;

u =
    (rhs[0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) -
     M[0][1]*(rhs[1]*M[2][2]-M[1][2]*rhs[2]) +
     M[0][2]*(rhs[1]*M[2][1]-M[1][1]*rhs[2])) / det;

v =
    (M[0][0]*(rhs[1]*M[2][2]-M[1][2]*rhs[2]) -
     rhs[0]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) +
     M[0][2]*(M[1][0]*rhs[2]-rhs[1]*M[2][0])) / det;

w =
    (M[0][0]*(M[1][1]*rhs[2]-rhs[1]*M[2][1]) -
     M[0][1]*(M[1][0]*rhs[2]-rhs[1]*M[2][0]) +
     rhs[0]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])) / det;

const double tol = 1.0e-10;

if(u < -tol || u > 1.0 + tol) return 0;
if(v < -tol || v > 1.0 + tol) return 0;
if(w < -tol || w > 1.0 + tol) return 0;

J[0] = s->vec_val[0];
J[1] = s->vec_val[1];
J[2] = s->vec_val[2];

return 1;

}

static double normalize_angle_rad_team21a0(double a){while(a >  M_PI) a -= 2.0 * M_PI;while(a < -M_PI) a += 2.0 * M_PI;return a;}

static int angle_in_sweep_team21a0(double theta, double start, double sweep){theta = normalize_angle_rad_team21a0(theta);start = normalize_angle_rad_team21a0(start);

double rel = normalize_angle_rad_team21a0(theta - start);

if(sweep > 0.0){
    return (rel >= -1.0e-12 && rel <= sweep + 1.0e-12);
} else {
    return (rel <=  1.0e-12 && rel >= sweep - 1.0e-12);
}

}

static int point_in_dsc_source(const TEAM7_SOURCE_DEF* s,const double x[3],double J[3]){const int ia = s->axis;const int ir = (ia + 1) % 3;const int is = (ia + 2) % 3;

double a0 = s->base[ia];
double a1 = s->base[ia] + s->height;

if(a0 > a1){
    double tmp = a0;
    a0 = a1;
    a1 = tmp;
}

if(x[ia] < a0 - 1.0e-12 || x[ia] > a1 + 1.0e-12) return 0;

const double dr = x[ir] - s->base[ir];
const double ds = x[is] - s->base[is];
const double r  = sqrt(dr*dr + ds*ds);

if(r < s->r_in  - 1.0e-12) return 0;
if(r > s->r_out + 1.0e-12) return 0;
if(r < 1.0e-20) return 0;

const double theta = atan2(ds, dr);
const double start = s->start_deg * M_PI / 180.0;
const double sweep = s->sweep_deg * M_PI / 180.0;

if(!angle_in_sweep_team21a0(theta, start, sweep)) return 0;

J[0] = 0.0;
J[1] = 0.0;
J[2] = 0.0;

/* Positive tangential direction in the local cross-section plane */
J[ir] = -s->value * ds / r;
J[is] =  s->value * dr / r;

return 1;

}

static int get_team21a0_known_Js_by_coord(int prop,const double x_ip[3],double Jmag,double Js[3]){
    Js[0] = 0.0;
    Js[1] = 0.0;
    Js[2] = 0.0;

    if(!is_team21a0_coil_prop(prop)) return 0;

    /* TEAM Problem 21a-2:
       prop 1 = COIL_1 (lower z), prop 2 = COIL_2 (upper z).
       The two coils carry source currents in opposite directions.
       If your measured B sign is reversed, swap these two signs. */
    const double coil_sign = (prop == 1) ? +1.0 : -1.0;
    const double J = coil_sign * Jmag;

    TEAM7_SOURCE_DEF srcs[8];
    int nsrc = 0;

    /* Geometry from Problem21a-2 .geo, all in metres. */
    const double ox0 = -0.135;
    const double ox1 =  0.135;
    const double oy0 = -0.135;
    const double oy1 =  0.135;

    const double ix0 = -0.100;
    const double ix1 =  0.100;
    const double iy0 = -0.100;
    const double iy1 =  0.100;

    const double rIn  = 0.010;
    const double rOut = 0.045;

    const double z0 = (prop == 1) ? -(0.024/2.0 + 0.217) : +(0.024/2.0);
    const double z1 = (prop == 1) ? -(0.024/2.0)         : +(0.024/2.0 + 0.217);
    const double h  = z1 - z0;

    /* The inner and outer rounded rectangles have the same corner centres. */
    const double cxBL = ix0 + rIn;  /* -0.090 */
    const double cyBL = iy0 + rIn;  /* -0.090 */
    const double cxBR = ix1 - rIn;  /* +0.090 */
    const double cyBR = iy0 + rIn;  /* -0.090 */
    const double cxTR = ix1 - rIn;  /* +0.090 */
    const double cyTR = iy1 - rIn;  /* +0.090 */
    const double cxTL = ix0 + rIn;  /* -0.090 */
    const double cyTL = iy1 - rIn;  /* +0.090 */

    /* bottom straight part: +x for prop 2, -x for prop 3 */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_PP,
        .value = 0.0,
        .p0 = {ox0 + rOut, oy0, z0},
        .p1 = {ox1 - rOut, oy0, z0},
        .p2 = {ox0 + rOut, iy0, z0},
        .p3 = {ox0 + rOut, oy0, z1},
        .vec_val = {+J, 0.0, 0.0}
    };

    /* right straight part: +y for prop 2, -y for prop 3 */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_PP,
        .value = 0.0,
        .p0 = {ix1, oy0 + rOut, z0},
        .p1 = {ox1, oy0 + rOut, z0},
        .p2 = {ix1, oy1 - rOut, z0},
        .p3 = {ix1, oy0 + rOut, z1},
        .vec_val = {0.0, +J, 0.0}
    };

    /* top straight part: -x for prop 2, +x for prop 3 */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_PP,
        .value = 0.0,
        .p0 = {ox0 + rOut, iy1, z0},
        .p1 = {ox1 - rOut, iy1, z0},
        .p2 = {ox0 + rOut, oy1, z0},
        .p3 = {ox0 + rOut, iy1, z1},
        .vec_val = {-J, 0.0, 0.0}
    };

    /* left straight part: -y for prop 2, +y for prop 3 */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_PP,
        .value = 0.0,
        .p0 = {ox0, oy0 + rOut, z0},
        .p1 = {ix0, oy0 + rOut, z0},
        .p2 = {ox0, oy1 - rOut, z0},
        .p3 = {ox0, oy0 + rOut, z1},
        .vec_val = {0.0, -J, 0.0}
    };

    /* bottom-left corner */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_DSC,
        .value = J,
        .base = {cxBL, cyBL, z0},
        .axis = 2,
        .height = h,
        .start_deg = 180.0,
        .sweep_deg = 90.0,
        .r_in = rIn,
        .r_out = rOut
    };

    /* bottom-right corner */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_DSC,
        .value = J,
        .base = {cxBR, cyBR, z0},
        .axis = 2,
        .height = h,
        .start_deg = -90.0,
        .sweep_deg = 90.0,
        .r_in = rIn,
        .r_out = rOut
    };

    /* top-right corner */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_DSC,
        .value = J,
        .base = {cxTR, cyTR, z0},
        .axis = 2,
        .height = h,
        .start_deg = 0.0,
        .sweep_deg = 90.0,
        .r_in = rIn,
        .r_out = rOut
    };

    /* top-left corner */
    srcs[nsrc++] = (TEAM7_SOURCE_DEF){
        .type = TEAM7_SRC_DSC,
        .value = J,
        .base = {cxTL, cyTL, z0},
        .axis = 2,
        .height = h,
        .start_deg = 90.0,
        .sweep_deg = 90.0,
        .r_in = rIn,
        .r_out = rOut
    };

    for(int k = 0; k < nsrc; ++k){
        int ok;

        if(srcs[k].type == TEAM7_SRC_PP){
            ok = point_in_parallelepiped_source(&srcs[k], x_ip, Js);
        } else {
            ok = point_in_dsc_source(&srcs[k], x_ip, Js);
        }

        if(ok) return 1;
    }

    Js[0] = 0.0;
    Js[1] = 0.0;
    Js[2] = 0.0;

    return 0;
}

/* =========================================================RHS assembly for A-phi========================================================= */

void set_element_vec_nedelec_Aphi_team21a0(MONOLIS*     monolis,BBFE_DATA*   fe,BBFE_BASIS*  basis,BBFE_BC*     bc,NEDELEC*     ned){
    (void)bc;

const int np = basis->num_integ_points;

double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
double* val_ip      = BB_std_calloc_1d_double(val_ip, np);

/* TEAM Problem 21a-2 coil cross-sectional area normal to current:
   radial build 35 mm x coil height 217 mm = 7595 mm^2. */
const double coil_area = (35.0 * 217.0) * 1.0e-6;

/* Problem 21a-2: 50 Hz sinusoidal RMS current. Keep the phase used by your solver convention. */
const double phase_rad = 90.0 * M_PI / 180.0;
//const double phase_rad = 0.0;
const double _Complex phase_factor = cos(phase_rad) + I * sin(phase_rad);

long long coil_elem_count = 0;
long long src_hit_count   = 0;

for(int e = 0; e < fe->total_num_elems; ++e){
    const int prop = ned->elem_prop[e];

    if(!is_team21a0_coil_prop(prop)) continue;

    coil_elem_count++;

    const double NI = get_coil_current_team21a0(prop);
    const double J_mag = NI / coil_area;

    BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

    for(int i = 0; i < ned->local_num_edges; ++i){
        const int gi = ned->nedelec_conn[e][i];
        const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

        for(int p = 0; p < np; ++p){
            double x_ip[3];
            double Js[3];

            get_interp_coords(e, p, fe, basis, x_ip);

            if(get_team21a0_known_Js_by_coord(prop, x_ip, J_mag, Js)){
                src_hit_count++;
                val_ip[p] = dot3_team21a0(Js, ned->N_edge[e][p][i]);
            } else {
                val_ip[p] = 0.0;
            }
        }

        {
            const double integ =
                BBFE_std_integ_calc(
                    np,
                    val_ip,
                    basis->integ_weight,
                    Jacobian_ip
                );

            /* Sign convention:
               Use += for curl(nu curl A) ... = Js.
               If your whole formulation uses the opposite RHS sign,
               change only this line consistently. */
            monolis->mat.C.B[gi] += (double)si * integ * phase_factor;
        }
    }
}

printf("[TEAM21a0 coil debug] coil_elem_count = %lld\n", coil_elem_count);
printf("[TEAM21a0 coil debug] src_hit_count   = %lld\n", src_hit_count);

BB_std_free_1d_double(Jacobian_ip, np);
BB_std_free_1d_double(val_ip, np);

}

/* =========================================================Matrix assembly for A-phi========================================================= */
void set_element_mat_nedelec_Aphi_team21a0(MONOLIS*     monolis,BBFE_DATA*   fe,BBFE_BASIS*  basis,BBFE_BC*     bc,NEDELEC*     ned){(void)bc;

const int np = basis->num_integ_points;

double* J_ip = BB_std_calloc_1d_double(J_ip, np);
double _Complex* val_ip_C = BB_std_calloc_1d_double_C(val_ip_C, np);

for(int e=0; e<fe->total_num_elems; ++e){

    int prop = ned->elem_prop[e];
    double nu;
    double sigma;
    get_material_for_prop_team21a0_Aphi(prop, &nu, &sigma);

    int nonlinear_mu = 0;
    (void)nonlinear_mu;

    BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

    /* curl-curl block (A-A) */
    for(int i=0; i<ned->local_num_edges; ++i){
        const int gi = ned->nedelec_conn[e][i];
        const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

        for(int j=0; j<ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int p=0; p<np; ++p){
                val_ip_C[p] = 0.0 + 0.0*I;
                val_ip_C[p] += nu * BBFE_elemmat_mag_mat_curl(
                    ned->curl_N_edge[e][p][i],
                    ned->curl_N_edge[e][p][j],
                    1.0
                );
            }

            double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
            v *= (double)(si*sj);

            monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
        }
    }

    /* sigma * iω * A mass block */
    for(int i=0; i<ned->local_num_edges; ++i){
        const int gi = ned->nedelec_conn[e][i];
        const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

        for(int j=0; j<ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int p=0; p<np; ++p){
                val_ip_C[p] = 0.0 + 0.0*I;
                val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                    ned->N_edge[e][p][i],
                    ned->N_edge[e][p][j],
                    sigma
                ) * omega_team21a0 * I;
            }

            double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
            v *= (double)(si*sj);
            monolis_add_scalar_to_sparse_matrix_C(monolis, gi, gj, 0, 0, v);
        }
    }

    /* phi-phi block */
    for(int m=0; m<fe->local_num_nodes; ++m){
        const int gm = fe->conn[e][m];

        for(int n=0; n<fe->local_num_nodes; ++n){
            const int gn = fe->conn[e][n];

            for(int p=0; p<np; ++p){
                val_ip_C[p] = 0.0 + 0.0*I;
                val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                    fe->geo[e][p].grad_N[m],
                    fe->geo[e][p].grad_N[n],
                    sigma
                ) * I * omega_team21a0;
            }

            double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
            monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gn, 0, 0, v);
        }
    }

    /* A-phi block */
    for(int j=0; j<ned->local_num_edges; ++j){
        const int gj = ned->nedelec_conn[e][j];
        const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

        for(int n=0; n<fe->local_num_nodes; ++n){
            const int gn = fe->conn[e][n];

            for(int p=0; p<np; ++p){
                val_ip_C[p] = 0.0 + 0.0*I;
                val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                    fe->geo[e][p].grad_N[n],
                    ned->N_edge[e][p][j],
                    sigma
                ) * I * omega_team21a0;
            }

            double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
            v *= (double)sj;

            monolis_add_scalar_to_sparse_matrix_C(monolis, gj, gn, 0, 0, v);
        }
    }

    /* phi-A block */
    for(int m=0; m<fe->local_num_nodes; ++m){
        const int gm = fe->conn[e][m];

        for(int i=0; i<ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int p=0; p<np; ++p){
                val_ip_C[p] = 0.0 + 0.0*I;
                val_ip_C[p] += BBFE_elemmat_mag_mat_mass(
                    ned->N_edge[e][p][i],
                    fe->geo[e][p].grad_N[m],
                    sigma
                ) * I * omega_team21a0;
            }

            double _Complex v = BBFE_std_integ_calc_C(np, val_ip_C, basis->integ_weight, J_ip);
            v *= (double)si;

            monolis_add_scalar_to_sparse_matrix_C(monolis, gm, gi, 0, 0, v);
        }
    }
}

BB_std_free_1d_double(J_ip, np);
BB_std_free_1d_double_C(val_ip_C, np);

}

/* =========================================================BC for A and phi========================================================= */

void apply_dirichlet_bc_for_A_and_phi_team21a0(MONOLIS* monolis,BBFE_DATA* fe,BBFE_BC* bc,NEDELEC* ned){const int nen = fe->local_num_nodes;

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

    if(prop == 4){
        for(int k=0; k<fe->local_num_nodes; ++k){
            int gn = fe->conn[e][k];
            node_is_conductor[gn] = 4;
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

for (int i = 0; i < num_nodes; ++i){
    if (node_is_conductor[i] == 1 || node_is_conductor[i] == 2 || node_is_conductor[i] == 4) {
    //if (node_is_conductor[i] == 1 || node_is_conductor[i] == 2 || node_is_conductor[i] == 3 || node_is_conductor[i] == 4) {

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