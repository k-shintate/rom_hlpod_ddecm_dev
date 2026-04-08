/* ============================================================
 * TEAM 21a-2 per-turn current constraint
 * penalty-method version
 * ============================================================ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================
 * constants / macros
 * ============================================================ */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEAM21A2_NRAD             17
#define TEAM21A2_NAX              18
#define TEAM21A2_NVOID             6
#define TEAM21A2_NTURNS_PER_COIL 300
#define TEAM21A2_NCOILS            2
#define TEAM21A2_NTURNS_TOTAL    (TEAM21A2_NCOILS * TEAM21A2_NTURNS_PER_COIL)

#define TEAM21A2_COIL1_ID          1
#define TEAM21A2_COIL2_ID          2

#define TEAM21A2_MM               (1.0e-3)

#define TEAM21A2_OUTER_X          (270.0 * TEAM21A2_MM)
#define TEAM21A2_OUTER_Y          (270.0 * TEAM21A2_MM)
#define TEAM21A2_INNER_X          (200.0 * TEAM21A2_MM)
#define TEAM21A2_INNER_Y          (200.0 * TEAM21A2_MM)
#define TEAM21A2_ROUTER           ( 45.0 * TEAM21A2_MM)
#define TEAM21A2_RINNER           ( 10.0 * TEAM21A2_MM)

#define TEAM21A2_COIL_HEIGHT      (217.0 * TEAM21A2_MM)
#define TEAM21A2_COIL_GAP_Z       ( 24.0 * TEAM21A2_MM)

#define TEAM21A2_TURN_RADIAL      ( 2.0 * TEAM21A2_MM)
#define TEAM21A2_TURN_AXIAL       ( 6.7 * TEAM21A2_MM)

#define TEAM21A2_RADIAL_AVAIL     (0.5 * (TEAM21A2_OUTER_X - TEAM21A2_INNER_X))
#define TEAM21A2_AXIAL_AVAIL      (TEAM21A2_COIL_HEIGHT)

#define TEAM21A2_GAP_RAD \
    ((TEAM21A2_RADIAL_AVAIL - TEAM21A2_NRAD * TEAM21A2_TURN_RADIAL) / (TEAM21A2_NRAD + 1))

#define TEAM21A2_GAP_AX \
    ((TEAM21A2_AXIAL_AVAIL - TEAM21A2_NAX * TEAM21A2_TURN_AXIAL) / (TEAM21A2_NAX + 1))

#define TEAM21A2_COIL1_Z0 (-(TEAM21A2_COIL_GAP_Z/2.0 + TEAM21A2_COIL_HEIGHT))
#define TEAM21A2_COIL2_Z0 ( +(TEAM21A2_COIL_GAP_Z/2.0))

#define TEAM21A2_I_RMS    (10.0)
#define TEAM21A2_FREQ_HZ  (50.0)
#define TEAM21A2_SIGMA_CU (5.8e7)

/* ------------------------------------------------------------
 * prop numbering assumption
 *   coil1 turns : 101 .. 400
 *   coil2 turns : 401 .. 700
 * ------------------------------------------------------------ */
#define TEAM21A2_PROP_COIL1_BEGIN 101
#define TEAM21A2_PROP_COIL2_BEGIN 401

/* ============================================================
 * structures
 * ============================================================ */
typedef struct {
    int    global_turn_id;   /* 0..599 */
    int    coil_id;          /* 1 or 2 */
    int    turn_id_in_coil;  /* 0..299 */

    int    i_rad;
    int    j_ax;

    double d0;
    double z0;
    double h;

    double area;
    double length_center;
    double center[3];
} TURN_INFO;

typedef struct {
    int  num_elems;
    int* elem_to_turn_gid;   /* -1 for non-turn elems */
} TURN_ELEM_MAP;

typedef struct {
    int    n_turns_total;
    double alpha[TEAM21A2_NTURNS_TOTAL];
    double target_current[TEAM21A2_NTURNS_TOTAL];
    double current_value[TEAM21A2_NTURNS_TOTAL];
} TURN_CURRENT_PENALTY_INFO;

/* ============================================================
 * local utility
 * ============================================================ */
static inline double clamp_double(double x, double a, double b)
{
    return (x < a) ? a : ((x > b) ? b : x);
}

static inline double dot3_local(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* ============================================================
 * slot / turn table
 * ============================================================ */
int team21a2_is_void_slot(int i_rad, int j_ax)
{
    if((j_ax == 0 || j_ax == TEAM21A2_NAX - 1) && (i_rad <= 2)) return 1;
    return 0;
}

void team21a2_compute_turn_center(const TURN_INFO* turn, double xc[3])
{
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = turn->z0 + 0.5 * turn->h;
}

double team21a2_compute_turn_centerline_length(const TURN_INFO* turn)
{
    const double dmid = turn->d0 + 0.5 * TEAM21A2_TURN_RADIAL;

    const double ox0 = -0.5 * TEAM21A2_OUTER_X;
    const double ox1 =  0.5 * TEAM21A2_OUTER_X;
    const double oy0 = -0.5 * TEAM21A2_OUTER_Y;
    const double oy1 =  0.5 * TEAM21A2_OUTER_Y;

    const double x0 = ox0 + dmid;
    const double x1 = ox1 - dmid;
    const double y0 = oy0 + dmid;
    const double y1 = oy1 - dmid;
    const double rr = TEAM21A2_ROUTER - dmid;

    const double Lx = (x1 - rr) - (x0 + rr);
    const double Ly = (y1 - rr) - (y0 + rr);

    return 2.0 * (Lx + Ly) + 2.0 * M_PI * rr;
}

int team21a2_build_turn_table(TURN_INFO turns[TEAM21A2_NTURNS_TOTAL])
{
    int k1 = 0;
    int k2 = 0;

    for(int j = 0; j < TEAM21A2_NAX; ++j){
        const double zBase1 = TEAM21A2_COIL1_Z0
            + TEAM21A2_GAP_AX
            + j * (TEAM21A2_TURN_AXIAL + TEAM21A2_GAP_AX);

        const double zBase2 = TEAM21A2_COIL2_Z0
            + TEAM21A2_GAP_AX
            + j * (TEAM21A2_TURN_AXIAL + TEAM21A2_GAP_AX);

        for(int i = 0; i < TEAM21A2_NRAD; ++i){
            if(team21a2_is_void_slot(i, j)) continue;

            const double d0 = TEAM21A2_GAP_RAD
                + i * (TEAM21A2_TURN_RADIAL + TEAM21A2_GAP_RAD);

            {
                TURN_INFO* t = &turns[k1];
                t->global_turn_id  = k1;
                t->coil_id         = TEAM21A2_COIL1_ID;
                t->turn_id_in_coil = k1;
                t->i_rad           = i;
                t->j_ax            = j;
                t->d0              = d0;
                t->z0              = zBase1;
                t->h               = TEAM21A2_TURN_AXIAL;
                t->area            = TEAM21A2_TURN_RADIAL * TEAM21A2_TURN_AXIAL;
                t->length_center   = team21a2_compute_turn_centerline_length(t);
                team21a2_compute_turn_center(t, t->center);
                ++k1;
            }

            {
                const int gid = TEAM21A2_NTURNS_PER_COIL + k2;
                TURN_INFO* t = &turns[gid];
                t->global_turn_id  = gid;
                t->coil_id         = TEAM21A2_COIL2_ID;
                t->turn_id_in_coil = k2;
                t->i_rad           = i;
                t->j_ax            = j;
                t->d0              = d0;
                t->z0              = zBase2;
                t->h               = TEAM21A2_TURN_AXIAL;
                t->area            = TEAM21A2_TURN_RADIAL * TEAM21A2_TURN_AXIAL;
                t->length_center   = team21a2_compute_turn_centerline_length(t);
                team21a2_compute_turn_center(t, t->center);
                ++k2;
            }
        }
    }

    if(k1 != TEAM21A2_NTURNS_PER_COIL || k2 != TEAM21A2_NTURNS_PER_COIL){
        fprintf(stderr,
            "[team21a2_build_turn_table] ERROR: k1=%d k2=%d\n", k1, k2);
        return 0;
    }

    return 1;
}

/* ============================================================
 * current target
 * ============================================================ */
double team21a2_get_turn_current_rms(int coil_id)
{
    if(coil_id == TEAM21A2_COIL1_ID) return +TEAM21A2_I_RMS;
    if(coil_id == TEAM21A2_COIL2_ID) return -TEAM21A2_I_RMS;
    return 0.0;
}

double team21a2_get_turn_current_peak(int coil_id, double t, int use_time_harmonic)
{
    const double Iamp  = sqrt(2.0) * TEAM21A2_I_RMS;
    const double omega = 2.0 * M_PI * TEAM21A2_FREQ_HZ;

    if(use_time_harmonic){
        if(coil_id == TEAM21A2_COIL1_ID) return +Iamp * sin(omega * t);
        if(coil_id == TEAM21A2_COIL2_ID) return -Iamp * sin(omega * t);
    }else{
        if(coil_id == TEAM21A2_COIL1_ID) return +Iamp;
        if(coil_id == TEAM21A2_COIL2_ID) return -Iamp;
    }
    return 0.0;
}

/* ============================================================
 * turn tangent
 * ============================================================ */
int team21a2_get_turn_tangent(
    const TURN_INFO* turn,
    const double x_ip[3],
    double tdir[3])
{
    const double ox0 = -0.5 * TEAM21A2_OUTER_X;
    const double ox1 =  0.5 * TEAM21A2_OUTER_X;
    const double oy0 = -0.5 * TEAM21A2_OUTER_Y;
    const double oy1 =  0.5 * TEAM21A2_OUTER_Y;

    const double dmid = turn->d0 + 0.5 * TEAM21A2_TURN_RADIAL;

    const double x0 = ox0 + dmid;
    const double x1 = ox1 - dmid;
    const double y0 = oy0 + dmid;
    const double y1 = oy1 - dmid;
    const double rr = TEAM21A2_ROUTER - dmid;

    const double xs = x1 - rr;
    const double ys = y1 - rr;

    const double x = x_ip[0];
    const double y = x_ip[1];

    double best_d2 = 1.0e300;
    double best_tx = 0.0, best_ty = 0.0;

    #define UPDATE_BEST(px, py, tx_, ty_) do { \
        const double dx_ = x - (px); \
        const double dy_ = y - (py); \
        const double d2_ = dx_*dx_ + dy_*dy_; \
        if(d2_ < best_d2){ \
            best_d2 = d2_; \
            best_tx = (tx_); \
            best_ty = (ty_); \
        } \
    } while(0)

    {
        const double px = clamp_double(x, -xs, xs);
        const double py = y1;
        UPDATE_BEST(px, py, -1.0, 0.0);
    }

    {
        const double px = clamp_double(x, -xs, xs);
        const double py = y0;
        UPDATE_BEST(px, py, +1.0, 0.0);
    }

    {
        const double px = x1;
        const double py = clamp_double(y, -ys, ys);
        UPDATE_BEST(px, py, 0.0, +1.0);
    }

    {
        const double px = x0;
        const double py = clamp_double(y, -ys, ys);
        UPDATE_BEST(px, py, 0.0, -1.0);
    }

    {
        const double cxy[4][2] = {
            { +xs, +ys },
            { -xs, +ys },
            { -xs, -ys },
            { +xs, -ys }
        };

        for(int k = 0; k < 4; ++k){
            const double cx = cxy[k][0];
            const double cy = cxy[k][1];
            const double vx = x - cx;
            const double vy = y - cy;
            const double vn = sqrt(vx*vx + vy*vy);
            if(vn < 1.0e-20) continue;

            const double ux = vx / vn;
            const double uy = vy / vn;

            int ok = 0;
            if(k == 0) ok = (ux >= 0.0 && uy >= 0.0);
            if(k == 1) ok = (ux <= 0.0 && uy >= 0.0);
            if(k == 2) ok = (ux <= 0.0 && uy <= 0.0);
            if(k == 3) ok = (ux >= 0.0 && uy <= 0.0);
            if(!ok) continue;

            const double px = cx + rr * ux;
            const double py = cy + rr * uy;

            const double tx = -uy;
            const double ty =  ux;
            UPDATE_BEST(px, py, tx, ty);
        }
    }

    #undef UPDATE_BEST

    {
        const double nt = sqrt(best_tx*best_tx + best_ty*best_ty);
        if(nt < 1.0e-20){
            tdir[0] = 0.0;
            tdir[1] = 0.0;
            tdir[2] = 0.0;
            return 0;
        }
        tdir[0] = best_tx / nt;
        tdir[1] = best_ty / nt;
        tdir[2] = 0.0;
        return 1;
    }
}

/* ============================================================
 * prop -> turn map
 * ============================================================ */
static inline int team21a2_prop_is_turn(int prop)
{
    return ((TEAM21A2_PROP_COIL1_BEGIN <= prop &&
             prop < TEAM21A2_PROP_COIL1_BEGIN + TEAM21A2_NTURNS_PER_COIL) ||
            (TEAM21A2_PROP_COIL2_BEGIN <= prop &&
             prop < TEAM21A2_PROP_COIL2_BEGIN + TEAM21A2_NTURNS_PER_COIL));
}

static inline int team21a2_prop_to_turn_gid(int prop)
{
    if(TEAM21A2_PROP_COIL1_BEGIN <= prop &&
       prop < TEAM21A2_PROP_COIL1_BEGIN + TEAM21A2_NTURNS_PER_COIL){
        return prop - TEAM21A2_PROP_COIL1_BEGIN;
    }

    if(TEAM21A2_PROP_COIL2_BEGIN <= prop &&
       prop < TEAM21A2_PROP_COIL2_BEGIN + TEAM21A2_NTURNS_PER_COIL){
        return TEAM21A2_NTURNS_PER_COIL + (prop - TEAM21A2_PROP_COIL2_BEGIN);
    }

    return -1;
}

void team21a2_build_elem_to_turn_map(
    const BBFE_DATA* fe,
    const NEDELEC* ned,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    TURN_ELEM_MAP* map)
{
    (void)turns;

    map->num_elems = fe->total_num_elems;
    map->elem_to_turn_gid = (int*)calloc((size_t)map->num_elems, sizeof(int));

    for(int e = 0; e < map->num_elems; ++e){
        const int prop = ned->elem_prop[e];
        map->elem_to_turn_gid[e] = team21a2_prop_to_turn_gid(prop);
    }
}

/* ============================================================
 * penalty info
 * ============================================================ */
void team21a2_init_turn_penalty_info(
    TURN_CURRENT_PENALTY_INFO* pen,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    double alpha_default,
    double current_time,
    int use_time_harmonic)
{
    pen->n_turns_total = TEAM21A2_NTURNS_TOTAL;

    for(int k = 0; k < TEAM21A2_NTURNS_TOTAL; ++k){
        pen->alpha[k] = alpha_default;
        pen->target_current[k] = team21a2_get_turn_current_peak(
            turns[k].coil_id, current_time, use_time_harmonic);
        pen->current_value[k] = 0.0;
    }
}

/* ============================================================
 * evaluate current per turn:
 *   I_k(x) = (1/L_k) * ∫_Ωk J·t dV
 * with
 *   J = sigma * ( -(A_curr - A_prev)/dt - grad(phi) )
 * ============================================================ */
void team21a2_eval_turn_currents(
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    const TURN_ELEM_MAP* map,
    TURN_CURRENT_PENALTY_INFO* pen)
{
    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    for(int k = 0; k < TEAM21A2_NTURNS_TOTAL; ++k){
        pen->current_value[k] = 0.0;
    }

    for(int e = 0; e < fe->total_num_elems; ++e){
        const int gid_turn = map->elem_to_turn_gid[e];
        if(gid_turn < 0) continue;

        const TURN_INFO* turn = &turns[gid_turn];
        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double x_ip[3], tdir[3];
            double A_curr[3], A_prev[3], grad_phi[3], Jc[3];

            get_interp_coords(e, p, fe, basis, x_ip);

            if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                val_ip_C[p] = 0.0;
                continue;
            }

            compute_A_ip(ned, e, p, x_curr, fe->total_num_nodes, A_curr);
            compute_A_ip(ned, e, p, x_prev, fe->total_num_nodes, A_prev);
            compute_grad_phi_ip(fe, e, p, x_curr, grad_phi);

            Jc[0] = TEAM21A2_SIGMA_CU * (-(A_curr[0] - A_prev[0]) * inv_dt - grad_phi[0]);
            Jc[1] = TEAM21A2_SIGMA_CU * (-(A_curr[1] - A_prev[1]) * inv_dt - grad_phi[1]);
            Jc[2] = TEAM21A2_SIGMA_CU * (-(A_curr[2] - A_prev[2]) * inv_dt - grad_phi[2]);

            val_ip_C[p] = dot3_local(Jc, tdir) / turn->length_center;
        }

        pen->current_value[gid_turn] +=
            BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

/* ============================================================
 * penalty residual:
 *   R += alpha_k * g_k * dI_k/dx
 * where
 *   g_k = I_k(x) - I_target
 * ============================================================ */
void assemble_turn_current_penalty_residual_team21a2(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    const TURN_ELEM_MAP* map,
    const TURN_CURRENT_PENALTY_INFO* pen,
    double dt)
{
    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_A       = BB_std_calloc_1d_double(val_A, np);
    double* val_P       = BB_std_calloc_1d_double(val_P, np);

    for(int e = 0; e < fe->total_num_elems; ++e){
        const int gid_turn = map->elem_to_turn_gid[e];
        if(gid_turn < 0) continue;

        const TURN_INFO* turn = &turns[gid_turn];
        const double alpha = pen->alpha[gid_turn];
        const double gk = pen->current_value[gid_turn] - pen->target_current[gid_turn];

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        for(int j = 0; j < ned->local_num_edges; ++j){
            const int gj = ned->nedelec_conn[e][j];
            const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

            for(int p = 0; p < np; ++p){
                double x_ip[3], tdir[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                    val_A[p] = 0.0;
                    continue;
                }

                val_A[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                         * inv_dt
                         * dot3_local(ned->N_edge[e][p][j], tdir);
            }

            {
                const double dIk_dAj =
                    BBFE_std_integ_calc(np, val_A, basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[gj] += (double)sj * alpha * gk * dIk_dAj;
            }
        }

        for(int n = 0; n < fe->local_num_nodes; ++n){
            const int gn = fe->conn[e][n];

            for(int p = 0; p < np; ++p){
                double x_ip[3], tdir[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                    val_P[p] = 0.0;
                    continue;
                }

                val_P[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                         * dot3_local(fe->geo[e][p].grad_N[n], tdir);
            }

            {
                const double dIk_dPhin =
                    BBFE_std_integ_calc(np, val_P, basis->integ_weight, Jacobian_ip);

                monolis->mat.R.B[gn] += alpha * gk * dIk_dPhin;
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_A, np);
    BB_std_free_1d_double(val_P, np);
}

/* ============================================================
 * penalty Jacobian (Gauss-Newton approximation):
 *   K += alpha_k * (dI_k/dx)^T (dI_k/dx)
 * ============================================================ */
void assemble_turn_current_penalty_jacobian_team21a2(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    const TURN_ELEM_MAP* map,
    const TURN_CURRENT_PENALTY_INFO* pen,
    double dt)
{
    const int np = basis->num_integ_points;
    const double inv_dt = 1.0 / dt;

    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_iA      = BB_std_calloc_1d_double(val_iA, np);
    double* val_jA      = BB_std_calloc_1d_double(val_jA, np);
    double* val_nP      = BB_std_calloc_1d_double(val_nP, np);
    double* val_mP      = BB_std_calloc_1d_double(val_mP, np);

    for(int e = 0; e < fe->total_num_elems; ++e){
        const int gid_turn = map->elem_to_turn_gid[e];
        if(gid_turn < 0) continue;

        const TURN_INFO* turn = &turns[gid_turn];
        const double alpha = pen->alpha[gid_turn];

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* A-A */
        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int p = 0; p < np; ++p){
                double x_ip[3], tdir[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                    val_iA[p] = 0.0;
                    continue;
                }

                val_iA[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                          * inv_dt
                          * dot3_local(ned->N_edge[e][p][i], tdir);
            }

            const double dIk_dAi =
                BBFE_std_integ_calc(np, val_iA, basis->integ_weight, Jacobian_ip);

            for(int j = 0; j < ned->local_num_edges; ++j){
                const int gj = ned->nedelec_conn[e][j];
                const int sj = (ned->edge_sign ? ned->edge_sign[e][j] : 1);

                for(int p = 0; p < np; ++p){
                    double x_ip[3], tdir[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                        val_jA[p] = 0.0;
                        continue;
                    }

                    val_jA[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                              * inv_dt
                              * dot3_local(ned->N_edge[e][p][j], tdir);
                }

                const double dIk_dAj =
                    BBFE_std_integ_calc(np, val_jA, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis, gi, gj, 0, 0,
                    (double)(si * sj) * alpha * dIk_dAi * dIk_dAj);
            }
        }

        /* phi-phi */
        for(int n = 0; n < fe->local_num_nodes; ++n){
            const int gn = fe->conn[e][n];

            for(int p = 0; p < np; ++p){
                double x_ip[3], tdir[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                    val_nP[p] = 0.0;
                    continue;
                }

                val_nP[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                          * dot3_local(fe->geo[e][p].grad_N[n], tdir);
            }

            const double dIk_dPhin =
                BBFE_std_integ_calc(np, val_nP, basis->integ_weight, Jacobian_ip);

            for(int m = 0; m < fe->local_num_nodes; ++m){
                const int gm = fe->conn[e][m];

                for(int p = 0; p < np; ++p){
                    double x_ip[3], tdir[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                        val_mP[p] = 0.0;
                        continue;
                    }

                    val_mP[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                              * dot3_local(fe->geo[e][p].grad_N[m], tdir);
                }

                const double dIk_dPhim =
                    BBFE_std_integ_calc(np, val_mP, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis, gn, gm, 0, 0,
                    alpha * dIk_dPhin * dIk_dPhim);
            }
        }

        /* A-phi and phi-A */
        for(int i = 0; i < ned->local_num_edges; ++i){
            const int gi = ned->nedelec_conn[e][i];
            const int si = (ned->edge_sign ? ned->edge_sign[e][i] : 1);

            for(int p = 0; p < np; ++p){
                double x_ip[3], tdir[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                    val_iA[p] = 0.0;
                    continue;
                }

                val_iA[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                          * inv_dt
                          * dot3_local(ned->N_edge[e][p][i], tdir);
            }

            const double dIk_dAi =
                BBFE_std_integ_calc(np, val_iA, basis->integ_weight, Jacobian_ip);

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];

                for(int p = 0; p < np; ++p){
                    double x_ip[3], tdir[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    if(!team21a2_get_turn_tangent(turn, x_ip, tdir)){
                        val_nP[p] = 0.0;
                        continue;
                    }

                    val_nP[p] = -(TEAM21A2_SIGMA_CU / turn->length_center)
                              * dot3_local(fe->geo[e][p].grad_N[n], tdir);
                }

                const double dIk_dPhin =
                    BBFE_std_integ_calc(np, val_nP, basis->integ_weight, Jacobian_ip);

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis, gi, gn, 0, 0,
                    (double)si * alpha * dIk_dAi * dIk_dPhin);

                monolis_add_scalar_to_sparse_matrix_R(
                    monolis, gn, gi, 0, 0,
                    (double)si * alpha * dIk_dPhin * dIk_dAi);
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_iA, np);
    BB_std_free_1d_double(val_jA, np);
    BB_std_free_1d_double(val_nP, np);
    BB_std_free_1d_double(val_mP, np);
}

/* ============================================================
 * example driver usage
 * ============================================================ */
void team21a2_apply_turn_current_penalty(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    int use_time_harmonic)
{
    static int initialized = 0;
    static TURN_INFO turns[TEAM21A2_NTURNS_TOTAL];
    static TURN_ELEM_MAP turn_map;
    static TURN_CURRENT_PENALTY_INFO pen;

    if(!initialized){
        team21a2_build_turn_table(turns);
        team21a2_build_elem_to_turn_map(fe, ned, turns, &turn_map);

        team21a2_init_turn_penalty_info(
            &pen,
            turns,
            1.0e8,             /* penalty coefficient: tune */
            current_time,
            use_time_harmonic);

        initialized = 1;
    }

    for(int k = 0; k < TEAM21A2_NTURNS_TOTAL; ++k){
        pen.target_current[k] = team21a2_get_turn_current_peak(
            turns[k].coil_id, current_time, use_time_harmonic);
    }

    team21a2_eval_turn_currents(
        fe, basis, ned, x_prev, x_curr, dt, turns, &turn_map, &pen);

    assemble_turn_current_penalty_jacobian_team21a2(
        monolis, fe, basis, ned, turns, &turn_map, &pen, dt);

    assemble_turn_current_penalty_residual_team21a2(
        monolis, fe, basis, ned, turns, &turn_map, &pen, dt);
}