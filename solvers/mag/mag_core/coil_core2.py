/* ============================================================
 * team21a2_turn_current.h
 * Per-turn total-current constraint for TEAM 21a-2
 * ============================================================ */
#ifndef TEAM21A2_TURN_CURRENT_H
#define TEAM21A2_TURN_CURRENT_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------- fixed benchmark settings ---------- */
#define TEAM21A2_NRAD            17
#define TEAM21A2_NAX             18
#define TEAM21A2_NVOID            6
#define TEAM21A2_NTURNS_PER_COIL 300
#define TEAM21A2_NCOILS            2
#define TEAM21A2_NTURNS_TOTAL   (TEAM21A2_NCOILS * TEAM21A2_NTURNS_PER_COIL)

#define TEAM21A2_COIL1_ID         1
#define TEAM21A2_COIL2_ID         2

/* ---------- geometry in SI ---------- */
#define TEAM21A2_MM              (1.0e-3)
#define TEAM21A2_OUTER_X         (270.0 * TEAM21A2_MM)
#define TEAM21A2_OUTER_Y         (270.0 * TEAM21A2_MM)
#define TEAM21A2_INNER_X         (200.0 * TEAM21A2_MM)
#define TEAM21A2_INNER_Y         (200.0 * TEAM21A2_MM)
#define TEAM21A2_ROUTER          ( 45.0 * TEAM21A2_MM)
#define TEAM21A2_RINNER          ( 10.0 * TEAM21A2_MM)

#define TEAM21A2_COIL_HEIGHT     (217.0 * TEAM21A2_MM)
#define TEAM21A2_COIL_GAP_Z      ( 24.0 * TEAM21A2_MM)

#define TEAM21A2_TURN_RADIAL     ( 2.0 * TEAM21A2_MM)
#define TEAM21A2_TURN_AXIAL      ( 6.7 * TEAM21A2_MM)

#define TEAM21A2_RADIAL_AVAIL    (0.5 * (TEAM21A2_OUTER_X - TEAM21A2_INNER_X))
#define TEAM21A2_AXIAL_AVAIL     (TEAM21A2_COIL_HEIGHT)

#define TEAM21A2_GAP_RAD \
    ((TEAM21A2_RADIAL_AVAIL - TEAM21A2_NRAD * TEAM21A2_TURN_RADIAL) / (TEAM21A2_NRAD + 1))
#define TEAM21A2_GAP_AX \
    ((TEAM21A2_AXIAL_AVAIL  - TEAM21A2_NAX  * TEAM21A2_TURN_AXIAL ) / (TEAM21A2_NAX  + 1))

#define TEAM21A2_COIL1_Z0 (-(TEAM21A2_COIL_GAP_Z/2.0 + TEAM21A2_COIL_HEIGHT))
#define TEAM21A2_COIL2_Z0 ( +(TEAM21A2_COIL_GAP_Z/2.0))

#define TEAM21A2_I_RMS     (10.0)
#define TEAM21A2_FREQ_HZ   (50.0)
#define TEAM21A2_SIGMA_CU  (5.8e7)

/* ---------- per-turn metadata ---------- */
typedef struct {
    int    global_turn_id;   /* 0..599 */
    int    coil_id;          /* 1 or 2 */
    int    turn_id_in_coil;  /* 0..299 */

    int    i_rad;            /* 0..16 candidate column */
    int    j_ax;             /* 0..17 candidate row */

    double d0;               /* inward offset of this turn outer boundary */
    double z0;               /* z start of this turn */
    double h;                /* turn axial height = 6.7 mm */

    double area;             /* conductor cross-section area */
    double length_center;    /* centerline length for diagnostics */

    double center[3];        /* rough geometric center */
} TURN_INFO;

/* ---------- extra unknown block for lambda ---------- */
typedef struct {
    int n_turns_total;       /* 600 */
    int row0_lambda;         /* first global row index of lambda block */
} TURN_CURRENT_SYSTEM_INFO;

/* ---------- map FE element -> turn ---------- */
typedef struct {
    int  num_elems;
    int* elem_to_turn_gid;   /* size = total_num_elems, -1 for non-turn elems */
} TURN_ELEM_MAP;

int team21a2_is_void_slot(int i_rad, int j_ax);
int team21a2_build_turn_table(TURN_INFO turns[TEAM21A2_NTURNS_TOTAL]);

double team21a2_get_turn_current_rms(int coil_id);
double team21a2_get_turn_current_peak(int coil_id, double t, int use_time_harmonic);

/* tangent on one resolved turn ring */
int team21a2_get_turn_tangent(
    const TURN_INFO* turn,
    const double x_ip[3],
    double tdir[3]);

/* utilities */
void team21a2_compute_turn_center(const TURN_INFO* turn, double xc[3]);
double team21a2_compute_turn_centerline_length(const TURN_INFO* turn);

/* FE mapping */
void team21a2_build_elem_to_turn_map(
    const BBFE_DATA* fe,
    const NEDELEC* ned,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    TURN_ELEM_MAP* map);

/* block-system helpers */
void team21a2_init_turn_system_info(
    const BBFE_DATA* fe,
    const NEDELEC* ned,
    TURN_CURRENT_SYSTEM_INFO* sys);

/* Jacobian/residual for constraints */
void assemble_turn_current_constraint_jacobian_team21a2(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    const TURN_ELEM_MAP* map,
    const TURN_CURRENT_SYSTEM_INFO* sys);

void assemble_turn_current_constraint_residual_team21a2(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,
    const double* x_curr,
    double dt,
    double current_time,
    const TURN_INFO turns[TEAM21A2_NTURNS_TOTAL],
    const TURN_ELEM_MAP* map,
    const TURN_CURRENT_SYSTEM_INFO* sys,
    int use_time_harmonic);

#endif


