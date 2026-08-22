#include "core_ROM.h"
//#include "core_Nitsche.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

const int BUFFER_SIZE = 1024;


/* TET4/TRI3 wall-mode configuration is owned by the core solver. */
#define T4_WALL_STRONG_NOSLIP              0
#define T4_WALL_NITSCHE_NOSLIP             1
#define T4_WALL_NITSCHE_FREESLIP           2
#define T4_WALL_NITSCHE_SPALDING           3
#define T4_WALL_STRONG_NORMAL_FREESLIP     4
#define T4_WALL_STRONG_NORMAL_SPALDING     5

#define T4_WALL_EXCHANGE_WALL_Q            0
#define T4_WALL_EXCHANGE_OPPOSITE_NODE     1

/*
 * These two runtime-parameter entry points are implemented with C linkage
 * in core_FOM.c, so keep C linkage only for them.
 */
#ifdef __cplusplus
extern "C" {
#endif

void BBFE_fluid_set_nitsche_runtime_parameters(
    double gamma_n,
    double dt_penalty_coeff);

void BBFE_fluid_set_tet4_wall_runtime_parameters(
    int wall_mode,
    int exchange_mode,
    double gamma_n,
    double dt_penalty_coeff,
    double kappa,
    double B,
    double ut_eps,
    double feature_angle_deg);

#ifdef __cplusplus
}
#endif


/*
 * core_FOM.h already declares these two functions with C++ linkage.
 * Define them here with the same linkage.
 *
 * Assumption required by the existing API:
 * surf->conn stores the same global node IDs used by the volume FE/BC data.
 */
int BBFE_fluid_sups_add_surface_velocity_Dirichlet(
    BBFE_BC* bc,
    const BBFE_DATA* surf,
    const double imposed_velocity[3])
{
    if (
        bc == NULL ||
        surf == NULL ||
        imposed_velocity == NULL ||
        bc->D_bc_exists == NULL ||
        bc->imposed_D_val == NULL
    ) {
        return -1;
    }

    if (
        bc->block_size < 3 ||
        bc->total_num_nodes <= 0
    ) {
        return -1;
    }

    for (int e = 0; e < surf->total_num_elems; ++e) {
        for (int a = 0; a < surf->local_num_nodes; ++a) {
            const int node_id = surf->conn[e][a];

            if (
                node_id < 0 ||
                node_id >= bc->total_num_nodes
            ) {
                fprintf(
                    stderr,
                    "%s ERROR: invalid surface node id %d "
                    "in BBFE_fluid_sups_add_surface_velocity_Dirichlet().\n",
                    CODENAME,
                    node_id);

                return -1;
            }

            for (int d = 0; d < 3; ++d) {
                const int index =
                    bc->block_size * node_id + d;

                /*
                 * A surface node can appear in several surface elements.
                 * Count a Dirichlet DOF only when it is newly constrained.
                 */
                if (!bc->D_bc_exists[index]) {
                    bc->num_D_bcs += 1;
                }

                bc->D_bc_exists[index] = true;
                bc->imposed_D_val[index] = imposed_velocity[d];
            }
        }
    }

    return 0;
}


void BBFE_fluid_sups_impose_surface_velocity_on_values(
    const BBFE_DATA* surf,
    double** v,
    const double imposed_velocity[3])
{
    if (
        surf == NULL ||
        v == NULL ||
        imposed_velocity == NULL
    ) {
        return;
    }

    for (int e = 0; e < surf->total_num_elems; ++e) {
        for (int a = 0; a < surf->local_num_nodes; ++a) {
            const int node_id = surf->conn[e][a];

            if (node_id < 0) {
                continue;
            }

            v[node_id][0] = imposed_velocity[0];
            v[node_id][1] = imposed_velocity[1];
            v[node_id][2] = imposed_velocity[2];
        }
    }
}

static double read_t4_nitsche_dt_penalty_coeff_from_env(void);

typedef struct T4WallOptions {
    int mode;
    int exchange_mode;
    double gamma_n;
    double c_dt;
    double kappa;
    double B;
    double ut_eps;
    double feature_angle_deg;
} T4WallOptions;

static const char* t4_wall_mode_name(int mode)
{
    switch (mode) {
    case T4_WALL_STRONG_NOSLIP:          return "strong-noslip";
    case T4_WALL_NITSCHE_NOSLIP:         return "nitsche-noslip";
    case T4_WALL_NITSCHE_FREESLIP:       return "nitsche-freeslip";
    case T4_WALL_NITSCHE_SPALDING:       return "nitsche-spalding";
    case T4_WALL_STRONG_NORMAL_FREESLIP: return "strong-normal-freeslip";
    case T4_WALL_STRONG_NORMAL_SPALDING: return "strong-normal-spalding";
    default:                              return "unknown";
    }
}

static int t4_parse_wall_mode_name(const char* text)
{
    if (text == NULL) return -1;
    if (strcmp(text, "strong-noslip") == 0)          return T4_WALL_STRONG_NOSLIP;
    if (strcmp(text, "nitsche-noslip") == 0)         return T4_WALL_NITSCHE_NOSLIP;
    if (strcmp(text, "nitsche-freeslip") == 0)       return T4_WALL_NITSCHE_FREESLIP;
    if (strcmp(text, "nitsche-spalding") == 0)       return T4_WALL_NITSCHE_SPALDING;
    if (strcmp(text, "strong-normal-freeslip") == 0 ||
        strcmp(text, "projected-strong-normal-freeslip") == 0)
        return T4_WALL_STRONG_NORMAL_FREESLIP;
    if (strcmp(text, "strong-normal-spalding") == 0 ||
        strcmp(text, "projected-strong-normal-spalding") == 0)
        return T4_WALL_STRONG_NORMAL_SPALDING;
    return -1;
}

static const char* t4_get_option_value(
    int argc,
    char* argv[],
    const char* name)
{
    const size_t nlen = strlen(name);
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], name) == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s ERROR: option %s requires a value.\n", CODENAME, name);
                exit(EXIT_FAILURE);
            }
            return argv[i + 1];
        }
        if (strncmp(argv[i], name, nlen) == 0 && argv[i][nlen] == '=') {
            return argv[i] + nlen + 1;
        }
    }
    return NULL;
}

static int t4_has_flag(int argc, char* argv[], const char* name)
{
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], name) == 0) return 1;
    }
    return 0;
}

static int t4_is_wall_option_name(const char* arg, int* consumes_next)
{
    static const char* names[] = {
        "--wall-mode",
        "--ahmed-wall-mode",
        "--wall-exchange",
        "--wall-gamma",
        "--wall-cdt",
        "--wall-kappa",
        "--wall-B",
        "--wall-ut-eps",
        "--wall-feature-angle"
    };

    if (consumes_next != NULL) *consumes_next = 0;
    if (arg == NULL) return 0;

    const int n = (int)(sizeof(names) / sizeof(names[0]));
    for (int i = 0; i < n; ++i) {
        const size_t len = strlen(names[i]);
        if (strcmp(arg, names[i]) == 0) {
            if (consumes_next != NULL) *consumes_next = 1;
            return 1;
        }
        if (strncmp(arg, names[i], len) == 0 && arg[len] == '=') {
            return 1;
        }
    }
    return 0;
}

/* Remove wall-specific options before passing argc/argv to the pre-existing
 * BBFE command-line parser.  This avoids relying on that parser to ignore
 * options it does not know about.  The strings themselves remain owned by the
 * original argv; only the pointer array is allocated here. */
static char** t4_build_framework_argv(int argc, char* argv[], int* framework_argc)
{
    char** out = (char**)calloc((size_t)argc + 1u, sizeof(char*));
    if (out == NULL) {
        fprintf(stderr, "%s ERROR: failed to allocate filtered argv.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    int nout = 0;
    if (argc > 0) out[nout++] = argv[0];

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--wall-help") == 0) continue;

        int consumes_next = 0;
        if (t4_is_wall_option_name(argv[i], &consumes_next)) {
            if (consumes_next && i + 1 < argc) ++i;
            continue;
        }

        out[nout++] = argv[i];
    }

    out[nout] = NULL;
    if (framework_argc != NULL) *framework_argc = nout;
    return out;
}

static void t4_print_wall_help_and_exit(void)
{
    printf(
        "TET4 Ahmed wall options:\n"
        "  --wall-mode MODE  (alias: --ahmed-wall-mode)\n"
        "    strong-noslip\n"
        "    nitsche-noslip\n"
        "    nitsche-freeslip\n"
        "    nitsche-spalding\n"
        "    strong-normal-freeslip\n"
        "    strong-normal-spalding\n"
        "  --wall-exchange opposite-node|wall-q\n"
        "  --wall-gamma VALUE\n"
        "  --wall-cdt VALUE\n"
        "  --wall-kappa VALUE\n"
        "  --wall-B VALUE\n"
        "  --wall-ut-eps VALUE\n"
        "  --wall-feature-angle DEG\n"
        "\n"
        "Default: --wall-mode nitsche-spalding --wall-exchange opposite-node\n"
        "Ground z=0 remains component-wise strong no-slip in this executable.\n"
        "Modes strong-normal-* use the experimental feature-aware projected-strong implementation.\n");
    exit(EXIT_SUCCESS);
}

static double t4_parse_double_option(
    int argc,
    char* argv[],
    const char* name,
    double default_value)
{
    const char* text = t4_get_option_value(argc, argv, name);
    if (text == NULL) return default_value;
    char* endptr = NULL;
    const double value = strtod(text, &endptr);
    if (endptr == text || endptr == NULL || *endptr != '\0' || !isfinite(value)) {
        fprintf(stderr, "%s ERROR: invalid %s=\"%s\".\n", CODENAME, name, text);
        exit(EXIT_FAILURE);
    }
    return value;
}

static T4WallOptions t4_read_wall_options(int argc, char* argv[])
{
    if (t4_has_flag(argc, argv, "--wall-help")) {
        t4_print_wall_help_and_exit();
    }

    T4WallOptions o;
    /* Keep the previous Case-C family as the default, but fix its exchange
     * velocity to the owner-TET opposite node. */
    o.mode = T4_WALL_NITSCHE_SPALDING;
    o.exchange_mode = T4_WALL_EXCHANGE_OPPOSITE_NODE;
    o.gamma_n = 100.0;
    o.c_dt = read_t4_nitsche_dt_penalty_coeff_from_env();
    o.kappa = 0.41;
    o.B = 5.2;
    o.ut_eps = 1.0e-12;
    o.feature_angle_deg = 30.0;

    const char* mode_text = t4_get_option_value(argc, argv, "--wall-mode");
    if (mode_text == NULL) {
        mode_text = t4_get_option_value(argc, argv, "--ahmed-wall-mode");
    }
    if (mode_text != NULL) {
        o.mode = t4_parse_wall_mode_name(mode_text);
        if (o.mode < 0) {
            fprintf(stderr,
                "%s ERROR: unsupported --wall-mode=%s.\n"
                "  allowed: strong-noslip, nitsche-noslip, nitsche-freeslip, "
                "nitsche-spalding, strong-normal-freeslip, strong-normal-spalding\n",
                CODENAME,
                mode_text);
            exit(EXIT_FAILURE);
        }
    }

    const char* exch = t4_get_option_value(argc, argv, "--wall-exchange");
    if (exch != NULL) {
        if (strcmp(exch, "opposite-node") == 0) {
            o.exchange_mode = T4_WALL_EXCHANGE_OPPOSITE_NODE;
        }
        else if (strcmp(exch, "wall-q") == 0) {
            o.exchange_mode = T4_WALL_EXCHANGE_WALL_Q;
        }
        else {
            fprintf(stderr, "%s ERROR: --wall-exchange must be opposite-node or wall-q.\n", CODENAME);
            exit(EXIT_FAILURE);
        }
    }

    o.gamma_n = t4_parse_double_option(argc, argv, "--wall-gamma", o.gamma_n);
    o.c_dt = t4_parse_double_option(argc, argv, "--wall-cdt", o.c_dt);
    o.kappa = t4_parse_double_option(argc, argv, "--wall-kappa", o.kappa);
    o.B = t4_parse_double_option(argc, argv, "--wall-B", o.B);
    o.ut_eps = t4_parse_double_option(argc, argv, "--wall-ut-eps", o.ut_eps);
    o.feature_angle_deg = t4_parse_double_option(
        argc, argv, "--wall-feature-angle", o.feature_angle_deg);

    if (o.gamma_n <= 0.0 || o.c_dt < 0.0 || o.kappa <= 0.0 || o.ut_eps <= 0.0 ||
        o.feature_angle_deg <= 0.0 || o.feature_angle_deg >= 90.0) {
        fprintf(stderr, "%s ERROR: invalid wall option range.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    return o;
}

#ifndef T4_AHMED_SPLIT_NUM_REGIONS
#define T4_AHMED_SPLIT_NUM_REGIONS 2
#endif

typedef struct T4AhmedBodySupportDiagnostics {
    double area[T4_AHMED_SPLIT_NUM_REGIONS];
    double projected_area[T4_AHMED_SPLIT_NUM_REGIONS];
    double int_n[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_pressure[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_molecular_viscous[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_normal_molecular_viscous[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_walllaw[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double force_normal_penalty_reaction[T4_AHMED_SPLIT_NUM_REGIONS][3];
    double face_count[T4_AHMED_SPLIT_NUM_REGIONS];
    double support_z_cut;
} T4AhmedBodySupportDiagnostics;

int calc_ahmed_body_support_diagnostics_tet4_tri3_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw_in[3],
    double gamma_n,
    double support_z_cut,
    const double eU_in[3],
    T4AhmedBodySupportDiagnostics* diag);

void ahmed_body_support_diagnostics_allreduce(
    T4AhmedBodySupportDiagnostics* diag,
    MONOLIS_COM* monolis_com);

void output_ahmed_body_support_drag_diagnostics(
    const T4AhmedBodySupportDiagnostics* diag,
    double t,
    double qA,
    double Aref,
    const double eU_in[3],
    const char* directory);

static double read_t4_nitsche_dt_penalty_coeff_from_env(void)
{
    const double default_value = 0.10;
    const char* text = getenv("T4_NITSCHE_DT_PENALTY_COEFF");

    if (text == NULL || text[0] == '\0') {
        return default_value;
    }

    char* endptr = NULL;
    const double value = strtod(text, &endptr);

    if (
        endptr == text ||
        endptr == NULL ||
        *endptr != '\0' ||
        !isfinite(value) ||
        value < 0.0
    ) {
        fprintf(
            stderr,
            "%s ERROR: invalid T4_NITSCHE_DT_PENALTY_COEFF=\"%s\"\n",
            CODENAME,
            text);
        exit(EXIT_FAILURE);
    }

    return value;
}


static double read_t4_ahmed_support_z_cut_from_env(void)
{
    const double default_value = 0.050;
    const char* text = getenv("T4_AHMED_SUPPORT_Z_CUT");

    if (text == NULL || text[0] == '\0') {
        return default_value;
    }

    char* endptr = NULL;
    const double value = strtod(text, &endptr);

    if (
        endptr == text ||
        endptr == NULL ||
        *endptr != '\0' ||
        !isfinite(value)
    ) {
        fprintf(
            stderr,
            "%s ERROR: invalid T4_AHMED_SUPPORT_Z_CUT=\"%s\"\n",
            CODENAME,
            text);
        exit(EXIT_FAILURE);
    }

    return value;
}

double hot_start_read_initialize_val(
    double*     int_val,
    const char* input_fname,
    const char* directory)
{
    int BUFFER_SIZE = 1024;
    int total_num_nodes;
    int ndof;
    double t = 0.0;
        FILE* fp;
        char fname[BUFFER_SIZE];
        char id[BUFFER_SIZE];

        fp = BBFE_sys_read_fopen(fp, "hot_start/start_time.dat", directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%lf", &(t));
    fclose(fp);

        fp = BBFE_sys_read_fopen(fp, input_fname, directory);
        fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes), &(ndof));
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fscanf(fp, "%lf", &(int_val[i * ndof + j]));
        }
    }
        fclose(fp);

    return t;
}

void hot_start_write_initialize_val(
    double*         int_val,
    const int       total_num_nodes,
    const int       ndof,
    const double    time,
    const char*     output_fname,
    const char*     directory)
{
        FILE* fp;
        char fname[BUFFER_SIZE];
        char id[BUFFER_SIZE];

        fp = BBFE_sys_write_fopen(fp, output_fname, directory);
        fprintf(fp, "initialization\n");
    fprintf(fp, "%d %d\n", total_num_nodes, ndof);
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fprintf(fp, "%e ", int_val[i * ndof + j]);
        }
        fprintf(fp, "\n");
    }
        fclose(fp);

    if(monolis_mpi_get_global_my_rank()==0){
        fp = BBFE_sys_write_fopen(fp, "hot_start/start_time.dat", directory);
        fprintf(fp, "start_time\n");
        fprintf(fp, "%lf", time);
        fclose(fp);
    }
}

static int t4_diag_trace_is_enabled(void)
{
    static int initialized = 0;
    static int enabled = 0;

    if (!initialized) {
        const char* env =
            getenv("T4_DIAG_TRACE");

        enabled =
            (
                env != NULL &&
                env[0] != '\0' &&
                strcmp(env, "0") != 0
            );

        initialized = 1;
    }

    return enabled;
}

static int t4_all_ranks_status_ok(
    int local_status,
    MONOLIS_COM* monolis_com,
    const char* label)
{
    if (monolis_com == NULL) {
        fprintf(
            stderr,
            "[rank %d][%s] monolis_com is NULL\n",
            monolis_mpi_get_global_my_rank(),
            label != NULL ? label : "diagnostic");

        return 0;
    }

    /*
     * local_status == 0 : success
     * local_status != 0 : failure
     */
    double failure_count =
        (local_status == 0)
        ? 0.0
        : 1.0;

    /*
     * 成功・失敗にかかわらず全rankが呼ぶ。
     */
    monolis_allreduce_R(
        1,
        &failure_count,
        MONOLIS_MPI_SUM,
        monolis_com->comm);

    /*
     * failure_countは全rankで同じ値になる。
     * したがって、この後の分岐も全rankで一致する。
     */
    if (failure_count > 0.5) {
        if (local_status != 0) {
            fprintf(
                stderr,
                "[rank %d][%s] local calculation failed: status=%d\n",
                monolis_mpi_get_global_my_rank(),
                label != NULL ? label : "diagnostic",
                local_status);
        }

        if (monolis_mpi_get_global_my_rank() == 0) {
            fprintf(
                stderr,
                "[%s] skipped because %.0f rank(s) failed\n",
                label != NULL ? label : "diagnostic",
                failure_count);
        }

        return 0;
    }

    return 1;
}

static void t4_diag_trace(
    const char* label,
    int step)
{
    if (!t4_diag_trace_is_enabled()) {
        return;
    }

    fprintf(
        stderr,
        "[rank %d] step=%d %s\n",
        monolis_mpi_get_global_my_rank(),
        step,
        label != NULL ? label : "(null)");

    fflush(stderr);
}

static int run_wall_and_continuity_diagnostics_parallel_safe(
    FE_SYSTEM* sys,
    int step,
    double t,
    double U_ref,
    const double Uw[3],
    const double eU[3],
    const double eP[3],
    double gamma_n,
    int write_nodal_yplus)
{
    if (sys == NULL) {
        fprintf(
            stderr,
            "[rank %d][diagnostics] sys is NULL\n",
            monolis_mpi_get_global_my_rank());

        return -1;
    }

    const int rank =
        monolis_mpi_get_global_my_rank();

    /*
     * Surface ownership policy:
     *   sys->surf          : overlap surface; use only for distributed Nitsche assembly.
     *   sys->surf_internal : non-duplicated/owned surface; use for diagnostics and forces.
     *
     * Keeping these aliases here prevents post-processing from accidentally
     * integrating the overlap surface more than once across MPI ranks.
     */
    const BBFE_DATA* wall_diag_surf =
        &(sys->surf_internal);

    const BBFE_BASIS* wall_diag_basis =
        &(sys->basis_surf_internal);

    int wall_mesh_global_ok = 0;
    int wall_velocity_global_ok = 0;
    int continuity_global_ok = 0;

    /*
     * 局所計算が失敗した場合にも、
     * 未初期化の値を使わないようにする。
     */
    T4WallMeshDiagnostics wall_mesh_diag;
    T4WallDiagnostics wall_velocity_diag;
    T4DomainContinuityDiagnostics continuity_diag;

    memset(
        &wall_mesh_diag,
        0,
        sizeof(wall_mesh_diag));

    memset(
        &wall_velocity_diag,
        0,
        sizeof(wall_velocity_diag));

    memset(
        &continuity_diag,
        0,
        sizeof(continuity_diag));

    /* ================================================================ */
    /* 1. Wall mesh / y+ diagnostics                                    */
    /* ================================================================ */

    t4_diag_trace(
        "before wall mesh local calculation",
        step);

    const int wall_mesh_status_local =
        calc_wall_mesh_diagnostics_tet4_tri3(
            wall_diag_surf,
            &(sys->fe),
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            U_ref,
            Uw,
            eU,
            gamma_n,
            write_nodal_yplus,
            &wall_mesh_diag);

    t4_diag_trace(
        "after wall mesh local calculation",
        step);

    wall_mesh_global_ok =
        t4_all_ranks_status_ok(
            wall_mesh_status_local,
            &(sys->mono_com),
            "wall mesh diagnostics");

    if (wall_mesh_global_ok) {
        t4_diag_trace(
            "before wall mesh diagnostics allreduce",
            step);

        wall_mesh_diagnostics_allreduce(
            &wall_mesh_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after wall mesh diagnostics allreduce",
            step);
    }

    /* ================================================================ */
    /* 2. Wall normal-BC / slip-speed / pressure diagnostics           */
    /* ================================================================ */

    const double p_ref_diag = 0.0;

    t4_diag_trace(
        "before wall velocity local calculation",
        step);

    const int wall_velocity_status_local =
        calc_wall_diagnostics_tet4_tri3_with_Uw(
            wall_diag_surf,
            &(sys->fe),
            wall_diag_basis,
            &(sys->vals),
            sys->vals.density,
            U_ref,
            p_ref_diag,
            eU,
            eP,
            Uw,
            &wall_velocity_diag);

    t4_diag_trace(
        "after wall velocity local calculation",
        step);

    wall_velocity_global_ok =
        t4_all_ranks_status_ok(
            wall_velocity_status_local,
            &(sys->mono_com),
            "wall normal/slip diagnostics");

    if (wall_velocity_global_ok) {
        t4_diag_trace(
            "before wall velocity diagnostics allreduce",
            step);

        wall_diagnostics_allreduce(
            &wall_velocity_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after wall velocity diagnostics allreduce",
            step);
    }

    /* ================================================================ */
    /* 3. Domain continuity diagnostics                                 */
    /* ================================================================ */

    t4_diag_trace(
        "before continuity local calculation",
        step);

    const int continuity_status_local =
        calc_domain_continuity_diagnostics_tet4(
            &(sys->fe),
            &(sys->vals),
            &continuity_diag);

    t4_diag_trace(
        "after continuity local calculation",
        step);

    continuity_global_ok =
        t4_all_ranks_status_ok(
            continuity_status_local,
            &(sys->mono_com),
            "domain continuity diagnostics");

    if (continuity_global_ok) {
        t4_diag_trace(
            "before continuity diagnostics allreduce",
            step);

        domain_continuity_diagnostics_allreduce(
            &continuity_diag,
            &(sys->mono_com));

        t4_diag_trace(
            "after continuity diagnostics allreduce",
            step);
    }

    if (wall_velocity_global_ok) {
        const int wall_location_status =
            output_wall_max_error_locations_tet4_tri3(
                wall_diag_surf,
                &(sys->fe),
                wall_diag_basis,
                &(sys->vals),
                U_ref,
                Uw,
                &(sys->mono_com),
                t,
                sys->cond.directory);

        if (
            wall_location_status != 0 &&
            rank == 0
        ) {
            /*
             * 最大誤差位置は補助診断であり、壁面統計本体とは独立。
             * 位置診断が失敗しても、すでにallreduce済みの
             * wall_velocity_diagは出力する。
             */
            fprintf(
                stderr,
                "[wall max location] failed; "
                "wall summary diagnostics will still be written.\n");
        }
    }

    if (continuity_global_ok) {
        const int div_location_status =
            output_domain_max_div_u_location_tet4(
                &(sys->fe),
                &(sys->vals),
                &(sys->mono_com),
                t,
                sys->cond.directory);

        if (
            div_location_status != 0 &&
            rank == 0
        ) {
            /*
             * 最大発散位置も補助診断。
             * 位置出力が失敗しても、領域発散統計本体は出力する。
             */
            fprintf(
                stderr,
                "[div-u max location] failed; "
                "continuity summary diagnostics will still be written.\n");
        }
    }

    /* ================================================================ */
    /* 4. Text output: rank 0 only                                      */
    /* ================================================================ */

    if (rank == 0) {
        if (wall_mesh_global_ok) {
            output_wall_mesh_diagnostics_tet4_tri3(
                &wall_mesh_diag,
                t,
                sys->cond.directory);
        }

        if (wall_velocity_global_ok) {
            fprintf(
                stderr,
                "[wall summary] writing at step=%d t=%.15e directory=%s\n",
                step,
                t,
                sys->cond.directory != NULL
                    ? sys->cond.directory
                    : "(null)");

            output_wall_diagnostics_tet4_tri3(
                &wall_velocity_diag,
                t,
                sys->cond.directory);
        } else {
            fprintf(
                stderr,
                "[wall summary] NOT written: "
                "wall_velocity_global_ok=0 at step=%d t=%.15e\n",
                step,
                t);
        }

        if (continuity_global_ok) {
            output_domain_continuity_diagnostics_tet4(
                &continuity_diag,
                t,
                sys->cond.directory);
        }
    }

    /*
     * どれか1つでも失敗した場合は、mainへ失敗を返す。
     * allreduceの順序自体は全rankで一致している。
     */
    if (
        !wall_mesh_global_ok ||
        !wall_velocity_global_ok ||
        !continuity_global_ok) {

        return -1;
    }

    return 0;
}

/*
 * Update the call in main(): add eP between eU and gamma_n.
 *
 *     run_wall_and_continuity_diagnostics_parallel_safe(
 *         &sys,
 *         step,
 *         t,
 *         U_ref,
 *         Uw,
 *         eU,
 *         eP,
 *         gamma_n,
 *         1);
 */



int main (
                int argc,
                char* argv[])
{
        printf("\n");

        FE_SYSTEM sys = {0};

        monolis_global_initialize();

    double t1 = monolis_get_time();
        double FOM_t1 = monolis_get_time();

        const T4WallOptions wall_options = t4_read_wall_options(argc, argv);

        int framework_argc = 0;
        char** framework_argv =
            t4_build_framework_argv(argc, argv, &framework_argc);

        sys.cond.directory =
            BBFE_fluid_get_directory_name(framework_argc, framework_argv, CODENAME);

        const double Uw[3] = {0.0, 0.0, 0.0};

        BBFE_fluid_set_tet4_wall_runtime_parameters(
            wall_options.mode,
            wall_options.exchange_mode,
            wall_options.gamma_n,
            wall_options.c_dt,
            wall_options.kappa,
            wall_options.B,
            wall_options.ut_eps,
            wall_options.feature_angle_deg);

        read_calc_conditions(&(sys.vals), sys.cond.directory);

        BBFE_fluid_pre(
                        &(sys.fe), &(sys.basis),
                        framework_argc, framework_argv, sys.cond.directory,
                        sys.vals.num_ip_each_axis);

        free(framework_argv);
        framework_argv = NULL;

        const char* filename;

        memory_allocation_nodal_values(
                        &(sys.vals),
                        sys.fe.total_num_nodes);

        filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
        BBFE_fluid_sups_read_Dirichlet_bc(
                        &(sys.bc),
                        filename,
                        sys.cond.directory,
                        sys.fe.total_num_nodes,
                        4);
/*
const int ground_bc_status =
    BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
        &(sys.bc),
        &(sys.fe),
        0.0,
        1.0e-10
    );

if (ground_bc_status != 0) {
    fprintf(
        stderr,
        "%s ERROR: failed to add ground no-slip condition.\n",
        CODENAME
    );

    exit(EXIT_FAILURE);
}
*/
    BBFE_fluid_sups_read_Dirichlet_bc_NR(
            &(sys.bc_NR),
            filename,
            sys.cond.directory,
            sys.fe.total_num_nodes,
            4);

            /*
    const int ground_bc_NR_status =
    BBFE_fluid_sups_add_zero_velocity_Dirichlet_on_z_plane(
        &(sys.bc_NR),
        &(sys.fe),
        0.0,
        1.0e-10
    );

    if (ground_bc_NR_status != 0) {
        fprintf(
            stderr,
            "%s ERROR: failed to add ground no-slip condition "
            "to Newton-increment BC.\n",
            CODENAME
        );

        exit(EXIT_FAILURE);
    }
    */

    pre_surface(
            //&(sys.monolis_com_surf),
            &(sys.surf),
            &(sys.basis_surf),
            sys.cond.directory,
            "surf_graph.dat",
            sys.vals.num_ip_each_axis);

    pre_surface_internal(
        //&(sys.monolis_com_surf),
        &(sys.surf_internal),
        &(sys.basis_surf_internal),
        sys.cond.directory,
        "surf_graph.dat",
        sys.vals.num_ip_each_axis);

    if (wall_options.mode == T4_WALL_STRONG_NOSLIP) {
        const double zero_increment[3] = {0.0, 0.0, 0.0};
        if (BBFE_fluid_sups_add_surface_velocity_Dirichlet(
                &(sys.bc), &(sys.surf), Uw) != 0 ||
            BBFE_fluid_sups_add_surface_velocity_Dirichlet(
                &(sys.bc_NR), &(sys.surf), zero_increment) != 0) {
            fprintf(stderr, "%s ERROR: failed to add strong Ahmed no-slip BC.\n", CODENAME);
            exit(EXIT_FAILURE);
        }
    }

    FILE* fp;

   if(monolis_mpi_get_global_my_rank() == 0){
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_velocity.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_velocity.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_Nitsche.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_drag_coeff_pressure.txt", sys.cond.directory);
            fclose(fp);
            fp = BBFE_sys_write_fopen(fp, "cylinder_lift_coeff_pressure.txt", sys.cond.directory);
            fclose(fp);

            const char* wall_mode_diag_files[] = {
                "cylinder_drag_coeff_walllaw.txt",
                "cylinder_drag_coeff_pressure_plus_walllaw.txt",
                "cylinder_drag_coeff_modelled_total.txt",
                "cylinder_drag_coeff_normal_molecular_viscous.txt",
                "cylinder_drag_coeff_molecular_viscous_splitcheck.txt",
                "cylinder_drag_coeff_normal_penalty_reaction.txt",
                "ahmed_support_z_cut.txt",

                "ahmed_mainbody_area.txt",
                "ahmed_mainbody_face_count.txt",
                "ahmed_mainbody_Aproj.txt",
                "ahmed_mainbody_Aproj_over_Aref.txt",
                "ahmed_mainbody_int_nx.txt",
                "ahmed_mainbody_int_ny.txt",
                "ahmed_mainbody_int_nz.txt",
                "ahmed_mainbody_Cd_pressure.txt",
                "ahmed_mainbody_Cd_molecular_viscous_diagnostic.txt",
                "ahmed_mainbody_Cd_normal_molecular_viscous.txt",
                "ahmed_mainbody_Cd_walllaw.txt",
                "ahmed_mainbody_Cd_normal_penalty_reaction.txt",
                "ahmed_mainbody_Cd_pressure_plus_walllaw.txt",
                "ahmed_mainbody_Cd_modelled_total.txt",

                "ahmed_support_area.txt",
                "ahmed_support_face_count.txt",
                "ahmed_support_Aproj.txt",
                "ahmed_support_Aproj_over_Aref.txt",
                "ahmed_support_int_nx.txt",
                "ahmed_support_int_ny.txt",
                "ahmed_support_int_nz.txt",
                "ahmed_support_Cd_pressure.txt",
                "ahmed_support_Cd_molecular_viscous_diagnostic.txt",
                "ahmed_support_Cd_normal_molecular_viscous.txt",
                "ahmed_support_Cd_walllaw.txt",
                "ahmed_support_Cd_normal_penalty_reaction.txt",
                "ahmed_support_Cd_pressure_plus_walllaw.txt",
                "ahmed_support_Cd_modelled_total.txt"
            };

            const int n_wall_mode_diag_files =
                (int)(sizeof(wall_mode_diag_files) / sizeof(wall_mode_diag_files[0]));

            for (int i = 0; i < n_wall_mode_diag_files; ++i) {
                fp = BBFE_sys_write_fopen(
                    fp,
                    wall_mode_diag_files[i],
                    sys.cond.directory);
                fclose(fp);
            }

            //fclose(fp);
    }

    BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

        BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

        //intialize for velocity and pressure
    initialize_velocity_pressure_karman_vortex(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);
    //initialize_velocity_pressure_ahmedbody(&(sys.surf), sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);
    double* vec = BB_std_calloc_1d_double(vec, sys.fe.total_num_nodes*4);

    BBFE_fluid_sups_add_velocity_pressure(
            sys.vals.v,
            sys.vals.p,
            vec,
            sys.fe.total_num_nodes);

    ROM_sys_hlpod_fe_add_Dbc(
            vec,
            &(sys.bc),
            sys.fe.total_num_nodes,
            4);

    BBFE_fluid_sups_renew_velocity(
            sys.vals.v,
            vec,
            sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
            sys.vals.p,
            vec,
            sys.fe.total_num_nodes);

    ROM_BB_vec_copy_2d(
            sys.vals.v,
            sys.vals.v_old,
            sys.fe.total_num_nodes,
            3);


output_files(&sys, 0, 0);
        /****************** solver ********************/
        double t = 0.0;
        int step = 0;
        int count = 0;  //for ROM

        t = 0.0; step = 0;

        /*
        //if(sys.rom_prm_p.hot_start == 1){
            char fname[BUFFER_SIZE];
            snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
            double* val = BB_std_calloc_1d_double(val, 4*sys.fe.total_num_nodes);
            double t_hs = hot_start_read_initialize_val(val, fname, sys.cond.directory);
            int step_hs = 0;

            printf("Hot start time: %lf\n", t);
            //printf("Hot start step: %d\n", step_rom);
            printf("sys.vals.finish_time - t = %lf\n", ((double)sys.vals.finish_time - t));

            BBFE_fluid_sups_renew_velocity(sys.vals.v, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_velocity(sys.vals.v_old, val, sys.fe.total_num_nodes);
            BBFE_fluid_sups_renew_pressure(sys.vals.p, val, sys.fe.total_num_nodes);

            //BBFE_fluid_sups_renew_velocity(sys.vals_rom.v, val, sys.fe.total_num_nodes);
            //BBFE_fluid_sups_renew_pressure(sys.vals_rom.p, val, sys.fe.total_num_nodes);

            BB_std_free_1d_double(val, 4*sys.fe.total_num_nodes);
        //}
        */

if (wall_options.mode == T4_WALL_STRONG_NOSLIP) {
    BBFE_fluid_sups_impose_surface_velocity_on_values(
        &(sys.surf), sys.vals.v, Uw);
    BBFE_fluid_sups_impose_surface_velocity_on_values(
        &(sys.surf), sys.vals.v_old, Uw);
}

//const double U_ref = 60.0;
//const double A_ref = 0.112032;

const double U_ref = 1.0;
const double A_ref = 0.08;

/*
 * Geometry-specific split plane: nominal Ahmed underbody / support-leg top.
 * Faces with centroid z < 0.05 m are classified as support legs.
 */
const double support_z_cut =
    read_t4_ahmed_support_z_cut_from_env();

const double eU[3] = {
    1.0,
    0.0,
    0.0
};

const double eP[3] = {
    0.0,
    1.0,
    0.0
};

/* Solver and force diagnostics use the same run-time wall configuration. */
const double gamma_n = wall_options.gamma_n;

if (monolis_mpi_get_global_my_rank() == 0) {
    printf("%s Ahmed wall mode: %s\n", CODENAME, t4_wall_mode_name(wall_options.mode));
    printf("%s   exchange=%s gamma=%.6g c_dt=%.6g kappa=%.6g B=%.6g feature_angle=%.3f deg\n",
        CODENAME,
        wall_options.exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE ? "opposite-node" : "wall-q",
        wall_options.gamma_n,
        wall_options.c_dt,
        wall_options.kappa,
        wall_options.B,
        wall_options.feature_angle_deg);

    if (wall_options.mode == T4_WALL_STRONG_NORMAL_FREESLIP ||
        wall_options.mode == T4_WALL_STRONG_NORMAL_SPALDING) {
        printf("%s   strong-normal implementation: feature-aware nodal state/update projection\n", CODENAME);
    }

    printf(
        "%s Ahmed force split: face-centroid z < %.6f m -> supports; otherwise main body\n",
        CODENAME,
        support_z_cut);

    FILE* wall_cfg_fp = NULL;
    wall_cfg_fp = BBFE_sys_write_fopen(
        wall_cfg_fp,
        "wall_runtime_config.txt",
        sys.cond.directory);
    fprintf(wall_cfg_fp, "wall_mode %s\n", t4_wall_mode_name(wall_options.mode));
    fprintf(wall_cfg_fp, "wall_mode_id %d\n", wall_options.mode);
    fprintf(wall_cfg_fp, "exchange %s\n",
        wall_options.exchange_mode == T4_WALL_EXCHANGE_OPPOSITE_NODE ? "opposite-node" : "wall-q");
    fprintf(wall_cfg_fp, "gamma_n %.17e\n", wall_options.gamma_n);
    fprintf(wall_cfg_fp, "c_dt %.17e\n", wall_options.c_dt);
    fprintf(wall_cfg_fp, "kappa %.17e\n", wall_options.kappa);
    fprintf(wall_cfg_fp, "B %.17e\n", wall_options.B);
    fprintf(wall_cfg_fp, "ut_eps %.17e\n", wall_options.ut_eps);
    fprintf(wall_cfg_fp, "feature_angle_deg %.17e\n", wall_options.feature_angle_deg);
    fprintf(wall_cfg_fp, "ground_wall strong-noslip\n");
    fclose(wall_cfg_fp);
}

/*
 * デバッグ中は100ステップごと程度にする。
 * 毎ステップowner mapを再構築すると非常に重い。
 */
const int diagnostic_interval = 1;
sys.vals.C_vms = 0.0;
sys.vals.vms_cap_coeff = 0.0;
    while (t < sys.vals.finish_time) {
        t += sys.vals.dt;
        step += 1;

        if (monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "\n%s ----------------- step %d ----------------\n",
                CODENAME,
                step);

            fflush(stdout);
        }

        /* ================================================================ */
        /* Solver                                                           */
        /* ================================================================ */

        t4_diag_trace(
            "before solver_fom_VMS",
            step);

        solver_fom_VMS(
            sys,
            &(sys.surf),
            &(sys.basis_surf),
            t,
            count);
printf(
    "[main after solver] C_vms=%.17g cap=%.17g\n",
    sys.vals.C_vms,
    sys.vals.vms_cap_coeff);

t4_diag_trace(
            "after solver_fom_VMS",
            step);

        count++;

        /* ================================================================ */
        /* Heavy diagnostics                                                */
        /* ================================================================ */

        const int do_diagnostics =
            (
                step == 1 ||
                step % diagnostic_interval == 0
            );

        if (do_diagnostics) {
            t4_diag_trace(
                "before parallel diagnostics wrapper",
                step);

            const int diagnostic_status =
                run_wall_and_continuity_diagnostics_parallel_safe(
                    &sys,
                    step,
                    t,
                    U_ref,
                    Uw,
                    eU,
                    eP,
                    gamma_n,
                    1);

            t4_diag_trace(
                "after parallel diagnostics wrapper",
                step);

            if (
                diagnostic_status != 0 &&
                monolis_mpi_get_global_my_rank() == 0) {

                fprintf(
                    stderr,
                    "[main] diagnostics failed at step %d\n",
                    step);
            }
        }

        /* ================================================================ */
        /* VTK output                                                       */
        /* ================================================================ */

        /*
        * output_filesはpartitionごとのVTK出力を行うので、
        * 原則として全rankが呼ぶ。
        *
        * yPlusを更新してから呼ぶ。
        */
        if (step % sys.vals.output_interval == 0) {
            t4_diag_trace(
                "before output_files",
                step);

            output_files(
                &sys,
                step,
                t);

            t4_diag_trace(
                "after output_files",
                step);

        }

        /* ================================================================ */
        /* Existing force / gamma / wall diagnostics                        */
        /* ================================================================ */

        /*
        * 既存の重い診断も毎ステップではなく、
        * diagnostic_intervalごとに変更する。
        */
        if (do_diagnostics) {
            double D_out = 0.0;
            double L_out = 0.0;
            double Cd_out = 0.0;
            double Cl_out = 0.0;

            double D_nitsche_out = 0.0;
            double L_nitsche_out = 0.0;
            double Cd_nitsche_out = 0.0;
            double Cl_nitsche_out = 0.0;

            double Cd_out_p = 0.0;
            double Cl_out_p = 0.0;

            const double qA =
                0.5
                * sys.vals.density
                * U_ref
                * U_ref
                * A_ref;

            /* ------------------------------------------------------------ */
            /* Local force calculation                                      */
            /* ------------------------------------------------------------ */

            const int force_status_local =
                calc_Cd_v_tet4_tri3_nitsche(
                    &(sys.surf_internal),
                    &(sys.fe),
                    &(sys.basis_surf_internal),
                    &(sys.mono_com),
                    &(sys.vals),
                    sys.vals.density,
                    sys.vals.viscosity,
                    U_ref,
                    A_ref,
                    eU,
                    eP,
                    Uw,
                    gamma_n,
                    &D_out,
                    &L_out,
                    &Cd_out,
                    &Cl_out,
                    &D_nitsche_out,
                    &L_nitsche_out,
                    &Cd_nitsche_out,
                    &Cl_nitsche_out);

            /*
            * 全rankでstatusを統一する。
            */
            const int force_global_ok =
                t4_all_ranks_status_ok(
                    force_status_local,
                    &(sys.mono_com),
                    "viscous/Nitsche force diagnostics");

            if (force_global_ok) {
                /*
                * 必ず全rankが同じ順番で呼ぶ。
                */
                monolis_allreduce_R(
                    1,
                    &Cd_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cl_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cd_nitsche_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);

                monolis_allreduce_R(
                    1,
                    &Cl_nitsche_out,
                    MONOLIS_MPI_SUM,
                    sys.mono_com.comm);
            }

            /* ------------------------------------------------------------ */
            /* Pressure force                                               */
            /* ------------------------------------------------------------ */

            calc_Cd_p_hex_or_tet(
                &(sys.surf_internal),
                &(sys.fe),
                &(sys.basis_surf_internal),
                &(sys.vals),
                eU,
                eP,
                &Cd_out_p,
                &Cl_out_p);

            monolis_allreduce_R(
                1,
                &Cd_out_p,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);

            monolis_allreduce_R(
                1,
                &Cl_out_p,
                MONOLIS_MPI_SUM,
                sys.mono_com.comm);

            /* ================================================================ */
            /* Case C: main-body / support-leg force split                      */
            /* ================================================================ */

            T4AhmedBodySupportDiagnostics split_diag;
            memset(&split_diag, 0, sizeof(split_diag));

            const int split_status_local =
                calc_ahmed_body_support_diagnostics_tet4_tri3_local(
                    &(sys.surf_internal),
                    &(sys.fe),
                    &(sys.basis_surf_internal),
                    &(sys.vals),
                    sys.vals.density,
                    sys.vals.viscosity,
                    Uw,
                    gamma_n,
                    support_z_cut,
                    eU,
                    &split_diag);

            const int split_global_ok =
                t4_all_ranks_status_ok(
                    split_status_local,
                    &(sys.mono_com),
                    "Ahmed body/support force split");

            if (split_global_ok) {
                ahmed_body_support_diagnostics_allreduce(
                    &split_diag,
                    &(sys.mono_com));
            }

            /* ================================================================ */
            /* Nitsche residual-row reaction diagnostic                         */
            /* ここへ t4n_nitsche_row_sum_main_patch.c の内容を追加              */
            /* ================================================================ */

            T4NitscheRowSumDiagnostics row_diag;

            const int row_status_local =
                calc_nitsche_row_sum_reaction_tet4_tri3_local(
                    &(sys.surf_internal),
                    &(sys.fe),
                    &(sys.basis_surf_internal),
                    &(sys.vals),
                    sys.vals.density,
                    sys.vals.viscosity,
                    Uw,
                    gamma_n,
                    &row_diag);

            const int row_global_ok =
                t4_all_ranks_status_ok(
                    row_status_local,
                    &(sys.mono_com),
                    "Nitsche residual-row reaction");

            if (row_global_ok) {
                nitsche_row_sum_diagnostics_allreduce(
                    &row_diag,
                    &(sys.mono_com));
            }

            /* ------------------------------------------------------------ */
            /* Shared text output: rank 0 only                              */
            /* ------------------------------------------------------------ */

            if (monolis_mpi_get_global_my_rank() == 0) {
                if (force_global_ok) {
                    ROM_std_hlpod_output_add_calc_time(
                        Cd_out / qA,
                        t,
                        "cylinder_drag_coeff_velocity.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cl_out / qA,
                        t,
                        "cylinder_lift_coeff_velocity.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cd_nitsche_out / qA,
                        t,
                        "cylinder_drag_coeff_Nitsche.txt",
                        sys.cond.directory);

                    ROM_std_hlpod_output_add_calc_time(
                        Cl_nitsche_out / qA,
                        t,
                        "cylinder_lift_coeff_Nitsche.txt",
                        sys.cond.directory);
                }

                ROM_std_hlpod_output_add_calc_time(
                    Cd_out_p / qA,
                    t,
                    "cylinder_drag_coeff_pressure.txt",
                    sys.cond.directory);

                ROM_std_hlpod_output_add_calc_time(
                    Cl_out_p / qA,
                    t,
                    "cylinder_lift_coeff_pressure.txt",
                    sys.cond.directory);

                if (split_global_ok) {
                    output_ahmed_body_support_drag_diagnostics(
                        &split_diag,
                        t,
                        qA,
                        A_ref,
                        eU,
                        sys.cond.directory);
                }
            }
        }

        /* ================================================================ */
        /* Hot-start output, when enabled                                   */
        /* ================================================================ */

        BBFE_fluid_sups_add_velocity_pressure(
                        sys.vals.v,
                        sys.vals.p,
                        sys.monolis.mat.R.X,
                        sys.fe.total_num_nodes);

        if(step%1 == 0){
                char fname[BUFFER_SIZE];
                snprintf(fname, BUFFER_SIZE, "hot_start/%s.%lf.%d.dat", "velosity_pressure", sys.vals.density, monolis_mpi_get_global_my_rank());
                //hot_start_write_initialize_val(sys.monolis.mat.R.X, sys.fe.total_num_nodes, 4, t, fname, sys.cond.directory);
            }
        else{
        }
    }

        BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
        BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), sys.fe.total_num_nodes, 4);
        monolis_finalize(&(sys.monolis));

        double t2 = monolis_get_time();
        int myrank = monolis_mpi_get_global_my_rank();

        if(myrank == 0) {
                printf("** Total time: %f\n", t2 - t1);
        }

        monolis_global_finalize();

        printf("\n");

        return 0;
}
