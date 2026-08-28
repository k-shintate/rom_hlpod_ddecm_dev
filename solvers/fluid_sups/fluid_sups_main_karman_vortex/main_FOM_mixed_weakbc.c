#include "core_FOM_mixed.h"
#include "core_Nitsche_prism.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Existing TET4 force diagnostic uses these runtime Nitsche parameters.
 * The implementation in core_FOM.c has C linkage, so the declaration here
 * must use the same linkage when this .c file is compiled as C++. */
#ifdef __cplusplus
extern "C" {
#endif
void BBFE_fluid_set_nitsche_runtime_parameters(
    double gamma_n, double dt_penalty_coeff);
#ifdef __cplusplus
}
#endif

/* -------------------------------------------------------------------------- */
/* Existing TET4/TRI3 physical-force routines from core_NR.c.                 */
/* These return LOCAL force contributions; main performs MPI_SUM.              */
/* -------------------------------------------------------------------------- */
void calc_Cd_p_hex_or_tet(
    BBFE_DATA* surf,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    VALUES* vals,
    const double eU_in[3],
    const double eP_in[3],
    double* D_out,
    double* L_out);

int calc_Cd_v_hex_or_tet_nitsche(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis,
    MONOLIS_COM* monolis_com,
    const VALUES* vals,
    double rho,
    double mu,
    double Uinf,
    double Aref,
    const double eU_in[3],
    const double eP_in[3],
    const double Uw[3],
    double gamma_n,
    double* D_out,
    double* L_out,
    double* Cd_out,
    double* Cl_out,
    double* D_nitsche_out,
    double* L_nitsche_out,
    double* Cd_nitsche_out,
    double* Cl_nitsche_out);

static const char* INPUT_FILENAME_D_BC_V = "D_bc_v.dat";

/*
 * Assembly surfaces: MPI overlap/duplicate copies are allowed.
 * These are used ONLY for distributed Nitsche matrix/RHS assembly.
 */
static const char* INPUT_FILENAME_SURF_PRI = "surf_prism_graph.dat";
static const char* INPUT_FILENAME_SURF_TET = "surf_tet_graph.dat";
static const char* INPUT_FILENAME_SURF_SINGLE = "surf_graph.dat";

/*
 * The SAME partitioned surface file is read in two ways:
 *   assembly : all local + overlap copies
 *   internal : only the first N_internal entries described by
 *              <filename>.n_internal.<rank>
 */

/*
 * No mesh-specific Ahmed wall-face count is hard-coded here.
 * A valid mesh is accepted when every supplied TRI3 wall face can be mapped
 * to the requested owner element family. Physical/internal surfaces are
 * taken from the partitioner's *.n_internal.<rank> metadata.
 */
static const double AHMED_INITIAL_U = 60.0;
static const double KARMAN_INITIAL_U_DEFAULT = 1.0;

typedef enum {
    FLUID_CASE_AHMED = 0,
    FLUID_CASE_KARMAN = 1
} FluidCaseMode;

static FluidCaseMode fluid_case_mode_from_env(void)
{
    const char* text = getenv("FLUID_CASE");
    if(text == NULL || text[0] == '\0') return FLUID_CASE_AHMED;

    if(strcmp(text, "ahmed") == 0 || strcmp(text, "AHMED") == 0 ||
       strcmp(text, "ahmedbody") == 0 || strcmp(text, "AHMEDBODY") == 0) {
        return FLUID_CASE_AHMED;
    }

    if(strcmp(text, "karman") == 0 || strcmp(text, "KARMAN") == 0 ||
       strcmp(text, "karman-vortex") == 0 || strcmp(text, "KARMAN-VORTEX") == 0) {
        return FLUID_CASE_KARMAN;
    }

    fprintf(stderr,
        "%s ERROR: unsupported FLUID_CASE=\"%s\". Use ahmed or karman.\n",
        CODENAME, text);
    exit(EXIT_FAILURE);
}

static const char* fluid_case_mode_name(FluidCaseMode mode)
{
    return (mode == FLUID_CASE_KARMAN) ? "karman" : "ahmed";
}

typedef int (*WallOwnerMapFn)(
    const BBFE_DATA*,
    const BBFE_DATA*,
    int**,
    int**,
    int*);

static int initialize_velocity_pressure_uniform(
    double** v,
    double* p,
    int nnode,
    double initial_u)
{
    if(v == NULL || p == NULL || nnode <= 0) return -1;

    for(int i=0; i<nnode; ++i) {
        v[i][0] = initial_u;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
        p[i] = 0.0;
    }
    return 0;
}

static int initialize_velocity_pressure_ahmed(
    double** v,
    double* p,
    int nnode,
    double inlet_u,
    const BBFE_DATA* surf_pri,
    const BBFE_DATA* surf_tet,
    int zero_pri_wall,
    int zero_tet_wall)
{
    if(v == NULL || p == NULL || nnode <= 0) {
        fprintf(stderr,
            "%s ERROR: invalid Ahmed initial-condition storage.\n",
            CODENAME);
        return -1;
    }

    /*
     * Far-field initial condition.  Strong Dirichlet rows are applied later.
     */
    for(int i=0; i<nnode; ++i) {
        v[i][0] = inlet_u;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
        p[i] = 0.0;
    }

    if(!zero_pri_wall && !zero_tet_wall) {
        return 0;
    }

    /*
     * Ahmed wall is weak/Nitsche, so D_bc_v.dat does not overwrite its
     * initial velocity.  Zero only the wall families whose Nitsche treatment
     * is enabled.  This is important for split A/B tests: an OFF family must
     * not be initialized to zero and then left unconstrained.
     *
     * Assembly surfaces may contain MPI overlap/duplicate faces.  That is
     * harmless here because this operation is idempotent assignment, not a
     * physical SUM reduction.
     */
    unsigned char* wall_node =
        (unsigned char*)calloc((size_t)nnode, sizeof(unsigned char));

    if(wall_node == NULL) {
        fprintf(stderr,
            "%s ERROR: failed to allocate Ahmed-wall initial-condition mask.\n",
            CODENAME);
        return -1;
    }

    const BBFE_DATA* surfaces[2] = {surf_pri, surf_tet};
    const int enabled[2] = {zero_pri_wall, zero_tet_wall};

    for(int is=0; is<2; ++is) {
        if(!enabled[is]) continue;
        const BBFE_DATA* surf = surfaces[is];
        if(surf == NULL || surf->total_num_elems <= 0) continue;

        for(int e=0; e<surf->total_num_elems; ++e) {
            for(int a=0; a<surf->local_num_nodes; ++a) {
                const int nid = surf->conn[e][a];

                if(nid < 0 || nid >= nnode) {
                    fprintf(stderr,
                        "%s ERROR: Ahmed initial-condition surface node "
                        "out of range: family=%s elem=%d local_node=%d "
                        "nid=%d nnode=%d.\n",
                        CODENAME,
                        (is == 0 ? "PRI6" : "TET4"),
                        e, a, nid, nnode);
                    free(wall_node);
                    return -1;
                }

                wall_node[nid] = 1;
            }
        }
    }

    int num_wall_nodes = 0;
    for(int i=0; i<nnode; ++i) {
        if(!wall_node[i]) continue;

        v[i][0] = 0.0;
        v[i][1] = 0.0;
        v[i][2] = 0.0;
        ++num_wall_nodes;
    }

    free(wall_node);
    return num_wall_nodes;
}

/*
 * Remove the constant-pressure null space when no pressure Dirichlet row
 * already exists.  Only owned/internal pressure rows are counted so MPI
 * overlap copies do not make the test ambiguous.
 *
 * Policy:
 *   - If an owned pressure Dirichlet row already exists anywhere, keep it.
 *   - Otherwise rank 0 pins the pressure of its first internal node to zero.
 *   - The absolute/state BC and Newton-increment BC are pinned together.
 */

static const char* fluid_env_string_or_default(
    const char* name,
    const char* default_value)
{
    const char* text = getenv(name);
    return (text != NULL && text[0] != '\0') ? text : default_value;
}

static int ahmed_env_flag(const char* name, int default_value)
{
    const char* text = getenv(name);
    if(text == NULL || text[0] == '\0') return default_value;
    if(strcmp(text, "0") == 0 || strcmp(text, "false") == 0 ||
       strcmp(text, "FALSE") == 0 || strcmp(text, "no") == 0 ||
       strcmp(text, "NO") == 0) return 0;
    return 1;
}

static double ahmed_env_positive_double(
    const char* name,
    double default_value,
    int allow_zero)
{
    const char* text = getenv(name);
    if(text == NULL || text[0] == '\0') return default_value;

    char* endptr = NULL;
    const double value = strtod(text, &endptr);
    const int bad =
        (endptr == text || endptr == NULL || *endptr != '\0' ||
         !isfinite(value) ||
         (allow_zero ? value < 0.0 : value <= 0.0));

    if(bad) {
        fprintf(stderr,
            "%s ERROR: invalid %s=\"%s\".\n",
            CODENAME, name, text);
        exit(EXIT_FAILURE);
    }
    return value;
}

static int ahmed_wall_mode_setting(
    const char* family_name,
    const char* common_name,
    int legacy_enable,
    int* enabled_out)
{
    const char* text = getenv(family_name);
    if(text == NULL || text[0] == '\0') text = getenv(common_name);

    if(text == NULL || text[0] == '\0') {
        if(enabled_out != NULL) *enabled_out = legacy_enable ? 1 : 0;
        return PRI6_WALL_NITSCHE_NOSLIP;
    }

    if(strcmp(text,"off") == 0 || strcmp(text,"none") == 0 || strcmp(text,"disabled") == 0) {
        if(enabled_out != NULL) *enabled_out = 0;
        return PRI6_WALL_NITSCHE_NOSLIP;
    }

    int mode = -1;
    if(BBFE_pri6_wall_mode_from_string(text,&mode) != 0) {
        fprintf(stderr,
            "%s ERROR: invalid %s=\"%s\". Allowed: off, strong-noslip, "
            "nitsche-noslip, nitsche-freeslip, nitsche-spalding, "
            "strong-normal-freeslip, strong-normal-spalding.\n",
            CODENAME,family_name,text);
        exit(EXIT_FAILURE);
    }
    if(enabled_out != NULL) *enabled_out = 1;
    return mode;
}

static int ahmed_wall_exchange_setting(
    const char* family_name,
    const char* common_name,
    int default_value)
{
    const char* text = getenv(family_name);
    if(text == NULL || text[0] == '\0') text = getenv(common_name);
    if(text == NULL || text[0] == '\0') return default_value;

    int mode = -1;
    if(BBFE_pri6_wall_exchange_from_string(text,&mode) != 0) {
        fprintf(stderr,
            "%s ERROR: invalid %s=\"%s\". Allowed: wall-q or opposite-face.\n",
            CODENAME,family_name,text);
        exit(EXIT_FAILURE);
    }
    return mode;
}

static double ahmed_env_family_double(
    const char* family_name,
    const char* common_name,
    double default_value,
    int allow_zero)
{
    const char* family = getenv(family_name);
    if(family != NULL && family[0] != '\0') {
        return ahmed_env_positive_double(family_name,default_value,allow_zero);
    }
    return ahmed_env_positive_double(common_name,default_value,allow_zero);
}

static void ahmed_add_strong_noslip_surface_bc(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf,
    const char* label)
{
    if(sys == NULL || surf == NULL) return;
    const int nnode = sys->fe_tet.total_num_nodes;
    int local_new_rows = 0;
    int local_nodes = 0;
    unsigned char* seen = (unsigned char*)calloc((size_t)(nnode > 0 ? nnode : 1),1);
    if(seen == NULL) {
        fprintf(stderr,"%s ERROR: strong wall BC allocation failed.\n",CODENAME);
        exit(EXIT_FAILURE);
    }

    for(int e=0; e<surf->total_num_elems; ++e) {
        for(int a=0; a<surf->local_num_nodes; ++a) {
            const int gid = surf->conn[e][a];
            if(gid < 0 || gid >= nnode) {
                free(seen);
                fprintf(stderr,"%s ERROR: invalid strong wall node %d.\n",CODENAME,gid);
                exit(EXIT_FAILURE);
            }
            if(!seen[gid]) { seen[gid]=1; ++local_nodes; }
            for(int d=0; d<3; ++d) {
                const int dof = 4*gid+d;
                if(sys->bc.D_bc_exists[dof] && fabs(sys->bc.imposed_D_val[dof]) > 1.0e-12) {
                    free(seen);
                    fprintf(stderr,
                        "%s ERROR: %s strong-noslip conflicts with nonzero strong BC "
                        "at node=%d component=%d value=%e.\n",
                        CODENAME,label ? label : "Ahmed",gid,d,sys->bc.imposed_D_val[dof]);
                    exit(EXIT_FAILURE);
                }
                if(!sys->bc.D_bc_exists[dof]) { ++sys->bc.num_D_bcs; ++local_new_rows; }
                sys->bc.D_bc_exists[dof] = 1;
                sys->bc.imposed_D_val[dof] = 0.0;
                if(!sys->bc_NR.D_bc_exists[dof]) ++sys->bc_NR.num_D_bcs;
                sys->bc_NR.D_bc_exists[dof] = 1;
                sys->bc_NR.imposed_D_val[dof] = 0.0;
            }
        }
    }
    free(seen);

    int global_nodes = local_nodes;
    int global_rows = local_new_rows;
    monolis_allreduce_I(1,&global_nodes,MONOLIS_MPI_SUM,sys->mono_com.comm);
    monolis_allreduce_I(1,&global_rows,MONOLIS_MPI_SUM,sys->mono_com.comm);
    if(monolis_mpi_get_global_my_rank() == 0) {
        printf("%s %s runtime strong-noslip: node_copies=%d new_velocity_rows=%d.\n",
               CODENAME,label ? label : "Ahmed",global_nodes,global_rows);
    }
}

static int ahmed_env_positive_int(
    const char* name,
    int default_value)
{
    const char* text = getenv(name);
    if(text == NULL || text[0] == '\0') return default_value;

    char* endptr = NULL;
    const long value = strtol(text, &endptr, 10);
    if(endptr == text || endptr == NULL || *endptr != '\0' || value <= 0 || value > 100000000L) {
        fprintf(stderr,
            "%s ERROR: invalid %s=\"%s\".\n",
            CODENAME, name, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static void scale_absolute_velocity_dirichlet_bc(
    BBFE_BC* bc,
    int nnode,
    double scale,
    const char* label)
{
    if(bc == NULL || bc->D_bc_exists == NULL ||
       bc->imposed_D_val == NULL || nnode <= 0) {
        fprintf(stderr,
            "%s ERROR: invalid velocity Dirichlet storage for %s scaling.\n",
            CODENAME, label ? label : "unknown");
        exit(EXIT_FAILURE);
    }

    if(!isfinite(scale) || scale <= 0.0) {
        fprintf(stderr,
            "%s ERROR: invalid velocity scale %.15e for %s.\n",
            CODENAME, scale, label ? label : "unknown");
        exit(EXIT_FAILURE);
    }

    int local_rows = 0;
    int local_nonzero = 0;
    double local_before_max = 0.0;
    double local_after_max = 0.0;

    for(int i=0; i<nnode; ++i) {
        for(int d=0; d<3; ++d) {
            const int dof = 4*i + d;
            if(!bc->D_bc_exists[dof]) continue;

            ++local_rows;
            const double before = bc->imposed_D_val[dof];
            if(fabs(before) > 0.0) ++local_nonzero;
            if(fabs(before) > local_before_max) local_before_max = fabs(before);

            bc->imposed_D_val[dof] = scale * before;
            if(fabs(bc->imposed_D_val[dof]) > local_after_max) {
                local_after_max = fabs(bc->imposed_D_val[dof]);
            }
        }
    }

    int global_rows = local_rows;
    int global_nonzero = local_nonzero;
    double global_before_max = local_before_max;
    double global_after_max = local_after_max;

    monolis_allreduce_I(1, &global_rows, MONOLIS_MPI_SUM,
                        monolis_mpi_get_global_comm());
    monolis_allreduce_I(1, &global_nonzero, MONOLIS_MPI_SUM,
                        monolis_mpi_get_global_comm());
    monolis_allreduce_R(1, &global_before_max, MONOLIS_MPI_MAX,
                        monolis_mpi_get_global_comm());
    monolis_allreduce_R(1, &global_after_max, MONOLIS_MPI_MAX,
                        monolis_mpi_get_global_comm());

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s %s velocity BC scale=%.9g: rows=%d nonzero_rows=%d "
            "max|u_D| %.9g -> %.9g.\n",
            CODENAME,
            label ? label : "absolute",
            scale,
            global_rows,
            global_nonzero,
            global_before_max,
            global_after_max);
    }
}

static void ensure_pressure_gauge(FE_SYSTEM_FLUID_MIXED* sys)
{
    if(sys == NULL) {
        fprintf(stderr, "%s ERROR: pressure-gauge system is NULL.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    const int rank = monolis_mpi_get_global_my_rank();
    const int nint = sys->mono_com.n_internal_vertex;

    int local_pressure_bc = 0;

    if(sys->bc_NR.D_bc_exists != NULL) {
        for(int i = 0; i < nint; ++i) {
            if(sys->bc_NR.D_bc_exists[4*i + 3]) {
                ++local_pressure_bc;
            }
        }
    }

    int global_pressure_bc = local_pressure_bc;
    monolis_allreduce_I(
        1,
        &global_pressure_bc,
        MONOLIS_MPI_SUM,
        sys->mono_com.comm);

    if(global_pressure_bc > 0) {
        if(rank == 0) {
            printf(
                "%s pressure gauge: existing owned pressure Dirichlet rows=%d; "
                "no extra gauge added.\n",
                CODENAME,
                global_pressure_bc);
        }
        return;
    }

    if(rank == 0) {
        if(nint <= 0) {
            fprintf(stderr,
                "%s ERROR: rank 0 has no internal node for pressure gauge.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }

        if(sys->bc.D_bc_exists == NULL ||
           sys->bc.imposed_D_val == NULL ||
           sys->bc_NR.D_bc_exists == NULL ||
           sys->bc_NR.imposed_D_val == NULL) {
            fprintf(stderr,
                "%s ERROR: Dirichlet arrays are not allocated for pressure gauge.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }

        const int node = 0;
        const int dof = 4*node + 3;

        if(!sys->bc.D_bc_exists[dof]) {
            sys->bc.num_D_bcs += 1;
        }
        sys->bc.D_bc_exists[dof] = 1;
        sys->bc.imposed_D_val[dof] = 0.0;

        if(!sys->bc_NR.D_bc_exists[dof]) {
            sys->bc_NR.num_D_bcs += 1;
        }
        sys->bc_NR.D_bc_exists[dof] = 1;
        sys->bc_NR.imposed_D_val[dof] = 0.0;

        printf(
            "%s pressure gauge: rank=0 local_internal_node=%d p=0 delta_p=0.\n",
            CODENAME,
            node);
    }
}

static void precheck_wall_owner_map(
    const char* label,
    const BBFE_DATA* surf,
    const BBFE_DATA* volume,
    WallOwnerMapFn make_owner_map,
    int is_physical_internal_surface,
    MONOLIS_COM* mono_com)
{
    int* owner = NULL;
    int* lface = NULL;
    int local_unmapped = 0;

    if(make_owner_map(
        surf,
        volume,
        &owner,
        &lface,
        &local_unmapped) != 0) {

        fprintf(stderr,
            "%s ERROR: %s owner-map precheck failed on rank %d.\n",
            CODENAME,
            label,
            monolis_mpi_get_global_my_rank());
        free(owner);
        free(lface);
        exit(EXIT_FAILURE);
    }

    const int local_surface = surf->total_num_elems;
    const int local_mapped = local_surface - local_unmapped;

    printf(
        "%s rank=%d %s owner-map: surf=%d volume=%d mapped=%d unmapped=%d\n",
        CODENAME,
        monolis_mpi_get_global_my_rank(),
        label,
        local_surface,
        volume->total_num_elems,
        local_mapped,
        local_unmapped);

    int global_surface = local_surface;
    int global_mapped = local_mapped;
    int global_unmapped = local_unmapped;

    monolis_allreduce_I(
        1, &global_surface, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(
        1, &global_mapped, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(
        1, &global_unmapped, MONOLIS_MPI_SUM, mono_com->comm);

    free(owner);
    free(lface);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s global %s owner-map: surf=%d mapped=%d unmapped=%d\n",
            CODENAME,
            label,
            global_surface,
            global_mapped,
            global_unmapped);
    }

    if(global_unmapped != 0) {
        if(is_physical_internal_surface) {
            /*
             * Physical/internal faces are the non-duplicated wall set used
             * for physical diagnostics.  Every one of these faces must have
             * a local owner volume element.
             */
            if(monolis_mpi_get_global_my_rank() == 0) {
                fprintf(stderr,
                    "%s ERROR: %s has %d physical/internal TRI3 face(s) "
                    "without a local owner.\n",
                    CODENAME,
                    label,
                    global_unmapped);
            }
            exit(EXIT_FAILURE);
        }
        else {
            /*
             * Assembly surfaces contain MPI overlap copies.  A copy can be
             * present on a rank that does not own the adjacent volume cell.
             * The TET4/PRI6 Nitsche assembly already treats owner=-1 as a
             * redundant copy and skips it, so this is not a mesh error.
             */
            if(monolis_mpi_get_global_my_rank() == 0) {
                fprintf(stderr,
                    "%s WARNING: %s has %d overlap TRI3 copy/copies without "
                    "a local owner; these redundant copies will be skipped.\n",
                    CODENAME,
                    label,
                    global_unmapped);
            }
        }
    }

    /*
     * Mesh-independent policy:
     *   - assembly surface: overlap/duplicate copies are allowed; every local
     *     face copy only needs a valid local owner.
     *   - physical/internal surface: the partitioner's *.n_internal.<rank>
     *     subset is used for SUM-type diagnostics.
     *
     * No fixed global face count is compared here. The wall resolution is a
     * property of the input mesh, not of the solver executable.
     */
    if(monolis_mpi_get_global_my_rank() == 0) {
        if(is_physical_internal_surface) {
            printf(
                "%s %s physical/internal surface: faces=%d "
                "(mesh-dependent count accepted)\n",
                CODENAME,
                label,
                global_surface);
        }
        else {
            printf(
                "%s %s assembly surface: copies=%d "
                "(mesh-dependent count; overlap copies allowed)\n",
                CODENAME,
                label,
                global_surface);
        }
    }
}

static int velocity_strong_mask(
    const BBFE_BC* bc_nr,
    int node_id)
{
    int mask = 0;
    if(bc_nr == NULL || bc_nr->D_bc_exists == NULL || node_id < 0) return 0;

    for(int d=0; d<3; ++d) {
        if(bc_nr->D_bc_exists[4*node_id+d]) mask |= (1 << d);
    }
    return mask;
}

/*
 * Strong/Nitsche junction policy
 * ------------------------------
 * A Nitsche wall can legitimately meet a strongly constrained boundary along
 * an edge or at a vertex (for example an Ahmed support touching the strong
 * ground boundary).  Therefore a wall node that also has one or more strong
 * velocity rows is NOT, by itself, an error.
 *
 * What is rejected is a complete TRI3 Nitsche face whose three vertices all
 * carry at least one strong velocity row.  Such a face is effectively covered
 * by another strong boundary definition and should not also be integrated as
 * an Ahmed Nitsche face.
 *
 * Use the INTERNAL/non-duplicated surfaces here.  This makes the face counts
 * physical counts rather than MPI assembly-copy counts.
 */
static void check_nitsche_wall_strong_bc_junctions(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_BC* bc_nr,
    int nnode,
    MONOLIS_COM* mono_com)
{
    unsigned char* wall_node =
        (unsigned char*)calloc((size_t)(nnode > 0 ? nnode : 1), sizeof(unsigned char));
    if(wall_node == NULL) {
        fprintf(stderr, "%s ERROR: wall-BC junction check allocation failed.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    const BBFE_DATA* surfaces[2] = {surf_pri_internal, surf_tet_internal};

    int local_faces_touching_strong = 0;
    int local_faces_all_nodes_strong = 0;
    int local_faces_all_nodes_full_velocity = 0;

    for(int s=0; s<2; ++s) {
        const BBFE_DATA* surf = surfaces[s];
        if(surf == NULL) continue;

        for(int e=0; e<surf->total_num_elems; ++e) {
            int face_any_strong = 0;
            int face_all_nodes_strong = 1;
            int face_all_nodes_full_velocity = 1;

            for(int a=0; a<surf->local_num_nodes; ++a) {
                const int i = surf->conn[e][a];
                if(i < 0 || i >= nnode) {
                    fprintf(stderr,
                        "%s ERROR: invalid wall node id=%d in strong-BC junction check.\n",
                        CODENAME, i);
                    free(wall_node);
                    exit(EXIT_FAILURE);
                }

                wall_node[i] = 1;
                const int mask = velocity_strong_mask(bc_nr, i);

                if(mask != 0) face_any_strong = 1;
                else          face_all_nodes_strong = 0;

                if(mask != 0x7) face_all_nodes_full_velocity = 0;
            }

            if(face_any_strong) ++local_faces_touching_strong;
            if(face_all_nodes_strong) ++local_faces_all_nodes_strong;
            if(face_all_nodes_full_velocity) ++local_faces_all_nodes_full_velocity;
        }
    }

    int local_intersection_node_copies = 0;
    int local_full_velocity_node_copies = 0;
    int local_component_node_copies[3] = {0,0,0};

    for(int i=0; i<nnode; ++i) {
        if(!wall_node[i]) continue;

        const int mask = velocity_strong_mask(bc_nr, i);
        if(mask == 0) continue;

        ++local_intersection_node_copies;
        if(mask == 0x7) ++local_full_velocity_node_copies;
        for(int d=0; d<3; ++d) {
            if(mask & (1 << d)) ++local_component_node_copies[d];
        }
    }

    free(wall_node);

    int global_intersection_node_copies = local_intersection_node_copies;
    int global_full_velocity_node_copies = local_full_velocity_node_copies;
    int global_component_node_copies[3] = {
        local_component_node_copies[0],
        local_component_node_copies[1],
        local_component_node_copies[2]
    };
    int global_faces_touching_strong = local_faces_touching_strong;
    int global_faces_all_nodes_strong = local_faces_all_nodes_strong;
    int global_faces_all_nodes_full_velocity = local_faces_all_nodes_full_velocity;

    monolis_allreduce_I(
        1, &global_intersection_node_copies, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(
        1, &global_full_velocity_node_copies, MONOLIS_MPI_SUM, mono_com->comm);
    for(int d=0; d<3; ++d) {
        monolis_allreduce_I(
            1, &global_component_node_copies[d], MONOLIS_MPI_SUM, mono_com->comm);
    }
    monolis_allreduce_I(
        1, &global_faces_touching_strong, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(
        1, &global_faces_all_nodes_strong, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(
        1, &global_faces_all_nodes_full_velocity, MONOLIS_MPI_SUM, mono_com->comm);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s Ahmed wall / strong-BC junctions: "
            "node_copies=%d full_velocity_node_copies=%d "
            "ux=%d uy=%d uz=%d touching_faces=%d\n",
            CODENAME,
            global_intersection_node_copies,
            global_full_velocity_node_copies,
            global_component_node_copies[0],
            global_component_node_copies[1],
            global_component_node_copies[2],
            global_faces_touching_strong);
    }

    if(global_faces_all_nodes_strong != 0) {
        fprintf(stderr,
            "%s ERROR: %d physical Ahmed TRI3 face(s) have strong velocity "
            "rows on all three vertices.\n"
            "%s ERROR: this is a face-level strong/Nitsche boundary overlap, "
            "not merely an allowed edge/vertex junction.\n",
            CODENAME,
            global_faces_all_nodes_strong,
            CODENAME);
        exit(EXIT_FAILURE);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        if(global_intersection_node_copies != 0) {
            printf(
                "%s Ahmed wall / strong-BC junction check: OK "
                "(edge/vertex sharing allowed; no fully overlapped TRI3 faces).\n",
                CODENAME);
        }
        else {
            printf("%s Ahmed wall / strong-BC junction check: OK (no shared nodes).\n",
                   CODENAME);
        }

        if(global_faces_all_nodes_full_velocity != 0) {
            /* This is redundant with the stricter all-nodes-strong check above,
             * but keep the value visible for debugging if the policy changes. */
            printf("%s note: fully strong-velocity overlapped faces=%d\n",
                   CODENAME, global_faces_all_nodes_full_velocity);
        }
    }
}

/*
 * Runtime-selectable benchmark flow.
 *
 * FLUID_CASE=ahmed (default): mixed TET4+PRI6 Ahmed-body flow.
 * FLUID_CASE=karman          : mixed TET4+PRI6 Karman-vortex volume flow,
 *                              with a single PRI6-owned TRI3 wall surface.
 *
 * Ahmed wall treatment:
 *   surf_prism_graph.dat / surf_tet_graph.dat
 *       -> overlap-capable assembly surfaces for Nitsche matrix/RHS
 *   the SAME files, truncated to N_internal via *.n_internal.<rank>
 *       -> non-duplicated physical surfaces for area/force/Cd/Cl SUMs
 *
 * Strong/Nitsche junction policy:
 *   Inlet/slip/ground outer-boundary rows remain strong.  Shared Ahmed-wall
 *   edge/vertex nodes are allowed.  What is forbidden is a complete Ahmed TRI3
 *   face that is simultaneously covered by strong velocity rows.
 */

typedef struct {
    double pressure[3];
    double viscous[3];
    /* Numerical wall-penalty reaction: reported separately, not added to total. */
    double nitsche_penalty[3];
    double total[3];
} AhmedPhysicalForce;

static void ahmed_force_zero(AhmedPhysicalForce* f)
{
    if(f != NULL) memset(f, 0, sizeof(*f));
}

static int ahmed_calc_physical_force_global(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_BASIS* basis_pri_internal,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_BASIS* basis_tet_internal,
    FE_SYSTEM_FLUID_MIXED* sys,
    int enable_pri_wall,
    int enable_tet_wall,
    const PRI6WallOptions* wall_pri,
    const PRI6WallOptions* wall_tet,
    double Uref,
    double Aref,
    AhmedPhysicalForce* force)
{
    if(force == NULL || sys == NULL || wall_pri == NULL || wall_tet == NULL) {
        return -1;
    }

    ahmed_force_zero(force);

    const double ex[3] = {1.0, 0.0, 0.0};
    const double ey[3] = {0.0, 1.0, 0.0};
    const double ez[3] = {0.0, 0.0, 1.0};
    const double Uw[3] = {0.0, 0.0, 0.0};

    /* PRI6 routine already performs its MPI reductions internally. */
    PRI6ForceDiagnostics pri;
    memset(&pri, 0, sizeof(pri));
    if(BBFE_pri6_force_diagnostics_global(
        surf_pri_internal,
        &(sys->fe_pri),
        basis_pri_internal,
        &(sys->vals),
        sys->vals.density,
        sys->vals.viscosity,
        Uw,
        wall_pri,
        &(sys->mono_com),
        &pri) != 0) {
        return -1;
    }

    for(int d=0; d<3; ++d) {
        force->pressure[d] += pri.force_pressure[d];
        force->viscous[d] += pri.force_viscous[d];
        if(enable_pri_wall && BBFE_pri6_wall_mode_uses_nitsche(wall_pri->wall_mode)) {
            force->nitsche_penalty[d] += pri.force_nitsche[d];
        }
    }

    /*
     * TET4 routines return LOCAL dimensional body-force projections.
     * Evaluate x/y together and z separately, then MPI_SUM the three axes.
     */
    double tet_p[3] = {0.0, 0.0, 0.0};
    double tet_v[3] = {0.0, 0.0, 0.0};
    double tet_nit[3] = {0.0, 0.0, 0.0};

    /* Keep the legacy TET4 force diagnostic penalty identical to the
     * split-Nitsche runtime options used by the mixed solver. */
    BBFE_fluid_set_nitsche_runtime_parameters(
        wall_tet->gamma_n, wall_tet->dt_penalty_coeff);

    calc_Cd_p_hex_or_tet(
        (BBFE_DATA*)surf_tet_internal,
        &(sys->fe_tet),
        (BBFE_BASIS*)basis_tet_internal,
        &(sys->vals),
        ex, ey,
        &tet_p[0], &tet_p[1]);

    {
        double dummy = 0.0;
        calc_Cd_p_hex_or_tet(
            (BBFE_DATA*)surf_tet_internal,
            &(sys->fe_tet),
            (BBFE_BASIS*)basis_tet_internal,
            &(sys->vals),
            ez, ex,
            &tet_p[2], &dummy);
    }

    {
        double dummy_cd = 0.0, dummy_cl = 0.0;
        double dummy_cdn = 0.0, dummy_cln = 0.0;
        if(calc_Cd_v_hex_or_tet_nitsche(
            surf_tet_internal,
            &(sys->fe_tet),
            basis_tet_internal,
            &(sys->mono_com),
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            Uref, Aref,
            ex, ey, Uw, wall_tet->gamma_n,
            &tet_v[0], &tet_v[1],
            &dummy_cd, &dummy_cl,
            &tet_nit[0], &tet_nit[1],
            &dummy_cdn, &dummy_cln) != 0) {
            return -1;
        }
    }

    {
        double dummy_l = 0.0;
        double dummy_cd = 0.0, dummy_cl = 0.0;
        double dummy_ln = 0.0;
        double dummy_cdn = 0.0, dummy_cln = 0.0;
        if(calc_Cd_v_hex_or_tet_nitsche(
            surf_tet_internal,
            &(sys->fe_tet),
            basis_tet_internal,
            &(sys->mono_com),
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            Uref, Aref,
            ez, ex, Uw, wall_tet->gamma_n,
            &tet_v[2], &dummy_l,
            &dummy_cd, &dummy_cl,
            &tet_nit[2], &dummy_ln,
            &dummy_cdn, &dummy_cln) != 0) {
            return -1;
        }
    }

    monolis_allreduce_R(3, tet_p, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_R(3, tet_v, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_R(3, tet_nit, MONOLIS_MPI_SUM, sys->mono_com.comm);

    if(!enable_tet_wall || !BBFE_pri6_wall_mode_uses_nitsche(wall_tet->wall_mode)) {
        tet_nit[0]=tet_nit[1]=tet_nit[2]=0.0;
    }

    for(int d=0; d<3; ++d) {
        force->pressure[d] += tet_p[d];
        force->viscous[d] += tet_v[d];
        force->nitsche_penalty[d] += tet_nit[d];
        force->total[d] = force->pressure[d] + force->viscous[d];
    }

    return 0;
}

static void ahmed_write_wall_runtime_config(
    const char* directory,
    int enable_pri,
    const PRI6WallOptions* pri,
    int enable_tet,
    const PRI6WallOptions* tet)
{
    if(monolis_mpi_get_global_my_rank() != 0 || pri == NULL || tet == NULL) return;
    char path[4096];
    if(directory == NULL || directory[0] == '\0' || strcmp(directory,".") == 0) {
        snprintf(path,sizeof(path),"ahmed_wall_runtime_config.txt");
    }
    else {
        const size_t n = strlen(directory);
        snprintf(path,sizeof(path),(n > 0 && directory[n-1] == '/') ? "%s%s" : "%s/%s",
                 directory,"ahmed_wall_runtime_config.txt");
    }
    FILE* fp = fopen(path,"w");
    if(fp == NULL) {
        fprintf(stderr,"%s WARNING: cannot write %s\n",CODENAME,path);
        return;
    }
    fprintf(fp,"PRI6 enabled %d\n",enable_pri ? 1 : 0);
    fprintf(fp,"PRI6 mode %s\n",enable_pri ? BBFE_pri6_wall_mode_name(pri->wall_mode) : "off");
    fprintf(fp,"PRI6 exchange %s\n",BBFE_pri6_wall_exchange_name(pri->exchange_mode));
    fprintf(fp,"PRI6 gamma %.17e\n",pri->gamma_n);
    fprintf(fp,"PRI6 c_dt %.17e\n",pri->dt_penalty_coeff);
    fprintf(fp,"PRI6 kappa %.17e\n",pri->kappa);
    fprintf(fp,"PRI6 B %.17e\n",pri->B);
    fprintf(fp,"PRI6 ut_eps %.17e\n",pri->ut_eps);
    fprintf(fp,"TET4 enabled %d\n",enable_tet ? 1 : 0);
    fprintf(fp,"TET4 mode %s\n",enable_tet ? BBFE_pri6_wall_mode_name(tet->wall_mode) : "off");
    fprintf(fp,"TET4 exchange %s\n",BBFE_pri6_wall_exchange_name(tet->exchange_mode));
    fprintf(fp,"TET4 gamma %.17e\n",tet->gamma_n);
    fprintf(fp,"TET4 c_dt %.17e\n",tet->dt_penalty_coeff);
    fprintf(fp,"TET4 kappa %.17e\n",tet->kappa);
    fprintf(fp,"TET4 B %.17e\n",tet->B);
    fprintf(fp,"TET4 ut_eps %.17e\n",tet->ut_eps);
    fprintf(fp,"projected_strong_feature_angle_deg %.17e\n",
            BBFE_mixed_wall_get_feature_angle_deg());
    fclose(fp);
}

static void ahmed_force_output_path(
    char* path, size_t path_size, const char* directory, const char* filename)
{
    if(path == NULL || path_size == 0) return;
    if(directory == NULL || directory[0] == '\0' || strcmp(directory, ".") == 0) {
        snprintf(path, path_size, "%s", filename);
        return;
    }
    const size_t n = strlen(directory);
    snprintf(path, path_size,
        (n > 0 && directory[n-1] == '/') ? "%s%s" : "%s/%s",
        directory, filename);
}

static void ahmed_reset_force_coefficients_file(const char* directory)
{
    char path[4096];
    ahmed_force_output_path(path, sizeof(path), directory,
                            "ahmed_force_coefficients.csv");
    FILE* fp = fopen(path, "w");
    if(fp == NULL) {
        fprintf(stderr, "%s WARNING: cannot reset %s\n", CODENAME, path);
        return;
    }
    fprintf(fp,
        "time,step,Uref,Aref,qA,"
        "Cd,Cy,Cl,"
        "Cd_pressure,Cy_pressure,Cl_pressure,"
        "Cd_viscous,Cy_viscous,Cl_viscous,"
        "Cd_nitsche_penalty,Cy_nitsche_penalty,Cl_nitsche_penalty,"
        "Fx,Fy,Fz,Fx_pressure,Fy_pressure,Fz_pressure,"
        "Fx_viscous,Fy_viscous,Fz_viscous,"
        "Fx_nitsche_penalty,Fy_nitsche_penalty,Fz_nitsche_penalty\n");
    fclose(fp);
}

static void ahmed_output_force_coefficients(
    const AhmedPhysicalForce* f,
    double t, int step, double rho, double Uref, double Aref,
    const char* directory)
{
    if(f == NULL) return;
    const double qA = 0.5 * rho * Uref * Uref * Aref;
    if(!isfinite(qA) || qA <= 0.0) {
        fprintf(stderr, "%s ERROR: invalid Ahmed force reference qA=%e\n",
                CODENAME, qA);
        return;
    }

    char path[4096];
    ahmed_force_output_path(path, sizeof(path), directory,
                            "ahmed_force_coefficients.csv");
    FILE* fp = fopen(path, "a");
    if(fp == NULL) {
        fprintf(stderr, "%s WARNING: cannot append %s\n", CODENAME, path);
        return;
    }

    fprintf(fp,
        "%.17e,%d,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e\n",
        t, step, Uref, Aref, qA,
        f->total[0]/qA, f->total[1]/qA, f->total[2]/qA,
        f->pressure[0]/qA, f->pressure[1]/qA, f->pressure[2]/qA,
        f->viscous[0]/qA, f->viscous[1]/qA, f->viscous[2]/qA,
        f->nitsche_penalty[0]/qA,
        f->nitsche_penalty[1]/qA,
        f->nitsche_penalty[2]/qA,
        f->total[0], f->total[1], f->total[2],
        f->pressure[0], f->pressure[1], f->pressure[2],
        f->viscous[0], f->viscous[1], f->viscous[2],
        f->nitsche_penalty[0],
        f->nitsche_penalty[1],
        f->nitsche_penalty[2]);
    fclose(fp);
}

static void ahmed_reset_wall_summary_csv(const char* directory)
{
    char path[4096];
    ahmed_force_output_path(path, sizeof(path), directory,
                            "ahmed_wall_summary.csv");
    FILE* fp = fopen(path, "w");
    if(fp == NULL) {
        fprintf(stderr, "%s WARNING: cannot reset %s\n", CODENAME, path);
        return;
    }
    fprintf(fp,
        "time,step,family,faces,qps,area,"
        "h_min,h_mean,h_max,"
        "tauw_min,tauw_mean,tauw_max,"
        "utau_min,utau_mean,utau_max,"
        "yplus_min,yplus_mean,yplus_max,"
        "abs_un_mean,abs_un_max,ut_mean,ut_max,"
        "p_ref,p_min,p_mean,p_max,cp_min,cp_mean,cp_max\n");
    fclose(fp);
}

static void ahmed_write_wall_summary_row(
    FILE* fp,
    double t,
    int step,
    const char* family,
    const MixedWallFamilyDiagnostics* d,
    double p_ref,
    double q_ref)
{
    if(fp == NULL || family == NULL || d == NULL) return;
    const double cp_min = (q_ref > 0.0) ? (d->pressure_min-p_ref)/q_ref : 0.0;
    const double cp_mean = (q_ref > 0.0) ? (d->pressure_mean-p_ref)/q_ref : 0.0;
    const double cp_max = (q_ref > 0.0) ? (d->pressure_max-p_ref)/q_ref : 0.0;
    fprintf(fp,
        "%.17e,%d,%s,%lld,%lld,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
        t, step, family, d->faces, d->quadrature_points, d->area,
        d->h_n_min, d->h_n_mean, d->h_n_max,
        d->tau_w_min, d->tau_w_mean, d->tau_w_max,
        d->u_tau_min, d->u_tau_mean, d->u_tau_max,
        d->y_plus_min, d->y_plus_mean, d->y_plus_max,
        d->abs_u_normal_mean, d->abs_u_normal_max,
        d->u_tangent_mean, d->u_tangent_max,
        p_ref, d->pressure_min, d->pressure_mean, d->pressure_max,
        cp_min, cp_mean, cp_max);
}

static void ahmed_output_wall_summary_csv(
    const MixedWallDiagnostics* diag,
    double t,
    int step,
    double p_ref,
    double rho,
    double Uref,
    const char* directory)
{
    if(diag == NULL) return;

    char path[4096];
    ahmed_force_output_path(path, sizeof(path), directory,
                            "ahmed_wall_summary.csv");
    FILE* fp = fopen(path, "a");
    if(fp == NULL) {
        fprintf(stderr, "%s WARNING: cannot append %s\n", CODENAME, path);
        return;
    }

    const double q_ref = 0.5*rho*Uref*Uref;
    ahmed_write_wall_summary_row(
        fp, t, step, "PRI6", &diag->pri6, p_ref, q_ref);
    ahmed_write_wall_summary_row(
        fp, t, step, "TET4", &diag->tet4, p_ref, q_ref);
    ahmed_write_wall_summary_row(
        fp, t, step, "TOTAL", &diag->total, p_ref, q_ref);
    fclose(fp);
}


static double ahmed_pressure_reference_global(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int* inlet_node_count_out)
{
    if(inlet_node_count_out != NULL) *inlet_node_count_out = 0;
    if(sys == NULL) return 0.0;

    const int nnode = sys->fe_tet.total_num_nodes;
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int nint = (comm_size == 1) ? nnode : sys->mono_com.n_internal_vertex;
    if(nint < 0 || nint > nnode) return 0.0;

    double psum = 0.0;
    int count = 0;
    for(int i=0; i<nint; ++i) {
        const int dof = 4*i;
        if(sys->bc.D_bc_exists != NULL &&
           sys->bc.imposed_D_val != NULL &&
           sys->bc.D_bc_exists[dof] &&
           fabs(sys->bc.imposed_D_val[dof]) > 1.0e-12) {
            psum += sys->vals.p[i];
            ++count;
        }
    }

    if(comm_size != 1) {
        monolis_allreduce_R(1,&psum,MONOLIS_MPI_SUM,sys->mono_com.comm);
        monolis_allreduce_I(1,&count,MONOLIS_MPI_SUM,sys->mono_com.comm);
    }

    /* Fallback: global internal-node mean pressure if no nonzero-ux inlet row. */
    if(count == 0) {
        psum = 0.0;
        count = 0;
        for(int i=0; i<nint; ++i) {
            psum += sys->vals.p[i];
            ++count;
        }
        if(comm_size != 1) {
            monolis_allreduce_R(1,&psum,MONOLIS_MPI_SUM,sys->mono_com.comm);
            monolis_allreduce_I(1,&count,MONOLIS_MPI_SUM,sys->mono_com.comm);
        }
    }

    if(inlet_node_count_out != NULL) *inlet_node_count_out = count;
    return (count > 0) ? psum/(double)count : 0.0;
}

typedef struct {
    long long tet_elements;
    long long pri_elements;
    double volume;
    double int_div_u;
    double int_abs_div_u;
    double int_div_u2;
    double max_abs_div_u;
    double mean_abs_div_u;
    double rms_div_u;
    double invalid_elements;
} AhmedContinuityDiagnostics;

static int ahmed_read_internal_element_count(
    const char* directory,
    const char* filename_base,
    int* n_internal_out)
{
    if(filename_base == NULL || n_internal_out == NULL) return -1;
    *n_internal_out = -1;

    char rel[4096];
    snprintf(rel,sizeof(rel),"parted.0/%s.n_internal.%d",
             filename_base,monolis_mpi_get_global_my_rank());
    char path[4096];
    ahmed_force_output_path(path,sizeof(path),directory,rel);

    FILE* fp = fopen(path,"r");
    if(fp == NULL) return -1;

    char label[4096];
    int aux = 0;
    int n_internal = -1;
    const int ok =
        (fscanf(fp,"%4095s %d",label,&aux) == 2 &&
         fscanf(fp,"%d",&n_internal) == 1 && n_internal >= 0);
    fclose(fp);
    if(!ok) return -1;

    *n_internal_out = n_internal;
    return 0;
}

static int ahmed_continuity_accumulate_family(
    const BBFE_DATA* fe,
    const BBFE_BASIS* basis,
    const VALUES* vals,
    int n_internal,
    AhmedContinuityDiagnostics* d,
    int family_is_pri)
{
    if(fe == NULL || basis == NULL || vals == NULL || d == NULL) return -1;
    if(n_internal < 0 || n_internal > fe->total_num_elems) return -1;
    if(n_internal == 0) return 0;
    if(fe->local_num_nodes <= 0 || basis->num_integ_points <= 0) return -1;

    if(family_is_pri) d->pri_elements += n_internal;
    else d->tet_elements += n_internal;

    const int nl = fe->local_num_nodes;
    const int np = basis->num_integ_points;

    for(int e=0; e<n_internal; ++e) {
        int invalid_elem = 0;
        for(int p=0; p<np; ++p) {
            const double J = fabs(fe->geo[e][p].Jacobian);
            const double w = basis->integ_weight[p] * J;
            if(!isfinite(w) || w <= 0.0) {
                invalid_elem = 1;
                continue;
            }

            double div_u = 0.0;
            for(int a=0; a<nl; ++a) {
                const int gid = fe->conn[e][a];
                div_u +=
                    vals->v[gid][0]*fe->geo[e][p].grad_N[a][0] +
                    vals->v[gid][1]*fe->geo[e][p].grad_N[a][1] +
                    vals->v[gid][2]*fe->geo[e][p].grad_N[a][2];
            }

            if(!isfinite(div_u)) {
                invalid_elem = 1;
                continue;
            }

            d->volume += w;
            d->int_div_u += div_u*w;
            d->int_abs_div_u += fabs(div_u)*w;
            d->int_div_u2 += div_u*div_u*w;
            d->max_abs_div_u = fmax(d->max_abs_div_u,fabs(div_u));
        }
        if(invalid_elem) d->invalid_elements += 1.0;
    }
    return 0;
}

static int ahmed_continuity_global(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int n_internal_tet,
    int n_internal_pri,
    AhmedContinuityDiagnostics* d)
{
    if(sys == NULL || d == NULL) return -1;
    memset(d,0,sizeof(*d));

    if(ahmed_continuity_accumulate_family(
        &sys->fe_tet,&sys->basis_tet,&sys->vals,n_internal_tet,d,0) != 0) {
        return -1;
    }
    if(ahmed_continuity_accumulate_family(
        &sys->fe_pri,&sys->basis_pri,&sys->vals,n_internal_pri,d,1) != 0) {
        return -1;
    }

    int nt = (int)d->tet_elements;
    int np = (int)d->pri_elements;
    double sum[5] = {
        d->volume,d->int_div_u,d->int_abs_div_u,d->int_div_u2,d->invalid_elements
    };
    double max_div = d->max_abs_div_u;

    if(monolis_mpi_get_global_comm_size() != 1) {
        monolis_allreduce_I(1,&nt,MONOLIS_MPI_SUM,sys->mono_com.comm);
        monolis_allreduce_I(1,&np,MONOLIS_MPI_SUM,sys->mono_com.comm);
        monolis_allreduce_R(5,sum,MONOLIS_MPI_SUM,sys->mono_com.comm);
        monolis_allreduce_R(1,&max_div,MONOLIS_MPI_MAX,sys->mono_com.comm);
    }

    d->tet_elements = nt;
    d->pri_elements = np;
    d->volume = sum[0];
    d->int_div_u = sum[1];
    d->int_abs_div_u = sum[2];
    d->int_div_u2 = sum[3];
    d->invalid_elements = sum[4];
    d->max_abs_div_u = max_div;
    if(d->volume > 0.0) {
        d->mean_abs_div_u = d->int_abs_div_u/d->volume;
        d->rms_div_u = sqrt(fmax(d->int_div_u2/d->volume,0.0));
    }
    return 0;
}

static void ahmed_reset_continuity_csv(const char* directory)
{
    char path[4096];
    ahmed_force_output_path(path,sizeof(path),directory,"ahmed_continuity.csv");
    FILE* fp = fopen(path,"w");
    if(fp == NULL) return;
    fprintf(fp,
        "time,step,tet_internal_elements,pri_internal_elements,volume,"
        "int_div_u,int_abs_div_u,mean_abs_div_u,rms_div_u,max_abs_div_u,"
        "net_mass_flux,abs_net_mass_flux,net_mass_flux_over_rhoUrefAref,"
        "invalid_elements\n");
    fclose(fp);
}

static void ahmed_output_continuity_csv(
    const AhmedContinuityDiagnostics* d,
    double t,
    int step,
    double rho,
    double Uref,
    double Aref,
    const char* directory)
{
    if(d == NULL) return;
    const double mdot = rho*d->int_div_u;
    const double mdot_ref = rho*fabs(Uref)*Aref;
    const double rel = (mdot_ref > 0.0) ? fabs(mdot)/mdot_ref : 0.0;

    char path[4096];
    ahmed_force_output_path(path,sizeof(path),directory,"ahmed_continuity.csv");
    FILE* fp = fopen(path,"a");
    if(fp == NULL) return;
    fprintf(fp,
        "%.17e,%d,%lld,%lld,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e\n",
        t,step,d->tet_elements,d->pri_elements,d->volume,
        d->int_div_u,d->int_abs_div_u,d->mean_abs_div_u,
        d->rms_div_u,d->max_abs_div_u,
        mdot,fabs(mdot),rel,d->invalid_elements);
    fclose(fp);
}

static void reduce_min_positive(double local_value, double* global_value, MONOLIS_COM* mono_com)
{
    double neg = (isfinite(local_value) && local_value < DBL_MAX)
        ? -local_value : -DBL_MAX;
    monolis_allreduce_R(1, &neg, MONOLIS_MPI_MAX, mono_com->comm);
    *global_value = (neg > -DBL_MAX) ? -neg : 0.0;
}

static void report_volume_bc_baseline_diagnostics(
    FE_SYSTEM_FLUID_MIXED* sys,
    const char* case_name)
{
    if(sys == NULL) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    const int nint = sys->mono_com.n_internal_vertex;
    if(nnode <= 0 || nint < 0 || nint > nnode) return;

    unsigned char* used_tet = (unsigned char*)calloc((size_t)nnode, 1u);
    unsigned char* used_pri = (unsigned char*)calloc((size_t)nnode, 1u);
    if(used_tet == NULL || used_pri == NULL) {
        fprintf(stderr, "%s ERROR: baseline diagnostic allocation failed.\n", CODENAME);
        free(used_tet); free(used_pri);
        exit(EXIT_FAILURE);
    }

    for(int e=0; e<sys->fe_tet.total_num_elems; ++e) {
        for(int a=0; a<sys->fe_tet.local_num_nodes; ++a) {
            const int i = sys->fe_tet.conn[e][a];
            if(i >= 0 && i < nnode) used_tet[i] = 1;
        }
    }
    for(int e=0; e<sys->fe_pri.total_num_elems; ++e) {
        for(int a=0; a<sys->fe_pri.local_num_nodes; ++a) {
            const int i = sys->fe_pri.conn[e][a];
            if(i >= 0 && i < nnode) used_pri[i] = 1;
        }
    }

    int membership[4] = {0,0,0,0}; /* unused, tet-only, pri-only, shared */
    int bc_rows[4] = {0,0,0,0};
    int bc_nonzero[4] = {0,0,0,0};
    double bc_max[4] = {0.0,0.0,0.0,0.0};
    int velocity_mask_count[4] = {0,0,0,0}; /* 0,1,2,3 constrained components */

    for(int i=0; i<nint; ++i) {
        const int mt = used_tet[i] ? 1 : 0;
        const int mp = used_pri[i] ? 1 : 0;
        const int cls = mt ? (mp ? 3 : 1) : (mp ? 2 : 0);
        ++membership[cls];

        int nvc = 0;
        for(int d=0; d<4; ++d) {
            const int dof = 4*i + d;
            if(sys->bc_NR.D_bc_exists != NULL && sys->bc_NR.D_bc_exists[dof]) {
                ++bc_rows[d];
                if(d < 3) ++nvc;
                if(sys->bc.imposed_D_val != NULL) {
                    const double av = fabs(sys->bc.imposed_D_val[dof]);
                    if(av > 0.0) ++bc_nonzero[d];
                    if(av > bc_max[d]) bc_max[d] = av;
                }
            }
        }
        if(nvc >= 0 && nvc <= 3) ++velocity_mask_count[nvc];
    }

    monolis_allreduce_I(4, membership, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_I(4, bc_rows, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_I(4, bc_nonzero, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_R(4, bc_max, MONOLIS_MPI_MAX, sys->mono_com.comm);
    monolis_allreduce_I(4, velocity_mask_count, MONOLIS_MPI_SUM, sys->mono_com.comm);

    double local_tet_vmin = DBL_MAX, local_tet_vmax = 0.0;
    double local_pri_vmin = DBL_MAX, local_pri_vmax = 0.0;
    double local_tet_jmin = DBL_MAX, local_tet_jmax = 0.0;
    double local_pri_jmin = DBL_MAX, local_pri_jmax = 0.0;
    long long local_tet_count = 0, local_pri_count = 0;

    if(sys->fe_tet.total_num_elems > 0 && sys->basis_tet.num_integ_points > 0) {
        const int np = sys->basis_tet.num_integ_points;
        double* jac = (double*)malloc(sizeof(double)*(size_t)np);
        if(jac == NULL) exit(EXIT_FAILURE);
        for(int e=0; e<sys->fe_tet.total_num_elems; ++e) {
            for(int q=0; q<np; ++q) {
                const double J = sys->fe_tet.geo[e][q].Jacobian;
                jac[q] = J;
                if(isfinite(J)) {
                    const double aJ = fabs(J);
                    if(aJ < local_tet_jmin) local_tet_jmin = aJ;
                    if(aJ > local_tet_jmax) local_tet_jmax = aJ;
                }
            }
            const double vol = fabs(BBFE_std_integ_calc_volume(
                np, sys->basis_tet.integ_weight, jac));
            if(isfinite(vol) && vol > 0.0) {
                if(vol < local_tet_vmin) local_tet_vmin = vol;
                if(vol > local_tet_vmax) local_tet_vmax = vol;
            }
            ++local_tet_count;
        }
        free(jac);
    }

    if(sys->fe_pri.total_num_elems > 0 && sys->basis_pri.num_integ_points > 0) {
        const int np = sys->basis_pri.num_integ_points;
        double* jac = (double*)malloc(sizeof(double)*(size_t)np);
        if(jac == NULL) exit(EXIT_FAILURE);
        for(int e=0; e<sys->fe_pri.total_num_elems; ++e) {
            for(int q=0; q<np; ++q) {
                const double J = sys->fe_pri.geo[e][q].Jacobian;
                jac[q] = J;
                if(isfinite(J)) {
                    const double aJ = fabs(J);
                    if(aJ < local_pri_jmin) local_pri_jmin = aJ;
                    if(aJ > local_pri_jmax) local_pri_jmax = aJ;
                }
            }
            const double vol = fabs(BBFE_std_integ_calc_volume(
                np, sys->basis_pri.integ_weight, jac));
            if(isfinite(vol) && vol > 0.0) {
                if(vol < local_pri_vmin) local_pri_vmin = vol;
                if(vol > local_pri_vmax) local_pri_vmax = vol;
            }
            ++local_pri_count;
        }
        free(jac);
    }

    double tet_vmin=0.0, pri_vmin=0.0, tet_jmin=0.0, pri_jmin=0.0;
    reduce_min_positive(local_tet_vmin, &tet_vmin, &(sys->mono_com));
    reduce_min_positive(local_pri_vmin, &pri_vmin, &(sys->mono_com));
    reduce_min_positive(local_tet_jmin, &tet_jmin, &(sys->mono_com));
    reduce_min_positive(local_pri_jmin, &pri_jmin, &(sys->mono_com));

    double tet_vmax=local_tet_vmax, pri_vmax=local_pri_vmax;
    double tet_jmax=local_tet_jmax, pri_jmax=local_pri_jmax;
    monolis_allreduce_R(1, &tet_vmax, MONOLIS_MPI_MAX, sys->mono_com.comm);
    monolis_allreduce_R(1, &pri_vmax, MONOLIS_MPI_MAX, sys->mono_com.comm);
    monolis_allreduce_R(1, &tet_jmax, MONOLIS_MPI_MAX, sys->mono_com.comm);
    monolis_allreduce_R(1, &pri_jmax, MONOLIS_MPI_MAX, sys->mono_com.comm);

    /* Element counts are intentionally long long locally but reduced as doubles
     * only for a lightweight diagnostic; values here are exactly representable. */
    double tet_count = (double)local_tet_count;
    double pri_count = (double)local_pri_count;
    monolis_allreduce_R(1, &tet_count, MONOLIS_MPI_SUM, sys->mono_com.comm);
    monolis_allreduce_R(1, &pri_count, MONOLIS_MPI_SUM, sys->mono_com.comm);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s baseline conditions (%s): rho=%.9e mu=%.9e dt=%.9e "
            "finish=%.9e ip_axis=%d mat_max_iter=%d mat_eps=%.9e\n",
            CODENAME, case_name ? case_name : "unknown",
            sys->vals.density, sys->vals.viscosity, sys->vals.dt,
            sys->vals.finish_time, sys->vals.num_ip_each_axis,
            sys->vals.mat_max_iter, sys->vals.mat_epsilon);
        printf(
            "%s owned node-family connectivity: unused=%d TET4-only=%d "
            "PRI6-only=%d TET4+PRI6(shared)=%d\n",
            CODENAME, membership[0], membership[1], membership[2], membership[3]);
        printf(
            "%s owned strong-BC rows: ux=%d uy=%d uz=%d p=%d; "
            "nonzero=(%d,%d,%d,%d); max|value|=(%.6e,%.6e,%.6e,%.6e)\n",
            CODENAME,
            bc_rows[0], bc_rows[1], bc_rows[2], bc_rows[3],
            bc_nonzero[0], bc_nonzero[1], bc_nonzero[2], bc_nonzero[3],
            bc_max[0], bc_max[1], bc_max[2], bc_max[3]);
        printf(
            "%s owned velocity-constraint masks: 0comp=%d 1comp=%d 2comp=%d 3comp=%d\n",
            CODENAME,
            velocity_mask_count[0], velocity_mask_count[1],
            velocity_mask_count[2], velocity_mask_count[3]);
        printf(
            "%s TET4 scale: elem_copies=%.0f |J|=[%.6e,%.6e] volume=[%.6e,%.6e] "
            "h=cbrt(V)=[%.6e,%.6e] Vratio=%.6e Jratio=%.6e\n",
            CODENAME, tet_count,
            tet_jmin, tet_jmax, tet_vmin, tet_vmax,
            tet_vmin>0.0?cbrt(tet_vmin):0.0,
            tet_vmax>0.0?cbrt(tet_vmax):0.0,
            tet_vmin>0.0?tet_vmax/tet_vmin:0.0,
            tet_jmin>0.0?tet_jmax/tet_jmin:0.0);
        printf(
            "%s PRI6 scale: elem_copies=%.0f |J|=[%.6e,%.6e] volume=[%.6e,%.6e] "
            "h=cbrt(V)=[%.6e,%.6e] Vratio=%.6e Jratio=%.6e\n",
            CODENAME, pri_count,
            pri_jmin, pri_jmax, pri_vmin, pri_vmax,
            pri_vmin>0.0?cbrt(pri_vmin):0.0,
            pri_vmax>0.0?cbrt(pri_vmax):0.0,
            pri_vmin>0.0?pri_vmax/pri_vmin:0.0,
            pri_jmin>0.0?pri_jmax/pri_jmin:0.0);

        if(membership[3] == 0 && tet_count > 0.0 && pri_count > 0.0) {
            fprintf(stderr,
                "%s WARNING: no owned node is shared by TET4 and PRI6. "
                "The two volume families may be disconnected.\n", CODENAME);
        }
        if(bc_rows[3] > 0) {
            fprintf(stderr,
                "%s NOTE: pressure Dirichlet rows are present (owned p rows=%d).\n",
                CODENAME, bc_rows[3]);
        }
    }

    free(used_tet);
    free(used_pri);
}

int main(int argc, char* argv[])
{
    FE_SYSTEM_FLUID_MIXED sys;
    memset(&sys, 0, sizeof(sys));

    /* Assembly surfaces: overlap copies are allowed. */
    BBFE_DATA surf_pri_assembly;
    BBFE_BASIS basis_surf_pri_assembly;
    BBFE_DATA surf_tet_assembly;
    BBFE_BASIS basis_surf_tet_assembly;

    /* Physical/internal surfaces: globally non-duplicated. */
    BBFE_DATA surf_pri_internal;
    BBFE_BASIS basis_surf_pri_internal;
    BBFE_DATA surf_tet_internal;
    BBFE_BASIS basis_surf_tet_internal;

    memset(&surf_pri_assembly, 0, sizeof(surf_pri_assembly));
    memset(&basis_surf_pri_assembly, 0, sizeof(basis_surf_pri_assembly));
    memset(&surf_tet_assembly, 0, sizeof(surf_tet_assembly));
    memset(&basis_surf_tet_assembly, 0, sizeof(basis_surf_tet_assembly));
    memset(&surf_pri_internal, 0, sizeof(surf_pri_internal));
    memset(&basis_surf_pri_internal, 0, sizeof(basis_surf_pri_internal));
    memset(&surf_tet_internal, 0, sizeof(surf_tet_internal));
    memset(&basis_surf_tet_internal, 0, sizeof(basis_surf_tet_internal));

    monolis_global_initialize();

    const FluidCaseMode case_mode = fluid_case_mode_from_env();
    const int is_ahmed = (case_mode == FLUID_CASE_AHMED);

    /* Surface file names are runtime-selectable and mesh/pipeline independent. */
    const char* surf_pri_filename = fluid_env_string_or_default(
        "AHMED_SURF_PRI_FILE", INPUT_FILENAME_SURF_PRI);
    const char* surf_tet_filename = fluid_env_string_or_default(
        "AHMED_SURF_TET_FILE", INPUT_FILENAME_SURF_TET);
    const char* surf_single_filename = fluid_env_string_or_default(
        "KARMAN_SURF_FILE", INPUT_FILENAME_SURF_SINGLE);

    /*
     * Read split-wall switches before the initial condition is constructed.
     * The same values are reused later by the Nitsche solver.
     */
    int enable_pri6_nitsche = 1;
    int enable_tet4_nitsche = is_ahmed ? 1 : 0;
    int wall_mode_pri = PRI6_WALL_NITSCHE_NOSLIP;
    int wall_mode_tet = PRI6_WALL_NITSCHE_NOSLIP;
    double ahmed_velocity_scale = 1.0;

    if(is_ahmed) {
        const int legacy_pri = ahmed_env_flag("AHMED_NITSCHE_PRI", 1);
        const int legacy_tet = ahmed_env_flag("AHMED_NITSCHE_TET", 1);
        wall_mode_pri = ahmed_wall_mode_setting(
            "AHMED_WALL_MODE_PRI","AHMED_WALL_MODE",legacy_pri,&enable_pri6_nitsche);
        wall_mode_tet = ahmed_wall_mode_setting(
            "AHMED_WALL_MODE_TET","AHMED_WALL_MODE",legacy_tet,&enable_tet4_nitsche);
        ahmed_velocity_scale = ahmed_env_positive_double(
            "AHMED_VELOCITY_SCALE", 1.0, 0);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf("%s runtime case: FLUID_CASE=%s\n",
               CODENAME, fluid_case_mode_name(case_mode));
        if(is_ahmed) {
            printf("%s surface mode: Ahmed dual TRI3 surfaces "
                   "(%s + %s; mesh-dependent counts accepted).\n",
                   CODENAME, surf_pri_filename, surf_tet_filename);
        }
        else {
            printf("%s surface mode: Karman mixed-volume / PRI6-wall TRI3 surface "
                   "(%s only).\n", CODENAME, surf_single_filename);
        }
    }

    sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);
    read_calc_conditions(&(sys.vals), sys.cond.directory);

    /*
     * Match the established solver_fom_VMS baseline explicitly.  The original
     * Ahmed VMS main set both coefficients to zero before entering the time
     * loop.  Keep runtime overrides for controlled VMS sensitivity tests.
     */
    const char* env_C_vms =
        is_ahmed ? "AHMED_C_VMS" : "KARMAN_C_VMS";
    const char* env_vms_cap =
        is_ahmed ? "AHMED_VMS_CAP_COEFF" : "KARMAN_VMS_CAP_COEFF";

    sys.vals.C_vms = ahmed_env_positive_double(
        env_C_vms, 0.0, 1);
    sys.vals.vms_cap_coeff = ahmed_env_positive_double(
        env_vms_cap, 0.0, 1);

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s mixed VMS parameters: C_vms=%.6g vms_cap_coeff=%.6g\n",
            CODENAME,
            sys.vals.C_vms,
            sys.vals.vms_cap_coeff);

        /*
         * Current core_NR VMS semantics use
         *   nu_cap = cap_coeff * h_e * |u|
         * and clamp nu_vms to nu_cap.  Therefore C_vms>0 with cap_coeff=0
         * disables the added VMS viscosity rather than meaning "no cap".
         */
        if(sys.vals.C_vms > 0.0 && sys.vals.vms_cap_coeff <= 0.0) {
            fprintf(
                stderr,
                "%s WARNING: %s>0 but %s=0; "
                "the current VMS kernel clamps nu_vms to zero. "
                "Set a positive cap coefficient to activate VMS viscosity.\n",
                CODENAME, env_C_vms, env_vms_cap);
        }
    }

    BBFE_fluid_pre_mixed(
        &(sys.fe_tet), &(sys.basis_tet),
        &(sys.fe_pri), &(sys.basis_pri),
        argc, argv,
        sys.cond.directory,
        sys.vals.num_ip_each_axis);

    int global_tet_elems_case = sys.fe_tet.total_num_elems;
    int global_pri_elems_case = sys.fe_pri.total_num_elems;
    monolis_allreduce_I(
        1, &global_tet_elems_case, MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());
    monolis_allreduce_I(
        1, &global_pri_elems_case, MONOLIS_MPI_SUM,
        monolis_mpi_get_global_comm());

    if(case_mode == FLUID_CASE_KARMAN) {
        /*
         * Karman uses the SAME mixed TET4+PRI6 volume discretization as the
         * mixed solver.  The case distinction is only in the wall-surface
         * input: surf_graph.dat is the single PRI6-owned TRI3 wall surface.
         * Ranks are allowed to contain TET4 only, PRI6 only, or both.
         */
        if(global_tet_elems_case <= 0 && global_pri_elems_case <= 0) {
            if(monolis_mpi_get_global_my_rank() == 0) {
                fprintf(stderr,
                    "%s ERROR: FLUID_CASE=karman has no volume elements: "
                    "global TET4=%d PRI6=%d.\n",
                    CODENAME, global_tet_elems_case, global_pri_elems_case);
            }
            exit(EXIT_FAILURE);
        }

        if(global_pri_elems_case <= 0) {
            if(monolis_mpi_get_global_my_rank() == 0) {
                fprintf(stderr,
                    "%s ERROR: FLUID_CASE=karman requires PRI6 elements for "
                    "the surf_graph.dat wall ownership: global TET4=%d PRI6=%d.\n",
                    CODENAME, global_tet_elems_case, global_pri_elems_case);
            }
            exit(EXIT_FAILURE);
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "%s Karman volume check: mixed TET4+PRI6 enabled "
                "(global TET4=%d PRI6=%d); wall owner family=PRI6.\n",
                CODENAME, global_tet_elems_case, global_pri_elems_case);
        }
    }

    const int nnode = sys.fe_tet.total_num_nodes;

    memory_allocation_nodal_values_AB2(&(sys.vals), nnode);

    const char* filename = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        INPUT_FILENAME_D_BC_V);

    BBFE_fluid_sups_read_Dirichlet_bc_perturbation(
        &(sys.bc), &(sys.vals), filename, sys.cond.directory, nnode, 4);

    BBFE_fluid_sups_read_Dirichlet_bc_NR(
        &(sys.bc_NR), filename, sys.cond.directory, nnode, 4);

    /*
     * Ahmed-only Reynolds/conditioning A/B control.  Scale the absolute
     * strong velocity values from D_bc_v.dat; bc_NR intentionally remains
     * zero on constrained Newton increments.  Zero wall/slip values are
     * unchanged, while the 60 m/s inlet becomes 60*scale.
     */
    if(is_ahmed && fabs(ahmed_velocity_scale - 1.0) > 1.0e-15) {
        scale_absolute_velocity_dirichlet_bc(
            &(sys.bc), nnode, ahmed_velocity_scale, "Ahmed absolute");
    }
    else if(is_ahmed && monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s Ahmed velocity scale: AHMED_VELOCITY_SCALE=1 "
            "(D_bc_v.dat unchanged).\n",
            CODENAME);
    }

    /* ------------------------------------------------------------ */
    /* Volume geometry.                                             */
    /* ------------------------------------------------------------ */

    if(BBFE_fluid_mixed_has_elements(&(sys.fe_tet))) {
        BBFE_elemmat_set_Jacobi_mat(&(sys.fe_tet), &(sys.basis_tet));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_tet), &(sys.basis_tet));
    }

    if(BBFE_fluid_mixed_has_elements(&(sys.fe_pri))) {
        BBFE_elemmat_set_Jacobi_mat(&(sys.fe_pri), &(sys.basis_pri));
        BBFE_elemmat_set_shapefunc_derivative(&(sys.fe_pri), &(sys.basis_pri));
    }

    BBFE_fluid_mixed_check_positive_jacobian(
        &(sys.fe_tet), &(sys.basis_tet), "TET4", 0.0);

    BBFE_fluid_mixed_check_positive_jacobian(
        &(sys.fe_pri), &(sys.basis_pri), "PRI6", 0.0);

    /* ------------------------------------------------------------ */
    /* Unified TET4 + PRI6 Monolis graph.                           */
    /* ------------------------------------------------------------ */

    BBFE_fluid_init_monomat_mixed(
        &(sys.monolis), &(sys.mono_com),
        &(sys.fe_tet), &(sys.fe_pri), sys.cond.directory);

    /*
     * Pressure gauge is case-specific and OFF by default:
     *   AHMED_PRESSURE_GAUGE=0 or KARMAN_PRESSURE_GAUGE=0.
     * A manually injected rank-local pressure row changed
     * the NR=0 linear system in earlier tests; enable it only deliberately.
     */
    const char* env_pressure_gauge =
        is_ahmed ? "AHMED_PRESSURE_GAUGE" : "KARMAN_PRESSURE_GAUGE";
    const int use_pressure_gauge =
        ahmed_env_flag(env_pressure_gauge, 0);

    if(use_pressure_gauge) {
        ensure_pressure_gauge(&sys);
    }
    else if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s pressure gauge: DISABLED by %s=0.\n",
            CODENAME, env_pressure_gauge);
    }

    report_volume_bc_baseline_diagnostics(
        &sys, fluid_case_mode_name(case_mode));

    /* mono_com is valid here; report BCs only after communication init. */
    BBFE_fluid_mixed_report_prism_bc(&sys);

    BBFE_fluid_mixed_report_prism_partition(
        &(sys.fe_pri), &(sys.mono_com));

    /* ------------------------------------------------------------ */
    /* Case-specific wall surfaces.                                 */
    /* ------------------------------------------------------------ */

    if(is_ahmed) {
        /*
         * Ahmed: dual owner families.  Assembly copies may overlap; physical
         * internal surfaces are also loaded for non-duplicated diagnostics.
         */
        BBFE_fluid_pre_surface_mixed(
            &surf_pri_assembly,
            &basis_surf_pri_assembly,
            sys.cond.directory,
            surf_pri_filename,
            sys.vals.num_ip_each_axis);

        BBFE_fluid_pre_surface_mixed(
            &surf_tet_assembly,
            &basis_surf_tet_assembly,
            sys.cond.directory,
            surf_tet_filename,
            sys.vals.num_ip_each_axis);

        BBFE_fluid_pre_surface_internal_mixed(
            &surf_pri_internal,
            &basis_surf_pri_internal,
            sys.cond.directory,
            surf_pri_filename,
            sys.vals.num_ip_each_axis);

        BBFE_fluid_pre_surface_internal_mixed(
            &surf_tet_internal,
            &basis_surf_tet_internal,
            sys.cond.directory,
            surf_tet_filename,
            sys.vals.num_ip_each_axis);

        precheck_wall_owner_map(
            "PRI6-wall assembly",
            &surf_pri_assembly,
            &(sys.fe_pri),
            BBFE_pri6_tri3_make_owner_map,
            0,
            &(sys.mono_com));

        precheck_wall_owner_map(
            "TET4-wall assembly",
            &surf_tet_assembly,
            &(sys.fe_tet),
            BBFE_tet4_tri3_make_owner_map,
            0,
            &(sys.mono_com));

        precheck_wall_owner_map(
            "PRI6-wall internal",
            &surf_pri_internal,
            &(sys.fe_pri),
            BBFE_pri6_tri3_make_owner_map,
            1,
            &(sys.mono_com));

        precheck_wall_owner_map(
            "TET4-wall internal",
            &surf_tet_internal,
            &(sys.fe_tet),
            BBFE_tet4_tri3_make_owner_map,
            1,
            &(sys.mono_com));

        /* Strong no-slip modes are true component-wise Dirichlet rows. */
        if(enable_pri6_nitsche &&
           BBFE_pri6_wall_mode_is_strong_noslip(wall_mode_pri)) {
            ahmed_add_strong_noslip_surface_bc(
                &sys,&surf_pri_assembly,"Ahmed PRI6");
        }
        if(enable_tet4_nitsche &&
           BBFE_pri6_wall_mode_is_strong_noslip(wall_mode_tet)) {
            ahmed_add_strong_noslip_surface_bc(
                &sys,&surf_tet_assembly,"Ahmed TET4");
        }

        /* Only Nitsche-mode faces participate in the strong/Nitsche overlap check. */
        check_nitsche_wall_strong_bc_junctions(
            (enable_pri6_nitsche && BBFE_pri6_wall_mode_uses_nitsche(wall_mode_pri))
                ? &surf_pri_internal : NULL,
            (enable_tet4_nitsche && BBFE_pri6_wall_mode_uses_nitsche(wall_mode_tet))
                ? &surf_tet_internal : NULL,
            &(sys.bc_NR),
            nnode,
            &(sys.mono_com));

        MixedWallPhysicalSurfaceDiagnostics physical_surface_diag;
        if(BBFE_mixed_wall_physical_surface_diagnostics_global(
            &surf_pri_internal, &(sys.fe_pri),
            &surf_tet_internal, &(sys.fe_tet),
            &(sys.mono_com),
            &physical_surface_diag) != 0) {
            fprintf(stderr,
                "%s ERROR: physical/internal Ahmed-wall diagnostics failed.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "%s Ahmed physical wall: PRI6 faces=%lld area=%e, "
                "TET4 faces=%lld area=%e, total faces=%lld area=%e\n",
                CODENAME,
                physical_surface_diag.pri_faces,
                physical_surface_diag.pri_area,
                physical_surface_diag.tet_faces,
                physical_surface_diag.tet_area,
                physical_surface_diag.total_faces,
                physical_surface_diag.total_area);
        }
    }
    else {
        /*
         * Karman: mixed TET4+PRI6 volume, PRI6-owned wall only.
         * Only the selected KARMAN_SURF_FILE is required for the wall surface.
         * Do NOT attempt to open Ahmed-specific dual surface files or
         * any *.n_internal companions.
         */
        BBFE_fluid_pre_surface_mixed(
            &surf_pri_assembly,
            &basis_surf_pri_assembly,
            sys.cond.directory,
            surf_single_filename,
            sys.vals.num_ip_each_axis);

        precheck_wall_owner_map(
            "Karman PRI6-wall assembly",
            &surf_pri_assembly,
            &(sys.fe_pri),
            BBFE_pri6_tri3_make_owner_map,
            0,
            &(sys.mono_com));

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "%s Karman surface input: %s only (PRI6-owned wall); "
                "TET4 remains active in the volume, while Ahmed dual-surface/internal diagnostics are disabled.\n",
                CODENAME, surf_single_filename);
        }
    }

    /* ------------------------------------------------------------ */
    /* Case-specific initial condition; strong BCs overwrite later. */
    /* ------------------------------------------------------------ */

    if(is_ahmed) {
        const int init_wall_zero_requested =
            ahmed_env_flag("AHMED_INIT_WALL_ZERO", 1);

        /*
         * Only Nitsche-enabled wall families may be wall-zero initialized.
         * With both families OFF, force a uniform baseline even if the user
         * leaves AHMED_INIT_WALL_ZERO=1 in an older launch command.
         */
        const int zero_pri_wall =
            enable_pri6_nitsche &&
            (BBFE_pri6_wall_mode_is_strong_noslip(wall_mode_pri) ||
             (init_wall_zero_requested && wall_mode_pri == PRI6_WALL_NITSCHE_NOSLIP));
        const int zero_tet_wall =
            enable_tet4_nitsche &&
            (BBFE_pri6_wall_mode_is_strong_noslip(wall_mode_tet) ||
             (init_wall_zero_requested && wall_mode_tet == PRI6_WALL_NITSCHE_NOSLIP));
        const double ahmed_initial_u =
            AHMED_INITIAL_U * ahmed_velocity_scale;

        if(init_wall_zero_requested && !zero_pri_wall && !zero_tet_wall &&
           monolis_mpi_get_global_my_rank() == 0) {
            fprintf(
                stderr,
                "%s WARNING: AHMED_INIT_WALL_ZERO=1 but neither enabled family uses "
                "a no-slip mode; keeping the non-no-slip wall initialization unchanged.\n",
                CODENAME);
        }

        const int local_init_wall_nodes =
            initialize_velocity_pressure_ahmed(
                sys.vals.v,
                sys.vals.p,
                nnode,
                ahmed_initial_u,
                &surf_pri_assembly,
                &surf_tet_assembly,
                zero_pri_wall,
                zero_tet_wall);

        if(local_init_wall_nodes < 0) {
            fprintf(stderr,
                "%s ERROR: Ahmed initial-condition construction failed.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }

        int global_init_wall_node_copies = local_init_wall_nodes;
        monolis_allreduce_I(
            1,
            &global_init_wall_node_copies,
            MONOLIS_MPI_SUM,
            sys.mono_com.comm);

        if(monolis_mpi_get_global_my_rank() == 0) {
            if(zero_pri_wall || zero_tet_wall) {
                printf(
                    "%s Ahmed initial condition: far field=(%.6g,0,0), "
                    "wall-zero PRI6=%s TET4=%s, wall_node_copies=%d.\n",
                    CODENAME,
                    ahmed_initial_u,
                    zero_pri_wall ? "ON" : "OFF",
                    zero_tet_wall ? "ON" : "OFF",
                    global_init_wall_node_copies);
            }
            else {
                printf(
                    "%s Ahmed initial condition: uniform=(%.6g,0,0); "
                    "wall-zero initialization DISABLED for all wall families.\n",
                    CODENAME,
                    ahmed_initial_u);
            }
        }
    }
    else {
        const double karman_initial_u =
            ahmed_env_positive_double(
                "KARMAN_INITIAL_U", KARMAN_INITIAL_U_DEFAULT, 1);

        if(initialize_velocity_pressure_uniform(
            sys.vals.v, sys.vals.p, nnode, karman_initial_u) != 0) {
            fprintf(stderr,
                "%s ERROR: Karman initial-condition construction failed.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "%s Karman initial condition: uniform=(%.6g,0,0), p=0; "
                "strong D_bc_v rows overwrite this field.\n",
                CODENAME, karman_initial_u);
        }
    }

    double* initial_vec = BB_std_calloc_1d_double(NULL, 4*nnode);

    BBFE_fluid_sups_add_velocity_pressure(
        sys.vals.v,
        sys.vals.p,
        initial_vec,
        nnode);

    for(int i=0; i<4*nnode; ++i) {
        if(sys.bc.D_bc_exists[i]) {
            initial_vec[i] = sys.bc.imposed_D_val[i];
        }
    }

    BBFE_fluid_mixed_sync_solution_vector(
        &(sys.mono_com), nnode, initial_vec);

    BBFE_fluid_sups_renew_velocity(
        sys.vals.v, initial_vec, nnode);

    BBFE_fluid_sups_renew_pressure(
        sys.vals.p, initial_vec, nnode);

    BB_std_free_1d_double(initial_vec, 4*nnode);

    BBFE_fluid_mixed_sync_state(
        &(sys.mono_com), &(sys.vals), nnode);

    ROM_BB_vec_copy_2d(
        sys.vals.v,
        sys.vals.v_old,
        nnode,
        3);

    /* ------------------------------------------------------------ */
    /* Case wall: symmetric Nitsche no-slip; Ahmed may use both owners. */
    /* ------------------------------------------------------------ */

    const double Uw_diag[3] = {0.0, 0.0, 0.0};

    PRI6WallOptions wall_pri;
    PRI6WallOptions wall_tet;
    BBFE_pri6_wall_options_default(&wall_pri);
    BBFE_pri6_wall_options_default(&wall_tet);

    wall_pri.wall_mode = wall_mode_pri;
    wall_tet.wall_mode = wall_mode_tet;
    wall_pri.exchange_mode = PRI6_WALL_EXCHANGE_OPPOSITE_FACE;
    wall_tet.exchange_mode = PRI6_WALL_EXCHANGE_OPPOSITE_FACE;

    if(is_ahmed) {
        /*
         * Backward-compatible common defaults.  Family-specific variables
         * override common values.  Legacy AHMED_NITSCHE_GAMMA/CDT remain
         * accepted for the Nitsche penalty coefficients.
         */
        const double gamma_common = ahmed_env_positive_double(
            "AHMED_NITSCHE_GAMMA", 1.0e4, 0);
        const double cdt_common = ahmed_env_positive_double(
            "AHMED_NITSCHE_CDT", 1.0, 1);

        wall_pri.gamma_n = ahmed_env_positive_double(
            "AHMED_NITSCHE_GAMMA_PRI", gamma_common, 0);
        wall_pri.dt_penalty_coeff = ahmed_env_positive_double(
            "AHMED_NITSCHE_CDT_PRI", cdt_common, 1);
        wall_tet.gamma_n = ahmed_env_positive_double(
            "AHMED_NITSCHE_GAMMA_TET", gamma_common, 0);
        wall_tet.dt_penalty_coeff = ahmed_env_positive_double(
            "AHMED_NITSCHE_CDT_TET", cdt_common, 1);

        wall_pri.exchange_mode = ahmed_wall_exchange_setting(
            "AHMED_WALL_EXCHANGE_PRI","AHMED_WALL_EXCHANGE",
            PRI6_WALL_EXCHANGE_OPPOSITE_FACE);
        wall_tet.exchange_mode = ahmed_wall_exchange_setting(
            "AHMED_WALL_EXCHANGE_TET","AHMED_WALL_EXCHANGE",
            PRI6_WALL_EXCHANGE_OPPOSITE_FACE);

        wall_pri.kappa = ahmed_env_family_double(
            "AHMED_WALL_KAPPA_PRI","AHMED_WALL_KAPPA",0.41,0);
        wall_tet.kappa = ahmed_env_family_double(
            "AHMED_WALL_KAPPA_TET","AHMED_WALL_KAPPA",0.41,0);
        wall_pri.B = ahmed_env_family_double(
            "AHMED_WALL_B_PRI","AHMED_WALL_B",5.2,0);
        wall_tet.B = ahmed_env_family_double(
            "AHMED_WALL_B_TET","AHMED_WALL_B",5.2,0);
        wall_pri.ut_eps = ahmed_env_family_double(
            "AHMED_WALL_UT_EPS_PRI","AHMED_WALL_UT_EPS",1.0e-12,0);
        wall_tet.ut_eps = ahmed_env_family_double(
            "AHMED_WALL_UT_EPS_TET","AHMED_WALL_UT_EPS",1.0e-12,0);

        const double feature_angle = ahmed_env_positive_double(
            "AHMED_WALL_FEATURE_ANGLE",30.0,0);
        if(feature_angle >= 90.0 ||
           BBFE_mixed_wall_set_feature_angle_deg(feature_angle) != 0) {
            fprintf(stderr,
                "%s ERROR: AHMED_WALL_FEATURE_ANGLE must satisfy 0 < angle < 90 deg.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }
    }
    else {
        /* Karman keeps its established single PRI6-owned wall path. */
        enable_pri6_nitsche = 1;
        enable_tet4_nitsche = 0;

        wall_pri.gamma_n = ahmed_env_positive_double(
            "KARMAN_NITSCHE_GAMMA", 1.0e4, 0);
        wall_pri.dt_penalty_coeff = ahmed_env_positive_double(
            "KARMAN_NITSCHE_CDT", 1.0, 1);

        wall_tet = wall_pri;
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        if(is_ahmed) {
            printf(
                "%s Ahmed wall PRI6: %s mode=%s exchange=%s "
                "gamma=%.6g c_dt=%.6g kappa=%.6g B=%.6g ut_eps=%.3e\n",
                CODENAME,
                enable_pri6_nitsche ? "ON" : "OFF",
                enable_pri6_nitsche ? BBFE_pri6_wall_mode_name(wall_pri.wall_mode) : "off",
                enable_pri6_nitsche ? BBFE_pri6_wall_exchange_name(wall_pri.exchange_mode) : "-",
                wall_pri.gamma_n,wall_pri.dt_penalty_coeff,
                wall_pri.kappa,wall_pri.B,wall_pri.ut_eps);
            printf(
                "%s Ahmed wall TET4: %s mode=%s exchange=%s "
                "gamma=%.6g c_dt=%.6g kappa=%.6g B=%.6g ut_eps=%.3e\n",
                CODENAME,
                enable_tet4_nitsche ? "ON" : "OFF",
                enable_tet4_nitsche ? BBFE_pri6_wall_mode_name(wall_tet.wall_mode) : "off",
                enable_tet4_nitsche ? BBFE_pri6_wall_exchange_name(wall_tet.exchange_mode) : "-",
                wall_tet.gamma_n,wall_tet.dt_penalty_coeff,
                wall_tet.kappa,wall_tet.B,wall_tet.ut_eps);
            printf(
                "%s Ahmed wall controls: AHMED_WALL_MODE[_PRI/_TET], "
                "AHMED_WALL_EXCHANGE[_PRI/_TET], "
                "AHMED_NITSCHE_GAMMA_PRI/TET, AHMED_NITSCHE_CDT_PRI/TET, "
                "AHMED_WALL_KAPPA/B/UT_EPS[_PRI/_TET], "
                "AHMED_WALL_FEATURE_ANGLE=%.6g, AHMED_VELOCITY_SCALE=%.9g.\n",
                CODENAME,BBFE_mixed_wall_get_feature_angle_deg(),ahmed_velocity_scale);
        }
        else {
            printf(
                "%s Karman Nitsche: PRI6=ON gamma=%.6g c_dt=%.6g; "
                "TET4=OFF\n",
                CODENAME,
                wall_pri.gamma_n,
                wall_pri.dt_penalty_coeff);
        }

        if(is_ahmed) {
            ahmed_write_wall_runtime_config(
                sys.cond.directory,enable_pri6_nitsche,&wall_pri,
                enable_tet4_nitsche,&wall_tet);
        }

        printf(
            "%s pressure gauge runtime switch (%s): %s\n",
            CODENAME, env_pressure_gauge,
            use_pressure_gauge ? "ON" : "OFF");
    }

    /* ------------------------------------------------------------ */
    /* Physical wall diagnostics: internal/non-duplicated surfaces only. */
    /* ------------------------------------------------------------ */

    int diagnostic_interval = 0;
    if(is_ahmed) {
        diagnostic_interval =
            ahmed_env_positive_int("AHMED_DIAGNOSTIC_INTERVAL", 10);

        if(monolis_mpi_get_global_my_rank() == 0) {
            BBFE_mixed_reset_wall_diagnostics_file(sys.cond.directory);
            ahmed_reset_wall_summary_csv(sys.cond.directory);
            printf(
                "%s Ahmed wall diagnostics: step=1 and every %d step(s); "
                "source=internal/non-duplicated surfaces.\n",
                CODENAME,
                diagnostic_interval);
            printf(
                "%s y+ definition: tau_w=|Pt mu(grad(u)+grad(u)^T)n|, "
                "u_tau=sqrt(tau_w/rho), y+=rho*u_tau*h_n/mu.\n",
                CODENAME);
        }
    }
    else if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s Ahmed-only wall y+/force diagnostics: DISABLED for FLUID_CASE=karman.\n",
            CODENAME);
    }

    double ahmed_force_Uref = 0.0;
    double ahmed_force_Aref = 0.0;
    int ahmed_internal_tet_elems = -1;
    int ahmed_internal_pri_elems = -1;
    int ahmed_continuity_enabled = 0;
    if(is_ahmed) {
        ahmed_force_Uref = ahmed_env_positive_double(
            "AHMED_U_REF", AHMED_INITIAL_U * ahmed_velocity_scale, 0);
        ahmed_force_Aref = ahmed_env_positive_double(
            "AHMED_A_REF", 0.112032, 0);

        if(monolis_mpi_get_global_my_rank() == 0) {
            ahmed_reset_force_coefficients_file(sys.cond.directory);
            printf(
                "%s Ahmed force coefficients: U_ref=%.9g A_ref=%.9g "
                "qA=%.9e; output=ahmed_force_coefficients.csv\n",
                CODENAME,
                ahmed_force_Uref,
                ahmed_force_Aref,
                0.5 * sys.vals.density * ahmed_force_Uref * ahmed_force_Uref *
                    ahmed_force_Aref);
            printf(
                "%s Ahmed coefficient axes: Cd=x(streamwise), "
                "Cy=y(side), Cl=z(vertical).\n",
                CODENAME);
        }

        if(ahmed_read_internal_element_count(
               sys.cond.directory, FLUID_MIXED_ELEM_TET_FILE,
               &ahmed_internal_tet_elems) == 0 &&
           ahmed_read_internal_element_count(
               sys.cond.directory, FLUID_MIXED_ELEM_PRI_FILE,
               &ahmed_internal_pri_elems) == 0) {
            ahmed_continuity_enabled = 1;
            if(monolis_mpi_get_global_my_rank() == 0) {
                ahmed_reset_continuity_csv(sys.cond.directory);
                printf(
                    "%s Ahmed continuity diagnostics: physical/internal "
                    "TET4+PRI6 elements; output=ahmed_continuity.csv\n",
                    CODENAME);
            }
        }
        else if(monolis_mpi_get_global_my_rank() == 0) {
            fprintf(stderr,
                "%s WARNING: graph_elem_*.dat.n_internal metadata unavailable; "
                "continuity CSV disabled to avoid overlap double counting.\n",
                CODENAME);
        }
    }

    double t = 0.0;
    int step = 0;
    int file_num = 0;

    output_files_mixed(&sys, file_num, t);
    ++file_num;

    /*
     * Surface-copy policy:
     *   - surf_*_assembly: overlap/duplicate copies are allowed and assembled.
     *   - surf_*_internal: each physical face occurs once globally and is the
     *     ONLY surface class allowed for SUM-type area/force/Cd/Cl diagnostics.
     */

    while(t < sys.vals.finish_time) {
        t += sys.vals.dt;
        ++step;

        const int solver_status =
            solver_fom_VMS_mixed_dual_nitsche(
                &sys,
                &surf_pri_assembly,
                &basis_surf_pri_assembly,
                &surf_tet_assembly,
                &basis_surf_tet_assembly,
                enable_pri6_nitsche,
                &wall_pri,
                enable_tet4_nitsche,
                &wall_tet,
                t,
                step);

        if(solver_status != 0) {
            if(monolis_mpi_get_global_my_rank() == 0) {
                fprintf(stderr,
                    "%s ERROR: %s VMS Nitsche solver failed "
                    "at step=%d t=%.15e. Time stepping is stopped.\n",
                    CODENAME, fluid_case_mode_name(case_mode), step, t);
            }
            exit(EXIT_FAILURE);
        }

        /*
         * Physical/SUM diagnostics must use *_internal, never *_assembly.
         * The state is synchronized first so wall interpolation sees current
         * owner values on overlap nodes.
         */
        if(is_ahmed && (step == 1 || step % diagnostic_interval == 0)) {
            BBFE_fluid_mixed_sync_state(
                &(sys.mono_com), &(sys.vals), nnode);

            MixedWallDiagnostics wall_diag;
            if(BBFE_mixed_wall_diagnostics_global(
                &surf_pri_internal,
                &(sys.fe_pri),
                &basis_surf_pri_internal,
                &surf_tet_internal,
                &(sys.fe_tet),
                &basis_surf_tet_internal,
                &(sys.vals),
                sys.vals.density,
                sys.vals.viscosity,
                Uw_diag,
                &(sys.mono_com),
                &wall_diag) != 0) {

                fprintf(stderr,
                    "%s ERROR: mixed Ahmed-wall y+/state diagnostics failed "
                    "at step=%d.\n",
                    CODENAME,
                    step);
                exit(EXIT_FAILURE);
            }

            int p_ref_nodes = 0;
            const double p_ref =
                ahmed_pressure_reference_global(&sys, &p_ref_nodes);

            if(monolis_mpi_get_global_my_rank() == 0) {
                printf(
                    "[Ahmed wall diag step=%d] "
                    "y+ PRI6=[%.3e mean %.3e max %.3e] "
                    "TET4=[%.3e mean %.3e max %.3e] "
                    "TOTAL=[%.3e mean %.3e max %.3e]\n",
                    step,
                    wall_diag.pri6.y_plus_min,
                    wall_diag.pri6.y_plus_mean,
                    wall_diag.pri6.y_plus_max,
                    wall_diag.tet4.y_plus_min,
                    wall_diag.tet4.y_plus_mean,
                    wall_diag.tet4.y_plus_max,
                    wall_diag.total.y_plus_min,
                    wall_diag.total.y_plus_mean,
                    wall_diag.total.y_plus_max);

                printf(
                    "[Ahmed wall diag step=%d] "
                    "max|u.n|=%e max|u_t|=%e "
                    "tau_w(mean/max)=%e/%e h_n(mean/min/max)=%e/%e/%e\n",
                    step,
                    wall_diag.total.abs_u_normal_max,
                    wall_diag.total.u_tangent_max,
                    wall_diag.total.tau_w_mean,
                    wall_diag.total.tau_w_max,
                    wall_diag.total.h_n_mean,
                    wall_diag.total.h_n_min,
                    wall_diag.total.h_n_max);

                BBFE_mixed_output_wall_diagnostics(
                    &wall_diag,
                    t,
                    sys.cond.directory);

                printf(
                    "[Ahmed wall diag step=%d] p_ref=%e reference_nodes=%d\n",
                    step, p_ref, p_ref_nodes);

                ahmed_output_wall_summary_csv(
                    &wall_diag,
                    t,
                    step,
                    p_ref,
                    sys.vals.density,
                    ahmed_force_Uref,
                    sys.cond.directory);
            }
        }

        if(is_ahmed && ahmed_continuity_enabled &&
           (step == 1 || step % diagnostic_interval == 0)) {
            AhmedContinuityDiagnostics cdiag;
            if(ahmed_continuity_global(
                &sys,
                ahmed_internal_tet_elems,
                ahmed_internal_pri_elems,
                &cdiag) != 0) {
                fprintf(stderr,
                    "%s ERROR: Ahmed continuity diagnostics failed at step=%d.\n",
                    CODENAME,step);
                exit(EXIT_FAILURE);
            }

            if(monolis_mpi_get_global_my_rank() == 0) {
                const double net_mdot = sys.vals.density*cdiag.int_div_u;
                const double mdot_ref =
                    sys.vals.density*fabs(ahmed_force_Uref)*ahmed_force_Aref;
                printf(
                    "[Ahmed continuity step=%d] mean|divu|=%e rms=%e max=%e "
                    "net_mdot=%e rel_ref=%e\n",
                    step,
                    cdiag.mean_abs_div_u,
                    cdiag.rms_div_u,
                    cdiag.max_abs_div_u,
                    net_mdot,
                    mdot_ref > 0.0 ? fabs(net_mdot)/mdot_ref : 0.0);
                ahmed_output_continuity_csv(
                    &cdiag,t,step,sys.vals.density,
                    ahmed_force_Uref,ahmed_force_Aref,
                    sys.cond.directory);
            }
        }

        if(is_ahmed && (step == 1 || step % diagnostic_interval == 0)) {
            AhmedPhysicalForce force;
            if(ahmed_calc_physical_force_global(
                &surf_pri_internal,
                &basis_surf_pri_internal,
                &surf_tet_internal,
                &basis_surf_tet_internal,
                &sys,
                enable_pri6_nitsche,
                enable_tet4_nitsche,
                &wall_pri,
                &wall_tet,
                ahmed_force_Uref,
                ahmed_force_Aref,
                &force) != 0) {

                fprintf(stderr,
                    "%s ERROR: Ahmed physical-force/Cd diagnostics failed at step=%d.\n",
                    CODENAME, step);
                exit(EXIT_FAILURE);
            }

            if(monolis_mpi_get_global_my_rank() == 0) {
                const double qA =
                    0.5 * sys.vals.density * ahmed_force_Uref *
                    ahmed_force_Uref * ahmed_force_Aref;

                printf(
                    "[Ahmed force step=%d] Cd=%+.6e (p=%+.6e visc=%+.6e) "
                    "Cy=%+.6e Cl=%+.6e NitschePenalty=(%+.3e,%+.3e,%+.3e)\n",
                    step,
                    force.total[0]/qA,
                    force.pressure[0]/qA,
                    force.viscous[0]/qA,
                    force.total[1]/qA,
                    force.total[2]/qA,
                    force.nitsche_penalty[0]/qA,
                    force.nitsche_penalty[1]/qA,
                    force.nitsche_penalty[2]/qA);

                ahmed_output_force_coefficients(
                    &force, t, step, sys.vals.density,
                    ahmed_force_Uref, ahmed_force_Aref,
                    sys.cond.directory);
            }
        }

        if(step % sys.vals.output_interval == 0) {
            output_files_mixed(&sys, file_num, t);
            ++file_num;
        }
    }

    /*
     * Force/Cd/Cl rule:
     *   Any current or future physical force routine MUST receive
     *   surf_pri_internal / surf_tet_internal (and their matching bases), not
     *   the assembly surfaces.  The internal surfaces have already passed the
     *   mesh-independent owner-map checks above.
     *
     * The actual combined PRI6+TET4 Ahmed drag/lift formula is intentionally
     * kept separate from the assembly solver.  Existing PRI6 force diagnostics
     * and TET4 force routines can be called here, but only with *_internal.
     */

    BBFE_fluid_finalize_surface_mixed(
        &surf_pri_assembly, &basis_surf_pri_assembly);

    BBFE_fluid_finalize_surface_mixed(
        &surf_tet_assembly, &basis_surf_tet_assembly);

    BBFE_fluid_finalize_surface_mixed(
        &surf_pri_internal, &basis_surf_pri_internal);

    BBFE_fluid_finalize_surface_mixed(
        &surf_tet_internal, &basis_surf_tet_internal);

    BBFE_fluid_finalize_mixed(
        &(sys.fe_tet), &(sys.basis_tet),
        &(sys.fe_pri), &(sys.basis_pri));

    BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), nnode, 4);
    BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc_NR), nnode, 4);

    monolis_finalize(&(sys.monolis));
    monolis_global_finalize();

    return 0;
}
