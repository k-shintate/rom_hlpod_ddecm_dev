#include "core_Nitsche_prism.h"
#include "fluid_linalg_compat.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/* Small vector / matrix utilities                                             */
/* -------------------------------------------------------------------------- */

static double pri6_dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double pri6_norm3(const double a[3])
{
    return sqrt(pri6_dot3(a, a));
}

static void pri6_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static double pri6_det3(const double A[3][3])
{
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
      - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
      + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int pri6_inv3(const double A[3][3], double invA[3][3])
{
    const double det = pri6_det3(A);
    if(!isfinite(det) || fabs(det) < 1.0e-300) return 0;

    const double id = 1.0 / det;

    invA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) * id;
    invA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) * id;
    invA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * id;

    invA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) * id;
    invA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * id;
    invA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) * id;

    invA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) * id;
    invA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) * id;
    invA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * id;

    return 1;
}

/* -------------------------------------------------------------------------- */
/* Wall options                                                               */
/* -------------------------------------------------------------------------- */

void BBFE_pri6_wall_options_default(PRI6WallOptions* opt)
{
    if(opt == NULL) return;

    opt->wall_mode        = PRI6_WALL_NITSCHE_NOSLIP;
    opt->exchange_mode    = PRI6_WALL_EXCHANGE_OPPOSITE_FACE;
    opt->gamma_n          = 100.0;
    opt->dt_penalty_coeff = 0.10;
    opt->kappa            = 0.41;
    opt->B                = 5.2;
    opt->ut_eps           = 1.0e-12;
}

static double g_mixed_wall_feature_angle_deg = 30.0;

const char* BBFE_pri6_wall_mode_name(int mode)
{
    switch(mode) {
    case PRI6_WALL_STRONG_NOSLIP:          return "strong-noslip";
    case PRI6_WALL_NITSCHE_NOSLIP:         return "nitsche-noslip";
    case PRI6_WALL_NITSCHE_FREESLIP:       return "nitsche-freeslip";
    case PRI6_WALL_NITSCHE_SPALDING:       return "nitsche-spalding";
    case PRI6_WALL_STRONG_NORMAL_FREESLIP: return "strong-normal-freeslip";
    case PRI6_WALL_STRONG_NORMAL_SPALDING: return "strong-normal-spalding";
    default:                                return "unknown";
    }
}

int BBFE_pri6_wall_mode_from_string(const char* text, int* mode_out)
{
    if(text == NULL || mode_out == NULL) return -1;
    if(strcmp(text,"strong-noslip") == 0) {
        *mode_out = PRI6_WALL_STRONG_NOSLIP; return 0;
    }
    if(strcmp(text,"nitsche-noslip") == 0) {
        *mode_out = PRI6_WALL_NITSCHE_NOSLIP; return 0;
    }
    if(strcmp(text,"nitsche-freeslip") == 0 || strcmp(text,"free-slip") == 0) {
        *mode_out = PRI6_WALL_NITSCHE_FREESLIP; return 0;
    }
    if(strcmp(text,"nitsche-spalding") == 0 || strcmp(text,"spalding") == 0) {
        *mode_out = PRI6_WALL_NITSCHE_SPALDING; return 0;
    }
    if(strcmp(text,"strong-normal-freeslip") == 0 ||
       strcmp(text,"projected-strong-normal-freeslip") == 0) {
        *mode_out = PRI6_WALL_STRONG_NORMAL_FREESLIP; return 0;
    }
    if(strcmp(text,"strong-normal-spalding") == 0 ||
       strcmp(text,"projected-strong-normal-spalding") == 0) {
        *mode_out = PRI6_WALL_STRONG_NORMAL_SPALDING; return 0;
    }
    return -1;
}

const char* BBFE_pri6_wall_exchange_name(int mode)
{
    switch(mode) {
    case PRI6_WALL_EXCHANGE_WALL_Q:        return "wall-q";
    case PRI6_WALL_EXCHANGE_OPPOSITE_FACE: return "opposite-face";
    default:                                return "unknown";
    }
}

int BBFE_pri6_wall_exchange_from_string(const char* text, int* mode_out)
{
    if(text == NULL || mode_out == NULL) return -1;
    if(strcmp(text,"wall-q") == 0 || strcmp(text,"wall") == 0) {
        *mode_out = PRI6_WALL_EXCHANGE_WALL_Q; return 0;
    }
    if(strcmp(text,"opposite-face") == 0 || strcmp(text,"opposite") == 0 ||
       strcmp(text,"opposite-node") == 0) {
        *mode_out = PRI6_WALL_EXCHANGE_OPPOSITE_FACE; return 0;
    }
    return -1;
}

int BBFE_pri6_wall_mode_uses_nitsche(int mode)
{
    return mode == PRI6_WALL_NITSCHE_NOSLIP ||
           mode == PRI6_WALL_NITSCHE_FREESLIP ||
           mode == PRI6_WALL_NITSCHE_SPALDING;
}

int BBFE_pri6_wall_mode_uses_spalding(int mode)
{
    return mode == PRI6_WALL_NITSCHE_SPALDING ||
           mode == PRI6_WALL_STRONG_NORMAL_SPALDING;
}

int BBFE_pri6_wall_mode_uses_projected_strong_normal(int mode)
{
    return mode == PRI6_WALL_STRONG_NORMAL_FREESLIP ||
           mode == PRI6_WALL_STRONG_NORMAL_SPALDING;
}

int BBFE_pri6_wall_mode_is_strong_noslip(int mode)
{
    return mode == PRI6_WALL_STRONG_NOSLIP;
}

int BBFE_mixed_wall_set_feature_angle_deg(double feature_angle_deg)
{
    if(!isfinite(feature_angle_deg) || feature_angle_deg <= 0.0 ||
       feature_angle_deg >= 90.0) return -1;
    g_mixed_wall_feature_angle_deg = feature_angle_deg;
    return 0;
}

double BBFE_mixed_wall_get_feature_angle_deg(void)
{
    return g_mixed_wall_feature_angle_deg;
}

static int pri6_validate_wall_options(const PRI6WallOptions* opt)
{
    if(opt == NULL) return 0;

    if(opt->wall_mode < PRI6_WALL_STRONG_NOSLIP ||
       opt->wall_mode > PRI6_WALL_STRONG_NORMAL_SPALDING) return 0;

    if(opt->exchange_mode != PRI6_WALL_EXCHANGE_WALL_Q &&
       opt->exchange_mode != PRI6_WALL_EXCHANGE_OPPOSITE_FACE) return 0;

    if(!isfinite(opt->gamma_n) || opt->gamma_n <= 0.0) return 0;
    if(!isfinite(opt->dt_penalty_coeff) || opt->dt_penalty_coeff < 0.0) return 0;
    if(!isfinite(opt->kappa) || opt->kappa <= 0.0) return 0;
    if(!isfinite(opt->B)) return 0;
    if(!isfinite(opt->ut_eps) || opt->ut_eps <= 0.0) return 0;

    return 1;
}

/* -------------------------------------------------------------------------- */
/* TRI3 -> PRI6 owner map                                                     */
/* -------------------------------------------------------------------------- */

typedef long long pri6_node_id_t;

static const int PRI6_TRI_FACE[2][3] = {
    {0, 1, 2},
    {3, 4, 5}
};

typedef struct {
    pri6_node_id_t key[3];
    int owner_elem;
    int owner_lface;
} PRI6FaceEntry;

typedef struct {
    pri6_node_id_t key[3];
    int surf_face_id;
} PRI6SurfEntry;

static void pri6_sort3(pri6_node_id_t k[3])
{
    if(k[1] < k[0]) { pri6_node_id_t t=k[0]; k[0]=k[1]; k[1]=t; }
    if(k[2] < k[1]) { pri6_node_id_t t=k[1]; k[1]=k[2]; k[2]=t; }
    if(k[1] < k[0]) { pri6_node_id_t t=k[0]; k[0]=k[1]; k[1]=t; }
}

static int pri6_cmp3(const pri6_node_id_t a[3], const pri6_node_id_t b[3])
{
    if(a[0] != b[0]) return a[0] < b[0] ? -1 : +1;
    if(a[1] != b[1]) return a[1] < b[1] ? -1 : +1;
    if(a[2] != b[2]) return a[2] < b[2] ? -1 : +1;
    return 0;
}

static int pri6_cmp_face_entry(const void* A, const void* B)
{
    const PRI6FaceEntry* a = (const PRI6FaceEntry*)A;
    const PRI6FaceEntry* b = (const PRI6FaceEntry*)B;
    const int c = pri6_cmp3(a->key, b->key);
    if(c) return c;
    if(a->owner_elem != b->owner_elem) return a->owner_elem < b->owner_elem ? -1 : +1;
    if(a->owner_lface != b->owner_lface) return a->owner_lface < b->owner_lface ? -1 : +1;
    return 0;
}

static int pri6_cmp_surf_entry(const void* A, const void* B)
{
    const PRI6SurfEntry* a = (const PRI6SurfEntry*)A;
    const PRI6SurfEntry* b = (const PRI6SurfEntry*)B;
    const int c = pri6_cmp3(a->key, b->key);
    if(c) return c;
    if(a->surf_face_id != b->surf_face_id) return a->surf_face_id < b->surf_face_id ? -1 : +1;
    return 0;
}

int BBFE_pri6_tri3_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    int** owner_elem_out,
    int** owner_lface_out,
    int* num_unmapped_local)
{
    if(owner_elem_out == NULL || owner_lface_out == NULL) return -1;

    *owner_elem_out = NULL;
    *owner_lface_out = NULL;
    if(num_unmapped_local) *num_unmapped_local = 0;

    if(surf == NULL || fe_pri == NULL) return -1;

    const int nsurf = surf->total_num_elems;
    if(nsurf < 0) return -1;

    /* Empty surface copy is valid on MPI ranks with no wall/prism cells. */
    if(nsurf == 0) return 0;

    if(surf->local_num_nodes != 3) {
        fprintf(stderr,
            "%s ERROR: active wall surface must be TRI3; got NEN=%d\n",
            CODENAME, surf->local_num_nodes);
        return -1;
    }

    /* A rank can theoretically have a wall copy but no local owner prism.
       Mark those faces unmapped instead of crashing; the global precheck will
       catch such a partitioning inconsistency. */
    if(!BBFE_fluid_mixed_has_elements(fe_pri)) {
        int* owner = (int*)malloc(sizeof(int)*(size_t)nsurf);
        int* lface = (int*)malloc(sizeof(int)*(size_t)nsurf);
        if(owner == NULL || lface == NULL) {
            free(owner); free(lface);
            return -1;
        }
        for(int es=0; es<nsurf; ++es) {
            owner[es] = -1;
            lface[es] = -1;
        }
        *owner_elem_out = owner;
        *owner_lface_out = lface;
        if(num_unmapped_local) *num_unmapped_local = nsurf;
        return 0;
    }

    if(fe_pri->local_num_nodes != 6) {
        fprintf(stderr,
            "%s ERROR: active prism family must be PRI6; got NEN=%d\n",
            CODENAME, fe_pri->local_num_nodes);
        return -1;
    }

    const int nvolface = 2 * fe_pri->total_num_elems;

    PRI6FaceEntry* vol = (PRI6FaceEntry*)malloc(sizeof(PRI6FaceEntry) * (size_t)nvolface);
    PRI6SurfEntry* srf = (PRI6SurfEntry*)malloc(sizeof(PRI6SurfEntry) * (size_t)nsurf);
    int* owner = (int*)malloc(sizeof(int) * (size_t)nsurf);
    int* lface = (int*)malloc(sizeof(int) * (size_t)nsurf);

    if(vol == NULL || srf == NULL || owner == NULL || lface == NULL) {
        free(vol); free(srf); free(owner); free(lface);
        return -1;
    }

    int k = 0;
    for(int e=0; e<fe_pri->total_num_elems; ++e) {
        for(int f=0; f<2; ++f) {
            pri6_node_id_t key[3] = {
                (pri6_node_id_t)fe_pri->conn[e][PRI6_TRI_FACE[f][0]],
                (pri6_node_id_t)fe_pri->conn[e][PRI6_TRI_FACE[f][1]],
                (pri6_node_id_t)fe_pri->conn[e][PRI6_TRI_FACE[f][2]]
            };
            pri6_sort3(key);
            for(int a=0; a<3; ++a) vol[k].key[a] = key[a];
            vol[k].owner_elem = e;
            vol[k].owner_lface = f;
            ++k;
        }
    }

    for(int es=0; es<nsurf; ++es) {
        pri6_node_id_t key[3] = {
            (pri6_node_id_t)surf->conn[es][0],
            (pri6_node_id_t)surf->conn[es][1],
            (pri6_node_id_t)surf->conn[es][2]
        };
        pri6_sort3(key);
        for(int a=0; a<3; ++a) srf[es].key[a] = key[a];
        srf[es].surf_face_id = es;
        owner[es] = -1;
        lface[es] = -1;
    }

    qsort(vol, (size_t)nvolface, sizeof(PRI6FaceEntry), pri6_cmp_face_entry);
    qsort(srf, (size_t)nsurf, sizeof(PRI6SurfEntry), pri6_cmp_surf_entry);

    int i=0, j=0;
    while(i<nvolface && j<nsurf) {
        const int c = pri6_cmp3(vol[i].key, srf[j].key);
        if(c < 0) {
            pri6_node_id_t key0[3] = {vol[i].key[0],vol[i].key[1],vol[i].key[2]};
            do { ++i; } while(i<nvolface && pri6_cmp3(key0,vol[i].key)==0);
        }
        else if(c > 0) {
            pri6_node_id_t key0[3] = {srf[j].key[0],srf[j].key[1],srf[j].key[2]};
            do { ++j; } while(j<nsurf && pri6_cmp3(key0,srf[j].key)==0);
        }
        else {
            pri6_node_id_t key0[3] = {srf[j].key[0],srf[j].key[1],srf[j].key[2]};
            const int oe = vol[i].owner_elem;
            const int of = vol[i].owner_lface;

            /* A duplicated surface face can occur in overlap data. Map every
               copy to the same local owner when it exists on this rank. */
            do {
                const int es = srf[j].surf_face_id;
                owner[es] = oe;
                lface[es] = of;
                ++j;
            } while(j<nsurf && pri6_cmp3(key0,srf[j].key)==0);

            do { ++i; } while(i<nvolface && pri6_cmp3(key0,vol[i].key)==0);
        }
    }

    int unmapped = 0;
    for(int es=0; es<nsurf; ++es) {
        if(owner[es] < 0) ++unmapped;
    }

    free(vol);
    free(srf);

    *owner_elem_out = owner;
    *owner_lface_out = lface;
    if(num_unmapped_local) *num_unmapped_local = unmapped;

    return 0;
}

/* -------------------------------------------------------------------------- */
/* Surface / prism geometry                                                   */
/* -------------------------------------------------------------------------- */

static int pri6_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for(int a=0; a<3; ++a) {
        s2v[a] = -1;
        const int gid = surf->conn[es][a];
        for(int b=0; b<6; ++b) {
            if(fe->conn[ke][b] == gid) {
                s2v[a] = b;
                break;
            }
        }
        if(s2v[a] < 0) return 0;
    }
    return 1;
}

static void pri6_tri3_normal(
    const double xtri[3][3],
    const double* dN_dr,
    const double* dN_ds,
    double normal[3])
{
    double xr[3] = {0.0,0.0,0.0};
    double xs[3] = {0.0,0.0,0.0};

    for(int a=0; a<3; ++a) {
        for(int d=0; d<3; ++d) {
            xr[d] += dN_dr[a] * xtri[a][d];
            xs[d] += dN_ds[a] * xtri[a][d];
        }
    }
    pri6_cross3(xr, xs, normal);
}

static double pri6_tri_area_factor(const BBFE_BASIS* basis_surf)
{
    double sum = 0.0;
    for(int p=0; p<basis_surf->num_integ_points; ++p) sum += basis_surf->integ_weight[p];

    if(fabs(sum - 0.5) < 1.0e-10) return 1.0;
    if(fabs(sum - 1.0) < 1.0e-10) return 0.5;

    /* Normalize arbitrary TRI3 rule to reference-triangle area 1/2. */
    return (fabs(sum) > 1.0e-30) ? 0.5/sum : 0.0;
}

static double pri6_tet_volume4(
    const double a[3], const double b[3], const double c[3], const double d[3])
{
    double ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    double ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    double ad[3] = {d[0]-a[0], d[1]-a[1], d[2]-a[2]};
    double cr[3];
    pri6_cross3(ac, ad, cr);
    return fabs(pri6_dot3(ab, cr)) / 6.0;
}

static double pri6_physical_volume(const double x[6][3])
{
    /* Standard wedge decomposition: (0,1,2,3), (1,2,4,3), (2,4,5,3). */
    return
        pri6_tet_volume4(x[0],x[1],x[2],x[3]) +
        pri6_tet_volume4(x[1],x[2],x[4],x[3]) +
        pri6_tet_volume4(x[2],x[4],x[5],x[3]);
}

/* Build volume N and physical grad(N) at a TRI3 quadrature point. */
static int pri6_eval_on_tri_face(
    const double xvol[6][3],
    int lface,
    const double Nv_from_surface[6],
    double Nv[6],
    double gradN[6][3],
    double J_inv[3][3],
    double* detJ_abs)
{
    if(lface != 0 && lface != 1) return 0;

    double xi[3];
    if(lface == 0) {
        xi[0] = Nv_from_surface[1];
        xi[1] = Nv_from_surface[2];
        xi[2] = -1.0;
    }
    else {
        xi[0] = Nv_from_surface[4];
        xi[1] = Nv_from_surface[5];
        xi[2] = +1.0;
    }

    BBFE_fluid_mixed_shapefunc_prism1st_get_val(xi, Nv);

    double dr[6], ds[6], dz[6];
    BBFE_fluid_mixed_shapefunc_prism1st_get_derivative(xi, dr, ds, dz);

    double J[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    for(int a=0; a<6; ++a) {
        const double dref[3] = {dr[a],ds[a],dz[a]};
        for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
                J[i][j] += xvol[a][i] * dref[j];
            }
        }
    }

    const double detJ = pri6_det3(J);
    if(!isfinite(detJ) || fabs(detJ) < 1.0e-300) return 0;
    if(!pri6_inv3(J, J_inv)) return 0;

    for(int a=0; a<6; ++a) {
        const double dref[3] = {dr[a],ds[a],dz[a]};
        for(int i=0; i<3; ++i) {
            gradN[a][i] = 0.0;
            for(int j=0; j<3; ++j) gradN[a][i] += dref[j] * J_inv[j][i];
        }
    }

    *detJ_abs = fabs(detJ);
    return 1;
}

static void pri6_opposite_face_weights(int lface, const double Nv[6], double Nopp[6])
{
    for(int a=0; a<6; ++a) Nopp[a] = 0.0;

    if(lface == 0) {
        Nopp[3] = Nv[0];
        Nopp[4] = Nv[1];
        Nopp[5] = Nv[2];
    }
    else {
        Nopp[0] = Nv[3];
        Nopp[1] = Nv[4];
        Nopp[2] = Nv[5];
    }
}

/* -------------------------------------------------------------------------- */
/* Nitsche / wall-law local operators                                         */
/* -------------------------------------------------------------------------- */

static void pri6_get_penalty_components(
    const PRI6WallOptions* opt,
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double* beta,
    double* beta_mu,
    double* beta_dt)
{
    const double h = fmax(h_n, 1.0e-12);
    const double dte = fmax(dt, 1.0e-30);
    const double b_mu = opt->gamma_n * mu_eff / h;
    const double b_dt =
        opt->gamma_n * opt->dt_penalty_coeff * rho * h / dte;

    if(beta != NULL) *beta = b_mu + b_dt;
    if(beta_mu != NULL) *beta_mu = b_mu;
    if(beta_dt != NULL) *beta_dt = b_dt;
}

static void pri6_get_penalty(
    const PRI6WallOptions* opt,
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double* beta)
{
    pri6_get_penalty_components(
        opt, rho, mu_eff, h_n, dt, beta, NULL, NULL);
}

static double pri6_spalding_residual(
    double uplus,
    double Rey,
    double kappa,
    double B)
{
    const double A = exp(-kappa*B);
    const double x = kappa*uplus;
    if(x > 50.0) return DBL_MAX*0.25;

    const double yplus = uplus + A*(exp(x)-1.0-x-0.5*x*x-x*x*x/6.0);
    return uplus*yplus - Rey;
}

static double pri6_compute_utau(
    const PRI6WallOptions* opt,
    double ut_norm,
    double y,
    double rho,
    double mu)
{
    const double eps=1.0e-30;
    if(ut_norm<=eps || y<=eps || rho<=eps || mu<=eps) return 0.0;
    const double nu=mu/rho;
    const double Rey=ut_norm*y/nu;
    if(Rey<=eps) return 0.0;

    double lo=1.0e-12, hi=1.0;
    while(pri6_spalding_residual(hi,Rey,opt->kappa,opt->B)<0.0 && hi<1.0e6) hi*=2.0;

    for(int it=0; it<80; ++it) {
        const double mid=0.5*(lo+hi);
        if(pri6_spalding_residual(mid,Rey,opt->kappa,opt->B)>0.0) hi=mid;
        else lo=mid;
    }

    const double uplus=0.5*(lo+hi);
    return uplus>eps ? ut_norm/uplus : 0.0;
}

static void pri6_spalding_coefficients(
    const PRI6WallOptions* opt,
    const double u_exchange[3],
    const double Uw[3],
    const double n[3],
    double rho,
    double mu,
    double y,
    double ut[3],
    double* beta_wall)
{
    const double eps=1.0e-30;
    const double r[3] = {
        u_exchange[0]-Uw[0],
        u_exchange[1]-Uw[1],
        u_exchange[2]-Uw[2]
    };
    const double rn=pri6_dot3(r,n);
    for(int a=0;a<3;++a) ut[a]=r[a]-rn*n[a];

    const double ut_norm=pri6_norm3(ut);
    const double y_eff=fmax(y,1.0e-12);
    double bw=fmax(mu,eps)/y_eff;

    if(ut_norm > opt->ut_eps) {
        const double utau=pri6_compute_utau(opt,ut_norm,y_eff,fmax(rho,eps),fmax(mu,eps));
        if(isfinite(utau) && utau>0.0) bw=fmax(rho,eps)*utau*utau/fmax(ut_norm,opt->ut_eps);
    }

    if(!isfinite(bw) || bw<0.0) bw=fmax(mu,eps)/y_eff;
    *beta_wall=bw;
}

static void pri6_fullvec_nitsche_local(
    double vec[4], double mat[4][4],
    double Ni, double Nj,
    const double gradNi[3], const double gradNj[3],
    const double n[3], const double u[3], double p,
    const double grad_u[3][3], const double Uw[3],
    double rho, double mu, double h_n, double dt,
    const PRI6WallOptions* opt)
{
    double beta=0.0;
    pri6_get_penalty(opt,rho,mu,h_n,dt,&beta);

    double r[3]={u[0]-Uw[0],u[1]-Uw[1],u[2]-Uw[2]};
    const double rn=pri6_dot3(r,n);

    /*
     * Symmetric Newtonian Cauchy stress:
     *
     *   sigma = -p I + mu (grad(u) + grad(u)^T)
     *   traction = sigma n
     */
    double gun[3]  ={0.0,0.0,0.0};  /* grad(u)   n */
    double gTun[3] ={0.0,0.0,0.0};  /* grad(u)^T n */

    for(int a=0;a<3;++a) {
        for(int b=0;b<3;++b) {
            gun[a]  += grad_u[a][b]*n[b];
            gTun[a] += grad_u[b][a]*n[b];
        }
    }

    double traction[3];
    for(int a=0;a<3;++a) {
        traction[a] = -p*n[a] + mu*(gun[a] + gTun[a]);
    }

    const double gni=pri6_dot3(gradNi,n);
    const double gnj=pri6_dot3(gradNj,n);
    const double gradNi_dot_r=pri6_dot3(gradNi,r);

    for(int a=0;a<3;++a) {
        /* consistency */
        vec[a] += dt*(-Ni*traction[a]);

        /*
         * adjoint consistency:
         * -(sigma(v,0)n) . (u-Uw)
         */
        vec[a] += dt*(-mu*(gni*r[a] + n[a]*gradNi_dot_r));

        /* penalty */
        vec[a] += dt*(Ni*beta*r[a]);
    }

    /* pressure-test part of adjoint consistency */
    vec[3] += dt*Ni*rn;

    for(int a=0;a<3;++a) {
        for(int b=0;b<3;++b) {
            const double dab=(a==b)?1.0:0.0;

            /*
             * d[(sigma(u,p)n)_a]/d U_{j,b}
             * = mu [delta_ab (gradNj.n) + gradNj[a] n[b]]
             */
            const double dtr =
                mu*(dab*gnj + gradNj[a]*n[b]);

            mat[a][b] += dt*(-Ni*dtr);

            /*
             * derivative of adjoint-consistency term
             */
            mat[a][b] += dt*(-mu*(
                gni*Nj*dab
                + n[a]*gradNi[b]*Nj));

            mat[a][b] += dt*(Ni*Nj*beta*dab);
        }

        /* derivative of consistency traction wrt pressure */
        mat[a][3] += dt*(Ni*Nj*n[a]);
    }

    for(int b=0;b<3;++b) {
        mat[3][b] += dt*(Ni*Nj*n[b]);
    }
}

static void pri6_normal_nitsche_local(
    double vec[4], double mat[4][4],
    double Ni, double Nj,
    const double gradNi[3], const double gradNj[3],
    const double n[3], const double u[3], double p,
    const double grad_u[3][3], const double Uw[3],
    double rho, double mu, double h_n, double dt,
    const PRI6WallOptions* opt)
{
    double beta=0.0;
    pri6_get_penalty(opt,rho,mu,h_n,dt,&beta);

    const double r[3]={u[0]-Uw[0],u[1]-Uw[1],u[2]-Uw[2]};
    const double rn=pri6_dot3(r,n);

    /*
     * Normal component of symmetric Cauchy traction:
     *
     *   n.sigma.n = -p + 2 mu n.grad(u).n
     */
    double gun[3]={0.0,0.0,0.0};
    for(int a=0;a<3;++a) {
        for(int b=0;b<3;++b) {
            gun[a]+=grad_u[a][b]*n[b];
        }
    }

    const double traction_n =
        -p + 2.0*mu*pri6_dot3(gun,n);

    const double gni=pri6_dot3(gradNi,n);
    const double gnj=pri6_dot3(gradNj,n);

    for(int a=0;a<3;++a) {
        const double wn = Ni*n[a];

        /*
         * n.sigma(v,0).n for v = Ni e_a:
         *   2 mu n[a] (gradNi.n)
         */
        const double twn =
            2.0*mu*n[a]*gni;

        vec[a] += dt*(
            -wn*traction_n
            -twn*rn
            +beta*wn*rn);
    }

    vec[3] += dt*Ni*rn;

    for(int a=0;a<3;++a) {
        for(int b=0;b<3;++b) {
            const double dtn =
                2.0*mu*n[b]*gnj;

            mat[a][b] += dt*(-Ni*n[a]*dtn);
            mat[a][b] += dt*(-2.0*mu*n[a]*gni*Nj*n[b]);
            mat[a][b] += dt*(beta*Ni*n[a]*Nj*n[b]);
        }

        mat[a][3] += dt*(Ni*Nj*n[a]);
    }

    for(int b=0;b<3;++b) {
        mat[3][b] += dt*(Ni*Nj*n[b]);
    }
}

static void pri6_add_spalding_local(
    double vec[4], double mat[4][4],
    double Ni, double exchange_Nj,
    const double n[3], const double u_exchange[3], const double Uw[3],
    double rho, double mu, double y, double dt,
    const PRI6WallOptions* opt)
{
    double ut[3], beta_wall=0.0;
    pri6_spalding_coefficients(opt,u_exchange,Uw,n,rho,mu,y,ut,&beta_wall);

    for(int a=0;a<3;++a) vec[a] += dt*Ni*beta_wall*ut[a];
    for(int a=0;a<3;++a) {
        for(int b=0;b<3;++b) {
            const double Pt=(a==b?1.0:0.0)-n[a]*n[b];
            mat[a][b] += dt*Ni*beta_wall*Pt*exchange_Nj;
        }
    }
}

static void pri6_wall_local_vecmat(
    double vec[4], double mat[4][4],
    double Ni, double Nj, double exchange_Nj,
    const double gradNi[3], const double gradNj[3],
    const double n[3], const double u[3], const double u_exchange[3],
    double p, const double grad_u[3][3], const double Uw[3],
    double rho, double mu_eff, double mu_wall,
    double h_n, double y_exchange, double dt,
    const PRI6WallOptions* opt)
{
    for(int a=0;a<4;++a) {
        vec[a]=0.0;
        for(int b=0;b<4;++b) mat[a][b]=0.0;
    }

    switch(opt->wall_mode) {
    case PRI6_WALL_STRONG_NOSLIP:
        return;
    case PRI6_WALL_NITSCHE_NOSLIP:
        pri6_fullvec_nitsche_local(vec,mat,Ni,Nj,gradNi,gradNj,n,u,p,grad_u,Uw,
                                  rho,mu_eff,h_n,dt,opt);
        return;
    case PRI6_WALL_NITSCHE_FREESLIP:
        pri6_normal_nitsche_local(vec,mat,Ni,Nj,gradNi,gradNj,n,u,p,grad_u,Uw,
                                 rho,mu_eff,h_n,dt,opt);
        return;
    case PRI6_WALL_NITSCHE_SPALDING:
        pri6_normal_nitsche_local(vec,mat,Ni,Nj,gradNi,gradNj,n,u,p,grad_u,Uw,
                                 rho,mu_eff,h_n,dt,opt);
        pri6_add_spalding_local(vec,mat,Ni,exchange_Nj,n,u_exchange,Uw,
                               rho,mu_wall,y_exchange,dt,opt);
        return;
    case PRI6_WALL_STRONG_NORMAL_FREESLIP:
        return;
    case PRI6_WALL_STRONG_NORMAL_SPALDING:
        pri6_add_spalding_local(vec,mat,Ni,exchange_Nj,n,u_exchange,Uw,
                               rho,mu_wall,y_exchange,dt,opt);
        return;
    default:
        return;
    }
}

/* -------------------------------------------------------------------------- */
/* PRI6/TRI3 Nitsche assembly                                                 */
/* -------------------------------------------------------------------------- */

int set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_pri,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    PRI6NitscheDiagnostics* diag)
{
    if(diag) {
        memset(diag,0,sizeof(*diag));
        diag->h_n_min=HUGE_VAL;
        diag->beta_min=HUGE_VAL;
        diag->beta_mu_min=HUGE_VAL;
        diag->beta_dt_min=HUGE_VAL;
    }

    if(monolis==NULL || fe_pri==NULL || basis_pri==NULL || surf==NULL ||
       basis_surf==NULL || vals==NULL || Uw==NULL || !pri6_validate_wall_options(opt)) {
        fprintf(stderr,"%s ERROR: invalid input to PRI6/TRI3 Nitsche assembly.\n",CODENAME);
        return -1;
    }

    if(surf->total_num_elems==0) return 0;
    if(fe_pri->local_num_nodes!=6 || surf->local_num_nodes!=3) return -1;

    int *owner=NULL,*lface=NULL,unmapped=0;
    if(BBFE_pri6_tri3_make_owner_map(surf,fe_pri,&owner,&lface,&unmapped)!=0) {
        free(owner); free(lface); return -1;
    }

    const int np=basis_surf->num_integ_points;
    const double area_factor=pri6_tri_area_factor(basis_surf);
    if(np<=0 || area_factor<=0.0) { free(owner); free(lface); return -1; }

    if(diag) {
        diag->surface_faces=surf->total_num_elems;
        diag->skipped_no_local_owner=unmapped;
    }

    for(int es=0; es<surf->total_num_elems; ++es) {
        const int ke=owner[es];
        const int lf=lface[es];
        if(ke<0) continue;
        if(ke>=fe_pri->total_num_elems || (lf!=0 && lf!=1)) {
            free(owner); free(lface); return -1;
        }

        int s2v[3];
        if(!pri6_surface_to_volume_nodes(surf,fe_pri,es,ke,s2v)) {
            free(owner); free(lface); return -1;
        }

        for(int a=0;a<3;++a) {
            int onface=0;
            for(int q=0;q<3;++q) if(s2v[a]==PRI6_TRI_FACE[lf][q]) onface=1;
            if(!onface) { free(owner); free(lface); return -1; }
        }

        double xvol[6][3],xsurf[3][3];
        for(int a=0;a<6;++a) {
            const int gid=fe_pri->conn[ke][a];
            for(int d=0;d<3;++d) xvol[a][d]=fe_pri->x[gid][d];
        }
        for(int a=0;a<3;++a) {
            const int gid=surf->conn[es][a];
            for(int d=0;d<3;++d) xsurf[a][d]=fe_pri->x[gid][d];
        }

        double xc[3]={0.0,0.0,0.0},xf[3]={0.0,0.0,0.0};
        for(int a=0;a<6;++a) for(int d=0;d<3;++d) xc[d]+=xvol[a][d]/6.0;
        for(int a=0;a<3;++a) for(int d=0;d<3;++d) xf[d]+=xsurf[a][d]/3.0;
        const double outdir[3]={xf[0]-xc[0],xf[1]-xc[1],xf[2]-xc[2]};

        double e1[3]={xsurf[1][0]-xsurf[0][0],xsurf[1][1]-xsurf[0][1],xsurf[1][2]-xsurf[0][2]};
        double e2[3]={xsurf[2][0]-xsurf[0][0],xsurf[2][1]-xsurf[0][1],xsurf[2][2]-xsurf[0][2]};
        double cr[3]; pri6_cross3(e1,e2,cr);
        const double face_area=0.5*pri6_norm3(cr);
        const double volume=pri6_physical_volume(xvol);
        const double h_n=fmax(volume/fmax(face_area,1.0e-30),1.0e-12);
        const double h_e_vms=cbrt(fmax(volume,1.0e-36));

        if(!isfinite(face_area)||face_area<=0.0||!isfinite(volume)||volume<=0.0) {
            free(owner); free(lface); return -1;
        }

        int qp_count=0;
        for(int p=0;p<np;++p) {
            double nraw[3];
            pri6_tri3_normal(xsurf,basis_surf->dN_dxi[p],basis_surf->dN_det[p],nraw);
            const double Jface=pri6_norm3(nraw);
            if(!isfinite(Jface)||Jface<=1.0e-300) { free(owner);free(lface);return -1; }

            double n[3]={nraw[0]/Jface,nraw[1]/Jface,nraw[2]/Jface};
            if(pri6_dot3(n,outdir)<0.0) { n[0]*=-1.0;n[1]*=-1.0;n[2]*=-1.0; }

            double Nv_surface[6]={0.0,0.0,0.0,0.0,0.0,0.0};
            for(int a=0;a<3;++a) Nv_surface[s2v[a]]=basis_surf->N[p][a];

            double Nv[6],gradN[6][3],Jinv[3][3],detJ=0.0;
            if(!pri6_eval_on_tri_face(xvol,lf,Nv_surface,Nv,gradN,Jinv,&detJ)) {
                free(owner);free(lface);return -1;
            }

            /* The surface-interpolated and canonical prism shape values should
               be identical on the triangular end face. */
            for(int a=0;a<6;++a) {
                if(fabs(Nv[a]-Nv_surface[a])>1.0e-9) {
                    fprintf(stderr,
                        "%s ERROR: PRI6 face interpolation mismatch es=%d owner=%d lf=%d qp=%d node=%d Nv=%e Ns=%e\n",
                        CODENAME,es,ke,lf,p,a,Nv[a],Nv_surface[a]);
                    free(owner);free(lface);return -1;
                }
            }

            double u[3]={0.0,0.0,0.0},uold[3]={0.0,0.0,0.0};
            double grad_u[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
            double p_q=0.0,grad_p[3]={0.0,0.0,0.0};

            for(int a=0;a<6;++a) {
                const int gid=fe_pri->conn[ke][a];
                for(int c=0;c<3;++c) {
                    u[c]+=Nv[a]*vals->v[gid][c];
                    uold[c]+=Nv[a]*vals->v_old[gid][c];
                    for(int d=0;d<3;++d) grad_u[c][d]+=gradN[a][d]*vals->v[gid][c];
                }
                p_q+=Nv[a]*vals->p[gid];
                for(int d=0;d<3;++d) grad_p[d]+=gradN[a][d]*vals->p[gid];
            }

            double Nexchange[6];
            if(opt->exchange_mode==PRI6_WALL_EXCHANGE_OPPOSITE_FACE)
                pri6_opposite_face_weights(lf,Nv,Nexchange);
            else
                for(int a=0;a<6;++a) Nexchange[a]=Nv[a];

            double u_exchange[3]={0.0,0.0,0.0};
            for(int a=0;a<6;++a) {
                const int gid=fe_pri->conn[ke][a];
                for(int c=0;c<3;++c) u_exchange[c]+=Nexchange[a]*vals->v[gid][c];
            }

            double* grad_uptr[3]={grad_u[0],grad_u[1],grad_u[2]};
            double mu_eff=mu_molecular,tau=0.0,tau_c=0.0;
            BBFE_vms_mu_eff_tau(&mu_eff,&tau,&tau_c,Jinv,detJ,h_e_vms,
                                u,uold,grad_uptr,grad_p,rho,mu_molecular,
                                vals->dt,vals->C_vms,vals->vms_cap_coeff);
            (void)tau;(void)tau_c;
            if(!isfinite(mu_eff)||mu_eff<=0.0) mu_eff=mu_molecular;

            const double w=basis_surf->integ_weight[p]*area_factor*Jface;
            double beta=0.0, beta_mu=0.0, beta_dt=0.0;
            if(BBFE_pri6_wall_mode_uses_nitsche(opt->wall_mode)) {
                pri6_get_penalty_components(
                    opt,rho,mu_eff,h_n,vals->dt,&beta,&beta_mu,&beta_dt);
            }
            const double beta_ratio =
                (beta_mu > 0.0) ? beta_dt/beta_mu : 0.0;

            if(diag) {
                diag->area+=w;
                diag->h_n_min=fmin(diag->h_n_min,h_n);
                diag->h_n_max=fmax(diag->h_n_max,h_n);
                diag->beta_min=fmin(diag->beta_min,beta);
                diag->beta_max=fmax(diag->beta_max,beta);
                diag->beta_mu_min=fmin(diag->beta_mu_min,beta_mu);
                diag->beta_mu_max=fmax(diag->beta_mu_max,beta_mu);
                diag->beta_dt_min=fmin(diag->beta_dt_min,beta_dt);
                diag->beta_dt_max=fmax(diag->beta_dt_max,beta_dt);
                diag->beta_dt_over_beta_mu_max=fmax(
                    diag->beta_dt_over_beta_mu_max,beta_ratio);
                diag->beta_area_sum += w*beta;
                diag->beta_mu_area_sum += w*beta_mu;
                diag->beta_dt_area_sum += w*beta_dt;
                diag->beta_ratio_area_sum += w*beta_ratio;
                const double rr[3]={u[0]-Uw[0],u[1]-Uw[1],u[2]-Uw[2]};
                diag->max_abs_normal_residual=fmax(diag->max_abs_normal_residual,
                                                   fabs(pri6_dot3(rr,n)));
            }

            for(int i=0;i<6;++i) {
                const int gi=fe_pri->conn[ke][i];
                double vec[4],dummy[4][4];
                const double zero[3]={0.0,0.0,0.0};
                pri6_wall_local_vecmat(vec,dummy,Nv[i],0.0,0.0,gradN[i],zero,n,
                                      u,u_exchange,p_q,grad_u,Uw,rho,mu_eff,mu_molecular,
                                      h_n,h_n,vals->dt,opt);

                for(int a=0;a<4;++a) {
                    const double val=w*vec[a];
                    monolis->mat.R.B[4*gi+a]-=val;
                    if(diag) diag->rhs_abs_sum+=fabs(val);
                }

                for(int j=0;j<6;++j) {
                    const int gj=fe_pri->conn[ke][j];
                    double dvec[4],mat[4][4];
                    pri6_wall_local_vecmat(dvec,mat,Nv[i],Nv[j],Nexchange[j],
                                          gradN[i],gradN[j],n,u,u_exchange,p_q,grad_u,Uw,
                                          rho,mu_eff,mu_molecular,h_n,h_n,vals->dt,opt);
                    for(int a=0;a<4;++a) for(int b=0;b<4;++b) {
                        const double val=w*mat[a][b];
                        monolis_add_scalar_to_sparse_matrix_R(monolis,gi,gj,a,b,val);
                        if(diag) diag->matrix_abs_sum+=fabs(val);
                    }
                }
            }

            ++qp_count;
            if(diag) ++diag->quadrature_points;
        }

        if(qp_count!=np) { free(owner);free(lface);return -1; }
        if(diag) { ++diag->mapped_faces; ++diag->assembled_faces; }
    }

    if(diag) {
        if(!isfinite(diag->h_n_min)) diag->h_n_min=0.0;
        if(!isfinite(diag->beta_min)) diag->beta_min=0.0;
        if(!isfinite(diag->beta_mu_min)) diag->beta_mu_min=0.0;
        if(!isfinite(diag->beta_dt_min)) diag->beta_dt_min=0.0;
    }

    free(owner);
    free(lface);
    return 0;
}


/* -------------------------------------------------------------------------- */
/* TET4/TRI3 owner map + Nitsche assembly                                     */
/* -------------------------------------------------------------------------- */

/*
 * Face numbering is only used as a local identifier.  Surface orientation is
 * taken from surf->conn and then flipped to point away from the owner cell, so
 * these triples do not need an outward orientation here.
 */
static const int TET4_TRI_FACE[4][3] = {
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
};

int BBFE_tet4_tri3_make_owner_map(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_tet,
    int** owner_elem_out,
    int** owner_lface_out,
    int* num_unmapped_local)
{
    if(owner_elem_out == NULL || owner_lface_out == NULL) return -1;

    *owner_elem_out = NULL;
    *owner_lface_out = NULL;
    if(num_unmapped_local) *num_unmapped_local = 0;

    if(surf == NULL || fe_tet == NULL) return -1;

    const int nsurf = surf->total_num_elems;
    if(nsurf < 0) return -1;
    if(nsurf == 0) return 0;

    if(surf->local_num_nodes != 3) {
        fprintf(stderr,
            "%s ERROR: active TET wall surface must be TRI3; got NEN=%d\n",
            CODENAME, surf->local_num_nodes);
        return -1;
    }

    if(!BBFE_fluid_mixed_has_elements(fe_tet)) {
        int* owner = (int*)malloc(sizeof(int)*(size_t)nsurf);
        int* lface = (int*)malloc(sizeof(int)*(size_t)nsurf);
        if(owner == NULL || lface == NULL) {
            free(owner); free(lface);
            return -1;
        }
        for(int es=0; es<nsurf; ++es) {
            owner[es] = -1;
            lface[es] = -1;
        }
        *owner_elem_out = owner;
        *owner_lface_out = lface;
        if(num_unmapped_local) *num_unmapped_local = nsurf;
        return 0;
    }

    if(fe_tet->local_num_nodes != 4) {
        fprintf(stderr,
            "%s ERROR: active tetra family must be TET4; got NEN=%d\n",
            CODENAME, fe_tet->local_num_nodes);
        return -1;
    }

    const int nvolface = 4 * fe_tet->total_num_elems;

    PRI6FaceEntry* vol =
        (PRI6FaceEntry*)malloc(sizeof(PRI6FaceEntry)*(size_t)nvolface);
    PRI6SurfEntry* srf =
        (PRI6SurfEntry*)malloc(sizeof(PRI6SurfEntry)*(size_t)nsurf);
    int* owner = (int*)malloc(sizeof(int)*(size_t)nsurf);
    int* lface = (int*)malloc(sizeof(int)*(size_t)nsurf);

    if(vol == NULL || srf == NULL || owner == NULL || lface == NULL) {
        free(vol); free(srf); free(owner); free(lface);
        return -1;
    }

    int k = 0;
    for(int e=0; e<fe_tet->total_num_elems; ++e) {
        for(int f=0; f<4; ++f) {
            pri6_node_id_t key[3] = {
                (pri6_node_id_t)fe_tet->conn[e][TET4_TRI_FACE[f][0]],
                (pri6_node_id_t)fe_tet->conn[e][TET4_TRI_FACE[f][1]],
                (pri6_node_id_t)fe_tet->conn[e][TET4_TRI_FACE[f][2]]
            };
            pri6_sort3(key);
            for(int a=0; a<3; ++a) vol[k].key[a] = key[a];
            vol[k].owner_elem = e;
            vol[k].owner_lface = f;
            ++k;
        }
    }

    for(int es=0; es<nsurf; ++es) {
        pri6_node_id_t key[3] = {
            (pri6_node_id_t)surf->conn[es][0],
            (pri6_node_id_t)surf->conn[es][1],
            (pri6_node_id_t)surf->conn[es][2]
        };
        pri6_sort3(key);
        for(int a=0; a<3; ++a) srf[es].key[a] = key[a];
        srf[es].surf_face_id = es;
        owner[es] = -1;
        lface[es] = -1;
    }

    qsort(vol, (size_t)nvolface, sizeof(PRI6FaceEntry), pri6_cmp_face_entry);
    qsort(srf, (size_t)nsurf, sizeof(PRI6SurfEntry), pri6_cmp_surf_entry);

    int i=0, j=0;
    while(i<nvolface && j<nsurf) {
        const int c = pri6_cmp3(vol[i].key, srf[j].key);
        if(c < 0) {
            pri6_node_id_t key0[3] = {vol[i].key[0],vol[i].key[1],vol[i].key[2]};
            do { ++i; } while(i<nvolface && pri6_cmp3(key0,vol[i].key)==0);
        }
        else if(c > 0) {
            pri6_node_id_t key0[3] = {srf[j].key[0],srf[j].key[1],srf[j].key[2]};
            do { ++j; } while(j<nsurf && pri6_cmp3(key0,srf[j].key)==0);
        }
        else {
            pri6_node_id_t key0[3] = {srf[j].key[0],srf[j].key[1],srf[j].key[2]};
            const int oe = vol[i].owner_elem;
            const int of = vol[i].owner_lface;
            do {
                const int es = srf[j].surf_face_id;
                owner[es] = oe;
                lface[es] = of;
                ++j;
            } while(j<nsurf && pri6_cmp3(key0,srf[j].key)==0);
            do { ++i; } while(i<nvolface && pri6_cmp3(key0,vol[i].key)==0);
        }
    }

    int unmapped = 0;
    for(int es=0; es<nsurf; ++es) {
        if(owner[es] < 0) ++unmapped;
    }

    free(vol);
    free(srf);

    *owner_elem_out = owner;
    *owner_lface_out = lface;
    if(num_unmapped_local) *num_unmapped_local = unmapped;

    return 0;
}

static int tet4_surface_to_volume_nodes(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    int es,
    int ke,
    int s2v[3])
{
    for(int a=0; a<3; ++a) {
        s2v[a] = -1;
        const int gid = surf->conn[es][a];
        for(int b=0; b<4; ++b) {
            if(fe->conn[ke][b] == gid) {
                s2v[a] = b;
                break;
            }
        }
        if(s2v[a] < 0) return 0;
    }
    return 1;
}

static int tet4_geometry(
    const double x[4][3],
    double gradN[4][3],
    double J_inv[3][3],
    double* detJ_abs,
    double* volume)
{
    static const double dNref[4][3] = {
        {-1.0, -1.0, -1.0},
        {+1.0,  0.0,  0.0},
        { 0.0, +1.0,  0.0},
        { 0.0,  0.0, +1.0}
    };

    double J[3][3] = {
        {0.0,0.0,0.0},
        {0.0,0.0,0.0},
        {0.0,0.0,0.0}
    };

    for(int a=0; a<4; ++a) {
        for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
                J[i][j] += x[a][i] * dNref[a][j];
            }
        }
    }

    const double detJ = pri6_det3(J);
    if(!isfinite(detJ) || fabs(detJ) < 1.0e-300) return 0;
    if(!pri6_inv3(J, J_inv)) return 0;

    for(int a=0; a<4; ++a) {
        for(int i=0; i<3; ++i) {
            gradN[a][i] = 0.0;
            for(int j=0; j<3; ++j) {
                gradN[a][i] += dNref[a][j] * J_inv[j][i];
            }
        }
    }

    *detJ_abs = fabs(detJ);
    *volume = fabs(detJ) / 6.0;
    return 1;
}

static void tet4_opposite_vertex_weights(int lface, double Nopp[4])
{
    for(int a=0; a<4; ++a) Nopp[a] = 0.0;
    if(lface < 0 || lface >= 4) return;

    for(int a=0; a<4; ++a) {
        int on_face = 0;
        for(int q=0; q<3; ++q) {
            if(TET4_TRI_FACE[lface][q] == a) {
                on_face = 1;
                break;
            }
        }
        if(!on_face) {
            Nopp[a] = 1.0;
            return;
        }
    }
}

int set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
    MONOLIS* monolis,
    const BBFE_DATA* fe_tet,
    const BBFE_BASIS* basis_tet,
    const BBFE_DATA* surf,
    const BBFE_BASIS* basis_surf,
    VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    const PRI6WallOptions* opt,
    PRI6NitscheDiagnostics* diag)
{
    (void)basis_tet;

    if(diag) {
        memset(diag, 0, sizeof(*diag));
        diag->h_n_min = HUGE_VAL;
        diag->beta_min = HUGE_VAL;
        diag->beta_mu_min = HUGE_VAL;
        diag->beta_dt_min = HUGE_VAL;
    }

    if(monolis==NULL || fe_tet==NULL || surf==NULL || basis_surf==NULL ||
       vals==NULL || Uw==NULL || !pri6_validate_wall_options(opt)) {
        fprintf(stderr,
            "%s ERROR: invalid input to TET4/TRI3 Nitsche assembly.\n",
            CODENAME);
        return -1;
    }

    if(surf->total_num_elems == 0) return 0;
    if(fe_tet->local_num_nodes != 4 || surf->local_num_nodes != 3) return -1;

    int *owner=NULL, *lface=NULL, unmapped=0;
    if(BBFE_tet4_tri3_make_owner_map(
            surf, fe_tet, &owner, &lface, &unmapped) != 0) {
        free(owner); free(lface);
        return -1;
    }

    const int np = basis_surf->num_integ_points;
    const double area_factor = pri6_tri_area_factor(basis_surf);
    if(np <= 0 || area_factor <= 0.0) {
        free(owner); free(lface);
        return -1;
    }

    if(diag) {
        diag->surface_faces = surf->total_num_elems;
        diag->skipped_no_local_owner = unmapped;
    }

    for(int es=0; es<surf->total_num_elems; ++es) {
        const int ke = owner[es];
        const int lf = lface[es];
        if(ke < 0) continue;
        if(ke >= fe_tet->total_num_elems || lf < 0 || lf >= 4) {
            free(owner); free(lface);
            return -1;
        }

        int s2v[3];
        if(!tet4_surface_to_volume_nodes(surf, fe_tet, es, ke, s2v)) {
            free(owner); free(lface);
            return -1;
        }

        for(int a=0; a<3; ++a) {
            int on_face = 0;
            for(int q=0; q<3; ++q) {
                if(s2v[a] == TET4_TRI_FACE[lf][q]) on_face = 1;
            }
            if(!on_face) {
                free(owner); free(lface);
                return -1;
            }
        }

        double xvol[4][3], xsurf[3][3];
        for(int a=0; a<4; ++a) {
            const int gid = fe_tet->conn[ke][a];
            for(int d=0; d<3; ++d) xvol[a][d] = fe_tet->x[gid][d];
        }
        for(int a=0; a<3; ++a) {
            const int gid = surf->conn[es][a];
            for(int d=0; d<3; ++d) xsurf[a][d] = fe_tet->x[gid][d];
        }

        double gradN[4][3], Jinv[3][3], detJ_abs=0.0, volume=0.0;
        if(!tet4_geometry(xvol, gradN, Jinv, &detJ_abs, &volume)) {
            free(owner); free(lface);
            return -1;
        }

        double xc[3]={0.0,0.0,0.0}, xf[3]={0.0,0.0,0.0};
        for(int a=0; a<4; ++a) for(int d=0; d<3; ++d) xc[d] += xvol[a][d]/4.0;
        for(int a=0; a<3; ++a) for(int d=0; d<3; ++d) xf[d] += xsurf[a][d]/3.0;
        const double outdir[3] = {xf[0]-xc[0], xf[1]-xc[1], xf[2]-xc[2]};

        double e1[3] = {
            xsurf[1][0]-xsurf[0][0],
            xsurf[1][1]-xsurf[0][1],
            xsurf[1][2]-xsurf[0][2]
        };
        double e2[3] = {
            xsurf[2][0]-xsurf[0][0],
            xsurf[2][1]-xsurf[0][1],
            xsurf[2][2]-xsurf[0][2]
        };
        double cr[3];
        pri6_cross3(e1, e2, cr);
        const double face_area = 0.5 * pri6_norm3(cr);

        if(!isfinite(face_area) || face_area <= 0.0 ||
           !isfinite(volume) || volume <= 0.0) {
            free(owner); free(lface);
            return -1;
        }

        /* For a tetrahedron V = A_face * h_normal / 3. */
        const double h_n =
            fmax(3.0*volume/fmax(face_area,1.0e-30), 1.0e-12);
        const double h_e_vms = cbrt(fmax(volume,1.0e-36));

        int qp_count = 0;
        for(int p=0; p<np; ++p) {
            double nraw[3];
            pri6_tri3_normal(
                xsurf,
                basis_surf->dN_dxi[p],
                basis_surf->dN_det[p],
                nraw);
            const double Jface = pri6_norm3(nraw);
            if(!isfinite(Jface) || Jface <= 1.0e-300) {
                free(owner); free(lface);
                return -1;
            }

            double n[3] = {
                nraw[0]/Jface,
                nraw[1]/Jface,
                nraw[2]/Jface
            };
            if(pri6_dot3(n,outdir) < 0.0) {
                n[0]*=-1.0; n[1]*=-1.0; n[2]*=-1.0;
            }

            double Nv[4] = {0.0,0.0,0.0,0.0};
            for(int a=0; a<3; ++a) {
                Nv[s2v[a]] = basis_surf->N[p][a];
            }

            double u[3]={0.0,0.0,0.0}, uold[3]={0.0,0.0,0.0};
            double grad_u[3][3] = {
                {0.0,0.0,0.0},
                {0.0,0.0,0.0},
                {0.0,0.0,0.0}
            };
            double p_q=0.0, grad_p[3]={0.0,0.0,0.0};

            for(int a=0; a<4; ++a) {
                const int gid = fe_tet->conn[ke][a];
                for(int c=0; c<3; ++c) {
                    u[c] += Nv[a]*vals->v[gid][c];
                    uold[c] += Nv[a]*vals->v_old[gid][c];
                    for(int d=0; d<3; ++d) {
                        grad_u[c][d] += gradN[a][d]*vals->v[gid][c];
                    }
                }
                p_q += Nv[a]*vals->p[gid];
                for(int d=0; d<3; ++d) {
                    grad_p[d] += gradN[a][d]*vals->p[gid];
                }
            }

            double Nexchange[4];
            if(opt->exchange_mode == PRI6_WALL_EXCHANGE_OPPOSITE_FACE) {
                /* TET4 analogue of the PRI6 opposite face: the opposite vertex. */
                tet4_opposite_vertex_weights(lf, Nexchange);
            }
            else {
                for(int a=0; a<4; ++a) Nexchange[a] = Nv[a];
            }

            double u_exchange[3]={0.0,0.0,0.0};
            for(int a=0; a<4; ++a) {
                const int gid = fe_tet->conn[ke][a];
                for(int c=0; c<3; ++c) {
                    u_exchange[c] += Nexchange[a]*vals->v[gid][c];
                }
            }

            double* grad_uptr[3] = {grad_u[0],grad_u[1],grad_u[2]};
            double mu_eff=mu_molecular, tau=0.0, tau_c=0.0;
            BBFE_vms_mu_eff_tau(
                &mu_eff, &tau, &tau_c,
                Jinv, detJ_abs, h_e_vms,
                u, uold, grad_uptr, grad_p,
                rho, mu_molecular,
                vals->dt, vals->C_vms, vals->vms_cap_coeff);
            (void)tau; (void)tau_c;
            if(!isfinite(mu_eff) || mu_eff <= 0.0) mu_eff = mu_molecular;

            const double w =
                basis_surf->integ_weight[p] * area_factor * Jface;
            double beta=0.0, beta_mu=0.0, beta_dt=0.0;
            if(BBFE_pri6_wall_mode_uses_nitsche(opt->wall_mode)) {
                pri6_get_penalty_components(
                    opt, rho, mu_eff, h_n, vals->dt,
                    &beta, &beta_mu, &beta_dt);
            }
            const double beta_ratio =
                (beta_mu > 0.0) ? beta_dt/beta_mu : 0.0;

            if(diag) {
                diag->area += w;
                diag->h_n_min = fmin(diag->h_n_min,h_n);
                diag->h_n_max = fmax(diag->h_n_max,h_n);
                diag->beta_min = fmin(diag->beta_min,beta);
                diag->beta_max = fmax(diag->beta_max,beta);
                diag->beta_mu_min = fmin(diag->beta_mu_min,beta_mu);
                diag->beta_mu_max = fmax(diag->beta_mu_max,beta_mu);
                diag->beta_dt_min = fmin(diag->beta_dt_min,beta_dt);
                diag->beta_dt_max = fmax(diag->beta_dt_max,beta_dt);
                diag->beta_dt_over_beta_mu_max = fmax(
                    diag->beta_dt_over_beta_mu_max,beta_ratio);
                diag->beta_area_sum += w*beta;
                diag->beta_mu_area_sum += w*beta_mu;
                diag->beta_dt_area_sum += w*beta_dt;
                diag->beta_ratio_area_sum += w*beta_ratio;
                const double rr[3] = {
                    u[0]-Uw[0],u[1]-Uw[1],u[2]-Uw[2]
                };
                diag->max_abs_normal_residual = fmax(
                    diag->max_abs_normal_residual,
                    fabs(pri6_dot3(rr,n)));
            }

            for(int i=0; i<4; ++i) {
                const int gi = fe_tet->conn[ke][i];
                double vec[4], dummy[4][4];
                const double zero[3]={0.0,0.0,0.0};
                pri6_wall_local_vecmat(
                    vec,dummy,
                    Nv[i],0.0,0.0,
                    gradN[i],zero,n,
                    u,u_exchange,p_q,grad_u,Uw,
                    rho,mu_eff,mu_molecular,
                    h_n,h_n,vals->dt,opt);

                for(int a=0; a<4; ++a) {
                    const double val = w*vec[a];
                    monolis->mat.R.B[4*gi+a] -= val;
                    if(diag) diag->rhs_abs_sum += fabs(val);
                }

                for(int j=0; j<4; ++j) {
                    const int gj = fe_tet->conn[ke][j];
                    double dvec[4], mat[4][4];
                    pri6_wall_local_vecmat(
                        dvec,mat,
                        Nv[i],Nv[j],Nexchange[j],
                        gradN[i],gradN[j],n,
                        u,u_exchange,p_q,grad_u,Uw,
                        rho,mu_eff,mu_molecular,
                        h_n,h_n,vals->dt,opt);

                    for(int a=0; a<4; ++a) {
                        for(int b=0; b<4; ++b) {
                            const double val = w*mat[a][b];
                            monolis_add_scalar_to_sparse_matrix_R(
                                monolis,gi,gj,a,b,val);
                            if(diag) diag->matrix_abs_sum += fabs(val);
                        }
                    }
                }
            }

            ++qp_count;
            if(diag) ++diag->quadrature_points;
        }

        if(qp_count != np) {
            free(owner); free(lface);
            return -1;
        }
        if(diag) {
            ++diag->mapped_faces;
            ++diag->assembled_faces;
        }
    }

    if(diag) {
        if(!isfinite(diag->h_n_min)) diag->h_n_min = 0.0;
        if(!isfinite(diag->beta_min)) diag->beta_min = 0.0;
        if(!isfinite(diag->beta_mu_min)) diag->beta_mu_min = 0.0;
        if(!isfinite(diag->beta_dt_min)) diag->beta_dt_min = 0.0;
    }

    free(owner);
    free(lface);
    return 0;
}

/* -------------------------------------------------------------------------- */
/* Mixed PRI6/TET4 projected-strong wall support                              */
/* -------------------------------------------------------------------------- */

typedef struct {
    int is_wall;
    int ncluster;
    int rank;
    double cluster_sum[3][3];
    double qn[3][3];
    double Pt[3][3];
} MixedWallNodeFeature;

static MixedWallNodeFeature* g_mixed_wall_feature = NULL;
static int g_mixed_wall_feature_nnode = 0;
static const BBFE_DATA* g_mixed_wall_projector_pri = NULL;
static const BBFE_DATA* g_mixed_wall_projector_tet = NULL;
static int g_mixed_wall_projector_use_pri = 0;
static int g_mixed_wall_projector_use_tet = 0;
static double g_mixed_wall_projector_angle = -1.0;

static void mixed_wall_projector_release(void)
{
    free(g_mixed_wall_feature);
    g_mixed_wall_feature = NULL;
    g_mixed_wall_feature_nnode = 0;
    g_mixed_wall_projector_pri = NULL;
    g_mixed_wall_projector_tet = NULL;
    g_mixed_wall_projector_use_pri = 0;
    g_mixed_wall_projector_use_tet = 0;
    g_mixed_wall_projector_angle = -1.0;
}

static void mixed_wall_normalize3(double a[3])
{
    const double nrm = pri6_norm3(a);
    if(nrm > 1.0e-30) {
        a[0] /= nrm;
        a[1] /= nrm;
        a[2] /= nrm;
    }
}

static int mixed_wall_add_surface_normals(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_ref,
    double cos_threshold)
{
    if(surf == NULL || surf->total_num_elems == 0) return 0;
    if(surf->local_num_nodes != 3 || fe_ref == NULL || fe_ref->x == NULL) return -1;

    for(int es=0; es<surf->total_num_elems; ++es) {
        const int g0 = surf->conn[es][0];
        const int g1 = surf->conn[es][1];
        const int g2 = surf->conn[es][2];
        if(g0 < 0 || g0 >= g_mixed_wall_feature_nnode ||
           g1 < 0 || g1 >= g_mixed_wall_feature_nnode ||
           g2 < 0 || g2 >= g_mixed_wall_feature_nnode) return -1;

        const double e1[3] = {
            fe_ref->x[g1][0]-fe_ref->x[g0][0],
            fe_ref->x[g1][1]-fe_ref->x[g0][1],
            fe_ref->x[g1][2]-fe_ref->x[g0][2]
        };
        const double e2[3] = {
            fe_ref->x[g2][0]-fe_ref->x[g0][0],
            fe_ref->x[g2][1]-fe_ref->x[g0][1],
            fe_ref->x[g2][2]-fe_ref->x[g0][2]
        };
        double n[3];
        pri6_cross3(e1,e2,n);
        const double nraw = pri6_norm3(n);
        if(!isfinite(nraw) || nraw <= 1.0e-30) continue;
        const double area = 0.5*nraw;
        n[0] /= nraw; n[1] /= nraw; n[2] /= nraw;

        const int gids[3] = {g0,g1,g2};
        for(int aa=0; aa<3; ++aa) {
            MixedWallNodeFeature* f = &g_mixed_wall_feature[gids[aa]];
            f->is_wall = 1;
            int best = -1;
            double best_abs_dot = -1.0;
            double best_signed_dot = 1.0;
            for(int k=0; k<f->ncluster; ++k) {
                double mean[3] = {
                    f->cluster_sum[k][0],
                    f->cluster_sum[k][1],
                    f->cluster_sum[k][2]
                };
                mixed_wall_normalize3(mean);
                const double d = pri6_dot3(mean,n);
                const double ad = fabs(d);
                if(ad > best_abs_dot) {
                    best_abs_dot = ad;
                    best_signed_dot = d;
                    best = k;
                }
            }
            if(best >= 0 && best_abs_dot >= cos_threshold) {
                const double sgn = best_signed_dot < 0.0 ? -1.0 : 1.0;
                for(int d=0; d<3; ++d) f->cluster_sum[best][d] += area*sgn*n[d];
            }
            else if(f->ncluster < 3) {
                const int k = f->ncluster++;
                for(int d=0; d<3; ++d) f->cluster_sum[k][d] = area*n[d];
            }
            else {
                if(best < 0) best = 0;
                const double sgn = best_signed_dot < 0.0 ? -1.0 : 1.0;
                for(int d=0; d<3; ++d) f->cluster_sum[best][d] += area*sgn*n[d];
            }
        }
    }
    return 0;
}

static int mixed_wall_projector_prepare(
    const BBFE_DATA* surf_pri,
    const BBFE_DATA* surf_tet,
    const FE_SYSTEM_FLUID_MIXED* sys,
    int use_pri,
    int use_tet)
{
    if(sys == NULL || sys->fe_tet.total_num_nodes <= 0 || sys->fe_tet.x == NULL) return -1;
    const int nnode = sys->fe_tet.total_num_nodes;
    const double angle = g_mixed_wall_feature_angle_deg;

    if(g_mixed_wall_feature != NULL &&
       g_mixed_wall_feature_nnode == nnode &&
       g_mixed_wall_projector_pri == surf_pri &&
       g_mixed_wall_projector_tet == surf_tet &&
       g_mixed_wall_projector_use_pri == use_pri &&
       g_mixed_wall_projector_use_tet == use_tet &&
       fabs(g_mixed_wall_projector_angle-angle) < 1.0e-12) return 0;

    mixed_wall_projector_release();
    g_mixed_wall_feature = (MixedWallNodeFeature*)calloc(
        (size_t)nnode, sizeof(MixedWallNodeFeature));
    if(g_mixed_wall_feature == NULL) return -1;
    g_mixed_wall_feature_nnode = nnode;
    g_mixed_wall_projector_pri = surf_pri;
    g_mixed_wall_projector_tet = surf_tet;
    g_mixed_wall_projector_use_pri = use_pri;
    g_mixed_wall_projector_use_tet = use_tet;
    g_mixed_wall_projector_angle = angle;

    const double cos_threshold = cos(angle*(3.14159265358979323846/180.0));
    if(use_pri && mixed_wall_add_surface_normals(surf_pri,&(sys->fe_tet),cos_threshold) != 0) {
        mixed_wall_projector_release(); return -1;
    }
    if(use_tet && mixed_wall_add_surface_normals(surf_tet,&(sys->fe_tet),cos_threshold) != 0) {
        mixed_wall_projector_release(); return -1;
    }

    long long rank_count[4] = {0,0,0,0};
    for(int gid=0; gid<nnode; ++gid) {
        MixedWallNodeFeature* f = &g_mixed_wall_feature[gid];
        for(int a=0; a<3; ++a) {
            for(int b=0; b<3; ++b) {
                f->Pt[a][b] = (a==b) ? 1.0 : 0.0;
                f->qn[a][b] = 0.0;
            }
        }
        if(!f->is_wall) continue;

        int rank = 0;
        for(int k=0; k<f->ncluster && rank<3; ++k) {
            double q[3] = {
                f->cluster_sum[k][0],f->cluster_sum[k][1],f->cluster_sum[k][2]
            };
            mixed_wall_normalize3(q);
            for(int j=0; j<rank; ++j) {
                const double alpha = pri6_dot3(q,f->qn[j]);
                for(int d=0; d<3; ++d) q[d] -= alpha*f->qn[j][d];
            }
            const double qnrm = pri6_norm3(q);
            if(qnrm <= 1.0e-6) continue;
            for(int d=0; d<3; ++d) f->qn[rank][d] = q[d]/qnrm;
            ++rank;
        }

        /* Existing component-wise strong rows are additional zero-increment
         * constraints at Ahmed/outer-boundary junction nodes. */
        if(sys->bc_NR.D_bc_exists != NULL) {
            for(int d=0; d<3 && rank<3; ++d) {
                if(!sys->bc_NR.D_bc_exists[4*gid+d]) continue;
                if(sys->bc.D_bc_exists != NULL && sys->bc.D_bc_exists[4*gid+d] &&
                   fabs(sys->bc.imposed_D_val[4*gid+d]) > 1.0e-12) {
                    fprintf(stderr,
                        "%s ERROR: projected-strong Ahmed wall intersects a nonzero "
                        "strong velocity row at node=%d component=%d value=%e.\n",
                        CODENAME,gid,d,sys->bc.imposed_D_val[4*gid+d]);
                    mixed_wall_projector_release(); return -1;
                }
                double q[3] = {0.0,0.0,0.0};
                q[d] = 1.0;
                for(int j=0; j<rank; ++j) {
                    const double alpha = pri6_dot3(q,f->qn[j]);
                    for(int c=0; c<3; ++c) q[c] -= alpha*f->qn[j][c];
                }
                const double qnrm = pri6_norm3(q);
                if(qnrm <= 1.0e-6) continue;
                for(int c=0; c<3; ++c) f->qn[rank][c] = q[c]/qnrm;
                ++rank;
            }
        }

        f->rank = rank;
        if(rank >= 0 && rank <= 3) ++rank_count[rank];
        for(int j=0; j<rank; ++j) {
            for(int a=0; a<3; ++a) {
                for(int b=0; b<3; ++b) {
                    f->Pt[a][b] -= f->qn[j][a]*f->qn[j][b];
                }
            }
        }
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s projected-strong mixed wall projector: feature_angle=%.3f deg "
            "rank1=%lld rank2=%lld rank3=%lld (rank-0 local copies)\n",
            CODENAME,angle,rank_count[1],rank_count[2],rank_count[3]);
    }
    return 0;
}

static void mixed_wall_project_relative_vector(int gid, double q[3])
{
    if(g_mixed_wall_feature == NULL || gid < 0 || gid >= g_mixed_wall_feature_nnode ||
       !g_mixed_wall_feature[gid].is_wall) return;
    const MixedWallNodeFeature* f = &g_mixed_wall_feature[gid];
    const double in[3] = {q[0],q[1],q[2]};
    for(int a=0; a<3; ++a) {
        q[a] = 0.0;
        for(int b=0; b<3; ++b) q[a] += f->Pt[a][b]*in[b];
    }
}

static void mixed_wall_project_velocity_state(double** v, const double Uw[3])
{
    if(v == NULL || Uw == NULL || g_mixed_wall_feature == NULL) return;
    for(int gid=0; gid<g_mixed_wall_feature_nnode; ++gid) {
        if(!g_mixed_wall_feature[gid].is_wall) continue;
        double q[3] = {v[gid][0]-Uw[0],v[gid][1]-Uw[1],v[gid][2]-Uw[2]};
        mixed_wall_project_relative_vector(gid,q);
        for(int d=0; d<3; ++d) v[gid][d] = Uw[d]+q[d];
    }
}

static void mixed_wall_project_velocity_increment(double** dv)
{
    if(dv == NULL || g_mixed_wall_feature == NULL) return;
    for(int gid=0; gid<g_mixed_wall_feature_nnode; ++gid) {
        if(!g_mixed_wall_feature[gid].is_wall) continue;
        double q[3] = {dv[gid][0],dv[gid][1],dv[gid][2]};
        mixed_wall_project_relative_vector(gid,q);
        for(int d=0; d<3; ++d) dv[gid][d] = q[d];
    }
}

static void mixed_wall_project_velocity_rhs(double* rhs)
{
    if(rhs == NULL || g_mixed_wall_feature == NULL) return;
    for(int gid=0; gid<g_mixed_wall_feature_nnode; ++gid) {
        if(!g_mixed_wall_feature[gid].is_wall) continue;
        double q[3] = {rhs[4*gid+0],rhs[4*gid+1],rhs[4*gid+2]};
        mixed_wall_project_relative_vector(gid,q);
        rhs[4*gid+0]=q[0]; rhs[4*gid+1]=q[1]; rhs[4*gid+2]=q[2];
    }
}

static int mixed_wall_apply_strong_noslip_family(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf,
    const double Uw[3])
{
    if(sys == NULL || surf == NULL || Uw == NULL) return -1;
    const int nnode = sys->fe_tet.total_num_nodes;
    if(sys->bc.D_bc_exists == NULL || sys->bc.imposed_D_val == NULL ||
       sys->bc_NR.D_bc_exists == NULL || sys->bc_NR.imposed_D_val == NULL) return -1;

    for(int es=0; es<surf->total_num_elems; ++es) {
        for(int a=0; a<surf->local_num_nodes; ++a) {
            const int gid = surf->conn[es][a];
            if(gid < 0 || gid >= nnode) return -1;
            for(int d=0; d<3; ++d) {
                const int dof = 4*gid+d;
                if(sys->bc.D_bc_exists[dof] &&
                   fabs(sys->bc.imposed_D_val[dof]-Uw[d]) > 1.0e-12) {
                    fprintf(stderr,
                        "%s ERROR: strong-noslip Ahmed wall conflicts with existing "
                        "velocity BC at node=%d component=%d old=%e wall=%e.\n",
                        CODENAME,gid,d,sys->bc.imposed_D_val[dof],Uw[d]);
                    return -1;
                }
                if(!sys->bc.D_bc_exists[dof]) ++sys->bc.num_D_bcs;
                sys->bc.D_bc_exists[dof] = 1;
                sys->bc.imposed_D_val[dof] = Uw[d];
                if(!sys->bc_NR.D_bc_exists[dof]) ++sys->bc_NR.num_D_bcs;
                sys->bc_NR.D_bc_exists[dof] = 1;
                sys->bc_NR.imposed_D_val[dof] = 0.0;
            }
            for(int d=0; d<3; ++d) {
                sys->vals.v[gid][d] = Uw[d];
                sys->vals.v_old[gid][d] = Uw[d];
            }
        }
    }
    return 0;
}

static int mixed_wall_mode_uses_surface_operator(int mode)
{
    return BBFE_pri6_wall_mode_uses_nitsche(mode) ||
           BBFE_pri6_wall_mode_uses_spalding(mode);
}

/* -------------------------------------------------------------------------- */
/* Mixed volume + prism-wall Newton/VMS driver                                */
/* -------------------------------------------------------------------------- */

typedef struct {
    int    local_prism_elements;

    /* PRI6-only RHS increment, separated by [ux,uy,uz,p]. */
    double rhs_l1[4];
    double rhs_linf[4];
} PRI6VolumeAssemblyDiagnostics;

/*
 * IMPORTANT:
 * Use exactly the same SUPS/PSPG/LSIC Newton operator as the established
 * single-element FOM.  set_element_mat_NR_Tezuer() already uses
 * fe->local_num_nodes dynamically, so the same implementation is valid for
 * TET4 and PRI6 as long as basis/geometry are correct.
 *
 * The previous mixed driver called the experimental split VMS routines here.
 * That changed the governing discrete operator only in the mixed executable
 * and made component-by-component debugging very difficult.  Do not do that.
 *
 * set_element_mat_NR_Tezuer() assembles BOTH the tangent matrix and the
 * residual vector.  Therefore this routine is used both for the Newton
 * assembly and for residual re-evaluation.  On residual-only re-evaluation
 * the matrix entries are harmless because the matrix is discarded afterwards.
 */
static void pri6_assemble_volume_sups_mixed(
    FE_SYSTEM_FLUID_MIXED*         sys,
    PRI6VolumeAssemblyDiagnostics* diag)
{
    const int has_tet = BBFE_fluid_mixed_has_elements(&(sys->fe_tet));
    const int has_pri = BBFE_fluid_mixed_has_elements(&(sys->fe_pri));
    const int nnode   = sys->fe_tet.total_num_nodes;
    const int ndof    = 4*nnode;

    if(diag != NULL) {
        memset(diag, 0, sizeof(*diag));
        diag->local_prism_elements =
            has_pri ? sys->fe_pri.total_num_elems : 0;
    }

    /* First assemble TET4 using the established FOM operator. */
    if(has_tet) {
        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe_tet),
            &(sys->basis_tet),
            &(sys->vals));
    }

    if(has_pri) {
        double* before = NULL;

        if(diag != NULL) {
            before = (double*)malloc(sizeof(double)*(size_t)ndof);
            if(before == NULL) {
                fprintf(stderr,
                    "%s ERROR: PRI6 RHS diagnostic allocation failed.\n",
                    CODENAME);
                exit(EXIT_FAILURE);
            }
            memcpy(before,
                   sys->monolis.mat.R.B,
                   sizeof(double)*(size_t)ndof);
        }

        /* Same SUPS/PSPG/LSIC Newton operator, now with PRI6 basis. */
        set_element_mat_NR_Tezuer(
            &(sys->monolis),
            &(sys->fe_pri),
            &(sys->basis_pri),
            &(sys->vals));

        if(diag != NULL) {
            for(int i = 0; i < nnode; ++i) {
                for(int d = 0; d < 4; ++d) {
                    const int id = 4*i + d;
                    const double a = fabs(
                        sys->monolis.mat.R.B[id] - before[id]);

                    diag->rhs_l1[d] += a;
                    if(a > diag->rhs_linf[d]) {
                        diag->rhs_linf[d] = a;
                    }
                }
            }
            free(before);
        }
    }
}

/* -------------------------------------------------------------------------- */
/* Physical/internal wall-surface SUM diagnostics                             */
/* -------------------------------------------------------------------------- */

static int mixed_wall_tri3_area_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe,
    double* area_out)
{
    if(area_out == NULL || surf == NULL || fe == NULL) return -1;

    *area_out = 0.0;

    if(surf->total_num_elems == 0) return 0;
    if(surf->local_num_nodes != 3) return -1;

    for(int e=0; e<surf->total_num_elems; ++e) {
        const int i0 = surf->conn[e][0];
        const int i1 = surf->conn[e][1];
        const int i2 = surf->conn[e][2];

        if(i0 < 0 || i0 >= fe->total_num_nodes ||
           i1 < 0 || i1 >= fe->total_num_nodes ||
           i2 < 0 || i2 >= fe->total_num_nodes) {
            return -1;
        }

        const double e1[3] = {
            fe->x[i1][0] - fe->x[i0][0],
            fe->x[i1][1] - fe->x[i0][1],
            fe->x[i1][2] - fe->x[i0][2]
        };
        const double e2[3] = {
            fe->x[i2][0] - fe->x[i0][0],
            fe->x[i2][1] - fe->x[i0][1],
            fe->x[i2][2] - fe->x[i0][2]
        };
        double cr[3];
        pri6_cross3(e1, e2, cr);
        const double a = 0.5 * pri6_norm3(cr);

        if(!isfinite(a) || a <= 0.0) return -1;
        *area_out += a;
    }

    return 0;
}

int BBFE_mixed_wall_physical_surface_diagnostics_global(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_DATA* fe_pri,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_DATA* fe_tet,
    MONOLIS_COM* mono_com,
    MixedWallPhysicalSurfaceDiagnostics* diag)
{
    if(surf_pri_internal == NULL || fe_pri == NULL ||
       surf_tet_internal == NULL || fe_tet == NULL ||
       mono_com == NULL || diag == NULL) {
        return -1;
    }

    memset(diag, 0, sizeof(*diag));

    double pri_area_local = 0.0;
    double tet_area_local = 0.0;
    int local_bad = 0;

    if(mixed_wall_tri3_area_local(
            surf_pri_internal, fe_pri, &pri_area_local) != 0) {
        local_bad = 1;
    }
    if(mixed_wall_tri3_area_local(
            surf_tet_internal, fe_tet, &tet_area_local) != 0) {
        local_bad = 1;
    }

    /* Every rank participates before any rank returns. */
    int global_bad = local_bad;
    monolis_allreduce_I(1, &global_bad, MONOLIS_MPI_SUM, mono_com->comm);
    if(global_bad != 0) return -1;

    int pri_faces = surf_pri_internal->total_num_elems;
    int tet_faces = surf_tet_internal->total_num_elems;

    monolis_allreduce_I(1, &pri_faces, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &tet_faces, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &pri_area_local, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &tet_area_local, MONOLIS_MPI_SUM, mono_com->comm);

    diag->pri_faces = pri_faces;
    diag->tet_faces = tet_faces;
    diag->total_faces = (long long)pri_faces + (long long)tet_faces;
    diag->pri_area = pri_area_local;
    diag->tet_area = tet_area_local;
    diag->total_area = pri_area_local + tet_area_local;

    return 0;
}


/* -------------------------------------------------------------------------- */
/* Mixed Ahmed-wall y+ / wall-state diagnostics                              */
/* -------------------------------------------------------------------------- */

typedef struct {
    long long faces;
    long long quadrature_points;

    double area;

    double int_h_n;
    double h_n_min;
    double h_n_max;

    double int_tau_w;
    double tau_w_min;
    double tau_w_max;

    double int_u_tau;
    double u_tau_min;
    double u_tau_max;

    double int_y_plus;
    double y_plus_min;
    double y_plus_max;

    double int_abs_u_normal;
    double abs_u_normal_max;

    double int_u_tangent;
    double u_tangent_max;

    double int_pressure;
    double pressure_min;
    double pressure_max;
} MixedWallDiagAccum;

static void mixed_wall_diag_accum_init(MixedWallDiagAccum* a)
{
    if(a == NULL) return;
    memset(a, 0, sizeof(*a));
    a->h_n_min = DBL_MAX;
    a->tau_w_min = DBL_MAX;
    a->u_tau_min = DBL_MAX;
    a->y_plus_min = DBL_MAX;
    a->pressure_min = DBL_MAX;
    a->pressure_max = -DBL_MAX;
}

static void mixed_wall_diag_accum_sample(
    MixedWallDiagAccum* a,
    double w,
    double h_n,
    const double u[3],
    double p,
    const double grad_u[3][3],
    const double n[3],
    double rho,
    double mu_molecular,
    const double Uw[3])
{
    if(a == NULL || w <= 0.0) return;

    const double rel[3] = {
        u[0] - Uw[0],
        u[1] - Uw[1],
        u[2] - Uw[2]
    };

    const double un = pri6_dot3(rel, n);
    const double ut[3] = {
        rel[0] - un*n[0],
        rel[1] - un*n[1],
        rel[2] - un*n[2]
    };
    const double ut_mag = pri6_norm3(ut);

    /* Molecular viscous Cauchy traction t_mu = mu (grad u + grad u^T) n. */
    double gun[3] = {0.0,0.0,0.0};
    double gTun[3] = {0.0,0.0,0.0};
    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            gun[i]  += grad_u[i][j] * n[j];
            gTun[i] += grad_u[j][i] * n[j];
        }
    }

    double traction_mu[3];
    for(int i=0; i<3; ++i) {
        traction_mu[i] = mu_molecular * (gun[i] + gTun[i]);
    }

    const double tn = pri6_dot3(traction_mu, n);
    const double tau_vec[3] = {
        traction_mu[0] - tn*n[0],
        traction_mu[1] - tn*n[1],
        traction_mu[2] - tn*n[2]
    };

    const double tau_w = pri6_norm3(tau_vec);
    const double u_tau =
        (rho > 0.0 && tau_w >= 0.0) ? sqrt(tau_w / rho) : 0.0;
    const double y_plus =
        (mu_molecular > 0.0) ? rho*u_tau*h_n/mu_molecular : 0.0;

    a->area += w;
    a->quadrature_points += 1;

    a->int_h_n += w*h_n;
    a->h_n_min = fmin(a->h_n_min, h_n);
    a->h_n_max = fmax(a->h_n_max, h_n);

    a->int_tau_w += w*tau_w;
    a->tau_w_min = fmin(a->tau_w_min, tau_w);
    a->tau_w_max = fmax(a->tau_w_max, tau_w);

    a->int_u_tau += w*u_tau;
    a->u_tau_min = fmin(a->u_tau_min, u_tau);
    a->u_tau_max = fmax(a->u_tau_max, u_tau);

    a->int_y_plus += w*y_plus;
    a->y_plus_min = fmin(a->y_plus_min, y_plus);
    a->y_plus_max = fmax(a->y_plus_max, y_plus);

    a->int_abs_u_normal += w*fabs(un);
    a->abs_u_normal_max = fmax(a->abs_u_normal_max, fabs(un));

    a->int_u_tangent += w*ut_mag;
    a->u_tangent_max = fmax(a->u_tangent_max, ut_mag);

    a->int_pressure += w*p;
    a->pressure_min = fmin(a->pressure_min, p);
    a->pressure_max = fmax(a->pressure_max, p);
}

static int mixed_wall_diag_pri_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    MixedWallDiagAccum* a)
{
    if(surf == NULL || fe_pri == NULL || basis_surf == NULL ||
       vals == NULL || Uw == NULL || a == NULL ||
       rho <= 0.0 || mu_molecular <= 0.0) return -1;

    mixed_wall_diag_accum_init(a);
    if(surf->total_num_elems == 0) return 0;
    if(surf->local_num_nodes != 3 || fe_pri->local_num_nodes != 6) return -1;

    int *owner=NULL, *lface=NULL, unmapped=0;
    if(BBFE_pri6_tri3_make_owner_map(
        surf, fe_pri, &owner, &lface, &unmapped) != 0) {
        free(owner); free(lface); return -1;
    }
    if(unmapped != 0) {
        free(owner); free(lface); return -1;
    }

    const int np = basis_surf->num_integ_points;
    const double area_factor = pri6_tri_area_factor(basis_surf);
    if(np <= 0 || area_factor <= 0.0) {
        free(owner); free(lface); return -1;
    }

    for(int es=0; es<surf->total_num_elems; ++es) {
        const int ke = owner[es];
        const int lf = lface[es];
        if(ke < 0 || ke >= fe_pri->total_num_elems || (lf != 0 && lf != 1)) {
            free(owner); free(lface); return -1;
        }

        int s2v[3];
        if(!pri6_surface_to_volume_nodes(surf, fe_pri, es, ke, s2v)) {
            free(owner); free(lface); return -1;
        }

        double xvol[6][3], xsurf[3][3];
        for(int ia=0; ia<6; ++ia) {
            const int gid = fe_pri->conn[ke][ia];
            for(int d=0; d<3; ++d) xvol[ia][d] = fe_pri->x[gid][d];
        }
        for(int ia=0; ia<3; ++ia) {
            const int gid = surf->conn[es][ia];
            for(int d=0; d<3; ++d) xsurf[ia][d] = fe_pri->x[gid][d];
        }

        double xc[3]={0.0,0.0,0.0}, xf[3]={0.0,0.0,0.0};
        for(int ia=0; ia<6; ++ia) for(int d=0; d<3; ++d) xc[d] += xvol[ia][d]/6.0;
        for(int ia=0; ia<3; ++ia) for(int d=0; d<3; ++d) xf[d] += xsurf[ia][d]/3.0;
        const double outdir[3] = {xf[0]-xc[0], xf[1]-xc[1], xf[2]-xc[2]};

        double e1[3] = {
            xsurf[1][0]-xsurf[0][0],
            xsurf[1][1]-xsurf[0][1],
            xsurf[1][2]-xsurf[0][2]
        };
        double e2[3] = {
            xsurf[2][0]-xsurf[0][0],
            xsurf[2][1]-xsurf[0][1],
            xsurf[2][2]-xsurf[0][2]
        };
        double cr[3];
        pri6_cross3(e1, e2, cr);
        const double face_area = 0.5*pri6_norm3(cr);
        const double volume = pri6_physical_volume(xvol);
        if(!isfinite(face_area) || face_area <= 0.0 ||
           !isfinite(volume) || volume <= 0.0) {
            free(owner); free(lface); return -1;
        }
        const double h_n = fmax(volume/fmax(face_area,1.0e-30), 1.0e-12);

        for(int p=0; p<np; ++p) {
            double nraw[3];
            pri6_tri3_normal(
                xsurf, basis_surf->dN_dxi[p], basis_surf->dN_det[p], nraw);
            const double Jface = pri6_norm3(nraw);
            if(!isfinite(Jface) || Jface <= 1.0e-300) {
                free(owner); free(lface); return -1;
            }

            double n[3] = {nraw[0]/Jface,nraw[1]/Jface,nraw[2]/Jface};
            if(pri6_dot3(n,outdir) < 0.0) {
                n[0]*=-1.0; n[1]*=-1.0; n[2]*=-1.0;
            }

            double Nv_surface[6] = {0,0,0,0,0,0};
            for(int ia=0; ia<3; ++ia) {
                Nv_surface[s2v[ia]] = basis_surf->N[p][ia];
            }

            double Nv[6], gradN[6][3], Jinv[3][3], detJ_abs=0.0;
            if(!pri6_eval_on_tri_face(
                xvol, lf, Nv_surface, Nv, gradN, Jinv, &detJ_abs)) {
                free(owner); free(lface); return -1;
            }

            double u[3]={0.0,0.0,0.0};
            double grad_u[3][3] = {
                {0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}
            };
            double p_q = 0.0;
            for(int ia=0; ia<6; ++ia) {
                const int gid = fe_pri->conn[ke][ia];
                for(int c=0; c<3; ++c) {
                    u[c] += Nv[ia]*vals->v[gid][c];
                    for(int d=0; d<3; ++d) {
                        grad_u[c][d] += gradN[ia][d]*vals->v[gid][c];
                    }
                }
                p_q += Nv[ia]*vals->p[gid];
            }

            const double w = basis_surf->integ_weight[p]*area_factor*Jface;
            mixed_wall_diag_accum_sample(
                a, w, h_n, u, p_q, grad_u, n,
                rho, mu_molecular, Uw);
        }
        a->faces += 1;
    }

    free(owner);
    free(lface);
    return 0;
}

static int mixed_wall_diag_tet_local(
    const BBFE_DATA* surf,
    const BBFE_DATA* fe_tet,
    const BBFE_BASIS* basis_surf,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    MixedWallDiagAccum* a)
{
    if(surf == NULL || fe_tet == NULL || basis_surf == NULL ||
       vals == NULL || Uw == NULL || a == NULL ||
       rho <= 0.0 || mu_molecular <= 0.0) return -1;

    mixed_wall_diag_accum_init(a);
    if(surf->total_num_elems == 0) return 0;
    if(surf->local_num_nodes != 3 || fe_tet->local_num_nodes != 4) return -1;

    int *owner=NULL, *lface=NULL, unmapped=0;
    if(BBFE_tet4_tri3_make_owner_map(
        surf, fe_tet, &owner, &lface, &unmapped) != 0) {
        free(owner); free(lface); return -1;
    }
    if(unmapped != 0) {
        free(owner); free(lface); return -1;
    }

    const int np = basis_surf->num_integ_points;
    const double area_factor = pri6_tri_area_factor(basis_surf);
    if(np <= 0 || area_factor <= 0.0) {
        free(owner); free(lface); return -1;
    }

    for(int es=0; es<surf->total_num_elems; ++es) {
        const int ke = owner[es];
        const int lf = lface[es];
        if(ke < 0 || ke >= fe_tet->total_num_elems || lf < 0 || lf >= 4) {
            free(owner); free(lface); return -1;
        }

        int s2v[3];
        if(!tet4_surface_to_volume_nodes(surf, fe_tet, es, ke, s2v)) {
            free(owner); free(lface); return -1;
        }

        double xvol[4][3], xsurf[3][3];
        for(int ia=0; ia<4; ++ia) {
            const int gid = fe_tet->conn[ke][ia];
            for(int d=0; d<3; ++d) xvol[ia][d] = fe_tet->x[gid][d];
        }
        for(int ia=0; ia<3; ++ia) {
            const int gid = surf->conn[es][ia];
            for(int d=0; d<3; ++d) xsurf[ia][d] = fe_tet->x[gid][d];
        }

        double gradN[4][3], Jinv[3][3], detJ_abs=0.0, volume=0.0;
        if(!tet4_geometry(xvol, gradN, Jinv, &detJ_abs, &volume)) {
            free(owner); free(lface); return -1;
        }

        double xc[3]={0.0,0.0,0.0}, xf[3]={0.0,0.0,0.0};
        for(int ia=0; ia<4; ++ia) for(int d=0; d<3; ++d) xc[d] += xvol[ia][d]/4.0;
        for(int ia=0; ia<3; ++ia) for(int d=0; d<3; ++d) xf[d] += xsurf[ia][d]/3.0;
        const double outdir[3] = {xf[0]-xc[0],xf[1]-xc[1],xf[2]-xc[2]};

        double e1[3] = {
            xsurf[1][0]-xsurf[0][0],
            xsurf[1][1]-xsurf[0][1],
            xsurf[1][2]-xsurf[0][2]
        };
        double e2[3] = {
            xsurf[2][0]-xsurf[0][0],
            xsurf[2][1]-xsurf[0][1],
            xsurf[2][2]-xsurf[0][2]
        };
        double cr[3];
        pri6_cross3(e1, e2, cr);
        const double face_area = 0.5*pri6_norm3(cr);
        if(!isfinite(face_area) || face_area <= 0.0 ||
           !isfinite(volume) || volume <= 0.0) {
            free(owner); free(lface); return -1;
        }
        const double h_n = fmax(3.0*volume/fmax(face_area,1.0e-30), 1.0e-12);

        for(int p=0; p<np; ++p) {
            double nraw[3];
            pri6_tri3_normal(
                xsurf, basis_surf->dN_dxi[p], basis_surf->dN_det[p], nraw);
            const double Jface = pri6_norm3(nraw);
            if(!isfinite(Jface) || Jface <= 1.0e-300) {
                free(owner); free(lface); return -1;
            }

            double n[3] = {nraw[0]/Jface,nraw[1]/Jface,nraw[2]/Jface};
            if(pri6_dot3(n,outdir) < 0.0) {
                n[0]*=-1.0; n[1]*=-1.0; n[2]*=-1.0;
            }

            double Nv[4] = {0,0,0,0};
            for(int ia=0; ia<3; ++ia) Nv[s2v[ia]] = basis_surf->N[p][ia];

            double u[3]={0.0,0.0,0.0};
            double grad_u[3][3] = {
                {0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}
            };
            double p_q = 0.0;
            for(int ia=0; ia<4; ++ia) {
                const int gid = fe_tet->conn[ke][ia];
                for(int c=0; c<3; ++c) {
                    u[c] += Nv[ia]*vals->v[gid][c];
                    for(int d=0; d<3; ++d) {
                        grad_u[c][d] += gradN[ia][d]*vals->v[gid][c];
                    }
                }
                p_q += Nv[ia]*vals->p[gid];
            }

            const double w = basis_surf->integ_weight[p]*area_factor*Jface;
            mixed_wall_diag_accum_sample(
                a, w, h_n, u, p_q, grad_u, n,
                rho, mu_molecular, Uw);
        }
        a->faces += 1;
    }

    free(owner);
    free(lface);
    return 0;
}

static void mixed_wall_diag_reduce_family(
    const MixedWallDiagAccum* local,
    MONOLIS_COM* mono_com,
    MixedWallFamilyDiagnostics* out)
{
    memset(out, 0, sizeof(*out));

    int faces = (int)local->faces;
    int qps = (int)local->quadrature_points;
    monolis_allreduce_I(1, &faces, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &qps, MONOLIS_MPI_SUM, mono_com->comm);

    double area = local->area;
    double int_h = local->int_h_n;
    double int_tau = local->int_tau_w;
    double int_utau = local->int_u_tau;
    double int_yp = local->int_y_plus;
    double int_un = local->int_abs_u_normal;
    double int_ut = local->int_u_tangent;
    double int_p = local->int_pressure;

    monolis_allreduce_R(1, &area, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_h, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_tau, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_utau, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_yp, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_un, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_ut, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &int_p, MONOLIS_MPI_SUM, mono_com->comm);

    double hmax = local->h_n_max;
    double taumax = local->tau_w_max;
    double utaumax = local->u_tau_max;
    double ypmax = local->y_plus_max;
    double unmax = local->abs_u_normal_max;
    double utmax = local->u_tangent_max;
    double pmax = (local->quadrature_points > 0) ? local->pressure_max : -DBL_MAX;

    monolis_allreduce_R(1, &hmax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &taumax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &utaumax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &ypmax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &unmax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &utmax, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &pmax, MONOLIS_MPI_MAX, mono_com->comm);

    double neg_hmin = (local->quadrature_points > 0) ? -local->h_n_min : -DBL_MAX;
    double neg_taumin = (local->quadrature_points > 0) ? -local->tau_w_min : -DBL_MAX;
    double neg_utaumin = (local->quadrature_points > 0) ? -local->u_tau_min : -DBL_MAX;
    double neg_ypmin = (local->quadrature_points > 0) ? -local->y_plus_min : -DBL_MAX;
    double neg_pmin = (local->quadrature_points > 0) ? -local->pressure_min : -DBL_MAX;

    monolis_allreduce_R(1, &neg_hmin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_taumin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_utaumin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_ypmin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_pmin, MONOLIS_MPI_MAX, mono_com->comm);

    out->faces = faces;
    out->quadrature_points = qps;
    out->area = area;

    if(area > 0.0) {
        out->h_n_mean = int_h/area;
        out->tau_w_mean = int_tau/area;
        out->u_tau_mean = int_utau/area;
        out->y_plus_mean = int_yp/area;
        out->abs_u_normal_mean = int_un/area;
        out->u_tangent_mean = int_ut/area;
        out->pressure_mean = int_p/area;
    }

    out->h_n_min = (neg_hmin > -DBL_MAX) ? -neg_hmin : 0.0;
    out->h_n_max = hmax;
    out->tau_w_min = (neg_taumin > -DBL_MAX) ? -neg_taumin : 0.0;
    out->tau_w_max = taumax;
    out->u_tau_min = (neg_utaumin > -DBL_MAX) ? -neg_utaumin : 0.0;
    out->u_tau_max = utaumax;
    out->y_plus_min = (neg_ypmin > -DBL_MAX) ? -neg_ypmin : 0.0;
    out->y_plus_max = ypmax;
    out->abs_u_normal_max = unmax;
    out->u_tangent_max = utmax;
    out->pressure_min = (neg_pmin > -DBL_MAX) ? -neg_pmin : 0.0;
    out->pressure_max = (pmax > -DBL_MAX) ? pmax : 0.0;
}

static void mixed_wall_diag_combine(
    const MixedWallFamilyDiagnostics* a,
    const MixedWallFamilyDiagnostics* b,
    MixedWallFamilyDiagnostics* t)
{
    memset(t, 0, sizeof(*t));
    t->faces = a->faces + b->faces;
    t->quadrature_points = a->quadrature_points + b->quadrature_points;
    t->area = a->area + b->area;

    if(t->area > 0.0) {
        t->h_n_mean = (a->h_n_mean*a->area + b->h_n_mean*b->area)/t->area;
        t->tau_w_mean = (a->tau_w_mean*a->area + b->tau_w_mean*b->area)/t->area;
        t->u_tau_mean = (a->u_tau_mean*a->area + b->u_tau_mean*b->area)/t->area;
        t->y_plus_mean = (a->y_plus_mean*a->area + b->y_plus_mean*b->area)/t->area;
        t->abs_u_normal_mean =
            (a->abs_u_normal_mean*a->area + b->abs_u_normal_mean*b->area)/t->area;
        t->u_tangent_mean =
            (a->u_tangent_mean*a->area + b->u_tangent_mean*b->area)/t->area;
        t->pressure_mean =
            (a->pressure_mean*a->area + b->pressure_mean*b->area)/t->area;
    }

    const int ha = a->quadrature_points > 0;
    const int hb = b->quadrature_points > 0;
    if(ha && hb) {
        t->h_n_min = fmin(a->h_n_min,b->h_n_min);
        t->tau_w_min = fmin(a->tau_w_min,b->tau_w_min);
        t->u_tau_min = fmin(a->u_tau_min,b->u_tau_min);
        t->y_plus_min = fmin(a->y_plus_min,b->y_plus_min);
        t->pressure_min = fmin(a->pressure_min,b->pressure_min);
    }
    else if(ha) {
        t->h_n_min=a->h_n_min; t->tau_w_min=a->tau_w_min;
        t->u_tau_min=a->u_tau_min; t->y_plus_min=a->y_plus_min;
        t->pressure_min=a->pressure_min;
    }
    else if(hb) {
        t->h_n_min=b->h_n_min; t->tau_w_min=b->tau_w_min;
        t->u_tau_min=b->u_tau_min; t->y_plus_min=b->y_plus_min;
        t->pressure_min=b->pressure_min;
    }

    t->h_n_max = fmax(a->h_n_max,b->h_n_max);
    t->tau_w_max = fmax(a->tau_w_max,b->tau_w_max);
    t->u_tau_max = fmax(a->u_tau_max,b->u_tau_max);
    t->y_plus_max = fmax(a->y_plus_max,b->y_plus_max);
    t->abs_u_normal_max = fmax(a->abs_u_normal_max,b->abs_u_normal_max);
    t->u_tangent_max = fmax(a->u_tangent_max,b->u_tangent_max);
    if(ha && hb) t->pressure_max = fmax(a->pressure_max,b->pressure_max);
    else if(ha) t->pressure_max = a->pressure_max;
    else if(hb) t->pressure_max = b->pressure_max;
}

int BBFE_mixed_wall_diagnostics_global(
    const BBFE_DATA* surf_pri_internal,
    const BBFE_DATA* fe_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet_internal,
    const BBFE_DATA* fe_tet,
    const BBFE_BASIS* basis_surf_tet,
    const VALUES* vals,
    double rho,
    double mu_molecular,
    const double Uw[3],
    MONOLIS_COM* mono_com,
    MixedWallDiagnostics* diag)
{
    if(surf_pri_internal == NULL || fe_pri == NULL || basis_surf_pri == NULL ||
       surf_tet_internal == NULL || fe_tet == NULL || basis_surf_tet == NULL ||
       vals == NULL || Uw == NULL || mono_com == NULL || diag == NULL ||
       rho <= 0.0 || mu_molecular <= 0.0) return -1;

    MixedWallDiagAccum pri_local, tet_local;
    int local_bad = 0;

    if(mixed_wall_diag_pri_local(
        surf_pri_internal, fe_pri, basis_surf_pri,
        vals, rho, mu_molecular, Uw, &pri_local) != 0) local_bad = 1;

    if(mixed_wall_diag_tet_local(
        surf_tet_internal, fe_tet, basis_surf_tet,
        vals, rho, mu_molecular, Uw, &tet_local) != 0) local_bad = 1;

    int global_bad = local_bad;
    monolis_allreduce_I(1, &global_bad, MONOLIS_MPI_SUM, mono_com->comm);
    if(global_bad != 0) return -1;

    memset(diag, 0, sizeof(*diag));
    mixed_wall_diag_reduce_family(&pri_local, mono_com, &diag->pri6);
    mixed_wall_diag_reduce_family(&tet_local, mono_com, &diag->tet4);
    mixed_wall_diag_combine(&diag->pri6, &diag->tet4, &diag->total);
    return 0;
}

static int mixed_wall_diag_path(
    char* out,
    size_t n,
    const char* directory,
    const char* filename)
{
    if(out == NULL || n == 0 || filename == NULL) return -1;
    if(directory == NULL || directory[0] == '\0' || strcmp(directory,".") == 0) {
        return snprintf(out,n,"%s",filename) < (int)n ? 0 : -1;
    }
    const size_t len = strlen(directory);
    const char* sep = (len > 0 && directory[len-1] == '/') ? "" : "/";
    return snprintf(out,n,"%s%s%s",directory,sep,filename) < (int)n ? 0 : -1;
}

void BBFE_mixed_reset_wall_diagnostics_file(const char* directory)
{
    char path[4096];
    if(mixed_wall_diag_path(
        path,sizeof(path),directory,"ahmed_wall_diagnostics.txt") != 0) return;

    FILE* fp = fopen(path,"w");
    if(fp == NULL) return;

    fprintf(fp,
        "# t family faces qps area "
        "h_min h_mean h_max "
        "tauw_min tauw_mean tauw_max "
        "utau_min utau_mean utau_max "
        "yplus_min yplus_mean yplus_max "
        "abs_un_mean abs_un_max ut_mean ut_max "
        "p_min p_mean p_max\n");
    fclose(fp);
}

static void mixed_wall_diag_write_family(
    FILE* fp,
    double t,
    const char* family,
    const MixedWallFamilyDiagnostics* d)
{
    fprintf(fp,
        "%.15e %s %lld %lld %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e "
        "%.15e %.15e %.15e %.15e "
        "%.15e %.15e %.15e\n",
        t,family,d->faces,d->quadrature_points,d->area,
        d->h_n_min,d->h_n_mean,d->h_n_max,
        d->tau_w_min,d->tau_w_mean,d->tau_w_max,
        d->u_tau_min,d->u_tau_mean,d->u_tau_max,
        d->y_plus_min,d->y_plus_mean,d->y_plus_max,
        d->abs_u_normal_mean,d->abs_u_normal_max,
        d->u_tangent_mean,d->u_tangent_max,
        d->pressure_min,d->pressure_mean,d->pressure_max);
}

void BBFE_mixed_output_wall_diagnostics(
    const MixedWallDiagnostics* diag,
    double t,
    const char* directory)
{
    if(diag == NULL) return;

    char path[4096];
    if(mixed_wall_diag_path(
        path,sizeof(path),directory,"ahmed_wall_diagnostics.txt") != 0) return;

    FILE* fp = fopen(path,"a");
    if(fp == NULL) return;

    mixed_wall_diag_write_family(fp,t,"PRI6",&diag->pri6);
    mixed_wall_diag_write_family(fp,t,"TET4",&diag->tet4);
    mixed_wall_diag_write_family(fp,t,"TOTAL",&diag->total);
    fclose(fp);
}

/*
 * Assembly diagnostics reduction.
 *
 * IMPORTANT:
 *   SUM-reduced fields (surface_faces, mapped_faces, assembled_faces, area,
 *   rhs_abs_sum, matrix_abs_sum, quadrature_points) count MPI surface overlap
 *   copies when assembly surfaces are duplicated.  They describe the actual
 *   distributed assembly work and MUST NOT be interpreted as unique physical
 *   wall area/force.
 *
 *   MIN/MAX fields (h_n_min/max, beta_min/max, max_abs_normal_residual) are
 *   duplicate-invariant and remain valid diagnostics.
 */
static void pri6_reduce_nitsche_diagnostics(
    const PRI6NitscheDiagnostics* local,
    PRI6NitscheDiagnostics*       global,
    MONOLIS_COM*                  mono_com)
{
    memset(global, 0, sizeof(*global));

    int surface_faces = (int)local->surface_faces;
    int mapped_faces = (int)local->mapped_faces;
    int skipped_faces = (int)local->skipped_no_local_owner;
    int assembled_faces = (int)local->assembled_faces;
    int qps = (int)local->quadrature_points;

    monolis_allreduce_I(1, &surface_faces,   MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &mapped_faces,    MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &skipped_faces,   MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &assembled_faces, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_I(1, &qps,             MONOLIS_MPI_SUM, mono_com->comm);

    global->surface_faces = surface_faces;
    global->mapped_faces = mapped_faces;
    global->skipped_no_local_owner = skipped_faces;
    global->assembled_faces = assembled_faces;
    global->quadrature_points = qps;

    double area = local->area;
    double rhs_l1 = local->rhs_abs_sum;
    double mat_l1 = local->matrix_abs_sum;
    double beta_sum = local->beta_area_sum;
    double beta_mu_sum = local->beta_mu_area_sum;
    double beta_dt_sum = local->beta_dt_area_sum;
    double beta_ratio_sum = local->beta_ratio_area_sum;

    monolis_allreduce_R(1, &area,   MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &rhs_l1, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &mat_l1, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &beta_sum, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &beta_mu_sum, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &beta_dt_sum, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &beta_ratio_sum, MONOLIS_MPI_SUM, mono_com->comm);

    global->area = area;
    global->rhs_abs_sum = rhs_l1;
    global->matrix_abs_sum = mat_l1;
    global->beta_area_sum = beta_sum;
    global->beta_mu_area_sum = beta_mu_sum;
    global->beta_dt_area_sum = beta_dt_sum;
    global->beta_ratio_area_sum = beta_ratio_sum;
    if(area > 0.0) {
        global->beta_mean = beta_sum/area;
        global->beta_mu_mean = beta_mu_sum/area;
        global->beta_dt_mean = beta_dt_sum/area;
        global->beta_dt_over_beta_mu_mean = beta_ratio_sum/area;
    }

    double hmax = local->h_n_max;
    double bmax = local->beta_max;
    double bmu_max = local->beta_mu_max;
    double bdt_max = local->beta_dt_max;
    double bratio_max = local->beta_dt_over_beta_mu_max;
    double unmax = local->max_abs_normal_residual;

    monolis_allreduce_R(1, &hmax,  MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &bmax,  MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &bmu_max, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &bdt_max, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &bratio_max, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &unmax, MONOLIS_MPI_MAX, mono_com->comm);

    global->h_n_max = hmax;
    global->beta_max = bmax;
    global->beta_mu_max = bmu_max;
    global->beta_dt_max = bdt_max;
    global->beta_dt_over_beta_mu_max = bratio_max;
    global->max_abs_normal_residual = unmax;

    /* No MIN reduction wrapper is assumed; reduce negative minima by MAX. */
    double neg_hmin = (local->assembled_faces > 0) ? -local->h_n_min : -DBL_MAX;
    double neg_bmin = (local->assembled_faces > 0) ? -local->beta_min : -DBL_MAX;
    double neg_bmu_min = (local->assembled_faces > 0) ? -local->beta_mu_min : -DBL_MAX;
    double neg_bdt_min = (local->assembled_faces > 0) ? -local->beta_dt_min : -DBL_MAX;

    monolis_allreduce_R(1, &neg_hmin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_bmin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_bmu_min, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_bdt_min, MONOLIS_MPI_MAX, mono_com->comm);

    global->h_n_min = (neg_hmin > -DBL_MAX) ? -neg_hmin : 0.0;
    global->beta_min = (neg_bmin > -DBL_MAX) ? -neg_bmin : 0.0;
    global->beta_mu_min = (neg_bmu_min > -DBL_MAX) ? -neg_bmu_min : 0.0;
    global->beta_dt_min = (neg_bdt_min > -DBL_MAX) ? -neg_bdt_min : 0.0;
}

void solver_fom_NR_mixed_prism_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_wall,
    const BBFE_BASIS* basis_wall,
    const PRI6WallOptions* wall_opt,
    double t,
    int step)
{
    (void)t;

    if(sys == NULL || surf_wall == NULL || basis_wall == NULL ||
       !pri6_validate_wall_options(wall_opt)) {
        fprintf(stderr, "%s ERROR: invalid mixed prism Nitsche solver input.\n", CODENAME);
        exit(EXIT_FAILURE);
    }

    const int nnode = sys->fe_tet.total_num_nodes;
    const int ndof = 4*nnode;
    const int nint = sys->mono_com.n_internal_vertex;
    const int max_iter_NR = 20;
    const double eps_NR = 1.0e-6;
    const double tiny = 1.0e-30;
    const double Uw[3] = {0.0, 0.0, 0.0};

    double* r0 = (double*)calloc((size_t)ndof, sizeof(double));
    if(r0 == NULL) exit(EXIT_FAILURE);

    for(int it = 0; it < max_iter_NR; ++it) {
        if(monolis_mpi_get_global_my_rank() == 0) {
            printf("\n%s ---- mixed PRI6 Nitsche step %d NR %d ----\n",
                   CODENAME, step, it);
        }

        /*
         * Critical for mixed MPI meshes:
         * elements on rank 3 can use overlap/interface nodes whose owner is a
         * different rank.  Always refresh v,p before element interpolation.
         */
        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i = 0; i < ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        PRI6VolumeAssemblyDiagnostics volume_diag;
        pri6_assemble_volume_sups_mixed(sys, &volume_diag);

        PRI6NitscheDiagnostics local_wall_diag;
        if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
            &(sys->monolis),
            &(sys->fe_pri),
            &(sys->basis_pri),
            surf_wall,
            basis_wall,
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            Uw,
            wall_opt,
            &local_wall_diag) != 0) {

            fprintf(stderr, "%s ERROR: PRI6/TRI3 Nitsche assembly failed.\n", CODENAME);
            free(r0);
            exit(EXIT_FAILURE);
        }

        int global_pri_elems = volume_diag.local_prism_elements;
        double global_rhs_l1[4];
        double global_rhs_linf[4];

        for(int d = 0; d < 4; ++d) {
            global_rhs_l1[d]   = volume_diag.rhs_l1[d];
            global_rhs_linf[d] = volume_diag.rhs_linf[d];
        }

        monolis_allreduce_I(
            1, &global_pri_elems, MONOLIS_MPI_SUM, sys->mono_com.comm);

        for(int d = 0; d < 4; ++d) {
            monolis_allreduce_R(
                1, &global_rhs_l1[d],
                MONOLIS_MPI_SUM, sys->mono_com.comm);
            monolis_allreduce_R(
                1, &global_rhs_linf[d],
                MONOLIS_MPI_MAX, sys->mono_com.comm);
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "[PRI6 SUPS step=%d NR=%d] elems=%d\n"
                "  ux: L1=%e Linf=%e\n"
                "  uy: L1=%e Linf=%e\n"
                "  uz: L1=%e Linf=%e\n"
                "   p: L1=%e Linf=%e\n",
                step, it, global_pri_elems,
                global_rhs_l1[0], global_rhs_linf[0],
                global_rhs_l1[1], global_rhs_linf[1],
                global_rhs_l1[2], global_rhs_linf[2],
                global_rhs_l1[3], global_rhs_linf[3]);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        if(it == 0) {
            memcpy(r0, sys->monolis.mat.R.B, sizeof(double)*(size_t)ndof);
        }

        BBFE_sys_monowrap_solve(
            &(sys->monolis),
            &(sys->mono_com),
            sys->monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        /*
         * Do not read overlap entries of X before synchronizing them.
         * This is especially important because all PRI6 cells are on one rank
         * but some TET/PRI interface nodes can be owned by neighbouring ranks.
         */
        BBFE_fluid_mixed_sync_solution_vector(
            &(sys->mono_com), nnode, sys->monolis.mat.R.X);

        BBFE_fluid_sups_renew_velocity(
            sys->vals.delta_v,
            sys->monolis.mat.R.X,
            nnode);

        BBFE_fluid_sups_renew_pressure(
            sys->vals.delta_p,
            sys->monolis.mat.R.X,
            nnode);

        BBFE_fluid_mixed_report_newton_diagnostics(sys, step, it);
        BBFE_fluid_mixed_report_prism_update(sys, step, it);

        update_velocity_pressure_NR(
            sys->vals.v,
            sys->vals.delta_v,
            sys->vals.p,
            sys->vals.delta_p,
            nnode);

        /* Make the updated state consistent before residual re-evaluation. */
        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        /*
         * Report the UPDATED state.  In v8 this was called before
         * update_velocity_pressure_NR(), which made NR=0 look like
         * uy=uz=p=0 even when the Newton increment was nonzero.
         */
        BBFE_fluid_mixed_report_prism_components(sys, step, it);

        /* Re-evaluate residual at the updated state. */
        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i = 0; i < ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        pri6_assemble_volume_sups_mixed(sys, NULL);

        if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
            &(sys->monolis),
            &(sys->fe_pri),
            &(sys->basis_pri),
            surf_wall,
            basis_wall,
            &(sys->vals),
            sys->vals.density,
            sys->vals.viscosity,
            Uw,
            wall_opt,
            NULL) != 0) {

            free(r0);
            exit(EXIT_FAILURE);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        double n0 = 0.0;
        double nr = 0.0;

        for(int i = 0; i < 4*nint; ++i) {
            n0 += r0[i]*r0[i];
            nr += sys->monolis.mat.R.B[i]*sys->monolis.mat.R.B[i];
        }

        monolis_allreduce_R(1, &n0, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &nr, MONOLIS_MPI_SUM, sys->mono_com.comm);

        const double rel = sqrt(nr)/fmax(sqrt(n0), tiny);

        PRI6NitscheDiagnostics global_wall_diag;
        pri6_reduce_nitsche_diagnostics(
            &local_wall_diag, &global_wall_diag, &(sys->mono_com));

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "[mixed-pri6-nitsche] residual rel=%e copy_area=%e "
                "h=[%e,%e] beta=[%e,%e] mapped=%lld skipped=%lld "
                "|A_wall|_1=%e |R_wall|_1=%e max|u.n|=%e\n",
                rel,
                global_wall_diag.area,
                global_wall_diag.h_n_min,
                global_wall_diag.h_n_max,
                global_wall_diag.beta_min,
                global_wall_diag.beta_max,
                global_wall_diag.mapped_faces,
                global_wall_diag.skipped_no_local_owner,
                global_wall_diag.matrix_abs_sum,
                global_wall_diag.rhs_abs_sum,
                global_wall_diag.max_abs_normal_residual);
        }

        if(rel < eps_NR) break;
    }

    /* v is synchronized at the end of every Newton iteration. */
    ROM_BB_vec_copy_2d(sys->vals.v, sys->vals.v_old, nnode, 3);
    free(r0);
}



/* -------------------------------------------------------------------------- */
/* Mixed TET4 + PRI6 VMS/Newton + dual TRI3 Nitsche driver                    */
/*                                                                            */
/* This follows the established solver_fom_VMS() algorithm:                    */
/*   1) VMS linear matrix/residual                                              */
/*   2) VMS nonlinear matrix/residual                                           */
/*   3) wall terms                                                              */
/*   4) strong Dirichlet                                                        */
/*   5) BiCGSAFE Newton correction                                              */
/*   6) updated-state residual re-evaluation using residual kernels only        */
/*                                                                            */
/* TET4 and PRI6 volume contributions are accumulated into the same Monolis     */
/* matrix by the *_VMS_mixed wrappers in core_FOM_mixed.c.                     */
/* -------------------------------------------------------------------------- */

static void mixed_vms_assemble_full(
    FE_SYSTEM_FLUID_MIXED* sys)
{
    set_element_mat_NR_linear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));

    set_element_vec_NR_linear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));

    set_element_mat_NR_nonlinear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));

    set_element_vec_NR_nonlinear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));
}

static void mixed_vms_assemble_residual_only(
    FE_SYSTEM_FLUID_MIXED* sys)
{
    set_element_vec_NR_linear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));

    set_element_vec_NR_nonlinear_VMS_mixed(
        &(sys->monolis),
        &(sys->fe_tet), &(sys->basis_tet),
        &(sys->fe_pri), &(sys->basis_pri),
        &(sys->vals));
}

static int mixed_vms_vector_is_finite_global(
    const double* vec,
    int n,
    MONOLIS_COM* mono_com,
    const char* label)
{
    int local_bad = 0;

    if(vec == NULL || n < 0 || mono_com == NULL) {
        local_bad = 1;
    }
    else {
        for(int i=0; i<n; ++i) {
            if(!isfinite(vec[i])) {
                local_bad = 1;
                break;
            }
        }
    }

    int global_bad = local_bad;
    if(monolis_mpi_get_global_comm_size() != 1) {
        monolis_allreduce_I(
            1, &global_bad, MONOLIS_MPI_SUM, mono_com->comm);
    }

    if(global_bad != 0 && monolis_mpi_get_global_my_rank() == 0) {
        fprintf(stderr,
            "%s ERROR: non-finite value detected in %s on %d rank(s).\n",
            CODENAME,
            label != NULL ? label : "vector",
            global_bad);
    }

    return global_bad == 0;
}


static void mixed_vms_report_rhs_components(
    const FE_SYSTEM_FLUID_MIXED* sys,
    int step,
    int it)
{
    if(sys == NULL || sys->monolis.mat.R.B == NULL) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int nint = (comm_size == 1) ? nnode : sys->mono_com.n_internal_vertex;
    if(nint <= 0) return;

    double sumsq[4] = {0.0,0.0,0.0,0.0};
    double sumabs[4] = {0.0,0.0,0.0,0.0};
    double linf[4] = {0.0,0.0,0.0,0.0};

    for(int i=0; i<nint; ++i) {
        for(int d=0; d<4; ++d) {
            const double x = sys->monolis.mat.R.B[4*i+d];
            const double ax = fabs(x);
            sumsq[d] += x*x;
            sumabs[d] += ax;
            if(ax > linf[d]) linf[d] = ax;
        }
    }

    if(comm_size != 1) {
        monolis_allreduce_R(4, sumsq, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(4, sumabs, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(4, linf, MONOLIS_MPI_MAX, sys->mono_com.comm);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "[Newton RHS step=%d NR=%d] "
            "ux:L1=%e L2=%e Linf=%e  "
            "uy:L1=%e L2=%e Linf=%e  "
            "uz:L1=%e L2=%e Linf=%e  "
            "p:L1=%e L2=%e Linf=%e\n",
            step, it,
            sumabs[0], sqrt(fmax(sumsq[0],0.0)), linf[0],
            sumabs[1], sqrt(fmax(sumsq[1],0.0)), linf[1],
            sumabs[2], sqrt(fmax(sumsq[2],0.0)), linf[2],
            sumabs[3], sqrt(fmax(sumsq[3],0.0)), linf[3]);
    }
}


static int mixed_vms_is_ahmed_case(void)
{
    const char* c = getenv("FLUID_CASE");
    if(c == NULL || c[0] == '\0') return 1;
    return strcmp(c,"ahmed") == 0 || strcmp(c,"AHMED") == 0 ||
           strcmp(c,"ahmedbody") == 0 || strcmp(c,"AHMEDBODY") == 0;
}

static void mixed_vms_history_path(
    char* path, size_t n, const char* directory)
{
    if(path == NULL || n == 0) return;
    if(directory == NULL || directory[0] == '\0' || strcmp(directory,".") == 0) {
        snprintf(path,n,"ahmed_solver_history.csv");
        return;
    }
    const size_t len = strlen(directory);
    snprintf(path,n,
        (len > 0 && directory[len-1] == '/') ?
            "%sahmed_solver_history.csv" : "%s/ahmed_solver_history.csv",
        directory);
}

static void mixed_vms_append_solver_history(
    FE_SYSTEM_FLUID_MIXED* sys,
    double t,
    int step,
    int it,
    double nrmr,
    double nrm0,
    double rel,
    int converged,
    const PRI6NitscheDiagnostics* pri,
    const PRI6NitscheDiagnostics* tet)
{
    if(sys == NULL || !mixed_vms_is_ahmed_case()) return;

    const int nnode = sys->fe_tet.total_num_nodes;
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int nint = (comm_size == 1) ? nnode : sys->mono_com.n_internal_vertex;
    if(nint < 0 || nint > nnode) return;

    double max_du = 0.0;
    double max_dp = 0.0;
    double max_abs_u[3] = {0.0,0.0,0.0};
    double pmin = DBL_MAX;
    double pmax = -DBL_MAX;
    double pabs = 0.0;

    for(int i=0; i<nint; ++i) {
        const double du2 =
            sys->vals.delta_v[i][0]*sys->vals.delta_v[i][0] +
            sys->vals.delta_v[i][1]*sys->vals.delta_v[i][1] +
            sys->vals.delta_v[i][2]*sys->vals.delta_v[i][2];
        max_du = fmax(max_du, sqrt(fmax(du2,0.0)));
        max_dp = fmax(max_dp, fabs(sys->vals.delta_p[i]));
        for(int d=0; d<3; ++d) {
            max_abs_u[d] = fmax(max_abs_u[d], fabs(sys->vals.v[i][d]));
        }
        pmin = fmin(pmin, sys->vals.p[i]);
        pmax = fmax(pmax, sys->vals.p[i]);
        pabs = fmax(pabs, fabs(sys->vals.p[i]));
    }

    if(comm_size != 1) {
        monolis_allreduce_R(1,&max_du,MONOLIS_MPI_MAX,sys->mono_com.comm);
        monolis_allreduce_R(1,&max_dp,MONOLIS_MPI_MAX,sys->mono_com.comm);
        for(int d=0; d<3; ++d)
            monolis_allreduce_R(1,&max_abs_u[d],MONOLIS_MPI_MAX,sys->mono_com.comm);
        monolis_allreduce_R(1,&pabs,MONOLIS_MPI_MAX,sys->mono_com.comm);

        double neg_pmin = (nint > 0) ? -pmin : -DBL_MAX;
        double pmax_red = (nint > 0) ? pmax : -DBL_MAX;
        monolis_allreduce_R(1,&neg_pmin,MONOLIS_MPI_MAX,sys->mono_com.comm);
        monolis_allreduce_R(1,&pmax_red,MONOLIS_MPI_MAX,sys->mono_com.comm);
        pmin = (neg_pmin > -DBL_MAX) ? -neg_pmin : 0.0;
        pmax = (pmax_red > -DBL_MAX) ? pmax_red : 0.0;
    }
    else if(nint == 0) {
        pmin = pmax = 0.0;
    }

    if(monolis_mpi_get_global_my_rank() != 0) return;

    char path[4096];
    mixed_vms_history_path(path,sizeof(path),sys->cond.directory);
    FILE* fp = fopen(path, (step == 1 && it == 0) ? "w" : "a");
    if(fp == NULL) {
        fprintf(stderr,"%s WARNING: cannot write %s\n",CODENAME,path);
        return;
    }

    if(step == 1 && it == 0) {
        fprintf(fp,
            "time,step,NR,converged,residual_nrm,residual_initial,residual_rel,"
            "max_delta_u,max_delta_p,max_abs_ux,max_abs_uy,max_abs_uz,"
            "p_min,p_max,max_abs_p,linear_max_iter,linear_tolerance,"
            "pri_beta_min,pri_beta_mean,pri_beta_max,"
            "pri_beta_mu_min,pri_beta_mu_mean,pri_beta_mu_max,"
            "pri_beta_dt_min,pri_beta_dt_mean,pri_beta_dt_max,"
            "pri_beta_dt_over_mu_mean,pri_beta_dt_over_mu_max,"
            "tet_beta_min,tet_beta_mean,tet_beta_max,"
            "tet_beta_mu_min,tet_beta_mu_mean,tet_beta_mu_max,"
            "tet_beta_dt_min,tet_beta_dt_mean,tet_beta_dt_max,"
            "tet_beta_dt_over_mu_mean,tet_beta_dt_over_mu_max\n");
    }

    PRI6NitscheDiagnostics z;
    memset(&z,0,sizeof(z));
    if(pri == NULL) pri = &z;
    if(tet == NULL) tet = &z;

    fprintf(fp,
        "%.17e,%d,%d,%d,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%d,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
        t,step,it,converged,nrmr,nrm0,rel,
        max_du,max_dp,max_abs_u[0],max_abs_u[1],max_abs_u[2],
        pmin,pmax,pabs,sys->vals.mat_max_iter,sys->vals.mat_epsilon,
        pri->beta_min,pri->beta_mean,pri->beta_max,
        pri->beta_mu_min,pri->beta_mu_mean,pri->beta_mu_max,
        pri->beta_dt_min,pri->beta_dt_mean,pri->beta_dt_max,
        pri->beta_dt_over_beta_mu_mean,pri->beta_dt_over_beta_mu_max,
        tet->beta_min,tet->beta_mean,tet->beta_max,
        tet->beta_mu_min,tet->beta_mu_mean,tet->beta_mu_max,
        tet->beta_dt_min,tet->beta_dt_mean,tet->beta_dt_max,
        tet->beta_dt_over_beta_mu_mean,tet->beta_dt_over_beta_mu_max);
    fclose(fp);
}

int solver_fom_VMS_mixed_dual_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet,
    const BBFE_BASIS* basis_surf_tet,
    int enable_pri6_nitsche,
    const PRI6WallOptions* wall_opt_pri,
    int enable_tet4_nitsche,
    const PRI6WallOptions* wall_opt_tet,
    double t,
    int step)
{
    if(sys == NULL ||
       surf_pri == NULL || basis_surf_pri == NULL ||
       surf_tet == NULL || basis_surf_tet == NULL ||
       (enable_pri6_nitsche &&
        (wall_opt_pri == NULL || !pri6_validate_wall_options(wall_opt_pri))) ||
       (enable_tet4_nitsche &&
        (wall_opt_tet == NULL || !pri6_validate_wall_options(wall_opt_tet)))) {
        fprintf(stderr,
            "%s ERROR: invalid mixed VMS split-Nitsche solver input.\n",
            CODENAME);
        return -1;
    }

    const int nnode = sys->fe_tet.total_num_nodes;
    const int ndof = 4*nnode;
    const int comm_size = monolis_mpi_get_global_comm_size();
    const int nint =
        (comm_size == 1)
        ? nnode
        : sys->mono_com.n_internal_vertex;

    const int max_iter_NR = 20;
    const double eps_NR = 1.0e-6;
    const double tiny = 1.0e-30;
    const double Uw[3] = {0.0, 0.0, 0.0};

    int global_pri_wall_faces = surf_pri->total_num_elems;
    int global_tet_wall_faces = surf_tet->total_num_elems;
    if(comm_size != 1) {
        monolis_allreduce_I(
            1, &global_pri_wall_faces, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_I(
            1, &global_tet_wall_faces, MONOLIS_MPI_SUM, sys->mono_com.comm);
    }

    const int use_pri_wall =
        (enable_pri6_nitsche != 0) && (global_pri_wall_faces > 0);
    const int use_tet_wall =
        (enable_tet4_nitsche != 0) && (global_tet_wall_faces > 0);

    const int pri_strong_noslip = use_pri_wall &&
        BBFE_pri6_wall_mode_is_strong_noslip(wall_opt_pri->wall_mode);
    const int tet_strong_noslip = use_tet_wall &&
        BBFE_pri6_wall_mode_is_strong_noslip(wall_opt_tet->wall_mode);
    const int pri_projected = use_pri_wall &&
        BBFE_pri6_wall_mode_uses_projected_strong_normal(wall_opt_pri->wall_mode);
    const int tet_projected = use_tet_wall &&
        BBFE_pri6_wall_mode_uses_projected_strong_normal(wall_opt_tet->wall_mode);
    const int assemble_pri_wall = use_pri_wall &&
        mixed_wall_mode_uses_surface_operator(wall_opt_pri->wall_mode);
    const int assemble_tet_wall = use_tet_wall &&
        mixed_wall_mode_uses_surface_operator(wall_opt_tet->wall_mode);
    const int use_projected = pri_projected || tet_projected;

    if(pri_strong_noslip &&
       mixed_wall_apply_strong_noslip_family(sys,surf_pri,Uw) != 0) return -1;
    if(tet_strong_noslip &&
       mixed_wall_apply_strong_noslip_family(sys,surf_tet,Uw) != 0) return -1;

    if(use_projected) {
        if(mixed_wall_projector_prepare(
            surf_pri,surf_tet,sys,pri_projected,tet_projected) != 0) {
            fprintf(stderr,
                "%s ERROR: failed to prepare mixed projected-strong wall projector.\n",
                CODENAME);
            return -1;
        }
        mixed_wall_project_velocity_state(sys->vals.v,Uw);
        mixed_wall_project_velocity_state(sys->vals.v_old,Uw);
        BBFE_fluid_mixed_sync_state(&(sys->mono_com),&(sys->vals),nnode);
    }

    if(monolis_mpi_get_global_my_rank() == 0) {
        printf(
            "%s wall runtime: PRI6=%s mode=%s faces=%d exchange=%s "
            "gamma=%.6g c_dt=%.6g; "
            "TET4=%s mode=%s faces=%d exchange=%s gamma=%.6g c_dt=%.6g\n",
            CODENAME,
            use_pri_wall ? "ON" : "OFF",
            use_pri_wall ? BBFE_pri6_wall_mode_name(wall_opt_pri->wall_mode) : "off",
            global_pri_wall_faces,
            use_pri_wall ? BBFE_pri6_wall_exchange_name(wall_opt_pri->exchange_mode) : "-",
            wall_opt_pri != NULL ? wall_opt_pri->gamma_n : 0.0,
            wall_opt_pri != NULL ? wall_opt_pri->dt_penalty_coeff : 0.0,
            use_tet_wall ? "ON" : "OFF",
            use_tet_wall ? BBFE_pri6_wall_mode_name(wall_opt_tet->wall_mode) : "off",
            global_tet_wall_faces,
            use_tet_wall ? BBFE_pri6_wall_exchange_name(wall_opt_tet->exchange_mode) : "-",
            wall_opt_tet != NULL ? wall_opt_tet->gamma_n : 0.0,
            wall_opt_tet != NULL ? wall_opt_tet->dt_penalty_coeff : 0.0);
    }

    if(nnode <= 0 || ndof <= 0 || nint < 0 || nint > nnode) {
        if(monolis_mpi_get_global_my_rank() == 0) {
            fprintf(stderr,
                "%s ERROR: invalid mixed VMS dimensions: "
                "nnode=%d nint=%d ndof=%d.\n",
                CODENAME, nnode, nint, ndof);
        }
        return -1;
    }

    double* r0 = (double*)calloc((size_t)ndof, sizeof(double));
    if(r0 == NULL) {
        fprintf(stderr, "%s ERROR: VMS residual allocation failed.\n", CODENAME);
        return -1;
    }

    int converged = 0;

    for(int it=0; it<max_iter_NR; ++it) {
        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "\n%s ---- mixed VMS TET4+PRI6 / runtime wall "
                "PRI6=%s(%s) TET4=%s(%s) step %d NR %d ----\n",
                CODENAME,
                use_pri_wall ? "ON" : "OFF",
                use_pri_wall ? BBFE_pri6_wall_mode_name(wall_opt_pri->wall_mode) : "off",
                use_tet_wall ? "ON" : "OFF",
                use_tet_wall ? BBFE_pri6_wall_mode_name(wall_opt_tet->wall_mode) : "off",
                step, it);
        }

        /* The VMS kernels evaluate both current and previous-time states. */
        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i=0; i<ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        /* Same four-stage VMS volume operator as solver_fom_VMS(). */
        mixed_vms_assemble_full(sys);

        if(assemble_pri_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
                &(sys->monolis),
                &(sys->fe_pri),
                &(sys->basis_pri),
                surf_pri,
                basis_surf_pri,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_pri,
                NULL) != 0) {

                fprintf(stderr,
                    "%s ERROR: PRI6/TRI3 Nitsche assembly failed in VMS solver.\n",
                    CODENAME);
                free(r0);
                return -1;
            }
        }

        if(assemble_tet_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
                &(sys->monolis),
                &(sys->fe_tet),
                &(sys->basis_tet),
                surf_tet,
                basis_surf_tet,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_tet,
                NULL) != 0) {

                fprintf(stderr,
                    "%s ERROR: TET4/TRI3 Nitsche assembly failed in VMS solver.\n",
                    CODENAME);
                free(r0);
                return -1;
            }
        }

        if(use_projected) {
            mixed_wall_project_velocity_rhs(sys->monolis.mat.R.B);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        if(!mixed_vms_vector_is_finite_global(
            sys->monolis.mat.R.B,
            ndof,
            &(sys->mono_com),
            "VMS Newton RHS before solve")) {
            free(r0);
            return -1;
        }

        if(it == 0) {
            memcpy(r0, sys->monolis.mat.R.B, sizeof(double)*(size_t)ndof);
        }

        mixed_vms_report_rhs_components(sys, step, it);

        /* Use the same linear-solver family as the established VMS solver. */
        ROM_monowrap_solve(
            &(sys->monolis),
            &(sys->mono_com),
            sys->monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        /* MONOLIS owns only internal X entries; synchronize overlap copies
         * before interpreting X as nodal Newton increments. */
        BBFE_fluid_mixed_sync_solution_vector(
            &(sys->mono_com), nnode, sys->monolis.mat.R.X);

        if(!mixed_vms_vector_is_finite_global(
            sys->monolis.mat.R.X,
            ndof,
            &(sys->mono_com),
            "VMS Newton correction")) {
            free(r0);
            return -1;
        }

        BBFE_fluid_sups_renew_velocity(
            sys->vals.delta_v,
            sys->monolis.mat.R.X,
            nnode);

        BBFE_fluid_sups_renew_pressure(
            sys->vals.delta_p,
            sys->monolis.mat.R.X,
            nnode);

        if(use_projected) {
            mixed_wall_project_velocity_increment(sys->vals.delta_v);
        }

        BBFE_fluid_mixed_report_newton_diagnostics(sys, step, it);
        BBFE_fluid_mixed_report_prism_update(sys, step, it);

        update_velocity_pressure_NR(
            sys->vals.v,
            sys->vals.delta_v,
            sys->vals.p,
            sys->vals.delta_p,
            nnode);

        if(use_projected) {
            mixed_wall_project_velocity_state(sys->vals.v,Uw);
        }

        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        BBFE_fluid_mixed_report_prism_components(sys, step, it);

        /* Updated-state residual.  Match solver_fom_VMS(): only residual
         * volume kernels are needed; the wall helper currently also adds a
         * matrix, which is harmless because the matrix is discarded. */
        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i=0; i<ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        mixed_vms_assemble_residual_only(sys);

        PRI6NitscheDiagnostics eval_pri_diag;
        PRI6NitscheDiagnostics eval_tet_diag;
        memset(&eval_pri_diag, 0, sizeof(eval_pri_diag));
        memset(&eval_tet_diag, 0, sizeof(eval_tet_diag));

        if(assemble_pri_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
                &(sys->monolis),
                &(sys->fe_pri),
                &(sys->basis_pri),
                surf_pri,
                basis_surf_pri,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_pri,
                &eval_pri_diag) != 0) {
                free(r0);
                return -1;
            }
        }

        if(assemble_tet_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
                &(sys->monolis),
                &(sys->fe_tet),
                &(sys->basis_tet),
                surf_tet,
                basis_surf_tet,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_tet,
                &eval_tet_diag) != 0) {
                free(r0);
                return -1;
            }
        }

        if(use_projected) {
            mixed_wall_project_velocity_rhs(sys->monolis.mat.R.B);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        if(!mixed_vms_vector_is_finite_global(
            sys->monolis.mat.R.B,
            ndof,
            &(sys->mono_com),
            "VMS updated-state residual")) {
            free(r0);
            return -1;
        }

        double n0 = 0.0;
        double nr = 0.0;

        for(int i=0; i<4*nint; ++i) {
            n0 += r0[i]*r0[i];
            nr += sys->monolis.mat.R.B[i]*sys->monolis.mat.R.B[i];
        }

        if(comm_size != 1) {
            monolis_allreduce_R(1, &n0, MONOLIS_MPI_SUM, sys->mono_com.comm);
            monolis_allreduce_R(1, &nr, MONOLIS_MPI_SUM, sys->mono_com.comm);
        }

        const double nrm0 = sqrt(fmax(n0, 0.0));
        const double nrmr = sqrt(fmax(nr, 0.0));
        const double rel = nrmr/fmax(nrm0, tiny);

        PRI6NitscheDiagnostics global_pri_diag;
        PRI6NitscheDiagnostics global_tet_diag;
        memset(&global_pri_diag, 0, sizeof(global_pri_diag));
        memset(&global_tet_diag, 0, sizeof(global_tet_diag));
        if(assemble_pri_wall) {
            pri6_reduce_nitsche_diagnostics(
                &eval_pri_diag, &global_pri_diag, &(sys->mono_com));
        }
        if(assemble_tet_wall) {
            pri6_reduce_nitsche_diagnostics(
                &eval_tet_diag, &global_tet_diag, &(sys->mono_com));
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "[mixed-VMS-split-nitsche] residual nrm=%e old=%e rel=%e "
                "PRI6=%s TET4=%s\n",
                nrmr, nrm0, rel,
                use_pri_wall ? "ON" : "OFF",
                use_tet_wall ? "ON" : "OFF");

            if(assemble_pri_wall) {
                printf(
                    "[PRI6-wall updated-state] copy_area=%e h=[%e,%e] "
                    "beta=[%e mean %e max %e] "
                    "beta_mu=[%e mean %e max %e] "
                    "beta_dt=[%e mean %e max %e] dt/mu(mean/max)=%e/%e "
                    "mapped=%lld skipped=%lld max|u.n|=%e\n",
                    global_pri_diag.area,
                    global_pri_diag.h_n_min,
                    global_pri_diag.h_n_max,
                    global_pri_diag.beta_min,
                    global_pri_diag.beta_mean,
                    global_pri_diag.beta_max,
                    global_pri_diag.beta_mu_min,
                    global_pri_diag.beta_mu_mean,
                    global_pri_diag.beta_mu_max,
                    global_pri_diag.beta_dt_min,
                    global_pri_diag.beta_dt_mean,
                    global_pri_diag.beta_dt_max,
                    global_pri_diag.beta_dt_over_beta_mu_mean,
                    global_pri_diag.beta_dt_over_beta_mu_max,
                    global_pri_diag.mapped_faces,
                    global_pri_diag.skipped_no_local_owner,
                    global_pri_diag.max_abs_normal_residual);
            }
            else {
                printf("[PRI6-wall updated-state] surface operator OFF (mode=%s)\n",
                       use_pri_wall ? BBFE_pri6_wall_mode_name(wall_opt_pri->wall_mode) : "off");
            }

            if(assemble_tet_wall) {
                printf(
                    "[TET4-wall updated-state] copy_area=%e h=[%e,%e] "
                    "beta=[%e mean %e max %e] "
                    "beta_mu=[%e mean %e max %e] "
                    "beta_dt=[%e mean %e max %e] dt/mu(mean/max)=%e/%e "
                    "mapped=%lld skipped=%lld max|u.n|=%e\n",
                    global_tet_diag.area,
                    global_tet_diag.h_n_min,
                    global_tet_diag.h_n_max,
                    global_tet_diag.beta_min,
                    global_tet_diag.beta_mean,
                    global_tet_diag.beta_max,
                    global_tet_diag.beta_mu_min,
                    global_tet_diag.beta_mu_mean,
                    global_tet_diag.beta_mu_max,
                    global_tet_diag.beta_dt_min,
                    global_tet_diag.beta_dt_mean,
                    global_tet_diag.beta_dt_max,
                    global_tet_diag.beta_dt_over_beta_mu_mean,
                    global_tet_diag.beta_dt_over_beta_mu_max,
                    global_tet_diag.mapped_faces,
                    global_tet_diag.skipped_no_local_owner,
                    global_tet_diag.max_abs_normal_residual);
            }
            else {
                printf("[TET4-wall updated-state] surface operator OFF (mode=%s)\n",
                       use_tet_wall ? BBFE_pri6_wall_mode_name(wall_opt_tet->wall_mode) : "off");
            }
        }

        mixed_vms_append_solver_history(
            sys,t,step,it,nrmr,nrm0,rel,(rel < eps_NR),
            assemble_pri_wall ? &global_pri_diag : NULL,
            assemble_tet_wall ? &global_tet_diag : NULL);

        if(rel < eps_NR) {
            converged = 1;
            ROM_BB_vec_copy_2d(
                sys->vals.v,
                sys->vals.v_old,
                nnode,
                3);

            if(monolis_mpi_get_global_my_rank() == 0) {
                printf(
                    "%s mixed VMS Newton converged: step=%d it=%d rel=%e\n",
                    CODENAME, step, it, rel);
            }
            break;
        }
    }

    free(r0);

    if(!converged) {
        if(monolis_mpi_get_global_my_rank() == 0) {
            fprintf(stderr,
                "%s ERROR: mixed VMS Newton did not converge at step=%d "
                "within %d iterations; v_old is NOT advanced.\n",
                CODENAME, step, max_iter_NR);
        }
        return -1;
    }

    return 0;
}

void solver_fom_NR_mixed_dual_nitsche(
    FE_SYSTEM_FLUID_MIXED* sys,
    const BBFE_DATA* surf_pri,
    const BBFE_BASIS* basis_surf_pri,
    const BBFE_DATA* surf_tet,
    const BBFE_BASIS* basis_surf_tet,
    int enable_pri6_nitsche,
    const PRI6WallOptions* wall_opt_pri,
    int enable_tet4_nitsche,
    const PRI6WallOptions* wall_opt_tet,
    double t,
    int step)
{
    (void)t;

    if(sys == NULL ||
       surf_pri == NULL || basis_surf_pri == NULL ||
       surf_tet == NULL || basis_surf_tet == NULL ||
       (enable_pri6_nitsche &&
        (wall_opt_pri == NULL || !pri6_validate_wall_options(wall_opt_pri))) ||
       (enable_tet4_nitsche &&
        (wall_opt_tet == NULL || !pri6_validate_wall_options(wall_opt_tet)))) {
        fprintf(stderr,
            "%s ERROR: invalid mixed split-Nitsche solver input.\n",
            CODENAME);
        exit(EXIT_FAILURE);
    }

    const int nnode = sys->fe_tet.total_num_nodes;
    const int ndof = 4*nnode;
    const int nint = sys->mono_com.n_internal_vertex;
    const int max_iter_NR = 20;
    const double eps_NR = 1.0e-6;
    const double tiny = 1.0e-30;
    const double Uw[3] = {0.0, 0.0, 0.0};

    const int use_pri_wall = (enable_pri6_nitsche != 0);
    const int use_tet_wall = (enable_tet4_nitsche != 0);
    const int pri_strong_noslip = use_pri_wall &&
        BBFE_pri6_wall_mode_is_strong_noslip(wall_opt_pri->wall_mode);
    const int tet_strong_noslip = use_tet_wall &&
        BBFE_pri6_wall_mode_is_strong_noslip(wall_opt_tet->wall_mode);
    const int pri_projected = use_pri_wall &&
        BBFE_pri6_wall_mode_uses_projected_strong_normal(wall_opt_pri->wall_mode);
    const int tet_projected = use_tet_wall &&
        BBFE_pri6_wall_mode_uses_projected_strong_normal(wall_opt_tet->wall_mode);
    const int assemble_pri_wall = use_pri_wall &&
        mixed_wall_mode_uses_surface_operator(wall_opt_pri->wall_mode);
    const int assemble_tet_wall = use_tet_wall &&
        mixed_wall_mode_uses_surface_operator(wall_opt_tet->wall_mode);
    const int use_projected = pri_projected || tet_projected;

    if(pri_strong_noslip &&
       mixed_wall_apply_strong_noslip_family(sys,surf_pri,Uw) != 0) exit(EXIT_FAILURE);
    if(tet_strong_noslip &&
       mixed_wall_apply_strong_noslip_family(sys,surf_tet,Uw) != 0) exit(EXIT_FAILURE);
    if(use_projected) {
        if(mixed_wall_projector_prepare(
            surf_pri,surf_tet,sys,pri_projected,tet_projected) != 0) {
            fprintf(stderr,
                "%s ERROR: failed to prepare projected-strong wall projector.\n",
                CODENAME);
            exit(EXIT_FAILURE);
        }
        mixed_wall_project_velocity_state(sys->vals.v,Uw);
        mixed_wall_project_velocity_state(sys->vals.v_old,Uw);
        BBFE_fluid_mixed_sync_state(&(sys->mono_com),&(sys->vals),nnode);
    }

    double* r0 = (double*)calloc((size_t)ndof, sizeof(double));
    if(r0 == NULL) exit(EXIT_FAILURE);

    for(int it=0; it<max_iter_NR; ++it) {
        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "\n%s ---- mixed runtime wall PRI6=%s(%s) TET4=%s(%s) step %d NR %d ----\n",
                CODENAME,
                use_pri_wall ? "ON" : "OFF",
                use_pri_wall ? BBFE_pri6_wall_mode_name(wall_opt_pri->wall_mode) : "off",
                use_tet_wall ? "ON" : "OFF",
                use_tet_wall ? BBFE_pri6_wall_mode_name(wall_opt_tet->wall_mode) : "off",
                step, it);
        }

        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i=0; i<ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        PRI6VolumeAssemblyDiagnostics volume_diag;
        pri6_assemble_volume_sups_mixed(sys, &volume_diag);

        PRI6NitscheDiagnostics local_pri_diag;
        PRI6NitscheDiagnostics local_tet_diag;
        memset(&local_pri_diag, 0, sizeof(local_pri_diag));
        memset(&local_tet_diag, 0, sizeof(local_tet_diag));

        if(assemble_pri_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
                &(sys->monolis),
                &(sys->fe_pri),
                &(sys->basis_pri),
                surf_pri,
                basis_surf_pri,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_pri,
                &local_pri_diag) != 0) {

                fprintf(stderr,
                    "%s ERROR: PRI6/TRI3 Nitsche assembly failed.\n",
                    CODENAME);
                free(r0);
                exit(EXIT_FAILURE);
            }
        }

        if(assemble_tet_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
                &(sys->monolis),
                &(sys->fe_tet),
                &(sys->basis_tet),
                surf_tet,
                basis_surf_tet,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_tet,
                &local_tet_diag) != 0) {

                fprintf(stderr,
                    "%s ERROR: TET4/TRI3 Nitsche assembly failed.\n",
                    CODENAME);
                free(r0);
                exit(EXIT_FAILURE);
            }
        }

        int global_pri_elems = volume_diag.local_prism_elements;
        double global_rhs_l1[4];
        double global_rhs_linf[4];
        for(int d=0; d<4; ++d) {
            global_rhs_l1[d] = volume_diag.rhs_l1[d];
            global_rhs_linf[d] = volume_diag.rhs_linf[d];
        }

        monolis_allreduce_I(
            1, &global_pri_elems, MONOLIS_MPI_SUM, sys->mono_com.comm);
        for(int d=0; d<4; ++d) {
            monolis_allreduce_R(
                1, &global_rhs_l1[d], MONOLIS_MPI_SUM, sys->mono_com.comm);
            monolis_allreduce_R(
                1, &global_rhs_linf[d], MONOLIS_MPI_MAX, sys->mono_com.comm);
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            printf(
                "[PRI6 SUPS step=%d NR=%d] elems=%d\n"
                "  ux: L1=%e Linf=%e\n"
                "  uy: L1=%e Linf=%e\n"
                "  uz: L1=%e Linf=%e\n"
                "   p: L1=%e Linf=%e\n",
                step, it, global_pri_elems,
                global_rhs_l1[0], global_rhs_linf[0],
                global_rhs_l1[1], global_rhs_linf[1],
                global_rhs_l1[2], global_rhs_linf[2],
                global_rhs_l1[3], global_rhs_linf[3]);
        }

        /* Ahmed wall DOFs are intentionally absent from D_bc_v.dat. */
        if(use_projected) {
            mixed_wall_project_velocity_rhs(sys->monolis.mat.R.B);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        if(it == 0) {
            memcpy(r0, sys->monolis.mat.R.B, sizeof(double)*(size_t)ndof);
        }

        BBFE_sys_monowrap_solve(
            &(sys->monolis),
            &(sys->mono_com),
            sys->monolis.mat.R.X,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        BBFE_fluid_mixed_sync_solution_vector(
            &(sys->mono_com), nnode, sys->monolis.mat.R.X);

        BBFE_fluid_sups_renew_velocity(
            sys->vals.delta_v,
            sys->monolis.mat.R.X,
            nnode);

        BBFE_fluid_sups_renew_pressure(
            sys->vals.delta_p,
            sys->monolis.mat.R.X,
            nnode);

        if(use_projected) {
            mixed_wall_project_velocity_increment(sys->vals.delta_v);
        }

        BBFE_fluid_mixed_report_newton_diagnostics(sys, step, it);
        BBFE_fluid_mixed_report_prism_update(sys, step, it);

        update_velocity_pressure_NR(
            sys->vals.v,
            sys->vals.delta_v,
            sys->vals.p,
            sys->vals.delta_p,
            nnode);

        if(use_projected) {
            mixed_wall_project_velocity_state(sys->vals.v,Uw);
        }

        BBFE_fluid_mixed_sync_state(
            &(sys->mono_com), &(sys->vals), nnode);

        BBFE_fluid_mixed_report_prism_components(sys, step, it);

        /* Residual re-evaluation at the updated state. */
        monolis_clear_mat_value_R(&(sys->monolis));
        for(int i=0; i<ndof; ++i) {
            sys->monolis.mat.R.B[i] = 0.0;
            sys->monolis.mat.R.X[i] = 0.0;
        }

        pri6_assemble_volume_sups_mixed(sys, NULL);

        if(assemble_pri_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_pri6_tri3(
                &(sys->monolis),
                &(sys->fe_pri),
                &(sys->basis_pri),
                surf_pri,
                basis_surf_pri,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_pri,
                NULL) != 0) {
                free(r0);
                exit(EXIT_FAILURE);
            }
        }

        if(assemble_tet_wall) {
            if(set_wall_face_vecmat_symmetric_nitsche_tet4_tri3(
                &(sys->monolis),
                &(sys->fe_tet),
                &(sys->basis_tet),
                surf_tet,
                basis_surf_tet,
                &(sys->vals),
                sys->vals.density,
                sys->vals.viscosity,
                Uw,
                wall_opt_tet,
                NULL) != 0) {
                free(r0);
                exit(EXIT_FAILURE);
            }
        }

        if(use_projected) {
            mixed_wall_project_velocity_rhs(sys->monolis.mat.R.B);
        }

        BBFE_sys_monowrap_set_Dirichlet_bc(
            &(sys->monolis),
            nnode,
            4,
            &(sys->bc_NR),
            sys->monolis.mat.R.B);

        double n0 = 0.0;
        double nr = 0.0;
        for(int i=0; i<4*nint; ++i) {
            n0 += r0[i]*r0[i];
            nr += sys->monolis.mat.R.B[i]*sys->monolis.mat.R.B[i];
        }

        monolis_allreduce_R(1, &n0, MONOLIS_MPI_SUM, sys->mono_com.comm);
        monolis_allreduce_R(1, &nr, MONOLIS_MPI_SUM, sys->mono_com.comm);
        const double rel = sqrt(nr)/fmax(sqrt(n0), tiny);

        PRI6NitscheDiagnostics global_pri_diag;
        PRI6NitscheDiagnostics global_tet_diag;
        memset(&global_pri_diag, 0, sizeof(global_pri_diag));
        memset(&global_tet_diag, 0, sizeof(global_tet_diag));
        if(assemble_pri_wall) {
            pri6_reduce_nitsche_diagnostics(
                &local_pri_diag, &global_pri_diag, &(sys->mono_com));
        }
        if(assemble_tet_wall) {
            pri6_reduce_nitsche_diagnostics(
                &local_tet_diag, &global_tet_diag, &(sys->mono_com));
        }

        if(monolis_mpi_get_global_my_rank() == 0) {
            if(assemble_pri_wall) {
                printf(
                    "[PRI6-wall Nitsche assembly] copy_area=%e h=[%e,%e] beta=[%e,%e] "
                    "mapped_copies=%lld skipped=%lld |A|_1(copy-sum)=%e |R|_1(copy-sum)=%e max|u.n|=%e\n",
                    global_pri_diag.area,
                    global_pri_diag.h_n_min,
                    global_pri_diag.h_n_max,
                    global_pri_diag.beta_min,
                    global_pri_diag.beta_max,
                    global_pri_diag.mapped_faces,
                    global_pri_diag.skipped_no_local_owner,
                    global_pri_diag.matrix_abs_sum,
                    global_pri_diag.rhs_abs_sum,
                    global_pri_diag.max_abs_normal_residual);
            }
            else {
                printf("[PRI6-wall] surface operator OFF (mode=%s)\n",
                       use_pri_wall ? BBFE_pri6_wall_mode_name(wall_opt_pri->wall_mode) : "off");
            }

            if(assemble_tet_wall) {
                printf(
                    "[TET4-wall Nitsche assembly] copy_area=%e h=[%e,%e] beta=[%e,%e] "
                    "mapped_copies=%lld skipped=%lld |A|_1(copy-sum)=%e |R|_1(copy-sum)=%e max|u.n|=%e\n",
                    global_tet_diag.area,
                    global_tet_diag.h_n_min,
                    global_tet_diag.h_n_max,
                    global_tet_diag.beta_min,
                    global_tet_diag.beta_max,
                    global_tet_diag.mapped_faces,
                    global_tet_diag.skipped_no_local_owner,
                    global_tet_diag.matrix_abs_sum,
                    global_tet_diag.rhs_abs_sum,
                    global_tet_diag.max_abs_normal_residual);
            }
            else {
                printf("[TET4-wall] surface operator OFF (mode=%s)\n",
                       use_tet_wall ? BBFE_pri6_wall_mode_name(wall_opt_tet->wall_mode) : "off");
            }

            printf(
                "[mixed-dual-nitsche] residual rel=%e mapped_total=%lld "
                "skipped_total=%lld max|u.n|=%e\n",
                rel,
                global_pri_diag.mapped_faces + global_tet_diag.mapped_faces,
                global_pri_diag.skipped_no_local_owner +
                    global_tet_diag.skipped_no_local_owner,
                fmax(global_pri_diag.max_abs_normal_residual,
                     global_tet_diag.max_abs_normal_residual));
        }

        if(rel < eps_NR) break;
    }

    ROM_BB_vec_copy_2d(sys->vals.v, sys->vals.v_old, nnode, 3);
    free(r0);
}

/* PRI6 Karman force diagnostics: keep this include at EOF. */
#include "core_Nitsche_prism_diagnostics.inc"
