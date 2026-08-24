#include "core_Nitsche_prism.h"

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

static void pri6_get_penalty(
    const PRI6WallOptions* opt,
    double rho,
    double mu_eff,
    double h_n,
    double dt,
    double* beta)
{
    const double h = fmax(h_n, 1.0e-12);
    const double dte = fmax(dt, 1.0e-30);
    *beta = opt->gamma_n * mu_eff / h
          + opt->gamma_n * opt->dt_penalty_coeff * rho * h / dte;
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
            double beta=0.0; pri6_get_penalty(opt,rho,mu_eff,h_n,vals->dt,&beta);

            if(diag) {
                diag->area+=w;
                diag->h_n_min=fmin(diag->h_n_min,h_n);
                diag->h_n_max=fmax(diag->h_n_max,h_n);
                diag->beta_min=fmin(diag->beta_min,beta);
                diag->beta_max=fmax(diag->beta_max,beta);
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
    }

    free(owner);
    free(lface);
    return 0;
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

    monolis_allreduce_R(1, &area,   MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &rhs_l1, MONOLIS_MPI_SUM, mono_com->comm);
    monolis_allreduce_R(1, &mat_l1, MONOLIS_MPI_SUM, mono_com->comm);

    global->area = area;
    global->rhs_abs_sum = rhs_l1;
    global->matrix_abs_sum = mat_l1;

    double hmax = local->h_n_max;
    double bmax = local->beta_max;
    double unmax = local->max_abs_normal_residual;

    monolis_allreduce_R(1, &hmax,  MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &bmax,  MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &unmax, MONOLIS_MPI_MAX, mono_com->comm);

    global->h_n_max = hmax;
    global->beta_max = bmax;
    global->max_abs_normal_residual = unmax;

    /* No MIN reduction wrapper is assumed; reduce negative minima by MAX. */
    double neg_hmin = (local->assembled_faces > 0) ? -local->h_n_min : -DBL_MAX;
    double neg_bmin = (local->assembled_faces > 0) ? -local->beta_min : -DBL_MAX;

    monolis_allreduce_R(1, &neg_hmin, MONOLIS_MPI_MAX, mono_com->comm);
    monolis_allreduce_R(1, &neg_bmin, MONOLIS_MPI_MAX, mono_com->comm);

    global->h_n_min = (neg_hmin > -DBL_MAX) ? -neg_hmin : 0.0;
    global->beta_min = (neg_bmin > -DBL_MAX) ? -neg_bmin : 0.0;
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
                "[mixed-pri6-nitsche] residual rel=%e area=%e "
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

/* PRI6 Karman force diagnostics: keep this include at EOF. */
#include "core_Nitsche_prism_diagnostics.inc"
