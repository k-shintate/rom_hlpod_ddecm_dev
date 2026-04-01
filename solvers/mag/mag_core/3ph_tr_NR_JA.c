#include "3ph_tr_NR.h"
#include "3ph_tr_NR_JA.h"

/* ============================================================
 * Voltage excitation settings
 * ============================================================ */
static const double V_AMP  = 100.0;   /* [V] phase voltage amplitude (temporary) */
static const double R_COIL = 0.1;     /* [Ohm] equivalent phase winding resistance */

/* 各相電流（時間ステップ更新） */
double g_phase_current[3]      = {0.0, 0.0, 0.0};
double g_phase_current_prev[3] = {0.0, 0.0, 0.0};

/* ============================================================
 * Applied phase voltage
 * ============================================================ */
static inline double get_phase_voltage(int phase_idx, double t)
{
    double omega = 2.0 * M_PI * FREQ_HZ;

    if(phase_idx == 0) return V_AMP * sin(omega * t);
    if(phase_idx == 1) return V_AMP * sin(omega * t - 2.0 * M_PI / 3.0);
    if(phase_idx == 2) return V_AMP * sin(omega * t + 2.0 * M_PI / 3.0);

    return 0.0;
}

/* ============================================================
 * Jiles-Atherton parameter / state
 * ============================================================ */

typedef struct {
    double Ms;       /* saturation magnetization */
    double a;        /* Langevin parameter */
    double k;        /* pinning parameter */
    double c;        /* reversible coefficient */
    double alpha_j;  /* inter-domain coupling */
} JA_Param;

typedef struct {
    double b_prev[3];      /* accepted B at time n */
    double h_prev[3];      /* accepted H at time n */

    double b_iter_in[3];   /* input B for current NR material update */
    double h_iter_out[3];  /* output H for current NR material update */

    double dh_db[3][3];    /* tangent reluctivity tensor dH/dB */
} JA_IPState;

/* ============================================================
 * Global / module-level JA storage
 * One state per (element, integration point)
 * ============================================================ */

static JA_Param g_ja_param = {
    .Ms      = 1.6e6,   /* placeholder */
    .a       = 400.0,   /* placeholder */
    .k       = 250.0,   /* placeholder */
    .c       = 0.2,     /* placeholder */
    .alpha_j = 1.0e-4   /* placeholder */
};

static JA_IPState** g_ja_state = NULL;   /* [elem][ip] */
static int g_ja_num_elem = 0;
static int g_ja_num_ip   = 0;

/* material update substeps inside one FE evaluation */
static const int JA_SUBSTEP = 10;

/* ============================================================
 * Material coefficient policy (prop: 1..3 coil, 4 core, 5 air)
 * ============================================================ */

const double Sigma_coil = 5.77e7;                  // 無次元化（典型導電率）
const double Sigma_core =  3.72e3;                  // 無次元化（典型導電率）


void get_sigmas_for_prop_JA(
    int prop,
    double* sigma_mass_A,   /* used for (sigma/dt) * M_A in Jacobian */
    double* sigma_cpl,      /* used for coupling C and C^T/dt */
    double* sigma_phi       /* used for phi-phi Laplace term */
){
    if(prop == 1 || prop == 2 || prop == 3){
        /* coil (conductor) */
        *sigma_mass_A = Sigma_coil;
        *sigma_cpl    = Sigma_coil;
        *sigma_phi    = Sigma_coil;
    } else if(prop == 4){
        /* core (A-only by default; phi is not solved here) */
        *sigma_mass_A = Sigma_core;
        *sigma_cpl    = Sigma_core;
        *sigma_phi    = Sigma_core;
    } else if(prop == 5){
        /* air: coupling cut; add small penalties for solvability */
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    } else {
        *sigma_mass_A = 0.0;
        *sigma_cpl    = 0.0;
        *sigma_phi    = 0.0;
    }
}


/* ============================================================
 * Small vector / matrix helpers
 * ============================================================ */

static inline void vec3_zero(double v[3]){
    v[0] = v[1] = v[2] = 0.0;
}

static inline void vec3_copy(double dst[3], const double src[3]){
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}

static inline void vec3_add(double out[3], const double a[3], const double b[3]){
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
}

static inline void vec3_sub(double out[3], const double a[3], const double b[3]){
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

static inline void vec3_scale(double out[3], const double a[3], double s){
    out[0] = a[0] * s;
    out[1] = a[1] * s;
    out[2] = a[2] * s;
}

static inline double vec3_dot(const double a[3], const double b[3]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double vec3_norm(const double a[3]){
    return sqrt(vec3_dot(a, a));
}

static inline void vec3_axpy(double y[3], double alpha, const double x[3]){
    y[0] += alpha * x[0];
    y[1] += alpha * x[1];
    y[2] += alpha * x[2];
}

static inline void mat3_eye(double A[3][3]){
    memset(A, 0, sizeof(double)*9);
    A[0][0] = A[1][1] = A[2][2] = 1.0;
}

static inline void mat3_copy(double dst[3][3], const double src[3][3]){
    memcpy(dst, src, sizeof(double)*9);
}

static inline void mat3_vec(double y[3], const double A[3][3], const double x[3]){
    y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
    y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
    y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

static inline void mat3_add(double C[3][3], const double A[3][3], const double B[3][3]){
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

static inline void mat3_sub(double C[3][3], const double A[3][3], const double B[3][3]){
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

static inline void mat3_scale(double C[3][3], const double A[3][3], double s){
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            C[i][j] = A[i][j] * s;
        }
    }
}

static inline void mat3_outer(double A[3][3], const double u[3], const double v[3]){
    for(int i=0;i<3;++i){
        for(int j=0;j<3;++j){
            A[i][j] = u[i]*v[j];
        }
    }
}

static inline double mat3_det(const double A[3][3]){
    return
        A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
        A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
        A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

static int mat3_inv(double Ainv[3][3], const double A[3][3]){
    const double det = mat3_det(A);
    const double eps = 1.0e-20;

    if(fabs(det) < eps){
        return 0;
    }

    Ainv[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / det;
    Ainv[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) / det;
    Ainv[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / det;

    Ainv[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) / det;
    Ainv[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / det;
    Ainv[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) / det;

    Ainv[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / det;
    Ainv[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) / det;
    Ainv[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / det;

    return 1;
}

/* ============================================================
 * JA state allocation / initialization / commit
 * ============================================================ */

void ja_allocate_states(int n_elem, int n_ip)
{
    g_ja_num_elem = n_elem;
    g_ja_num_ip   = n_ip;

    g_ja_state = (JA_IPState**)calloc(n_elem, sizeof(JA_IPState*));
    for(int e=0; e<n_elem; ++e){
        g_ja_state[e] = (JA_IPState*)calloc(n_ip, sizeof(JA_IPState));
    }
}

void ja_free_states(void)
{
    if(g_ja_state == NULL) return;
    for(int e=0; e<g_ja_num_elem; ++e){
        free(g_ja_state[e]);
    }
    free(g_ja_state);
    g_ja_state = NULL;
    g_ja_num_elem = 0;
    g_ja_num_ip   = 0;
}

void ja_init_states_zero(void)
{
    if(g_ja_state == NULL) return;

    for(int e=0; e<g_ja_num_elem; ++e){
        for(int p=0; p<g_ja_num_ip; ++p){
            vec3_zero(g_ja_state[e][p].b_prev);
            vec3_zero(g_ja_state[e][p].h_prev);
            vec3_zero(g_ja_state[e][p].b_iter_in);
            vec3_zero(g_ja_state[e][p].h_iter_out);
            mat3_eye(g_ja_state[e][p].dh_db);
        }
    }
}

void ja_commit_step_states(BBFE_DATA* fe, BBFE_BASIS* basis, NEDELEC* ned, const double* x_curr)
{
    const int np = basis->num_integ_points;
    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){
        if(ned->elem_prop[e] != 4) continue; /* core only */

        for(int p=0; p<np; ++p){
            JA_IPState* s = &g_ja_state[e][p];

            double B[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B);

            vec3_copy(s->b_prev, B);
            vec3_copy(s->h_prev, s->h_iter_out);
        }
    }
}

/* ============================================================
 * Langevin / anhysteretic utilities
 * Based on the paper's inverse JA formulation
 * ============================================================ */

static double coth_safe(double x)
{
    const double ax = fabs(x);
    if(ax < 1.0e-8){
        /* coth(x) ~ 1/x + x/3 */
        return 1.0/x + x/3.0;
    }
    return cosh(x)/sinh(x);
}

/* scalar anhysteretic magnitude */
static double ja_man_scalar(double he_mag, const JA_Param* prm)
{
    const double a = prm->a;
    const double Ms = prm->Ms;
    const double eps = 1.0e-14;

    const double x = he_mag / fmax(a, eps);

    if(fabs(x) < 1.0e-8){
        /* Langevin small-x expansion: coth(x) - 1/x ~ x/3 */
        return Ms * (x / 3.0);
    }

    return Ms * (coth_safe(x) - 1.0/x);
}

static double ja_dman_dhe_scalar(double he_mag, const JA_Param* prm)
{
    const double a = prm->a;
    const double Ms = prm->Ms;
    const double eps = 1.0e-14;

    const double x = he_mag / fmax(a, eps);

    if(fabs(x) < 1.0e-6){
        /* derivative of Ms * L(x), L'(0)=1/3 */
        return Ms / (3.0 * fmax(a, eps));
    }

    /* d/dhe [Ms(coth(x)-1/x)] with x=he/a */
    const double csch = 1.0 / sinh(x);
    return Ms * ( (-csch*csch) + 1.0/(x*x) ) / fmax(a, eps);
}

/* vector anhysteretic magnetization and its tensor derivative dman/dhe */
static void ja_compute_man_and_tangent(
    const double he[3],
    const JA_Param* prm,
    double man[3],
    double dman_dhe[3][3]
){
    const double eps = 1.0e-14;
    const double he_mag = fmax(vec3_norm(he), eps);

    const double man_mag = ja_man_scalar(he_mag, prm);
    const double dman_mag = ja_dman_dhe_scalar(he_mag, prm);

    double u[3] = { he[0]/he_mag, he[1]/he_mag, he[2]/he_mag };

    man[0] = man_mag * u[0];
    man[1] = man_mag * u[1];
    man[2] = man_mag * u[2];

    /* d(man)/d(he) = (man_mag/he_mag)(I - uu^T) + dman_mag (uu^T) */
    double II[3][3], uu[3][3], term1[3][3], term2[3][3], tmp[3][3];
    mat3_eye(II);
    mat3_outer(uu, u, u);
    mat3_sub(tmp, II, uu);
    mat3_scale(term1, tmp, man_mag / he_mag);
    mat3_scale(term2, uu, dman_mag);
    mat3_add(dman_dhe, term1, term2);
}

/* ============================================================
 * Approximate inverse vector JA material update
 *
 * Input:
 *   b_prev, h_prev : accepted state at previous time step
 *   b_curr         : trial B at current Newton iterate
 *
 * Output:
 *   h_curr         : updated H
 *   dh_db          : consistent-ish tangent dH/dB
 *
 * Notes:
 * - This is an implementation scaffold based on the paper.
 * - The paper stores the magnetic field per Gauss point and uses it
 *   during NR iterations. :contentReference[oaicite:3]{index=3}
 * - For robustness, we substep B from b_prev to b_curr.
 * ============================================================ */

void ja_inverse_update(
    const JA_Param* prm,
    const double b_prev[3],
    const double h_prev[3],
    const double b_curr[3],
    double h_curr[3],
    double dh_db[3][3]
){
    const double mu0 = 4.0 * M_PI * 1.0e-7;
    const double eps = 1.0e-14;

    double b_k[3], h_k[3];
    vec3_copy(b_k, b_prev);
    vec3_copy(h_k, h_prev);

    mat3_eye(dh_db);

    for(int sub=0; sub<JA_SUBSTEP; ++sub){
        const double t0 = (double)sub / (double)JA_SUBSTEP;
        const double t1 = (double)(sub+1) / (double)JA_SUBSTEP;

        double b0[3], b1[3], db[3];
        for(int i=0;i<3;++i){
            b0[i] = b_prev[i] + (b_curr[i] - b_prev[i]) * t0;
            b1[i] = b_prev[i] + (b_curr[i] - b_prev[i]) * t1;
            db[i] = b1[i] - b0[i];
        }

        /* recover previous total magnetization from b = mu0 (h + m) */
        double m_k[3];
        for(int i=0;i<3;++i){
            m_k[i] = b_k[i]/mu0 - h_k[i];
        }

        /* effective field */
        double he[3];
        for(int i=0;i<3;++i){
            he[i] = h_k[i] + prm->alpha_j * m_k[i];
        }

        double man[3], dman_dhe[3][3];
        ja_compute_man_and_tangent(he, prm, man, dman_dhe);

        /* irreversible direction delta = (man - mirr)/|man - mirr| */
        double mirr_k[3];
        {
            /* m = (1-c) mirr + c man => mirr = (m - c man)/(1-c) */
            const double denom = fmax(1.0 - prm->c, 1.0e-12);
            for(int i=0;i<3;++i){
                mirr_k[i] = (m_k[i] - prm->c * man[i]) / denom;
            }
        }

        double diff[3];
        vec3_sub(diff, man, mirr_k);

        double diff_mag = vec3_norm(diff);
        double delta_dir[3] = {0.0, 0.0, 0.0};
        if(diff_mag > eps){
            delta_dir[0] = diff[0]/diff_mag;
            delta_dir[1] = diff[1]/diff_mag;
            delta_dir[2] = diff[2]/diff_mag;
        }

        /* local linearization:
           db = mu0 (dh + dm)
           dm = ((1-c) dmirr/dhe + c dman/dhe) dhe
           dhe = dh + alpha_j dm

           We approximate dmirr/dhe using the paper's directional rule.
        */
        double proj = 0.0;
        {
            /* crude estimate of dhe direction from previous tangent step */
            proj = vec3_dot(delta_dir, db);
        }

        double dmirr_dhe[3][3];
        memset(dmirr_dhe, 0, sizeof(dmirr_dhe));

        if(proj > 0.0 && diff_mag > eps){
            double outer_dd[3][3];
            mat3_outer(outer_dd, delta_dir, delta_dir);
            /* simplified version of directional irreversible slope */
            const double slope = diff_mag / fmax(prm->k, eps);
            mat3_scale(dmirr_dhe, outer_dd, slope);
        }

        /* A = dman/dhe, G = (1-c) dmirr/dhe + c A */
        double G[3][3];
        {
            double t1_[3][3], t2_[3][3];
            mat3_scale(t1_, dmirr_dhe, (1.0 - prm->c));
            mat3_scale(t2_, dman_dhe, prm->c);
            mat3_add(G, t1_, t2_);
        }

        /* dm = G dhe, dhe = dh + alpha_j dm = dh + alpha_j G dhe
           => (I - alpha_j G) dhe = dh
           => dm = G (I - alpha_j G)^(-1) dh
           => db = mu0 (I + G (I - alpha_j G)^(-1)) dh
        */
        double II[3][3], alphaG[3][3], M1[3][3], M1_inv[3][3];
        mat3_eye(II);
        mat3_scale(alphaG, G, prm->alpha_j);
        mat3_sub(M1, II, alphaG);

        if(!mat3_inv(M1_inv, M1)){
            mat3_eye(dh_db);
            vec3_copy(h_curr, h_k);
            return;
        }

        double Tm[3][3]; /* dm/dh */
        {
            double col[3];
            for(int j=0;j<3;++j){
                double ej[3] = {0.0,0.0,0.0};
                ej[j] = 1.0;
                double z[3], gcol[3];
                mat3_vec(z, M1_inv, ej);
                mat3_vec(gcol, G, z);
                Tm[0][j] = gcol[0];
                Tm[1][j] = gcol[1];
                Tm[2][j] = gcol[2];
            }
        }

        double db_dh[3][3];
        {
            double tmp[3][3];
            mat3_add(tmp, II, Tm);
            mat3_scale(db_dh, tmp, mu0);
        }

        if(!mat3_inv(dh_db, db_dh)){
            mat3_eye(dh_db);
            vec3_copy(h_curr, h_k);
            return;
        }

        /* update dh = (dH/dB) dB */
        double dh[3];
        mat3_vec(dh, dh_db, db);

        h_k[0] += dh[0];
        h_k[1] += dh[1];
        h_k[2] += dh[2];

        b_k[0] = b1[0];
        b_k[1] = b1[1];
        b_k[2] = b1[2];
    }

    vec3_copy(h_curr, h_k);
}

/* ============================================================
 * Residual-norm utility for line search
 * ============================================================ */

static double local_sum_sq(const double* v, int n) {
    double s = 0.0;
    for(int i=0; i<n; ++i) s += v[i]*v[i];
    return s;
}

static double local_max_abs(const double* v, int n) {
    double m = 0.0;
    for(int i=0; i<n; ++i){
        double a = fabs(v[i]);
        if(a > m) m = a;
    }
    return m;
}

/* ============================================================
 * Residual Assembly (JA version)
 *
 * Important:
 * Replace only the core magnetic constitutive part.
 * Keep your existing coil / phi / external circuit terms.
 * ============================================================ */
/* ============================================================
 * Residual Assembly (Newton) : B = -F(x)
 * JA inverse model in core region (prop == 4)
 * Voltage excitation version
 * ============================================================ */
void initialize_g_phase_current(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_R    = BB_std_calloc_1d_double(val_ip_R, np);

    const double inv_dt = 1.0 / dt;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    /* initial phase current guess */
    if(current_time == dt){
        printf("initialize g_phase_current\n");
        for(int k = 0; k < 3; k++){
            double Vapp = get_phase_voltage(k, current_time);
            g_phase_current[k] = Vapp / R_COIL;
        }
    }
}

void set_element_vec_NR_Aphi_JA(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_prev,   /* time n */
    const double* x_curr,   /* time n+1 (Newton iter) */
    double dt,
    double current_time     /* time n+1 */
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_R    = BB_std_calloc_1d_double(val_ip_R, np);

    const double inv_dt = 1.0 / dt;
    const double eps_r  = 1.0e-14;

    const int nEdge = fe->total_num_nodes;

    for(int e = 0; e < fe->total_num_elems; ++e){
        const int prop = ned->elem_prop[e];

        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_JA(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        /* coil info */
        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ==========================================================
         * [1] Source term (coil excitation)
         * ========================================================== */
        if(is_coil){
            int phase_idx = prop - 1; /* prop 1,2,3 -> U,V,W */
            double I_t = 0.0;
            if(0 <= phase_idx && phase_idx < 3){
                I_t = g_phase_current[phase_idx];
            }

            double J_mag = (coil.turns * I_t) / coil.area;

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int p = 0; p < np; ++p){
                    double x_ip[3];
                    get_interp_coords(e, p, fe, basis, x_ip);

                    double d[3] = {
                        x_ip[0] - coil.center[0],
                        x_ip[1] - coil.center[1],
                        x_ip[2] - coil.center[2]
                    };
                    double da = dot3(d, coil.axis);

                    double r[3] = {
                        d[0] - da * coil.axis[0],
                        d[1] - da * coil.axis[1],
                        d[2] - da * coil.axis[2]
                    };

                    double tdir[3] = {
                        coil.axis[1] * r[2] - coil.axis[2] * r[1],
                        coil.axis[2] * r[0] - coil.axis[0] * r[2],
                        coil.axis[0] * r[1] - coil.axis[1] * r[0]
                    };
                    double n_tdir = norm3(tdir);

                    double Js[3] = {0.0, 0.0, 0.0};
                    if(n_tdir > eps_r){
                        double inv_n = 1.0 / n_tdir;
                        Js[0] = J_mag * tdir[0] * inv_n;
                        Js[1] = J_mag * tdir[1] * inv_n;
                        Js[2] = J_mag * tdir[2] * inv_n;
                    }

                    val_ip_R[p] = dot3(Js, ned->N_edge[e][p][i]);
                }

                double integ = BBFE_std_integ_calc(
                    np, val_ip_R, basis->integ_weight, Jacobian_ip
                );
                monolis->mat.R.B[gi] += (double)si * integ;
            }
        }

        /* ==========================================================
         * [2] A-mass term
         * ========================================================== */
        if(sigma_mass_A > 0.0){
            int use_history = (prop == 1 || prop == 2 || prop == 3 || prop == 5);

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                double acc = 0.0;
                for(int j = 0; j < ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    double coeffA;
                    if(use_history){
                        coeffA = (x_curr[gj] - x_prev[gj]);
                    }else{
                        coeffA = (x_curr[gj] - x_prev[gj]);
                    }

                    for(int p = 0; p < np; ++p){
                        val_ip_R[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            ned->N_edge[e][p][j],
                            sigma_mass_A
                        );
                    }

                    double mij = BBFE_std_integ_calc(
                        np, val_ip_R, basis->integ_weight, Jacobian_ip
                    );
                    acc += (double)(si * sj) * mij * coeffA * inv_dt;
                }

                monolis->mat.R.B[gi] -= acc;
            }
        }

        /* ==========================================================
         * [3] Magnetic term in A-equation:
         *     \int curl(W_i) · H(B) dΩ
         *     core(prop==4): JA inverse
         *     others       : linear H = nu B
         * ========================================================== */
        for(int i = 0; i < ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int p = 0; p < np; ++p){
                const double* ci = ned->curl_N_edge[e][p][i];

                double B[3];
                double H[3];

                compute_B_ip(ned, e, p, x_curr, nEdge, B);

                if(prop == 4){
                    JA_IPState* s = &g_ja_state[e][p];

                    vec3_copy(s->b_iter_in, B);

                    ja_inverse_update(
                        &g_ja_param,
                        s->b_prev,      /* previous converged B */
                        s->h_prev,      /* previous converged H */
                        s->b_iter_in,   /* current Newton/FPI input B */
                        s->h_iter_out,  /* output H(B) */
                        s->dh_db        /* tangent if needed elsewhere */
                    );

                    vec3_copy(H, s->h_iter_out);
                }else{
                    H[0] = NU_LIN * B[0];
                    H[1] = NU_LIN * B[1];
                    H[2] = NU_LIN * B[2];
                }

                val_ip_R[p] = vec3_dot(ci, H);
            }

            double v = BBFE_std_integ_calc(
                np, val_ip_R, basis->integ_weight, Jacobian_ip
            );

            /* residual B = -F(x) */
            monolis->mat.R.B[gi] -= (double)si * v;
        }

        /* ==========================================================
         * [4] A-Phi coupling in A-equation
         * ========================================================== */
        if(sigma_cpl > 0.0){
            for(int j = 0; j < ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                double acc = 0.0;
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double phi_val = x_curr[gn];

                    for(int p = 0; p < np; ++p){
                        val_ip_R[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n],
                            ned->N_edge[e][p][j],
                            sigma_cpl
                        );
                    }

                    double cjn = BBFE_std_integ_calc(
                        np, val_ip_R, basis->integ_weight, Jacobian_ip
                    );
                    acc += (double)sj * cjn * phi_val;
                }

                monolis->mat.R.B[gj] -= acc;
            }
        }

        /* ==========================================================
         * [5] Phi-equation
         * ========================================================== */
        if(sigma_cpl > 0.0){
            for(int n = 0; n < fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];

                double acc = 0.0;
                for(int i = 0; i < ned->local_num_edges; ++i){
                    int gi = ned->nedelec_conn[e][i];
                    int si = ned->edge_sign[e][i];
                    double da = x_curr[gi] - x_prev[gi];

                    for(int p = 0; p < np; ++p){
                        val_ip_R[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i],
                            fe->geo[e][p].grad_N[n],
                            sigma_cpl
                        );
                    }

                    double gin = BBFE_std_integ_calc(
                        np, val_ip_R, basis->integ_weight, Jacobian_ip
                    );
                    acc += (double)si * gin * da * inv_dt;
                }

                monolis->mat.R.B[gn] -= acc;
            }
        }

        /* ==========================================================
         * [6] Phi Laplace term
         * ========================================================== */
        {
            double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;

            if(sigma_laplace > 0.0){
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];

                    double acc = 0.0;
                    for(int m = 0; m < fe->local_num_nodes; ++m){
                        int gm = fe->conn[e][m];
                        double phi_m = x_curr[gm];

                        for(int p = 0; p < np; ++p){
                            val_ip_R[p] = BBFE_elemmat_mag_mat_mass(
                                fe->geo[e][p].grad_N[n],
                                fe->geo[e][p].grad_N[m],
                                sigma_laplace
                            );
                        }

                        double knm = BBFE_std_integ_calc(
                            np, val_ip_R, basis->integ_weight, Jacobian_ip
                        );
                        acc += knm * phi_m;
                    }

                    monolis->mat.R.B[gn] -= acc;
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_R, np);
}
/* ============================================================
 * Jacobian Assembly (JA version)
 *
 * This replaces the core tangent:
 *   nu * ci.cj + alpha (ci.B)(cj.B)
 * with:
 *   ci^T (dH/dB) cj
 *
 * Your current Jacobian uses the former non-hysteretic law. :contentReference[oaicite:4]{index=4}
 * ============================================================ */

void set_element_mat_NR_Aphi_JA(
    MONOLIS* monolis,
    BBFE_DATA* fe,
    BBFE_BASIS* basis,
    NEDELEC* ned,
    const double* x_curr,
    double dt
){
    const int np = basis->num_integ_points;
    double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
    double* val_ip_C    = BB_std_calloc_1d_double(val_ip_C, np);

    const double inv_dt = 1.0/dt;
    const int nEdge = fe->total_num_nodes;

    for(int e=0; e<fe->total_num_elems; ++e){

        int prop = ned->elem_prop[e];
        double sigma_mass_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_JA(prop, &sigma_mass_A, &sigma_cpl, &sigma_phi);

        BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

        /* ========= [1] A-A : Curl-Curl tangent ========= */
        for(int i=0; i<ned->local_num_edges; ++i){
            int gi = ned->nedelec_conn[e][i];
            int si = ned->edge_sign[e][i];

            for(int j=0; j<ned->local_num_edges; ++j){
                int gj = ned->nedelec_conn[e][j];
                int sj = ned->edge_sign[e][j];

                for(int p=0; p<np; ++p){
                    const double* ci = ned->curl_N_edge[e][p][i];
                    const double* cj = ned->curl_N_edge[e][p][j];

                    if(prop == 4){
                        //printf("g_ja_num_elem=%d fe->total_num_elems=%d\n", g_ja_num_elem, fe->total_num_elems);
                        //printf("g_ja_num_ip=%d basis->num_integ_points=%d\n", g_ja_num_ip, basis->num_integ_points);

                        double t = monolis_mpi_get_global_my_rank();
                        JA_IPState* s = &g_ja_state[e][p];
                        /* ensure material state at current iterate exists */
                        double B[3];
                        compute_B_ip(ned, e, p, x_curr, nEdge, B);

                        vec3_copy(s->b_iter_in, B);

                        ja_inverse_update(
                            &g_ja_param,
                            s->b_prev,
                            s->h_prev,
                            s->b_iter_in,
                            s->h_iter_out,
                            s->dh_db
                        );

                        double tmp[3];
                        mat3_vec(tmp, s->dh_db, cj);
                        val_ip_C[p] = vec3_dot(ci, tmp);
                    }else{
                        val_ip_C[p] = NU_LIN * vec3_dot(ci, cj);
                    }
                }

                double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, (double)(si*sj) * v);
            }
        }

        /* ========= [2] A-A : Mass (same as original) ========= */
        if(sigma_mass_A > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], ned->N_edge[e][p][j], sigma_mass_A);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0,
                        (double)(si*sj) * v * inv_dt);
                }
            }
        }

        /* ========= [3] Phi-Phi : same as original ========= */
        double sigma_laplace = (sigma_cpl > 0.0) ? sigma_cpl : sigma_phi;

        if(sigma_laplace > 0.0){
            for(int i=0; i<fe->local_num_nodes; ++i){
                int gi = fe->conn[e][i];
                for(int j=0; j<fe->local_num_nodes; ++j){
                    int gj = fe->conn[e][j];
                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], sigma_laplace);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gi, gj, 0, 0, v);
                }
            }
        }

        /* ========= [4] A-Phi : same as original ========= */
        if(sigma_cpl > 0.0){
            for(int n=0; n<fe->local_num_nodes; ++n){
                int gn = fe->conn[e][n];
                for(int j=0; j<ned->local_num_edges; ++j){
                    int gj = ned->nedelec_conn[e][j];
                    int sj = ned->edge_sign[e][j];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            fe->geo[e][p].grad_N[n], ned->N_edge[e][p][j], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gj, gn, 0, 0, (double)sj * v);
                }
            }
        }

        /* ========= [5] Phi-A : same as original ========= */
        if(sigma_cpl > 0.0){
            for(int i=0; i<ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];

                for(int n=0; n<fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];

                    for(int p=0; p<np; ++p){
                        val_ip_C[p] = BBFE_elemmat_mag_mat_mass(
                            ned->N_edge[e][p][i], fe->geo[e][p].grad_N[n], sigma_cpl);
                    }
                    double v = BBFE_std_integ_calc(np, val_ip_C, basis->integ_weight, Jacobian_ip);
                    monolis_add_scalar_to_sparse_matrix_R(monolis, gn, gi, 0, 0, (double)si * v * inv_dt);
                }
            }
        }
    }

    BB_std_free_1d_double(Jacobian_ip, np);
    BB_std_free_1d_double(val_ip_C, np);
}

/* ============================================================
 * Residual norm evaluation for line search
 *
 * The paper uses relaxation coefficient search based on minimizing
 * residual norm with trial alphas 1, 1/2, 1/4, ... :contentReference[oaicite:5]{index=5}
 * ============================================================ */

double compute_total_residual_norm_sq(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_trial,
    double dt,
    double current_time
){
    const int n_dof_total = sys->fe.total_num_nodes; /* adjust if total dof differs */
    double local = 0.0;
    double global = 0.0;

    monolis_clear_mat_value_R(&(sys->monolis));

    set_element_vec_NR_Aphi_JA(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->ned),
        x_prev,
        x_trial,
        dt,
        current_time
    );

    apply_dirichlet_bc_for_A_and_phi(
        &(sys->monolis),
        &(sys->fe),
        &(sys->bc),
        &(sys->ned)
    );

    local = local_sum_sq(sys->monolis.mat.R.B, n_dof_total);
    global = local;

    monolis_allreduce_R(1, &global, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    return global;
}

/* ============================================================
 * Relaxation coefficient search
 *
 * alpha = 1, 1/2, 1/4, ...
 * stop when W(alpha_next) > W(alpha_current)
 * based on the paper's strategy. :contentReference[oaicite:6]{index=6}
 * ============================================================ */

double search_relaxation_factor_JA(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    const double* delta,
    int n_dof_total,
    double dt,
    double current_time
){
    const int m_max = 8;

    double* x_trial = (double*)calloc(n_dof_total, sizeof(double));
    double best_alpha = 1.0;
    double prev_W = 1.0e+30;

    for(int m=0; m<=m_max; ++m){
        const double alpha = pow(0.5, (double)m);

        memcpy(x_trial, x_curr, sizeof(double)*n_dof_total);
        update_Aphi_NR(x_trial, delta, n_dof_total, alpha);

        double W = compute_total_residual_norm_sq(sys, x_prev, x_trial, dt, current_time);

        if(m == 0){
            best_alpha = alpha;
            prev_W = W;
            continue;
        }

        if(W > prev_W){
            break;
        }else{
            best_alpha = alpha;
            prev_W = W;
        }
    }

    free(x_trial);
    return best_alpha;
}

/* ============================================================
 * Example Newton loop skeleton
 * ============================================================ */

int solver_fom_NR_Aphi_JA(
    FE_SYSTEM* sys,
    double* x_prev,
    double* x_curr,
    int n_dof_total,
    double dt,
    double current_time)
{
    const int max_iter = 20;
    const double    tol_res = 1e-6;
    const double tol_dx = 1e-6;
    double* dx   = (double*)calloc((size_t)n_dof_total, sizeof(double));

    double res_local_old = local_sum_sq(sys->monolis.mat.R.B, n_dof_total);
    double res_global_old = res_local_old;

    /* 初期値 */
    for(int i=0; i<n_dof_total; ++i){
        x_curr[i] = x_prev[i];
    }


    if(dt==current_time){
        printf("initialize g_phase_current");
        
        initialize_g_phase_current(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->ned),
            x_prev,
            x_curr,
            dt,
            current_time
        );

    }

    for(int it=0; it<max_iter; ++it){

        monolis_clear_mat_value_R(&(sys->monolis));

        set_element_mat_NR_Aphi_JA(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->ned),
            x_curr,
            dt
        );

        set_element_vec_NR_Aphi_JA(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->ned),
            x_prev,
            x_curr,
            dt,
            current_time
        );

        apply_dirichlet_bc_for_A_and_phi(
            &(sys->monolis),
            &(sys->fe),
            &(sys->bc),
            &(sys->ned)
        );

        if(it == 0){
            res_local_old = local_sum_sq(sys->monolis.mat.R.B, n_dof_total);
            res_global_old = res_local_old;

            monolis_allreduce_R(1, &res_global_old, MONOLIS_MPI_SUM, sys->monolis_com.comm);

            res_global_old = sqrt(res_global_old);
        }

        monowrap_solve_R(
            &(sys->monolis),
            &(sys->monolis_com),
            dx,
            MONOLIS_ITER_BICGSAFE,
            MONOLIS_PREC_DIAG,
            sys->vals.mat_max_iter,
            sys->vals.mat_epsilon);

        monolis_clear_mat_value_R(&(sys->monolis));

        set_element_mat_NR_Aphi_JA(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->ned),
            x_curr,
            dt
        );

        set_element_vec_NR_Aphi_JA(
            &(sys->monolis),
            &(sys->fe),
            &(sys->basis),
            &(sys->ned),
            x_prev,
            x_curr,
            dt,
            current_time
        );

        apply_dirichlet_bc_for_A_and_phi(
            &(sys->monolis),
            &(sys->fe),
            &(sys->bc),
            &(sys->ned)
        );

        double alpha = search_relaxation_factor_JA(
            sys, x_prev, x_curr, dx, n_dof_total, dt, current_time
        );

        update_Aphi_NR(x_curr, dx, n_dof_total, alpha);

        double res_local = local_sum_sq(sys->monolis.mat.R.B, n_dof_total);
        double dx_local  = local_max_abs(dx, n_dof_total);

        double res_global = res_local;
        double dx_global  = dx_local;

        monolis_allreduce_R(1, &res_global, MONOLIS_MPI_SUM, sys->monolis_com.comm);
        monolis_allreduce_R(1, &dx_global,  MONOLIS_MPI_MAX, sys->monolis_com.comm);

        res_global = sqrt(res_global);

        if(res_global/res_global_old < tol_res){
            ja_commit_step_states(&(sys->fe), &(sys->basis), &(sys->ned), x_curr);
            return 1;
        }

        if(monolis_mpi_get_global_my_rank()==0){
            printf("res_global = %e, dx_global = %e", res_global/res_global_old, dx_global);
        }
        
    }

    return 0;
}


/* ============================================================
 * log_accuracy_metrics
 * Voltage-excitation version
 * ============================================================ */
void log_accuracy_metrics_JA(
    FE_SYSTEM* sys,
    const double* x_prev,
    const double* x_curr,
    int step,
    double t,
    double dt,
    double P_input
){
    BBFE_DATA* fe     = &(sys->fe);
    BBFE_BASIS* basis = &(sys->basis);
    NEDELEC* ned      = &(sys->ned);

    const int np    = basis->num_integ_points;
    const int nEdge = fe->total_num_nodes;
    const double inv_dt = 1.0 / dt;

    const double B_SAT_THRESHOLD = 1.5;
    const double T_period = 1.0 / FREQ_HZ;

    const double search_radius = 2.0e-3;
    const double target_pos1[3] = { 0.0025, 0.0350, 0.0200 };
    const double target_pos2[3] = { 0.0350, 0.0350, 0.0200 };
    const double target_pos3[3] = { 0.0675, 0.0350, 0.0200 };

    const int N_LOG = 29;
    double local_vals[29];
    for(int i = 0; i < N_LOG; i++) local_vals[i] = 0.0;

    const double NU0 = 1.0 / (4.0 * M_PI * 1.0e-7);

    /* applied phase voltages */
    double V_app[3];
    for(int k = 0; k < 3; k++){
        V_app[k] = get_phase_voltage(k, t);
    }

    /* current values used in this step */
    double I_src[3] = {
        g_phase_current[0],
        g_phase_current[1],
        g_phase_current[2]
    };

    const char* fname;
    FILE* fp_in;
    char id[128];
    int tmp;
    int num_internal_elems = 0;

    fname = monolis_get_global_input_file_name(
        MONOLIS_DEFAULT_TOP_DIR,
        MONOLIS_DEFAULT_PART_DIR,
        "graph_nedelec_elem.dat.n_internal"
    );
    fp_in = BBFE_sys_read_fopen(fp_in, fname, sys->cond.directory);
    fscanf(fp_in, "%s %d", id, &(tmp));
    fscanf(fp_in, "%d", &(num_internal_elems));
    fclose(fp_in);

    for(int e = 0; e < num_internal_elems; ++e){
        int prop = ned->elem_prop[e];

        int region_type = 0;
        if(prop == 1 || prop == 2 || prop == 3) region_type = 1;
        else if(prop == 4) region_type = 2;
        else if(prop == 5) region_type = 3;

        double sigma_A, sigma_cpl, sigma_phi;
        get_sigmas_for_prop_JA(prop, &sigma_A, &sigma_cpl, &sigma_phi);

        COIL_INFO coil;
        int is_coil = get_coil_info(prop, &coil);
        int phase_idx = -1;
        if(prop == 1) phase_idx = 0;
        if(prop == 2) phase_idx = 1;
        if(prop == 3) phase_idx = 2;

        double center[3] = {0.0, 0.0, 0.0};
        for(int n = 0; n < fe->local_num_nodes; ++n){
            int node_id = fe->conn[e][n];
            center[0] += fe->x[node_id][0];
            center[1] += fe->x[node_id][1];
            center[2] += fe->x[node_id][2];
        }
        center[0] /= fe->local_num_nodes;
        center[1] /= fe->local_num_nodes;
        center[2] /= fe->local_num_nodes;

        double dist1 = sqrt(
            (center[0]-target_pos1[0])*(center[0]-target_pos1[0]) +
            (center[1]-target_pos1[1])*(center[1]-target_pos1[1]) +
            (center[2]-target_pos1[2])*(center[2]-target_pos1[2])
        );
        double dist2 = sqrt(
            (center[0]-target_pos2[0])*(center[0]-target_pos2[0]) +
            (center[1]-target_pos2[1])*(center[1]-target_pos2[1]) +
            (center[2]-target_pos2[2])*(center[2]-target_pos2[2])
        );
        double dist3 = sqrt(
            (center[0]-target_pos3[0])*(center[0]-target_pos3[0]) +
            (center[1]-target_pos3[1])*(center[1]-target_pos3[1]) +
            (center[2]-target_pos3[2])*(center[2]-target_pos3[2])
        );

        for(int p = 0; p < np; ++p){
            double w_detJ = basis->integ_weight[p] * fe->geo[e][p].Jacobian;

            if(region_type == 1) local_vals[11] += w_detJ;
            if(region_type == 2) local_vals[12] += w_detJ;
            if(region_type == 3) local_vals[13] += w_detJ;

            double B_curr[3], B_prev[3];
            compute_B_ip(ned, e, p, x_curr, nEdge, B_curr);
            compute_B_ip(ned, e, p, x_prev, nEdge, B_prev);

            double Bmag = norm3(B_curr);
            if(Bmag > local_vals[14]) local_vals[14] = Bmag;

            double current_nu;
            if(region_type == 2){
                current_nu = get_reluctivity_nu(fmax(Bmag, 1.0e-14));
            }else{
                current_nu = NU0;
            }

            if(region_type == 3){
                if(current_nu > local_vals[15]) local_vals[15] = current_nu;
            }

            if(region_type == 2 && Bmag > B_SAT_THRESHOLD){
                local_vals[19] += w_detJ;
            }

            if(region_type == 2){
                if(dist1 < search_radius){
                    local_vals[20] = current_nu;
                    local_vals[21] = dist1;
                    local_vals[22] = Bmag;
                }
                if(dist2 < search_radius){
                    local_vals[23] = current_nu;
                    local_vals[24] = dist2;
                    local_vals[25] = Bmag;
                }
                if(dist3 < search_radius){
                    local_vals[26] = current_nu;
                    local_vals[27] = dist3;
                    local_vals[28] = Bmag;
                }
            }

            double A_curr[3] = {0.0, 0.0, 0.0};
            double A_prev[3] = {0.0, 0.0, 0.0};
            double grad_phi[3] = {0.0, 0.0, 0.0};

            for(int i = 0; i < ned->local_num_edges; ++i){
                int gi = ned->nedelec_conn[e][i];
                int si = ned->edge_sign[e][i];
                const double* N = ned->N_edge[e][p][i];
                double ac = x_curr[gi] * si;
                double ap = x_prev[gi] * si;
                for(int d = 0; d < 3; d++){
                    A_curr[d] += ac * N[d];
                    A_prev[d] += ap * N[d];
                }
            }

            if(sigma_cpl > 0.0 || sigma_phi > 0.0){
                for(int n = 0; n < fe->local_num_nodes; ++n){
                    int gn = fe->conn[e][n];
                    double ph = x_curr[gn];
                    const double* dN = fe->geo[e][p].grad_N[n];
                    for(int d = 0; d < 3; d++){
                        grad_phi[d] += ph * dN[d];
                    }
                }
            }

            double E_vec[3];
            for(int d = 0; d < 3; d++){
                E_vec[d] = -(A_curr[d] - A_prev[d]) * inv_dt - grad_phi[d];
            }

            if(is_coil && phase_idx >= 0){
                double x_ip[3];
                get_interp_coords(e, p, fe, basis, x_ip);

                double d_vec[3] = {
                    x_ip[0] - coil.center[0],
                    x_ip[1] - coil.center[1],
                    x_ip[2] - coil.center[2]
                };
                double da = dot3(d_vec, coil.axis);
                double r_vec[3] = {
                    d_vec[0] - da * coil.axis[0],
                    d_vec[1] - da * coil.axis[1],
                    d_vec[2] - da * coil.axis[2]
                };
                double t_vec[3] = {
                    coil.axis[1]*r_vec[2] - coil.axis[2]*r_vec[1],
                    coil.axis[2]*r_vec[0] - coil.axis[0]*r_vec[2],
                    coil.axis[0]*r_vec[1] - coil.axis[1]*r_vec[0]
                };
                double t_norm = norm3(t_vec);

                if(t_norm > 1.0e-12){
                    double J0_mag = coil.turns / coil.area;
                    double J0[3] = {
                        J0_mag*t_vec[0]/t_norm,
                        J0_mag*t_vec[1]/t_norm,
                        J0_mag*t_vec[2]/t_norm
                    };
                    local_vals[3 + phase_idx] += dot3(E_vec, J0) * w_detJ;
                }
            }

            double dLoss = 0.0;
            if(sigma_A > 0.0){
                dLoss = sigma_A * dot3(E_vec, E_vec) * w_detJ;
            }

            local_vals[6] += dLoss;
            if(region_type == 1) local_vals[16] += dLoss;
            if(region_type == 2) local_vals[17] += dLoss;
            if(region_type == 3) local_vals[18] += dLoss;

            double H_vec[3] = {
                current_nu*B_curr[0],
                current_nu*B_curr[1],
                current_nu*B_curr[2]
            };
            double dB_dt[3] = {
                (B_curr[0]-B_prev[0])*inv_dt,
                (B_curr[1]-B_prev[1])*inv_dt,
                (B_curr[2]-B_prev[2])*inv_dt
            };
            double dW_val = dot3(H_vec, dB_dt) * w_detJ;

            local_vals[7] += dW_val;
            if(region_type == 1) local_vals[8]  += dW_val;
            if(region_type == 2) local_vals[9]  += dW_val;
            if(region_type == 3) local_vals[10] += dW_val;
        }
    }

    double global_vals[29];

    double val_sum[29];
    for(int i = 0; i < 29; i++) val_sum[i] = local_vals[i];
    val_sum[14] = 0.0;
    val_sum[15] = 0.0;
    for(int i = 20; i < 29; i++) val_sum[i] = 0.0;

    monolis_allreduce_R(29, val_sum, MONOLIS_MPI_SUM, sys->monolis_com.comm);

    double val_max[29];
    for(int i = 0; i < 29; i++) val_max[i] = 0.0;
    val_max[14] = local_vals[14];
    val_max[15] = local_vals[15];
    for(int i = 20; i < 29; i++) val_max[i] = local_vals[i];

    monolis_allreduce_R(29, val_max, MONOLIS_MPI_MAX, sys->monolis_com.comm);

    for(int i = 0; i < 29; i++) global_vals[i] = val_sum[i];
    global_vals[14] = val_max[14];
    global_vals[15] = val_max[15];
    for(int i = 20; i < 29; i++) global_vals[i] = val_max[i];

    if(sys->monolis_com.my_rank == 0){
        double V_fem[3] = { global_vals[3], global_vals[4], global_vals[5] };

        double V_term[3];
        for(int k = 0; k < 3; k++){
            /* applied terminal voltage */
            V_term[k] = V_app[k];
        }

        /* current update for next step
         * sign convention follows original V_term = R*I - V_fem
         * => I = (V_app + V_fem)/R
         */
        for(int k = 0; k < 3; k++){
            g_phase_current_prev[k] = g_phase_current[k];
            g_phase_current[k] = (V_app[k] + V_fem[k]) / R_COIL;
        }

        double p1_nu = global_vals[20], p1_B = global_vals[22];
        double p2_nu = global_vals[23], p2_B = global_vals[25];
        double p3_nu = global_vals[26], p3_B = global_vals[28];

        double Loss_Coil  = global_vals[16];
        double Loss_Core  = global_vals[17];
        double Loss_Air   = global_vals[18];
        double Loss_Total = Loss_Coil + Loss_Core + Loss_Air;

        double CoreVol = global_vals[12];
        double CoreSatVol = global_vals[19];
        double CoreSatRatio = 0.0;
        if(CoreVol > 1.0e-30) CoreSatRatio = CoreSatVol / CoreVol;

        static int initialized = 0;
        static double run_max_B = -1.0;
        static double run_max_B_time = 0.0;
        static int run_max_B_step = -1;

        static double cycle_start_time = 0.0;
        static double sum_Iu2_dt = 0.0;
        static double sum_Iv2_dt = 0.0;
        static double sum_Iw2_dt = 0.0;
        static double sum_Pin_dt = 0.0;
        static double sum_LossCore_dt = 0.0;
        static double sum_LossTotal_dt = 0.0;
        static double sum_CoreSatRatio_dt = 0.0;
        static double cycle_max_B = 0.0;

        if(step == 0 || !initialized){
            initialized = 1;
            run_max_B = -1.0;
            run_max_B_time = t;
            run_max_B_step = step;

            cycle_start_time = t;
            sum_Iu2_dt = sum_Iv2_dt = sum_Iw2_dt = 0.0;
            sum_Pin_dt = 0.0;
            sum_LossCore_dt = 0.0;
            sum_LossTotal_dt = 0.0;
            sum_CoreSatRatio_dt = 0.0;
            cycle_max_B = 0.0;
        }

        if(global_vals[14] > run_max_B){
            run_max_B = global_vals[14];
            run_max_B_time = t;
            run_max_B_step = step;
        }

        sum_Iu2_dt += I_src[0]*I_src[0]*dt;
        sum_Iv2_dt += I_src[1]*I_src[1]*dt;
        sum_Iw2_dt += I_src[2]*I_src[2]*dt;
        sum_Pin_dt += P_input * dt;
        sum_LossCore_dt += Loss_Core * dt;
        sum_LossTotal_dt += Loss_Total * dt;
        sum_CoreSatRatio_dt += CoreSatRatio * dt;

        if(global_vals[14] > cycle_max_B) cycle_max_B = global_vals[14];

        printf(
            "  [VoltageExc] Step:%d t=%.4e | "
            "Vapp=(%.4e, %.4e, %.4e) | "
            "I=(%.4e, %.4e, %.4e) -> next=(%.4e, %.4e, %.4e) | "
            "MaxB=%.4e | SatRatio=%.4e\n",
            step, t,
            V_app[0], V_app[1], V_app[2],
            I_src[0], I_src[1], I_src[2],
            g_phase_current[0], g_phase_current[1], g_phase_current[2],
            global_vals[14], CoreSatRatio
        );

        FILE* fp_out;
        fp_out = BBFE_sys_write_add_fopen(fp_out, "terminal_voltage_log.csv", sys->cond.directory);

        if(step == 0){
            fprintf(
                fp_out,
                "Step,Time,"
                "Vapp_u,Vapp_v,Vapp_w,"
                "I_u,I_v,I_w,"
                "I_next_u,I_next_v,I_next_w,"
                "V_t_u,V_t_v,V_t_w,"
                "EMF_u,EMF_v,EMF_w,"
                "Power_In,"
                "Max_B,Max_Nu,"
                "Loss_Coil,Loss_Core,Loss_Air,Loss_Total,"
                "Core_Vol,Core_SatVol,Core_SatRatio,"
                "P1_Nu,P1_B,P2_Nu,P2_B,P3_Nu,P3_B,"
                "Run_Max_B,Run_Max_B_Time,Run_Max_B_Step\n"
            );
        }

        fprintf(
            fp_out,
            "%d,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,"
            "%.6e,%.6e,"
            "%.6e,%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,"
            "%.6e,%.6e,%d\n",
            step, t,
            V_app[0], V_app[1], V_app[2],
            I_src[0], I_src[1], I_src[2],
            g_phase_current[0], g_phase_current[1], g_phase_current[2],
            V_term[0], V_term[1], V_term[2],
            -V_fem[0], -V_fem[1], -V_fem[2],
            P_input,
            global_vals[14], global_vals[15],
            Loss_Coil, Loss_Core, Loss_Air, Loss_Total,
            CoreVol, CoreSatVol, CoreSatRatio,
            p1_nu, p1_B, p2_nu, p2_B, p3_nu, p3_B,
            run_max_B, run_max_B_time, run_max_B_step
        );
        fclose(fp_out);

        if((t - cycle_start_time + dt) >= T_period - 1.0e-14){
            double cycle_len = t - cycle_start_time + dt;
            if(cycle_len < 1.0e-14) cycle_len = T_period;

            double Irms_u = sqrt(sum_Iu2_dt / cycle_len);
            double Irms_v = sqrt(sum_Iv2_dt / cycle_len);
            double Irms_w = sqrt(sum_Iw2_dt / cycle_len);

            double Pin_avg = sum_Pin_dt / cycle_len;
            double LossCore_avg = sum_LossCore_dt / cycle_len;
            double LossTotal_avg = sum_LossTotal_dt / cycle_len;
            double CoreSatRatio_avg = sum_CoreSatRatio_dt / cycle_len;

            FILE* fp_sum;
            fp_sum = BBFE_sys_write_add_fopen(fp_sum, "transmag_cycle_summary.csv", sys->cond.directory);

            if(step == 0 || fabs(cycle_start_time) < 1.0e-14){
                fprintf(
                    fp_sum,
                    "CycleStart,CycleEnd,CycleLength,"
                    "Irms_u,Irms_v,Irms_w,"
                    "Pin_avg,LossCore_avg,LossTotal_avg,"
                    "CoreSatRatio_avg,CycleMax_B,"
                    "Run_Max_B,Run_Max_B_Time,Run_Max_B_Step\n"
                );
            }

            fprintf(
                fp_sum,
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,%.6e,"
                "%.6e,%.6e,"
                "%.6e,%.6e,%d\n",
                cycle_start_time, t + dt, cycle_len,
                Irms_u, Irms_v, Irms_w,
                Pin_avg, LossCore_avg, LossTotal_avg,
                CoreSatRatio_avg, cycle_max_B,
                run_max_B, run_max_B_time, run_max_B_step
            );
            fclose(fp_sum);

            cycle_start_time = t + dt;
            sum_Iu2_dt = 0.0;
            sum_Iv2_dt = 0.0;
            sum_Iw2_dt = 0.0;
            sum_Pin_dt = 0.0;
            sum_LossCore_dt = 0.0;
            sum_LossTotal_dt = 0.0;
            sum_CoreSatRatio_dt = 0.0;
            cycle_max_B = 0.0;
        }
    }
}