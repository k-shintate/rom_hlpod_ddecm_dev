#include "team11_core.h"

enum { VTK_TETRA = 10, VTK_HEXAHEDRON = 12 };


/* ======= 球ベッセル ======= */
void s_sin_cos(double x,double *sx,double *cx){ *sx=sin(x); *cx=cos(x); }
static inline double j0_sph(double x){ if(fabs(x)<1e-12) return 1.0; double s,c; s_sin_cos(x,&s,&c); return s/x; }
static inline double y0_sph(double x){ if(fabs(x)<1e-12) return 0.0; double s,c; s_sin_cos(x,&s,&c); return -c/x; }
static inline double j1_sph(double x){
    if(fabs(x)<1e-6) return x/3.0;
    double s,c; s_sin_cos(x,&s,&c); return s/(x*x) - c/x;
}

static inline double y1_sph(double x){
    double s,c; s_sin_cos(x,&s,&c); return -c/(x*x) - s/x;
}

void sph_to_cart(double th,double ph,double Br,double Bt,double B[3]){
    double sth=sin(th), cth=cos(th), s=sin(ph), c=cos(ph);
    double erx=sth*c, ery=sth*s, erz=cth;
    double etx=cth*c, ety=cth*s, etz=-sth;
    B[0]=Br*erx + Bt*etx;
    B[1]=Br*ery + Bt*ety;
    B[2]=Br*erz + Bt*etz;
}

/* ======= TEAM#11 固有モード（9モード） =======
 * tau = mu0*sigma/lambda^2 で算出済み（単位 s）
 * A, alpha は B0=1T のステップに対する値
 */
typedef struct { double lambda, tau, A, alpha; } Mode;
static const Mode TEAM11_MODES[9] = {
  {109.442332, 52.4576935254833e-3, -0.146072, -7.186145},
  {646.869386,  1.501573766604363e-3,  0.025943, -2.360774},
  {1266.114873, 0.39195268824446033e-3,-0.022434, -0.728356},
  {1891.300433, 0.1756543156286042e-3, 0.016342, -0.450044},
  {2518.039735, 0.09909567728327648e-3,-0.012634, -0.328911},
  {3145.407739, 0.06350763884388948e-3, 0.010249, -0.260027},
  {3773.091600, 0.044135207303649365e-3,-0.008605, -0.215304},
  {4400.956395, 0.03244036539134097e-3, 0.007410, -0.183836},
  {5028.934437, 0.024844366140195407e-3,-0.006503, -0.160456}
};

/* ======= 物性・幾何（TEAM#11 既定） ======= */
typedef struct {
    double a, b;      /* inner/outer radius [m] */
    double mu, sigma; /* [H/m], [S/m] */
    double B0;        /* step amplitude [T] (+z) */
    double cx, cy, cz;/* sphere center [m] */
} Team11Param;

void team11_default_param(Team11Param *p){
    p->a=0.05; p->b=0.055;
    p->mu=4.0*M_PI*1e-7; p->sigma=5.0e8;
    p->B0=1.0;
    p->cx=p->cy=p->cz=0.0;
}

/* ======= 厳密解で一点の B を評価 =======
 * - 内部 r<a： a(r,t)=B0 r/2 - r Σ A alpha e^{-t/τ} ⇒ Br, Bt が一定係数×cos/sin
 * - 導体 a<r<b： φ=j1+β y1, β=-j0(λb)/y0(λb) を用いた式
 * - 外部 r>b： C(t)=Σ A b^2 φ(b) e^{-t/τ} からの双極子補正
 */
void team11_eval_B_at_point(const Team11Param *P, double x,double y,double z, double t, double B[3]){
    const double a=P->a, b=P->b, B0=P->B0;
    /* 球中心へ平行移動 */
    double X=x-P->cx, Y=y-P->cy, Z=z-P->cz;
    double r = sqrt(X*X+Y*Y+Z*Z);

    if(r<1e-14){
        double s=0.0; for(int i=0;i<9;++i) s += TEAM11_MODES[i].A*TEAM11_MODES[i].alpha*exp(-t/TEAM11_MODES[i].tau);
        B[0]=0.0; B[1]=0.0; B[2]=B0 - 2.0*s; return;
    }
    double th = acos(fmax(-1.0, fmin(1.0, Z/r)));
    double ph = atan2(Y,X);

    /* 先に C(t)=Σ A b^2 φ(b) e^{-t/τ} を計算 */
    double C=0.0;
    for(int i=0;i<9;++i){
        double lam=TEAM11_MODES[i].lambda, xb=lam*b;
        double beta = - j0_sph(xb)/y0_sph(xb);
        double phib = j1_sph(xb) + beta*y1_sph(xb);
        C += TEAM11_MODES[i].A * (b*b*phib) * exp(-t/TEAM11_MODES[i].tau);
    }

    if(r<a){
        double s=0.0; for(int i=0;i<9;++i) s += TEAM11_MODES[i].A*TEAM11_MODES[i].alpha*exp(-t/TEAM11_MODES[i].tau);
        double S=B0 - 2.0*s;
        double Br= S*cos(th), Bt= -S*sin(th);
        sph_to_cart(th,ph,Br,Bt,B); return;
    }else if(r<=b){
        double Sphi=0.0, Spr=0.0;
        for(int i=0;i<9;++i){
            double lam=TEAM11_MODES[i].lambda, xr=lam*r, xb=lam*b;
            double beta = - j0_sph(xb)/y0_sph(xb);
            double phi = j1_sph(xr) + beta*y1_sph(xr);
            double phiplus = -phi + r*lam*( j0_sph(xr) + beta*y0_sph(xr) );
            double e = exp(-t/TEAM11_MODES[i].tau);
            Sphi += TEAM11_MODES[i].A*phi*e;
            Spr  += TEAM11_MODES[i].A*phiplus*e;
        }
        double Br =  cos(th)*( B0 - 2.0*Sphi/r );
        double Bt = -(sin(th)/r)*( B0*r - Spr );
        sph_to_cart(th,ph,Br,Bt,B); return;
    }else{
        double invr3 = 1.0/(r*r*r);
        double Br =  cos(th)*( B0 - 2.0*C*invr3 );
        double Bt = -sin(th)*( B0 + C*invr3 );
        sph_to_cart(th,ph,Br,Bt,B); return;
    }
}

/* ======= VTK 出力（UnstructuredGrid, ASCII） ======= */
void write_unstructured_vtk_ascii(
    BBFE_DATA* fe,
    const double (*B_node)[3],
    const char* directory,
    const char* filename)
{
    char path[1024];
    if (directory && directory[0]) snprintf(path, sizeof(path), "%s/%s", directory, filename);
    else snprintf(path, sizeof(path), "%s", filename);

    FILE* fp = fopen(path, "w");
    if (!fp) { perror("fopen"); return; }

    const int N = fe->total_num_nodes;
    const int M = fe->total_num_elems;
    const int L = fe->local_num_nodes;

    /* 0/1 始まり自動判定 */
    int min_id = INT_MAX;
    for (int e=0; e<M; ++e)
        for (int a=0; a<L; ++a)
            if (fe->conn[e][a] < min_id) min_id = fe->conn[e][a];
    const int shift = (min_id == 1) ? 1 : 0;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "B field at nodes (TEAM#11 analytic)\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", N);
    for (int i=0; i<N; ++i)
        fprintf(fp, "%.9g %.9g %.9g\n", fe->x[i][0], fe->x[i][1], fe->x[i][2]);

    fprintf(fp, "CELLS %d %d\n", M, M*(L+1));
    for (int e=0; e<M; ++e) {
        if (L == 4) {
            fprintf(fp, "4 %d %d %d %d\n",
                fe->conn[e][0]-shift, fe->conn[e][1]-shift, fe->conn[e][2]-shift, fe->conn[e][3]-shift);
        } else if (L == 8) {
            fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",
                fe->conn[e][0]-shift, fe->conn[e][1]-shift, fe->conn[e][2]-shift, fe->conn[e][3]-shift,
                fe->conn[e][4]-shift, fe->conn[e][5]-shift, fe->conn[e][6]-shift, fe->conn[e][7]-shift);
        } else {
            fprintf(stderr, "Unsupported local_num_nodes=%d\n", L);
            fclose(fp); return;
        }
    }

    fprintf(fp, "CELL_TYPES %d\n", M);
    const int ctype = (L == 4) ? VTK_TETRA : VTK_HEXAHEDRON;
    for (int e=0; e<M; ++e) fprintf(fp, "%d\n", ctype);

    fprintf(fp, "POINT_DATA %d\n", N);
    fprintf(fp, "VECTORS B_node float\n");
    for (int i=0; i<N; ++i)
        fprintf(fp, "%.9g %.9g %.9g\n", B_node[i][0], B_node[i][1], B_node[i][2]);

    fclose(fp);
}

void fill_B_nodes_team11(BBFE_DATA* fe, double t_sec, const Team11Param *P, double (*B_node)[3]){
    for (int i=0; i<fe->total_num_nodes; ++i){
        team11_eval_B_at_point(P, fe->x[i][0], fe->x[i][1], fe->x[i][2], t_sec, B_node[i]);
    }
}

/* ======= 例：t=1,5,10,20 ms を出力 ======= */
void output_B_node_vtk_team11_series(BBFE_DATA* fe, const char* dir){
    Team11Param P; team11_default_param(&P);

    const double times_ms[] = {1,5,10,20};
    const int Nt = (int)(sizeof(times_ms)/sizeof(times_ms[0]));

    double (*B_node)[3] = (double(*)[3])calloc((size_t)fe->total_num_nodes, sizeof *B_node);
    if(!B_node){ perror("calloc"); return; }

    for(int it=0; it<Nt; ++it){
        const double t = times_ms[it] * 1e-3;
        fill_B_nodes_team11(fe, t, &P, B_node);

        char name[128];
        snprintf(name, sizeof(name), "B_team11_%04dms.vtk", (int)times_ms[it]);
        write_unstructured_vtk_ascii(fe, (const double(*)[3])B_node, dir, name);
    }
    free(B_node);
}

static double team11_C_unitstep(double t_sec, double b){
    double C = 0.0;
    for(int i=0;i<9;++i){
        double lam = TEAM11_MODES[i].lambda;
        double xb  = lam*b;
        double beta = - j0_sph(xb)/y0_sph(xb);
        double phib = j1_sph(xb) + beta*y1_sph(xb);
        C += TEAM11_MODES[i].A * (b*b*phib) * exp(-t_sec/TEAM11_MODES[i].tau);
    }
    return C; // B0=1T のステップ応答
}

void A_uniform(const double B[3], const double x[3], double A[3]){
    // 0.5 * (B x x)
    A[0] = 0.5*( B[1]*x[2] - B[2]*x[1] );
    A[1] = 0.5*( B[2]*x[0] - B[0]*x[2] );
    A[2] = 0.5*( B[0]*x[1] - B[1]*x[0] );
}

void A_dipole_from_C(double C, const double x[3], double A[3]){
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r  = sqrt(r2);
    double r3 = r2*r;
    if(r3 < 1e-30){ A[0]=A[1]=A[2]=0.0; return; }

    // A_dip = -C * (ez x x)/r^3 = -C*(-y, x, 0)/r^3 = (C*y/r^3, -C*x/r^3, 0)
    A[0] =  C * x[1] / r3;
    A[1] = -C * x[0] / r3;
    A[2] =  0.0;
}

static inline double edge_integral_A_dot_dl(
    const double B[3], double C,
    const double x1[3], const double x2[3])
{
    double dx[3] = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };

    // 2点Gauss on s in [0,1]
    const double gp[2] = { 0.5 - 0.2886751345948129, 0.5 + 0.2886751345948129 };
    const double gw[2] = { 0.5, 0.5 };

    double In = 0.0;
    for(int q=0;q<2;++q){
        double s = gp[q];
        double x[3] = { x1[0] + s*dx[0], x1[1] + s*dx[1], x1[2] + s*dx[2] };

        double Au[3], Ad[3], A[3];
        A_uniform(B, x, Au);
        A_dipole_from_C(C, x, Ad);
        A[0]=Au[0]+Ad[0]; A[1]=Au[1]+Ad[1]; A[2]=Au[2]+Ad[2];

        In += gw[q] * (A[0]*dx[0] + A[1]*dx[1] + A[2]*dx[2]);
    }
    return In;
}

void apply_nedelec_boundary_conditions_TEAM11_time(
    MONOLIS*    monolis,
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    double      Bz_t,
    double      t_sec
){
    const double b = 0.055;
    const double C = team11_C_unitstep(t_sec, b) * (Bz_t / 1.0);
    const double B_vec[3] = {0.0, 0.0, Bz_t};

    const int nen = fe->local_num_nodes;
    int num_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if (nen == 4) {
        num_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if (nen == 8) {
        num_edges = 12;
        edge_tbl = hex_edge_conn;
    } else {
        fprintf(stderr, "Unsupported element type: local_num_nodes=%d\n", nen);
        exit(EXIT_FAILURE);
    }

    for (int e = 0; e < fe->total_num_elems; ++e){
        for (int ed = 0; ed < num_edges; ++ed){

            const int ln1 = edge_tbl[ed][0];
            const int ln2 = edge_tbl[ed][1];
            const int gn1 = fe->conn[e][ln1];
            const int gn2 = fe->conn[e][ln2];

            const int gid = ned->nedelec_conn[e][ed];

            if (!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;

            const double x1[3] = { fe->x[gn1][0], fe->x[gn1][1], fe->x[gn1][2] };
            const double x2[3] = { fe->x[gn2][0], fe->x[gn2][1], fe->x[gn2][2] };

            double Aint = edge_integral_A_dot_dl(B_vec, C, x1, x2);

            if (ned->edge_sign) {
                Aint *= ned->edge_sign[e][ed];
            }

            const int ge = gid;

            monolis_set_Dirichlet_bc_R(monolis, monolis->mat.R.B, ge, 0, Aint);
        }
    }
}



void apply_nedelec_boundary_conditions_TEAM11_NR(
    MONOLIS*    monolis,
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    double      Bz_t,
    double      t_sec,
    double*      val)
{
    const double b = 0.055;
    const double C = team11_C_unitstep(t_sec, b) * (Bz_t / 1.0);
    const double B_vec[3] = {0.0, 0.0, Bz_t};

    const int nen = fe->local_num_nodes;
    int num_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if (nen == 4) {
        num_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if (nen == 8) {
        num_edges = 12;
        edge_tbl = hex_edge_conn;
    } else {
        fprintf(stderr, "Unsupported element type: local_num_nodes=%d\n", nen);
        exit(EXIT_FAILURE);
    }

    for (int e = 0; e < fe->total_num_elems; ++e){

        for (int ed = 0; ed < num_edges; ++ed){

            const int ln1 = edge_tbl[ed][0];
            const int ln2 = edge_tbl[ed][1];
            const int gn1 = fe->conn[e][ln1];
            const int gn2 = fe->conn[e][ln2];

            const int gid = ned->nedelec_conn[e][ed];

            if (!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;

            const double x1[3] = { fe->x[gn1][0], fe->x[gn1][1], fe->x[gn1][2] };
            const double x2[3] = { fe->x[gn2][0], fe->x[gn2][1], fe->x[gn2][2] };

            double Aint = edge_integral_A_dot_dl(B_vec, C, x1, x2);

            if (ned->edge_sign) {
                Aint *= ned->edge_sign[e][ed]; /* ±1 */
            }

            const int ge = gid;

            val[ge] = Aint;
        }
    }

}
