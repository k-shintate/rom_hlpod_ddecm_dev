
#include "convdiff_core.h"
#include "nedelec_core.h"
#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

static const int BUFFER_SIZE = 10000;

const double mu0 = 4.0*M_PI*1e-7; // H/m
//const double Nu  = 1.0 / mu0;     // ← ここを μ0 に合わせる

void compute_B_cell_average(
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    NEDELEC*     ned,
    const double* sol,
    double**     B_cell
){
    const int ne = fe->total_num_elems;
    const int np = basis->num_integ_points;

    if(!B_cell){
        fprintf(stderr, "compute_B_cell_average: B_cell must be pre-allocated\n");
        exit(EXIT_FAILURE);
    }

    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    for(int e = 0; e < ne; ++e){
        double num[3] = {0.0, 0.0, 0.0};
        double den = 0.0;

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        for(int p = 0; p < np; ++p){
            double Tp[3]      = {0.0, 0.0, 0.0};
            double gradPhi[3] = {0.0, 0.0, 0.0};
            double Hp[3];
            double Bp[3];

            if(ned->elem_prop==4){
                for(int j = 0; j < ned->local_num_edges; ++j){
                    const int ge = ned->nedelec_conn[e][j];
                    const double cj =
                        (ned->edge_sign ? ned->edge_sign[e][j] : 1) * sol[ge];

                    Tp[0] += cj * ned->N_edge[e][p][j][0];
                    Tp[1] += cj * ned->N_edge[e][p][j][1];
                    Tp[2] += cj * ned->N_edge[e][p][j][2];
                }
            }

            for(int n = 0; n < fe->local_num_nodes; ++n){
                const int gn = fe->conn[e][n];
                const double phi_n = sol[gn];

                gradPhi[0] += phi_n * fe->geo[e][p].grad_N[n][0];
                gradPhi[1] += phi_n * fe->geo[e][p].grad_N[n][1];
                gradPhi[2] += phi_n * fe->geo[e][p].grad_N[n][2];
            }

            Hp[0] = Tp[0] - gradPhi[0];
            Hp[1] = Tp[1] - gradPhi[1];
            Hp[2] = Tp[2] - gradPhi[2];

            Bp[0] = mu0 * Hp[0];
            Bp[1] = mu0 * Hp[1];
            Bp[2] = mu0 * Hp[2];

            const double w = basis->integ_weight[p] * J_ip[p];
            num[0] += w * Bp[0];
            num[1] += w * Bp[1];
            num[2] += w * Bp[2];
            den    += w;
        }

        if(den > 0.0){
            B_cell[e][0] = num[0] / den;
            B_cell[e][1] = num[1] / den;
            B_cell[e][2] = num[2] / den;
        } else {
            B_cell[e][0] = 0.0;
            B_cell[e][1] = 0.0;
            B_cell[e][2] = 0.0;
        }
    }

    BB_std_free_1d_double(J_ip, np);
}
/*
void compute_B_cell_average(
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    NEDELEC*     ned,
    const double* Aphi,
    double**      B_cell
){
    const int ne = fe->total_num_elems;
    const int np = basis->num_integ_points;

    if(!B_cell){
        fprintf(stderr, "compute_B_cell_average: B_cell must be pre-allocated\n");
        exit(EXIT_FAILURE);
    }

    double* J_ip = BB_std_calloc_1d_double(J_ip, np);

    for(int e=0; e<ne; ++e){
        double num[3] = {0,0,0};
        double den    = 0.0;

        BBFE_elemmat_set_Jacobian_array(J_ip, np, e, fe);

        for(int p=0; p<np; ++p){
            // B(x_ip) = curl A_h = Σ_j c_j * curl N_j
            double Bp[3] = {0,0,0};
            for(int j=0; j<ned->local_num_edges; ++j){
                const int ge = ned->nedelec_conn[e][j];
                const double cj =
                        ned->edge_sign[e][j] * Aphi[ge];

                Bp[0] += cj * ned->curl_N_edge[e][p][j][0];
                Bp[1] += cj * ned->curl_N_edge[e][p][j][1];
                Bp[2] += cj * ned->curl_N_edge[e][p][j][2];
            }

            const double w = basis->integ_weight[p] * J_ip[p];
            num[0] += w * Bp[0];
            num[1] += w * Bp[1];
            num[2] += w * Bp[2];
            den    += w;
        }

        if(den > 0){
            B_cell[e][0] = num[0]/den;
            B_cell[e][1] = num[1]/den;
            B_cell[e][2] = num[2]/den;
        }else{
            B_cell[e][0] = B_cell[e][1] = B_cell[e][2] = 0.0;
        }
    }

    BB_std_free_1d_double(J_ip, np);
}
*/

// B_cell -> 節点B
void accumulate_B_cell_to_nodes(
    BBFE_DATA* fe,
    double** B_cell,
    double** B_node_out
){
    int N = fe->total_num_nodes, ne = fe->total_num_elems;
    int* cnt = (int*)calloc(N, sizeof(int));
    for (int i=0;i<N;++i){ B_node_out[i][0]=B_node_out[i][1]=B_node_out[i][2]=0.0; }

    for (int e=0; e<ne; ++e){
        for (int a=0;a<fe->local_num_nodes;++a){
            int g = fe->conn[e][a];
            B_node_out[g][0] += B_cell[e][0];
            B_node_out[g][1] += B_cell[e][1];
            B_node_out[g][2] += B_cell[e][2];
            cnt[g] += 1;
        }
    }
    for (int g=0; g<N; ++g) if (cnt[g]>0){
        B_node_out[g][0] /= cnt[g];
        B_node_out[g][1] /= cnt[g];
        B_node_out[g][2] /= cnt[g];
    }
    free(cnt);
}

void set_elem_types(
    BBFE_DATA* fe,
    NEDELEC* ned,
    double*     elem_type
){
    for(int e=0; e<fe->total_num_elems; ++e){
        for(int n=0; n<fe->local_num_nodes; ++n){

            int prop = ned->elem_prop[e];
            if(prop==6){
                elem_type[fe->conn[e][n]] = prop;
            }
        }
    }

    for(int e=0; e<fe->total_num_elems; ++e){
        for(int n=0; n<fe->local_num_nodes; ++n){

            int prop = ned->elem_prop[e];
            if(prop!=6){
                elem_type[fe->conn[e][n]] = prop;
            }
        }
    }
}

void output_B_node_vtk(
    BBFE_DATA* fe, BBFE_BASIS* basis, NEDELEC* ned, const double* Aphi,
    const char* filename, const char* directory)
{
    double** B_cell = BB_std_calloc_2d_double(B_cell, fe->total_num_elems, 3);
    double** B_node = BB_std_calloc_2d_double(B_node, fe->total_num_nodes, 3);
    double* elem_type = BB_std_calloc_1d_double(elem_type, fe->total_num_nodes);

    compute_B_cell_average(fe, basis, ned, Aphi, B_cell);
    accumulate_B_cell_to_nodes(fe, B_cell, B_node);

    FILE* fp = BBFE_sys_write_fopen(fp, filename, directory);
    switch (fe->local_num_nodes){
        case 4: BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA); break;
        case 8: BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON); break;
    }
    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    BB_vtk_write_point_vals_vector(fp, B_node, fe->total_num_nodes, "B_node");

    set_elem_types(fe, ned, elem_type);
    
    BB_vtk_write_point_vals_scalar(fp, elem_type, fe->total_num_nodes, "elem_type");
    fclose(fp);

    BB_std_free_2d_double(B_cell, fe->total_num_elems, 3);
    BB_std_free_2d_double(B_node, fe->total_num_nodes, 3);
}

void copy_Aphi_to_V_phi_time2(
    BBFE_DATA* fe,
    NEDELEC* ned,
    double * Aphi,
    double * V,
    double * phi,
    const int total_num_elems)
{

    for(int i = 0; i < total_num_elems; i++){
        for(int j = 0; j < ned->local_num_edges; j++){
            V[ned->nedelec_conn[i][j]] = Aphi[ned->nedelec_conn[i][j]];
        }
    }

    for(int i = 0; i < total_num_elems; i++){
        for(int j = 0; j < fe->local_num_nodes; j++){
            phi[fe->conn[i][j]] = Aphi[fe->conn[i][j]];
        }
    }

}

void compute_nedelec_edge_coords(
    BBFE_DATA* fe,
    NEDELEC* ned,
    int total_elems)
{
    int num_edges_per_elem = 0;
    const int (*curr_edge_conn)[2] = NULL;

    if (fe->local_num_nodes == 8) {
        num_edges_per_elem = 12;
        curr_edge_conn = hex_edge_conn;
    } else if (fe->local_num_nodes == 4) {
        num_edges_per_elem = 6;
        curr_edge_conn = tet_edge_conn;
    } else {
        fprintf(stderr, "Error: Unsupported element type (nodes=%d)\n", fe->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    ned->nedelec_coords = BB_std_calloc_2d_double(ned->nedelec_coords, fe->total_num_nodes, 3);

    for (int e = 0; e < total_elems; e++) {
        for (int n = 0; n < num_edges_per_elem; n++) {
            int local_node1 = curr_edge_conn[n][0];
            int local_node2 = curr_edge_conn[n][1];

            int global_node1 = fe->conn[e][local_node1];
            int global_node2 = fe->conn[e][local_node2];

            // 中点を計算
            double midpoint[3];
            for (int d = 0; d < 3; d++) {
                midpoint[d] = 0.5 * (fe->x[global_node1][d] + fe->x[global_node2][d]);
            }

            int global_edge = ned->nedelec_conn[e][n];
            for (int d = 0; d < 3; d++) {
                ned->nedelec_coords[global_edge][d] = midpoint[d];
            }
        }
    }
}

typedef struct {
    int g0, g1, g2;   /* sorted */
    int e;            /* element id */
    int f;            /* local face id 0..3 */
} FaceRec;

static inline void sort3_int(int* a, int* b, int* c){
    int x=*a,y=*b,z=*c,t;
    if(x>y){t=x;x=y;y=t;}
    if(y>z){t=y;y=z;z=t;}
    if(x>y){t=x;x=y;y=t;}
    *a=x; *b=y; *c=z;
}

static int cmp_face(const void* A, const void* B){
    const FaceRec* a=(const FaceRec*)A;
    const FaceRec* b=(const FaceRec*)B;
    if(a->g0!=b->g0) return (a->g0<b->g0)?-1:1;
    if(a->g1!=b->g1) return (a->g1<b->g1)?-1:1;
    if(a->g2!=b->g2) return (a->g2<b->g2)?-1:1;
    return 0;
}

static int find_local_edge_id_tet(int lnA, int lnB){
    for(int ed=0; ed<6; ++ed){
        int a = tet_edge_conn[ed][0];
        int b = tet_edge_conn[ed][1];
        if((a==lnA && b==lnB) || (a==lnB && b==lnA)) return ed;
    }
    return -1;
}

void build_dirichlet_edge_mask_from_boundary_faces_tet(
    const BBFE_DATA* fe,
    const BBFE_BC* bc,
    const NEDELEC* ned,
    bool* is_dir_edge,          /* size = num_global_edges, zero-initialized */
    int  num_global_edges
){
    const int M = fe->total_num_elems;
    const int nen = fe->local_num_nodes;
    if(nen != 4){
        fprintf(stderr, "This helper is for Tet4 only (nen=%d)\n", nen);
        return;
    }

    /* 4 faces per tet */
    static const int face_nodes[4][3] = {
        {0,1,2},{0,1,3},{0,2,3},{1,2,3}
    };

    FaceRec* faces = (FaceRec*)malloc(sizeof(FaceRec) * (size_t)M * 4);
    if(!faces){ perror("malloc"); return; }

    int idx=0;
    for(int e=0; e<M; ++e){
        for(int f=0; f<4; ++f){
            int gA = fe->conn[e][ face_nodes[f][0] ];
            int gB = fe->conn[e][ face_nodes[f][1] ];
            int gC = fe->conn[e][ face_nodes[f][2] ];
            int a=gA,b=gB,c=gC;
            sort3_int(&a,&b,&c);
            faces[idx++] = (FaceRec){ .g0=a,.g1=b,.g2=c,.e=e,.f=f };
        }
    }

    qsort(faces, (size_t)idx, sizeof(FaceRec), cmp_face);

    /* 同一面の出現回数を見る */
    for(int i=0; i<idx; ){
        int j=i+1;
        while(j<idx && cmp_face(&faces[i], &faces[j])==0) j++;

        const int count = j - i;
        if(count == 1){
            /* 境界面。OuterBoundaryかどうか判定：3節点すべてDirichletノードなら採用 */
            const int e = faces[i].e;
            const int f = faces[i].f;

            int ln0 = face_nodes[f][0];
            int ln1 = face_nodes[f][1];
            int ln2 = face_nodes[f][2];

            int gn0 = fe->conn[e][ln0];
            int gn1 = fe->conn[e][ln1];
            int gn2 = fe->conn[e][ln2];

            if(bc->D_bc_exists[gn0] && bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2]){
                /* この境界三角形の3辺だけをDirichlet対象にする */
                const int lned01 = find_local_edge_id_tet(ln0, ln1);
                const int lned12 = find_local_edge_id_tet(ln1, ln2);
                const int lned20 = find_local_edge_id_tet(ln2, ln0);

                const int ledges[3] = { lned01, lned12, lned20 };
                for(int k=0;k<3;++k){
                    int led = ledges[k];
                    if(led < 0) continue;
                    int gid = ned->nedelec_conn[e][led];
                    if(0 <= gid && gid < num_global_edges){
                        is_dir_edge[gid] = 1;
                    }
                }
            }
        }

        i = j;
    }

    free(faces);
}


void BBFE_mag_set_basis(
        BBFE_DATA*     fe,
        BBFE_BASIS*    basis,
        NEDELEC*       ned,
        int            local_num_nodes,
        int            num_integ_points_each_axis)
{
    (void)num_integ_points_each_axis;

    int num_edges_per_elem = 0;
    if (local_num_nodes == 8) {
        num_edges_per_elem = 12; /* Hex */
        ned->local_num_edges = 12;
    } else if (local_num_nodes == 4) {
        num_edges_per_elem = 6;  /* Tet */
        ned->local_num_edges = 6;
    } else {
        fprintf(stderr, "Error: unsupported local_num_nodes=%d\n", local_num_nodes);
        exit(EXIT_FAILURE);
    }

    ned->N_edge      = BB_std_calloc_4d_double(ned->N_edge,
                           fe->total_num_elems, basis->num_integ_points,
                           num_edges_per_elem, 3);

    ned->curl_N_edge = BB_std_calloc_4d_double(ned->curl_N_edge,
                           fe->total_num_elems, basis->num_integ_points,
                           num_edges_per_elem, 3);

    double J_inv[3][3];

    for (int e = 0; e < fe->total_num_elems; ++e) {
        for (int ip = 0; ip < basis->num_integ_points; ++ip) {

            BB_calc_mat3d_inverse(fe->geo[e][ip].J, fe->geo[e][ip].Jacobian, J_inv);

            if (local_num_nodes == 8) {
                BBFE_std_shapefunc_hex1st_nedelec_get_val(
                    basis->integ_point[ip],
                    ned->N_edge[e][ip],
                    J_inv);

                BBFE_std_shapefunc_hex1st_nedelec_get_curl(
                    basis->integ_point[ip],
                    ned->curl_N_edge[e][ip],
                    fe->geo[e][ip].J,
                    fe->geo[e][ip].Jacobian);
            }
            else {
                BBFE_std_shapefunc_tet1st_nedelec_get_val(
                    basis->integ_point[ip],
                    ned->N_edge[e][ip],
                    J_inv);

                BBFE_std_shapefunc_tet1st_nedelec_get_curl(
                    basis->integ_point[ip],
                    ned->curl_N_edge[e][ip],
                    fe->geo[e][ip].J,
                    fe->geo[e][ip].Jacobian);
            }
        }
    }
}


void read_edge_sign(
    BBFE_DATA *fe,
    NEDELEC *ned,
    const char *filename,
    const char *directory,
    int num_integ_points)
{
    FILE* fp;
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    int graph_ndof;
    int tmp;
    int total_num_elems;
    int local_num_nodes;

    const char* fname;
	fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, filename);

    fp = BBFE_sys_read_fopen(fp, fname, directory);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%s", char_n_internal);
    fscanf(fp, "%d %d", &(total_num_elems), &(local_num_nodes));

    ned->edge_sign = BB_std_calloc_2d_int(ned->edge_sign, total_num_elems, local_num_nodes);

    for (int e = 0; e < total_num_elems; e++)
    {
        for (int i = 0; i < fe->local_num_nodes; i++)
        {
            fscanf(fp, "%d", &(tmp));
        }
        for (int i = 0; i < ned->local_num_edges; i++)
        {
            fscanf(fp, "%d", &(ned->edge_sign[e][i]));
        }
    }

    fclose(fp);
}


void ROM_std_hlpod_set_nonzero_pattern_bcsr(
    MONOLIS*     	monolis,
    const char*     label,
    const char*		directory)
{
    const char* fname;
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);

    fp = BBFE_sys_read_fopen(fp, fname, directory);

    fscanf(fp, "%d", &(num_nodes));

    for(int e = 0; e < num_nodes; e++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes) );

        for(int i = 0; i < num_adj_nodes; i++) {
            fscanf(fp, "%d", &(tmp));
        }
        sum += num_adj_nodes;
    }
    fclose(fp);
    int* index = BB_std_calloc_1d_int(index, num_nodes + 1);
    int* item = BB_std_calloc_1d_int(item, sum);

    sum = 0;
    int index_sum = 0;
    index[0] = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);
    fp = BBFE_sys_read_fopen(fp, fname, directory);
    fscanf(fp, "%d", &(num_nodes));

    for(int i = 0; i < num_nodes; i++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes));

        for(int j = 0; j < num_adj_nodes; j++) {
            fscanf(fp, "%d", &(tmp));
            item[index_sum] = tmp;
            index_sum++;
        }
        sum += num_adj_nodes;
        index[i + 1] = sum;
    }	
    fclose(fp);

    monolis_get_nonzero_pattern_by_nodal_graph_R(
        monolis,
        num_nodes,
        1,
        index,
        item);

}

void set_elem_prop(
		BBFE_DATA*      fe,
        NEDELEC*        ned,
		const char*     directory,
        const char*     filename)
{
    ned->elem_prop = BB_std_calloc_1d_int(ned->elem_prop, fe->total_num_elems);

    FILE* fp;
	fp = BBFE_sys_read_fopen(fp, filename, directory);

    int total_num_elems;
    int tmp;
    char id[BUFFER_SIZE];

    fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &total_num_elems, &tmp);

	for(int e=0; e< total_num_elems; e++) {
		fscanf(fp, "%d", &(ned->elem_prop[e]));
	}

	fclose(fp);
}

void read_elem_types(
    BBFE_DATA*    fe,
    NEDELEC *ned,
    const char*     label,
    const char*		directory)
{
    const char* fname;
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;


    fp = BBFE_sys_read_fopen(fp, label, directory);
    fscanf(fp, "%d %d", &(tmp), &(fe->local_num_nodes));
    printf("%s elem types: %d\n", CODENAME, fe->local_num_nodes);

    fclose(fp);

    if (fe->local_num_nodes == 8) {
        ned->local_num_edges = 12;
    } else if (fe->local_num_nodes == 4) {
        ned->local_num_edges = 6;
    } else {
        fprintf(stderr, "Error: unsupported local_num_nodes=%d\n", fe->local_num_nodes);
        exit(EXIT_FAILURE);
    }

}

void read_num_nodes(
    BBFE_DATA*    fe,
    const char*     label,
    const char*		directory)
{
    const char* fname;
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);

    fp = BBFE_sys_read_fopen(fp, fname, directory);

    fscanf(fp, "%d", &(fe->total_num_nodes));


    fclose(fp);

}

void read_connectivity_graph_lag_nedelec(
    BBFE_DATA *fe,
    NEDELEC *ned,
    const char *filename,
    const char *directory,
    int num_integ_points)
{
    FILE* fp;
    char fname_n_internal_graph[BUFFER_SIZE];
    char char_n_internal[BUFFER_SIZE];
    int graph_ndof;
    int tmp;
    int total_num_elems;

    const char* fname;
	fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, filename);

    fp = BBFE_sys_read_fopen(fp, fname, directory);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &(total_num_elems));
    printf("%s Num. elements: %d\n", CODENAME, total_num_elems);
    int local_num_nodes = ned->local_num_edges; // 1st-order rectangle
    fe->total_num_elems = total_num_elems;
    printf("%s Num. elements: %d\n", CODENAME, fe->total_num_elems);
    BBFE_sys_memory_allocation_elem(fe, num_integ_points, 3);
    ned->nedelec_conn = BB_std_calloc_2d_int(ned->nedelec_conn, total_num_elems, local_num_nodes);

    if (total_num_elems < 0)
    {
        exit(EXIT_FAILURE);
    }
    if (total_num_elems == 0)
    {
        return;
    }

    int id = 0; // 　elem id (not used)
    for (int e = 0; e < total_num_elems; e++)
    {
        if (fscanf(fp, "%d ", &id) != 1)
        {
            exit(EXIT_FAILURE);
        }
        if (fscanf(fp, "%d ", &(tmp)) != 1)
        {
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < (fe->local_num_nodes); i++)
        {
            if(fscanf(fp, "%d", &(fe->conn[e][i])) != 1)
            {
                exit(EXIT_FAILURE);
            }
        }
        for (int i = 0; i < (local_num_nodes); i++)
        {
            if(fscanf(fp, "%d", &(ned->nedelec_conn[e][i])) != 1)
            {
                exit(EXIT_FAILURE);
            }
        }
    }


    fclose(fp);
}


void monowrap_solve_C(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double _Complex*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon)
{
    monolis_set_method   (monolis, solver_type);
    monolis_set_precond  (monolis, precond_type);
    monolis_set_maxiter  (monolis, num_max_iters);
    monolis_set_tolerance(monolis, epsilon);
    //monolis_show_iterlog (monolis, false);

    monolis_solve_C(
            monolis,
            monolis_com,
            monolis->mat.C.B,
            sol_vec);
}

void monowrap_solve_R(
        MONOLIS*      monolis,
        MONOLIS_COM*  monolis_com,
        double*       sol_vec,
        const int     solver_type,
        const int     precond_type,
        const int     num_max_iters,
        const double  epsilon)
{
    monolis_set_method   (monolis, solver_type);
    monolis_set_precond  (monolis, precond_type);
    monolis_set_maxiter  (monolis, num_max_iters);
    monolis_set_tolerance(monolis, epsilon);
    monolis_show_iterlog (monolis, true);

    monolis_solve_R(
            monolis,
            monolis_com,
            monolis->mat.R.B,
            sol_vec);
}

void BBFE_fluid_sups_read_Dirichlet_bc_NR(
                BBFE_BC*     bc,
                const char*  filename,
                const char*  directory,
                const int    total_num_nodes,
                const int    block_size)
{
        bc->total_num_nodes = total_num_nodes;
        bc->block_size      = block_size;

        BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
        int n = total_num_nodes * bc->block_size;

        for(int i=0; i<n; i++) {
                bc->D_bc_exists[i]   = false;
                bc->imposed_D_val[i] = 0.0;
        }

        FILE* fp;
        fp = BBFE_sys_read_fopen(fp, filename, directory);
        if( fp == NULL ) {
                printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
                                CODENAME, filename);
                return;
        }

        int tmp;
        BB_std_scan_line(&fp, BUFFER_SIZE,
                        "%d %d", &(bc->num_D_bcs), &(tmp));
        printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

        for(int i=0; i<(bc->num_D_bcs); i++) {
                int node_id;  int block_id;  double val;
                BB_std_scan_line(&fp, BUFFER_SIZE,
                                "%d %d %lf", &node_id, &block_id, &val);

                int index = (bc->block_size)*node_id + block_id;
                bc->D_bc_exists[ index ]   = true;
                bc->imposed_D_val[ index ] = 0.0;
        }

        fclose(fp);
}


void ROM_std_hlpod_set_nonzero_pattern_bcsr_C(
    MONOLIS*     	monolis,
    const char*     label,
    const char*		directory)
{
    const char* fname;
    FILE* fp;

    int num_nodes;
    int num_adj_nodes;
    int tmp;
    int sum = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);

    fp = BBFE_sys_read_fopen(fp, fname, directory);

    fscanf(fp, "%d", &(num_nodes));

    for(int e = 0; e < num_nodes; e++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes) );

        for(int i = 0; i < num_adj_nodes; i++) {
            fscanf(fp, "%d", &(tmp));
        }
        sum += num_adj_nodes;
    }
    fclose(fp);
    int* index = BB_std_calloc_1d_int(index, num_nodes + 1);
    int* item = BB_std_calloc_1d_int(item, sum);

    sum = 0;
    int index_sum = 0;
    index[0] = 0;

    fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, label);
    fp = BBFE_sys_read_fopen(fp, fname, directory);
    fscanf(fp, "%d", &(num_nodes));

    for(int i = 0; i < num_nodes; i++) {
        fscanf(fp, "%d", &(tmp));
        fscanf(fp, "%d", &(num_adj_nodes));

        for(int j = 0; j < num_adj_nodes; j++) {
            fscanf(fp, "%d", &(tmp));
            item[index_sum] = tmp;
            index_sum++;
        }
        sum += num_adj_nodes;
        index[i + 1] = sum;
    }	
    fclose(fp);

    monolis_get_nonzero_pattern_by_nodal_graph_C(
        monolis,
        num_nodes,
        1,
        index,
        item);

}


void BBFE_mag_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
        NEDELEC*      ned,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution)
{
	BB_calc_void();

	int n_axis = num_integ_points_each_axis;
	const char* filename;

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_NODE);
	BBFE_sys_read_node(
			fe,
			filename,
			directory);
    
    read_num_nodes(
			fe,
            "graph.dat",
            directory);
    
    read_elem_types(
			fe,
        	ned,
            "elem.dat",
            directory);

    read_connectivity_graph_lag_nedelec(
			fe,
			ned,
            "graph_nedelec_elem.dat",
			directory,
			n_axis*n_axis*n_axis);

    read_edge_sign(
			fe,
			ned,
            "nedelec_edge_sign.dat",
			directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			fe->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_convdiff_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);

	monolis_initialize(monolis);

	monolis_com_initialize_by_parted_files(
			monolis_com,
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"graph.dat");

    ROM_std_hlpod_set_nonzero_pattern_bcsr(
            monolis,
            "graph.dat",
            directory);

}


void BBFE_mag_pre_C(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
        NEDELEC*      ned,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution)
{
	BB_calc_void();

	int n_axis = num_integ_points_each_axis;
	const char* filename;

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_NODE);
	BBFE_sys_read_node(
			fe,
			filename,
			directory);
    
    read_num_nodes(
			fe,
            "graph.dat",
            directory);
    
    read_elem_types(
			fe,
        	ned,
            "elem.dat",
            directory);

    read_connectivity_graph_lag_nedelec(
			fe,
			ned,
            "graph_nedelec_elem.dat",
			directory,
			n_axis*n_axis*n_axis);

    read_edge_sign(
			fe,
			ned,
            "nedelec_edge_sign.dat",
			directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			fe->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_convdiff_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);

	monolis_initialize(monolis);

	monolis_com_initialize_by_parted_files(
			monolis_com,
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"graph.dat");

    ROM_std_hlpod_set_nonzero_pattern_bcsr_C(
            monolis,
            "graph.dat",
            directory);

}