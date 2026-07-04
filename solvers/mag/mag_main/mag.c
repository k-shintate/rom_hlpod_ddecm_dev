
#include "math.h"

#include "./../mag_core/convdiff_core.h"
#include "./../mag_core/nedelec_core.h"

#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include "mag_dataset.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-10;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;


const double DELTA    = 1.0E-06;
const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";


#define TET2_NEDGE        6
#define TET2_NFACE        4
#define TET2_EDGE_MODE    2
#define TET2_FACE_MODE    2
#define TET2_NDOF         20
#define TET2_FACE_OFFSET  12

static const int tet_face_conn[4][3] = {
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
};

static const int tet_face_edge_lids[4][3] = {
    {0, 1, 2},
    {0, 4, 3},
    {2, 5, 3},
    {1, 5, 4}
};

typedef struct {
    int a, b, c;
    int elem;
    int lface;
} FaceCand2nd;

static inline void swap_int_2nd(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

static inline void sort3_int_2nd(int* a, int* b, int* c)
{
    if(*a > *b) swap_int_2nd(a, b);
    if(*b > *c) swap_int_2nd(b, c);
    if(*a > *b) swap_int_2nd(a, b);
}

static int cmp_facecand_2nd(const void* pa, const void* pb)
{
    const FaceCand2nd* A = (const FaceCand2nd*)pa;
    const FaceCand2nd* B = (const FaceCand2nd*)pb;

    if(A->a != B->a) return (A->a < B->a) ? -1 : 1;
    if(A->b != B->b) return (A->b < B->b) ? -1 : 1;
    if(A->c != B->c) return (A->c < B->c) ? -1 : 1;
    if(A->elem != B->elem) return (A->elem < B->elem) ? -1 : 1;
    return (A->lface < B->lface) ? -1 : (A->lface > B->lface);
}


static int cmp_edgecand(const void *pa, const void *pb) {
    const EdgeCand *A = (const EdgeCand*)pa;
    const EdgeCand *B = (const EdgeCand*)pb;
    if (A->a != B->a) return (A->a < B->a) ? -1 : 1;
    if (A->b != B->b) return (A->b < B->b) ? -1 : 1;
    if (A->elem != B->elem) return (A->elem < B->elem) ? -1 : 1;
    return (A->ledge < B->ledge) ? -1 : (A->ledge > B->ledge);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}


#define EDGE_ID(x)   (abs((x)) - 1)
#define EDGE_SIGN(x) ((x) > 0 ? +1 : -1)


void output_set_elems_nedelec_unstructured(
    BBFE_DATA   *fe,
    NEDELEC     *ned,
    const char* fname,
    const char* directory)
{
    const int num_elems = fe->total_num_elems;
    const int num_nodes = fe->total_num_nodes;

    (void)directory;

    /* --- 要素タイプの判定と設定 --- */
    int num_edges_per_elem = 0;
    const int (*curr_edge_conn)[2] = NULL;

    if (fe->local_num_nodes == 8) {
        /* Hexahedron */
        num_edges_per_elem = 12;
        curr_edge_conn = hex_edge_conn;
        printf("Generating Nédélec edge connectivity for Hexahedra (%d elems)...\n", num_elems);
    } 
    else if (fe->local_num_nodes == 4) {
        /* Tetrahedron */
        num_edges_per_elem = 6;
        curr_edge_conn = tet_edge_conn;
        printf("Generating Nédélec edge connectivity for Tetrahedra (%d elems)...\n", num_elems);
    } 
    else {
        fprintf(stderr, "Error: Unsupported element type (nodes=%d)\n", fe->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    /* 出力配列の確保：要素数 × 要素あたりのエッジ数 */
    ned->nedelec_conn = BB_std_calloc_2d_int(ned->nedelec_conn, num_elems, num_edges_per_elem);
    ned->edge_sign    = BB_std_calloc_2d_int(ned->edge_sign,    num_elems, num_edges_per_elem);

    /* 全要素の全ローカルエッジを候補として列挙 */
    const size_t M = (size_t)num_elems * num_edges_per_elem;
    EdgeCand *cand = (EdgeCand*)malloc(sizeof(EdgeCand) * M);
    if (!cand) { fprintf(stderr, "alloc failed (cand)\n"); exit(EXIT_FAILURE); }

    size_t t = 0;
    for (int e = 0; e < num_elems; ++e) {
        for (int le = 0; le < num_edges_per_elem; ++le) {
            /* 現在の要素タイプに対応したテーブルを参照 */
            int u = fe->conn[e][ curr_edge_conn[le][0] ];
            int v = fe->conn[e][ curr_edge_conn[le][1] ];
            
            int a = (u < v) ? u : v;
            int b = (u < v) ? v : u;
            
            cand[t++] = (EdgeCand){ .a=a, .b=b, .u=u, .v=v, .elem=e, .ledge=le };
        }
    }

    /* 無向キー(a,b)でソート */
    qsort(cand, M, sizeof(EdgeCand), cmp_edgecand);

    /* 最大 M 本のエッジ用に一旦確保 */
    int (*edge_nodes)[2] = (int (*)[2])malloc(sizeof(int[2]) * M);
    if (!edge_nodes) { fprintf(stderr, "alloc failed (edge_nodes)\n"); exit(EXIT_FAILURE); }

    int num_edges = 0;

    /* 走査しながらユニーク化 & グローバルID付与 */
    for (size_t i = 0; i < M; ) {
        size_t j = i + 1;
        while (j < M && cand[j].a == cand[i].a && cand[j].b == cand[i].b) ++j;

        const int eid = num_edges++;            /* 新しいグローバルエッジID */
        edge_nodes[eid][0] = cand[i].a;         /* min */
        edge_nodes[eid][1] = cand[i].b;         /* max */

        /* 同じキーに属する候補全てへ eid と符号を配布 */
        for (size_t k = i; k < j; ++k) {
            /* u->v が min->max (a->b) と一致すれば +1 */
            const int sign = (cand[k].u == cand[k].a) ? +1 : -1;
            
            ned->nedelec_conn[cand[k].elem][cand[k].ledge] = eid;
            ned->edge_sign[cand[k].elem][cand[k].ledge]    = sign;
        }

        i = j;
    }

    /* 最後に edge_nodes をジャストサイズに縮小して ned に保持 */
    ned->num_edges = num_edges;
    
    // エッジ数が0の場合(ありえないが念のため)の対応を含めるなら条件分岐推奨
    if (num_edges > 0) {
        int (*tmp_nodes)[2] = (int (*)[2])realloc(edge_nodes, sizeof(int[2]) * (size_t)num_edges);
        if (!tmp_nodes) {
            fprintf(stderr, "realloc failed (edge_nodes)\n");
            free(edge_nodes);
            free(cand);
            exit(EXIT_FAILURE);
        }
    } else {
        free(edge_nodes);
    }

    free(cand);

	FILE* fp;

	fp = fopen(fname, "w");
    printf("ned->part_num_nodes = %d fe->total_num_nodes = %d", ned->part_num_nodes, fe->total_num_nodes);

    fprintf(fp,"#elem_id\n");
    fprintf(fp,"%d 10\n", fe->total_num_elems);
    int num = 0;
    for(int e = 0; e < fe->total_num_elems; e++){
        for(int n = 0; n < fe->local_num_nodes; n++){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }
        for(int n = 0; n < ned->local_num_edges; n++){
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + fe->total_num_nodes);
        }
        fprintf(fp, "\n");
    }    
    fclose(fp);

    //#bool_elem
    //#78499 1
    //5
    fp = BBFE_sys_write_fopen(fp, "distval_elem_lagnedelec_elem.dat", directory);
    fprintf(fp,"#elem_id\n");
    fprintf(fp,"%d 10\n", fe->total_num_elems);
    num = 0;
    for(int e = 0; e < fe->total_num_elems; e++){
        //fprintf(fp, "%d %d ", e, 10);
        for(int n = 0; n < fe->local_num_nodes; n++){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }
        for(int n = 0; n < ned->local_num_edges; n++){
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + fe->total_num_nodes);
        }
        fprintf(fp, "\n");
    }
    
    ご対応いただきありがとうございます。
    根本くんと確認しましたが、マザボ側からは出力できるものの GPU からは出力されず、修理前と同様の状態でした。
    また、モニターを変える、ケーブルを変えるなどは一通り行いました。

    お手数ですが、
    動作テストの内容と初期設定の際に必要な設定やコマンド等があるかをお聞きできますでしょうか。








    fclose(fp);


    fp = BBFE_sys_write_fopen(fp, "node_coordinate_elem.dat", directory);

    fprintf(fp,"#elem_id\n");
    fprintf(fp,"%d %d\n", fe->total_num_elems, 3*(fe->local_num_nodes+ned->local_num_edges));
    //int num = 0;
    for(int e = 0; e < fe->total_num_elems; e++){
        for(int n = 0; n < fe->local_num_nodes; n++){
            fprintf(fp, "%lf %lf %lf ", fe->x[fe->conn[e][n]][0], fe->x[fe->conn[e][n]][1], fe->x[fe->conn[e][n]][2]);
        }
        for(int n = 0; n < ned->local_num_edges; n++){
            fprintf(fp, "0.0 0.0 0.0 ");
        }
        fprintf(fp, "\n");
    }    
    fclose(fp);


    fp = BBFE_sys_write_fopen(fp, "graph_nedelec_elem.dat", directory);
    fprintf(fp,"%d\n", fe->total_num_elems);

    num = 0;
    for(int e = 0; e < fe->total_num_elems; e++){
        // ノードID
        if(ned->elem_prop[e] == 3){
            fprintf(fp, "%d ", e);
            fprintf(fp, "%d ", 10);
            for(int n = 0; n < fe->local_num_nodes; n++){
                fprintf(fp, "%d ", ned->phi_conn[num][n]);
            }
            num++;
        }
        else{
            fprintf(fp, "%d ", e);
            fprintf(fp, "%d ", 6);
        }
        // エッジID (総ノード数をオフセットとして加算)
        for(int n = 0; n < ned->local_num_edges; n++){
            //fprintf(fp, "%d ", ned->nedelec_conn[e][n] + fe->total_num_nodes);
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + ned->part_num_nodes);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);


    fp = BBFE_sys_write_fopen(fp, "nedelec_edge_sign.dat", directory);
    fprintf(fp, "#edge_sign\n");
    fprintf(fp, "%d %d\n", num_elems, fe->local_num_nodes+ned->local_num_edges);

    for (int e = 0; e < num_elems; ++e) {
        for (int le = 0; le < fe->local_num_nodes; ++le) {
            fprintf(fp, "10 ");
        }
        for (int le = 0; le < ned->local_num_edges; ++le) {
            fprintf(fp, "%d ", ned->edge_sign[e][le]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);


    fp = BBFE_sys_write_fopen(fp, "nedelec_conn_all.dat", directory);
    fprintf(fp,"%d %d\n", fe->total_num_elems, ned->local_num_edges);
    for(int e = 0; e < fe->total_num_elems; e++){
        for(int n = 0; n < ned->local_num_edges; n++){
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + ned->part_num_nodes);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
}

void output_set_elems_nedelec_unstructured_2nd(
    BBFE_DATA*   fe,
    NEDELEC*     ned,
    const char*  fname,
    const char*  directory)
{
    if(fe->local_num_nodes != 4){
        fprintf(stderr,
            "output_set_elems_nedelec_unstructured_2nd: Tet4 only. local_num_nodes=%d\n",
            fe->local_num_nodes);
        exit(EXIT_FAILURE);
    }

    const int num_elems = fe->total_num_elems;
    const int edge_modes = TET2_EDGE_MODE;
    const int face_modes = TET2_FACE_MODE;

    ned->local_num_edges = TET2_NDOF;

    ned->nedelec_conn = BB_std_calloc_2d_int(
        ned->nedelec_conn,
        num_elems,
        TET2_NDOF
    );

    /*
     * 2次では edge_sign は計算本体では使わない。
     * ただしデバッグ確認用として edge_dir 6個分だけ出力する。
     */
    ned->edge_sign = BB_std_calloc_2d_int(
        ned->edge_sign,
        num_elems,
        TET2_NEDGE
    );

    /* ============================================================
     * 1. unique edge を作る
     * ============================================================ */
    const size_t Me = (size_t)num_elems * TET2_NEDGE;

    EdgeCand* ecand = (EdgeCand*)malloc(sizeof(EdgeCand) * Me);
    if(ecand == NULL){
        fprintf(stderr, "alloc failed: ecand\n");
        exit(EXIT_FAILURE);
    }

    size_t te = 0;

    for(int e = 0; e < num_elems; ++e){
        for(int le = 0; le < TET2_NEDGE; ++le){
            int u = fe->conn[e][tet_edge_conn[le][0]];
            int v = fe->conn[e][tet_edge_conn[le][1]];

            int a = (u < v) ? u : v;
            int b = (u < v) ? v : u;

            ecand[te++] = (EdgeCand){
                .a = a,
                .b = b,
                .u = u,
                .v = v,
                .elem = e,
                .ledge = le
            };
        }
    }

    qsort(ecand, Me, sizeof(EdgeCand), cmp_edgecand);

    int next_ned_dof = 0;
    int num_unique_edges = 0;

    for(size_t i = 0; i < Me; ){
        size_t j = i + 1;

        while(j < Me &&
              ecand[j].a == ecand[i].a &&
              ecand[j].b == ecand[i].b){
            ++j;
        }

        /*
         * この topological edge に 2 個の Nedelec DOF を割り当てる。
         */
        int gdof_edge_mode[2];
        for(int m = 0; m < edge_modes; ++m){
            gdof_edge_mode[m] = next_ned_dof++;
        }

        for(size_t k = i; k < j; ++k){
            const int e  = ecand[k].elem;
            const int le = ecand[k].ledge;

            const int sign = (ecand[k].u == ecand[k].a) ? +1 : -1;
            ned->edge_sign[e][le] = sign;

            for(int m = 0; m < edge_modes; ++m){
                const int lid = edge_modes * le + m;
                ned->nedelec_conn[e][lid] = gdof_edge_mode[m];
            }
        }

        ++num_unique_edges;
        i = j;
    }

    free(ecand);

    /* ============================================================
     * 2. unique face を作る
     * ============================================================ */
    const size_t Mf = (size_t)num_elems * TET2_NFACE;

    FaceCand2nd* fcand = (FaceCand2nd*)malloc(sizeof(FaceCand2nd) * Mf);
    if(fcand == NULL){
        fprintf(stderr, "alloc failed: fcand\n");
        exit(EXIT_FAILURE);
    }

    size_t tf = 0;

    for(int e = 0; e < num_elems; ++e){
        for(int lf = 0; lf < TET2_NFACE; ++lf){
            int a = fe->conn[e][tet_face_conn[lf][0]];
            int b = fe->conn[e][tet_face_conn[lf][1]];
            int c = fe->conn[e][tet_face_conn[lf][2]];

            sort3_int_2nd(&a, &b, &c);

            fcand[tf++] = (FaceCand2nd){
                .a = a,
                .b = b,
                .c = c,
                .elem = e,
                .lface = lf
            };
        }
    }

    qsort(fcand, Mf, sizeof(FaceCand2nd), cmp_facecand_2nd);

    int num_unique_faces = 0;

    for(size_t i = 0; i < Mf; ){
        size_t j = i + 1;

        while(j < Mf &&
              fcand[j].a == fcand[i].a &&
              fcand[j].b == fcand[i].b &&
              fcand[j].c == fcand[i].c){
            ++j;
        }

        /*
         * この topological face に 2 個の Nedelec DOF を割り当てる。
         */
        int gdof_face_mode[2];
        for(int m = 0; m < face_modes; ++m){
            gdof_face_mode[m] = next_ned_dof++;
        }

        for(size_t k = i; k < j; ++k){
            const int e  = fcand[k].elem;
            const int lf = fcand[k].lface;

            for(int m = 0; m < face_modes; ++m){
                const int lid = TET2_FACE_OFFSET + face_modes * lf + m;
                ned->nedelec_conn[e][lid] = gdof_face_mode[m];
            }
        }

        ++num_unique_faces;
        i = j;
    }

    free(fcand);

    /*
     * 既存コードの都合で num_edges を「Nedelec DOF 数」として使う。
     * 名前は不正確だが、構造体を変更しないならこれが最小変更。
     */
    ned->num_edges = next_ned_dof;

    printf("Tet2 Nedelec connectivity generated:\n");
    printf("  unique edges = %d\n", num_unique_edges);
    printf("  unique faces = %d\n", num_unique_faces);
    printf("  Nedelec DOFs = %d\n", next_ned_dof);

    /* ============================================================
     * 3. nedelec_elem_2nd.dat
     *
     * 1行:
     * node0 node1 node2 node3 ned0 ... ned19
     * graph_ndof = 4 + 20 = 24
     * ============================================================ */
    FILE* fp;

    fp = fopen(fname, "w");
    if(fp == NULL){
        perror(fname);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "#elem_id\n");
    fprintf(fp, "%d %d\n", fe->total_num_elems, fe->local_num_nodes + TET2_NDOF);

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int n = 0; n < fe->local_num_nodes; ++n){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }

        for(int j = 0; j < TET2_NDOF; ++j){
            fprintf(fp, "%d ", ned->nedelec_conn[e][j] + fe->total_num_nodes);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    /* ============================================================
     * 4. distval_elem_lagnedelec_elem_2nd.dat
     * ============================================================ */
    fp = BBFE_sys_write_fopen(fp, "distval_elem_lagnedelec_elem_2nd.dat", directory);

    fprintf(fp, "#elem_id\n");
    fprintf(fp, "%d %d\n", fe->total_num_elems, fe->local_num_nodes + TET2_NDOF);

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int n = 0; n < fe->local_num_nodes; ++n){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }

        for(int j = 0; j < TET2_NDOF; ++j){
            fprintf(fp, "%d ", ned->nedelec_conn[e][j] + fe->total_num_nodes);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    /* ============================================================
     * 5. node_coordinate_elem_2nd.dat
     *
     * 2次 Nedelec DOF は節点座標ではないので、ここでは 0 を出す。
     * 可視化用に edge/face 中心を出したい場合は別途拡張する。
     * ============================================================ */
    fp = BBFE_sys_write_fopen(fp, "node_coordinate_elem_2nd.dat", directory);

    fprintf(fp, "#elem_id\n");
    fprintf(fp, "%d %d\n",
        fe->total_num_elems,
        3 * (fe->local_num_nodes + TET2_NDOF));

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int n = 0; n < fe->local_num_nodes; ++n){
            const int gn = fe->conn[e][n];

            fprintf(fp, "%lf %lf %lf ",
                fe->x[gn][0],
                fe->x[gn][1],
                fe->x[gn][2]);
        }

        for(int j = 0; j < TET2_NDOF; ++j){
            fprintf(fp, "0.0 0.0 0.0 ");
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    /* ============================================================
     * 6. graph_nedelec_elem_2nd.dat
     *
     * elem_prop==4 のときは phi_conn 4個 + Nedelec 20個。
     * それ以外は Nedelec 20個だけ。
     * ============================================================ */
    fp = BBFE_sys_write_fopen(fp, "graph_nedelec_elem_2nd.dat", directory);

    fprintf(fp, "%d\n", fe->total_num_elems);

    int phi_elem_count = 0;

    for(int e = 0; e < fe->total_num_elems; ++e){
        fprintf(fp, "%d ", e);

        if(ned->elem_prop[e] == 4){
            fprintf(fp, "%d ", fe->local_num_nodes + TET2_NDOF);

            for(int n = 0; n < fe->local_num_nodes; ++n){
                fprintf(fp, "%d ", ned->phi_conn[phi_elem_count][n]);
            }

            ++phi_elem_count;
        }else{
            fprintf(fp, "%d ", TET2_NDOF);
        }

        for(int j = 0; j < TET2_NDOF; ++j){
            fprintf(fp, "%d ", ned->nedelec_conn[e][j] + ned->part_num_nodes);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    /* ============================================================
     * 7. nedelec_orient_2nd.dat
     *
     * 2次計算では edge_sign は掛けない方針だが、
     * デバッグ確認用に 6 edge の向きを保存する。
     * ============================================================ */
    fp = BBFE_sys_write_fopen(fp, "nedelec_orient_2nd.dat", directory);

    fprintf(fp, "#nedelec_orient_2nd\n");
    fprintf(fp, "%d %d %d %d\n",
        fe->total_num_elems,
        2,
        fe->local_num_nodes,
        TET2_NDOF);

    for(int e = 0; e < fe->total_num_elems; ++e){
        fprintf(fp, "%d ", e);

        for(int le = 0; le < TET2_NEDGE; ++le){
            fprintf(fp, "%d ", ned->edge_sign[e][le]);
        }

        /*
         * face permutation は今は基底評価側で global node order に揃える前提。
         * 将来 face transform を導入するなら 0..5 の順列IDを保存する。
         */
        for(int lf = 0; lf < TET2_NFACE; ++lf){
            fprintf(fp, "%d ", 0);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    /* ============================================================
     * 8. nedelec_conn_all_2nd.dat
     * ============================================================ */
    fp = BBFE_sys_write_fopen(fp, "nedelec_conn_all_2nd.dat", directory);

    fprintf(fp, "%d %d\n", fe->total_num_elems, TET2_NDOF);

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int j = 0; j < TET2_NDOF; ++j){
            fprintf(fp, "%d ", ned->nedelec_conn[e][j] + ned->part_num_nodes);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);

    fp = BBFE_sys_write_fopen(fp, "distval_elem_global_node_id_2nd.dat", directory);

    fprintf(fp, "#elem_global_node_id_2nd\n");
    fprintf(fp, "%d %d\n", fe->total_num_elems, fe->local_num_nodes);

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int n = 0; n < fe->local_num_nodes; ++n){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}


void compute_nedelec_edge_coords(
    BBFE_DATA* fe,
    NEDELEC* ned,
    int total_elems,
    int total_edges,
    const char* directory)
{
    ned->nedelec_coords = BB_std_calloc_2d_double(ned->nedelec_coords, total_edges, 3);

    for (int e = 0; e < total_elems; e++) {
        for (int n = 0; n < ned->local_num_edges; n++) {
            int local_node1 = hex_edge_conn[n][0];
            int local_node2 = hex_edge_conn[n][1];

            int global_node1 = fe->conn[e][local_node1];
            int global_node2 = fe->conn[e][local_node2];

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

void write_bc_ned(
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    const char* directory)
{

	FILE* fp;
    const int nen = fe->local_num_nodes;

    int is_dir_edge_n = ned->num_edges;
    bool* is_dir_edge = BB_std_calloc_1d_bool(is_dir_edge, is_dir_edge_n);
    build_dirichlet_edge_mask_from_boundary_faces_tet(fe, bc, ned, is_dir_edge, is_dir_edge_n);

    int n_local_edges = 0;
    const int (*edge_tbl)[2] = NULL;

    if(nen == 4){
        n_local_edges = 6;
        edge_tbl = tet_edge_conn;
    } else if(nen == 8){
        n_local_edges = 12;
        edge_tbl = hex_edge_conn;
    }

    int num = 0;

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int i = 0; i < n_local_edges; ++i){
            int gn1 = fe->conn[e][edge_tbl[i][0]];
            int gn2 = fe->conn[e][edge_tbl[i][1]];
            int ged = ned->nedelec_conn[e][i];

            if(!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;
            if(!is_dir_edge[ged]) continue;

            num++;
        }
    }

    fp = BBFE_sys_write_fopen(fp, "D_bc_ned.dat", directory);
    fprintf(fp, "%d 1\n", num);

    for(int e = 0; e < fe->total_num_elems; ++e){
        for(int i = 0; i < n_local_edges; ++i){
            int gn1 = fe->conn[e][edge_tbl[i][0]];
            int gn2 = fe->conn[e][edge_tbl[i][1]];
            int ged = ned->nedelec_conn[e][i];

            if(!(bc->D_bc_exists[gn1] && bc->D_bc_exists[gn2])) continue;
            if(!is_dir_edge[ged]) continue;

            int ged1 = ged + ned->part_num_nodes;

            fprintf(fp, "%d %d %lf\n", ged1, 0, 0.0);

        }
    }

    fclose(fp);

    BB_std_free_1d_bool(is_dir_edge, is_dir_edge_n);

}

void write_bc_ned_2nd(
    BBFE_DATA*  fe,
    BBFE_BC*    bc,
    NEDELEC*    ned,
    const char* directory)
{
    if(fe->local_num_nodes != 4){
        fprintf(stderr, "write_bc_ned_2nd: Tet4 only\n");
        exit(EXIT_FAILURE);
    }

    const int total_num_dof = ned->num_edges;

    bool* is_dir_dof = BB_std_calloc_1d_bool(
        is_dir_dof,
        total_num_dof
    );

    build_dirichlet_dof_mask_from_boundary_faces_tet_2nd(
        fe,
        bc,
        ned,
        is_dir_dof,
        total_num_dof
    );

    int num = 0;

    for(int gdof = 0; gdof < total_num_dof; ++gdof){
        if(is_dir_dof[gdof]) ++num;
    }

    FILE* fp = BBFE_sys_write_fopen(fp, "D_bc_ned_2nd.dat", directory);

    fprintf(fp, "%d 1\n", num);

    for(int gdof = 0; gdof < total_num_dof; ++gdof){
        if(is_dir_dof[gdof]){
            fprintf(fp, "%d %d %lf\n",
                gdof + ned->part_num_nodes,
                0,
                0.0);
        }
    }

    fclose(fp);

    BB_std_free_1d_bool(is_dir_dof, total_num_dof);
}



int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

    /*
     * Nedelec order switch for mesh/connectivity generation.
     *
     * Default is 1 because the 1st-order implementation has already been
     * verified. To test the 2nd-order implementation without editing this
     * file, run with an environment variable, for example:
     *
     *   NEDELEC_ORDER=2 ./your_executable <directory>
     *
     * You can also change DEFAULT_NEDELEC_ORDER below to 2 while developing.
     */
    const int DEFAULT_NEDELEC_ORDER = 2;
    int nedelec_order = DEFAULT_NEDELEC_ORDER;

    const char* env_order = getenv("NEDELEC_ORDER");
    if(env_order != NULL && env_order[0] != '\0'){
        nedelec_order = atoi(env_order);
    }

    if(nedelec_order != 1 && nedelec_order != 2){
        fprintf(stderr,
            "Unsupported NEDELEC_ORDER=%d. Use 1 or 2.\n",
            nedelec_order);
        exit(EXIT_FAILURE);
    }

    printf("Nedelec order for mesh generation: %d\n", nedelec_order);

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

    if(nedelec_order == 2 && sys.fe.local_num_nodes != 4){
        fprintf(stderr,
            "2nd-order Nedelec mesh generation supports Tet4 only. "
            "local_num_nodes=%d\n",
            sys.fe.local_num_nodes);
        exit(EXIT_FAILURE);
    }

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));

    if(nedelec_order == 1){
        /*
         * Verified 1st-order path. Keep this branch identical in behavior to
         * the existing implementation.
         */
        BBFE_mag_set_basis(
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            sys.fe.local_num_nodes,
            sys.vals.num_ip_each_axis);
    }else{
        /*
         * 2nd-order Tet Nedelec candidate basis.
         * This sets sys.ned.local_num_edges = TET2_NDOF (=20).
         */
        BBFE_mag_set_basis_2nd_pre(
            &(sys.fe),
            &(sys.basis),
            &(sys.ned),
            sys.fe.local_num_nodes,
            sys.vals.num_ip_each_axis);
    }

    const char* filename;
    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			BLOCK_SIZE);

    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem_bool.dat");

    set_elem_prop(
            &(sys.fe),
            &(sys.ned),
            sys.cond.directory,
            filename);

    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "material_node_mapping.dat");

    read_part_node_id(
            &(sys.fe),
            &(sys.ned),
            sys.cond.directory,
            filename);

    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "part_elem_conn.dat");

    read_part_elem_id(
            &(sys.fe),
            &(sys.ned),
            sys.cond.directory,
            filename);

    if(nedelec_order == 1){
        filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            "nedelec_elem.dat");

        output_set_elems_nedelec_unstructured(
            &(sys.fe),
            &(sys.ned),
            filename,
            sys.cond.directory);

        compute_nedelec_edge_coords(
            &(sys.fe),
            &(sys.ned),
            sys.fe.total_num_elems,
            sys.ned.num_edges,
            sys.cond.directory);

        write_bc_ned(
            &(sys.fe),
            &(sys.bc),
            &(sys.ned),
            sys.cond.directory);
    }else{
        filename = monolis_get_global_input_file_name(
            MONOLIS_DEFAULT_TOP_DIR,
            MONOLIS_DEFAULT_PART_DIR,
            "nedelec_elem_2nd.dat");

        output_set_elems_nedelec_unstructured_2nd(
            &(sys.fe),
            &(sys.ned),
            filename,
            sys.cond.directory);

        compute_nedelec_edge_coords_2nd(
            &(sys.fe),
            &(sys.ned),
            sys.fe.total_num_elems,
            sys.ned.num_edges,
            sys.cond.directory);

        write_bc_ned_2nd(
            &(sys.fe),
            &(sys.bc),
            &(sys.ned),
            sys.cond.directory);
    }

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
