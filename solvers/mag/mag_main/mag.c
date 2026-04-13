
#include "math.h"

#include "./../mag_core/convdiff_core.h"
#include "./../mag_core/nedelec_core.h"

#include "./../mag_core/elemmat.h"
#include "./../mag_core/shapefunc.h"
#include "./../mag_core/std.h"

#include "mag_dataset.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 3;
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


static int cmp_edgecand(const void *pa, const void *pb) {
    const EdgeCand *A = (const EdgeCand*)pa;
    const EdgeCand *B = (const EdgeCand*)pb;
    if (A->a != B->a) return (A->a < B->a) ? -1 : 1;
    if (A->b != B->b) return (A->b < B->b) ? -1 : 1;
    if (A->elem != B->elem) return (A->elem < B->elem) ? -1 : 1;
    return (A->ledge < B->ledge) ? -1 : (A->ledge > B->ledge);
}


void output_set_elems_nedelec_unstructured(
    BBFE_DATA   *fe,
    NEDELEC     *ned,
    const char* directory)
{
    const int num_elems = fe->total_num_elems;

    (void)directory;

    /* 念のため初期化（呼び出し側が安全になる） */
    ned->num_edges  = 0;

    /* --- 要素タイプの判定と設定 --- */
    int num_edges_per_elem = 0;
    const int (*curr_edge_conn)[2] = NULL;

    if (fe->local_num_nodes == 8) {
        num_edges_per_elem = 12;
        curr_edge_conn = hex_edge_conn;
        printf("Generating Nédélec edge connectivity for Hexahedra (%d elems)...\n", num_elems);
    }
    else if (fe->local_num_nodes == 4) {
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
    const size_t M = (size_t)num_elems * (size_t)num_edges_per_elem;
    EdgeCand *cand = (EdgeCand*)malloc(sizeof(EdgeCand) * M);
    if (!cand) { fprintf(stderr, "alloc failed (cand)\n"); exit(EXIT_FAILURE); }

    size_t t = 0;
    for (int e = 0; e < num_elems; ++e) {
        for (int le = 0; le < num_edges_per_elem; ++le) {
            const int u = fe->conn[e][ curr_edge_conn[le][0] ];
            const int v = fe->conn[e][ curr_edge_conn[le][1] ];

            const int a = (u < v) ? u : v;
            const int b = (u < v) ? v : u;

            cand[t++] = (EdgeCand){ .a=a, .b=b, .u=u, .v=v, .elem=e, .ledge=le };
        }
    }

    /* 無向キー(a,b)でソート */
    qsort(cand, M, sizeof(EdgeCand), cmp_edgecand);

    /* 最大 M 本のエッジ用に一旦確保 */
    int (*edge_nodes)[2] = (int (*)[2])malloc(sizeof(int[2]) * M);
    if (!edge_nodes) { fprintf(stderr, "alloc failed (edge_nodes)\n"); free(cand); exit(EXIT_FAILURE); }

    int num_edges = 0;

    /* 走査しながらユニーク化 & グローバルID付与 */
    for (size_t i = 0; i < M; ) {
        size_t j = i + 1;
        while (j < M && cand[j].a == cand[i].a && cand[j].b == cand[i].b) ++j;

        const int eid = num_edges++;
        edge_nodes[eid][0] = cand[i].a; /* min */
        edge_nodes[eid][1] = cand[i].b; /* max */

        for (size_t k = i; k < j; ++k) {
            const int sign = (cand[k].u == cand[k].a) ? +1 : -1;
            ned->nedelec_conn[cand[k].elem][cand[k].ledge] = eid;
            ned->edge_sign[cand[k].elem][cand[k].ledge]    = sign;
        }

        i = j;
    }

    /* edge_nodes を縮小して ned に保持 */
    ned->num_edges = num_edges;

    if (num_edges > 0) {
        int (*tmp_nodes)[2] = (int (*)[2])realloc(edge_nodes, sizeof(int[2]) * (size_t)num_edges);
        if (!tmp_nodes) {
            fprintf(stderr, "realloc failed (edge_nodes)\n");
            free(edge_nodes);
            free(cand);
            exit(EXIT_FAILURE);
        }
        edge_nodes = tmp_nodes;
    } else {
        free(edge_nodes);
        edge_nodes = NULL;
    }

    free(cand);

    printf("Done. Unique global edges = %d\n", ned->num_edges);
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
    fprintf(fp,"%d %d\n", fe->total_num_elems, fe->local_num_nodes + ned->local_num_edges);

    // データ書き出し
    for(int e = 0; e < fe->total_num_elems; e++){
        // ノードID
        for(int n = 0; n < fe->local_num_nodes; n++){
            fprintf(fp, "%d ", fe->conn[e][n]);
        }
        // エッジID (総ノード数をオフセットとして加算)
        for(int n = 0; n < ned->local_num_edges; n++){
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + fe->total_num_nodes);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);


    fp = BBFE_sys_write_add_fopen(fp, "nedelec_edge_sign.dat", directory);
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
}

/* * 関数2: 汎用版 (Edges Only 出力)
 * ファイル出力フォーマット:
 * [Edge ID 1 + Offset] ... [Edge ID M + Offset]
 */
void output_set_elems_nedelec_unstructured_only(
    BBFE_DATA   *fe,
    NEDELEC     *ned,
    const char* fname,
    const char* directory)
{
    const int num_elems = fe->total_num_elems;
    // const int num_nodes = fe->total_num_nodes;

    (void)directory;

    /* --- 要素タイプ判定 --- */
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

    ned->local_num_edges = num_edges_per_elem;

    /* --- メモリ確保 --- */
    ned->nedelec_conn = BB_std_calloc_2d_int(ned->nedelec_conn, num_elems, num_edges_per_elem);

    /* --- 候補リスト作成 --- */
    const size_t M = (size_t)num_elems * num_edges_per_elem;
    EdgeCand *cand = (EdgeCand*)malloc(sizeof(EdgeCand) * M);
    if (!cand) { fprintf(stderr, "alloc failed\n"); exit(EXIT_FAILURE); }

    size_t t = 0;
    for (int e = 0; e < num_elems; ++e) {
        for (int le = 0; le < num_edges_per_elem; ++le) {
            int u = fe->conn[e][ curr_edge_conn[le][0] ];
            int v = fe->conn[e][ curr_edge_conn[le][1] ];
            int a = (u < v) ? u : v;
            int b = (u < v) ? v : u;
            cand[t++] = (EdgeCand){ .a=a, .b=b, .u=u, .v=v, .elem=e, .ledge=le };
        }
    }

    qsort(cand, M, sizeof(EdgeCand), cmp_edgecand);

    int (*edge_nodes)[2] = (int (*)[2])malloc(sizeof(int[2]) * M);
    if (!edge_nodes) { fprintf(stderr, "alloc failed\n"); exit(EXIT_FAILURE); }

    int num_edges = 0;

    /* --- ユニーク化とID付与 --- */
    for (size_t i = 0; i < M; ) {
        size_t j = i + 1;
        while (j < M && cand[j].a == cand[i].a && cand[j].b == cand[i].b) ++j;

        const int eid = num_edges++;
        edge_nodes[eid][0] = cand[i].a;
        edge_nodes[eid][1] = cand[i].b;

        for (size_t k = i; k < j; ++k) {
            ned->nedelec_conn[cand[k].elem][cand[k].ledge] = eid;
        }
        i = j;
    }

    /* --- ファイル出力 (Edges Only) --- */
    FILE* fp = fopen(fname, "w");
    if (!fp) { perror("fopen"); exit(EXIT_FAILURE); }

    // ヘッダー: 要素数, 1要素あたりのエッジ数
    fprintf(fp, "%d %d\n", fe->total_num_elems, ned->local_num_edges);

    ned->num_edges = num_edges;
    
    if (num_edges > 0) {
        int (*tmp)[2] = (int (*)[2])realloc(edge_nodes, sizeof(int[2]) * num_edges);
        if (tmp) edge_nodes = tmp;
    } else {
        free(edge_nodes); edge_nodes = NULL;
    }

    if (edge_nodes) free(edge_nodes);
    free(cand);

    for(int e = 0; e < fe->total_num_elems; e++){
        for(int n = 0; n < ned->local_num_edges; n++){
            // ノード数オフセットを加算
            fprintf(fp, "%d ", ned->nedelec_conn[e][n] + fe->total_num_nodes);
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


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));
    
    BBFE_mag_set_basis(
        &(sys.fe),
        &(sys.basis),
        &(sys.ned),
        sys.fe.local_num_nodes,
        sys.vals.num_ip_each_axis);

    const char* filename;
/*
    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem_bool.dat");

    set_elem_prop(
            &(sys.fe),
            &(sys.ned),
            sys.cond.directory,
            filename);
*/

    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "nedelec_elem.dat");

    output_set_elems_nedelec_unstructured(
			&(sys.fe),
            &(sys.ned),
            filename,
            sys.cond.directory);

    filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "nedelec_elem_only.dat");

    output_set_elems_nedelec_unstructured_only(
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

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));

	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
