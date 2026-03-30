/* merge_test_0based.c
 *
 * Fortran版 merge_test を C に移植し、C内部では「節点番号を 0 スタート」で扱う版。
 *
 * 方針（0スタート化）:
 * - bc/load ファイルから読んだ節点番号（元は 1..nnode）を、読み込み直後に -1 して 0..nnode-1 に変換して保持する
 * - merged 後に復元して出力する ibc_out / iload_out の節点番号も 0..nnode-1 で出力する（=0スタート出力）
 * - node.dat は節点番号がファイルに無い（行順）ので、値配列の添字が 0 始まりになるだけ（データ自体は同じ）
 *
 * 注意:
 * - monolis/gedatsu の C API の正確なシグネチャは環境依存です。
 *   ここでは「典型的な C 側の引数形」を仮定して書いています。
 *   もしコンパイルが通らない場合は、該当 API のプロトタイプに合わせて引数型を調整してください。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "monolis_utils.h"
#include "gedatsu_def_graph_c.h"
#include "gedatsu_graph_merge_c.h"
#include "monolis_std_list_c.h"

#ifndef MONOLIS_CHARLEN
#define MONOLIS_CHARLEN 256
#endif

#define CHECK_ALLOC(p) do { if(!(p)) { \
  fprintf(stderr,"alloc failed %s:%d\n",__FILE__,__LINE__); exit(1);} } while(0)

typedef struct {
  int node_id;     /* 0-based node id */
  int old_index;   /* 元の j (0-based) */
} PairNodeIdx;

/* Fortran 側の kint/kdouble 相当（環境で違うならここを合わせる） */
#ifndef kint
#define kint int
#endif

#ifndef kdouble
#define kdouble double
#endif

#define CHECK_ALLOC(p) do { \
  if(!(p)) { \
    fprintf(stderr, "alloc failed %s:%d\n", __FILE__, __LINE__); \
    exit(1); \
  } \
} while(0)

/* ============================================================
 * graph フォーマット
 *
 * 入力:
 *  1行目: n_vertex
 *  以降 n_vertex 行:
 *    vertex_id  in  item_1 ... item_in
 *
 * 出力も同形式
 *
 * CSR:
 *  index は長さ n_vertex+1, 0-based, index[0]=0
 *  i 行目の隣接は item[index[i] ... index[i+1)-1]
 * ============================================================ */
void monolis_input_graph(const char* fname,
                         kint* n_vertex,
                         kint** vertex_id,
                         kint** index,
                         kint** item)
{
  FILE* fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s\n", fname);
    exit(1);
  }

  if(fscanf(fp, "%d", n_vertex) != 1){
    fprintf(stderr, "failed read n_vertex from %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(*n_vertex < 0){
    fprintf(stderr, "invalid n_vertex=%d in %s\n", *n_vertex, fname);
    fclose(fp);
    exit(1);
  }

  /* 1st pass: count nz */
  kint nz = 0;
  for(kint i=0; i<*n_vertex; ++i){
    kint tmp_id, in;
    if(fscanf(fp, "%d %d", &tmp_id, &in) != 2){
      fprintf(stderr, "failed read header of vertex line %d in %s\n", (int)i, fname);
      fclose(fp);
      exit(1);
    }
    if(in < 0){
      fprintf(stderr, "invalid degree=%d at line %d in %s\n", (int)in, (int)i, fname);
      fclose(fp);
      exit(1);
    }
    nz += in;

    /* skip items */
    for(kint j=0; j<in; ++j){
      kint dummy;
      if(fscanf(fp, "%d", &dummy) != 1){
        fprintf(stderr, "failed skip item at line %d in %s\n", (int)i, fname);
        fclose(fp);
        exit(1);
      }
    }
  }

  fclose(fp);

  *vertex_id = (kint*)malloc(sizeof(kint) * (size_t)(*n_vertex));
  *index     = (kint*)malloc(sizeof(kint) * (size_t)(*n_vertex + 1));
  *item      = (kint*)malloc(sizeof(kint) * (size_t)nz);
  CHECK_ALLOC(*vertex_id);
  CHECK_ALLOC(*index);
  CHECK_ALLOC(*item);

  /* 2nd pass: read actual */
  fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s (2nd pass)\n", fname);
    exit(1);
  }

  {
    kint nv2;
    if(fscanf(fp, "%d", &nv2) != 1 || nv2 != *n_vertex){
      fprintf(stderr, "n_vertex mismatch in %s\n", fname);
      fclose(fp);
      exit(1);
    }
  }

  (*index)[0] = 0;
  kint pos = 0;
  for(kint i=0; i<*n_vertex; ++i){
    kint in;
    if(fscanf(fp, "%d %d", &((*vertex_id)[i]), &in) != 2){
      fprintf(stderr, "failed read vertex_id/deg at i=%d in %s\n", (int)i, fname);
      fclose(fp);
      exit(1);
    }
    for(kint j=0; j<in; ++j){
      if(fscanf(fp, "%d", &((*item)[pos + j])) != 1){
        fprintf(stderr, "failed read item at i=%d in %s\n", (int)i, fname);
        fclose(fp);
        exit(1);
      }
    }
    pos += in;
    (*index)[i+1] = pos;
  }

  fclose(fp);
}

void monolis_output_graph(const char* fname,
                          kint n_vertex,
                          const kint* vertex_id,
                          const kint* index,
                          const kint* item)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  fprintf(fp, "%d\n", (int)n_vertex);

  for(kint i=0; i<n_vertex; ++i){
    kint jS = index[i];
    kint jE = index[i+1];
    kint in = jE - jS;

    fprintf(fp, "%d %d", (int)vertex_id[i], (int)in);
    for(kint j=0; j<in; ++j){
      fprintf(fp, " %d", (int)item[jS + j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

/* ============================================================
 * n_internal フォーマット
 *
 * Fortran版:
 *  write "#n_internal 1"
 *  write n_internal
 *
 * inputは:
 *  1行目: label i
 *  2行目: n_internal
 * ============================================================ */
void monolis_input_internal_vertex_number(const char* fname, kint* n_internal)
{
  FILE* fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s\n", fname);
    exit(1);
  }

  char label[MONOLIS_CHARLEN];
  kint dummy_i;

  if(fscanf(fp, "%255s %d", label, &dummy_i) != 2){
    fprintf(stderr, "failed read header in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(fscanf(fp, "%d", n_internal) != 1){
    fprintf(stderr, "failed read n_internal in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  fclose(fp);
}

void monolis_output_internal_vertex_number(const char* fname, kint n_internal)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  fprintf(fp, "#n_internal 1\n");
  fprintf(fp, "%d\n", (int)n_internal);

  fclose(fp);
}

/* ============================================================
 * global id フォーマット
 *
 * Fortran版 output:
 *  "#id"
 *  n_vertex 1
 *  vertex_id (n_vertex 行)
 *
 * input:
 *  label
 *  n_vertex n_dof
 *  vertex_id...
 * ============================================================ */
void monolis_input_global_id(const char* fname,
                             kint* n_vertex,
                             kint** vertex_id)
{
  FILE* fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s\n", fname);
    exit(1);
  }

  char label[MONOLIS_CHARLEN];
  kint n_dof;

  if(fscanf(fp, "%255s", label) != 1){
    fprintf(stderr, "failed read label in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(fscanf(fp, "%d %d", n_vertex, &n_dof) != 2){
    fprintf(stderr, "failed read n_vertex n_dof in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(*n_vertex < 0){
    fprintf(stderr, "invalid n_vertex=%d in %s\n", *n_vertex, fname);
    fclose(fp);
    exit(1);
  }

  *vertex_id = (kint*)malloc(sizeof(kint) * (size_t)(*n_vertex));
  CHECK_ALLOC(*vertex_id);

  for(kint i=0; i<*n_vertex; ++i){
    if(fscanf(fp, "%d", &((*vertex_id)[i])) != 1){
      fprintf(stderr, "failed read vertex_id[%d] in %s\n", (int)i, fname);
      fclose(fp);
      exit(1);
    }
  }

  fclose(fp);
}

void monolis_output_global_id(const char* fname,
                              kint n_vertex,
                              const kint* vertex_id)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  fprintf(fp, "#id\n");
  fprintf(fp, "%d %d\n", (int)n_vertex, 1);

  for(kint i=0; i<n_vertex; ++i){
    fprintf(fp, "%d\n", (int)vertex_id[i]);
  }

  fclose(fp);
}

/* ============================================================
 * bc フォーマット (R)
 *
 * Fortran版:
 *  1行目: n_bc n_dof
 *  以降 n_bc 行: node_id dir value
 *
 * C側の ibc は 2*n_bc の連続配列:
 *   ibc[2*i+0]=node_id, ibc[2*i+1]=dir
 * rbc は n_bc
 * ============================================================ */
void monolis_input_bc_R(const char* fname,
                        kint* n_bc,
                        kint* n_dof,
                        kint** i_bc,
                        kdouble** r_bc)
{
  FILE* fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s\n", fname);
    exit(1);
  }

  if(fscanf(fp, "%d %d", n_bc, n_dof) != 2){
    fprintf(stderr, "failed read n_bc n_dof in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(*n_bc < 0 || *n_dof <= 0){
    fprintf(stderr, "invalid header in %s (n_bc=%d, n_dof=%d)\n",
            fname, (int)*n_bc, (int)*n_dof);
    fclose(fp);
    exit(1);
  }

  *i_bc = (kint*)malloc(sizeof(kint) * (size_t)(2 * (*n_bc)));
  *r_bc = (kdouble*)malloc(sizeof(kdouble) * (size_t)(*n_bc));
  CHECK_ALLOC(*i_bc);
  CHECK_ALLOC(*r_bc);

  for(kint i=0; i<*n_bc; ++i){
    kint node_id, dir;
    kdouble val;
    if(fscanf(fp, "%d %d %lf", &node_id, &dir, &val) != 3){
      fprintf(stderr, "failed read bc line %d in %s\n", (int)i, fname);
      fclose(fp);
      exit(1);
    }
    (*i_bc)[2*i + 0] = node_id;
    (*i_bc)[2*i + 1] = dir;
    (*r_bc)[i] = val;
  }

  fclose(fp);
}

void monolis_output_bc_R(const char* fname,
                         kint n_bc,
                         kint n_dof,
                         const kint* i_bc,
                         const kdouble* r_bc)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  fprintf(fp, "%d %d\n", (int)n_bc, (int)n_dof);

  for(kint i=0; i<n_bc; ++i){
    fprintf(fp, "%d %d %.14e\n",
            (int)i_bc[2*i + 0],
            (int)i_bc[2*i + 1],
            (double)r_bc[i]);
  }

  fclose(fp);
}

/* ============================================================
 * distval_R 出力
 *
 * Fortran版:
 *  1行目: label
 *  2行目: n_node n_dof
 *  以降 n_node 行、各行 n_dof 個を 1pe22.14 で出力
 *
 * val は node-major を想定:
 *   val[node*n_dof + dof]
 * ============================================================ */
void monolis_output_distval_R(const char* fname,
                              const char* label,
                              kint n_node,
                              kint n_dof,
                              const kdouble* val)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  fprintf(fp, "%s\n", label);
  fprintf(fp, "%d %d\n", (int)n_node, (int)n_dof);

  for(kint i=0; i<n_node; ++i){
    for(kint j=0; j<n_dof; ++j){
      /* Fortranの 1pe22.14 相当っぽく指数表記で揃える */
      fprintf(fp, "%.14e", (double)val[(size_t)i*(size_t)n_dof + (size_t)j]);
      if(j != n_dof-1) fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}


static int cmp_pair_node_id(const void* a, const void* b)
{
  const PairNodeIdx* pa = (const PairNodeIdx*)a;
  const PairNodeIdx* pb = (const PairNodeIdx*)b;
  if(pa->node_id < pb->node_id) return -1;
  if(pa->node_id > pb->node_id) return  1;
  return 0;
}

/* node_id(0-based) 配列（長さ n）から node_id 昇順の perm を作る
 * 返り値 perm[j] = ソート後 j 番目に来る「元のインデックス」
 */
static int* make_perm_by_sort_node_id_0based(const int* node_id_0based, int n)
{
  PairNodeIdx* pairs = (PairNodeIdx*)malloc(sizeof(PairNodeIdx)* (size_t)n);
  CHECK_ALLOC(pairs);

  for(int j=0; j<n; ++j){
    pairs[j].node_id   = node_id_0based[j];
    pairs[j].old_index = j;
  }

  qsort(pairs, (size_t)n, sizeof(PairNodeIdx), cmp_pair_node_id);

  int* perm = (int*)malloc(sizeof(int)* (size_t)n);
  CHECK_ALLOC(perm);

  for(int j=0; j<n; ++j){
    perm[j] = pairs[j].old_index;
  }

  free(pairs);
  return perm;
}

/* node.dat 読み込み（0-based 配列に格納）
 * ファイル形式:
 *  1行目: label
 *  2行目: n_node n_dof
 *  以降: n_node 行、各行 n_dof 個
 *
 * val は node-major: val[node*n_dof + dof]
 */
static int monolis_input_distval_R_c_0based(
    const char* fname,
    char        label[MONOLIS_CHARLEN],
    int*        n_node,
    int*        n_dof,
    double**    val
){
  FILE* fp = fopen(fname, "r");
  if(!fp) return -1;

  if(fscanf(fp, "%255s", label) != 1){
    fclose(fp);
    return -2;
  }

  if(fscanf(fp, "%d %d", n_node, n_dof) != 2){
    fclose(fp);
    return -3;
  }

  if(*n_node <= 0 || *n_dof <= 0){
    fclose(fp);
    return -4;
  }

  *val = (double*)malloc((size_t)(*n_node) * (size_t)(*n_dof) * sizeof(double));
  if(!(*val)){
    fclose(fp);
    return -5;
  }

  for(int i=0; i<*n_node; ++i){         /* 0..n_node-1 */
    for(int j=0; j<*n_dof; ++j){        /* 0..n_dof-1  */
      if(fscanf(fp, "%lf", &((*val)[(size_t)i*(size_t)(*n_dof) + (size_t)j])) != 1){
        free(*val);
        *val = NULL;
        fclose(fp);
        return -6;
      }
    }
  }

  fclose(fp);
  return 0;
}

void monolis_output_com_table_main_c(const char* fname,
                                           int n_neib,
                                           const int* neib_pe,
                                           const int* index,
                                           const int* item)
{
  FILE* fp = fopen(fname, "w");
  if(!fp){
    fprintf(stderr, "cannot open %s for write\n", fname);
    exit(1);
  }

  if(n_neib == 0){
    fprintf(fp, "0 0\n");
    fclose(fp);
    return;
  }

  int nz = index[n_neib]; /* index は 0-based で長さ n_neib+1 を想定 */
  fprintf(fp, "%d %d\n", n_neib, nz);

  for(int i=0; i<n_neib; ++i){
    fprintf(fp, "%d\n", neib_pe[i]);
  }
  for(int i=0; i<n_neib+1; ++i){
    fprintf(fp, "%d\n", index[i]);
  }
  for(int i=0; i<nz; ++i){
    fprintf(fp, "%d\n", item[i]);
  }

  fclose(fp);
}

/* Fortran互換の公開名 */
void monolis_output_send_com_table(const char* fname, const MONOLIS_COM* COM)
{
  monolis_output_com_table_main_c(fname,
                                 COM->send_n_neib,
                                 COM->send_neib_pe,
                                 COM->send_index,
                                 COM->send_item);
}

void monolis_output_recv_com_table(const char* fname, const MONOLIS_COM* COM)
{
  monolis_output_com_table_main_c(fname,
                                 COM->recv_n_neib,
                                 COM->recv_neib_pe,
                                 COM->recv_index,
                                 COM->recv_item);
}

/* ---------- input side（必要なら） ---------- */

void monolis_input_com_table_main_c(const char* fname,
                                          int* n_neib,
                                          int** neib_pe,
                                          int** index,
                                          int** item)
{
  FILE* fp = fopen(fname, "r");
  if(!fp){
    fprintf(stderr, "cannot open %s\n", fname);
    exit(1);
  }

  int nz = 0;
  if(fscanf(fp, "%d %d", n_neib, &nz) != 2){
    fprintf(stderr, "failed read header in %s\n", fname);
    fclose(fp);
    exit(1);
  }

  if(*n_neib < 0 || nz < 0){
    fprintf(stderr, "invalid header in %s (n_neib=%d nz=%d)\n", fname, *n_neib, nz);
    fclose(fp);
    exit(1);
  }

  if(*n_neib == 0){
    /* Fortran版は pointer を長さ1で確保して return。
       C版も安全のため長さ1を確保しておく（NULLでもいいが呼び出し側次第）。 */
    *neib_pe = (int*)malloc(sizeof(int));
    *index   = (int*)malloc(sizeof(int));
    *item    = (int*)malloc(sizeof(int));
    CHECK_ALLOC(*neib_pe);
    CHECK_ALLOC(*index);
    CHECK_ALLOC(*item);
    (*neib_pe)[0] = 0;
    (*index)[0]   = 0;
    (*item)[0]    = 0;
    fclose(fp);
    return;
  }

  *neib_pe = (int*)malloc(sizeof(int) * (size_t)(*n_neib));
  *index   = (int*)malloc(sizeof(int) * (size_t)(*n_neib + 1));
  *item    = (int*)malloc(sizeof(int) * (size_t)nz);
  CHECK_ALLOC(*neib_pe);
  CHECK_ALLOC(*index);
  CHECK_ALLOC(*item);

  for(int i=0; i<*n_neib; ++i){
    if(fscanf(fp, "%d", &((*neib_pe)[i])) != 1){
      fprintf(stderr, "failed read neib_pe[%d] in %s\n", i, fname);
      fclose(fp);
      exit(1);
    }
  }
  for(int i=0; i<*n_neib+1; ++i){
    if(fscanf(fp, "%d", &((*index)[i])) != 1){
      fprintf(stderr, "failed read index[%d] in %s\n", i, fname);
      fclose(fp);
      exit(1);
    }
  }

  /* 簡易整合チェック：index[n_neib] が nz か */
  if((*index)[*n_neib] != nz){
    /* ファイル先頭の nz と index末尾が食い違う場合でも
       読めるようにはするが警告は出しておく */
    fprintf(stderr, "warning: nz mismatch in %s (header=%d index_last=%d)\n",
            fname, nz, (*index)[*n_neib]);
    /* index末尾に合わせるなら:
       nz = (*index)[*n_neib];
       reallocate が必要なのでここではしない */
  }

  for(int i=0; i<nz; ++i){
    if(fscanf(fp, "%d", &((*item)[i])) != 1){
      fprintf(stderr, "failed read item[%d] in %s\n", i, fname);
      fclose(fp);
      exit(1);
    }
  }

  fclose(fp);
}

void monolis_input_send_com_table(const char* fname, MONOLIS_COM* COM)
{
  monolis_input_com_table_main_c(fname,
                                &COM->send_n_neib,
                                &COM->send_neib_pe,
                                &COM->send_index,
                                &COM->send_item);
}

void monolis_input_recv_com_table(const char* fname, MONOLIS_COM* COM)
{
  monolis_input_com_table_main_c(fname,
                                &COM->recv_n_neib,
                                &COM->recv_neib_pe,
                                &COM->recv_index,
                                &COM->recv_item);
}

/* ------------------------------------------------------------
 * main
 * ------------------------------------------------------------ */
int main(int argc, char** argv)
{
  monolis_mpi_initialize();

  if(argc < 2){
    if(monolis_mpi_get_global_my_rank() == 0){
      fprintf(stderr, "Error : argv is need.\n");
    }
    monolis_mpi_finalize();
    return 1;
  }

  int ndomain_basis = atoi(argv[1]);
  int nprocs = monolis_mpi_get_global_comm_size();
  int myrank = monolis_mpi_get_global_my_rank();

  if(ndomain_basis < nprocs){
    if(myrank == 0) fprintf(stderr, "Error : ndomain_basis is invalid.\n");
    monolis_mpi_finalize();
    return 1;
  }

  if(ndomain_basis == 1 && nprocs == 1){
    monolis_mpi_finalize();
    return 0;
  }

  /* metagraph から n_graphs を読む */
  int n_graphs = 0;
  if(nprocs == 1){
    const char* fname = "./parted.0/parted.1/metagraph.dat";
    FILE* fp = fopen(fname, "r");
    if(!fp){
      fprintf(stderr, "cannot open %s\n", fname);
      monolis_mpi_finalize();
      return 1;
    }
    if(fscanf(fp, "%d", &n_graphs) != 1){
      fprintf(stderr, "failed read n_graphs from %s\n", fname);
      fclose(fp);
      monolis_mpi_finalize();
      return 1;
    }
    fclose(fp);
  } else {
    char fname[MONOLIS_CHARLEN];
    strcpy(fname, monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR,
                                                     MONOLIS_DEFAULT_PART_DIR,
                                                     "metagraph.dat.n_internal"));
    monolis_input_internal_vertex_number(fname, &n_graphs);
  }

  printf("my_rank %d n_graphs %d\n", myrank, n_graphs);

  /* merged structures */
  GEDATSU_GRAPH merged_nodal_graph;
  GEDATSU_GRAPH merged_conn_graph;
  MONOLIS_COM   merged_monoCOM;

  gedatsu_graph_initialize(&merged_nodal_graph);
  gedatsu_graph_initialize(&merged_conn_graph);
  monolis_com_initialize_by_self(&merged_monoCOM);

  /* arrays */
  GEDATSU_GRAPH* nodal_graphs = (GEDATSU_GRAPH*)malloc(sizeof(GEDATSU_GRAPH) * (size_t)n_graphs);
  GEDATSU_GRAPH* conn_graphs  = (GEDATSU_GRAPH*)malloc(sizeof(GEDATSU_GRAPH) * (size_t)n_graphs);
  MONOLIS_COM*   monoCOMs     = (MONOLIS_COM*)malloc(sizeof(MONOLIS_COM)   * (size_t)n_graphs);
  CHECK_ALLOC(nodal_graphs);
  CHECK_ALLOC(conn_graphs);
  CHECK_ALLOC(monoCOMs);

  /* list structs（id は 0..n_graphs-1） */
  MONOLIS_LIST_I n_dof_list_node, n_dof_list_bc, n_dof_list_load;
  MONOLIS_LIST_I list_ibc, list_iload;
  MONOLIS_LIST_R list_node, list_rbc, list_rload;

  monolis_list_initialize_I(&n_dof_list_node, n_graphs);
  monolis_list_initialize_I(&n_dof_list_bc,   n_graphs);
  monolis_list_initialize_I(&n_dof_list_load, n_graphs);
  monolis_list_initialize_I(&list_ibc,        n_graphs);
  monolis_list_initialize_I(&list_iload,      n_graphs);
  monolis_list_initialize_R(&list_node,       n_graphs);
  monolis_list_initialize_R(&list_rbc,        n_graphs);
  monolis_list_initialize_R(&list_rload,      n_graphs);

  /* domain_id (この rank が読む subgraph のID) */
  int* domain_id = (int*)malloc(sizeof(int) * (size_t)n_graphs);
  CHECK_ALLOC(domain_id);

  if(nprocs == 1){
    for(int i=0;i<n_graphs;i++){
      domain_id[i] = i; /* 0..n_graphs-1 */
    }
  } else {
    char cnum[64];
    snprintf(cnum, sizeof(cnum), "%d", myrank);
    char fname[MONOLIS_CHARLEN];
    snprintf(fname, sizeof(fname), "./parted.0/metagraph.dat.id.%s", cnum);

    FILE* fp = fopen(fname, "r");
    if(!fp){
      fprintf(stderr, "cannot open %s\n", fname);
      monolis_mpi_finalize();
      return 1;
    }

    /* 2行読み飛ばし */
    {
      char buf[1024];
      fgets(buf, sizeof(buf), fp);
      fgets(buf, sizeof(buf), fp);
    }

    for(int i=0;i<n_graphs;i++){
      if(fscanf(fp, "%d", &domain_id[i]) != 1){
        fprintf(stderr, "failed read domain_id from %s\n", fname);
        fclose(fp);
        monolis_mpi_finalize();
        return 1;
      }
      domain_id[i] -= 1; /* Fortran同様 -1（ここは 0-based の domain id にしたい） */
    }
    fclose(fp);
  }

  /* ------------------------------------------------------------
   * read graphs & distvals, convert to list (0-based)
   * ------------------------------------------------------------ */
  for(int gi=0; gi<n_graphs; ++gi){
    gedatsu_graph_initialize(&nodal_graphs[gi]);
    gedatsu_graph_initialize(&conn_graphs[gi]);
    monolis_com_initialize_by_self(&monoCOMs[gi]);

    char cnum[64];
    snprintf(cnum, sizeof(cnum), "%d", domain_id[gi]);

    /* ---------- nodal graph ---------- */
    {
      char fname[MONOLIS_CHARLEN];

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/graph.dat.%s", cnum);
      monolis_input_graph(fname,
        &nodal_graphs[gi].n_vertex,
        &nodal_graphs[gi].vertex_id,
        &nodal_graphs[gi].index,
        &nodal_graphs[gi].item);

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/graph.dat.n_internal.%s", cnum);
      monolis_input_internal_vertex_number(fname, &nodal_graphs[gi].n_internal_vertex);

      free(nodal_graphs[gi].vertex_id);
      nodal_graphs[gi].vertex_id = NULL;

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/graph.dat.id.%s", cnum);
      monolis_input_global_id(fname, nodal_graphs[gi].n_vertex, &nodal_graphs[gi].vertex_id);

      nodal_graphs[gi].vertex_domain_id = (int*)malloc(sizeof(int) * (size_t)nodal_graphs[gi].n_vertex);
      CHECK_ALLOC(nodal_graphs[gi].vertex_domain_id);
      for(int j=0; j<nodal_graphs[gi].n_vertex; ++j){
        nodal_graphs[gi].vertex_domain_id[j] = myrank;
      }

      monolis_com_input_comm_table(&monoCOMs[gi], "./", "parted.1", "graph.dat");
    }

    /* ---------- connectivity graph ---------- */
    {
      char fname[MONOLIS_CHARLEN];

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/connectivity.dat.%s", cnum);
      monolis_input_graph(fname,
        &conn_graphs[gi].n_vertex,
        &conn_graphs[gi].vertex_id,
        &conn_graphs[gi].index,
        &conn_graphs[gi].item);

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/connectivity.dat.n_internal.%s", cnum);
      monolis_input_internal_vertex_number(fname, &conn_graphs[gi].n_internal_vertex);

      free(conn_graphs[gi].vertex_id);
      conn_graphs[gi].vertex_id = NULL;

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/connectivity.dat.id.%s", cnum);
      monolis_input_global_id(fname, conn_graphs[gi].n_vertex, &conn_graphs[gi].vertex_id);

      conn_graphs[gi].vertex_domain_id = (int*)malloc(sizeof(int) * (size_t)conn_graphs[gi].n_vertex);
      CHECK_ALLOC(conn_graphs[gi].vertex_domain_id);
      for(int j=0; j<conn_graphs[gi].n_vertex; ++j){
        conn_graphs[gi].vertex_domain_id[j] = myrank;
      }
    }

    /* ---------- distval node ---------- */
    char label[MONOLIS_CHARLEN];
    int  nnode = 0, ndof = 0;
    double* node_val = NULL; /* node-major: node*n_dof + dof */

    {
      char fname[MONOLIS_CHARLEN];
      snprintf(fname, sizeof(fname), "./parted.0/parted.1/node.dat.%s", cnum);

      int ierr = monolis_input_distval_R_c_0based(fname, label, &nnode, &ndof, &node_val);
      if(ierr != 0){
        fprintf(stderr, "failed monolis_input_distval_R_c_0based %s (ierr=%d)\n", fname, ierr);
        monolis_mpi_finalize();
        return 1;
      }
      if(ndof != 3){
        fprintf(stderr, "ndof isn't 3\n");
        monolis_mpi_finalize();
        return 1;
      }

      /* n_dof_list_node: 各節点の dof 数（全部 3） */
      int* dof_cnt = (int*)malloc(sizeof(int) * (size_t)nnode);
      CHECK_ALLOC(dof_cnt);
      for(int j=0; j<nnode; ++j) dof_cnt[j] = ndof;
      monolis_list_set_I(&n_dof_list_node, gi, nnode, dof_cnt);
      free(dof_cnt);

      /* list_node: node_val をそのまま 1D 配列として登録（node-major） */
      monolis_list_set_R(&list_node, gi, nnode*ndof, node_val);

      free(node_val);
      node_val = NULL;
    }

    /* ---------- bc / load ---------- */
    int nbc=0, nload=0, ndof_bc=0, ndof_load=0;
    int* ibc = NULL;     /* 想定: 2*nbc 連続配列 [node_id, dir, node_id, dir, ...] */
    int* iload = NULL;   /* 同上 */
    double* rbc = NULL;  /* nbc */
    double* rload = NULL;/* nload */

    {
      char fname[MONOLIS_CHARLEN];

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/bc.dat.%s", cnum);
      monolis_input_bc_R(fname, &nbc, &ndof_bc, &ibc, &rbc);
      if(ndof_bc != 3){
        fprintf(stderr, "ndof isn't 3\n");
        monolis_mpi_finalize();
        return 1;
      }

      snprintf(fname, sizeof(fname), "./parted.0/parted.1/load.dat.%s", cnum);
      monolis_input_bc_R(fname, &nload, &ndof_load, &iload, &rload);
      if(ndof_load != 3){
        fprintf(stderr, "ndof isn't 3\n");
        monolis_mpi_finalize();
        return 1;
      }

      /* -----------------------
       * 0スタート化：node_id を -1 する
       * ----------------------- */
      for(int j=0; j<nbc; ++j){
        ibc[2*j + 0] -= 1; /* 1..nnode -> 0..nnode-1 */
      }
      for(int j=0; j<nload; ++j){
        iload[2*j + 0] -= 1;
      }

      /* --- bc を list に変換（node_id 昇順で並べる） --- */

      /* n_dof_list_bc: 各節点に付く Dirichlet 条件数（0-based node_id でカウント） */
      int* cnt = (int*)calloc((size_t)nnode, sizeof(int));
      CHECK_ALLOC(cnt);
      for(int j=0; j<nbc; ++j){
        int node0 = ibc[2*j + 0];
        if(node0 < 0 || node0 >= nnode){
          fprintf(stderr, "bc node_id out of range (node0=%d, nnode=%d)\n", node0, nnode);
          monolis_mpi_finalize();
          return 1;
        }
        cnt[node0] += 1;
      }
      monolis_list_set_I(&n_dof_list_bc, gi, nnode, cnt);
      free(cnt);

      /* perm: ibc の node_id(0-based) で昇順 */
      int* node_ids = (int*)malloc(sizeof(int) * (size_t)nbc);
      CHECK_ALLOC(node_ids);
      for(int j=0; j<nbc; ++j) node_ids[j] = ibc[2*j + 0];
      int* perm = make_perm_by_sort_node_id_0based(node_ids, nbc);
      free(node_ids);

      /* list_ibc: 拘束方向（ソート順） */
      int* bc_dir = (int*)malloc(sizeof(int) * (size_t)nbc);
      CHECK_ALLOC(bc_dir);
      for(int j=0; j<nbc; ++j){
        int old = perm[j];
        bc_dir[j] = ibc[2*old + 1];
      }
      monolis_list_set_I(&list_ibc, gi, nbc, bc_dir);
      free(bc_dir);

      /* list_rbc: 値（ソート順） */
      double* bc_val = (double*)malloc(sizeof(double) * (size_t)nbc);
      CHECK_ALLOC(bc_val);
      for(int j=0; j<nbc; ++j){
        int old = perm[j];
        bc_val[j] = rbc[old];
      }
      monolis_list_set_R(&list_rbc, gi, nbc, bc_val);
      free(bc_val);

      free(perm);

      /* --- load を list に変換（node_id 昇順で並べる） --- */

      cnt = (int*)calloc((size_t)nnode, sizeof(int));
      CHECK_ALLOC(cnt);
      for(int j=0; j<nload; ++j){
        int node0 = iload[2*j + 0];
        if(node0 < 0 || node0 >= nnode){
          fprintf(stderr, "load node_id out of range (node0=%d, nnode=%d)\n", node0, nnode);
          monolis_mpi_finalize();
          return 1;
        }
        cnt[node0] += 1;
      }
      monolis_list_set_I(&n_dof_list_load, gi, nnode, cnt);
      free(cnt);

      node_ids = (int*)malloc(sizeof(int) * (size_t)nload);
      CHECK_ALLOC(node_ids);
      for(int j=0; j<nload; ++j) node_ids[j] = iload[2*j + 0];
      perm = make_perm_by_sort_node_id_0based(node_ids, nload);
      free(node_ids);

      int* load_dir = (int*)malloc(sizeof(int) * (size_t)nload);
      CHECK_ALLOC(load_dir);
      for(int j=0; j<nload; ++j){
        int old = perm[j];
        load_dir[j] = iload[2*old + 1];
      }
      monolis_list_set_I(&list_iload, gi, nload, load_dir);
      free(load_dir);

      double* load_val = (double*)malloc(sizeof(double) * (size_t)nload);
      CHECK_ALLOC(load_val);
      for(int j=0; j<nload; ++j){
        int old = perm[j];
        load_val[j] = rload[old];
      }
      monolis_list_set_R(&list_rload, gi, nload, load_val);
      free(load_val);

      free(perm);

      free(ibc);
      free(iload);
      free(rbc);
      free(rload);
    }
  }

  /* ------------------------------------------------------------
   * merge
   * ------------------------------------------------------------ */
  gedatsu_merge_nodal_subgraphs(n_graphs, nodal_graphs, monoCOMs,
                               &merged_nodal_graph, &merged_monoCOM,
                               ORDER_NODAL_ID);

  gedatsu_merge_connectivity_subgraphs(n_graphs, nodal_graphs,
                                       &merged_nodal_graph, &merged_monoCOM,
                                       n_graphs, conn_graphs,
                                       &merged_conn_graph);

  int*    merged_n_dof_list_node = NULL;
  double* merged_array_node      = NULL;

  int*    merged_n_dof_list_bc   = NULL;
  int*    merged_array_ibc       = NULL;
  double* merged_array_rbc       = NULL;

  int*    merged_n_dof_list_load = NULL;
  int*    merged_array_iload     = NULL;
  double* merged_array_rload     = NULL;

  gedatsu_merge_distval_R(n_graphs, nodal_graphs, &merged_nodal_graph,
                          &n_dof_list_node, &list_node,
                          &merged_n_dof_list_node, &merged_array_node);

  gedatsu_merge_distval_I(n_graphs, nodal_graphs, &merged_nodal_graph,
                          &n_dof_list_bc, &list_ibc,
                          &merged_n_dof_list_bc, &merged_array_ibc);

  gedatsu_merge_distval_R(n_graphs, nodal_graphs, &merged_nodal_graph,
                          &n_dof_list_bc, &list_rbc,
                          &merged_n_dof_list_bc, &merged_array_rbc);

  gedatsu_merge_distval_I(n_graphs, nodal_graphs, &merged_nodal_graph,
                          &n_dof_list_load, &list_iload,
                          &merged_n_dof_list_load, &merged_array_iload);

  gedatsu_merge_distval_R(n_graphs, nodal_graphs, &merged_nodal_graph,
                          &n_dof_list_load, &list_rload,
                          &merged_n_dof_list_load, &merged_array_rload);

  /* ------------------------------------------------------------
   * list -> 出力用配列へ復元（0-based node id で bc/load を作る）
   * ------------------------------------------------------------ */
  int nnode_merged = merged_nodal_graph.n_vertex;
  int ndof = 3;

  /* node_out は node-major (node*n_dof + dof) */
  double* node_out = (double*)malloc(sizeof(double) * (size_t)nnode_merged * (size_t)ndof);
  CHECK_ALLOC(node_out);
  for(int i=0; i<nnode_merged; ++i){
    for(int j=0; j<ndof; ++j){
      node_out[(size_t)i*(size_t)ndof + (size_t)j] =
        merged_array_node[(size_t)i*(size_t)ndof + (size_t)j];
    }
  }

  /* bc の要素数 = sum(merged_n_dof_list_bc) */
  int nbc = 0;
  for(int i=0; i<nnode_merged; ++i) nbc += merged_n_dof_list_bc[i];

  int* ibc_out = (int*)malloc(sizeof(int) * (size_t)2 * (size_t)nbc);
  double* rbc_out = (double*)malloc(sizeof(double) * (size_t)nbc);
  CHECK_ALLOC(ibc_out);
  CHECK_ALLOC(rbc_out);

  int pos = 0;
  for(int node0=0; node0<nnode_merged; ++node0){
    int cnt = merged_n_dof_list_bc[node0];
    for(int t=0; t<cnt; ++t){
      /* 0-based node id で復元 */
      ibc_out[2*pos + 0] = node0;
      ibc_out[2*pos + 1] = merged_array_ibc[pos];
      rbc_out[pos]       = merged_array_rbc[pos];
      pos++;
    }
  }

  /* load の要素数 = sum(merged_n_dof_list_load) */
  int nload = 0;
  for(int i=0; i<nnode_merged; ++i) nload += merged_n_dof_list_load[i];

  int* iload_out = (int*)malloc(sizeof(int) * (size_t)2 * (size_t)nload);
  double* rload_out = (double*)malloc(sizeof(double) * (size_t)nload);
  CHECK_ALLOC(iload_out);
  CHECK_ALLOC(rload_out);

  pos = 0;
  for(int node0=0; node0<nnode_merged; ++node0){
    int cnt = merged_n_dof_list_load[node0];
    for(int t=0; t<cnt; ++t){
      iload_out[2*pos + 0] = node0;                 /* 0-based */
      iload_out[2*pos + 1] = merged_array_iload[pos];
      rload_out[pos]       = merged_array_rload[pos];
      pos++;
    }
  }

  /* ------------------------------------------------------------
   * output
   * ------------------------------------------------------------ */
  {
    char fname[MONOLIS_CHARLEN];

    /* graph */
    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "graph.dat"));
    monolis_output_graph(fname,
      merged_nodal_graph.n_vertex,
      merged_nodal_graph.vertex_id,
      merged_nodal_graph.index,
      merged_nodal_graph.item);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "graph.dat.n_internal"));
    monolis_output_internal_vertex_number(fname, merged_nodal_graph.n_internal_vertex);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "graph.dat.id"));
    monolis_output_global_id(fname, merged_nodal_graph.n_vertex, merged_nodal_graph.vertex_id);

    /* com */
    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "graph.dat.send"));
    monolis_output_send_com_table(fname, &merged_monoCOM);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "graph.dat.recv"));
    monolis_output_recv_com_table(fname, &merged_monoCOM);

    /* connectivity */
    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "connectivity.dat"));
    monolis_output_graph(fname,
      merged_conn_graph.n_vertex,
      merged_conn_graph.vertex_id,
      merged_conn_graph.index,
      merged_conn_graph.item);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "connectivity.dat.n_internal"));
    monolis_output_internal_vertex_number(fname, merged_conn_graph.n_internal_vertex);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "connectivity.dat.id"));
    monolis_output_global_id(fname, merged_conn_graph.n_vertex, merged_conn_graph.vertex_id);

    /* distval node */
    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "node.dat"));
    {
      char out_label[MONOLIS_CHARLEN];
      strcpy(out_label, "#node");
      /* ここも C API の期待する val の並びに合わせてください。
         今回は node-major の 1D を渡す想定。 */
      monolis_output_distval_R(fname, out_label, nnode_merged, ndof, node_out);
    }

    /* bc / load：節点番号は 0-based で出力する */
    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "bc.dat"));
    monolis_output_bc_R(fname, nbc, ndof, ibc_out, rbc_out);

    strcpy(fname, monolis_get_global_output_file_name("./parted.0/", "", "load.dat"));
    monolis_output_bc_R(fname, nload, ndof, iload_out, rload_out);
  }

  /* ------------------------------------------------------------
   * cleanup
   * ------------------------------------------------------------ */
  free(domain_id);

  free(node_out);
  free(ibc_out);
  free(rbc_out);
  free(iload_out);
  free(rload_out);

  free(merged_n_dof_list_node);
  free(merged_array_node);
  free(merged_n_dof_list_bc);
  free(merged_array_ibc);
  free(merged_array_rbc);
  free(merged_n_dof_list_load);
  free(merged_array_iload);
  free(merged_array_rload);

  monolis_list_finalize_I(&n_dof_list_node, n_graphs);
  monolis_list_finalize_I(&n_dof_list_bc,   n_graphs);
  monolis_list_finalize_I(&n_dof_list_load, n_graphs);
  monolis_list_finalize_I(&list_ibc,        n_graphs);
  monolis_list_finalize_I(&list_iload,      n_graphs);
  monolis_list_finalize_R(&list_node,       n_graphs);
  monolis_list_finalize_R(&list_rbc,        n_graphs);
  monolis_list_finalize_R(&list_rload,      n_graphs);

  free(nodal_graphs);
  free(conn_graphs);
  free(monoCOMs);

  monolis_mpi_finalize();
  return 0;
}
