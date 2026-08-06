#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static const char* CODENAME = "util/vtkddm/vtkddm_tet4 >";
static const char* VOIDNAME = "                         ";

#define BUFFER_SIZE 10000
#define TET4_NUM_NODE 4
#define VTK_TETRA 10
#define COORD_SCALE 1000000000000.0

static const char* DEF_DIRECTORY          = ".";
static const char* DEF_PART_DIR           = "parted.0";
static const char* DEF_PART_NODE_BODY     = "sorted_nodes.dat";
static const char* DEF_PART_ELEM_BODY     = "sorted_local_elem.dat";
static const char* DEF_OUTFILE_VTK        = "vtkddm_tet4.vtk";
static const char* VTK_CELL_DATA_LABEL    = "num_parted";

static const char* OPT_DIRECTORY          = "-d";
static const char* OPT_PART_DIR           = "-id";
static const char* OPT_PART_NODE_BODY     = "-np";
static const char* OPT_PART_ELEM_BODY     = "-ip";
static const char* OPT_OUTFILE_VTK        = "-ov";
static const char* OPT_INDEX_BASE         = "-base";

typedef enum
{
    INDEX_BASE_AUTO = -1,
    INDEX_BASE_0    = 0,
    INDEX_BASE_1    = 1
} INDEX_BASE_MODE;

typedef struct
{
    const char* directory;
    const char* part_dir;
    const char* part_node_body;
    const char* part_elem_body;
    const char* outfile_vtk;
    INDEX_BASE_MODE index_base_mode;
} SETTINGS;

typedef struct
{
    int count;
    int capacity;

    double* x;
    double* y;
    double* z;

    long long* qx;
    long long* qy;
    long long* qz;
} NODE_ARRAY;

typedef struct
{
    int size;
    int count;
    int* ids;
} NODE_HASH;

typedef struct
{
    int conn[TET4_NUM_NODE];
    int key[TET4_NUM_NODE];
    int rank;
} ELEM_RECORD;

typedef struct
{
    int count;
    int capacity;
    ELEM_RECORD* elem;
} ELEM_ARRAY;

typedef struct
{
    int count;
    int capacity;
    int* conn;
    double* cell_value;
} VTK_CELL_ARRAY;

/* -------------------------------------------------------------------------- */
/* Small utilities                                                             */
/* -------------------------------------------------------------------------- */

static void* checked_malloc(size_t nbytes)
{
    void* p = malloc(nbytes);

    if(p == NULL) {
        fprintf(stderr, "%s Failed to allocate %zu bytes.\n", CODENAME, nbytes);
        exit(EXIT_FAILURE);
    }

    return p;
}

static void* checked_realloc(void* ptr, size_t nbytes)
{
    void* p = realloc(ptr, nbytes);

    if(p == NULL) {
        fprintf(stderr, "%s Failed to reallocate %zu bytes.\n", CODENAME, nbytes);
        exit(EXIT_FAILURE);
    }

    return p;
}

static int find_option(int argc, char* argv[], const char* option)
{
    for(int i = 0; i < argc; i++) {
        if(strcmp(argv[i], option) == 0) {
            return i;
        }
    }

    return -1;
}

static const char* get_option_string(
    int         argc,
    char*       argv[],
    const char* option,
    const char* default_value)
{
    int pos = find_option(argc, argv, option);

    if(pos < 0) {
        return default_value;
    }

    if(pos + 1 >= argc) {
        fprintf(stderr, "%s Option %s requires a value.\n", CODENAME, option);
        exit(EXIT_FAILURE);
    }

    return argv[pos + 1];
}

static INDEX_BASE_MODE get_index_base_mode(int argc, char* argv[])
{
    const char* value = get_option_string(argc, argv, OPT_INDEX_BASE, "auto");

    if(strcmp(value, "auto") == 0) {
        return INDEX_BASE_AUTO;
    }

    if(strcmp(value, "0") == 0) {
        return INDEX_BASE_0;
    }

    if(strcmp(value, "1") == 0) {
        return INDEX_BASE_1;
    }

    fprintf(stderr,
            "%s Invalid %s value: %s. Use auto, 0, or 1.\n",
            CODENAME,
            OPT_INDEX_BASE,
            value);
    exit(EXIT_FAILURE);
}

static void print_usage(const char* argv0)
{
    printf("%s Please input parameters.\n", CODENAME);
    printf("%s Format:\n", CODENAME);
    printf("%s    %s [the num. of parallel]\n", VOIDNAME, argv0);
    printf("%s Options:\n", CODENAME);
    printf("%s    %s [input/output directory]\n", VOIDNAME, OPT_DIRECTORY);
    printf("%s    %s [partition directory name]\n", VOIDNAME, OPT_PART_DIR);
    printf("%s    %s [partition node filename body]\n", VOIDNAME, OPT_PART_NODE_BODY);
    printf("%s       Example: -np sorted_local_node.dat reads sorted_local_node.dat.0, ...\n", VOIDNAME);
    printf("%s    %s [partition element filename body]\n", VOIDNAME, OPT_PART_ELEM_BODY);
    printf("%s       Example: -ip sorted_local_elem.dat reads sorted_local_elem.dat.0, ...\n", VOIDNAME);
    printf("%s    %s [node index base: auto, 0, or 1]\n", VOIDNAME, OPT_INDEX_BASE);
    printf("%s    %s [output filename (VTK)]\n", VOIDNAME, OPT_OUTFILE_VTK);
    printf("\n");
}

static int args_manager(SETTINGS* set, int argc, char* argv[])
{
    if(argc < 2) {
        print_usage(argv[0]);
        exit(EXIT_SUCCESS);
    }

    int num_parallel = atoi(argv[1]);

    set->directory       = get_option_string(argc, argv, OPT_DIRECTORY,      DEF_DIRECTORY);
    set->part_dir        = get_option_string(argc, argv, OPT_PART_DIR,       DEF_PART_DIR);
    set->part_node_body  = get_option_string(argc, argv, OPT_PART_NODE_BODY, DEF_PART_NODE_BODY);
    set->part_elem_body  = get_option_string(argc, argv, OPT_PART_ELEM_BODY, DEF_PART_ELEM_BODY);
    set->outfile_vtk     = get_option_string(argc, argv, OPT_OUTFILE_VTK,    DEF_OUTFILE_VTK);
    set->index_base_mode = get_index_base_mode(argc, argv);

    printf("%s The num. of parallel                : %d\n", CODENAME, num_parallel);
    printf("%s Input/output directory              : %s\n", CODENAME, set->directory);
    printf("%s Input directory (partition files)   : %s\n", CODENAME, set->part_dir);
    printf("%s Input filename body (part nodes)    : %s\n", CODENAME, set->part_node_body);
    printf("%s Input filename body (part elements) : %s\n", CODENAME, set->part_elem_body);
    printf("%s Element node index base             : ", CODENAME);

    if(set->index_base_mode == INDEX_BASE_AUTO) {
        printf("auto\n");
    }
    else {
        printf("%d\n", (int)set->index_base_mode);
    }

    printf("%s Output filename (VTK)               : %s\n", CODENAME, set->outfile_vtk);
    printf("%s Global node.dat                     : not used\n", CODENAME);
    printf("%s Global elem.dat                     : not used\n", CODENAME);

    return num_parallel;
}

static void make_rank_filename(
    char*           output,
    int             output_size,
    const char*     body,
    const SETTINGS* set,
    int             rank)
{
    snprintf(output,
             output_size,
             "%s/%s.%d",
             set->part_dir,
             body,
             rank);
}

static FILE* open_read_file(const SETTINGS* set, const char* filename)
{
    char path[BUFFER_SIZE];

    snprintf(path, BUFFER_SIZE, "%s/%s", set->directory, filename);

    printf("%s Reading file \"%s\".\n", CODENAME, path);

    FILE* fp = fopen(path, "r");

    if(fp == NULL) {
        fprintf(stderr, "%s Failed to open file \"%s\".\n", CODENAME, path);
        exit(EXIT_FAILURE);
    }

    return fp;
}

static FILE* open_write_file(const SETTINGS* set, const char* filename)
{
    char path[BUFFER_SIZE];

    snprintf(path, BUFFER_SIZE, "%s/%s", set->directory, filename);

    printf("%s Writing file \"%s\".\n", CODENAME, path);

    FILE* fp = fopen(path, "w");

    if(fp == NULL) {
        fprintf(stderr, "%s Failed to open file \"%s\" for writing.\n", CODENAME, path);
        exit(EXIT_FAILURE);
    }

    return fp;
}

/* -------------------------------------------------------------------------- */
/* Global node merging by coordinate                                           */
/* -------------------------------------------------------------------------- */

static long long quantize_coord(double value)
{
    double scaled = value * COORD_SCALE;

    if(scaled >= 0.0) {
        return (long long)(scaled + 0.5);
    }

    return (long long)(scaled - 0.5);
}

static unsigned int hash_u64(uint64_t x)
{
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;

    return (unsigned int)x;
}

static unsigned int hash_node_key(long long qx, long long qy, long long qz)
{
    uint64_t h = (uint64_t)qx;

    h ^= ((uint64_t)qy + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2));
    h ^= ((uint64_t)qz + UINT64_C(0x9e3779b97f4a7c15) + (h << 6) + (h >> 2));

    return hash_u64(h);
}

static int next_power_of_two(int value)
{
    int n = 1;

    while(n < value) {
        n <<= 1;
    }

    return n;
}

static void node_array_init(NODE_ARRAY* nodes)
{
    nodes->count = 0;
    nodes->capacity = 1024;

    nodes->x  = (double*)checked_malloc(sizeof(double)    * (size_t)nodes->capacity);
    nodes->y  = (double*)checked_malloc(sizeof(double)    * (size_t)nodes->capacity);
    nodes->z  = (double*)checked_malloc(sizeof(double)    * (size_t)nodes->capacity);
    nodes->qx = (long long*)checked_malloc(sizeof(long long) * (size_t)nodes->capacity);
    nodes->qy = (long long*)checked_malloc(sizeof(long long) * (size_t)nodes->capacity);
    nodes->qz = (long long*)checked_malloc(sizeof(long long) * (size_t)nodes->capacity);
}

static void node_array_reserve(NODE_ARRAY* nodes, int min_capacity)
{
    if(nodes->capacity >= min_capacity) {
        return;
    }

    while(nodes->capacity < min_capacity) {
        nodes->capacity *= 2;
    }

    nodes->x  = (double*)checked_realloc(nodes->x,  sizeof(double)    * (size_t)nodes->capacity);
    nodes->y  = (double*)checked_realloc(nodes->y,  sizeof(double)    * (size_t)nodes->capacity);
    nodes->z  = (double*)checked_realloc(nodes->z,  sizeof(double)    * (size_t)nodes->capacity);
    nodes->qx = (long long*)checked_realloc(nodes->qx, sizeof(long long) * (size_t)nodes->capacity);
    nodes->qy = (long long*)checked_realloc(nodes->qy, sizeof(long long) * (size_t)nodes->capacity);
    nodes->qz = (long long*)checked_realloc(nodes->qz, sizeof(long long) * (size_t)nodes->capacity);
}

static int node_array_append(
    NODE_ARRAY* nodes,
    double      x,
    double      y,
    double      z,
    long long   qx,
    long long   qy,
    long long   qz)
{
    node_array_reserve(nodes, nodes->count + 1);

    int id = nodes->count;

    nodes->x[id]  = x;
    nodes->y[id]  = y;
    nodes->z[id]  = z;
    nodes->qx[id] = qx;
    nodes->qy[id] = qy;
    nodes->qz[id] = qz;

    nodes->count++;

    return id;
}

static void node_hash_init(NODE_HASH* hash, int initial_size)
{
    hash->size = next_power_of_two(initial_size);

    if(hash->size < 2048) {
        hash->size = 2048;
    }

    hash->count = 0;
    hash->ids = (int*)checked_malloc(sizeof(int) * (size_t)hash->size);

    for(int i = 0; i < hash->size; i++) {
        hash->ids[i] = -1;
    }
}

static void node_hash_free(NODE_HASH* hash)
{
    free(hash->ids);

    hash->ids = NULL;
    hash->size = 0;
    hash->count = 0;
}

static int node_key_equals(
    const NODE_ARRAY* nodes,
    int               id,
    long long         qx,
    long long         qy,
    long long         qz)
{
    return nodes->qx[id] == qx &&
           nodes->qy[id] == qy &&
           nodes->qz[id] == qz;
}

static void node_hash_insert_existing(
    NODE_HASH*        hash,
    const NODE_ARRAY* nodes,
    int               node_id)
{
    unsigned int h = hash_node_key(
        nodes->qx[node_id],
        nodes->qy[node_id],
        nodes->qz[node_id]);

    int mask = hash->size - 1;
    int pos = (int)(h & (unsigned int)mask);

    while(hash->ids[pos] >= 0) {
        pos = (pos + 1) & mask;
    }

    hash->ids[pos] = node_id;
    hash->count++;
}

static void node_hash_rebuild(NODE_HASH* hash, const NODE_ARRAY* nodes)
{
    NODE_HASH new_hash;

    node_hash_init(&new_hash, hash->size * 2);

    for(int i = 0; i < nodes->count; i++) {
        node_hash_insert_existing(&new_hash, nodes, i);
    }

    node_hash_free(hash);

    *hash = new_hash;
}

static int find_or_add_global_node(
    NODE_ARRAY* nodes,
    NODE_HASH*  hash,
    double      x,
    double      y,
    double      z)
{
    long long qx = quantize_coord(x);
    long long qy = quantize_coord(y);
    long long qz = quantize_coord(z);

    if((hash->count + 1) * 10 > hash->size * 7) {
        node_hash_rebuild(hash, nodes);
    }

    unsigned int h = hash_node_key(qx, qy, qz);
    int mask = hash->size - 1;
    int pos = (int)(h & (unsigned int)mask);

    while(hash->ids[pos] >= 0) {
        int id = hash->ids[pos];

        if(node_key_equals(nodes, id, qx, qy, qz)) {
            return id;
        }

        pos = (pos + 1) & mask;
    }

    int new_id = node_array_append(nodes, x, y, z, qx, qy, qz);

    hash->ids[pos] = new_id;
    hash->count++;

    return new_id;
}

/* -------------------------------------------------------------------------- */
/* Element records                                                             */
/* -------------------------------------------------------------------------- */

static void elem_array_init(ELEM_ARRAY* elems)
{
    elems->count = 0;
    elems->capacity = 1024;
    elems->elem = (ELEM_RECORD*)checked_malloc(
        sizeof(ELEM_RECORD) * (size_t)elems->capacity);
}

static void elem_array_reserve(ELEM_ARRAY* elems, int min_capacity)
{
    if(elems->capacity >= min_capacity) {
        return;
    }

    while(elems->capacity < min_capacity) {
        elems->capacity *= 2;
    }

    elems->elem = (ELEM_RECORD*)checked_realloc(
        elems->elem,
        sizeof(ELEM_RECORD) * (size_t)elems->capacity);
}

static void sort_four_ints(int value[TET4_NUM_NODE])
{
    for(int i = 1; i < TET4_NUM_NODE; i++) {
        int key = value[i];
        int j = i - 1;

        while(j >= 0 && value[j] > key) {
            value[j + 1] = value[j];
            j--;
        }

        value[j + 1] = key;
    }
}

static void elem_array_append(
    ELEM_ARRAY* elems,
    int         n0,
    int         n1,
    int         n2,
    int         n3,
    int         rank)
{
    elem_array_reserve(elems, elems->count + 1);

    ELEM_RECORD* e = &(elems->elem[elems->count]);

    e->conn[0] = n0;
    e->conn[1] = n1;
    e->conn[2] = n2;
    e->conn[3] = n3;
    e->rank = rank;

    for(int i = 0; i < TET4_NUM_NODE; i++) {
        e->key[i] = e->conn[i];
    }

    sort_four_ints(e->key);

    elems->count++;
}

static int compare_elem_record(const void* lhs, const void* rhs)
{
    const ELEM_RECORD* a = (const ELEM_RECORD*)lhs;
    const ELEM_RECORD* b = (const ELEM_RECORD*)rhs;

    for(int i = 0; i < TET4_NUM_NODE; i++) {
        if(a->key[i] < b->key[i]) {
            return -1;
        }

        if(a->key[i] > b->key[i]) {
            return 1;
        }
    }

    if(a->rank < b->rank) {
        return -1;
    }

    if(a->rank > b->rank) {
        return 1;
    }

    return 0;
}

static int same_elem_key(const ELEM_RECORD* a, const ELEM_RECORD* b)
{
    for(int i = 0; i < TET4_NUM_NODE; i++) {
        if(a->key[i] != b->key[i]) {
            return 0;
        }
    }

    return 1;
}

static void vtk_cell_array_init(VTK_CELL_ARRAY* cells)
{
    cells->count = 0;
    cells->capacity = 1024;

    cells->conn = (int*)checked_malloc(
        sizeof(int) * (size_t)cells->capacity * TET4_NUM_NODE);

    cells->cell_value = (double*)checked_malloc(
        sizeof(double) * (size_t)cells->capacity);
}

static void vtk_cell_array_reserve(VTK_CELL_ARRAY* cells, int min_capacity)
{
    if(cells->capacity >= min_capacity) {
        return;
    }

    while(cells->capacity < min_capacity) {
        cells->capacity *= 2;
    }

    cells->conn = (int*)checked_realloc(
        cells->conn,
        sizeof(int) * (size_t)cells->capacity * TET4_NUM_NODE);

    cells->cell_value = (double*)checked_realloc(
        cells->cell_value,
        sizeof(double) * (size_t)cells->capacity);
}

static void vtk_cell_array_append(
    VTK_CELL_ARRAY* cells,
    const int       conn[TET4_NUM_NODE],
    double          value)
{
    vtk_cell_array_reserve(cells, cells->count + 1);

    int offset = cells->count * TET4_NUM_NODE;

    for(int i = 0; i < TET4_NUM_NODE; i++) {
        cells->conn[offset + i] = conn[i];
    }

    cells->cell_value[cells->count] = value;
    cells->count++;
}

static void build_unique_cells(ELEM_ARRAY* records, VTK_CELL_ARRAY* cells)
{
    if(records->count == 0) {
        return;
    }

    qsort(records->elem,
          (size_t)records->count,
          sizeof(ELEM_RECORD),
          compare_elem_record);

    int i = 0;

    while(i < records->count) {
        int j = i + 1;
        int first_rank = records->elem[i].rank;
        int shared_by_multiple_ranks = 0;

        while(j < records->count &&
              same_elem_key(&(records->elem[i]), &(records->elem[j]))) {
            if(records->elem[j].rank != first_rank) {
                shared_by_multiple_ranks = 1;
            }

            j++;
        }

        double cell_value =
            shared_by_multiple_ranks ? -1.0 : (double)(first_rank + 1);

        vtk_cell_array_append(
            cells,
            records->elem[i].conn,
            cell_value);

        i = j;
    }
}

/* -------------------------------------------------------------------------- */
/* Partition files                                                             */
/* -------------------------------------------------------------------------- */

static int read_partition_nodes(
    const SETTINGS* set,
    const char*     filename,
    NODE_ARRAY*     global_nodes,
    NODE_HASH*      node_hash,
    int**           local_to_global_out)
{
    FILE* fp = open_read_file(set, filename);

    int num_local_nodes = 0;

    if(fscanf(fp, "%d", &num_local_nodes) != 1) {
        fprintf(stderr,
                "%s Failed to read node header in \"%s\". Expected first line: num_nodes.\n",
                CODENAME,
                filename);
        exit(EXIT_FAILURE);
    }

    if(num_local_nodes <= 0) {
        fprintf(stderr,
                "%s Invalid number of local nodes in \"%s\": %d.\n",
                CODENAME,
                filename,
                num_local_nodes);
        exit(EXIT_FAILURE);
    }

    int* local_to_global =
        (int*)checked_malloc(sizeof(int) * (size_t)num_local_nodes);

    for(int i = 0; i < num_local_nodes; i++) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        if(fscanf(fp, "%lf %lf %lf", &x, &y, &z) != 3) {
            fprintf(stderr,
                    "%s Failed to read coordinate in \"%s\" at row %d.\n",
                    CODENAME,
                    filename,
                    i);
            exit(EXIT_FAILURE);
        }

        local_to_global[i] = find_or_add_global_node(
            global_nodes,
            node_hash,
            x,
            y,
            z);
    }

    fclose(fp);

    *local_to_global_out = local_to_global;

    return num_local_nodes;
}

static int read_first_nonempty_line(FILE* fp, char* line, int line_size)
{
    while(fgets(line, line_size, fp) != NULL) {
        char* p = line;

        while(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') {
            p++;
        }

        if(*p != '\0') {
            return 1;
        }
    }

    return 0;
}

static void read_partition_header_allow_label(
    FILE*       fp,
    const char* filename,
    int*        num_elems,
    int*        num_local_nodes_per_elem)
{
    char line[BUFFER_SIZE];

    if(!read_first_nonempty_line(fp, line, BUFFER_SIZE)) {
        fprintf(stderr,
                "%s Failed to read element header in \"%s\".\n",
                CODENAME,
                filename);
        exit(EXIT_FAILURE);
    }

    if(sscanf(line, "%d %d", num_elems, num_local_nodes_per_elem) == 2) {
        return;
    }

    /*
     * 1行目が #elem_id などのラベルの場合は読み飛ばし、
     * 次の行を「要素数 節点数」として読む。
     */
    if(!read_first_nonempty_line(fp, line, BUFFER_SIZE)) {
        fprintf(stderr,
                "%s Failed to read element size line after label in \"%s\".\n",
                CODENAME,
                filename);
        exit(EXIT_FAILURE);
    }

    if(sscanf(line, "%d %d", num_elems, num_local_nodes_per_elem) != 2) {
        fprintf(stderr,
                "%s Invalid element size line in \"%s\": %s\n",
                CODENAME,
                filename,
                line);
        exit(EXIT_FAILURE);
    }
}

static int infer_index_base(
    INDEX_BASE_MODE mode,
    int             num_local_nodes,
    int             min_index,
    int             max_index,
    const char*     filename)
{
    if(mode == INDEX_BASE_0) {
        return 0;
    }

    if(mode == INDEX_BASE_1) {
        return 1;
    }

    if(min_index == 0) {
        return 0;
    }

    if(max_index == num_local_nodes) {
        return 1;
    }

    if(max_index > num_local_nodes) {
        fprintf(stderr,
                "%s Invalid node index in \"%s\": max index = %d, num_local_nodes = %d.\n",
                CODENAME,
                filename,
                max_index,
                num_local_nodes);
        exit(EXIT_FAILURE);
    }

    /*
     * 0-based/1-based の判定が曖昧な場合。
     * C系の出力を想定して0-basedを既定にする。
     * 1-basedなら実行時に -base 1 を指定する。
     */
    printf("%s Index base in \"%s\" is ambiguous; using 0-based. "
           "Use -base 1 if your local element files are 1-based.\n",
           CODENAME,
           filename);

    return 0;
}

static void append_partition_elements(
    const SETTINGS* set,
    const char*     filename,
    int             rank,
    const int*      local_to_global,
    int             num_local_nodes,
    ELEM_ARRAY*     elem_records)
{
    FILE* fp = open_read_file(set, filename);

    int num_elems = 0;
    int file_num_local_nodes_per_elem = 0;

    read_partition_header_allow_label(
        fp,
        filename,
        &num_elems,
        &file_num_local_nodes_per_elem);

    if(num_elems <= 0 || file_num_local_nodes_per_elem < TET4_NUM_NODE) {
        fprintf(stderr,
                "%s Invalid element header in \"%s\": num_elems = %d, num_nodes_per_elem = %d.\n",
                CODENAME,
                filename,
                num_elems,
                file_num_local_nodes_per_elem);
        exit(EXIT_FAILURE);
    }

    /*
     * 各行のうち先頭4つだけ保存する。
     * 5列目以降は読み捨てる。
     */
    int* raw = (int*)checked_malloc(
        sizeof(int) * (size_t)num_elems * TET4_NUM_NODE);

    int min_index = 2147483647;
    int max_index = -2147483647;

    for(int i = 0; i < num_elems; i++) {
        for(int j = 0; j < file_num_local_nodes_per_elem; j++) {
            int value = 0;

            if(fscanf(fp, "%d", &value) != 1) {
                fprintf(stderr,
                        "%s Failed to read element in \"%s\" at row %d, col %d.\n",
                        CODENAME,
                        filename,
                        i,
                        j);
                exit(EXIT_FAILURE);
            }

            if(j < TET4_NUM_NODE) {
                raw[i * TET4_NUM_NODE + j] = value;

                if(value < min_index) {
                    min_index = value;
                }

                if(value > max_index) {
                    max_index = value;
                }
            }
        }
    }

    fclose(fp);

    int index_base = infer_index_base(
        set->index_base_mode,
        num_local_nodes,
        min_index,
        max_index,
        filename);

    for(int i = 0; i < num_elems; i++) {
        int g[TET4_NUM_NODE];

        for(int j = 0; j < TET4_NUM_NODE; j++) {
            int local_index = raw[i * TET4_NUM_NODE + j] - index_base;

            if(local_index < 0 || local_index >= num_local_nodes) {
                fprintf(stderr,
                        "%s Local node index out of range in \"%s\": row = %d, col = %d, "
                        "value = %d, index_base = %d, num_local_nodes = %d.\n",
                        CODENAME,
                        filename,
                        i,
                        j,
                        raw[i * TET4_NUM_NODE + j],
                        index_base,
                        num_local_nodes);
                exit(EXIT_FAILURE);
            }

            g[j] = local_to_global[local_index];
        }

        elem_array_append(
            elem_records,
            g[0],
            g[1],
            g[2],
            g[3],
            rank);
    }

    free(raw);
}

static void read_all_partition_files(
    const SETTINGS* set,
    int             num_parallel,
    NODE_ARRAY*     global_nodes,
    NODE_HASH*      node_hash,
    ELEM_ARRAY*     elem_records)
{
    long long total_local_nodes_read = 0;
    long long total_local_elem_records_read = 0;

    for(int rank = 0; rank < num_parallel; rank++) {
        char node_filename[BUFFER_SIZE];
        char elem_filename[BUFFER_SIZE];

        make_rank_filename(
            node_filename,
            BUFFER_SIZE,
            set->part_node_body,
            set,
            rank);

        make_rank_filename(
            elem_filename,
            BUFFER_SIZE,
            set->part_elem_body,
            set,
            rank);

        int* local_to_global = NULL;

        int num_local_nodes = read_partition_nodes(
            set,
            node_filename,
            global_nodes,
            node_hash,
            &local_to_global);

        int before_elem_count = elem_records->count;

        append_partition_elements(
            set,
            elem_filename,
            rank,
            local_to_global,
            num_local_nodes,
            elem_records);

        total_local_nodes_read += (long long)num_local_nodes;
        total_local_elem_records_read +=
            (long long)(elem_records->count - before_elem_count);

        free(local_to_global);
    }

    printf("%s Number of local nodes read          : %lld\n",
           CODENAME,
           total_local_nodes_read);

    printf("%s Number of unique VTK points         : %d\n",
           CODENAME,
           global_nodes->count);

    printf("%s Number of local element records read: %lld\n",
           CODENAME,
           total_local_elem_records_read);
}

/* -------------------------------------------------------------------------- */
/* VTK output                                                                  */
/* -------------------------------------------------------------------------- */

static void write_vtk(
    const SETTINGS*       set,
    const NODE_ARRAY*     nodes,
    const VTK_CELL_ARRAY* cells)
{
    FILE* fp = open_write_file(set, set->outfile_vtk);

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtkddm_tet4 generated from partitioned node/element files\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d double\n", nodes->count);

    for(int i = 0; i < nodes->count; i++) {
        fprintf(fp,
                "%.15e %.15e %.15e\n",
                nodes->x[i],
                nodes->y[i],
                nodes->z[i]);
    }

    fprintf(fp,
            "CELLS %d %lld\n",
            cells->count,
            (long long)cells->count * (long long)(TET4_NUM_NODE + 1));

    for(int i = 0; i < cells->count; i++) {
        int offset = i * TET4_NUM_NODE;

        fprintf(fp,
                "%d %d %d %d %d\n",
                TET4_NUM_NODE,
                cells->conn[offset + 0],
                cells->conn[offset + 1],
                cells->conn[offset + 2],
                cells->conn[offset + 3]);
    }

    fprintf(fp, "CELL_TYPES %d\n", cells->count);

    for(int i = 0; i < cells->count; i++) {
        fprintf(fp, "%d\n", VTK_TETRA);
    }

    fprintf(fp, "CELL_DATA %d\n", cells->count);
    fprintf(fp, "SCALARS %s double 1\n", VTK_CELL_DATA_LABEL);
    fprintf(fp, "LOOKUP_TABLE default\n");

    for(int i = 0; i < cells->count; i++) {
        fprintf(fp, "%.15e\n", cells->cell_value[i]);
    }

    fclose(fp);
}

/* -------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
    printf("\n");

    SETTINGS set;
    int num_parallel = args_manager(&set, argc, argv);

    if(num_parallel <= 0) {
        fprintf(stderr,
                "%s Invalid num_parallel = %d. Please pass a positive number.\n",
                CODENAME,
                num_parallel);
        return EXIT_FAILURE;
    }

    NODE_ARRAY global_nodes;
    NODE_HASH node_hash;
    ELEM_ARRAY elem_records;
    VTK_CELL_ARRAY cells;

    node_array_init(&global_nodes);
    node_hash_init(&node_hash, 2048);
    elem_array_init(&elem_records);
    vtk_cell_array_init(&cells);

    read_all_partition_files(
        &set,
        num_parallel,
        &global_nodes,
        &node_hash,
        &elem_records);

    build_unique_cells(
        &elem_records,
        &cells);

    printf("%s Number of unique VTK cells          : %d\n",
           CODENAME,
           cells.count);

    write_vtk(
        &set,
        &global_nodes,
        &cells);

    free(global_nodes.x);
    free(global_nodes.y);
    free(global_nodes.z);
    free(global_nodes.qx);
    free(global_nodes.qy);
    free(global_nodes.qz);

    node_hash_free(&node_hash);

    free(elem_records.elem);
    free(cells.conn);
    free(cells.cell_value);

    return EXIT_SUCCESS;
}