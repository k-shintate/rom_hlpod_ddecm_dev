/*
 * st_partition_postprocessor.c
 *
 * グローバル時空間（ST）グラフの分割結果 owner[global_st_id] から、
 * 各 MPI ランク用の以下のファイルを生成する後処理プログラム。
 *
 *   rank_XXXXXX/st_node_map.dat
 *   rank_XXXXXX/graph.dat
 *   rank_XXXXXX/st_elem.dat
 *   rank_XXXXXX/st_elem_map.dat
 *   rank_XXXXXX/st_meta.dat
 *
 * ローカル ST 節点番号は、所有節点（internal）を先頭、ハロー節点を
 * 後方に並べる。ランク r のローカル ST セルは、セル内に rank r が
 * 所有する ST 節点を1個以上含むセルとする。そのため境界セルは複数
 * ランクに重複して配置され得る。
 *
 * owner.dat の対応形式（空行と # 以降は無視）:
 *
 *   A: N 行の owner ベクトル
 *      owner(gid=0)
 *      owner(gid=1)
 *      ...
 *
 *   B: ヘッダ付き owner ベクトル
 *      N
 *      owner(gid=0)
 *      ...
 *
 *   C: N 行の gid-owner ペア
 *      gid owner
 *      ...
 *
 *   D: ヘッダ付き gid-owner ペア
 *      N NUM_PARTS
 *      gid owner
 *      ...
 *
 * ビルド:
 *   cc -O2 -std=c11 -Wall -Wextra -pedantic \
 *      st_partition_postprocessor.c -o st_partition_postprocessor
 *
 * 実行:
 *   ./st_partition_postprocessor NUM_PARTS \
 *      st_meta.dat graph.dat st_elem.dat st_elem_map.dat \
 *      owner.dat OUTPUT_DIR
 */

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
    int degree;
    int num_time_dofs;
    int num_window_slabs;
    int num_global_space_nodes;
    int num_global_space_elements;
    int num_space_nodes_per_element;
    int num_st_nodes_per_cell;
    int64_t num_global_st_nodes;
    int64_t num_global_st_cells;
} ST_META;

typedef struct {
    int64_t num_nodes;
    int64_t* index;
    int* item;
    int64_t num_items;
} GLOBAL_GRAPH;

typedef struct {
    int64_t num_cells;
    int width;
    int64_t* conn;
} GLOBAL_ST_ELEM;

typedef struct {
    int64_t num_cells;
    int num_space_nodes;
    int64_t* global_cell;
    int* slab;
    int* global_space_elem;
    int* global_space_conn;
} GLOBAL_ST_ELEM_MAP;

typedef struct {
    int64_t* data;
    size_t size;
    size_t capacity;
} I64_VECTOR;

typedef struct {
    int* data;
    size_t size;
    size_t capacity;
} INT_VECTOR;

typedef struct {
    int64_t a;
    int64_t b;
    int count;
} OWNER_LINE;

typedef struct {
    OWNER_LINE* data;
    size_t size;
    size_t capacity;
} OWNER_LINES;

static void fail(const char* message)
{
    fprintf(stderr, "ERROR: %s\n", message);
    exit(EXIT_FAILURE);
}

static void fail_path(const char* action, const char* path)
{
    fprintf(stderr, "ERROR: %s: %s (%s)\n", action, path, strerror(errno));
    exit(EXIT_FAILURE);
}

static void* checked_calloc(size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        fail("allocation size overflow");
    }
    void* p = calloc(count, size);
    if(p == NULL) {
        fail("memory allocation failed");
    }
    return p;
}

static void* checked_malloc(size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        fail("allocation size overflow");
    }
    void* p = malloc(count * size);
    if(p == NULL) {
        fail("memory allocation failed");
    }
    return p;
}

static void* checked_realloc(void* ptr, size_t count, size_t size)
{
    if(size != 0 && count > SIZE_MAX / size) {
        fail("reallocation size overflow");
    }
    void* p = realloc(ptr, count * size);
    if(p == NULL) {
        fail("memory reallocation failed");
    }
    return p;
}

static int compare_i64(const void* lhs, const void* rhs)
{
    const int64_t a = *(const int64_t*)lhs;
    const int64_t b = *(const int64_t*)rhs;
    return (a > b) - (a < b);
}

static int compare_int(const void* lhs, const void* rhs)
{
    const int a = *(const int*)lhs;
    const int b = *(const int*)rhs;
    return (a > b) - (a < b);
}

static void i64_vector_push(I64_VECTOR* v, int64_t value)
{
    if(v->size == v->capacity) {
        const size_t next = (v->capacity == 0) ? 16 : 2 * v->capacity;
        if(next < v->capacity) {
            fail("int64 vector capacity overflow");
        }
        v->data = (int64_t*)checked_realloc(v->data, next, sizeof(int64_t));
        v->capacity = next;
    }
    v->data[v->size++] = value;
}

static void int_vector_push(INT_VECTOR* v, int value)
{
    if(v->size == v->capacity) {
        const size_t next = (v->capacity == 0) ? 16 : 2 * v->capacity;
        if(next < v->capacity) {
            fail("integer vector capacity overflow");
        }
        v->data = (int*)checked_realloc(v->data, next, sizeof(int));
        v->capacity = next;
    }
    v->data[v->size++] = value;
}

static void owner_lines_push(OWNER_LINES* v, OWNER_LINE value)
{
    if(v->size == v->capacity) {
        const size_t next = (v->capacity == 0) ? 64 : 2 * v->capacity;
        if(next < v->capacity) {
            fail("owner-line capacity overflow");
        }
        v->data = (OWNER_LINE*)checked_realloc(
            v->data, next, sizeof(OWNER_LINE));
        v->capacity = next;
    }
    v->data[v->size++] = value;
}

static void i64_vector_sort_unique(I64_VECTOR* v)
{
    if(v->size == 0) {
        return;
    }
    qsort(v->data, v->size, sizeof(int64_t), compare_i64);
    size_t write = 1;
    for(size_t read = 1; read < v->size; read++) {
        if(v->data[read] != v->data[write - 1]) {
            v->data[write++] = v->data[read];
        }
    }
    v->size = write;
}

static void int_vector_sort_unique(INT_VECTOR* v)
{
    if(v->size == 0) {
        return;
    }
    qsort(v->data, v->size, sizeof(int), compare_int);
    size_t write = 1;
    for(size_t read = 1; read < v->size; read++) {
        if(v->data[read] != v->data[write - 1]) {
            v->data[write++] = v->data[read];
        }
    }
    v->size = write;
}

static void make_directory(const char* path)
{
    if(mkdir(path, 0777) != 0 && errno != EEXIST) {
        fail_path("cannot create directory", path);
    }
}

static void make_directories_recursive(const char* path)
{
    char buffer[PATH_MAX];
    const size_t length = strlen(path);
    if(length == 0 || length >= sizeof(buffer)) {
        fail("output directory path is empty or too long");
    }
    memcpy(buffer, path, length + 1);

    for(char* p = buffer + 1; *p != '\0'; p++) {
        if(*p == '/') {
            *p = '\0';
            make_directory(buffer);
            *p = '/';
        }
    }
    make_directory(buffer);
}

static void make_path(
    char* path,
    size_t path_size,
    const char* directory,
    const char* filename)
{
    const int n = snprintf(path, path_size, "%s/%s", directory, filename);
    if(n < 0 || (size_t)n >= path_size) {
        fail("path is too long");
    }
}

static FILE* open_read(const char* path)
{
    FILE* fp = fopen(path, "r");
    if(fp == NULL) {
        fail_path("cannot open input file", path);
    }
    return fp;
}

static FILE* open_write(const char* path)
{
    FILE* fp = fopen(path, "w");
    if(fp == NULL) {
        fail_path("cannot open output file", path);
    }
    (void)setvbuf(fp, NULL, _IOFBF, 1U << 20);
    return fp;
}

static int64_t expected_st_gid(
    int slab,
    int alpha,
    int global_space_node,
    const ST_META* meta)
{
    return
        ((int64_t)slab * meta->num_time_dofs + alpha)
        * (int64_t)meta->num_global_space_nodes
        + global_space_node;
}

static void decode_st_gid(
    int64_t gid,
    const ST_META* meta,
    int* slab,
    int* alpha,
    int* global_space_node)
{
    if(gid < 0 || gid >= meta->num_global_st_nodes) {
        fail("global ST node ID is out of range while decoding");
    }
    const int64_t nx = meta->num_global_space_nodes;
    *global_space_node = (int)(gid % nx);
    const int64_t temporal = gid / nx;
    *alpha = (int)(temporal % meta->num_time_dofs);
    *slab = (int)(temporal / meta->num_time_dofs);
}

static void read_meta(const char* path, ST_META* meta)
{
    memset(meta, 0, sizeof(*meta));
    FILE* fp = open_read(path);

    char key[128];
    char value[512];
    while(fscanf(fp, "%127s %511s", key, value) == 2) {
        if(strcmp(key, "degree") == 0) {
            meta->degree = atoi(value);
        }
        else if(strcmp(key, "num_time_dofs") == 0) {
            meta->num_time_dofs = atoi(value);
        }
        else if(strcmp(key, "num_window_slabs") == 0) {
            meta->num_window_slabs = atoi(value);
        }
        else if(strcmp(key, "num_global_space_nodes") == 0) {
            meta->num_global_space_nodes = atoi(value);
        }
        else if(strcmp(key, "num_global_space_elements") == 0) {
            meta->num_global_space_elements = atoi(value);
        }
        else if(strcmp(key, "num_space_nodes_per_element") == 0) {
            meta->num_space_nodes_per_element = atoi(value);
        }
        else if(strcmp(key, "num_st_nodes_per_cell") == 0) {
            meta->num_st_nodes_per_cell = atoi(value);
        }
        else if(strcmp(key, "num_global_st_nodes") == 0) {
            meta->num_global_st_nodes = strtoll(value, NULL, 10);
        }
        else if(strcmp(key, "num_global_st_cells") == 0) {
            meta->num_global_st_cells = strtoll(value, NULL, 10);
        }
    }
    fclose(fp);

    if(meta->num_time_dofs <= 0 ||
       meta->num_window_slabs <= 0 ||
       meta->num_global_space_nodes <= 0 ||
       meta->num_global_space_elements <= 0 ||
       meta->num_space_nodes_per_element <= 0 ||
       meta->num_global_st_nodes <= 0 ||
       meta->num_global_st_cells <= 0) {
        fail("invalid or incomplete st_meta.dat");
    }

    const int expected_width =
        meta->num_time_dofs * meta->num_space_nodes_per_element;
    if(meta->num_st_nodes_per_cell == 0) {
        meta->num_st_nodes_per_cell = expected_width;
    }
    if(meta->num_st_nodes_per_cell != expected_width) {
        fail("num_st_nodes_per_cell does not equal nt * spatial element width");
    }

    const int64_t expected_nodes =
        (int64_t)meta->num_window_slabs
        * meta->num_time_dofs
        * meta->num_global_space_nodes;
    const int64_t expected_cells =
        (int64_t)meta->num_window_slabs
        * meta->num_global_space_elements;

    if(meta->num_global_st_nodes != expected_nodes) {
        fail("num_global_st_nodes is inconsistent with metadata");
    }
    if(meta->num_global_st_cells != expected_cells) {
        fail("num_global_st_cells is inconsistent with metadata");
    }
    if(meta->num_global_st_nodes > INT_MAX) {
        fail("num_global_st_nodes exceeds INT_MAX used by graph files");
    }
}

static void read_global_graph(
    const char* path,
    const ST_META* meta,
    GLOBAL_GRAPH* graph)
{
    memset(graph, 0, sizeof(*graph));
    FILE* fp = open_read(path);

    int64_t num_nodes = -1;
    if(fscanf(fp, "%" SCNd64, &num_nodes) != 1 ||
       num_nodes != meta->num_global_st_nodes) {
        fclose(fp);
        fail("graph.dat node count is inconsistent with st_meta.dat");
    }

    INT_VECTOR* rows = (INT_VECTOR*)checked_calloc(
        (size_t)num_nodes, sizeof(INT_VECTOR));
    unsigned char* seen = (unsigned char*)checked_calloc(
        (size_t)num_nodes, sizeof(unsigned char));

    for(int64_t k = 0; k < num_nodes; k++) {
        int64_t row_id = -1;
        int64_t degree = -1;
        if(fscanf(fp, "%" SCNd64 " %" SCNd64, &row_id, &degree) != 2 ||
           row_id < 0 || row_id >= num_nodes || degree < 0 ||
           degree > INT_MAX || seen[row_id]) {
            fclose(fp);
            fail("invalid graph.dat row header");
        }
        seen[row_id] = 1;

        for(int64_t j = 0; j < degree; j++) {
            int64_t col = -1;
            if(fscanf(fp, "%" SCNd64, &col) != 1 ||
               col < 0 || col >= num_nodes || col > INT_MAX) {
                fclose(fp);
                fail("invalid graph.dat column ID");
            }
            if(col != row_id) {
                int_vector_push(&rows[row_id], (int)col);
            }
        }
    }
    fclose(fp);

    int64_t total_items = 0;
    for(int64_t i = 0; i < num_nodes; i++) {
        if(!seen[i]) {
            fail("graph.dat has a missing row");
        }
        int_vector_sort_unique(&rows[i]);
        if((int64_t)rows[i].size > INT64_MAX - total_items) {
            fail("graph item count overflow");
        }
        total_items += (int64_t)rows[i].size;
    }

    graph->num_nodes = num_nodes;
    graph->num_items = total_items;
    graph->index = (int64_t*)checked_calloc(
        (size_t)num_nodes + 1, sizeof(int64_t));
    graph->item = (int*)checked_malloc(
        (size_t)total_items, sizeof(int));

    int64_t position = 0;
    for(int64_t i = 0; i < num_nodes; i++) {
        graph->index[i] = position;
        for(size_t j = 0; j < rows[i].size; j++) {
            graph->item[position++] = rows[i].data[j];
        }
        free(rows[i].data);
    }
    graph->index[num_nodes] = position;

    free(rows);
    free(seen);
}

static void read_global_st_elem(
    const char* path,
    const ST_META* meta,
    GLOBAL_ST_ELEM* elem)
{
    memset(elem, 0, sizeof(*elem));
    FILE* fp = open_read(path);

    int64_t num_cells = -1;
    int width = -1;
    if(fscanf(fp, "%" SCNd64 " %d", &num_cells, &width) != 2 ||
       num_cells != meta->num_global_st_cells ||
       width != meta->num_st_nodes_per_cell) {
        fclose(fp);
        fail("st_elem.dat header is inconsistent with st_meta.dat");
    }

    if((uint64_t)num_cells > SIZE_MAX / (size_t)width) {
        fclose(fp);
        fail("st_elem.dat connectivity size overflow");
    }
    const size_t total = (size_t)num_cells * (size_t)width;
    elem->conn = (int64_t*)checked_malloc(total, sizeof(int64_t));

    for(size_t i = 0; i < total; i++) {
        if(fscanf(fp, "%" SCNd64, &elem->conn[i]) != 1 ||
           elem->conn[i] < 0 ||
           elem->conn[i] >= meta->num_global_st_nodes) {
            fclose(fp);
            fail("invalid global ST node ID in st_elem.dat");
        }
    }
    fclose(fp);

    elem->num_cells = num_cells;
    elem->width = width;
}

static void read_global_st_elem_map(
    const char* path,
    const ST_META* meta,
    const GLOBAL_ST_ELEM* elem,
    GLOBAL_ST_ELEM_MAP* map)
{
    memset(map, 0, sizeof(*map));
    FILE* fp = open_read(path);

    int64_t num_cells = -1;
    int num_space_nodes = -1;
    if(fscanf(fp, "%" SCNd64 " %d", &num_cells, &num_space_nodes) != 2 ||
       num_cells != elem->num_cells ||
       num_space_nodes != meta->num_space_nodes_per_element) {
        fclose(fp);
        fail("st_elem_map.dat header is inconsistent");
    }

    map->num_cells = num_cells;
    map->num_space_nodes = num_space_nodes;
    map->global_cell = (int64_t*)checked_malloc(
        (size_t)num_cells, sizeof(int64_t));
    map->slab = (int*)checked_malloc((size_t)num_cells, sizeof(int));
    map->global_space_elem = (int*)checked_malloc(
        (size_t)num_cells, sizeof(int));
    map->global_space_conn = (int*)checked_malloc(
        (size_t)num_cells * (size_t)num_space_nodes, sizeof(int));

    unsigned char* local_seen = (unsigned char*)checked_calloc(
        (size_t)num_cells, sizeof(unsigned char));
    unsigned char* global_seen = (unsigned char*)checked_calloc(
        (size_t)num_cells, sizeof(unsigned char));

    for(int64_t row = 0; row < num_cells; row++) {
        int64_t local_cell = -1;
        int64_t global_cell = -1;
        int slab = -1;
        int global_elem = -1;

        if(fscanf(
               fp,
               "%" SCNd64 " %" SCNd64 " %d %d",
               &local_cell,
               &global_cell,
               &slab,
               &global_elem) != 4 ||
           local_cell < 0 || local_cell >= num_cells ||
           global_cell < 0 || global_cell >= num_cells ||
           slab < 0 || slab >= meta->num_window_slabs ||
           global_elem < 0 ||
           global_elem >= meta->num_global_space_elements ||
           local_seen[local_cell] || global_seen[global_cell]) {
            fclose(fp);
            fail("invalid or duplicated row in st_elem_map.dat");
        }

        const int64_t expected_cell =
            (int64_t)slab * meta->num_global_space_elements + global_elem;
        if(global_cell != expected_cell) {
            fclose(fp);
            fail("global ST cell ID does not match slab and spatial element");
        }

        local_seen[local_cell] = 1;
        global_seen[global_cell] = 1;
        map->global_cell[local_cell] = global_cell;
        map->slab[local_cell] = slab;
        map->global_space_elem[local_cell] = global_elem;

        for(int a = 0; a < num_space_nodes; a++) {
            int global_space_node = -1;
            if(fscanf(fp, "%d", &global_space_node) != 1 ||
               global_space_node < 0 ||
               global_space_node >= meta->num_global_space_nodes) {
                fclose(fp);
                fail("invalid spatial node in st_elem_map.dat");
            }
            map->global_space_conn[
                (size_t)local_cell * (size_t)num_space_nodes + (size_t)a
            ] = global_space_node;
        }
    }
    fclose(fp);

    for(int64_t c = 0; c < num_cells; c++) {
        if(!local_seen[c]) {
            fail("st_elem_map.dat has a missing local cell row");
        }
        const int slab = map->slab[c];
        for(int alpha = 0; alpha < meta->num_time_dofs; alpha++) {
            for(int a = 0; a < num_space_nodes; a++) {
                const int global_space_node = map->global_space_conn[
                    (size_t)c * (size_t)num_space_nodes + (size_t)a
                ];
                const int64_t expected = expected_st_gid(
                    slab, alpha, global_space_node, meta);
                const int64_t actual = elem->conn[
                    (size_t)c * (size_t)elem->width
                    + (size_t)alpha * (size_t)num_space_nodes
                    + (size_t)a
                ];
                if(actual != expected) {
                    fprintf(
                        stderr,
                        "ERROR: global elem/map mismatch: local_cell=%" PRId64
                        " alpha=%d node=%d actual=%" PRId64
                        " expected=%" PRId64 "\n",
                        c, alpha, a, actual, expected);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    free(local_seen);
    free(global_seen);
}

static char* trim_leading_space(char* s)
{
    while(*s != '\0' && isspace((unsigned char)*s)) {
        s++;
    }
    return s;
}

static void read_owner_lines(const char* path, OWNER_LINES* lines)
{
    memset(lines, 0, sizeof(*lines));
    FILE* fp = open_read(path);
    char buffer[4096];
    int line_number = 0;

    while(fgets(buffer, sizeof(buffer), fp) != NULL) {
        line_number++;
        char* comment = strchr(buffer, '#');
        if(comment != NULL) {
            *comment = '\0';
        }
        char* p = trim_leading_space(buffer);
        if(*p == '\0') {
            continue;
        }

        errno = 0;
        char* end = NULL;
        const int64_t a = strtoll(p, &end, 10);
        if(errno != 0 || end == p) {
            fprintf(stderr, "ERROR: invalid owner.dat line %d\n", line_number);
            exit(EXIT_FAILURE);
        }

        p = trim_leading_space(end);
        OWNER_LINE line;
        line.a = a;
        line.b = 0;
        line.count = 1;

        if(*p != '\0') {
            errno = 0;
            const int64_t b = strtoll(p, &end, 10);
            if(errno != 0 || end == p) {
                fprintf(stderr, "ERROR: invalid owner.dat line %d\n", line_number);
                exit(EXIT_FAILURE);
            }
            p = trim_leading_space(end);
            if(*p != '\0') {
                fprintf(
                    stderr,
                    "ERROR: owner.dat line %d has more than two integers\n",
                    line_number);
                exit(EXIT_FAILURE);
            }
            line.b = b;
            line.count = 2;
        }
        owner_lines_push(lines, line);
    }
    fclose(fp);
}

static int* read_owner_file(
    const char* path,
    const ST_META* meta,
    int num_parts)
{
    OWNER_LINES lines;
    read_owner_lines(path, &lines);

    const size_t n = (size_t)meta->num_global_st_nodes;
    size_t start = 0;
    int vector_format = 0;
    int pair_format = 0;

    if(lines.size == n) {
        int all_one = 1;
        int all_two = 1;
        for(size_t i = 0; i < lines.size; i++) {
            all_one = all_one && lines.data[i].count == 1;
            all_two = all_two && lines.data[i].count == 2;
        }
        vector_format = all_one;
        pair_format = all_two;
    }
    else if(lines.size == n + 1) {
        const OWNER_LINE header = lines.data[0];
        if(header.count == 1 && header.a == meta->num_global_st_nodes) {
            vector_format = 1;
            start = 1;
            for(size_t i = start; i < lines.size; i++) {
                if(lines.data[i].count != 1) {
                    vector_format = 0;
                    break;
                }
            }
        }
        else if(header.count == 2 &&
                header.a == meta->num_global_st_nodes &&
                header.b == num_parts) {
            pair_format = 1;
            start = 1;
            for(size_t i = start; i < lines.size; i++) {
                if(lines.data[i].count != 2) {
                    pair_format = 0;
                    break;
                }
            }
        }
    }

    if(!vector_format && !pair_format) {
        free(lines.data);
        fail(
            "owner.dat format is not recognized; expected N owners or N gid-owner pairs");
    }

    int* owner = (int*)checked_malloc(n, sizeof(int));
    for(size_t i = 0; i < n; i++) {
        owner[i] = -1;
    }

    if(vector_format) {
        for(size_t gid = 0; gid < n; gid++) {
            const int64_t rank64 = lines.data[start + gid].a;
            if(rank64 < 0 || rank64 >= num_parts) {
                free(lines.data);
                free(owner);
                fail("owner rank is out of range in owner.dat");
            }
            owner[gid] = (int)rank64;
        }
    }
    else {
        for(size_t i = start; i < lines.size; i++) {
            const int64_t gid64 = lines.data[i].a;
            const int64_t rank64 = lines.data[i].b;
            if(gid64 < 0 || gid64 >= meta->num_global_st_nodes ||
               rank64 < 0 || rank64 >= num_parts) {
                free(lines.data);
                free(owner);
                fail("gid or owner rank is out of range in owner.dat");
            }
            if(owner[gid64] != -1) {
                free(lines.data);
                free(owner);
                fail("duplicated gid in owner.dat");
            }
            owner[gid64] = (int)rank64;
        }
    }

    for(size_t gid = 0; gid < n; gid++) {
        if(owner[gid] < 0) {
            free(lines.data);
            free(owner);
            fail("owner.dat does not define every global ST node");
        }
    }

    free(lines.data);
    return owner;
}

static void copy_file(const char* source, const char* destination)
{
    FILE* in = fopen(source, "rb");
    if(in == NULL) {
        fail_path("cannot open source file for copying", source);
    }
    FILE* out = fopen(destination, "wb");
    if(out == NULL) {
        fclose(in);
        fail_path("cannot open destination file for copying", destination);
    }

    unsigned char buffer[1U << 16];
    size_t count;
    while((count = fread(buffer, 1, sizeof(buffer), in)) > 0) {
        if(fwrite(buffer, 1, count, out) != count) {
            fclose(in);
            fclose(out);
            fail_path("failed while copying to", destination);
        }
    }
    if(ferror(in)) {
        fclose(in);
        fclose(out);
        fail_path("failed while reading", source);
    }
    fclose(in);
    fclose(out);
}

static void build_rank_sets(
    int rank,
    const ST_META* meta,
    const GLOBAL_GRAPH* graph,
    const GLOBAL_ST_ELEM* elem,
    const int* owner,
    I64_VECTOR* internal,
    I64_VECTOR* halo,
    I64_VECTOR* local_cells,
    unsigned char* local_mark)
{
    memset(internal, 0, sizeof(*internal));
    memset(halo, 0, sizeof(*halo));
    memset(local_cells, 0, sizeof(*local_cells));
    memset(local_mark, 0, (size_t)meta->num_global_st_nodes);

    for(int64_t gid = 0; gid < meta->num_global_st_nodes; gid++) {
        if(owner[gid] == rank) {
            i64_vector_push(internal, gid);
            local_mark[gid] = 1;
        }
    }

    /*
     * 所有 ST 節点を1個以上含むセルを選び、そのセルの全 ST 節点を
     * ローカル節点集合へ追加する。
     */
    for(int64_t c = 0; c < elem->num_cells; c++) {
        int required = 0;
        const int64_t* conn = &elem->conn[(size_t)c * (size_t)elem->width];
        for(int a = 0; a < elem->width; a++) {
            if(owner[conn[a]] == rank) {
                required = 1;
                break;
            }
        }
        if(required) {
            i64_vector_push(local_cells, c);
            for(int a = 0; a < elem->width; a++) {
                local_mark[conn[a]] = 1;
            }
        }
    }

    /* 所有行が参照するグラフ隣接節点を全てハロー候補へ加える。 */
    for(size_t i = 0; i < internal->size; i++) {
        const int64_t row = internal->data[i];
        for(int64_t p = graph->index[row]; p < graph->index[row + 1]; p++) {
            local_mark[graph->item[p]] = 1;
        }
    }

    for(int64_t gid = 0; gid < meta->num_global_st_nodes; gid++) {
        if(local_mark[gid] && owner[gid] != rank) {
            i64_vector_push(halo, gid);
        }
    }

    i64_vector_sort_unique(internal);
    i64_vector_sort_unique(halo);
    i64_vector_sort_unique(local_cells);
}

static void write_st_node_map(
    const char* rank_dir,
    const ST_META* meta,
    const I64_VECTOR* internal,
    const I64_VECTOR* halo,
    int* global_to_local,
    int64_t* local_to_global)
{
    const size_t num_internal = internal->size;
    const size_t num_total = internal->size + halo->size;
    if(num_total > INT_MAX) {
        fail("rank-local ST node count exceeds INT_MAX");
    }

    for(int64_t gid = 0; gid < meta->num_global_st_nodes; gid++) {
        global_to_local[gid] = -1;
    }

    size_t position = 0;
    for(size_t i = 0; i < internal->size; i++) {
        const int64_t gid = internal->data[i];
        global_to_local[gid] = (int)position;
        local_to_global[position++] = gid;
    }
    for(size_t i = 0; i < halo->size; i++) {
        const int64_t gid = halo->data[i];
        if(global_to_local[gid] >= 0) {
            fail("internal and halo ST node sets overlap");
        }
        global_to_local[gid] = (int)position;
        local_to_global[position++] = gid;
    }

    char path[PATH_MAX];
    make_path(path, sizeof(path), rank_dir, "st_node_map.dat");
    FILE* fp = open_write(path);
    fprintf(fp, "%zu %zu\n", num_internal, num_total);

    for(size_t local = 0; local < num_total; local++) {
        const int64_t gid = local_to_global[local];
        int slab = -1;
        int alpha = -1;
        int global_space_node = -1;
        decode_st_gid(
            gid, meta, &slab, &alpha, &global_space_node);
        fprintf(
            fp,
            "%zu %" PRId64 " %d %d %d\n",
            local,
            gid,
            slab,
            alpha,
            global_space_node);
    }
    fclose(fp);
}

static void write_local_graph(
    const char* rank_dir,
    const GLOBAL_GRAPH* graph,
    size_t num_internal,
    size_t num_total,
    const int64_t* local_to_global,
    const int* global_to_local)
{
    char path[PATH_MAX];
    make_path(path, sizeof(path), rank_dir, "graph.dat");
    FILE* fp = open_write(path);
    fprintf(fp, "%zu\n", num_total);

    INT_VECTOR columns = {0};
    for(size_t local_row = 0; local_row < num_total; local_row++) {
        columns.size = 0;
        const int64_t global_row = local_to_global[local_row];
        for(int64_t p = graph->index[global_row];
            p < graph->index[global_row + 1];
            p++) {
            const int global_col = graph->item[p];
            const int local_col = global_to_local[global_col];
            if(local_col >= 0 && (size_t)local_col != local_row) {
                int_vector_push(&columns, local_col);
            }
            else if(local_row < num_internal && local_col < 0) {
                fprintf(
                    stderr,
                    "ERROR: internal graph row %zu (global=%" PRId64
                    ") references missing global column %d\n",
                    local_row,
                    global_row,
                    global_col);
                exit(EXIT_FAILURE);
            }
        }
        int_vector_sort_unique(&columns);

        fprintf(fp, "%zu %zu", local_row, columns.size);
        for(size_t j = 0; j < columns.size; j++) {
            fprintf(fp, " %d", columns.data[j]);
        }
        fputc('\n', fp);
    }

    free(columns.data);
    fclose(fp);
}

static void write_local_cells(
    const char* rank_dir,
    const ST_META* meta,
    const GLOBAL_ST_ELEM* elem,
    const GLOBAL_ST_ELEM_MAP* map,
    const I64_VECTOR* local_cells,
    const int* global_to_local)
{
    char elem_path[PATH_MAX];
    char map_path[PATH_MAX];
    make_path(elem_path, sizeof(elem_path), rank_dir, "st_elem.dat");
    make_path(map_path, sizeof(map_path), rank_dir, "st_elem_map.dat");

    FILE* elem_fp = open_write(elem_path);
    FILE* map_fp = open_write(map_path);

    fprintf(
        elem_fp,
        "%zu %d\n",
        local_cells->size,
        elem->width);
    fprintf(
        map_fp,
        "%zu %d\n",
        local_cells->size,
        map->num_space_nodes);

    for(size_t local_cell = 0;
        local_cell < local_cells->size;
        local_cell++) {
        const int64_t source_cell = local_cells->data[local_cell];
        const int64_t* source_conn = &elem->conn[
            (size_t)source_cell * (size_t)elem->width
        ];

        for(int a = 0; a < elem->width; a++) {
            const int local_node = global_to_local[source_conn[a]];
            if(local_node < 0) {
                fprintf(
                    stderr,
                    "ERROR: selected ST cell %" PRId64
                    " contains a node absent from the rank-local map\n",
                    source_cell);
                exit(EXIT_FAILURE);
            }
            if(a != 0) {
                fputc(' ', elem_fp);
            }
            fprintf(elem_fp, "%d", local_node);
        }
        fputc('\n', elem_fp);

        const int64_t global_cell = map->global_cell[source_cell];
        const int slab = map->slab[source_cell];
        const int global_elem = map->global_space_elem[source_cell];
        const int64_t expected_cell =
            (int64_t)slab * meta->num_global_space_elements + global_elem;
        if(global_cell != expected_cell) {
            fail("cell metadata became inconsistent during output");
        }

        fprintf(
            map_fp,
            "%zu %" PRId64 " %d %d",
            local_cell,
            global_cell,
            slab,
            global_elem);
        for(int a = 0; a < map->num_space_nodes; a++) {
            fprintf(
                map_fp,
                " %d",
                map->global_space_conn[
                    (size_t)source_cell * (size_t)map->num_space_nodes
                    + (size_t)a
                ]);
        }
        fputc('\n', map_fp);
    }

    fclose(elem_fp);
    fclose(map_fp);
}

static void write_rank_info(
    const char* rank_dir,
    int rank,
    int num_parts,
    size_t num_internal,
    size_t num_halo,
    size_t num_cells)
{
    char path[PATH_MAX];
    make_path(path, sizeof(path), rank_dir, "partition_info.txt");
    FILE* fp = open_write(path);
    fprintf(fp, "rank %d\n", rank);
    fprintf(fp, "num_parts %d\n", num_parts);
    fprintf(fp, "num_internal_st_nodes %zu\n", num_internal);
    fprintf(fp, "num_halo_st_nodes %zu\n", num_halo);
    fprintf(fp, "num_total_st_nodes %zu\n", num_internal + num_halo);
    fprintf(fp, "num_local_st_cells %zu\n", num_cells);
    fclose(fp);
}

static void free_global_graph(GLOBAL_GRAPH* graph)
{
    free(graph->index);
    free(graph->item);
    memset(graph, 0, sizeof(*graph));
}

static void free_global_st_elem(GLOBAL_ST_ELEM* elem)
{
    free(elem->conn);
    memset(elem, 0, sizeof(*elem));
}

static void free_global_st_elem_map(GLOBAL_ST_ELEM_MAP* map)
{
    free(map->global_cell);
    free(map->slab);
    free(map->global_space_elem);
    free(map->global_space_conn);
    memset(map, 0, sizeof(*map));
}

static void print_usage(const char* program)
{
    fprintf(
        stderr,
        "Usage:\n"
        "  %s NUM_PARTS st_meta.dat graph.dat st_elem.dat "
        "st_elem_map.dat owner.dat OUTPUT_DIR\n"
        "\n"
        "Example:\n"
        "  %s 4 ./st_meta.dat ./graph.dat ./st_elem.dat "
        "./st_elem_map.dat ./owner.dat ./st_parted\n",
        program,
        program);
}

int main(int argc, char* argv[])
{
    if(argc != 8) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    const int num_parts = atoi(argv[1]);
    const char* meta_path = argv[2];
    const char* graph_path = argv[3];
    const char* elem_path = argv[4];
    const char* elem_map_path = argv[5];
    const char* owner_path = argv[6];
    const char* output_dir = argv[7];

    if(num_parts <= 0) {
        fail("NUM_PARTS must be positive");
    }

    ST_META meta;
    GLOBAL_GRAPH graph;
    GLOBAL_ST_ELEM elem;
    GLOBAL_ST_ELEM_MAP elem_map;

    fprintf(stderr, "Reading metadata...\n");
    read_meta(meta_path, &meta);
    fprintf(stderr, "Reading global graph...\n");
    read_global_graph(graph_path, &meta, &graph);
    fprintf(stderr, "Reading global ST connectivity...\n");
    read_global_st_elem(elem_path, &meta, &elem);
    fprintf(stderr, "Reading global ST cell metadata...\n");
    read_global_st_elem_map(elem_map_path, &meta, &elem, &elem_map);
    fprintf(stderr, "Reading node ownership...\n");
    int* owner = read_owner_file(owner_path, &meta, num_parts);

    make_directories_recursive(output_dir);

    const size_t global_nodes = (size_t)meta.num_global_st_nodes;
    int* global_to_local = (int*)checked_malloc(global_nodes, sizeof(int));
    int64_t* local_to_global = (int64_t*)checked_malloc(
        global_nodes, sizeof(int64_t));
    unsigned char* local_mark = (unsigned char*)checked_calloc(
        global_nodes, sizeof(unsigned char));

    char summary_path[PATH_MAX];
    make_path(
        summary_path,
        sizeof(summary_path),
        output_dir,
        "partition_summary.csv");
    FILE* summary_fp = open_write(summary_path);
    fprintf(
        summary_fp,
        "rank,num_internal_st_nodes,num_halo_st_nodes,"
        "num_total_st_nodes,num_local_st_cells\n");

    for(int rank = 0; rank < num_parts; rank++) {
        fprintf(stderr, "Building rank %d / %d...\n", rank, num_parts - 1);

        I64_VECTOR internal = {0};
        I64_VECTOR halo = {0};
        I64_VECTOR local_cells = {0};
        build_rank_sets(
            rank,
            &meta,
            &graph,
            &elem,
            owner,
            &internal,
            &halo,
            &local_cells,
            local_mark);

        if(internal.size == 0) {
            fprintf(
                stderr,
                "ERROR: rank %d owns no ST nodes; reduce NUM_PARTS or fix owner.dat.\n",
                rank);
            exit(EXIT_FAILURE);
        }

        char rank_dir[PATH_MAX];
        const int n = snprintf(
            rank_dir,
            sizeof(rank_dir),
            "%s/rank_%06d",
            output_dir,
            rank);
        if(n < 0 || (size_t)n >= sizeof(rank_dir)) {
            fail("rank output directory path is too long");
        }
        make_directories_recursive(rank_dir);

        write_st_node_map(
            rank_dir,
            &meta,
            &internal,
            &halo,
            global_to_local,
            local_to_global);

        write_local_graph(
            rank_dir,
            &graph,
            internal.size,
            internal.size + halo.size,
            local_to_global,
            global_to_local);

        write_local_cells(
            rank_dir,
            &meta,
            &elem,
            &elem_map,
            &local_cells,
            global_to_local);

        char local_meta_path[PATH_MAX];
        make_path(
            local_meta_path,
            sizeof(local_meta_path),
            rank_dir,
            "st_meta.dat");
        copy_file(meta_path, local_meta_path);

        write_rank_info(
            rank_dir,
            rank,
            num_parts,
            internal.size,
            halo.size,
            local_cells.size);

        fprintf(
            summary_fp,
            "%d,%zu,%zu,%zu,%zu\n",
            rank,
            internal.size,
            halo.size,
            internal.size + halo.size,
            local_cells.size);

        free(internal.data);
        free(halo.data);
        free(local_cells.data);
    }

    fclose(summary_fp);

    free(global_to_local);
    free(local_to_global);
    free(local_mark);
    free(owner);
    free_global_st_elem_map(&elem_map);
    free_global_st_elem(&elem);
    free_global_graph(&graph);

    printf("Generated rank-local ST files in: %s\n", output_dir);
    printf("  ranks               : %d\n", num_parts);
    printf("  global ST nodes     : %" PRId64 "\n", meta.num_global_st_nodes);
    printf("  global ST cells     : %" PRId64 "\n", meta.num_global_st_cells);
    printf("  summary             : %s\n", summary_path);
    printf("\n");
    printf("Note: this program creates the solver-side rank-local graph and maps.\n");
    printf("MONOLIS/GEDATSU communication-table files, if required by your build,\n");
    printf("must still be generated with the corresponding MONOLIS/GEDATSU tool.\n");

    return EXIT_SUCCESS;
}
