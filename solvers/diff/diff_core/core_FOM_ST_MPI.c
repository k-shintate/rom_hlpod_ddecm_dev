#include "core_FOM_ST_MPI.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
	int* data;
	int size;
	int capacity;
} ST_INT_VECTOR;

static void ST_fail(const char* message)
{
	fprintf(stderr, "ERROR: %s\n", message);
	exit(EXIT_FAILURE);
}

static void* ST_calloc(size_t n, size_t size)
{
	if(size != 0 && n > SIZE_MAX / size) {
		ST_fail("allocation size overflow");
	}
	void* p = calloc(n, size);
	if(p == NULL) {
		ST_fail("memory allocation failed");
	}
	return p;
}

static void* ST_realloc(void* ptr, size_t n, size_t size)
{
	if(size != 0 && n > SIZE_MAX / size) {
		ST_fail("reallocation size overflow");
	}
	void* p = realloc(ptr, n * size);
	if(p == NULL) {
		ST_fail("memory reallocation failed");
	}
	return p;
}

static int ST_compare_int(const void* a, const void* b)
{
	const int ia = *(const int*)a;
	const int ib = *(const int*)b;
	return (ia > ib) - (ia < ib);
}

static void ST_vector_push(ST_INT_VECTOR* v, int value)
{
	if(v->size == v->capacity) {
		const int next = (v->capacity == 0) ? 8 : 2 * v->capacity;
		if(next < v->capacity) {
			ST_fail("integer-vector capacity overflow");
		}
		v->data = (int*)ST_realloc(v->data, (size_t)next, sizeof(int));
		v->capacity = next;
	}
	v->data[v->size++] = value;
}

static void ST_vector_sort_unique_without_self(ST_INT_VECTOR* v, int self)
{
	if(v->size == 0) {
		return;
	}
	qsort(v->data, (size_t)v->size, sizeof(int), ST_compare_int);
	int write = 0;
	for(int read=0; read<v->size; read++) {
		const int value = v->data[read];
		if(value == self) {
			continue;
		}
		if(write == 0 || value != v->data[write - 1]) {
			v->data[write++] = value;
		}
	}
	v->size = write;
}

static uint64_t ST_hash_u64(uint64_t x)
{
	x += UINT64_C(0x9e3779b97f4a7c15);
	x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
	x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
	return x ^ (x >> 31);
}

void ST_ID_MAP_initialize(ST_ID_MAP* map, size_t expected_size)
{
	if(map == NULL) {
		ST_fail("ST_ID_MAP_initialize received NULL");
	}
	memset(map, 0, sizeof(*map));

	size_t capacity = 16;
	const size_t target = (expected_size < SIZE_MAX / 2)
		? 2 * expected_size
		: SIZE_MAX;
	while(capacity < target) {
		if(capacity > SIZE_MAX / 2) {
			ST_fail("ID-map capacity overflow");
		}
		capacity *= 2;
	}

	map->capacity = capacity;
	map->keys = (int64_t*)ST_calloc(capacity, sizeof(int64_t));
	map->values = (int*)ST_calloc(capacity, sizeof(int));
	map->used = (unsigned char*)ST_calloc(capacity, sizeof(unsigned char));
}

void ST_ID_MAP_insert(ST_ID_MAP* map, int64_t key, int value)
{
	if(map == NULL || map->capacity == 0) {
		ST_fail("ST_ID_MAP_insert called before initialization");
	}
	if(map->size * 10 >= map->capacity * 7) {
		ST_fail("ID map load factor exceeded; increase initial capacity");
	}

	const size_t mask = map->capacity - 1;
	size_t pos = (size_t)ST_hash_u64((uint64_t)key) & mask;
	while(map->used[pos]) {
		if(map->keys[pos] == key) {
			fprintf(stderr, "ERROR: duplicated global ID: %" PRId64 "\n", key);
			exit(EXIT_FAILURE);
		}
		pos = (pos + 1) & mask;
	}
	map->used[pos] = 1;
	map->keys[pos] = key;
	map->values[pos] = value;
	map->size++;
}

int ST_ID_MAP_find(const ST_ID_MAP* map, int64_t key)
{
	if(map == NULL || map->capacity == 0) {
		return -1;
	}
	const size_t mask = map->capacity - 1;
	size_t pos = (size_t)ST_hash_u64((uint64_t)key) & mask;
	const size_t start = pos;
	while(map->used[pos]) {
		if(map->keys[pos] == key) {
			return map->values[pos];
		}
		pos = (pos + 1) & mask;
		if(pos == start) {
			break;
		}
	}
	return -1;
}

void ST_ID_MAP_finalize(ST_ID_MAP* map)
{
	if(map == NULL) {
		return;
	}
	free(map->keys);
	free(map->values);
	free(map->used);
	memset(map, 0, sizeof(*map));
}
static FILE* ST_open_parted_file(
        const char* directory,
        const char* top_dir,
        const char* part_dir,
        const char* label)
{
    if(directory == NULL ||
       top_dir == NULL ||
       part_dir == NULL ||
       label == NULL) {
        ST_fail("NULL argument in ST_open_parted_file");
    }

    const char* filename =
        monolis_get_global_input_file_name(
                top_dir,
                part_dir,
                label);

    FILE* fp = NULL;

    fp = BBFE_sys_read_fopen(
            fp,
            filename,
            directory);

    if(fp == NULL) {
        fprintf(
                stderr,
                "ERROR: cannot open partitioned ST file \"%s\"\n",
                filename);

        exit(EXIT_FAILURE);
    }

    return fp;
}

void ST_partition_initialize(ST_PARTITION* st)
{
	if(st == NULL) {
		ST_fail("ST_partition_initialize received NULL");
	}
	memset(st, 0, sizeof(*st));
}

static void ST_read_node_map(ST_PARTITION* st, const char* directory)
{
	FILE* fp = ST_open_parted_file(
			directory,
			ST_MPI_TOP_DIR,
			ST_MPI_PART_DIR,
			"st_node_map.dat");

	if(fscanf(
			fp,
			"%d %d",
			&st->num_internal_st_nodes,
			&st->num_total_st_nodes) != 2 ||
	   st->num_internal_st_nodes < 0 ||
	   st->num_total_st_nodes <= 0 ||
	   st->num_internal_st_nodes > st->num_total_st_nodes) {
		fclose(fp);
		ST_fail("invalid st_node_map.dat header");
	}

	const int np = st->num_total_st_nodes;
	st->local_to_global_st = (int64_t*)ST_calloc((size_t)np, sizeof(int64_t));
	st->slab_id = (int*)ST_calloc((size_t)np, sizeof(int));
	st->time_dof = (int*)ST_calloc((size_t)np, sizeof(int));
	st->global_space_node_id = (int*)ST_calloc((size_t)np, sizeof(int));
	ST_ID_MAP_initialize(&st->global_to_local_st, (size_t)np);

	for(int row=0; row<np; row++) {
		int local_id = -1;
		int64_t global_id = -1;
		int slab = -1;
		int alpha = -1;
		int global_space = -1;
		if(fscanf(
				fp,
				"%d %" SCNd64 " %d %d %d",
				&local_id,
				&global_id,
				&slab,
				&alpha,
				&global_space) != 5) {
			fclose(fp);
			ST_fail("failed to read st_node_map.dat row");
		}
		if(local_id < 0 || local_id >= np ||
		   global_id < 0 || global_id >= st->meta.num_global_st_nodes ||
		   slab < 0 || slab >= st->meta.num_window_slabs ||
		   alpha < 0 || alpha >= st->meta.num_time_dofs ||
		   global_space < 0 || global_space >= st->meta.num_global_space_nodes) {
			fclose(fp);
			ST_fail("out-of-range value in st_node_map.dat");
		}
		const int64_t expected = ST_global_gid_partitioned(
				slab,
				alpha,
				global_space,
				st->meta.num_time_dofs,
				st->meta.num_global_space_nodes);
		if(global_id != expected) {
			fprintf(
					stderr,
					"ERROR: ST global ID mismatch at local node %d: "
					"file=%" PRId64 " expected=%" PRId64 "\n",
					local_id,
					global_id,
					expected);
			exit(EXIT_FAILURE);
		}
		st->local_to_global_st[local_id] = global_id;
		st->slab_id[local_id] = slab;
		st->time_dof[local_id] = alpha;
		st->global_space_node_id[local_id] = global_space;
		ST_ID_MAP_insert(&st->global_to_local_st, global_id, local_id);
	}
	fclose(fp);
}

static void ST_read_graph(ST_PARTITION* st, const char* directory)
{
	FILE* fp = ST_open_parted_file(
			directory,
			ST_MPI_TOP_DIR,
			ST_MPI_PART_DIR,
			"graph.dat");

	int num_nodes = 0;
	if(fscanf(fp, "%d", &num_nodes) != 1 ||
	   num_nodes != st->num_total_st_nodes) {
		fclose(fp);
		ST_fail("graph.dat node count is inconsistent with st_node_map.dat");
	}

	ST_INT_VECTOR* rows = (ST_INT_VECTOR*)ST_calloc(
			(size_t)num_nodes,
			sizeof(ST_INT_VECTOR));
	unsigned char* row_seen = (unsigned char*)ST_calloc(
			(size_t)num_nodes,
			sizeof(unsigned char));

	for(int r=0; r<num_nodes; r++) {
		int row_id = -1;
		int num_adj = -1;
		if(fscanf(fp, "%d %d", &row_id, &num_adj) != 2 ||
		   row_id < 0 || row_id >= num_nodes ||
		   num_adj < 0 || row_seen[row_id]) {
			fclose(fp);
			ST_fail("invalid graph.dat row header");
		}
		row_seen[row_id] = 1;
		for(int j=0; j<num_adj; j++) {
			int col = -1;
			if(fscanf(fp, "%d", &col) != 1 || col < 0 || col >= num_nodes) {
				fclose(fp);
				ST_fail("invalid graph.dat column");
			}
			ST_vector_push(&rows[row_id], col);
		}
	}
	fclose(fp);

	int64_t total_items_64 = 0;
	for(int i=0; i<num_nodes; i++) {
		if(!row_seen[i]) {
			ST_fail("missing row in graph.dat");
		}
		ST_vector_sort_unique_without_self(&rows[i], i);
		total_items_64 += rows[i].size;
	}
	if(total_items_64 > INT_MAX) {
		ST_fail("local ST graph exceeds INT_MAX");
	}

	st->graph_index = (int*)ST_calloc((size_t)num_nodes + 1, sizeof(int));
	st->graph_item = (int*)ST_calloc((size_t)total_items_64, sizeof(int));
	int position = 0;
	for(int i=0; i<num_nodes; i++) {
		st->graph_index[i] = position;
		for(int j=0; j<rows[i].size; j++) {
			st->graph_item[position++] = rows[i].data[j];
		}
		free(rows[i].data);
	}
	st->graph_index[num_nodes] = position;
	free(rows);
	free(row_seen);
}

static void ST_read_elem(ST_PARTITION* st, const char* directory)
{
	FILE* fp = ST_open_parted_file(
			directory,
			ST_MPI_TOP_DIR,
			ST_MPI_PART_DIR,
			"st_elem.dat");

	if(fscanf(
			fp,
			"%d %d",
			&st->num_local_st_cells,
			&st->num_st_nodes_per_cell) != 2 ||
	st->num_local_st_cells < 0 ||
	st->num_st_nodes_per_cell <= 0) {
		fclose(fp);
		ST_fail("invalid partitioned st_elem.dat header");
	}

	const int expected_nodes_per_cell =
		st->meta.num_time_dofs *
		st->meta.num_space_nodes_per_element;

	if(st->num_st_nodes_per_cell != expected_nodes_per_cell) {
		fclose(fp);
		ST_fail(
				"st_elem.dat width does not equal "
				"nt * spatial element width");
	}

	const size_t total =
		(size_t)st->num_local_st_cells *
		(size_t)st->num_st_nodes_per_cell;

	st->cell_st_conn =
		(int*)ST_calloc(total, sizeof(int));

	for(size_t k=0; k<total; k++) {
		if(fscanf(fp, "%d", &st->cell_st_conn[k]) != 1 ||
		st->cell_st_conn[k] < 0 ||
		st->cell_st_conn[k] >= st->num_total_st_nodes) {
			fclose(fp);
			ST_fail("invalid local ST node ID in st_elem.dat");
		}
	}

	fclose(fp);
}

static void ST_read_cell_map(ST_PARTITION* st, const char* directory)
{
	FILE* fp = ST_open_parted_file(
			directory,
			ST_MPI_TOP_DIR,
			ST_MPI_PART_DIR,
			"st_elem_map.dat");

	int num_cells = -1;
	int num_space_nodes = -1;
	if(fscanf(fp, "%d %d", &num_cells, &num_space_nodes) != 2 ||
	   num_cells != st->num_local_st_cells ||
	   num_space_nodes != st->meta.num_space_nodes_per_element) {
		fclose(fp);
		ST_fail("st_elem_map.dat header is inconsistent");
	}

	st->global_st_cell_id = (int64_t*)ST_calloc(
			(size_t)num_cells,
			sizeof(int64_t));
	st->cell_slab_id = (int*)ST_calloc((size_t)num_cells, sizeof(int));
	st->cell_global_space_elem_id = (int*)ST_calloc(
			(size_t)num_cells,
			sizeof(int));
	st->cell_global_space_conn = (int*)ST_calloc(
			(size_t)num_cells * (size_t)num_space_nodes,
			sizeof(int));

	ST_ID_MAP cell_seen;
	ST_ID_MAP_initialize(&cell_seen, (size_t)num_cells);

	for(int row=0; row<num_cells; row++) {
		int local_cell = -1;
		int64_t global_cell = -1;
		int slab = -1;
		int global_elem = -1;
		if(fscanf(
				fp,
				"%d %" SCNd64 " %d %d",
				&local_cell,
				&global_cell,
				&slab,
				&global_elem) != 4 ||
		   local_cell < 0 || local_cell >= num_cells ||
		   global_cell < 0 || global_cell >= st->meta.num_global_st_cells ||
		   slab < 0 || slab >= st->meta.num_window_slabs ||
		   global_elem < 0 || global_elem >= st->meta.num_global_space_elements) {
			fclose(fp);
			ST_fail("invalid st_elem_map.dat row");
		}
		if(ST_ID_MAP_find(&cell_seen, global_cell) >= 0) {
			fclose(fp);
			ST_fail("duplicated ST cell in one rank partition");
		}
		ST_ID_MAP_insert(&cell_seen, global_cell, local_cell);

		const int64_t expected_cell =
			(int64_t)slab * st->meta.num_global_space_elements + global_elem;
		if(global_cell != expected_cell) {
			fclose(fp);
			ST_fail("global ST cell ID does not match slab and spatial element");
		}

		st->global_st_cell_id[local_cell] = global_cell;
		st->cell_slab_id[local_cell] = slab;
		st->cell_global_space_elem_id[local_cell] = global_elem;
		for(int a=0; a<num_space_nodes; a++) {
			int global_space_node = -1;
			if(fscanf(fp, "%d", &global_space_node) != 1 ||
			   global_space_node < 0 ||
			   global_space_node >= st->meta.num_global_space_nodes) {
				fclose(fp);
				ST_fail("invalid spatial node in st_elem_map.dat");
			}
			st->cell_global_space_conn[
				local_cell * num_space_nodes + a
			] = global_space_node;
		}
	}
	fclose(fp);
	ST_ID_MAP_finalize(&cell_seen);
}

static void ST_validate_elem_and_map(ST_PARTITION* st)
{
	const int nt = st->meta.num_time_dofs;
	const int nl = st->meta.num_space_nodes_per_element;
	for(int c=0; c<st->num_local_st_cells; c++) {
		const int slab = st->cell_slab_id[c];
		for(int alpha=0; alpha<nt; alpha++) {
			for(int a=0; a<nl; a++) {
				const int lid = st->cell_st_conn[
					c * st->num_st_nodes_per_cell + alpha * nl + a
				];
				const int global_space = st->cell_global_space_conn[c * nl + a];
				const int64_t expected_gid = ST_global_gid_partitioned(
						slab,
						alpha,
						global_space,
						nt,
						st->meta.num_global_space_nodes);
				if(st->local_to_global_st[lid] != expected_gid) {
					fprintf(
							stderr,
							"ERROR: elem/map mismatch: cell=%d alpha=%d node=%d "
							"local ST=%d global ST=%" PRId64 " expected=%" PRId64 "\n",
							c,
							alpha,
							a,
							lid,
							st->local_to_global_st[lid],
							expected_gid);
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

static void ST_read_global_meta(
        ST_GLOBAL_META* meta,
        const char* directory)
{
    if(meta == NULL || directory == NULL) {
        ST_fail("NULL argument in ST_read_global_meta");
    }

    char path[PATH_MAX];

    const int n = snprintf(
            path,
            sizeof(path),
            "%s/st_meta.dat",
            directory);

    if(n < 0 || (size_t)n >= sizeof(path)) {
        ST_fail("st_meta.dat path is too long");
    }

    FILE* fp = fopen(path, "r");

    if(fp == NULL) {
        fprintf(
                stderr,
                "ERROR: cannot open %s: %s\n",
                path,
                strerror(errno));

        exit(EXIT_FAILURE);
    }

    memset(meta, 0, sizeof(*meta));

    char key[128];
    char value[512];

    while(fscanf(
            fp,
            "%127s %511s",
            key,
            value) == 2) {

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
        else if(strcmp(key, "num_global_st_nodes") == 0) {
            meta->num_global_st_nodes =
                strtoll(value, NULL, 10);
        }
        else if(strcmp(key, "num_global_st_cells") == 0) {
            meta->num_global_st_cells =
                strtoll(value, NULL, 10);
        }
        /*
         * 現在のソルバーで使用しないメタデータ項目は
         * 読み飛ばす。
         */
    }

    if(ferror(fp)) {
        fprintf(
                stderr,
                "ERROR: failed while reading %s\n",
                path);

        fclose(fp);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    if(meta->degree < 0) {
        fprintf(
                stderr,
                "ERROR: invalid degree in %s: %d\n",
                path,
                meta->degree);

        exit(EXIT_FAILURE);
    }

    if(meta->num_time_dofs <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_time_dofs in %s: %d\n",
                path,
                meta->num_time_dofs);

        exit(EXIT_FAILURE);
    }

    if(meta->num_window_slabs <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_window_slabs in %s: %d\n",
                path,
                meta->num_window_slabs);

        exit(EXIT_FAILURE);
    }

    if(meta->num_global_space_nodes <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_global_space_nodes in %s: %d\n",
                path,
                meta->num_global_space_nodes);

        exit(EXIT_FAILURE);
    }

    if(meta->num_global_space_elements <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_global_space_elements in %s: %d\n",
                path,
                meta->num_global_space_elements);

        exit(EXIT_FAILURE);
    }

    if(meta->num_space_nodes_per_element <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_space_nodes_per_element in %s: %d\n",
                path,
                meta->num_space_nodes_per_element);

        exit(EXIT_FAILURE);
    }

    if(meta->num_global_st_nodes <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_global_st_nodes in %s: %" PRId64 "\n",
                path,
                meta->num_global_st_nodes);

        exit(EXIT_FAILURE);
    }

    if(meta->num_global_st_cells <= 0) {
        fprintf(
                stderr,
                "ERROR: invalid num_global_st_cells in %s: %" PRId64 "\n",
                path,
                meta->num_global_st_cells);

        exit(EXIT_FAILURE);
    }

    /*
     * グローバルST節点数を検証する。
     *
     * global ST nodes
     *   = num_window_slabs
     *   × num_time_dofs
     *   × num_global_space_nodes
     */
    const int64_t expected_num_global_st_nodes =
        (int64_t)meta->num_window_slabs
        * (int64_t)meta->num_time_dofs
        * (int64_t)meta->num_global_space_nodes;

    if(meta->num_global_st_nodes !=
       expected_num_global_st_nodes) {

        fprintf(
                stderr,
                "ERROR: inconsistent num_global_st_nodes in %s\n"
                "  file value     = %" PRId64 "\n"
                "  expected value = %" PRId64 "\n",
                path,
                meta->num_global_st_nodes,
                expected_num_global_st_nodes);

        exit(EXIT_FAILURE);
    }

    /*
     * グローバルSTセル数を検証する。
     *
     * global ST cells
     *   = num_window_slabs
     *   × num_global_space_elements
     */
    const int64_t expected_num_global_st_cells =
        (int64_t)meta->num_window_slabs
        * (int64_t)meta->num_global_space_elements;

    if(meta->num_global_st_cells !=
       expected_num_global_st_cells) {

        fprintf(
                stderr,
                "ERROR: inconsistent num_global_st_cells in %s\n"
                "  file value     = %" PRId64 "\n"
                "  expected value = %" PRId64 "\n",
                path,
                meta->num_global_st_cells,
                expected_num_global_st_cells);

        exit(EXIT_FAILURE);
    }
}

void ST_partition_read(
		ST_PARTITION* st,
		const char* directory,
		const ST_TIME_ELEMENT* te,
		int num_window_slabs,
		const MONOLIS_COM* parent_com)
{
	if(st == NULL || directory == NULL || te == NULL ||
	   parent_com == NULL) {
		ST_fail("NULL argument in ST_partition_read");
	}
	ST_read_global_meta(&st->meta, directory);
	if(st->meta.degree != te->degree ||
	   st->meta.num_time_dofs != te->n_dof ||
	   st->meta.num_window_slabs != num_window_slabs) {
		fprintf(
				stderr,
				"ERROR: ST preprocessing/runtime mismatch.\n"
				"  meta degree/slabs/nt = %d/%d/%d\n"
				"  runtime degree/slabs/nt = %d/%d/%d\n",
				st->meta.degree,
				st->meta.num_window_slabs,
				st->meta.num_time_dofs,
				te->degree,
				num_window_slabs,
				te->n_dof);
		exit(EXIT_FAILURE);
	}

	ST_read_node_map(st, directory);
	ST_read_graph(st, directory);
	ST_read_elem(st, directory);
	ST_read_cell_map(st, directory);
	ST_validate_elem_and_map(st);

	memset(&st->monolis_com_st, 0, sizeof(st->monolis_com_st));
	monolis_com_initialize_by_parted_files(
			&st->monolis_com_st,
			monolis_mpi_get_global_comm(),
			ST_MPI_TOP_DIR,
			ST_MPI_PART_DIR,
			"graph.dat");
}

void ST_partition_finalize(ST_PARTITION* st)
{
	if(st == NULL) {
		return;
	}
	free(st->local_to_global_st);
	free(st->slab_id);
	free(st->time_dof);
	free(st->global_space_node_id);
	ST_ID_MAP_finalize(&st->global_to_local_st);
	free(st->graph_index);
	free(st->graph_item);
	free(st->cell_st_conn);
	free(st->global_st_cell_id);
	free(st->cell_slab_id);
	free(st->cell_global_space_elem_id);
	free(st->cell_global_space_conn);
	memset(st, 0, sizeof(*st));
}

void ST_space_map_initialize(ST_SPACE_MAP* map)
{
	if(map == NULL) {
		ST_fail("ST_space_map_initialize received NULL");
	}
	memset(map, 0, sizeof(*map));
}


void ST_space_map_build_identity_from_mesh(
		ST_SPACE_MAP* map,
		const BBFE_DATA* fe)
{
	if(map == NULL || fe == NULL) {
		ST_fail(
				"NULL argument in "
				"ST_space_map_build_identity_from_mesh");
	}

	if(fe->total_num_nodes <= 0 ||
	   fe->total_num_elems <= 0) {
		ST_fail(
				"invalid BBFE mesh dimensions in "
				"ST_space_map_build_identity_from_mesh");
	}

	if(map->local_to_global_space_node != NULL ||
	   map->local_to_global_space_elem != NULL ||
	   map->global_to_local_space_node.capacity != 0 ||
	   map->global_to_local_space_elem.capacity != 0) {
		ST_fail(
				"ST_space_map_build_identity_from_mesh "
				"received an already initialized map");
	}

	/*
	 * Current execution model:
	 *
	 *   - every MPI rank reads the complete spatial node.dat;
	 *   - every MPI rank reads the complete spatial elem.dat;
	 *   - local spatial IDs are therefore identical to global spatial IDs.
	 *
	 * The ST graph and ST cells are distributed independently.
	 */
	map->num_internal_space_nodes = fe->total_num_nodes;
	map->num_total_space_nodes = fe->total_num_nodes;
	map->num_local_space_elements = fe->total_num_elems;

	map->local_to_global_space_node =
		(int64_t*)ST_calloc(
				(size_t)map->num_total_space_nodes,
				sizeof(int64_t));

	map->local_to_global_space_elem =
		(int64_t*)ST_calloc(
				(size_t)map->num_local_space_elements,
				sizeof(int64_t));

	ST_ID_MAP_initialize(
			&map->global_to_local_space_node,
			(size_t)map->num_total_space_nodes);

	ST_ID_MAP_initialize(
			&map->global_to_local_space_elem,
			(size_t)map->num_local_space_elements);

	for(int i=0; i<map->num_total_space_nodes; i++) {
		map->local_to_global_space_node[i] = (int64_t)i;
		ST_ID_MAP_insert(
				&map->global_to_local_space_node,
				(int64_t)i,
				i);
	}

	for(int e=0; e<map->num_local_space_elements; e++) {
		map->local_to_global_space_elem[e] = (int64_t)e;
		ST_ID_MAP_insert(
				&map->global_to_local_space_elem,
				(int64_t)e,
				e);
	}
}

void ST_space_map_read(
		ST_SPACE_MAP* map,
		const char* directory,
		const char* top_dir,
		const char* part_dir)
{
	if(map == NULL || directory == NULL || top_dir == NULL || part_dir == NULL) {
		ST_fail("NULL argument in ST_space_map_read");
	}

	FILE* fp = ST_open_parted_file(
			directory,
			top_dir,
			part_dir,
			"space_node_map.dat");
	if(fscanf(
			fp,
			"%d %d",
			&map->num_internal_space_nodes,
			&map->num_total_space_nodes) != 2 ||
	   map->num_internal_space_nodes < 0 ||
	   map->num_total_space_nodes <= 0 ||
	   map->num_internal_space_nodes > map->num_total_space_nodes) {
		fclose(fp);
		ST_fail("invalid space_node_map.dat header");
	}
	map->local_to_global_space_node = (int64_t*)ST_calloc(
			(size_t)map->num_total_space_nodes,
			sizeof(int64_t));
	ST_ID_MAP_initialize(
			&map->global_to_local_space_node,
			(size_t)map->num_total_space_nodes);
	for(int row=0; row<map->num_total_space_nodes; row++) {
		int local_id = -1;
		int64_t global_id = -1;
		if(fscanf(fp, "%d %" SCNd64, &local_id, &global_id) != 2 ||
		   local_id < 0 || local_id >= map->num_total_space_nodes ||
		   global_id < 0) {
			fclose(fp);
			ST_fail("invalid space_node_map.dat row");
		}
		map->local_to_global_space_node[local_id] = global_id;
		ST_ID_MAP_insert(&map->global_to_local_space_node, global_id, local_id);
	}
	fclose(fp);

	fp = ST_open_parted_file(
			directory,
			top_dir,
			part_dir,
			"space_elem_map.dat");
	if(fscanf(fp, "%d", &map->num_local_space_elements) != 1 ||
	   map->num_local_space_elements < 0) {
		fclose(fp);
		ST_fail("invalid space_elem_map.dat header");
	}
	map->local_to_global_space_elem = (int64_t*)ST_calloc(
			(size_t)map->num_local_space_elements,
			sizeof(int64_t));
	ST_ID_MAP_initialize(
			&map->global_to_local_space_elem,
			(size_t)map->num_local_space_elements);
	for(int row=0; row<map->num_local_space_elements; row++) {
		int local_id = -1;
		int64_t global_id = -1;
		if(fscanf(fp, "%d %" SCNd64, &local_id, &global_id) != 2 ||
		   local_id < 0 || local_id >= map->num_local_space_elements ||
		   global_id < 0) {
			fclose(fp);
			ST_fail("invalid space_elem_map.dat row");
		}
		map->local_to_global_space_elem[local_id] = global_id;
		ST_ID_MAP_insert(&map->global_to_local_space_elem, global_id, local_id);
	}
	fclose(fp);
}

void ST_space_map_finalize(ST_SPACE_MAP* map)
{
	if(map == NULL) {
		return;
	}
	free(map->local_to_global_space_node);
	free(map->local_to_global_space_elem);
	ST_ID_MAP_finalize(&map->global_to_local_space_node);
	ST_ID_MAP_finalize(&map->global_to_local_space_elem);
	memset(map, 0, sizeof(*map));
}

int64_t ST_global_gid_partitioned(
		int slab,
		int alpha,
		int global_space_node,
		int num_time_dofs,
		int num_global_space_nodes)
{
	return
		((int64_t)slab * num_time_dofs + alpha)
		* (int64_t)num_global_space_nodes
		+ global_space_node;
}

void ST_set_nonzero_pattern_from_partition(
		MONOLIS* monolis,
		const ST_PARTITION* st)
{
	if(monolis == NULL || st == NULL ||
	   st->graph_index == NULL || st->graph_item == NULL) {
		ST_fail("invalid data in ST_set_nonzero_pattern_from_partition");
	}

	/*
	 * The MONOLIS nodal-graph constructor allocates a square NP x NP local
	 * pattern and initially sets N = NP.  We pass all internal+halo nodes,
	 * then restore N to the number of owned rows.
	 */
	monolis_get_nonzero_pattern_by_nodal_graph_R(
			monolis,
			st->num_total_st_nodes,
			1,
			st->graph_index,
			st->graph_item);
	monolis->mat.N = st->num_internal_st_nodes;
	monolis->mat.NP = st->num_total_st_nodes;
}

static void ST_compute_element_matrices(
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		int e,
		double* Me,
		double* Ke)
{
	const int nl = fe->local_num_nodes;
	const int np = basis->num_integ_points;

	double* val_ip = BB_std_calloc_1d_double(val_ip, np);
	double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
	double** local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	double** x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	double* a_ip = BB_std_calloc_1d_double(a_ip, np);
	double* k_ip = BB_std_calloc_1d_double(k_ip, np);

	BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
	BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
	for(int p=0; p<np; p++) {
		BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
		a_ip[p] = manusol_get_mass_coef(x_ip[p]);
		k_ip[p] = manusol_get_diff_coef(x_ip[p]);
	}

	for(int i=0; i<nl; i++) {
		for(int j=0; j<nl; j++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = BBFE_elemmat_convdiff_mat_mass(
						basis->N[p][i],
						basis->N[p][j],
						a_ip[p]);
			}
			Me[i * nl + j] = BBFE_std_integ_calc(
					np,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);

			for(int p=0; p<np; p++) {
				val_ip[p] = -BBFE_elemmat_convdiff_mat_diff(
						fe->geo[e][p].grad_N[i],
						fe->geo[e][p].grad_N[j],
						k_ip[p]);
			}
			Ke[i * nl + j] = BBFE_std_integ_calc(
					np,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);
		}
	}

	BB_std_free_1d_double(val_ip, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}

static void ST_compute_element_source_vectors(
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		const ST_TIME_ELEMENT* te,
		int e,
		double t_old,
		double dt,
		double* Fe)
{
	const int nl = fe->local_num_nodes;
	const int np = basis->num_integ_points;
	const int nt = te->n_dof;

	double* val_ip = BB_std_calloc_1d_double(val_ip, np);
	double* Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
	double** local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	double** x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	double** v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	double* a_ip = BB_std_calloc_1d_double(a_ip, np);
	double* k_ip = BB_std_calloc_1d_double(k_ip, np);

	BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);
	BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
	for(int p=0; p<np; p++) {
		BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
		manusol_get_conv_vel(v_ip[p], x_ip[p]);
		a_ip[p] = manusol_get_mass_coef(x_ip[p]);
		k_ip[p] = manusol_get_diff_coef(x_ip[p]);
	}

	for(int alpha=0; alpha<nt; alpha++) {
		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;
				for(int q=0; q<te->num_time_ip; q++) {
					const double tq = t_old + te->tau[q] * dt;
					const double coef = te->weight[q] * te->psi[alpha][q];
					if(fabs(coef) < 1.0e-30) {
						continue;
					}
					const double f = manusol_get_source(
							x_ip[p],
							tq,
							a_ip[p],
							v_ip[p],
							k_ip[p]);
					val_ip[p] += coef * BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i],
							f);
				}
			}
			Fe[alpha * nl + i] = BBFE_std_integ_calc(
					np,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);
		}
	}

	BB_std_free_1d_double(val_ip, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_x, nl, 3);
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
}

void ST_validate_partition_against_space_mesh(
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		const BBFE_DATA* fe)
{
	if(st == NULL || space_map == NULL || fe == NULL) {
		ST_fail(
				"NULL argument in "
				"ST_validate_partition_against_space_mesh");
	}

	if(fe->total_num_nodes !=
	   space_map->num_total_space_nodes) {
		fprintf(
				stderr,
				"ERROR: spatial node map and BBFE mesh differ.\n"
				"  map total nodes = %d\n"
				"  BBFE total nodes = %d\n",
				space_map->num_total_space_nodes,
				fe->total_num_nodes);
		exit(EXIT_FAILURE);
	}

	if(fe->total_num_elems !=
	   space_map->num_local_space_elements) {
		fprintf(
				stderr,
				"ERROR: spatial element map and BBFE mesh differ.\n"
				"  map spatial elements = %d\n"
				"  BBFE spatial elements = %d\n"
				"  local ST cells = %d\n",
				space_map->num_local_space_elements,
				fe->total_num_elems,
				st->num_local_st_cells);
		exit(EXIT_FAILURE);
	}

	if(fe->local_num_nodes !=
	   st->meta.num_space_nodes_per_element) {
		fprintf(
				stderr,
				"ERROR: spatial element widths differ.\n"
				"  BBFE nodes/element = %d\n"
				"  ST metadata nodes/element = %d\n",
				fe->local_num_nodes,
				st->meta.num_space_nodes_per_element);
		exit(EXIT_FAILURE);
	}

	if(st->meta.num_global_space_nodes !=
	   fe->total_num_nodes ||
	   st->meta.num_global_space_elements !=
	   fe->total_num_elems) {
		fprintf(
				stderr,
				"ERROR: ST global metadata and replicated "
				"BBFE spatial mesh differ.\n"
				"  ST global nodes/elements = %d/%d\n"
				"  BBFE nodes/elements = %d/%d\n",
				st->meta.num_global_space_nodes,
				st->meta.num_global_space_elements,
				fe->total_num_nodes,
				fe->total_num_elems);
		exit(EXIT_FAILURE);
	}

	const int nl = fe->local_num_nodes;

	for(int c=0; c<st->num_local_st_cells; c++) {
		const int global_e =
			st->cell_global_space_elem_id[c];

		const int e =
			ST_ID_MAP_find(
					&space_map->global_to_local_space_elem,
					global_e);

		if(e < 0 || e >= fe->total_num_elems) {
			fprintf(
					stderr,
					"ERROR: spatial element %d required by "
					"ST cell %" PRId64
					" is absent from the BBFE mesh.\n",
					global_e,
					st->global_st_cell_id[c]);
			exit(EXIT_FAILURE);
		}

		for(int a=0; a<nl; a++) {
			const int local_space_node =
				fe->conn[e][a];

			if(local_space_node < 0 ||
			   local_space_node >=
			   space_map->num_total_space_nodes) {
				ST_fail(
						"invalid local spatial connectivity");
			}

			const int64_t global_from_mesh =
				space_map
				->local_to_global_space_node[
					local_space_node
				];

			const int global_from_cell =
				st->cell_global_space_conn[
					c * nl + a
				];

			if(global_from_mesh != global_from_cell) {
				fprintf(
						stderr,
						"ERROR: spatial connectivity mismatch.\n"
						"  ST cell             = %" PRId64 "\n"
						"  spatial element     = %d\n"
						"  local node position = %d\n"
						"  BBFE global node    = %" PRId64 "\n"
						"  ST-map global node  = %d\n",
						st->global_st_cell_id[c],
						e,
						a,
						global_from_mesh,
						global_from_cell);
				exit(EXIT_FAILURE);
			}
		}
	}
}

void set_element_mat_ST_dG_window_partitioned(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map)
{
	if(monolis == NULL || fe == NULL || basis == NULL || vals == NULL ||
	   te == NULL || st == NULL || space_map == NULL) {
		ST_fail("NULL argument in set_element_mat_ST_dG_window_partitioned");
	}
	const int nt = te->n_dof;
	const int nl = fe->local_num_nodes;
	double* Me = (double*)ST_calloc((size_t)nl * nl, sizeof(double));
	double* Ke = (double*)ST_calloc((size_t)nl * nl, sizeof(double));

	for(int c=0; c<st->num_local_st_cells; c++) {
		const int slab = st->cell_slab_id[c];
		const int global_e = st->cell_global_space_elem_id[c];
		const int e = ST_ID_MAP_find(&space_map->global_to_local_space_elem, global_e);
		if(e < 0) {
			ST_fail("required spatial element is absent during matrix assembly");
		}
		ST_compute_element_matrices(fe, basis, e, Me, Ke);

		for(int alpha=0; alpha<nt; alpha++) {
			for(int i=0; i<nl; i++) {
				const int row = st->cell_st_conn[
					c * st->num_st_nodes_per_cell + alpha * nl + i
				];
				if(row >= st->num_internal_st_nodes) {
					continue;
				}

				for(int beta=0; beta<nt; beta++) {
					const double coef_mass =
						(te->Dt[alpha][beta] + te->E_self[alpha][beta])
						/ vals->dt;
					const double coef_diff = te->Mt[alpha][beta];
					for(int j=0; j<nl; j++) {
						const int col = st->cell_st_conn[
							c * st->num_st_nodes_per_cell + beta * nl + j
						];
						const double value =
							coef_mass * Me[i * nl + j]
							+ coef_diff * Ke[i * nl + j];
						monolis_add_scalar_to_sparse_matrix_R(
								monolis, row, col, 0, 0, value);
					}
				}

				if(slab > 0) {
					for(int beta=0; beta<nt; beta++) {
						const double coef = -te->E_prev[alpha][beta] / vals->dt;
						if(fabs(coef) < 1.0e-30) {
							continue;
						}
						for(int j=0; j<nl; j++) {
							const int global_space_j =
								st->cell_global_space_conn[c * nl + j];
							const int64_t global_col = ST_global_gid_partitioned(
									slab - 1,
									beta,
									global_space_j,
									nt,
									st->meta.num_global_space_nodes);
							const int col = ST_ID_MAP_find(
									&st->global_to_local_st,
									global_col);
							if(col < 0) {
								fprintf(
										stderr,
										"ERROR: missing previous-slab halo ST node: global=%" PRId64 "\n",
										global_col);
								exit(EXIT_FAILURE);
							}
							monolis_add_scalar_to_sparse_matrix_R(
									monolis,
									row,
									col,
									0,
									0,
									coef * Me[i * nl + j]);
						}
					}
				}
			}
		}
	}
	free(Me);
	free(Ke);
}

void set_element_vec_ST_dG_window_partitioned(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BASIS* basis,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		const double* T_prev_trace_global,
		double t_window_start)
{
	if(monolis == NULL || fe == NULL || basis == NULL || vals == NULL ||
	   te == NULL || st == NULL || space_map == NULL ||
	   T_prev_trace_global == NULL) {
		ST_fail("NULL argument in set_element_vec_ST_dG_window_partitioned");
	}
	const int nt = te->n_dof;
	const int nl = fe->local_num_nodes;
	double* Me = (double*)ST_calloc((size_t)nl * nl, sizeof(double));
	double* Ke_dummy = (double*)ST_calloc((size_t)nl * nl, sizeof(double));
	double* Fe = (double*)ST_calloc((size_t)nt * nl, sizeof(double));

	for(int c=0; c<st->num_local_st_cells; c++) {
		const int slab = st->cell_slab_id[c];
		const int global_e = st->cell_global_space_elem_id[c];
		const int e = ST_ID_MAP_find(&space_map->global_to_local_space_elem, global_e);
		if(e < 0) {
			ST_fail("required spatial element is absent during RHS assembly");
		}
		ST_compute_element_matrices(fe, basis, e, Me, Ke_dummy);
		memset(Fe, 0, (size_t)nt * nl * sizeof(double));
		ST_compute_element_source_vectors(
				fe,
				basis,
				te,
				e,
				t_window_start + slab * vals->dt,
				vals->dt,
				Fe);

		for(int alpha=0; alpha<nt; alpha++) {
			for(int i=0; i<nl; i++) {
				const int row = st->cell_st_conn[
					c * st->num_st_nodes_per_cell + alpha * nl + i
				];
				if(row >= st->num_internal_st_nodes) {
					continue;
				}
				monolis->mat.R.B[row] += Fe[alpha * nl + i];

				if(slab == 0 && fabs(te->e_minus[alpha]) > 1.0e-30) {
					double inflow = 0.0;
					for(int j=0; j<nl; j++) {
						const int global_space_j =
							st->cell_global_space_conn[c * nl + j];
						inflow += Me[i * nl + j]
							* T_prev_trace_global[global_space_j];
					}
					monolis->mat.R.B[row] +=
						te->e_minus[alpha] * inflow / vals->dt;
				}
			}
		}
	}
	free(Me);
	free(Ke_dummy);
	free(Fe);
}

static double ST_partitioned_Dirichlet_value(
		BBFE_DATA* fe,
		BBFE_BC* bc,
		const ST_SPACE_MAP* space_map,
		int global_space_node,
		double time,
		int* exists)
{
	const int local_space = ST_ID_MAP_find(
			&space_map->global_to_local_space_node,
			global_space_node);
	if(local_space < 0 || local_space >= fe->total_num_nodes) {
		fprintf(
				stderr,
				"ERROR: missing spatial node %d required by an internal ST node\n",
				global_space_node);
		exit(EXIT_FAILURE);
	}
	*exists = bc->D_bc_exists[local_space] ? 1 : 0;
	if(!*exists) {
		return 0.0;
	}
	return manusol_get_sol(
			fe->x[local_space][0],
			fe->x[local_space][1],
			fe->x[local_space][2],
			time);
}

void BBFE_sys_monowrap_set_Dirichlet_bc_ST_partitioned_scalar(
		MONOLIS* monolis,
		BBFE_DATA* fe,
		BBFE_BC* bc,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		double* rhs)
{
	if(monolis == NULL || fe == NULL || bc == NULL || vals == NULL ||
	   te == NULL || st == NULL || space_map == NULL || rhs == NULL) {
		ST_fail("NULL argument in partitioned Dirichlet application");
	}
	for(int lid=0; lid<st->num_internal_st_nodes; lid++) {
		const double time = t_window_start
			+ ((double)st->slab_id[lid] + te->dof_tau[st->time_dof[lid]])
			* vals->dt;
		int exists = 0;
		const double value = ST_partitioned_Dirichlet_value(
				fe,
				bc,
				space_map,
				st->global_space_node_id[lid],
				time,
				&exists);
		if(exists) {
			monolis_set_Dirichlet_bc_R(monolis, rhs, lid, 0, value);
		}
	}
}

void ST_restore_Dirichlet_values_partitioned_scalar(
		BBFE_DATA* fe,
		BBFE_BC* bc,
		VALUES* vals,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		double* T_window)
{
	if(fe == NULL || bc == NULL || vals == NULL || te == NULL ||
	   st == NULL || space_map == NULL || T_window == NULL) {
		ST_fail("NULL argument in partitioned Dirichlet restoration");
	}
	for(int lid=0; lid<st->num_internal_st_nodes; lid++) {
		const double time = t_window_start
			+ ((double)st->slab_id[lid] + te->dof_tau[st->time_dof[lid]])
			* vals->dt;
		int exists = 0;
		const double value = ST_partitioned_Dirichlet_value(
				fe,
				bc,
				space_map,
				st->global_space_node_id[lid],
				time,
				&exists);
		if(exists) {
			T_window[lid] = value;
		}
	}
}

void solver_fom_space_time_partitioned(
		FE_SYSTEM* sys,
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		const ST_SPACE_MAP* space_map,
		double t_window_start,
		int window_step,
		const double* T_prev_trace_global,
		double* T_window)
{
	if(sys == NULL || te == NULL || st == NULL || space_map == NULL ||
	   T_prev_trace_global == NULL || T_window == NULL) {
		ST_fail("NULL argument in solver_fom_space_time_partitioned");
	}
	const double t_window_end = t_window_start
		+ st->meta.num_window_slabs * sys->vals.dt;
	if(monolis_mpi_get_global_my_rank() == 0) {
		printf(
				"\n%s ----------------- ST window %d : [%e, %e] ----------------\n",
				CODENAME,
				window_step,
				t_window_start,
				t_window_end);
	}

	monolis_clear_mat_value_rhs_R(&sys->monolis);
	monolis_copy_mat_value_R(&sys->monolis0, &sys->monolis);

	set_element_vec_ST_dG_window_partitioned(
			&sys->monolis,
			&sys->fe,
			&sys->basis,
			&sys->vals,
			te,
			st,
			space_map,
			T_prev_trace_global,
			t_window_start);

	BBFE_sys_monowrap_set_Dirichlet_bc_ST_partitioned_scalar(
			&sys->monolis,
			&sys->fe,
			&sys->bc,
			&sys->vals,
			te,
			st,
			space_map,
			t_window_start,
			sys->monolis.mat.R.B);

	ST_restore_Dirichlet_values_partitioned_scalar(
			&sys->fe,
			&sys->bc,
			&sys->vals,
			te,
			st,
			space_map,
			t_window_start,
			T_window);

	BBFE_sys_monowrap_solve(
			&sys->monolis,
			(MONOLIS_COM*)&st->monolis_com_st,
			T_window,
			MONOLIS_ITER_BICGSTAB,
			MONOLIS_PREC_DIAG,
			sys->vals.mat_max_iter,
			sys->vals.mat_epsilon);

	ST_restore_Dirichlet_values_partitioned_scalar(
			&sys->fe,
			&sys->bc,
			&sys->vals,
			te,
			st,
			space_map,
			t_window_start,
			T_window);
}

void ST_extract_slab_right_trace_global(
		const ST_TIME_ELEMENT* te,
		const ST_PARTITION* st,
		int slab,
		const double* T_window,
		double* T_global)
{
	if(te == NULL || st == NULL || T_window == NULL || T_global == NULL ||
	   slab < 0 || slab >= st->meta.num_window_slabs) {
		ST_fail("invalid argument in ST_extract_slab_right_trace_global");
	}

	const int nxg = st->meta.num_global_space_nodes;

	memset(
			T_global,
			0,
			(size_t)nxg * sizeof(double));

	for(int lid = 0; lid < st->num_internal_st_nodes; lid++) {
		if(st->slab_id[lid] != slab) {
			continue;
		}

		const int alpha = st->time_dof[lid];
		const int global_space = st->global_space_node_id[lid];

		if(global_space < 0 || global_space >= nxg) {
			ST_fail(
					"space global ID out of range while "
					"extracting ST trace");
		}

		T_global[global_space] +=
			te->e_plus[alpha] * T_window[lid];
	}

	monolis_allreduce_R(
			nxg,
			T_global,
			MONOLIS_MPI_SUM,
			st->monolis_com_st.comm);
}

void ST_copy_replicated_space_nodal_vector_global(
		const ST_SPACE_MAP* space_map,
		const double* local_vector,
		double* global_vector,
		int num_global_space_nodes)
{
	if(space_map == NULL ||
	   local_vector == NULL ||
	   global_vector == NULL ||
	   num_global_space_nodes <= 0) {
		ST_fail(
				"invalid argument in "
				"ST_copy_replicated_space_nodal_vector_global");
	}

	if(space_map->num_total_space_nodes !=
	   num_global_space_nodes) {
		fprintf(
				stderr,
				"ERROR: replicated spatial mesh size differs "
				"from global ST metadata.\n"
				"  spatial nodes = %d\n"
				"  ST global spatial nodes = %d\n",
				space_map->num_total_space_nodes,
				num_global_space_nodes);
		exit(EXIT_FAILURE);
	}

	memset(
			global_vector,
			0,
			(size_t)num_global_space_nodes * sizeof(double));

	for(int lid=0;
			lid<space_map->num_total_space_nodes;
			lid++) {

		const int64_t gid =
			space_map->local_to_global_space_node[lid];

		if(gid < 0 ||
		   gid >= num_global_space_nodes) {
			ST_fail(
					"space global ID out of range while "
					"copying replicated spatial vector");
		}

		global_vector[gid] = local_vector[lid];
	}

	/*
	 * No collective operation is performed here.
	 * All ranks already hold the same complete spatial field.
	 */
}

void ST_copy_global_vector_to_local_space(
		const ST_SPACE_MAP* space_map,
		const double* global_vector,
		double* local_vector)
{
	if(space_map == NULL || global_vector == NULL || local_vector == NULL) {
		ST_fail("invalid argument in ST_copy_global_vector_to_local_space");
	}
	for(int lid=0; lid<space_map->num_total_space_nodes; lid++) {
		const int64_t gid = space_map->local_to_global_space_node[lid];
		local_vector[lid] = global_vector[gid];
	}
}

void ST_collect_space_nodal_vector_global(
		const ST_SPACE_MAP* space_map,
		const FE_SYSTEM* sys,
		const double* local_vector,
		double* global_vector,
		int num_global_space_nodes)
{
	if(sys == NULL) {
		ST_fail(
				"NULL FE_SYSTEM in "
				"ST_collect_space_nodal_vector_global");
	}

	if(sys->fe.total_num_nodes !=
	   num_global_space_nodes) {
		fprintf(
				stderr,
				"ERROR: replicated BBFE spatial-node count "
				"differs from ST metadata.\n"
				"  BBFE nodes = %d\n"
				"  ST global spatial nodes = %d\n",
				sys->fe.total_num_nodes,
				num_global_space_nodes);
		exit(EXIT_FAILURE);
	}

	ST_copy_replicated_space_nodal_vector_global(
			space_map,
			local_vector,
			global_vector,
			num_global_space_nodes);
}