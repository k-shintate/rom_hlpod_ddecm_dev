#include "HROM_stage3_parallel_mode0.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "monolis_nnls_c.h"

namespace {

struct Stage3RowKey {
    int32_t snapshot;
    int64_t owner_global_id;
    int32_t mode_offset;

    bool operator<(const Stage3RowKey& other) const
    {
        if (snapshot != other.snapshot) return snapshot < other.snapshot;
        if (owner_global_id != other.owner_global_id) {
            return owner_global_id < other.owner_global_id;
        }
        return mode_offset < other.mode_offset;
    }
};

struct Stage3RankMeta {
    int source_rank;
    int nrow;
    int nselected;
    int nlocal;
    std::vector<Stage3RowKey> row_keys;
    std::vector<double> rhs;
    std::vector<int64_t> selected_global_ids;
    std::vector<int64_t> local_global_ids;
    long matrix_offset;
    std::string path;

    Stage3RankMeta()
        : source_rank(-1), nrow(0), nselected(0), nlocal(0),
          matrix_offset(0)
    {}
};


struct LegacyRouteEntry {
    int64_t original_global_id;
    int parallel_elem_id;
    int fine_subdomain_id;
    int has_D_bc;
};

struct LegacySelectedRecord {
    int elem_id;
    double weight;
};

static std::string join_path(const char* directory, const std::string& relative)
{
    if (directory == NULL || directory[0] == '\0' ||
        (directory[0] == '.' && directory[1] == '\0')) {
        return relative;
    }

    std::string out(directory);
    if (!out.empty() && out[out.size() - 1] != '/') out += '/';
    out += relative;
    return out;
}

[[noreturn]] static void fail_mpi(
    MPI_Comm comm,
    int rank,
    const std::string& message)
{
    std::fprintf(stderr, "[STAGE3-MODE0-MPI-ERROR] rank=%d %s\n",
                 rank, message.c_str());
    std::fflush(stderr);
    MPI_Abort(comm, EXIT_FAILURE);
    std::abort();
}

static bool file_exists(const std::string& path)
{
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (fp == NULL) return false;
    std::fclose(fp);
    return true;
}

static int detect_input_rank_count(const char* directory)
{
    int count = 0;
    for (;;) {
        char rel[256];
        std::snprintf(rel, sizeof(rel), "DDECM/stage3_rank.%d.bin", count);
        if (!file_exists(join_path(directory, rel))) break;
        ++count;
    }
    return count;
}

static bool read_exact(FILE* fp, void* ptr, size_t size, size_t count)
{
    return std::fread(ptr, size, count, fp) == count;
}

static Stage3RankMeta read_rank_meta(
    MPI_Comm comm,
    int solver_rank,
    const char* directory,
    int expected_source_rank)
{
    Stage3RankMeta meta;

    char rel[256];
    std::snprintf(rel, sizeof(rel),
                  "DDECM/stage3_rank.%d.bin", expected_source_rank);
    meta.path = join_path(directory, rel);

    FILE* fp = std::fopen(meta.path.c_str(), "rb");
    if (fp == NULL) {
        fail_mpi(comm, solver_rank,
                 std::string("cannot open ") + meta.path);
    }

    char magic[8];
    int32_t version = 0;
    int32_t source_rank32 = -1;
    int32_t nrow32 = -1;
    int32_t nselected32 = -1;
    int32_t nlocal32 = -1;

    if (!read_exact(fp, magic, sizeof(char), 8) ||
        !read_exact(fp, &version, sizeof(int32_t), 1) ||
        !read_exact(fp, &source_rank32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nrow32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nselected32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nlocal32, sizeof(int32_t), 1)) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("failed to read Stage-3 header: ") + meta.path);
    }

    const char expected_magic[8] = {'H','3','N','N','L','S','4','\0'};
    if (std::memcmp(magic, expected_magic, 8) != 0 || version != 4) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("unsupported Stage-3 rank-file format: ") +
                 meta.path);
    }

    if (source_rank32 != expected_source_rank ||
        nrow32 < 0 || nselected32 < 0 || nlocal32 < 0) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("invalid Stage-3 header values: ") + meta.path);
    }

    meta.source_rank = static_cast<int>(source_rank32);
    meta.nrow = static_cast<int>(nrow32);
    meta.nselected = static_cast<int>(nselected32);
    meta.nlocal = static_cast<int>(nlocal32);

    meta.row_keys.resize(static_cast<size_t>(meta.nrow));
    for (int r = 0; r < meta.nrow; ++r) {
        if (!read_exact(fp, &meta.row_keys[static_cast<size_t>(r)].snapshot,
                        sizeof(int32_t), 1) ||
            !read_exact(fp, &meta.row_keys[static_cast<size_t>(r)].owner_global_id,
                        sizeof(int64_t), 1) ||
            !read_exact(fp, &meta.row_keys[static_cast<size_t>(r)].mode_offset,
                        sizeof(int32_t), 1)) {
            std::fclose(fp);
            fail_mpi(comm, solver_rank,
                     std::string("failed to read row keys: ") + meta.path);
        }
    }

    meta.rhs.resize(static_cast<size_t>(meta.nrow), 0.0);
    if (meta.nrow > 0 &&
        !read_exact(fp, meta.rhs.data(), sizeof(double),
                    static_cast<size_t>(meta.nrow))) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("failed to read RHS: ") + meta.path);
    }

    meta.selected_global_ids.resize(
        static_cast<size_t>(meta.nselected), static_cast<int64_t>(-1));
    if (meta.nselected > 0 &&
        !read_exact(fp, meta.selected_global_ids.data(), sizeof(int64_t),
                    static_cast<size_t>(meta.nselected))) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("failed to read selected IDs: ") + meta.path);
    }

    /* Stage-2 weights are not an initial active set in parallel mode 0.
       Read-and-discard them only to advance to the next field. */
    if (meta.nselected > 0) {
        std::vector<double> discard_weights(
            static_cast<size_t>(meta.nselected), 0.0);
        if (!read_exact(fp, discard_weights.data(), sizeof(double),
                        static_cast<size_t>(meta.nselected))) {
            std::fclose(fp);
            fail_mpi(comm, solver_rank,
                     std::string("failed to read selected weights: ") +
                     meta.path);
        }
    }

    meta.local_global_ids.resize(
        static_cast<size_t>(meta.nlocal), static_cast<int64_t>(-1));
    if (meta.nlocal > 0 &&
        !read_exact(fp, meta.local_global_ids.data(), sizeof(int64_t),
                    static_cast<size_t>(meta.nlocal))) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("failed to read local element IDs: ") +
                 meta.path);
    }

    meta.matrix_offset = std::ftell(fp);
    if (meta.matrix_offset < 0) {
        std::fclose(fp);
        fail_mpi(comm, solver_rank,
                 std::string("ftell failed: ") + meta.path);
    }

    std::fclose(fp);
    return meta;
}

static int block_begin(int n, int rank, int nproc)
{
    return static_cast<int>((static_cast<long long>(n) * rank) / nproc);
}

static int block_end(int n, int rank, int nproc)
{
    return static_cast<int>((static_cast<long long>(n) * (rank + 1)) / nproc);
}


static bool read_stage3_route_file(
    const char* directory,
    int source_rank,
    std::vector<LegacyRouteEntry>& routes,
    std::string& error)
{
    char rel[256];
    std::snprintf(rel, sizeof(rel), "DDECM/stage3_route.%d.txt", source_rank);
    const std::string path = join_path(directory, rel);

    FILE* fp = std::fopen(path.c_str(), "r");
    if (fp == NULL) {
        error = std::string("cannot open legacy Stage-3 route file: ") + path;
        return false;
    }

    char magic[128];
    int nroute = -1;
    if (std::fscanf(fp, "%127s", magic) != 1 ||
        std::strcmp(magic, "HROM_STAGE3_ROUTE_V2") != 0 ||
        std::fscanf(fp, "%d", &nroute) != 1 ||
        nroute < 0) {
        std::fclose(fp);
        error = std::string("invalid legacy Stage-3 route header: ") + path;
        return false;
    }

    routes.clear();
    routes.reserve(static_cast<size_t>(nroute));

    for (int i = 0; i < nroute; ++i) {
        long long original = -1;
        LegacyRouteEntry e;
        e.original_global_id = -1;
        e.parallel_elem_id = -1;
        e.fine_subdomain_id = -1;
        e.has_D_bc = 0;

        if (std::fscanf(fp, "%lld %d %d %d",
                        &original,
                        &e.parallel_elem_id,
                        &e.fine_subdomain_id,
                        &e.has_D_bc) != 4) {
            std::fclose(fp);
            error = std::string("invalid legacy Stage-3 route record: ") + path;
            return false;
        }

        e.original_global_id = static_cast<int64_t>(original);
        routes.push_back(e);
    }

    std::fclose(fp);
    return true;
}


static bool file_exists_stage3(const std::string& path)
{
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (fp == NULL) return false;
    std::fclose(fp);
    return true;
}


static bool copy_file_if_missing_stage3(
    const std::string& src,
    const std::string& dst,
    std::string& error)
{
    if (file_exists_stage3(dst)) return true;

    FILE* in = std::fopen(src.c_str(), "rb");
    if (in == NULL) return true;

    FILE* out = std::fopen(dst.c_str(), "wb");
    if (out == NULL) {
        std::fclose(in);
        error = std::string("failed to create legacy backup file: ") + dst;
        return false;
    }

    char buf[65536];
    while (true) {
        const size_t n = std::fread(buf, 1, sizeof(buf), in);
        if (n > 0 && std::fwrite(buf, 1, n, out) != n) {
            std::fclose(in);
            std::fclose(out);
            error = std::string("failed to write legacy backup file: ") + dst;
            return false;
        }
        if (n < sizeof(buf)) {
            if (std::ferror(in)) {
                std::fclose(in);
                std::fclose(out);
                error = std::string("failed to read legacy backup source: ") + src;
                return false;
            }
            break;
        }
    }

    std::fclose(in);
    std::fclose(out);
    return true;
}


static std::string legacy_lb_path(
    const char* directory,
    int subdomain,
    bool dbc)
{
    char rel[256];
    std::snprintf(
        rel, sizeof(rel),
        dbc ? "DDECM/lb_selected_elem_D_bc.%d.txt"
            : "DDECM/lb_selected_elem.%d.txt",
        subdomain);
    return join_path(directory, rel);
}


/*
 * Exact legacy-output responsibility of HROM_stage3_serial_main.c:
 *
 *   final original-global Stage-3 selection
 *      -> stage3_route.<source-rank>.txt
 *      -> DDECM/lb_selected_elem[_D_bc].<fine-subdomain>.txt
 *
 * Do NOT create DDECM/selected_elem.<rank>.txt here.  In the conventional
 * workflow those rank files are produced later by the existing parallel-HROM
 * conversion path (HROM_ddecm_read_selected_elems_para()).
 */
static bool write_exact_serial_legacy_lb_files(
    const char* directory,
    int input_rank_files,
    const std::vector<int64_t>& candidate_ids,
    const std::vector<double>& x,
    double weight_tol,
    int* routed_selected_out,
    std::string& error)
{
    struct OutEntry {
        long long id;
        double w;
    };

    std::vector<std::pair<int64_t, double> > selected;
    selected.reserve(candidate_ids.size());

    for (size_t i = 0; i < candidate_ids.size(); ++i) {
        if (x[i] > weight_tol) {
            selected.push_back(
                std::make_pair(candidate_ids[i], x[i]));
        }
    }

    std::unordered_map<int64_t, std::vector<LegacyRouteEntry> >
        routes_by_original;
    std::unordered_set<int> subdomains;

    for (int src = 0; src < input_rank_files; ++src) {
        std::vector<LegacyRouteEntry> routes;
        if (!read_stage3_route_file(directory, src, routes, error)) {
            return false;
        }

        for (size_t i = 0; i < routes.size(); ++i) {
            const LegacyRouteEntry& e = routes[i];
            if (e.has_D_bc != 0 && e.has_D_bc != 1) {
                error = "Stage-3 route has invalid D_bc class";
                return false;
            }
            routes_by_original[e.original_global_id].push_back(e);
            subdomains.insert(e.fine_subdomain_id);
        }
    }

    std::unordered_map<int, std::vector<OutEntry> > no_out;
    std::unordered_map<int, std::vector<OutEntry> > dbc_out;
    std::unordered_map<int, std::unordered_set<long long> > seen_no;
    std::unordered_map<int, std::unordered_set<long long> > seen_dbc;

    size_t distributed_routes = 0;

    for (size_t p = 0; p < selected.size(); ++p) {
        const int64_t original = selected[p].first;
        const double w = selected[p].second;

        const std::unordered_map<
            int64_t,
            std::vector<LegacyRouteEntry> >::const_iterator it =
                routes_by_original.find(original);

        if (it == routes_by_original.end()) {
            std::fprintf(
                stderr,
                "ERROR: final original-global element %lld has no exact route\n",
                static_cast<long long>(original));
            error = "final Stage-3 element has no exact route";
            return false;
        }

        const std::vector<LegacyRouteEntry>& rr = it->second;
        for (size_t j = 0; j < rr.size(); ++j) {
            const LegacyRouteEntry& e = rr[j];
            const int s = e.fine_subdomain_id;
            const long long id = static_cast<long long>(e.parallel_elem_id);

            if (e.has_D_bc) {
                if (seen_dbc[s].insert(id).second) {
                    OutEntry oe;
                    oe.id = id;
                    oe.w = w;
                    dbc_out[s].push_back(oe);
                    distributed_routes++;
                }
            }
            else {
                if (seen_no[s].insert(id).second) {
                    OutEntry oe;
                    oe.id = id;
                    oe.w = w;
                    no_out[s].push_back(oe);
                    distributed_routes++;
                }
            }
        }
    }

    for (std::unordered_set<int>::const_iterator sit = subdomains.begin();
         sit != subdomains.end(); ++sit) {
        const int s = *sit;
        const std::string no_path = legacy_lb_path(directory, s, false);
        const std::string db_path = legacy_lb_path(directory, s, true);

        if (!copy_file_if_missing_stage3(
                no_path, no_path + ".pre_stage3_parallel.bak", error) ||
            !copy_file_if_missing_stage3(
                db_path, db_path + ".pre_stage3_parallel.bak", error)) {
            return false;
        }

        std::sort(
            no_out[s].begin(), no_out[s].end(),
            [](const OutEntry& a, const OutEntry& b) {
                return a.id < b.id;
            });
        std::sort(
            dbc_out[s].begin(), dbc_out[s].end(),
            [](const OutEntry& a, const OutEntry& b) {
                return a.id < b.id;
            });

        FILE* fp = std::fopen(no_path.c_str(), "w");
        if (fp == NULL) {
            error = std::string("cannot write legacy no-D_bc file: ") + no_path;
            return false;
        }
        std::fprintf(fp, "%zu\n", no_out[s].size());
        for (size_t i = 0; i < no_out[s].size(); ++i) {
            std::fprintf(
                fp, "%lld %.30e\n",
                no_out[s][i].id,
                no_out[s][i].w);
        }
        std::fclose(fp);

        fp = std::fopen(db_path.c_str(), "w");
        if (fp == NULL) {
            error = std::string("cannot write legacy D_bc file: ") + db_path;
            return false;
        }
        std::fprintf(fp, "%zu\n", dbc_out[s].size());
        for (size_t i = 0; i < dbc_out[s].size(); ++i) {
            std::fprintf(
                fp, "%lld %.30e\n",
                dbc_out[s][i].id,
                dbc_out[s][i].w);
        }
        std::fclose(fp);

        std::printf(
            "[STAGE3-MODE0-LEGACY-LB] subdomain=%d D_bc=%zu no_D_bc=%zu total=%zu\n",
            s,
            dbc_out[s].size(),
            no_out[s].size(),
            dbc_out[s].size() + no_out[s].size());
    }

    std::printf(
        "[STAGE3-MODE0-LEGACY-DONE] physical_selected=%zu distributed_exact_routes=%zu\n",
        selected.size(),
        distributed_routes);
    std::fflush(stdout);

    if (routed_selected_out != NULL) {
        *routed_selected_out = static_cast<int>(selected.size());
    }

    return true;
}


static void write_selected_output(
    MPI_Comm comm,
    int rank,
    const char* directory,
    const char* output_relative_path,
    const std::vector<int64_t>& candidate_ids,
    const std::vector<double>& x,
    double weight_tol)
{
    std::string relative =
        (output_relative_path != NULL && output_relative_path[0] != '\0')
        ? std::string(output_relative_path)
        : std::string("DDECM/stage3_selected_elem.txt");

    const std::string path = join_path(directory, relative);
    FILE* fp = std::fopen(path.c_str(), "w");
    if (fp == NULL) {
        fail_mpi(comm, rank, std::string("cannot write ") + path);
    }

    int selected = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] > weight_tol) ++selected;
    }

    std::fprintf(fp, "%d\n", selected);
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] <= weight_tol) continue;
        std::fprintf(fp, "%lld %.30e\n",
                     static_cast<long long>(candidate_ids[i]), x[i]);
    }
    std::fclose(fp);

    std::printf(
        "[STAGE3-MODE0-MPI-OUTPUT] selected=%d candidates=%zu file=%s\n",
        selected, candidate_ids.size(), path.c_str());
    std::fflush(stdout);
}

static void write_stage3_diagnostics(
    MPI_Comm comm,
    int rank,
    const char* directory,
    int input_rank_files,
    int solver_ranks,
    int rows,
    int candidates,
    int selected,
    int max_iter,
    double tol_outer,
    double tol_inner,
    double weight_tol,
    double scale,
    double residual_scaled,
    double residual_raw,
    double relative_residual)
{
    const std::string residual_path =
        join_path(directory, "DDECM/stage3_NNLS_residual.dat");
    FILE* fr = std::fopen(residual_path.c_str(), "w");
    if (fr == NULL) {
        fail_mpi(comm, rank, std::string("cannot write ") + residual_path);
    }
    /* Keep the residual file numeric-only for compatibility with simple readers. */
    std::fprintf(fr, "%.15e\n", residual_raw);
    std::fclose(fr);

    const std::string summary_path =
        join_path(directory, "DDECM/stage_reduction_summary.txt");
    FILE* fs = std::fopen(summary_path.c_str(), "w");
    if (fs == NULL) {
        fail_mpi(comm, rank, std::string("cannot write ") + summary_path);
    }
    std::fprintf(fs, "stage3_mode 0\n");
    std::fprintf(fs, "solver parallel_mpi_monolis_sparse_nnls\n");
    std::fprintf(fs, "input_rank_files %d\n", input_rank_files);
    std::fprintf(fs, "solver_ranks %d\n", solver_ranks);
    std::fprintf(fs, "online_output legacy_lb_selected_elem_files\n");
    std::fprintf(fs, "rows %d\n", rows);
    std::fprintf(fs, "candidates %d\n", candidates);
    std::fprintf(fs, "selected %d\n", selected);
    std::fprintf(fs, "max_iter %d\n", max_iter);
    std::fprintf(fs, "tol_outer %.15e\n", tol_outer);
    std::fprintf(fs, "tol_inner %.15e\n", tol_inner);
    std::fprintf(fs, "weight_tol %.15e\n", weight_tol);
    std::fprintf(fs, "normalization_scale %.15e\n", scale);
    std::fprintf(fs, "residual_scaled %.15e\n", residual_scaled);
    std::fprintf(fs, "residual_raw %.15e\n", residual_raw);
    std::fprintf(fs, "relative_residual %.15e\n", relative_residual);
    std::fclose(fs);
}

} // namespace

extern "C" int HROM_stage3_parallel_mode0_run(
    MPI_Comm comm,
    const char* directory,
    int input_rank_count,
    int max_iter,
    double tol_outer,
    double tol_inner,
    double weight_tol,
    const char* output_relative_path)
{
    int rank = 0;
    int nproc = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    if (max_iter <= 0) max_iter = 5000;
    if (!(tol_outer > 0.0) || !std::isfinite(tol_outer)) {
        tol_outer = 1.0e-10;
    }
    if (!std::isfinite(weight_tol) || weight_tol < 0.0) {
        weight_tol = 0.0;
    }

    int num_input_ranks = input_rank_count;
    if (rank == 0 && num_input_ranks <= 0) {
        num_input_ranks = detect_input_rank_count(directory);
    }
    MPI_Bcast(&num_input_ranks, 1, MPI_INT, 0, comm);

    if (num_input_ranks <= 0) {
        fail_mpi(comm, rank,
                 "no contiguous DDECM/stage3_rank.<rank>.bin files found");
    }

    /*
     * Metadata are intentionally read on every solver rank.  The metadata are
     * small relative to the matrices and this makes the global row/candidate
     * ordering deterministic without a custom serialization layer.
     */
    std::vector<Stage3RankMeta> metas;
    metas.reserve(static_cast<size_t>(num_input_ranks));

    std::set<Stage3RowKey> row_set;
    std::set<int64_t> candidate_set;

    for (int src = 0; src < num_input_ranks; ++src) {
        Stage3RankMeta meta = read_rank_meta(comm, rank, directory, src);
        for (size_t r = 0; r < meta.row_keys.size(); ++r) {
            row_set.insert(meta.row_keys[r]);
        }
        for (size_t i = 0; i < meta.selected_global_ids.size(); ++i) {
            candidate_set.insert(meta.selected_global_ids[i]);
        }
        metas.push_back(std::move(meta));
    }

    std::vector<Stage3RowKey> global_rows(row_set.begin(), row_set.end());
    std::vector<int64_t> candidate_ids(candidate_set.begin(), candidate_set.end());

    const int m = static_cast<int>(global_rows.size());
    const int n_global = static_cast<int>(candidate_ids.size());

    if (m <= 0) {
        fail_mpi(comm, rank, "global Stage-3 row set is empty");
    }
    if (n_global <= 0) {
        fail_mpi(comm, rank, "Stage-2 candidate union is empty");
    }
    std::map<Stage3RowKey, int> global_row_index;
    for (int r = 0; r < m; ++r) {
        global_row_index.emplace(global_rows[static_cast<size_t>(r)], r);
    }

    const int c_begin = block_begin(n_global, rank, nproc);
    const int c_end = block_end(n_global, rank, nproc);
    const int n_loc = c_end - c_begin;

    std::unordered_map<int64_t, int> owned_candidate_col;
    owned_candidate_col.reserve(static_cast<size_t>(2 * n_loc + 1));
    for (int c = c_begin; c < c_end; ++c) {
        owned_candidate_col.emplace(
            candidate_ids[static_cast<size_t>(c)], c - c_begin);
    }

    /* Row-major C storage; monolis_nnls_c wrapper accepts double** A[j][i]. */
    std::vector<double> A_storage(
        static_cast<size_t>(m) * static_cast<size_t>(n_loc), 0.0);
    std::vector<double*> A(static_cast<size_t>(m), NULL);
    for (int r = 0; r < m; ++r) {
        A[static_cast<size_t>(r)] =
            A_storage.data() + static_cast<size_t>(r) * n_loc;
    }

    std::vector<double> b(static_cast<size_t>(m), 0.0);

    /* Assemble the replicated global RHS. */
    for (size_t f = 0; f < metas.size(); ++f) {
        const Stage3RankMeta& meta = metas[f];
        for (int r = 0; r < meta.nrow; ++r) {
            const std::map<Stage3RowKey, int>::const_iterator it =
                global_row_index.find(meta.row_keys[static_cast<size_t>(r)]);
            if (it == global_row_index.end()) {
                fail_mpi(comm, rank, "internal global-row mapping failure");
            }
            b[static_cast<size_t>(it->second)] +=
                meta.rhs[static_cast<size_t>(r)];
        }
    }

    /*
     * Each solver rank reads all source rank files but keeps only its own
     * candidate columns.  Duplicate physical-element copies across source
     * rank files are accumulated into the same original-global candidate
     * column after row-key alignment.
     */
    for (size_t f = 0; f < metas.size(); ++f) {
        const Stage3RankMeta& meta = metas[f];

        std::vector<std::pair<int, int> > local_column_route;
        for (int c = 0; c < meta.nlocal; ++c) {
            const int64_t gid = meta.local_global_ids[static_cast<size_t>(c)];
            const std::unordered_map<int64_t, int>::const_iterator it =
                owned_candidate_col.find(gid);
            if (it != owned_candidate_col.end()) {
                local_column_route.push_back(std::make_pair(c, it->second));
            }
        }

        FILE* fp = std::fopen(meta.path.c_str(), "rb");
        if (fp == NULL) {
            fail_mpi(comm, rank,
                     std::string("cannot reopen ") + meta.path);
        }
        if (std::fseek(fp, meta.matrix_offset, SEEK_SET) != 0) {
            std::fclose(fp);
            fail_mpi(comm, rank,
                     std::string("cannot seek matrix in ") + meta.path);
        }

        std::vector<double> row_buffer(
            static_cast<size_t>(meta.nlocal), 0.0);

        for (int rr = 0; rr < meta.nrow; ++rr) {
            if (meta.nlocal > 0 &&
                !read_exact(fp, row_buffer.data(), sizeof(double),
                            static_cast<size_t>(meta.nlocal))) {
                std::fclose(fp);
                fail_mpi(comm, rank,
                         std::string("failed reading matrix in ") + meta.path);
            }

            if (local_column_route.empty()) continue;

            const int grow = global_row_index.find(
                meta.row_keys[static_cast<size_t>(rr)])->second;

            double* const dst = A[static_cast<size_t>(grow)];
            for (size_t q = 0; q < local_column_route.size(); ++q) {
                const int src_col = local_column_route[q].first;
                const int dst_col = local_column_route[q].second;
                dst[dst_col] += row_buffer[static_cast<size_t>(src_col)];
            }
        }

        std::fclose(fp);
    }

    /* Check dimensions are identical on all solver ranks. */
    int m_min = 0, m_max = 0, n_min = 0, n_max = 0;
    MPI_Allreduce(&m, &m_min, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(&m, &m_max, 1, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(&n_global, &n_min, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(&n_global, &n_max, 1, MPI_INT, MPI_MAX, comm);
    if (m_min != m_max || n_min != n_max) {
        fail_mpi(comm, rank,
                 "solver ranks constructed inconsistent global dimensions");
    }

    /* Common positive scaling leaves the mathematical NNLS minimizer intact. */
    double local_frob_sq = 0.0;
    for (size_t i = 0; i < A_storage.size(); ++i) {
        const double v = A_storage[i];
        if (!std::isfinite(v)) {
            fail_mpi(comm, rank, "non-finite value in assembled local A");
        }
        local_frob_sq += v * v;
    }

    double global_frob_sq = 0.0;
    MPI_Allreduce(&local_frob_sq, &global_frob_sq,
                  1, MPI_DOUBLE, MPI_SUM, comm);

    double bnorm_raw_sq = 0.0;
    for (int r = 0; r < m; ++r) {
        const double v = b[static_cast<size_t>(r)];
        if (!std::isfinite(v)) {
            fail_mpi(comm, rank, "non-finite value in assembled global b");
        }
        bnorm_raw_sq += v * v;
    }
    const double bnorm_raw = std::sqrt(std::max(0.0, bnorm_raw_sq));

    double scale = std::sqrt(std::max(0.0, global_frob_sq));
    if (!std::isfinite(scale) || scale <= DBL_MIN) scale = 1.0;

    for (size_t i = 0; i < A_storage.size(); ++i) A_storage[i] /= scale;
    for (int r = 0; r < m; ++r) b[static_cast<size_t>(r)] /= scale;

    if (!(tol_inner > 0.0) || !std::isfinite(tol_inner)) {
        tol_inner = 1.0e-14 * static_cast<double>(std::max(m, n_global));
    }

    if (rank == 0) {
        std::printf(
            "[STAGE3-MODE0-MPI-START] input_rank_files=%d solver_ranks=%d "
            "rows=%d candidates=%d max_iter=%d tol_outer=%.15e "
            "tol_inner=%.15e scale=%.15e bnorm_raw=%.15e\n",
            num_input_ranks, nproc, m, n_global, max_iter,
            tol_outer, tol_inner, scale, bnorm_raw);
        std::fflush(stdout);
    }

    std::printf(
        "[STAGE3-MODE0-MPI-BLOCK] rank=%d global_cols=[%d,%d) n_loc=%d\n",
        rank, c_begin, c_end, n_loc);
    std::fflush(stdout);

    std::vector<double> x_local(static_cast<size_t>(n_loc), 0.0);
    double residual_scaled = 0.0;

    const MPI_Fint fcomm = MPI_Comm_c2f(comm);
    const int monolis_comm = static_cast<int>(fcomm);

    monolis_optimize_nnls_R_with_sparse_solution(
        A.data(),
        b.data(),
        x_local.data(),
        m,
        n_loc,
        max_iter,
        tol_outer,
        tol_inner,
        &residual_scaled,
        monolis_comm);

    std::vector<int> recvcounts;
    std::vector<int> displs;
    std::vector<double> x_global;
    if (rank == 0) {
        recvcounts.resize(static_cast<size_t>(nproc), 0);
        displs.resize(static_cast<size_t>(nproc), 0);
        for (int p = 0; p < nproc; ++p) {
            const int pb = block_begin(n_global, p, nproc);
            const int pe = block_end(n_global, p, nproc);
            recvcounts[static_cast<size_t>(p)] = pe - pb;
            displs[static_cast<size_t>(p)] = pb;
        }
        x_global.resize(static_cast<size_t>(n_global), 0.0);
    }

    MPI_Gatherv(
        x_local.data(), n_loc, MPI_DOUBLE,
        rank == 0 ? x_global.data() : NULL,
        rank == 0 ? recvcounts.data() : NULL,
        rank == 0 ? displs.data() : NULL,
        MPI_DOUBLE, 0, comm);

    int legacy_output_ok = 1;

    if (rank == 0) {
        const double residual_raw = residual_scaled * scale;
        const double relative_residual =
            (bnorm_raw > DBL_MIN) ? residual_raw / bnorm_raw : residual_raw;

        int selected = 0;
        for (size_t i = 0; i < x_global.size(); ++i) {
            if (x_global[i] > weight_tol) ++selected;
        }

        std::printf(
            "[STAGE3-MODE0-MPI-DONE] selected=%d/%d "
            "residual_scaled=%.15e residual_raw=%.15e "
            "relative_residual=%.15e\n",
            selected, n_global, residual_scaled, residual_raw,
            relative_residual);
        std::fflush(stdout);

        /*
         * Legacy online-HROM output is authoritative.  The old one-file
         * original-global selection is now optional diagnostics only.
         */
        int routed_selected = 0;
        std::string legacy_error;
        const bool legacy_ok = write_exact_serial_legacy_lb_files(
            directory,
            num_input_ranks,
            candidate_ids,
            x_global,
            weight_tol,
            &routed_selected,
            legacy_error);

        if (!legacy_ok) {
            legacy_output_ok = 0;
            std::fprintf(
                stderr,
                "[STAGE3-MODE0-LEGACY-ERROR] %s\n",
                legacy_error.c_str());
            std::fflush(stderr);
        }

        if (output_relative_path != NULL &&
            output_relative_path[0] != '\0') {
            write_selected_output(
                comm, rank, directory, output_relative_path,
                candidate_ids, x_global, weight_tol);
        }

        write_stage3_diagnostics(
            comm, rank, directory, num_input_ranks, nproc, m, n_global,
            selected, max_iter, tol_outer, tol_inner, weight_tol, scale,
            residual_scaled, residual_raw, relative_residual);
    }

    MPI_Bcast(&legacy_output_ok, 1, MPI_INT, 0, comm);
    MPI_Barrier(comm);

    if (!legacy_output_ok) {
        return EXIT_FAILURE;
    }
    return 0;
}
