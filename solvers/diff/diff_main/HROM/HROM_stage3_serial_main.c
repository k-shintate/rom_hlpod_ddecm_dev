/*
 * Stage-3 serial NNLS version 7.
 *
 * Fixes two structural bugs from the previous version:
 *   1) Cross-rank rows are matched by a GLOBAL row key
 *      (snapshot, global mode-owner ID, mode offset), not rank-local j.
 *   2) RHS contributions AND element-column contributions are both summed
 *      across all MPI-rank files.  No adopted-rank/owner column selection.
 *
 * After solving Stage 3, this program uses stage3_route.<rank>.txt written by
 * the MPI stage to return the final weights to the exact legacy
 * DDECM/lb_selected_elem*.txt files used by the parallel online calculation.
 */
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "monolis_nnls_c.h"

struct RowKey {
    int32_t snapshot = 0;
    int64_t owner_global_id = 0;
    int32_t mode_offset = 0;

    bool operator==(const RowKey& other) const {
        return snapshot == other.snapshot &&
               owner_global_id == other.owner_global_id &&
               mode_offset == other.mode_offset;
    }
};

struct RowKeyHash {
    size_t operator()(const RowKey& k) const {
        size_t h = std::hash<long long>()((long long)k.owner_global_id);
        h ^= std::hash<int>()((int)k.snapshot) + 0x9e3779b9u + (h<<6) + (h>>2);
        h ^= std::hash<int>()((int)k.mode_offset) + 0x9e3779b9u + (h<<6) + (h>>2);
        return h;
    }
};

struct Stage3RankData {
    int rank = -1;
    std::vector<RowKey> row_key;
    std::vector<double> rhs;
    std::vector<int64_t> selected_elem_id;
    std::vector<double> selected_weight;
    std::vector<int64_t> local_elem_id;
    std::vector<double> local_matrix; /* row-major */
};

struct RouteEntry {
    int rank = -1;
    int fine_subdomain = -1;
    int64_t original_global_id = -1;
    int64_t parallel_elem_id = -1;
    int has_D_bc = -1;
};

struct InterfaceKey {
    int64_t source_global_id = -1;
    int64_t target_global_id = -1;

    bool operator==(const InterfaceKey& other) const {
        return source_global_id == other.source_global_id &&
               target_global_id == other.target_global_id;
    }
};

struct InterfaceKeyHash {
    size_t operator()(const InterfaceKey& k) const {
        size_t h = std::hash<long long>()(
            (long long)k.source_global_id);
        h ^= std::hash<long long>()(
            (long long)k.target_global_id)
            + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};

struct InterfaceGroup {
    InterfaceKey key;
    std::vector<int64_t> member_original_id;
};


struct Stage1SubdomainStats {
    int rank = -1;
    int fine_subdomain = -1;
    int local_unique = 0;
    int local_candidates = 0;
    int initial_active = 0;
    double residual = 0.0;
    std::vector<int64_t> candidate_parallel_ids;
};

struct Stage1RouteKey {
    int rank = -1;
    int fine_subdomain = -1;
    int64_t parallel_elem_id = -1;

    bool operator==(const Stage1RouteKey& other) const {
        return rank == other.rank &&
               fine_subdomain == other.fine_subdomain &&
               parallel_elem_id == other.parallel_elem_id;
    }
};

struct Stage1RouteKeyHash {
    size_t operator()(const Stage1RouteKey& k) const {
        size_t h = std::hash<int>()(k.rank);
        h ^= std::hash<int>()(k.fine_subdomain)
            + 0x9e3779b9u + (h<<6) + (h>>2);
        h ^= std::hash<long long>()((long long)k.parallel_elem_id)
            + 0x9e3779b9u + (h<<6) + (h>>2);
        return h;
    }
};

static std::string join_path(const char* directory, const std::string& relative)
{
    std::string base = (directory != NULL) ? directory : ".";
    if (!base.empty() && base.back() != '/') base.push_back('/');
    return base + relative;
}

static void die(const char* msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}

static bool file_exists(const std::string& path)
{
    FILE* fp = fopen(path.c_str(), "rb");
    if (fp == NULL) return false;
    fclose(fp);
    return true;
}

static bool read_exact(FILE* fp, void* p, size_t size, size_t count)
{
    return fread(p, size, count, fp) == count;
}

static void copy_file_if_missing(const std::string& src, const std::string& dst)
{
    if (file_exists(dst)) return;
    FILE* in = fopen(src.c_str(), "rb");
    if (in == NULL) return;
    FILE* out = fopen(dst.c_str(), "wb");
    if (out == NULL) {
        fclose(in);
        die("failed to create backup file");
    }
    char buf[65536];
    while (true) {
        size_t n = fread(buf, 1, sizeof(buf), in);
        if (n > 0 && fwrite(buf, 1, n, out) != n) {
            fclose(in); fclose(out); die("failed to write backup file");
        }
        if (n < sizeof(buf)) {
            if (ferror(in)) { fclose(in); fclose(out); die("failed to read backup source"); }
            break;
        }
    }
    fclose(in);
    fclose(out);
}

static Stage3RankData read_rank_file(const char* directory, int expected_rank)
{
    char rel[256];
    snprintf(rel, sizeof(rel), "DDECM/stage3_rank.%d.bin", expected_rank);
    const std::string path = join_path(directory, rel);
    FILE* fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    char magic[8];
    int32_t version=0, rank32=-1, nrow32=0, nselected32=0, nlocal32=0;
    if (!read_exact(fp, magic, 1, 8) ||
        !read_exact(fp, &version, sizeof(int32_t), 1) ||
        !read_exact(fp, &rank32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nrow32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nselected32, sizeof(int32_t), 1) ||
        !read_exact(fp, &nlocal32, sizeof(int32_t), 1)) {
        fclose(fp); die("failed to read Stage-3 v4 header");
    }

    const char expected_magic[8] = {'H','3','N','N','L','S','4','\0'};
    if (memcmp(magic, expected_magic, 8) != 0 || version != 4) {
        fclose(fp);
        die("Stage-3 rank files are not version 4; rerun MPI with the global-row-key writer");
    }
    if (rank32 != expected_rank || nrow32 < 0 || nselected32 < 0 || nlocal32 < 0) {
        fclose(fp); die("invalid Stage-3 v4 dimensions/rank");
    }

    Stage3RankData d;
    d.rank = rank32;
    d.row_key.resize((size_t)nrow32);
    d.rhs.resize((size_t)nrow32);
    d.selected_elem_id.resize((size_t)nselected32);
    d.selected_weight.resize((size_t)nselected32);
    d.local_elem_id.resize((size_t)nlocal32);
    d.local_matrix.resize((size_t)nrow32 * (size_t)nlocal32);

    for (int r=0; r<nrow32; ++r) {
        if (!read_exact(fp, &d.row_key[(size_t)r].snapshot, sizeof(int32_t), 1) ||
            !read_exact(fp, &d.row_key[(size_t)r].owner_global_id, sizeof(int64_t), 1) ||
            !read_exact(fp, &d.row_key[(size_t)r].mode_offset, sizeof(int32_t), 1)) {
            fclose(fp); die("failed to read global row keys");
        }
    }

    if (nrow32 > 0 && !read_exact(fp, d.rhs.data(), sizeof(double), (size_t)nrow32)) {
        fclose(fp); die("failed to read Stage-3 RHS");
    }
    if (nselected32 > 0 && !read_exact(fp, d.selected_elem_id.data(), sizeof(int64_t), (size_t)nselected32)) {
        fclose(fp); die("failed to read Stage-2 selected IDs");
    }
    if (nselected32 > 0 && !read_exact(fp, d.selected_weight.data(), sizeof(double), (size_t)nselected32)) {
        fclose(fp); die("failed to read Stage-2 selected weights");
    }
    if (nlocal32 > 0 && !read_exact(fp, d.local_elem_id.data(), sizeof(int64_t), (size_t)nlocal32)) {
        fclose(fp); die("failed to read rank-local original-global IDs");
    }
    const size_t matrix_count=(size_t)nrow32*(size_t)nlocal32;
    if (matrix_count > 0 && !read_exact(fp, d.local_matrix.data(), sizeof(double), matrix_count)) {
        fclose(fp); die("failed to read Stage-3 rank matrix");
    }
    fclose(fp);

    printf("[read] rank=%d rows=%d stage2_selected=%d local_columns=%d\n",
        expected_rank, nrow32, nselected32, nlocal32);
    return d;
}

static std::vector<RouteEntry> read_all_routes(const char* directory, int num_ranks)
{
    std::vector<RouteEntry> routes;
    for (int rank=0; rank<num_ranks; ++rank) {
        char rel[256];
        snprintf(rel, sizeof(rel), "DDECM/stage3_route.%d.txt", rank);
        const std::string path=join_path(directory, rel);
        FILE* fp=fopen(path.c_str(), "r");
        if (fp==NULL) {
            fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }
        char route_magic[64];
        if (fscanf(fp, "%63s", route_magic) != 1 ||
            strcmp(route_magic, "HROM_STAGE3_ROUTE_V2") != 0) {
            fclose(fp);
            fprintf(stderr,
                "ERROR: %s is an old Stage-3 route file.\n"
                "Re-run the MPI Stage-1/Stage-2 preprocessing with the "
                "route-D_bc writer enabled.\n",
                path.c_str());
            exit(EXIT_FAILURE);
        }

        int n=-1;
        if (fscanf(fp, "%d", &n)!=1 || n<0) { fclose(fp); die("invalid route header"); }
        for (int i=0; i<n; ++i) {
            long long original=-1, parallel=-1;
            int s=-1, dbc=-1;
            if (fscanf(fp, "%lld %lld %d %d",
                    &original, &parallel, &s, &dbc)!=4 ||
                (dbc != 0 && dbc != 1)) {
                fclose(fp); die("invalid Stage-3 route V2 entry");
            }
            RouteEntry e;
            e.rank=rank;
            e.original_global_id=(int64_t)original;
            e.parallel_elem_id=(int64_t)parallel;
            e.fine_subdomain=s;
            e.has_D_bc=dbc;
            routes.push_back(e);
        }
        fclose(fp);
    }
    printf("[route] total exact routes=%zu\n", routes.size());
    return routes;
}



static std::vector<InterfaceGroup> read_all_interface_groups(
    const char* directory,
    int num_ranks)
{
    /*
     * Interface coverage is imposed on a PHYSICAL (undirected) interface.
     *
     * stage3_interface.* is written from rank-local directed adjacency and may
     * contain asymmetric membership information:
     *
     *     (k,l): { ... }
     *     (l,k): { ... }
     *
     * A physical FE interface is the unordered pair {k,l}.  Therefore merge
     * both directions using the canonical key
     *
     *     (min(k,l), max(k,l))
     *
     * and take the union of their eligible element IDs.
     *
     * This is also important numerically: a directed group can be empty even
     * when the reverse direction contains valid interface-supporting columns.
     * The hard constraint is then imposed once per physical interface.
     */
    std::unordered_map<
        InterfaceKey,
        std::unordered_set<int64_t>,
        InterfaceKeyHash> merged;

    size_t directed_records = 0;
    size_t directed_empty_records = 0;
    size_t directed_memberships = 0;

    for (int rank = 0; rank < num_ranks; ++rank) {
        char rel[256];
        snprintf(
            rel,
            sizeof(rel),
            "DDECM/stage3_interface.%d.txt",
            rank);

        const std::string path =
            join_path(directory, rel);

        FILE* fp = fopen(path.c_str(), "r");

        if (fp == NULL) {
            fprintf(
                stderr,
                "ERROR: cannot open %s\n",
                path.c_str());
            exit(EXIT_FAILURE);
        }

        char magic[64];
        int ng = -1;

        if (fscanf(fp, "%63s", magic) != 1 ||
            strcmp(magic, "HROM_STAGE3_INTERFACE_V1") != 0 ||
            fscanf(fp, "%d", &ng) != 1 ||
            ng < 0) {
            fclose(fp);
            die("invalid Stage-3 interface-group header");
        }

        for (int g = 0; g < ng; ++g) {
            long long source = -1;
            long long target = -1;
            int nm = -1;

            if (fscanf(
                    fp,
                    "%lld %lld %d",
                    &source,
                    &target,
                    &nm) != 3 ||
                nm < 0) {
                fclose(fp);
                die("invalid Stage-3 interface-group row");
            }

            directed_records++;
            directed_memberships += (size_t)nm;

            if (nm == 0) {
                directed_empty_records++;
            }

            const int64_t a =
                std::min((int64_t)source, (int64_t)target);
            const int64_t b =
                std::max((int64_t)source, (int64_t)target);

            if (a == b) {
                fclose(fp);
                fprintf(
                    stderr,
                    "ERROR: self-interface (%lld,%lld) in %s\n",
                    source,
                    target,
                    path.c_str());
                exit(EXIT_FAILURE);
            }

            InterfaceKey key;
            key.source_global_id = a;
            key.target_global_id = b;

            /*
             * Create the physical group even when this directed record has
             * zero members.  The reverse direction or another rank may supply
             * members later.
             */
            std::unordered_set<int64_t>& members =
                merged[key];

            for (int k = 0; k < nm; ++k) {
                long long elem = -1;

                if (fscanf(fp, "%lld", &elem) != 1) {
                    fclose(fp);
                    die("invalid Stage-3 interface member");
                }

                members.insert((int64_t)elem);
            }
        }

        fclose(fp);
    }

    std::vector<InterfaceGroup> groups;
    groups.reserve(merged.size());

    size_t total_membership = 0;
    size_t empty_physical_groups = 0;

    for (const auto& kv : merged) {
        InterfaceGroup g;
        g.key = kv.first;
        g.member_original_id.assign(
            kv.second.begin(),
            kv.second.end());

        std::sort(
            g.member_original_id.begin(),
            g.member_original_id.end());

        if (g.member_original_id.empty()) {
            empty_physical_groups++;
        }

        total_membership +=
            g.member_original_id.size();

        groups.push_back(std::move(g));
    }

    std::sort(
        groups.begin(),
        groups.end(),
        [](const InterfaceGroup& a,
           const InterfaceGroup& b) {
            if (a.key.source_global_id !=
                b.key.source_global_id) {
                return a.key.source_global_id <
                       b.key.source_global_id;
            }

            return a.key.target_global_id <
                   b.key.target_global_id;
        });

    printf(
        "[interface-groups] directed_records=%zu "
        "directed_empty=%zu directed_memberships=%zu "
        "physical=%zu physical_memberships=%zu "
        "physical_empty=%zu\n",
        directed_records,
        directed_empty_records,
        directed_memberships,
        groups.size(),
        total_membership,
        empty_physical_groups);

    /*
     * A physical interface with no member in either direction cannot satisfy
     * sum_{i in I_g} x_i >= epsilon_g.  Fail here with a precise diagnostic
     * instead of later reporting "no available Stage-3 element column".
     */
    if (empty_physical_groups > 0) {
        for (const InterfaceGroup& g : groups) {
            if (g.member_original_id.empty()) {
                fprintf(
                    stderr,
                    "ERROR: physical interface {%lld,%lld} has no eligible "
                    "element in either directed record.  The interface-member "
                    "construction in the MPI preprocessing must be corrected "
                    "for this interface.\n",
                    (long long)g.key.source_global_id,
                    (long long)g.key.target_global_id);
            }
        }
        exit(EXIT_FAILURE);
    }

    return groups;
}


static std::vector<Stage1SubdomainStats> read_all_stage1_stats(
    const char* directory,
    int num_ranks)
{
    std::vector<Stage1SubdomainStats> stats;

    for (int rank = 0; rank < num_ranks; ++rank) {
        char rel[256];
        snprintf(
            rel,
            sizeof(rel),
            "DDECM/stage1_stats.%d.txt",
            rank);

        const std::string path =
            join_path(directory, rel);

        FILE* fp = fopen(path.c_str(), "r");
        if (fp == NULL) {
            fprintf(
                stderr,
                "ERROR: cannot open %s\n"
                "Re-run the MPI Stage-1/Stage-2 calculation with the "
                "Stage-1 statistics writer enabled.\n",
                path.c_str());
            exit(EXIT_FAILURE);
        }

        char magic[64];
        if (fscanf(fp, "%63s", magic) != 1 ||
            strcmp(magic, "HROM_STAGE1_STATS_V1") != 0) {
            fclose(fp);
            die("invalid Stage-1 stats magic");
        }

        int nsub = -1;
        if (fscanf(fp, "%d", &nsub) != 1 || nsub < 0) {
            fclose(fp);
            die("invalid Stage-1 stats subdomain count");
        }

        for (int m = 0; m < nsub; ++m) {
            Stage1SubdomainStats st;
            st.rank = rank;

            if (fscanf(
                    fp,
                    "%d %d %d %d %lf",
                    &st.fine_subdomain,
                    &st.local_unique,
                    &st.local_candidates,
                    &st.initial_active,
                    &st.residual) != 5) {
                fclose(fp);
                die("invalid Stage-1 subdomain stats row");
            }

            if (st.local_unique < 0 ||
                st.local_candidates < 0 ||
                st.local_candidates > st.local_unique) {
                fclose(fp);
                die("invalid Stage-1 subdomain counts");
            }

            st.candidate_parallel_ids.resize(
                (size_t)st.local_candidates);

            for (int k = 0;
                 k < st.local_candidates;
                 ++k) {

                long long id = -1;

                if (fscanf(fp, "%lld", &id) != 1) {
                    fclose(fp);
                    die("invalid Stage-1 candidate ID");
                }

                st.candidate_parallel_ids[(size_t)k] =
                    (int64_t)id;
            }

            stats.push_back(st);
        }

        fclose(fp);
    }

    std::sort(
        stats.begin(),
        stats.end(),
        [](const Stage1SubdomainStats& a,
           const Stage1SubdomainStats& b) {
            if (a.fine_subdomain != b.fine_subdomain) {
                return a.fine_subdomain < b.fine_subdomain;
            }
            return a.rank < b.rank;
        });

    printf(
        "[stage1-stats] total fine-subdomain records=%zu\n",
        stats.size());

    return stats;
}

static std::string legacy_path(const char* directory, int subdomain, bool dbc)
{
    char rel[256];
    snprintf(rel, sizeof(rel),
        dbc ? "DDECM/lb_selected_elem_D_bc.%d.txt" : "DDECM/lb_selected_elem.%d.txt",
        subdomain);
    return join_path(directory, rel);
}

/* 0=no-D_bc, 1=D_bc.  Prefer the pristine Stage-2 backup if available. */
static void read_class_file(
    const char* directory,
    int subdomain,
    bool dbc,
    std::unordered_map<long long,int>& class_by_sub_parallel,
    std::unordered_map<long long,int>& global_class_by_parallel)
{
    const std::string current=legacy_path(directory, subdomain, dbc);
    const std::string stage2=current+".stage2.bak";
    const std::string path=file_exists(stage2) ? stage2 : current;
    FILE* fp=fopen(path.c_str(), "r");
    if (fp==NULL) return;
    int n=0;
    if (fscanf(fp, "%d", &n)!=1 || n<0) { fclose(fp); die("invalid legacy selection header"); }
    for (int i=0; i<n; ++i) {
        long long id=-1;
        double w=0.0;
        if (fscanf(fp, "%lld %lf", &id, &w)!=2) { fclose(fp); die("invalid legacy selection row"); }
        (void)w;
        /* pack subdomain + parallel id (parallel IDs are assumed nonnegative int range in legacy files) */
        const unsigned long long packed=
            ((unsigned long long)(unsigned int)subdomain<<32) |
            (unsigned long long)(unsigned int)id;
        class_by_sub_parallel[(long long)packed]=dbc?1:0;
        auto it=global_class_by_parallel.find(id);
        if (it==global_class_by_parallel.end()) global_class_by_parallel.emplace(id, dbc?1:0);
        else if (dbc) it->second=1;
    }
    fclose(fp);
}

static void write_parallel_legacy_files(
    const char* directory,
    const std::vector<RouteEntry>& routes,
    const std::vector<std::pair<int64_t,double>>& selected)
{
    std::unordered_map<int64_t,std::vector<RouteEntry>> routes_by_original;
    std::unordered_set<int> subdomains;

    for (const RouteEntry& e: routes) {
        if (e.has_D_bc != 0 && e.has_D_bc != 1) {
            die("Stage-3 route has invalid D_bc class");
        }
        routes_by_original[e.original_global_id].push_back(e);
        subdomains.insert(e.fine_subdomain);
    }

    struct OutEntry { long long id; double w; };
    std::unordered_map<int,std::vector<OutEntry>> no_out, dbc_out;
    std::unordered_map<int,std::unordered_set<long long>> seen_no, seen_dbc;

    size_t distributed_routes = 0;

    for (const auto& p: selected) {
        const int64_t original=p.first;
        const double w=p.second;

        auto it=routes_by_original.find(original);
        if (it==routes_by_original.end()) {
            fprintf(stderr,
                "ERROR: final original-global element %lld has no exact route\n",
                (long long)original);
            exit(EXIT_FAILURE);
        }

        for (const RouteEntry& e: it->second) {
            const int s=e.fine_subdomain;
            const long long id=(long long)e.parallel_elem_id;

            if (e.has_D_bc) {
                if (seen_dbc[s].insert(id).second) {
                    dbc_out[s].push_back({id,w});
                    distributed_routes++;
                }
            }
            else {
                if (seen_no[s].insert(id).second) {
                    no_out[s].push_back({id,w});
                    distributed_routes++;
                }
            }
        }
    }

    for (int s: subdomains) {
        const std::string no_path=legacy_path(directory,s,false);
        const std::string db_path=legacy_path(directory,s,true);
        copy_file_if_missing(no_path, no_path+".pre_stage3_v10.bak");
        copy_file_if_missing(db_path, db_path+".pre_stage3_v10.bak");

        std::sort(no_out[s].begin(), no_out[s].end(),
            [](const OutEntry&a,const OutEntry&b){return a.id<b.id;});
        std::sort(dbc_out[s].begin(), dbc_out[s].end(),
            [](const OutEntry&a,const OutEntry&b){return a.id<b.id;});

        FILE* fp=fopen(no_path.c_str(),"w");
        if (fp==NULL) die("cannot write no-D_bc legacy file");
        fprintf(fp,"%zu\n",no_out[s].size());
        for (const OutEntry&e:no_out[s]) fprintf(fp,"%lld %.30e\n",e.id,e.w);
        fclose(fp);

        fp=fopen(db_path.c_str(),"w");
        if (fp==NULL) die("cannot write D_bc legacy file");
        fprintf(fp,"%zu\n",dbc_out[s].size());
        for (const OutEntry&e:dbc_out[s]) fprintf(fp,"%lld %.30e\n",e.id,e.w);
        fclose(fp);

        printf("[legacy-write-v10] subdomain=%d D_bc=%zu no_D_bc=%zu total=%zu\n",
            s, dbc_out[s].size(), no_out[s].size(),
            dbc_out[s].size()+no_out[s].size());
    }

    printf("[legacy-write-v10] unique physical Stage-3 elements=%zu distributed exact routes=%zu\n",
        selected.size(), distributed_routes);
}

static double** allocate_matrix(int m,int n,double** block_out)
{
    if (m<=0 || n<=0) return NULL;
    double** A=(double**)malloc(sizeof(double*)*(size_t)m);
    double* block=(double*)calloc((size_t)m*(size_t)n,sizeof(double));
    if (A==NULL || block==NULL) { free(A); free(block); die("failed to allocate matrix"); }
    for (int r=0;r<m;++r) A[r]=block+(size_t)r*(size_t)n;
    *block_out=block;
    return A;
}

int main(int argc,char** argv)
{
    if (argc < 3) {
        fprintf(
            stderr,
            "usage: %s <directory> <num_rank_files> "
            "[inner_max_iter=1000] [inner_tol=1e-8] "
            "[weight_tol=1e-8] [coverage_tau=weight_tol] "
            "[admm_max_iter=100] [rho=1.0] "
            "[admm_abs_tol=1e-10] [admm_rel_tol=1e-6] "
            "[verbose=1] [coverage_mode=2]\n"
            "  coverage_mode=0 : baseline Stage-3 (Stage-2 candidates only, no coverage)\n"
            "  coverage_mode=1 : coverage candidate expansion only, no hard constraint\n"
            "  coverage_mode=2 : proposed method (candidate expansion + hard coverage)\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    const char* directory = argv[1];
    const int num_ranks = atoi(argv[2]);

    const int requested_iter =
        (argc > 3) ? atoi(argv[3]) : 1000;

    const int max_iter =
        (requested_iter > 0)
        ? std::min(requested_iter, 1000)
        : 1000;

    const double tol =
        (argc > 4) ? atof(argv[4]) : 1.0e-8;

    const double weight_tol =
        (argc > 5) ? atof(argv[5]) : tol;

    const double coverage_tau =
        (argc > 6) ? atof(argv[6]) : weight_tol;

    const int admm_max_iter =
        (argc > 7) ? atoi(argv[7]) : 100;

    const double admm_rho =
        (argc > 8) ? atof(argv[8]) : 1.0;

    const double admm_abs_tol =
        (argc > 9) ? atof(argv[9]) : 1.0e-10;

    const double admm_rel_tol =
        (argc > 10) ? atof(argv[10]) : 1.0e-6;

    const int admm_verbose =
        (argc > 11) ? atoi(argv[11]) : 1;

    /*
     * coverage_mode:
     *   0 = exact comparison baseline:
     *       use only the Stage-2 selected union as Stage-3 candidates and
     *       solve the original unconstrained NNLS.
     *
     *   1 = candidate-expansion ablation:
     *       add one eligible candidate to an interface that has no Stage-2
     *       candidate, but solve the original unconstrained NNLS.
     *
     *   2 = proposed method:
     *       candidate expansion + Cx >= epsilon via coverage ADMM.
     */
    const int coverage_mode =
        (argc > 12) ? atoi(argv[12]) : 2;

    if (num_ranks <= 0) {
        die("num_rank_files must be positive");
    }

    if (coverage_mode < 0 || coverage_mode > 2) {
        die("coverage_mode must be 0, 1, or 2");
    }

    if (coverage_tau <= 0.0 ||
        admm_max_iter <= 0 ||
        admm_rho <= 0.0 ||
        admm_abs_tol <= 0.0 ||
        admm_rel_tol <= 0.0) {
        die("invalid interface-coverage ADMM option");
    }

    std::vector<Stage3RankData> data;
    data.reserve((size_t)num_ranks);
    std::unordered_set<RowKey,RowKeyHash> row_set;
    std::unordered_map<int64_t,double> candidate_seed_weight;

    for (int rank=0;rank<num_ranks;++rank) {
        Stage3RankData d=read_rank_file(directory,rank);
        for (const RowKey& k:d.row_key) row_set.insert(k);
        for (size_t i=0;i<d.selected_elem_id.size();++i) {
            const int64_t e=d.selected_elem_id[i];
            const double w=d.selected_weight[i];
            auto it=candidate_seed_weight.find(e);
            if (it==candidate_seed_weight.end()) candidate_seed_weight.emplace(e,w);
            else if (w>it->second) it->second=w;
        }
        data.push_back(std::move(d));
    }

    std::vector<RowKey> rows(row_set.begin(),row_set.end());
    std::sort(rows.begin(),rows.end(),[](const RowKey&a,const RowKey&b){
        if (a.snapshot!=b.snapshot) return a.snapshot<b.snapshot;
        if (a.owner_global_id!=b.owner_global_id) return a.owner_global_id<b.owner_global_id;
        return a.mode_offset<b.mode_offset;
    });

    const int stage2_selected_union_count =
        (int)candidate_seed_weight.size();

    std::vector<InterfaceGroup> interface_groups;

    /*
     * coverage_mode=0 is the true pre-coverage comparison baseline.
     * Do not even read/use the interface files in this mode.
     */
    if (coverage_mode > 0) {
        interface_groups =
            read_all_interface_groups(
                directory,
                num_ranks);

        if (interface_groups.empty()) {
            die("no physical interface groups were read");
        }
    }

    std::unordered_set<int64_t> coverage_added_element;
    int empty_interface_before_expansion = 0;

    /*
     * Candidate expansion belongs to modes 1 and 2 only.
     * In mode 0 the candidate pool remains exactly the Stage-2 union.
     */
    if (coverage_mode > 0) {
        /*
         * Compute a global column norm for every available physical element.
         * If Stage 2 removed every candidate from one interface, reintroduce
         * the strongest eligible column.
         */
        std::unordered_map<int64_t,double> all_column_norm_sq;

        for (const Stage3RankData& d : data) {
            const int nr = (int)d.row_key.size();
            const int nl = (int)d.local_elem_id.size();

            for (int lc = 0; lc < nl; ++lc) {
                double value = 0.0;

                for (int r = 0; r < nr; ++r) {
                    const double a =
                        d.local_matrix[
                            (size_t)r * (size_t)nl +
                            (size_t)lc];

                    value += a * a;
                }

                all_column_norm_sq[
                    d.local_elem_id[(size_t)lc]] += value;
            }
        }

        for (const InterfaceGroup& g : interface_groups) {
            bool covered = false;

            for (int64_t elem : g.member_original_id) {
                if (candidate_seed_weight.find(elem) !=
                    candidate_seed_weight.end()) {
                    covered = true;
                    break;
                }
            }

            if (covered) {
                continue;
            }

            empty_interface_before_expansion++;

            int64_t best_elem = -1;
            double best_norm = -1.0;

            for (int64_t elem : g.member_original_id) {
                const auto it =
                    all_column_norm_sq.find(elem);

                if (it == all_column_norm_sq.end()) {
                    continue;
                }

                if (it->second > best_norm) {
                    best_norm = it->second;
                    best_elem = elem;
                }
            }

            if (best_elem < 0) {
                fprintf(
                    stderr,
                    "ERROR: physical interface {%lld,%lld} has no "
                    "available Stage-3 element column\n",
                    (long long)g.key.source_global_id,
                    (long long)g.key.target_global_id);
                exit(EXIT_FAILURE);
            }

            candidate_seed_weight.emplace(best_elem, 0.0);
            coverage_added_element.insert(best_elem);
        }
    }

    std::vector<int64_t> candidates;
    candidates.reserve(candidate_seed_weight.size());

    for (const auto& kv : candidate_seed_weight) {
        candidates.push_back(kv.first);
    }

    std::sort(candidates.begin(), candidates.end());

    const int m = (int)rows.size();
    const int n = (int)candidates.size();

    if (m <= 0 || n <= 0) {
        die("empty Stage-3 system");
    }

    std::unordered_map<RowKey,int,RowKeyHash> row_index;
    row_index.reserve(rows.size()*2+1);
    for (int r=0;r<m;++r) row_index.emplace(rows[(size_t)r],r);
    std::unordered_map<int64_t,int> elem_index;
    elem_index.reserve(candidates.size()*2+1);
    for (int c=0;c<n;++c) elem_index.emplace(candidates[(size_t)c],c);

    std::vector<int> group_ptr;
    std::vector<int> group_item;
    std::vector<double> group_epsilon;

    group_ptr.reserve(interface_groups.size() + 1);
    group_epsilon.reserve(interface_groups.size());
    group_ptr.push_back(0);

    if (coverage_mode > 0) {
        for (const InterfaceGroup& g : interface_groups) {
            std::vector<int> local_item;

            for (int64_t elem : g.member_original_id) {
                const auto it = elem_index.find(elem);

                if (it != elem_index.end()) {
                    local_item.push_back(it->second);
                }
            }

            std::sort(local_item.begin(), local_item.end());
            local_item.erase(
                std::unique(
                    local_item.begin(),
                    local_item.end()),
                local_item.end());

            if (local_item.empty()) {
                fprintf(
                    stderr,
                    "ERROR: physical interface {%lld,%lld} has no "
                    "candidate after coverage expansion\n",
                    (long long)g.key.source_global_id,
                    (long long)g.key.target_global_id);
                exit(EXIT_FAILURE);
            }

            group_item.insert(
                group_item.end(),
                local_item.begin(),
                local_item.end());

            group_ptr.push_back((int)group_item.size());

            group_epsilon.push_back(
                2.0 * coverage_tau *
                (double)local_item.size());
        }
    }

    const int num_interface_group =
        (int)group_epsilon.size();

    double* Ablock=NULL;
    double** A=allocate_matrix(m,n,&Ablock);
    std::vector<double>b((size_t)m,0.0),x((size_t)n,0.0);
    std::vector<int> active((size_t)n,0);

    /* Sum BOTH RHS and columns over every rank: symmetric global assembly. */
    for (const Stage3RankData& d:data) {
        const int nr=(int)d.row_key.size();
        const int nl=(int)d.local_elem_id.size();
        for (int r=0;r<nr;++r) {
            const int gr=row_index.at(d.row_key[(size_t)r]);
            b[(size_t)gr]+=d.rhs[(size_t)r];
            for (int lc=0;lc<nl;++lc) {
                auto ei=elem_index.find(d.local_elem_id[(size_t)lc]);
                if (ei==elem_index.end()) continue;
                A[gr][ei->second]+=d.local_matrix[(size_t)r*(size_t)nl+(size_t)lc];
            }
        }

        int best=-1; double bestw=-1.0;
        for (int k=0;k<(int)d.selected_elem_id.size();++k) {
            if (d.selected_weight[(size_t)k]>bestw) { bestw=d.selected_weight[(size_t)k]; best=k; }
        }
        if (best>=0) {
            auto ei=elem_index.find(d.selected_elem_id[(size_t)best]);
            if (ei!=elem_index.end()) active[(size_t)ei->second]=1;
        }
    }

    /*
     * In mode 1 the added columns are merely available candidates; they are
     * NOT force-seeded into the active set.  This makes mode 1 a clean
     * "candidate pool only" ablation.
     *
     * Mode 2 uses them as warm-start columns for the constrained solve.
     */
    if (coverage_mode == 2) {
        for (int64_t elem : coverage_added_element) {
            const auto it = elem_index.find(elem);

            if (it != elem_index.end()) {
                active[(size_t)it->second] = 1;
            }
        }
    }

    double scale_sq=0.0;
    for (int r=0;r<m;++r) for (int c=0;c<n;++c) scale_sq+=A[r][c]*A[r][c];
    double scale=sqrt(scale_sq);
    if (scale<=DBL_MIN) scale=1.0;
    for (int r=0;r<m;++r) {
        b[(size_t)r]/=scale;
        for (int c=0;c<n;++c) A[r][c]/=scale;
    }

    /*
     * ------------------------------------------------------------
     * Coverage-anchor warm start
     *
     * The hard constraint is still Cx >= epsilon.  This block does NOT
     * replace that constraint.  It only ensures that at least one column
     * from every physical interface is presented to the inner sparse NNLS
     * at the beginning of every ADMM x-subproblem.
     *
     * This is important for tiny singleton constraints such as
     * epsilon = 2e-8: without an anchor, the sparse column-selection
     * heuristic can fail to discover that column before the inner
     * tolerance/iteration limit is reached.
     *
     * Choose the member having the largest norm in the normalized
     * Stage-3 matrix.  For a group whose objective columns are all zero,
     * fall back to its first member; the augmented coverage row still
     * gives that column a nonzero role in the ADMM subproblem.
     * ------------------------------------------------------------
     */
    std::vector<double> objective_col_norm_sq(
        (size_t)n,
        0.0);

    for (int r = 0; r < m; ++r) {
        for (int c = 0; c < n; ++c) {
            const double v = A[r][c];
            objective_col_norm_sq[(size_t)c] += v * v;
        }
    }

    std::vector<int> coverage_anchor(
        (size_t)num_interface_group,
        -1);

    std::unordered_set<int> unique_coverage_anchor;
    unique_coverage_anchor.reserve(
        (size_t)num_interface_group * 2 + 1);

    int singleton_interface_groups = 0;

    if (coverage_mode == 2) {
        for (int g = 0; g < num_interface_group; ++g) {
            const int kS = group_ptr[(size_t)g];
            const int kE = group_ptr[(size_t)g + 1];

            if (kE - kS == 1) {
                singleton_interface_groups++;
            }

            int best_c = -1;
            double best_norm_sq = -1.0;

            for (int k = kS; k < kE; ++k) {
                const int c = group_item[(size_t)k];

                if (c < 0 || c >= n) {
                    fprintf(
                        stderr,
                        "ERROR: invalid coverage group column: "
                        "g=%d c=%d n=%d\n",
                        g,
                        c,
                        n);
                    exit(EXIT_FAILURE);
                }

                const double cn =
                    objective_col_norm_sq[(size_t)c];

                if (best_c < 0 || cn > best_norm_sq) {
                    best_c = c;
                    best_norm_sq = cn;
                }
            }

            if (best_c < 0) {
                fprintf(
                    stderr,
                    "ERROR: interface group %d has no coverage anchor\n",
                    g);
                exit(EXIT_FAILURE);
            }

            coverage_anchor[(size_t)g] = best_c;
            active[(size_t)best_c] = 1;
            unique_coverage_anchor.insert(best_c);
        }
    }

    int nactive=0;
    for(int v:active) if(v)++nactive;

    const char* coverage_mode_name =
        (coverage_mode == 0)
        ? "baseline: Stage-2 candidates only, unconstrained NNLS"
        : (coverage_mode == 1)
          ? "candidate-expanded ablation, unconstrained NNLS"
          : "proposed: candidate expansion + hard interface coverage";

    printf("============================================================\n");
    printf("[Stage 3 v14 selectable coverage mode]\n");
    printf("coverage mode                      = %d\n", coverage_mode);
    printf("coverage mode description          = %s\n", coverage_mode_name);
    printf("rank files                         = %d\n", num_ranks);
    printf("global keyed rows                  = %d\n", m);
    printf("Stage-2 selected union elements    = %d\n",
        stage2_selected_union_count);
    printf("coverage-added candidate elements  = %zu\n",
        coverage_added_element.size());
    printf("Stage-3 candidate pool             = %d\n", n);
    printf("directed interface constraints     = %d\n",
        num_interface_group);
    printf("interfaces empty before expansion  = %d\n",
        empty_interface_before_expansion);
    printf("initial active elements            = %d\n", nactive);
    printf("unique coverage anchors            = %zu\n",
        unique_coverage_anchor.size());
    printf("singleton interface constraints    = %d\n",
        singleton_interface_groups);
    printf("inner NNLS max iterations          = %d\n", max_iter);
    printf("inner NNLS tolerance               = %.6e\n", tol);
    printf("coverage selection threshold tau   = %.6e\n",
        coverage_tau);
    printf("ADMM max iterations                = %d\n",
        admm_max_iter);
    printf("ADMM initial rho                   = %.6e\n",
        admm_rho);
    printf("global normalization scale         = %.15e\n",
        scale);
    printf("============================================================\n");

    double residual = 0.0;
    double max_coverage_violation = 0.0;
    int admm_outer_iter = 0;
    int coverage_status = 0;

    if (coverage_mode == 2) {
        coverage_status =
            monolis_optimize_nnls_R_with_sparse_solution_interface_coverage_admm(
                A,
                b.data(),
                x.data(),
                m,
                n,
                max_iter,
                tol,
                active.data(),
                num_interface_group,
                group_ptr.data(),
                group_item.data(),
                (int)group_item.size(),
                group_epsilon.data(),
                admm_max_iter,
                admm_rho,
                admm_abs_tol,
                admm_rel_tol,
                admm_verbose,
                &residual,
                &max_coverage_violation,
                &admm_outer_iter);
    }
    else {
        monolis_optimize_nnls_R_with_sparse_solution_initial_set(
            A,
            b.data(),
            x.data(),
            m,
            n,
            max_iter,
            tol,
            active.data(),
            &residual);
    }

    int uncovered_after_solve = 0;
    double minimum_group_sum = 0.0;
    double minimum_group_margin = 0.0;

    int minimum_margin_group = -1;
    int minimum_sum_group = -1;

    if (coverage_mode == 2) {
        minimum_group_sum = DBL_MAX;
        minimum_group_margin = DBL_MAX;

        for (int g = 0; g < num_interface_group; ++g) {
            double sum = 0.0;
            bool selected_member = false;

            for (int k = group_ptr[(size_t)g];
                 k < group_ptr[(size_t)g + 1];
                 ++k) {
                const int c = group_item[(size_t)k];
                sum += x[(size_t)c];

                if (x[(size_t)c] > weight_tol) {
                    selected_member = true;
                }
            }

            if (sum < minimum_group_sum) {
                minimum_group_sum = sum;
                minimum_sum_group = g;
            }

            const double margin =
                sum - group_epsilon[(size_t)g];

            if (margin < minimum_group_margin) {
                minimum_group_margin = margin;
                minimum_margin_group = g;
            }

            if (!selected_member) {
                uncovered_after_solve++;
            }
        }

        printf(
            "[COVERAGE-RESULT] status=%d outer_iter=%d "
            "max_violation=%.15e min_group_sum=%.15e "
            "min_margin=%.15e uncovered_after_threshold=%d\n",
            coverage_status,
            admm_outer_iter,
            max_coverage_violation,
            minimum_group_sum,
            minimum_group_margin,
            uncovered_after_solve);

        int printed_bad_group = 0;

        for (int g = 0; g < num_interface_group; ++g) {
            double sum = 0.0;
            double max_member_x = 0.0;
            int max_member_c = -1;
            int selected_count = 0;

            for (int k = group_ptr[(size_t)g];
                 k < group_ptr[(size_t)g + 1];
                 ++k) {

                const int c =
                    group_item[(size_t)k];

                const double xc =
                    x[(size_t)c];

                sum += xc;

                if (xc > max_member_x ||
                    max_member_c < 0) {

                    max_member_x = xc;
                    max_member_c = c;
                }

                if (xc > weight_tol) {
                    selected_count++;
                }
            }

            const double margin =
                sum - group_epsilon[(size_t)g];

            if (margin < -admm_abs_tol ||
                selected_count == 0) {

                const InterfaceGroup& ig =
                    interface_groups[(size_t)g];

                const int anchor_c =
                    coverage_anchor[(size_t)g];

                printf(
                    "[COVERAGE-BAD-GROUP] g=%d interface={%lld,%lld} "
                    "members=%d epsilon=%.15e sum=%.15e margin=%.15e "
                    "selected=%d max_x=%.15e "
                    "anchor_col=%d anchor_elem=%lld anchor_x=%.15e\n",
                    g,
                    (long long)ig.key.source_global_id,
                    (long long)ig.key.target_global_id,
                    group_ptr[(size_t)g + 1] -
                        group_ptr[(size_t)g],
                    group_epsilon[(size_t)g],
                    sum,
                    margin,
                    selected_count,
                    max_member_x,
                    anchor_c,
                    (anchor_c >= 0)
                        ? (long long)candidates[(size_t)anchor_c]
                        : -1LL,
                    (anchor_c >= 0)
                        ? x[(size_t)anchor_c]
                        : 0.0);

                if (++printed_bad_group >= 20) {
                    break;
                }
            }
        }

        if (minimum_margin_group >= 0) {
            const InterfaceGroup& ig =
                interface_groups[
                    (size_t)minimum_margin_group];

            printf(
                "[COVERAGE-WORST] g=%d interface={%lld,%lld} "
                "members=%d epsilon=%.15e margin=%.15e\n",
                minimum_margin_group,
                (long long)ig.key.source_global_id,
                (long long)ig.key.target_global_id,
                group_ptr[(size_t)minimum_margin_group + 1] -
                    group_ptr[(size_t)minimum_margin_group],
                group_epsilon[(size_t)minimum_margin_group],
                minimum_group_margin);
        }

        if (coverage_status != 0 ||
            uncovered_after_solve != 0 ||
            max_coverage_violation >
                10.0 * std::max(
                    admm_abs_tol,
                    admm_rel_tol *
                        *std::max_element(
                            group_epsilon.begin(),
                            group_epsilon.end()))) {
            die("interface-coverage constrained NNLS did not satisfy "
                "the requested hard constraints");
        }
    }
    else {
        printf(
            "[COVERAGE-RESULT] skipped (coverage_mode=%d: %s)\n",
            coverage_mode,
            coverage_mode_name);
    }

    std::vector<std::pair<int64_t,double>> selected;

    for (int c = 0; c < n; ++c) {
        if (x[(size_t)c] > weight_tol) {
            selected.emplace_back(
                candidates[(size_t)c],
                x[(size_t)c]);
        }
    }

    const double raw_residual = residual * scale;

    /*
     * Read the exact Stage-1/2 route table before printing reduction summaries.
     * The unique original-global IDs in this table are the physical elements
     * before Stage-1 reduction (ghost/rank copies are removed by the set).
     */
    const std::vector<RouteEntry> routes=read_all_routes(directory,num_ranks);
    std::unordered_set<int64_t> original_physical_elem_set;
    original_physical_elem_set.reserve(routes.size()*2+1);
    for (const RouteEntry& e : routes) {
        original_physical_elem_set.insert(e.original_global_id);
    }

    const int cumulative_num_original =
        (int)original_physical_elem_set.size();

    /*
     * Recover the exact Stage-1 selected physical-element union.
     * stage1_stats.* stores legacy parallel element IDs per fine subdomain;
     * stage3_route.* converts each exact occurrence to original-global ID.
     */
    const std::vector<Stage1SubdomainStats> stage1_stats =
        read_all_stage1_stats(directory, num_ranks);

    std::unordered_map<
        Stage1RouteKey,
        int64_t,
        Stage1RouteKeyHash> stage1_route_lookup;

    stage1_route_lookup.reserve(
        routes.size() * 2 + 1);

    for (const RouteEntry& e : routes) {
        Stage1RouteKey key;
        key.rank = e.rank;
        key.fine_subdomain = e.fine_subdomain;
        key.parallel_elem_id = e.parallel_elem_id;

        const auto old =
            stage1_route_lookup.find(key);

        if (old == stage1_route_lookup.end()) {
            stage1_route_lookup.emplace(
                key,
                e.original_global_id);
        }
        else if (old->second != e.original_global_id) {
            die("ambiguous exact Stage-1 route");
        }
    }

    std::unordered_set<int64_t> stage1_selected_original_set;
    size_t stage1_sum_local_candidates = 0;
    size_t stage1_sum_local_unique = 0;

    for (const Stage1SubdomainStats& st : stage1_stats) {
        stage1_sum_local_candidates +=
            (size_t)st.local_candidates;
        stage1_sum_local_unique +=
            (size_t)st.local_unique;

        for (int64_t parallel_id :
             st.candidate_parallel_ids) {

            Stage1RouteKey key;
            key.rank = st.rank;
            key.fine_subdomain =
                st.fine_subdomain;
            key.parallel_elem_id =
                parallel_id;

            const auto it =
                stage1_route_lookup.find(key);

            if (it == stage1_route_lookup.end()) {
                fprintf(
                    stderr,
                    "ERROR: Stage-1 candidate has no exact route: "
                    "rank=%d subdomain=%d parallel_elem=%lld\n",
                    st.rank,
                    st.fine_subdomain,
                    (long long)parallel_id);
                exit(EXIT_FAILURE);
            }

            stage1_selected_original_set.insert(
                it->second);
        }
    }

    const int stage1_num_selected_union =
        (int)stage1_selected_original_set.size();

    const int stage1_num_removed_global =
        cumulative_num_original -
        stage1_num_selected_union;

    const double stage1_retained_global_percent =
        (cumulative_num_original > 0)
        ? 100.0
            * (double)stage1_num_selected_union
            / (double)cumulative_num_original
        : 0.0;

    const double stage1_reduction_global_percent =
        (cumulative_num_original > 0)
        ? 100.0
            * (double)stage1_num_removed_global
            / (double)cumulative_num_original
        : 0.0;

    const int stage2_num_selected_union =
        stage2_selected_union_count;
    const int stage3_num_candidates = n;
    const int stage3_num_selected = (int)selected.size();
    const int stage3_num_removed =
        stage3_num_candidates - stage3_num_selected;

    const double stage3_retained_percent =
        (stage3_num_candidates > 0)
        ? 100.0 * (double)stage3_num_selected
            / (double)stage3_num_candidates
        : 0.0;

    const double stage3_reduction_percent =
        (stage3_num_candidates > 0)
        ? 100.0 * (double)stage3_num_removed
            / (double)stage3_num_candidates
        : 0.0;

    printf("\n");
    printf("============================================================\n");
    printf("[Stage 1 per-fine-subdomain reduction]\n");
    printf("subdomain  rank  original  selected  removed  retained(%%)  reduction(%%)  residual\n");

    for (const Stage1SubdomainStats& st : stage1_stats) {
        const int removed =
            st.local_unique - st.local_candidates;

        const double retained =
            (st.local_unique > 0)
            ? 100.0
                * (double)st.local_candidates
                / (double)st.local_unique
            : 0.0;

        const double reduction =
            (st.local_unique > 0)
            ? 100.0
                * (double)removed
                / (double)st.local_unique
            : 0.0;

        printf(
            "%9d  %4d  %8d  %8d  %7d  %11.3f  %12.3f  %.6e\n",
            st.fine_subdomain,
            st.rank,
            st.local_unique,
            st.local_candidates,
            removed,
            retained,
            reduction,
            st.residual);
    }

    printf("------------------------------------------------------------\n");
    printf("sum local original counts          = %zu\n",
        stage1_sum_local_unique);
    printf("sum local Stage-1 selected counts  = %zu\n",
        stage1_sum_local_candidates);
    printf("global original physical elements  = %d\n",
        cumulative_num_original);
    printf("global Stage-1 selected union       = %d\n",
        stage1_num_selected_union);
    printf("global Stage-1 removed              = %d\n",
        stage1_num_removed_global);
    printf("global Stage-1 retained             = %.3f %%\n",
        stage1_retained_global_percent);
    printf("global Stage-1 reduction            = %.3f %%\n",
        stage1_reduction_global_percent);
    printf("============================================================\n");

    printf("\n");
    printf("============================================================\n");
    printf("[Stage 3 reduction summary]\n");
    printf("Stage-3 candidate elements        = %d\n",
        stage3_num_candidates);
    printf("Stage-3 selected elements         = %d\n",
        stage3_num_selected);
    printf("Stage-3 removed elements          = %d\n",
        stage3_num_removed);
    printf("Stage-3 retained                  = %.3f %%\n",
        stage3_retained_percent);
    printf("Stage-3 reduction                 = %.3f %%\n",
        stage3_reduction_percent);
    printf("normalized residual               = %.15e\n",
        residual);
    printf("raw residual                      = %.15e\n",
        raw_residual);
    printf("============================================================\n");

    const int cumulative_num_selected = stage3_num_selected;
    const int cumulative_num_removed =
        cumulative_num_original - cumulative_num_selected;

    const double cumulative_retained_percent =
        (cumulative_num_original > 0)
        ? 100.0 * (double)cumulative_num_selected
            / (double)cumulative_num_original
        : 0.0;

    const double cumulative_reduction_percent =
        (cumulative_num_original > 0)
        ? 100.0 * (double)cumulative_num_removed
            / (double)cumulative_num_original
        : 0.0;

    const double stage2_union_retained_percent =
        (cumulative_num_original > 0)
        ? 100.0 * (double)stage2_num_selected_union
            / (double)cumulative_num_original
        : 0.0;

    const double stage2_union_reduction_percent =
        (cumulative_num_original > 0)
        ? 100.0 * (double)(cumulative_num_original-stage2_num_selected_union)
            / (double)cumulative_num_original
        : 0.0;

    const int stage2_removed_from_stage1 =
        stage1_num_selected_union -
        stage2_num_selected_union;

    const double stage2_retained_vs_stage1_percent =
        (stage1_num_selected_union > 0)
        ? 100.0
            * (double)stage2_num_selected_union
            / (double)stage1_num_selected_union
        : 0.0;

    const double stage2_reduction_vs_stage1_percent =
        (stage1_num_selected_union > 0)
        ? 100.0
            * (double)stage2_removed_from_stage1
            / (double)stage1_num_selected_union
        : 0.0;

    printf("\n");
    printf("============================================================\n");
    printf("[Stage 1 -> Stage 2 -> Stage 3 cumulative reduction]\n");
    printf("original physical elements        = %d\n",
        cumulative_num_original);
    printf("Stage-1 selected union elements   = %d\n",
        stage1_num_selected_union);
    printf("Stage-1 retained vs. original     = %.3f %%\n",
        stage1_retained_global_percent);
    printf("Stage-1 reduction vs. original    = %.3f %%\n",
        stage1_reduction_global_percent);
    printf("Stage-2 selected union elements   = %d\n",
        stage2_num_selected_union);
    printf("coverage-added Stage-3 candidates= %zu\n",
        coverage_added_element.size());
    printf("Stage-3 candidate pool            = %d\n",
        stage3_num_candidates);
    printf("Stage-2 removed from Stage-1      = %d\n",
        stage2_removed_from_stage1);
    printf("Stage-2 retained vs. Stage-1      = %.3f %%\n",
        stage2_retained_vs_stage1_percent);
    printf("Stage-2 reduction vs. Stage-1     = %.3f %%\n",
        stage2_reduction_vs_stage1_percent);
    printf("Stage-2 retained vs. original     = %.3f %%\n",
        stage2_union_retained_percent);
    printf("Stage-2 reduction vs. original    = %.3f %%\n",
        stage2_union_reduction_percent);
    printf("final Stage-3 selected elements   = %d\n",
        cumulative_num_selected);
    printf("Stage-3 retained vs. candidate pool= %.3f %%\n",
        stage3_retained_percent);
    printf("Stage-3 reduction vs. candidate pool= %.3f %%\n",
        stage3_reduction_percent);
    printf("cumulative removed elements       = %d\n",
        cumulative_num_removed);
    printf("cumulative retained               = %.3f %%\n",
        cumulative_retained_percent);
    printf("cumulative reduction              = %.3f %%\n",
        cumulative_reduction_percent);
    printf("============================================================\n");

    /*
     * Save the same reduction information in one text file for later
     * plotting, paper tables, and regression comparisons.
     */
    {
        const std::string report_path =
            join_path(
                directory,
                "DDECM/stage_reduction_summary.txt");

        FILE* fp_report =
            fopen(report_path.c_str(), "w");

        if (fp_report == NULL) {
            die("cannot write stage_reduction_summary.txt");
        }

        fprintf(
            fp_report,
            "coverage_mode %d\n",
            coverage_mode);
        fprintf(
            fp_report,
            "coverage_mode_description %s\n\n",
            coverage_mode_name);

        fprintf(
            fp_report,
            "[Stage 1 per-fine-subdomain reduction]\n");
        fprintf(
            fp_report,
            "subdomain rank original selected removed retained_percent reduction_percent residual\n");

        for (const Stage1SubdomainStats& st : stage1_stats) {
            const int removed =
                st.local_unique -
                st.local_candidates;

            const double retained =
                (st.local_unique > 0)
                ? 100.0
                    * (double)st.local_candidates
                    / (double)st.local_unique
                : 0.0;

            const double reduction =
                (st.local_unique > 0)
                ? 100.0
                    * (double)removed
                    / (double)st.local_unique
                : 0.0;

            fprintf(
                fp_report,
                "%d %d %d %d %d %.15e %.15e %.15e\n",
                st.fine_subdomain,
                st.rank,
                st.local_unique,
                st.local_candidates,
                removed,
                retained,
                reduction,
                st.residual);
        }

        fprintf(fp_report, "\n[Global reduction]\n");
        fprintf(
            fp_report,
            "original_physical_elements %d\n",
            cumulative_num_original);
        fprintf(
            fp_report,
            "stage1_selected_union %d\n",
            stage1_num_selected_union);
        fprintf(
            fp_report,
            "stage1_retained_percent %.15e\n",
            stage1_retained_global_percent);
        fprintf(
            fp_report,
            "stage1_reduction_percent %.15e\n",
            stage1_reduction_global_percent);
        fprintf(
            fp_report,
            "stage2_selected_union %d\n",
            stage2_num_selected_union);
        fprintf(
            fp_report,
            "coverage_added_candidates %zu\n",
            coverage_added_element.size());
        fprintf(
            fp_report,
            "stage3_candidate_pool %d\n",
            stage3_num_candidates);
        fprintf(
            fp_report,
            "stage2_retained_vs_stage1_percent %.15e\n",
            stage2_retained_vs_stage1_percent);
        fprintf(
            fp_report,
            "stage2_reduction_vs_stage1_percent %.15e\n",
            stage2_reduction_vs_stage1_percent);
        fprintf(
            fp_report,
            "stage2_retained_vs_original_percent %.15e\n",
            stage2_union_retained_percent);
        fprintf(
            fp_report,
            "stage2_reduction_vs_original_percent %.15e\n",
            stage2_union_reduction_percent);
        fprintf(
            fp_report,
            "stage3_selected %d\n",
            stage3_num_selected);
        fprintf(
            fp_report,
            "stage3_retained_vs_stage2_percent %.15e\n",
            stage3_retained_percent);
        fprintf(
            fp_report,
            "stage3_reduction_vs_candidate_pool_percent %.15e\n",
            stage3_reduction_percent);
        fprintf(
            fp_report,
            "coverage_directed_groups %d\n",
            num_interface_group);
        fprintf(
            fp_report,
            "coverage_max_violation %.30e\n",
            max_coverage_violation);
        fprintf(
            fp_report,
            "coverage_min_group_sum %.30e\n",
            minimum_group_sum);
        fprintf(
            fp_report,
            "coverage_uncovered_after_threshold %d\n",
            uncovered_after_solve);
        fprintf(
            fp_report,
            "coverage_admm_status %d\n",
            coverage_status);
        fprintf(
            fp_report,
            "coverage_admm_outer_iterations %d\n",
            admm_outer_iter);
        fprintf(
            fp_report,
            "cumulative_retained_percent %.15e\n",
            cumulative_retained_percent);
        fprintf(
            fp_report,
            "cumulative_reduction_percent %.15e\n",
            cumulative_reduction_percent);
        fprintf(
            fp_report,
            "normalized_residual %.30e\n",
            residual);
        fprintf(
            fp_report,
            "raw_residual %.30e\n",
            raw_residual);

        fclose(fp_report);

        printf(
            "[report] wrote %s\n",
            report_path.c_str());
    }

    printf("[Stage 3 v14 result] mode=%d selected=%d / %d "
           "reduction=%.3f%% normalized_residual=%.15e "
           "raw_residual=%.15e\n",
        coverage_mode,
        stage3_num_selected,
        stage3_num_candidates,
        stage3_reduction_percent,
        residual,
        raw_residual);

    {
        const std::string path=join_path(directory,"DDECM/stage3_selected_elem.txt");
        FILE* fp=fopen(path.c_str(),"w");
        if(fp==NULL) die("cannot write stage3_selected_elem.txt");
        fprintf(fp,"%zu\n",selected.size());
        for(const auto&p:selected) fprintf(fp,"%lld %.30e\n",(long long)p.first,p.second);
        fclose(fp);
    }
    {
        const std::string path=join_path(directory,"DDECM/stage3_NNLS_residual.dat");
        FILE* fp=fopen(path.c_str(),"w");
        if(fp==NULL) die("cannot write stage3_NNLS_residual.dat");
        fprintf(fp,"%.30e\n",raw_residual);
        fclose(fp);
    }

    write_parallel_legacy_files(directory,routes,selected);

    free(Ablock);
    free(A);
    return EXIT_SUCCESS;
}
