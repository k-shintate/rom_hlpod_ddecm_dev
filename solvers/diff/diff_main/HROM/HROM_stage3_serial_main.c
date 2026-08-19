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
#include <chrono>
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

using HROMClock = std::chrono::steady_clock;

static double HROM_elapsed_sec(
    const HROMClock::time_point& a,
    const HROMClock::time_point& b)
{
    return std::chrono::duration<double>(b - a).count();
}

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

/*
 * Mode-4 internal operator group.
 *
 * member_original_id contains strictly internal owned elements of one fine
 * subdomain.  Elements that belong to any physical interface group are removed
 * so that an internal operator constraint cannot be satisfied only by an
 * interface-supporting element.
 */
struct InternalOperatorGroup {
    int fine_subdomain = -1;
    std::vector<int64_t> member_original_id;
};

/*
 * Mode-5 adaptive enrichment diagnostics.
 *
 * Mode 5 starts from the exact mode-0 Stage-3 candidate set (the Stage-2
 * selected union), solves an unconstrained sparse NNLS with a deliberately
 * coarse seed tolerance, measures operator directions that are not represented
 * by the current output support, and reintroduces only the highest-priority
 * missing columns.  Hard singleton coverage is optional and is used only as a
 * final fallback when candidate enrichment alone cannot retain a required
 * direction.
 */
struct Mode5AdaptiveResult {
    int rounds = 0;
    int added_candidates = 0;
    int deficient_groups_final = 0;
    double final_max_defect = 0.0;
    std::vector<int64_t> added_element;
    std::vector<int64_t> hard_anchor_element;
};

struct Mode5System {
    int m = 0;
    int n = 0;
    std::vector<int64_t> candidates;
    std::vector<double> A; /* row-major, normalized */
    std::vector<double> b; /* normalized */
    std::vector<double> x;
    double scale = 1.0;
    double residual = 0.0;
};

struct Mode5OperatorGroupDesc {
    int kind = 0; /* 0=interface, 1=internal */
    int64_t owner_a = -1;
    int64_t owner_b = -1;
    std::vector<int> rows;
    std::vector<int64_t> members;
};

struct Mode5DefectCandidate {
    double defect = 0.0;
    double priority = 0.0;
    double importance = 0.0;
    double novelty = 0.0;
    double residual_alignment = 0.0;
    int group_index = -1;
    int64_t elem = -1;
    bool already_candidate = false;
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


/*
 * Build strictly-internal fine-subdomain element groups for coverage_mode=4.
 *
 * stage3_route.* contains exact routes for the owned/internal elements of each
 * fine subdomain.  Remove every element that appears in any physical interface
 * group; the remaining members are the true interior set used for internal
 * operator-signature coverage.
 */
static std::vector<InternalOperatorGroup> build_internal_operator_groups(
    const std::vector<RouteEntry>& routes,
    const std::vector<InterfaceGroup>& interface_groups)
{
    std::unordered_set<int64_t> interface_member;
    for (const InterfaceGroup& g : interface_groups) {
        for (int64_t elem : g.member_original_id) {
            interface_member.insert(elem);
        }
    }

    std::unordered_map<int,std::unordered_set<int64_t>> by_subdomain;
    for (const RouteEntry& e : routes) {
        if (interface_member.find(e.original_global_id) != interface_member.end()) {
            continue;
        }
        by_subdomain[e.fine_subdomain].insert(e.original_global_id);
    }

    std::vector<InternalOperatorGroup> groups;
    groups.reserve(by_subdomain.size());

    size_t total_membership = 0;
    for (const auto& kv : by_subdomain) {
        InternalOperatorGroup g;
        g.fine_subdomain = kv.first;
        g.member_original_id.assign(kv.second.begin(),kv.second.end());
        std::sort(g.member_original_id.begin(),g.member_original_id.end());
        total_membership += g.member_original_id.size();
        groups.push_back(std::move(g));
    }

    std::sort(
        groups.begin(),groups.end(),
        [](const InternalOperatorGroup& a,const InternalOperatorGroup& b) {
            return a.fine_subdomain < b.fine_subdomain;
        });

    printf(
        "[internal-groups] fine_subdomains=%zu strictly_internal_memberships=%zu "
        "excluded_interface_elements=%zu\n",
        groups.size(),total_membership,interface_member.size());

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


struct ThresholdRefitResult {
    int enabled = 0;
    int passes = 0;
    int initial_support = 0;
    int final_support = 0;
    int stable = 1;
    double residual_before = 0.0;
    double residual_after_initial_cut = 0.0;
    double residual_final = 0.0;
    int coverage_status = 0;
    double max_coverage_violation = 0.0;
    int coverage_outer_iter = 0;
};

static double compute_stage3_residual(
    double** A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    int m,
    int n)
{
    long double r2 = 0.0L;
    for (int r = 0; r < m; ++r) {
        long double ax = 0.0L;
        for (int c = 0; c < n; ++c) {
            ax += (long double)A[r][c] * (long double)x[(size_t)c];
        }
        const long double d = ax - (long double)b[(size_t)r];
        r2 += d*d;
    }
    return std::sqrt((double)r2);
}

static std::vector<int> threshold_support(
    std::vector<double>& x,
    double weight_tol)
{
    std::vector<int> support;
    support.reserve(x.size());
    for (int c = 0; c < (int)x.size(); ++c) {
        if (x[(size_t)c] > weight_tol) {
            support.push_back(c);
        }
        else {
            x[(size_t)c] = 0.0;
        }
    }
    return support;
}

/*
 * Threshold-consistent restricted refit.
 *
 * Keep the continuous NNLS itself unchanged.  The online output threshold
 * first defines S={c | x_c>weight_tol}; then the same least-squares problem
 * is solved again on columns S only.  In coverage modes, every physical and
 * operator-aware coverage row is restricted to S and enforced again with the
 * original epsilon.  If refitting creates new coefficients <=weight_tol,
 * those columns are removed and the restricted problem is solved again.
 * Since columns can only disappear, the support shrinks monotonically.
 *
 * The final residual is therefore evaluated for exactly the coefficient
 * vector written to the online HROM files, without replacing the convex NNLS
 * by an implicit nonconvex condition x_i=0 or x_i>weight_tol.
 */
static ThresholdRefitResult threshold_consistent_refit(
    double** A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int m,
    int n,
    int max_iter,
    double tol,
    double weight_tol,
    int coverage_mode,
    const std::vector<int>& solver_group_ptr,
    const std::vector<int>& solver_group_item,
    const std::vector<double>& solver_group_epsilon,
    int admm_max_iter,
    double admm_rho,
    double admm_abs_tol,
    double admm_rel_tol,
    int admm_verbose,
    int max_pass)
{
    ThresholdRefitResult out;
    out.enabled = 1;
    out.residual_before = compute_stage3_residual(A,b,x,m,n);

    std::vector<int> support = threshold_support(x,weight_tol);
    out.initial_support = (int)support.size();
    out.residual_after_initial_cut = compute_stage3_residual(A,b,x,m,n);

    printf(
        "[THRESHOLD-REFIT-START] weight_tol=%.15e support=%d/%d "
        "residual_before=%.15e residual_after_cut=%.15e cut_delta=%.15e\n",
        weight_tol,
        out.initial_support,
        n,
        out.residual_before,
        out.residual_after_initial_cut,
        out.residual_after_initial_cut - out.residual_before);

    if (support.empty()) {
        die("threshold refit removed every Stage-3 candidate");
    }

    if (max_pass <= 0) max_pass = n;
    out.stable = 0;

    for (int pass = 1; pass <= max_pass; ++pass) {
        const int ns = (int)support.size();
        double* Asub_block = NULL;
        double** Asub = allocate_matrix(m,ns,&Asub_block);
        std::vector<double> xsub((size_t)ns,0.0);
        std::vector<int> active_sub((size_t)ns,1);
        std::vector<int> old_to_sub((size_t)n,-1);

        for (int j = 0; j < ns; ++j) {
            old_to_sub[(size_t)support[(size_t)j]] = j;
        }
        for (int r = 0; r < m; ++r) {
            for (int j = 0; j < ns; ++j) {
                Asub[r][j] = A[r][support[(size_t)j]];
            }
        }

        double refit_residual = 0.0;
        int refit_status = 0;
        double refit_max_violation = 0.0;
        int refit_outer_iter = 0;

        if (coverage_mode >= 2) {
            const int ng = (int)solver_group_epsilon.size();
            std::vector<int> sub_group_ptr;
            std::vector<int> sub_group_item;
            std::vector<double> sub_group_epsilon = solver_group_epsilon;
            sub_group_ptr.reserve((size_t)ng + 1);
            sub_group_ptr.push_back(0);

            for (int g = 0; g < ng; ++g) {
                std::vector<int> local;
                for (int k = solver_group_ptr[(size_t)g];
                     k < solver_group_ptr[(size_t)g + 1]; ++k) {
                    const int old_col = solver_group_item[(size_t)k];
                    const int sub_col = old_to_sub[(size_t)old_col];
                    if (sub_col >= 0) local.push_back(sub_col);
                }
                std::sort(local.begin(),local.end());
                local.erase(std::unique(local.begin(),local.end()),local.end());
                if (local.empty()) {
                    fprintf(
                        stderr,
                        "ERROR: threshold refit emptied coverage group %d at pass %d; "
                        "weight_tol=%.15e\n",
                        g,pass,weight_tol);
                    free(Asub);
                    free(Asub_block);
                    exit(EXIT_FAILURE);
                }
                sub_group_item.insert(
                    sub_group_item.end(),local.begin(),local.end());
                sub_group_ptr.push_back((int)sub_group_item.size());
            }

            refit_status =
                monolis_optimize_nnls_R_with_sparse_solution_interface_coverage_admm(
                    Asub,
                    b.data(),
                    xsub.data(),
                    m,
                    ns,
                    max_iter,
                    tol,
                    active_sub.data(),
                    ng,
                    sub_group_ptr.data(),
                    sub_group_item.data(),
                    (int)sub_group_item.size(),
                    sub_group_epsilon.data(),
                    admm_max_iter,
                    admm_rho,
                    admm_abs_tol,
                    admm_rel_tol,
                    admm_verbose,
                    &refit_residual,
                    &refit_max_violation,
                    &refit_outer_iter);

            if (refit_status != 0) {
                fprintf(
                    stderr,
                    "ERROR: threshold-restricted coverage refit failed: "
                    "pass=%d status=%d max_violation=%.15e\n",
                    pass,refit_status,refit_max_violation);
                free(Asub);
                free(Asub_block);
                exit(EXIT_FAILURE);
            }
        }
        else {
            monolis_optimize_nnls_R_with_sparse_solution_initial_set(
                Asub,
                const_cast<double*>(b.data()),
                xsub.data(),
                m,
                ns,
                max_iter,
                tol,
                active_sub.data(),
                &refit_residual);
        }

        std::fill(x.begin(),x.end(),0.0);
        for (int j = 0; j < ns; ++j) {
            x[(size_t)support[(size_t)j]] = xsub[(size_t)j];
        }

        std::vector<int> new_support = threshold_support(x,weight_tol);
        const double residual_after_cut = compute_stage3_residual(A,b,x,m,n);
        const int removed = ns - (int)new_support.size();

        printf(
            "[THRESHOLD-REFIT-PASS] pass=%d restricted_cols=%d kept=%zu "
            "removed=%d solver_residual=%.15e residual_after_cut=%.15e "
            "coverage_status=%d max_violation=%.15e\n",
            pass,
            ns,
            new_support.size(),
            removed,
            refit_residual,
            residual_after_cut,
            refit_status,
            refit_max_violation);

        out.passes = pass;
        out.coverage_status = refit_status;
        out.max_coverage_violation = refit_max_violation;
        out.coverage_outer_iter = refit_outer_iter;

        free(Asub);
        free(Asub_block);

        if (new_support.size() == support.size()) {
            out.stable = 1;
            support.swap(new_support);
            break;
        }

        if (new_support.empty()) {
            die("threshold refit produced an empty support");
        }
        support.swap(new_support);
    }

    if (!out.stable) {
        die("threshold refit support did not stabilize; increase threshold_refit_max_pass");
    }

    out.final_support = (int)support.size();
    out.residual_final = compute_stage3_residual(A,b,x,m,n);

    int small_positive = 0;
    for (double v : x) {
        if (v > 0.0 && v <= weight_tol) small_positive++;
    }
    if (small_positive != 0) {
        die("threshold refit ended with positive coefficients below weight_tol");
    }

    printf(
        "[THRESHOLD-REFIT-DONE] passes=%d support=%d/%d stable=%d "
        "residual_before=%.15e residual_final=%.15e small_positive=%d\n",
        out.passes,
        out.final_support,
        n,
        out.stable,
        out.residual_before,
        out.residual_final,
        small_positive);

    return out;
}


static Mode5System mode5_solve_unconstrained_system(
    const std::vector<Stage3RankData>& data,
    const std::vector<RowKey>& rows,
    const std::unordered_map<int64_t,double>& candidate_seed_weight,
    int max_iter,
    double solve_tol)
{
    Mode5System sys;
    sys.m = (int)rows.size();
    sys.candidates.reserve(candidate_seed_weight.size());
    for (const auto& kv : candidate_seed_weight) {
        sys.candidates.push_back(kv.first);
    }
    std::sort(sys.candidates.begin(),sys.candidates.end());
    sys.n = (int)sys.candidates.size();

    if (sys.m <= 0 || sys.n <= 0) {
        die("mode-5 temporary system is empty");
    }

    std::unordered_map<RowKey,int,RowKeyHash> row_index;
    row_index.reserve(rows.size()*2 + 1);
    for (int r = 0; r < sys.m; ++r) {
        row_index.emplace(rows[(size_t)r],r);
    }

    std::unordered_map<int64_t,int> elem_index;
    elem_index.reserve(sys.candidates.size()*2 + 1);
    for (int c = 0; c < sys.n; ++c) {
        elem_index.emplace(sys.candidates[(size_t)c],c);
    }

    sys.A.assign((size_t)sys.m*(size_t)sys.n,0.0);
    sys.b.assign((size_t)sys.m,0.0);
    sys.x.assign((size_t)sys.n,0.0);
    std::vector<int> active((size_t)sys.n,0);

    for (const Stage3RankData& d : data) {
        const int nr = (int)d.row_key.size();
        const int nl = (int)d.local_elem_id.size();

        for (int r = 0; r < nr; ++r) {
            const int gr = row_index.at(d.row_key[(size_t)r]);
            sys.b[(size_t)gr] += d.rhs[(size_t)r];

            for (int lc = 0; lc < nl; ++lc) {
                const auto ei = elem_index.find(d.local_elem_id[(size_t)lc]);
                if (ei == elem_index.end()) continue;
                sys.A[(size_t)gr*(size_t)sys.n + (size_t)ei->second] +=
                    d.local_matrix[(size_t)r*(size_t)nl + (size_t)lc];
            }
        }

        int best = -1;
        double bestw = -1.0;
        for (int k = 0; k < (int)d.selected_elem_id.size(); ++k) {
            if (d.selected_weight[(size_t)k] > bestw) {
                bestw = d.selected_weight[(size_t)k];
                best = k;
            }
        }
        if (best >= 0) {
            const auto ei = elem_index.find(d.selected_elem_id[(size_t)best]);
            if (ei != elem_index.end()) active[(size_t)ei->second] = 1;
        }
    }

    long double scale_sq = 0.0L;
    for (double v : sys.A) scale_sq += (long double)v*(long double)v;
    sys.scale = std::sqrt((double)scale_sq);
    if (sys.scale <= DBL_MIN) sys.scale = 1.0;

    for (double& v : sys.A) v /= sys.scale;
    for (double& v : sys.b) v /= sys.scale;

    std::vector<double*> Arow((size_t)sys.m,NULL);
    for (int r = 0; r < sys.m; ++r) {
        Arow[(size_t)r] =
            sys.A.data() + (size_t)r*(size_t)sys.n;
    }

    monolis_optimize_nnls_R_with_sparse_solution_initial_set(
        Arow.data(),
        sys.b.data(),
        sys.x.data(),
        sys.m,
        sys.n,
        max_iter,
        solve_tol,
        active.data(),
        &sys.residual);

    return sys;
}


static void mode5_add_basis_vector(
    std::vector<std::vector<double>>& basis,
    const std::vector<double>& input)
{
    std::vector<double> v = input;
    for (const std::vector<double>& q : basis) {
        long double dot = 0.0L;
        for (size_t i = 0; i < v.size(); ++i) {
            dot += (long double)q[i]*(long double)v[i];
        }
        for (size_t i = 0; i < v.size(); ++i) {
            v[i] -= (double)dot*q[i];
        }
    }

    long double n2 = 0.0L;
    for (double a : v) n2 += (long double)a*(long double)a;
    const double nrm = std::sqrt((double)n2);
    if (nrm <= 1.0e-12) return;
    for (double& a : v) a /= nrm;
    basis.push_back(std::move(v));
}


static double mode5_signature_norm(const std::vector<double>& s)
{
    long double n2 = 0.0L;
    for (double v : s) n2 += (long double)v*(long double)v;
    return std::sqrt((double)n2);
}


static std::vector<double> mode5_assemble_signature(
    const std::vector<Stage3RankData>& data,
    const std::unordered_map<int64_t,std::vector<std::pair<int,int>>>& elem_locations,
    const std::unordered_map<RowKey,int,RowKeyHash>& row_index,
    int global_m,
    int64_t elem,
    const std::vector<int>& group_rows)
{
    std::vector<double> sig(group_rows.size(),0.0);
    if (group_rows.empty()) return sig;

    std::vector<int> row_to_local((size_t)global_m,-1);
    for (int j = 0; j < (int)group_rows.size(); ++j) {
        const int gr = group_rows[(size_t)j];
        if (gr >= 0 && gr < global_m) row_to_local[(size_t)gr] = j;
    }

    const auto loc_it = elem_locations.find(elem);
    if (loc_it == elem_locations.end()) return sig;

    for (const auto& loc : loc_it->second) {
        const int di = loc.first;
        const int lc = loc.second;
        const Stage3RankData& d = data[(size_t)di];
        const int nr = (int)d.row_key.size();
        const int nl = (int)d.local_elem_id.size();

        for (int r = 0; r < nr; ++r) {
            const auto ri = row_index.find(d.row_key[(size_t)r]);
            if (ri == row_index.end()) continue;
            const int gr = ri->second;
            if (gr < 0 || gr >= global_m) continue;
            const int j = row_to_local[(size_t)gr];
            if (j < 0) continue;
            sig[(size_t)j] +=
                d.local_matrix[(size_t)r*(size_t)nl + (size_t)lc];
        }
    }

    return sig;
}


static std::vector<Mode5OperatorGroupDesc> mode5_build_operator_groups(
    const std::vector<RowKey>& rows,
    const std::vector<InterfaceGroup>& interface_groups,
    const std::vector<InternalOperatorGroup>& internal_groups)
{
    std::unordered_map<int64_t,std::vector<int>> owner_rows;
    owner_rows.reserve(rows.size()/4 + 1);
    for (int r = 0; r < (int)rows.size(); ++r) {
        owner_rows[rows[(size_t)r].owner_global_id].push_back(r);
    }

    std::vector<Mode5OperatorGroupDesc> out;
    out.reserve(interface_groups.size() + internal_groups.size());

    for (const InterfaceGroup& ig : interface_groups) {
        Mode5OperatorGroupDesc g;
        g.kind = 0;
        g.owner_a = ig.key.source_global_id;
        g.owner_b = ig.key.target_global_id;
        g.members = ig.member_original_id;

        const auto a = owner_rows.find(g.owner_a);
        if (a != owner_rows.end()) {
            g.rows.insert(g.rows.end(),a->second.begin(),a->second.end());
        }
        const auto b = owner_rows.find(g.owner_b);
        if (b != owner_rows.end()) {
            g.rows.insert(g.rows.end(),b->second.begin(),b->second.end());
        }
        std::sort(g.rows.begin(),g.rows.end());
        g.rows.erase(std::unique(g.rows.begin(),g.rows.end()),g.rows.end());

        if (!g.rows.empty() && !g.members.empty()) out.push_back(std::move(g));
    }

    for (const InternalOperatorGroup& ig : internal_groups) {
        Mode5OperatorGroupDesc g;
        g.kind = 1;
        g.owner_a = (int64_t)ig.fine_subdomain;
        g.owner_b = -1;
        g.members = ig.member_original_id;
        const auto a = owner_rows.find(g.owner_a);
        if (a != owner_rows.end()) g.rows = a->second;
        if (!g.rows.empty() && !g.members.empty()) out.push_back(std::move(g));
    }

    return out;
}


static std::vector<Mode5DefectCandidate> mode5_evaluate_operator_defects(
    const Mode5System& sys,
    const std::vector<Mode5OperatorGroupDesc>& groups,
    const std::vector<std::vector<std::vector<double>>>& signatures,
    const std::vector<std::vector<double>>& signature_norm,
    double weight_tol,
    double signature_rtol,
    double defect_tol,
    bool return_new_only,
    int* deficient_group_count,
    double* max_defect_out)
{
    std::unordered_map<int64_t,int> cand_index;
    cand_index.reserve(sys.candidates.size()*2 + 1);
    for (int c = 0; c < sys.n; ++c) {
        cand_index.emplace(sys.candidates[(size_t)c],c);
    }

    std::vector<Mode5DefectCandidate> group_best;
    group_best.reserve(groups.size());

    int deficient = 0;
    double global_max_defect = 0.0;

    for (int gi = 0; gi < (int)groups.size(); ++gi) {
        const Mode5OperatorGroupDesc& g = groups[(size_t)gi];
        if (g.members.empty() || g.rows.empty()) continue;

        double max_norm = 0.0;
        for (double v : signature_norm[(size_t)gi]) {
            max_norm = std::max(max_norm,v);
        }
        if (max_norm <= DBL_MIN) continue;

        std::vector<std::vector<double>> basis;
        for (int j = 0; j < (int)g.members.size(); ++j) {
            const auto ci = cand_index.find(g.members[(size_t)j]);
            if (ci == cand_index.end()) continue;
            if (sys.x[(size_t)ci->second] <= weight_tol) continue;
            if (signature_norm[(size_t)gi][(size_t)j] <
                signature_rtol*max_norm) {
                continue;
            }
            mode5_add_basis_vector(
                basis,
                signatures[(size_t)gi][(size_t)j]);
        }

        std::vector<double> local_residual(g.rows.size(),0.0);
        long double local_r2 = 0.0L;
        for (int rr = 0; rr < (int)g.rows.size(); ++rr) {
            const int gr = g.rows[(size_t)rr];
            long double ax = 0.0L;
            for (int c = 0; c < sys.n; ++c) {
                ax += (long double)sys.A[
                    (size_t)gr*(size_t)sys.n + (size_t)c]
                    * (long double)sys.x[(size_t)c];
            }
            local_residual[(size_t)rr] =
                sys.b[(size_t)gr] - (double)ax;
            local_r2 +=
                (long double)local_residual[(size_t)rr]
                * (long double)local_residual[(size_t)rr];
        }
        const double local_rnorm = std::sqrt((double)local_r2);

        Mode5DefectCandidate best_any;
        best_any.group_index = gi;
        best_any.defect = -1.0;
        best_any.priority = -1.0;

        Mode5DefectCandidate best_new;
        best_new.group_index = gi;
        best_new.defect = -1.0;
        best_new.priority = -1.0;

        for (int j = 0; j < (int)g.members.size(); ++j) {
            const double sn = signature_norm[(size_t)gi][(size_t)j];
            if (sn < signature_rtol*max_norm || sn <= DBL_MIN) continue;

            const std::vector<double>& sig =
                signatures[(size_t)gi][(size_t)j];

            std::vector<double> v = sig;
            for (const std::vector<double>& q : basis) {
                long double dot = 0.0L;
                for (size_t k = 0; k < v.size(); ++k) {
                    dot += (long double)q[k]*(long double)v[k];
                }
                for (size_t k = 0; k < v.size(); ++k) {
                    v[k] -= (double)dot*q[k];
                }
            }

            const double novelty =
                mode5_signature_norm(v) / std::max(sn,DBL_MIN);
            const double importance = sn / max_norm;
            const double defect = importance*novelty;

            double alignment = 0.0;
            if (local_rnorm > DBL_MIN) {
                long double dot = 0.0L;
                for (size_t k = 0; k < sig.size(); ++k) {
                    dot += (long double)sig[k]
                         * (long double)local_residual[k];
                }
                alignment = std::max(
                    0.0,
                    std::min(
                        1.0,
                        (double)(dot /
                            ((long double)sn*(long double)local_rnorm))));
            }

            /*
             * Defect determines whether a direction is missing.  Residual
             * alignment only ranks equally important missing directions; it
             * cannot hide a geometrically missing operator direction.
             */
            const double priority =
                defect*(0.5 + 0.5*alignment);

            Mode5DefectCandidate cur;
            cur.defect = defect;
            cur.priority = priority;
            cur.importance = importance;
            cur.novelty = novelty;
            cur.residual_alignment = alignment;
            cur.group_index = gi;
            cur.elem = g.members[(size_t)j];
            cur.already_candidate =
                (cand_index.find(cur.elem) != cand_index.end());

            if (cur.defect > best_any.defect ||
                (cur.defect == best_any.defect &&
                 cur.priority > best_any.priority)) {
                best_any = cur;
            }

            if (!cur.already_candidate &&
                (cur.defect > best_new.defect ||
                 (cur.defect == best_new.defect &&
                  cur.priority > best_new.priority))) {
                best_new = cur;
            }
        }

        if (best_any.elem >= 0) {
            global_max_defect =
                std::max(global_max_defect,best_any.defect);
            if (best_any.defect > defect_tol) {
                deficient++;
                if (return_new_only) {
                    if (best_new.elem >= 0) group_best.push_back(best_new);
                }
                else {
                    group_best.push_back(best_any);
                }
            }
        }
    }

    if (deficient_group_count != NULL) *deficient_group_count = deficient;
    if (max_defect_out != NULL) *max_defect_out = global_max_defect;
    return group_best;
}


static Mode5AdaptiveResult mode5_adaptive_enrichment(
    const std::vector<Stage3RankData>& data,
    const std::vector<RowKey>& rows,
    const std::vector<InterfaceGroup>& interface_groups,
    const std::vector<InternalOperatorGroup>& internal_groups,
    std::unordered_map<int64_t,double>& candidate_seed_weight,
    int max_iter,
    double seed_tol,
    double weight_tol,
    double signature_rtol,
    double defect_tol,
    int max_round,
    int add_per_round,
    int hard_fallback,
    int verbose)
{
    Mode5AdaptiveResult out;

    std::vector<Mode5OperatorGroupDesc> groups =
        mode5_build_operator_groups(rows,interface_groups,internal_groups);

    if (groups.empty()) {
        fprintf(stderr,
            "WARNING: mode 5 found no operator groups; using mode-0 candidate set\n");
        return out;
    }

    std::unordered_map<RowKey,int,RowKeyHash> row_index;
    row_index.reserve(rows.size()*2 + 1);
    for (int r = 0; r < (int)rows.size(); ++r) {
        row_index.emplace(rows[(size_t)r],r);
    }

    std::unordered_map<int64_t,std::vector<std::pair<int,int>>> elem_locations;
    for (int di = 0; di < (int)data.size(); ++di) {
        const Stage3RankData& d = data[(size_t)di];
        for (int lc = 0; lc < (int)d.local_elem_id.size(); ++lc) {
            elem_locations[d.local_elem_id[(size_t)lc]].push_back(
                std::make_pair(di,lc));
        }
    }

    std::vector<std::vector<std::vector<double>>> signatures(groups.size());
    std::vector<std::vector<double>> signature_norm(groups.size());
    for (int gi = 0; gi < (int)groups.size(); ++gi) {
        const Mode5OperatorGroupDesc& g = groups[(size_t)gi];
        signatures[(size_t)gi].resize(g.members.size());
        signature_norm[(size_t)gi].resize(g.members.size(),0.0);
        for (int j = 0; j < (int)g.members.size(); ++j) {
            signatures[(size_t)gi][(size_t)j] =
                mode5_assemble_signature(
                    data,
                    elem_locations,
                    row_index,
                    (int)rows.size(),
                    g.members[(size_t)j],
                    g.rows);
            signature_norm[(size_t)gi][(size_t)j] =
                mode5_signature_norm(signatures[(size_t)gi][(size_t)j]);
        }
    }

    printf(
        "[MODE5-START] seed_candidates=%zu groups=%zu seed_tol=%.6e "
        "defect_tol=%.6e max_round=%d add_per_round=%d hard_fallback=%d\n",
        candidate_seed_weight.size(),
        groups.size(),
        seed_tol,
        defect_tol,
        max_round,
        add_per_round,
        hard_fallback);

    for (int round = 0; round < max_round; ++round) {
        Mode5System sys =
            mode5_solve_unconstrained_system(
                data,rows,candidate_seed_weight,max_iter,seed_tol);

        int deficient = 0;
        double max_defect = 0.0;
        std::vector<Mode5DefectCandidate> best =
            mode5_evaluate_operator_defects(
                sys,
                groups,
                signatures,
                signature_norm,
                weight_tol,
                signature_rtol,
                defect_tol,
                true,
                &deficient,
                &max_defect);

        int support = 0;
        for (double v : sys.x) if (v > weight_tol) support++;

        printf(
            "[MODE5-ROUND] round=%d candidates=%d support=%d residual=%.15e "
            "deficient_groups=%d max_defect=%.6e\n",
            round,
            sys.n,
            support,
            sys.residual,
            deficient,
            max_defect);

        out.rounds = round + 1;
        out.final_max_defect = max_defect;
        out.deficient_groups_final = deficient;

        if (deficient == 0 || max_defect <= defect_tol) break;

        std::vector<Mode5DefectCandidate> addable;
        for (const Mode5DefectCandidate& c : best) {
            if (!c.already_candidate && c.elem >= 0) addable.push_back(c);
        }

        std::sort(
            addable.begin(),addable.end(),
            [](const Mode5DefectCandidate& a,const Mode5DefectCandidate& b) {
                if (a.priority != b.priority) return a.priority > b.priority;
                if (a.defect != b.defect) return a.defect > b.defect;
                return a.elem < b.elem;
            });

        int added_this_round = 0;
        std::unordered_set<int64_t> added_now;
        for (const Mode5DefectCandidate& c : addable) {
            if (added_this_round >= add_per_round) break;
            if (!added_now.insert(c.elem).second) continue;
            if (candidate_seed_weight.find(c.elem) != candidate_seed_weight.end()) {
                continue;
            }

            candidate_seed_weight.emplace(c.elem,0.0);
            out.added_element.push_back(c.elem);
            out.added_candidates++;
            added_this_round++;

            if (verbose > 0) {
                const Mode5OperatorGroupDesc& g = groups[(size_t)c.group_index];
                printf(
                    "[MODE5-ADD] round=%d kind=%s owner={%lld,%lld} elem=%lld "
                    "defect=%.6e importance=%.6e novelty=%.6e alignment=%.6e\n",
                    round,
                    (g.kind == 0) ? "interface" : "internal",
                    (long long)g.owner_a,
                    (long long)g.owner_b,
                    (long long)c.elem,
                    c.defect,
                    c.importance,
                    c.novelty,
                    c.residual_alignment);
            }
        }

        if (added_this_round == 0) {
            printf(
                "[MODE5-STALL] round=%d no new column can be added; "
                "remaining deficient groups=%d\n",
                round,
                deficient);
            break;
        }
    }

    /* Final unconstrained diagnostic on the enriched pool. */
    Mode5System final_sys =
        mode5_solve_unconstrained_system(
            data,rows,candidate_seed_weight,max_iter,seed_tol);

    int final_deficient = 0;
    double final_max_defect = 0.0;
    std::vector<Mode5DefectCandidate> final_best =
        mode5_evaluate_operator_defects(
            final_sys,
            groups,
            signatures,
            signature_norm,
            weight_tol,
            signature_rtol,
            defect_tol,
            false,
            &final_deficient,
            &final_max_defect);

    out.deficient_groups_final = final_deficient;
    out.final_max_defect = final_max_defect;

    if (hard_fallback != 0 && final_deficient > 0) {
        std::unordered_set<int64_t> hard_seen;
        for (const Mode5DefectCandidate& c : final_best) {
            if (c.elem < 0 || c.defect <= defect_tol) continue;

            if (candidate_seed_weight.find(c.elem) == candidate_seed_weight.end()) {
                candidate_seed_weight.emplace(c.elem,0.0);
                out.added_element.push_back(c.elem);
                out.added_candidates++;
            }

            if (hard_seen.insert(c.elem).second) {
                out.hard_anchor_element.push_back(c.elem);
            }
        }
    }

    printf(
        "[MODE5-DONE] rounds=%d added_candidates=%d final_candidates=%zu "
        "deficient_groups=%d max_defect=%.6e hard_anchors=%zu\n",
        out.rounds,
        out.added_candidates,
        candidate_seed_weight.size(),
        out.deficient_groups_final,
        out.final_max_defect,
        out.hard_anchor_element.size());

    return out;
}

int main(int argc,char** argv)
{
    const auto timer_total_begin = HROMClock::now();

    if (argc < 3) {
        fprintf(
            stderr,
            "usage: %s <directory> <num_rank_files> "
            "[inner_max_iter=1000] [inner_tol=1e-8] "
            "[weight_tol=1e-8] [coverage_tau=weight_tol] "
            "[admm_max_iter=100] [rho=1.0] "
            "[admm_abs_tol=1e-10] [admm_rel_tol=1e-6] "
            "[verbose=1] [coverage_mode=2] "
            "[operator_clusters=2] [operator_signature_rtol=1e-8] "
            "[threshold_refit=1] [threshold_refit_max_pass=32] "
            "[mode5_seed_tol=max(inner_tol,1e-4)] [mode5_defect_tol=0.10] "
            "[mode5_max_round=32] [mode5_add_per_round=1] [mode5_hard_fallback=0]\n"
            "  coverage_mode=0 : baseline Stage-3 (Stage-2 candidates only, no coverage)\n"
            "  coverage_mode=1 : coverage candidate expansion only, no hard constraint\n"
            "  coverage_mode=2 : physical hard interface coverage\n"
            "  coverage_mode=3 : operator-aware hard interface coverage (physical + signature clusters)\n"
            "  coverage_mode=4 : domain-wide operator-aware coverage (interface + internal signatures)\n"
            "  coverage_mode=5 : mode-0-seeded adaptive operator enrichment (interface + internal defects)\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    const char* directory = argv[1];
    const int num_ranks = atoi(argv[2]);

    const int requested_iter =
        (argc > 3) ? atoi(argv[3]) : 1000;

    const int max_iter =
        (requested_iter > 0)
        ? requested_iter
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
     *   2 = physical hard coverage:
     *       candidate expansion + Cx >= epsilon.
     *
     *   3 = operator-aware hard interface coverage:
     *       physical coverage + interface operator-signature subgroups.
     *
     *   4 = domain-wide operator-aware hard coverage:
     *       mode 3 + internal operator-signature subgroups for every fine
     *       subdomain.  Internal signatures use A(R_s,c), where R_s are rows
     *       owned by fine subdomain s and c is a strictly internal element.
     *       Interface-supporting elements are excluded from the internal sets.
     *
     *   5 = mode-0-seeded adaptive operator enrichment:
     *       start from the Stage-2 selected union with an unconstrained sparse
     *       NNLS, detect missing interface/internal operator directions, add
     *       only the highest-priority missing columns, and repeat.  Optional
     *       singleton hard constraints are applied only as a final fallback.
     */
    const int coverage_mode =
        (argc > 12) ? atoi(argv[12]) : 2;

    const int operator_clusters =
        (argc > 13) ? atoi(argv[13]) : 2;

    const double operator_signature_rtol =
        (argc > 14) ? atof(argv[14]) : 1.0e-8;

    const int threshold_refit =
        (argc > 15) ? atoi(argv[15]) : 1;

    const int threshold_refit_max_pass =
        (argc > 16) ? atoi(argv[16]) : 32;

    const double mode5_seed_tol =
        (argc > 17) ? atof(argv[17]) : std::max(tol,1.0e-4);

    const double mode5_defect_tol =
        (argc > 18) ? atof(argv[18]) : 1.0e-1;

    const int mode5_max_round =
        (argc > 19) ? atoi(argv[19]) : 32;

    const int mode5_add_per_round =
        (argc > 20) ? atoi(argv[20]) : 1;

    const int mode5_hard_fallback =
        (argc > 21) ? atoi(argv[21]) : 0;

    if (num_ranks <= 0) {
        die("num_rank_files must be positive");
    }

    if (coverage_mode < 0 || coverage_mode > 5) {
        die("coverage_mode must be 0, 1, 2, 3, 4, or 5");
    }

    if (operator_clusters <= 0 || operator_clusters > 8 ||
        operator_signature_rtol <= 0.0 || operator_signature_rtol >= 1.0) {
        die("invalid operator-aware coverage option");
    }

    if (mode5_seed_tol <= 0.0 || mode5_defect_tol <= 0.0 ||
        mode5_defect_tol >= 1.0 || mode5_max_round <= 0 ||
        mode5_add_per_round <= 0 ||
        (mode5_hard_fallback != 0 && mode5_hard_fallback != 1)) {
        die("invalid mode-5 adaptive-enrichment option");
    }

    if ((threshold_refit != 0 && threshold_refit != 1) ||
        threshold_refit_max_pass <= 0 || weight_tol < 0.0) {
        die("invalid threshold-refit option");
    }

    if (threshold_refit != 0 &&
        ((coverage_mode >= 2 && coverage_mode <= 4) ||
         (coverage_mode == 5 && mode5_hard_fallback != 0)) &&
        !(2.0 * coverage_tau > weight_tol)) {
        die("threshold refit requires 2*coverage_tau > weight_tol so every coverage group is guaranteed to retain an output-weight member");
    }

    if (coverage_tau <= 0.0 ||
        admm_max_iter <= 0 ||
        admm_rho <= 0.0 ||
        admm_abs_tol <= 0.0 ||
        admm_rel_tol <= 0.0) {
        die("invalid interface-coverage ADMM option");
    }

    double time_read_rank = 0.0;
    double time_interface = 0.0;
    double time_candidate_expand = 0.0;
    double time_index_group = 0.0;
    double time_matrix_assembly = 0.0;
    double time_normalize_anchor = 0.0;
    double time_solver = 0.0;
    double time_postprocess = 0.0;
    double time_route_stats = 0.0;
    double time_legacy_write = 0.0;

    const auto timer_read_begin = HROMClock::now();

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

    time_read_rank =
        HROM_elapsed_sec(timer_read_begin, HROMClock::now());

    std::vector<RowKey> rows(row_set.begin(),row_set.end());
    std::sort(rows.begin(),rows.end(),[](const RowKey&a,const RowKey&b){
        if (a.snapshot!=b.snapshot) return a.snapshot<b.snapshot;
        if (a.owner_global_id!=b.owner_global_id) return a.owner_global_id<b.owner_global_id;
        return a.mode_offset<b.mode_offset;
    });

    const int stage2_selected_union_count =
        (int)candidate_seed_weight.size();

    std::vector<InterfaceGroup> interface_groups;
    std::vector<RouteEntry> routes;
    std::vector<InternalOperatorGroup> internal_operator_groups;
    Mode5AdaptiveResult mode5_result;

    /*
     * coverage_mode=0 is the true pre-coverage comparison baseline.
     * Do not even read/use the interface files in this mode.
     */
    const auto timer_interface_begin = HROMClock::now();

    if (coverage_mode > 0) {
        interface_groups =
            read_all_interface_groups(
                directory,
                num_ranks);

        if (interface_groups.empty()) {
            die("no physical interface groups were read");
        }
    }

    /*
     * Mode 4 needs the exact owned-element routes before candidate expansion
     * so that strictly-internal operator groups can reintroduce candidates.
     * The same route table is reused later for the legacy output/report.
     */
    if (coverage_mode == 4 || coverage_mode == 5) {
        routes = read_all_routes(directory,num_ranks);
        internal_operator_groups =
            build_internal_operator_groups(routes,interface_groups);

        if (internal_operator_groups.empty()) {
            fprintf(
                stderr,
                "WARNING: mode 4/5 found no strictly-internal fine-subdomain groups; "
                "the calculation reduces to mode-3 interface coverage\n");
        }
    }

    time_interface =
        HROM_elapsed_sec(timer_interface_begin, HROMClock::now());

    if (coverage_mode == 5) {
        mode5_result = mode5_adaptive_enrichment(
            data,
            rows,
            interface_groups,
            internal_operator_groups,
            candidate_seed_weight,
            max_iter,
            mode5_seed_tol,
            weight_tol,
            operator_signature_rtol,
            mode5_defect_tol,
            mode5_max_round,
            mode5_add_per_round,
            mode5_hard_fallback,
            admm_verbose);
    }

    std::unordered_set<int64_t> coverage_added_element;
    std::unordered_set<int64_t> internal_coverage_added_element;
    if (coverage_mode == 5) {
        for (int64_t elem : mode5_result.added_element) {
            coverage_added_element.insert(elem);
        }
    }
    int empty_interface_before_expansion = 0;
    int empty_internal_before_expansion = 0;

    /*
     * Candidate expansion is disabled only in mode 0.
     * Modes 1/2 guarantee physical interface availability; modes 3/4 make
     * enough interface candidates available for operator clustering; mode 4
     * additionally does the same for strictly-internal fine-subdomain groups.
     */
    const auto timer_candidate_begin = HROMClock::now();

    if (coverage_mode > 0 && coverage_mode != 5) {
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
            int current_count = 0;
            for (int64_t elem : g.member_original_id) {
                if (candidate_seed_weight.find(elem) !=
                    candidate_seed_weight.end()) {
                    current_count++;
                }
            }

            if (current_count == 0) {
                empty_interface_before_expansion++;
            }

            const int target_count =
                (coverage_mode >= 3)
                ? std::min(operator_clusters, (int)g.member_original_id.size())
                : 1;

            if (current_count >= target_count) {
                continue;
            }

            std::vector<std::pair<double,int64_t>> ranked;
            ranked.reserve(g.member_original_id.size());

            for (int64_t elem : g.member_original_id) {
                if (candidate_seed_weight.find(elem) !=
                    candidate_seed_weight.end()) {
                    continue;
                }

                const auto it = all_column_norm_sq.find(elem);
                if (it == all_column_norm_sq.end()) continue;
                ranked.emplace_back(it->second, elem);
            }

            std::sort(
                ranked.begin(), ranked.end(),
                [](const std::pair<double,int64_t>& a,
                   const std::pair<double,int64_t>& b) {
                    if (a.first != b.first) return a.first > b.first;
                    return a.second < b.second;
                });

            for (const auto& ne : ranked) {
                if (current_count >= target_count) break;
                candidate_seed_weight.emplace(ne.second, 0.0);
                coverage_added_element.insert(ne.second);
                current_count++;
            }

            if (current_count < target_count) {
                fprintf(
                    stderr,
                    "ERROR: physical interface {%lld,%lld} has only %d/%d "
                    "available Stage-3 candidate columns\n",
                    (long long)g.key.source_global_id,
                    (long long)g.key.target_global_id,
                    current_count,
                    target_count);
                exit(EXIT_FAILURE);
            }
        }

        /*
         * Mode 4: provide at least operator_clusters candidates in every
         * strictly-internal fine-subdomain group, mirroring the interface
         * operator candidate expansion of mode 3.  The actual operator
         * directions are determined after the normalized Stage-3 matrix A is
         * assembled.
         */
        if (coverage_mode == 4) {
            for (const InternalOperatorGroup& g : internal_operator_groups) {
                if (g.member_original_id.empty()) continue;

                int current_count = 0;
                for (int64_t elem : g.member_original_id) {
                    if (candidate_seed_weight.find(elem) !=
                        candidate_seed_weight.end()) {
                        current_count++;
                    }
                }

                if (current_count == 0) {
                    empty_internal_before_expansion++;
                }

                const int target_count =
                    std::min(operator_clusters,
                             (int)g.member_original_id.size());

                if (current_count >= target_count) continue;

                std::vector<std::pair<double,int64_t>> ranked;
                ranked.reserve(g.member_original_id.size());

                for (int64_t elem : g.member_original_id) {
                    if (candidate_seed_weight.find(elem) !=
                        candidate_seed_weight.end()) {
                        continue;
                    }

                    const auto it = all_column_norm_sq.find(elem);
                    if (it == all_column_norm_sq.end()) continue;
                    ranked.emplace_back(it->second,elem);
                }

                std::sort(
                    ranked.begin(),ranked.end(),
                    [](const std::pair<double,int64_t>& a,
                       const std::pair<double,int64_t>& b) {
                        if (a.first != b.first) return a.first > b.first;
                        return a.second < b.second;
                    });

                for (const auto& ne : ranked) {
                    if (current_count >= target_count) break;
                    candidate_seed_weight.emplace(ne.second,0.0);
                    coverage_added_element.insert(ne.second);
                    internal_coverage_added_element.insert(ne.second);
                    current_count++;
                }

                if (current_count < target_count) {
                    fprintf(
                        stderr,
                        "ERROR: internal subdomain %d has only %d/%d "
                        "available Stage-3 candidate columns\n",
                        g.fine_subdomain,current_count,target_count);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    time_candidate_expand =
        HROM_elapsed_sec(timer_candidate_begin, HROMClock::now());

    const auto timer_index_begin = HROMClock::now();

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

    if (coverage_mode > 0 && coverage_mode != 5) {
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

    int num_mode5_hard_group = 0;
    if (coverage_mode == 5 && !mode5_result.hard_anchor_element.empty()) {
        for (int64_t elem : mode5_result.hard_anchor_element) {
            const auto it = elem_index.find(elem);
            if (it == elem_index.end()) {
                fprintf(stderr,
                    "ERROR: mode-5 hard anchor %lld is absent from final candidate pool\n",
                    (long long)elem);
                exit(EXIT_FAILURE);
            }
            group_item.push_back(it->second);
            group_ptr.push_back((int)group_item.size());
            group_epsilon.push_back(2.0*coverage_tau);
            num_mode5_hard_group++;
        }
    }

    const int num_interface_group =
        (coverage_mode == 5) ? 0 : (int)group_epsilon.size();

    time_index_group =
        HROM_elapsed_sec(timer_index_begin, HROMClock::now());

    const auto timer_assembly_begin = HROMClock::now();

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

    time_matrix_assembly =
        HROM_elapsed_sec(timer_assembly_begin, HROMClock::now());

    const auto timer_norm_begin = HROMClock::now();

    /*
     * In mode 1 the added columns are merely available candidates; they are
     * NOT force-seeded into the active set.  This makes mode 1 a clean
     * "candidate pool only" ablation.
     *
     * Mode 2 uses them as warm-start columns for the constrained solve.
     */
    if (coverage_mode >= 2) {
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


    /*
     * ------------------------------------------------------------
     * Operator-aware interface coverage (coverage_mode=3)
     *
     * The physical groups above only state that some positive cubature
     * weight remains on each interface.  For mode 3, augment them with
     * subgroups obtained by clustering Stage-3 element operator signatures.
     *
     * For interface {p,q}, the signature of candidate c is
     *
     *     s_{g,c} = A(R_p union R_q, c),
     *
     * where R_p/R_q are rows whose owner_global_id is p/q.  Thus the
     * signature contains all training snapshots and retained mode offsets
     * associated with the two sides of the interface.  Greedy farthest-point
     * selection in cosine distance chooses up to operator_clusters distinct
     * directions, and every non-negligible member is assigned to its nearest
     * anchor.  Each cluster receives the same pigeon-hole coverage rule
     *
     *     sum_{c in cluster} x_c >= 2*tau*|cluster|,
     *
     * so at least one output-weight element survives in every represented
     * operator direction.  The original physical constraints are retained.
     * ------------------------------------------------------------
     */
    std::vector<int> solver_group_ptr = group_ptr;
    std::vector<int> solver_group_item = group_item;
    std::vector<double> solver_group_epsilon = group_epsilon;
    std::vector<int> operator_group_parent;
    std::vector<int> operator_group_cluster;
    std::vector<int> operator_group_anchor;
    std::vector<double> operator_group_anchor_norm;
    std::vector<int> internal_operator_group_anchor;

    int num_interface_operator_subgroups = 0;
    int num_internal_operator_subgroups = 0;

    std::unordered_map<int64_t,std::vector<int>> owner_rows;
    if (coverage_mode == 3 || coverage_mode == 4) {
        owner_rows.reserve(rows.size()/4 + 1);
        for (int r = 0; r < m; ++r) {
            owner_rows[rows[(size_t)r].owner_global_id].push_back(r);
        }
    }

    if (coverage_mode == 3 || coverage_mode == 4) {

        int physical_with_operator_rows = 0;
        int physical_without_operator_rows = 0;
        int physical_zero_signature = 0;
        int total_operator_clusters = 0;

        for (int g = 0; g < num_interface_group; ++g) {
            const InterfaceGroup& ig = interface_groups[(size_t)g];
            std::vector<int> op_rows;

            const auto rp = owner_rows.find(ig.key.source_global_id);
            if (rp != owner_rows.end()) {
                op_rows.insert(op_rows.end(), rp->second.begin(), rp->second.end());
            }
            const auto rq = owner_rows.find(ig.key.target_global_id);
            if (rq != owner_rows.end()) {
                op_rows.insert(op_rows.end(), rq->second.begin(), rq->second.end());
            }

            std::sort(op_rows.begin(), op_rows.end());
            op_rows.erase(std::unique(op_rows.begin(), op_rows.end()), op_rows.end());

            if (op_rows.empty()) {
                physical_without_operator_rows++;
                if (admm_verbose > 0) {
                    printf(
                        "[OPERATOR-COVERAGE-SKIP] g=%d interface={%lld,%lld} reason=no_owner_rows\n",
                        g,
                        (long long)ig.key.source_global_id,
                        (long long)ig.key.target_global_id);
                }
                continue;
            }
            physical_with_operator_rows++;

            const int kS = group_ptr[(size_t)g];
            const int kE = group_ptr[(size_t)g + 1];
            const int nm = kE - kS;
            std::vector<double> sig_norm((size_t)nm, 0.0);
            double max_sig_norm = 0.0;

            for (int j = 0; j < nm; ++j) {
                const int c = group_item[(size_t)(kS + j)];
                long double n2 = 0.0L;
                for (int rr : op_rows) {
                    const long double v = (long double)A[rr][c];
                    n2 += v*v;
                }
                sig_norm[(size_t)j] = std::sqrt((double)n2);
                max_sig_norm = std::max(max_sig_norm, sig_norm[(size_t)j]);
            }

            if (max_sig_norm <= DBL_MIN) {
                physical_zero_signature++;
                if (admm_verbose > 0) {
                    printf(
                        "[OPERATOR-COVERAGE-SKIP] g=%d interface={%lld,%lld} reason=zero_signature members=%d rows=%zu\n",
                        g,
                        (long long)ig.key.source_global_id,
                        (long long)ig.key.target_global_id,
                        nm,
                        op_rows.size());
                }
                continue;
            }

            std::vector<int> eligible;
            for (int j = 0; j < nm; ++j) {
                if (sig_norm[(size_t)j] >=
                    operator_signature_rtol * max_sig_norm) {
                    eligible.push_back(j);
                }
            }
            if (eligible.empty()) {
                physical_zero_signature++;
                continue;
            }

            auto cosine = [&](int ja, int jb) -> double {
                const int ca = group_item[(size_t)(kS + ja)];
                const int cb = group_item[(size_t)(kS + jb)];
                long double dot = 0.0L;
                for (int rr : op_rows) {
                    dot += (long double)A[rr][ca] * (long double)A[rr][cb];
                }
                const double den = sig_norm[(size_t)ja] * sig_norm[(size_t)jb];
                if (den <= DBL_MIN) return 1.0;
                double cs = (double)(dot / (long double)den);
                return std::max(-1.0, std::min(1.0, cs));
            };

            std::vector<int> anchors;
            int first = eligible.front();
            for (int j : eligible) {
                if (sig_norm[(size_t)j] > sig_norm[(size_t)first]) first = j;
            }
            anchors.push_back(first);

            const int requested = std::min(operator_clusters, (int)eligible.size());
            const double diversity_floor = 1.0e-4;

            while ((int)anchors.size() < requested) {
                int best = -1;
                double best_score = -1.0;
                for (int j : eligible) {
                    if (std::find(anchors.begin(), anchors.end(), j) != anchors.end()) continue;
                    double min_distance = 2.0;
                    for (int a : anchors) {
                        min_distance = std::min(min_distance, 1.0 - cosine(j,a));
                    }
                    const double importance = sig_norm[(size_t)j] / max_sig_norm;
                    const double score = importance * min_distance;
                    if (score > best_score) {
                        best_score = score;
                        best = j;
                    }
                }
                if (best < 0 || best_score < diversity_floor) break;
                anchors.push_back(best);
            }

            std::vector<std::vector<int>> clusters(anchors.size());
            for (int j : eligible) {
                int best_a = 0;
                double best_cos = -2.0;
                for (int aidx = 0; aidx < (int)anchors.size(); ++aidx) {
                    const double cs = cosine(j, anchors[(size_t)aidx]);
                    if (cs > best_cos) {
                        best_cos = cs;
                        best_a = aidx;
                    }
                }
                clusters[(size_t)best_a].push_back(
                    group_item[(size_t)(kS + j)]);
            }

            for (int aidx = 0; aidx < (int)anchors.size(); ++aidx) {
                std::vector<int>& cl = clusters[(size_t)aidx];
                if (cl.empty()) continue;
                std::sort(cl.begin(), cl.end());
                cl.erase(std::unique(cl.begin(), cl.end()), cl.end());

                solver_group_item.insert(
                    solver_group_item.end(), cl.begin(), cl.end());
                solver_group_ptr.push_back((int)solver_group_item.size());
                solver_group_epsilon.push_back(
                    2.0 * coverage_tau * (double)cl.size());

                const int aj = anchors[(size_t)aidx];
                const int anchor_col = group_item[(size_t)(kS + aj)];
                operator_group_parent.push_back(g);
                operator_group_cluster.push_back(aidx);
                operator_group_anchor.push_back(anchor_col);
                operator_group_anchor_norm.push_back(sig_norm[(size_t)aj]);
                total_operator_clusters++;

                if (admm_verbose > 0) {
                    printf(
                        "[OPERATOR-COVERAGE-GROUP] physical_g=%d interface={%lld,%lld} cluster=%d/%zu members=%zu anchor_col=%d anchor_elem=%lld anchor_norm=%.6e rows=%zu epsilon=%.6e\n",
                        g,
                        (long long)ig.key.source_global_id,
                        (long long)ig.key.target_global_id,
                        aidx + 1,
                        anchors.size(),
                        cl.size(),
                        anchor_col,
                        (long long)candidates[(size_t)anchor_col],
                        sig_norm[(size_t)aj],
                        op_rows.size(),
                        solver_group_epsilon.back());
                }
            }
        }

        printf(
            "[OPERATOR-COVERAGE-SUMMARY] physical=%d with_owner_rows=%d no_owner_rows=%d zero_signature=%d operator_subgroups=%d total_constraints=%zu clusters_requested=%d signature_rtol=%.3e\n",
            num_interface_group,
            physical_with_operator_rows,
            physical_without_operator_rows,
            physical_zero_signature,
            total_operator_clusters,
            solver_group_epsilon.size(),
            operator_clusters,
            operator_signature_rtol);

        num_interface_operator_subgroups = total_operator_clusters;
    }

    const int first_internal_solver_group =
        num_interface_group + num_interface_operator_subgroups;

    /*
     * ------------------------------------------------------------
     * Internal operator-aware coverage (coverage_mode=4)
     *
     * For fine subdomain s, use strictly-internal candidate elements only and
     * define
     *
     *     s_{s,c} = A(R_s,c),
     *
     * where R_s contains all Stage-3 rows whose owner_global_id is s.  The
     * same cosine-distance clustering used for interface signatures is applied
     * to these internal signatures.  Every represented internal operator
     * direction receives an independent hard coverage constraint.
     * ------------------------------------------------------------
     */
    if (coverage_mode == 4) {
        int internal_with_operator_rows = 0;
        int internal_without_operator_rows = 0;
        int internal_zero_signature = 0;
        int total_internal_operator_clusters = 0;

        for (const InternalOperatorGroup& ig : internal_operator_groups) {
            const auto rr_it = owner_rows.find((int64_t)ig.fine_subdomain);
            if (rr_it == owner_rows.end() || rr_it->second.empty()) {
                internal_without_operator_rows++;
                if (admm_verbose > 0) {
                    printf(
                        "[INTERNAL-OPERATOR-COVERAGE-SKIP] subdomain=%d reason=no_owner_rows\n",
                        ig.fine_subdomain);
                }
                continue;
            }

            const std::vector<int>& op_rows = rr_it->second;
            internal_with_operator_rows++;

            std::vector<int> member_cols;
            member_cols.reserve(ig.member_original_id.size());
            for (int64_t elem : ig.member_original_id) {
                const auto ei = elem_index.find(elem);
                if (ei != elem_index.end()) member_cols.push_back(ei->second);
            }
            std::sort(member_cols.begin(),member_cols.end());
            member_cols.erase(
                std::unique(member_cols.begin(),member_cols.end()),
                member_cols.end());

            if (member_cols.empty()) {
                fprintf(
                    stderr,
                    "ERROR: mode-4 internal subdomain %d has no candidate "
                    "column after internal candidate expansion\n",
                    ig.fine_subdomain);
                exit(EXIT_FAILURE);
            }

            const int nm = (int)member_cols.size();
            std::vector<double> sig_norm((size_t)nm,0.0);
            double max_sig_norm = 0.0;

            for (int j = 0; j < nm; ++j) {
                const int c = member_cols[(size_t)j];
                long double n2 = 0.0L;
                for (int rr : op_rows) {
                    const long double v = (long double)A[rr][c];
                    n2 += v*v;
                }
                sig_norm[(size_t)j] = std::sqrt((double)n2);
                max_sig_norm = std::max(max_sig_norm,sig_norm[(size_t)j]);
            }

            if (max_sig_norm <= DBL_MIN) {
                internal_zero_signature++;
                if (admm_verbose > 0) {
                    printf(
                        "[INTERNAL-OPERATOR-COVERAGE-SKIP] subdomain=%d "
                        "reason=zero_signature members=%d rows=%zu\n",
                        ig.fine_subdomain,nm,op_rows.size());
                }
                continue;
            }

            std::vector<int> eligible;
            for (int j = 0; j < nm; ++j) {
                if (sig_norm[(size_t)j] >=
                    operator_signature_rtol * max_sig_norm) {
                    eligible.push_back(j);
                }
            }
            if (eligible.empty()) {
                internal_zero_signature++;
                continue;
            }

            auto cosine = [&](int ja,int jb) -> double {
                const int ca = member_cols[(size_t)ja];
                const int cb = member_cols[(size_t)jb];
                long double dot = 0.0L;
                for (int rr : op_rows) {
                    dot += (long double)A[rr][ca] * (long double)A[rr][cb];
                }
                const double den =
                    sig_norm[(size_t)ja] * sig_norm[(size_t)jb];
                if (den <= DBL_MIN) return 1.0;
                double cs = (double)(dot / (long double)den);
                return std::max(-1.0,std::min(1.0,cs));
            };

            std::vector<int> anchors;
            int first = eligible.front();
            for (int j : eligible) {
                if (sig_norm[(size_t)j] > sig_norm[(size_t)first]) first = j;
            }
            anchors.push_back(first);

            const int requested =
                std::min(operator_clusters,(int)eligible.size());
            const double diversity_floor = 1.0e-4;

            while ((int)anchors.size() < requested) {
                int best = -1;
                double best_score = -1.0;
                for (int j : eligible) {
                    if (std::find(anchors.begin(),anchors.end(),j) !=
                        anchors.end()) {
                        continue;
                    }
                    double min_distance = 2.0;
                    for (int a : anchors) {
                        min_distance =
                            std::min(min_distance,1.0 - cosine(j,a));
                    }
                    const double importance =
                        sig_norm[(size_t)j] / max_sig_norm;
                    const double score = importance * min_distance;
                    if (score > best_score) {
                        best_score = score;
                        best = j;
                    }
                }
                if (best < 0 || best_score < diversity_floor) break;
                anchors.push_back(best);
            }

            std::vector<std::vector<int>> clusters(anchors.size());
            for (int j : eligible) {
                int best_a = 0;
                double best_cos = -2.0;
                for (int aidx = 0;
                     aidx < (int)anchors.size();
                     ++aidx) {
                    const double cs = cosine(j,anchors[(size_t)aidx]);
                    if (cs > best_cos) {
                        best_cos = cs;
                        best_a = aidx;
                    }
                }
                clusters[(size_t)best_a].push_back(
                    member_cols[(size_t)j]);
            }

            for (int aidx = 0;
                 aidx < (int)anchors.size();
                 ++aidx) {
                std::vector<int>& cl = clusters[(size_t)aidx];
                if (cl.empty()) continue;
                std::sort(cl.begin(),cl.end());
                cl.erase(std::unique(cl.begin(),cl.end()),cl.end());

                solver_group_item.insert(
                    solver_group_item.end(),cl.begin(),cl.end());
                solver_group_ptr.push_back((int)solver_group_item.size());
                solver_group_epsilon.push_back(
                    2.0 * coverage_tau * (double)cl.size());

                const int aj = anchors[(size_t)aidx];
                const int anchor_col = member_cols[(size_t)aj];
                internal_operator_group_anchor.push_back(anchor_col);
                total_internal_operator_clusters++;

                if (admm_verbose > 0) {
                    printf(
                        "[INTERNAL-OPERATOR-COVERAGE-GROUP] subdomain=%d "
                        "cluster=%d/%zu members=%zu anchor_col=%d "
                        "anchor_elem=%lld anchor_norm=%.6e rows=%zu epsilon=%.6e\n",
                        ig.fine_subdomain,
                        aidx + 1,
                        anchors.size(),
                        cl.size(),
                        anchor_col,
                        (long long)candidates[(size_t)anchor_col],
                        sig_norm[(size_t)aj],
                        op_rows.size(),
                        solver_group_epsilon.back());
                }
            }
        }

        num_internal_operator_subgroups = total_internal_operator_clusters;

        printf(
            "[INTERNAL-OPERATOR-COVERAGE-SUMMARY] groups=%zu "
            "with_owner_rows=%d no_owner_rows=%d zero_signature=%d "
            "operator_subgroups=%d total_constraints=%zu clusters_requested=%d "
            "signature_rtol=%.3e\n",
            internal_operator_groups.size(),
            internal_with_operator_rows,
            internal_without_operator_rows,
            internal_zero_signature,
            total_internal_operator_clusters,
            solver_group_epsilon.size(),
            operator_clusters,
            operator_signature_rtol);
    }

    const int num_solver_group = (int)solver_group_epsilon.size();
    const bool use_hard_coverage_solver =
        ((coverage_mode >= 2 && coverage_mode <= 4) ||
         (coverage_mode == 5 && num_solver_group > 0));

    std::vector<int> coverage_anchor(
        (size_t)num_interface_group,
        -1);

    std::unordered_set<int> unique_coverage_anchor;
    unique_coverage_anchor.reserve(
        (size_t)num_interface_group * 2 + 1);

    int singleton_interface_groups = 0;

    if (coverage_mode >= 2 && coverage_mode <= 4) {
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


    if (coverage_mode == 3 || coverage_mode == 4) {
        for (int c : operator_group_anchor) {
            if (c >= 0 && c < n) {
                active[(size_t)c] = 1;
                unique_coverage_anchor.insert(c);
            }
        }
    }

    if (coverage_mode == 4) {
        for (int c : internal_operator_group_anchor) {
            if (c >= 0 && c < n) {
                active[(size_t)c] = 1;
                unique_coverage_anchor.insert(c);
            }
        }
    }

    int nactive=0;
    for(int v:active) if(v)++nactive;

    const char* coverage_mode_name =
        (coverage_mode == 0)
        ? "baseline: Stage-2 candidates only, unconstrained NNLS"
        : (coverage_mode == 1)
          ? "candidate-expanded ablation, unconstrained NNLS"
        : (coverage_mode == 2)
          ? "physical hard interface coverage"
        : (coverage_mode == 3)
          ? "operator-aware interface coverage: physical + interface signatures"
        : (coverage_mode == 4)
          ? "domain-wide operator-aware coverage: interface + internal signatures"
          : "mode-0-seeded adaptive operator enrichment";

    printf("============================================================\n");
    printf("[Stage 3 v16 mode-5 adaptive operator enrichment]\n");
    printf("coverage mode                      = %d\n", coverage_mode);
    printf("coverage mode description          = %s\n", coverage_mode_name);
    printf("rank files                         = %d\n", num_ranks);
    printf("global keyed rows                  = %d\n", m);
    printf("Stage-2 selected union elements    = %d\n",
        stage2_selected_union_count);
    printf("coverage-added candidate elements  = %zu\n",
        coverage_added_element.size());
    printf("Stage-3 candidate pool             = %d\n", n);
    printf("physical interface constraints     = %d\n",
        num_interface_group);
    printf("interface operator subconstraints  = %d\n",
        num_interface_operator_subgroups);
    printf("internal operator subconstraints   = %d\n",
        num_internal_operator_subgroups);
    printf("total solver coverage constraints  = %d\n",
        num_solver_group);
    if (coverage_mode == 3 || coverage_mode == 4) {
        printf("operator clusters requested        = %d\n", operator_clusters);
        printf("operator signature relative floor  = %.6e\n", operator_signature_rtol);
    }
    if (coverage_mode == 4) {
        printf("strictly-internal groups           = %zu\n",
            internal_operator_groups.size());
        printf("internal-added candidate elements  = %zu\n",
            internal_coverage_added_element.size());
        printf("internal groups empty before expansion = %d\n",
            empty_internal_before_expansion);
    }
    if (coverage_mode == 5) {
        printf("mode5 seed tolerance               = %.6e\n", mode5_seed_tol);
        printf("mode5 operator defect tolerance    = %.6e\n", mode5_defect_tol);
        printf("mode5 adaptive rounds              = %d\n", mode5_result.rounds);
        printf("mode5 added candidates             = %d\n", mode5_result.added_candidates);
        printf("mode5 final maximum defect         = %.6e\n", mode5_result.final_max_defect);
        printf("mode5 final deficient groups       = %d\n", mode5_result.deficient_groups_final);
        printf("mode5 hard fallback groups         = %d\n", num_mode5_hard_group);
    }
    printf("interfaces empty before expansion  = %d\n",
        empty_interface_before_expansion);
    printf("initial active elements            = %d\n", nactive);
    printf("unique coverage anchors            = %zu\n",
        unique_coverage_anchor.size());
    printf("singleton interface constraints    = %d\n",
        singleton_interface_groups);
    printf("inner NNLS max iterations          = %d\n", max_iter);
    printf("inner NNLS tolerance               = %.6e\n", tol);
    printf("online output weight threshold     = %.6e\n", weight_tol);
    printf("coverage selection threshold tau   = %.6e\n",
        coverage_tau);
    printf("threshold-consistent refit         = %s\n",
        threshold_refit ? "enabled" : "disabled");
    printf("threshold refit max passes         = %d\n",
        threshold_refit_max_pass);
    printf("ADMM max iterations                = %d\n",
        admm_max_iter);
    printf("ADMM initial rho                   = %.6e\n",
        admm_rho);
    printf("global normalization scale         = %.15e\n",
        scale);
    printf("============================================================\n");

    time_normalize_anchor =
        HROM_elapsed_sec(timer_norm_begin, HROMClock::now());

    const auto timer_solver_begin = HROMClock::now();

    double residual = 0.0;
    double max_coverage_violation = 0.0;
    int admm_outer_iter = 0;
    int coverage_status = 0;

    if (use_hard_coverage_solver) {
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
                num_solver_group,
                solver_group_ptr.data(),
                solver_group_item.data(),
                (int)solver_group_item.size(),
                solver_group_epsilon.data(),
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

    ThresholdRefitResult threshold_refit_result;
    if (threshold_refit != 0 && weight_tol > 0.0) {
        const int refit_coverage_mode =
            use_hard_coverage_solver ? 2 : 0;
        threshold_refit_result = threshold_consistent_refit(
            A, b, x, m, n, max_iter, tol, weight_tol, refit_coverage_mode,
            solver_group_ptr, solver_group_item, solver_group_epsilon,
            admm_max_iter, admm_rho, admm_abs_tol, admm_rel_tol,
            admm_verbose, threshold_refit_max_pass);

        residual = threshold_refit_result.residual_final;
        if (use_hard_coverage_solver) {
            coverage_status = threshold_refit_result.coverage_status;
            max_coverage_violation = threshold_refit_result.max_coverage_violation;
            admm_outer_iter = threshold_refit_result.coverage_outer_iter;
        }
    }
    else {
        threshold_refit_result.enabled = 0;
        threshold_refit_result.initial_support =
            (int)std::count_if(
                x.begin(),x.end(),
                [weight_tol](double v){ return v > weight_tol; });
        threshold_refit_result.final_support = threshold_refit_result.initial_support;
        threshold_refit_result.residual_before = residual;
        threshold_refit_result.residual_after_initial_cut = residual;
        threshold_refit_result.residual_final = residual;
    }

    time_solver =
        HROM_elapsed_sec(timer_solver_begin, HROMClock::now());

    const auto timer_post_begin = HROMClock::now();

    int uncovered_after_solve = 0;
    double minimum_group_sum = 0.0;
    double minimum_group_margin = 0.0;

    int minimum_margin_group = -1;

    if (use_hard_coverage_solver) {
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


        auto check_operator_group_range = [&] (
            const int begin_group,
            const int end_group,
            const char* label) -> int {

            int operator_uncovered = 0;
            double operator_min_margin = DBL_MAX;
            int operator_worst_solver_g = -1;

            for (int sg = begin_group; sg < end_group; ++sg) {
                double sum = 0.0;
                bool selected_member = false;

                for (int k = solver_group_ptr[(size_t)sg];
                     k < solver_group_ptr[(size_t)sg + 1]; ++k) {
                    const int c = solver_group_item[(size_t)k];
                    sum += x[(size_t)c];
                    if (x[(size_t)c] > weight_tol) selected_member = true;
                }

                const double margin =
                    sum - solver_group_epsilon[(size_t)sg];

                if (margin < operator_min_margin) {
                    operator_min_margin = margin;
                    operator_worst_solver_g = sg;
                }
                if (!selected_member) operator_uncovered++;
            }

            printf(
                "[%s-OPERATOR-COVERAGE-RESULT] subgroups=%d "
                "uncovered_after_threshold=%d min_margin=%.15e "
                "worst_solver_group=%d\n",
                label,
                end_group - begin_group,
                operator_uncovered,
                (operator_min_margin == DBL_MAX) ? 0.0 : operator_min_margin,
                operator_worst_solver_g);

            return operator_uncovered;
        };

        if (coverage_mode == 3 || coverage_mode == 4) {
            const int interface_operator_uncovered =
                check_operator_group_range(
                    num_interface_group,
                    first_internal_solver_group,
                    "INTERFACE");

            if (interface_operator_uncovered != 0) {
                die("operator-aware interface coverage lost one or more signature clusters");
            }
        }

        if (coverage_mode == 4) {
            const int internal_operator_uncovered =
                check_operator_group_range(
                    first_internal_solver_group,
                    num_solver_group,
                    "INTERNAL");

            if (internal_operator_uncovered != 0) {
                die("internal operator-aware coverage lost one or more signature clusters");
            }
        }

        if (coverage_mode == 5 && num_mode5_hard_group > 0) {
            int mode5_uncovered = 0;
            double mode5_min_margin = DBL_MAX;
            for (int sg = 0; sg < num_solver_group; ++sg) {
                double sum = 0.0;
                bool selected_member = false;
                for (int k = solver_group_ptr[(size_t)sg];
                     k < solver_group_ptr[(size_t)sg + 1]; ++k) {
                    const int c = solver_group_item[(size_t)k];
                    sum += x[(size_t)c];
                    if (x[(size_t)c] > weight_tol) selected_member = true;
                }
                mode5_min_margin = std::min(
                    mode5_min_margin,
                    sum - solver_group_epsilon[(size_t)sg]);
                if (!selected_member) mode5_uncovered++;
            }
            printf(
                "[MODE5-HARD-COVERAGE-RESULT] groups=%d uncovered=%d "
                "min_margin=%.15e\n",
                num_mode5_hard_group,
                mode5_uncovered,
                (mode5_min_margin == DBL_MAX) ? 0.0 : mode5_min_margin);
            if (mode5_uncovered != 0) {
                die("mode-5 hard fallback lost one or more adaptive anchors");
            }
        }

        if (coverage_status != 0 ||
            uncovered_after_solve != 0 ||
            max_coverage_violation >
                10.0 * std::max(
                    admm_abs_tol,
                    admm_rel_tol *
                        *std::max_element(
                            solver_group_epsilon.begin(),
                            solver_group_epsilon.end()))) {
            die("coverage-constrained NNLS did not satisfy "
                "the requested hard constraints");
        }
    }
    else {
        printf(
            "[COVERAGE-RESULT] skipped (coverage_mode=%d: %s)\n",
            coverage_mode,
            coverage_mode_name);
    }

    time_postprocess =
        HROM_elapsed_sec(timer_post_begin, HROMClock::now());

    const auto timer_route_begin = HROMClock::now();

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
    if (routes.empty()) {
        routes = read_all_routes(directory,num_ranks);
    }
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
            "internal_coverage_added_candidates %zu\n",
            internal_coverage_added_element.size());
        fprintf(
            fp_report,
            "internal_operator_groups %zu\n",
            internal_operator_groups.size());
        fprintf(
            fp_report,
            "interface_operator_subconstraints %d\n",
            num_interface_operator_subgroups);
        fprintf(
            fp_report,
            "internal_operator_subconstraints %d\n",
            num_internal_operator_subgroups);
        fprintf(
            fp_report,
            "total_solver_coverage_constraints %d\n",
            num_solver_group);
        fprintf(
            fp_report,
            "mode5_seed_tol %.30e\n",
            mode5_seed_tol);
        fprintf(
            fp_report,
            "mode5_defect_tol %.30e\n",
            mode5_defect_tol);
        fprintf(
            fp_report,
            "mode5_rounds %d\n",
            mode5_result.rounds);
        fprintf(
            fp_report,
            "mode5_added_candidates %d\n",
            mode5_result.added_candidates);
        fprintf(
            fp_report,
            "mode5_final_max_defect %.30e\n",
            mode5_result.final_max_defect);
        fprintf(
            fp_report,
            "mode5_final_deficient_groups %d\n",
            mode5_result.deficient_groups_final);
        fprintf(
            fp_report,
            "mode5_hard_fallback_groups %d\n",
            num_mode5_hard_group);
        fprintf(
            fp_report,
            "stage3_candidate_pool %d\n",
            stage3_num_candidates);
        fprintf(
            fp_report,
            "threshold_refit_enabled %d\n",
            threshold_refit_result.enabled);
        fprintf(
            fp_report,
            "threshold_refit_passes %d\n",
            threshold_refit_result.passes);
        fprintf(
            fp_report,
            "threshold_refit_initial_support %d\n",
            threshold_refit_result.initial_support);
        fprintf(
            fp_report,
            "threshold_refit_final_support %d\n",
            threshold_refit_result.final_support);
        fprintf(
            fp_report,
            "threshold_refit_residual_before %.30e\n",
            threshold_refit_result.residual_before);
        fprintf(
            fp_report,
            "threshold_refit_residual_after_initial_cut %.30e\n",
            threshold_refit_result.residual_after_initial_cut);
        fprintf(
            fp_report,
            "threshold_refit_residual_final %.30e\n",
            threshold_refit_result.residual_final);
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

    printf("[Stage 3 v15 result] mode=%d selected=%d / %d "
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

    time_route_stats =
        HROM_elapsed_sec(timer_route_begin, HROMClock::now());

    const auto timer_legacy_begin = HROMClock::now();

    write_parallel_legacy_files(directory,routes,selected);

    time_legacy_write =
        HROM_elapsed_sec(timer_legacy_begin, HROMClock::now());

    const double time_total =
        HROM_elapsed_sec(timer_total_begin, HROMClock::now());

    printf("\n"
           "============================================================\n"
           "[STAGE3-TIMING]\n"
           "read rank files                  = %12.6f s  %6.2f %%\n"
           "read/build interface groups      = %12.6f s  %6.2f %%\n"
           "coverage candidate expansion     = %12.6f s  %6.2f %%\n"
           "candidate/group indexing         = %12.6f s  %6.2f %%\n"
           "global A,b assembly              = %12.6f s  %6.2f %%\n"
           "normalization + anchors          = %12.6f s  %6.2f %%\n"
           "NNLS / coverage solver           = %12.6f s  %6.2f %%\n"
           "coverage post-check              = %12.6f s  %6.2f %%\n"
           "route/stats/report               = %12.6f s  %6.2f %%\n"
           "legacy file write                = %12.6f s  %6.2f %%\n"
           "TOTAL                            = %12.6f s\n"
           "============================================================\n",
           time_read_rank,               100.0*time_read_rank/std::max(time_total,1.0e-300),
           time_interface,               100.0*time_interface/std::max(time_total,1.0e-300),
           time_candidate_expand,        100.0*time_candidate_expand/std::max(time_total,1.0e-300),
           time_index_group,             100.0*time_index_group/std::max(time_total,1.0e-300),
           time_matrix_assembly,         100.0*time_matrix_assembly/std::max(time_total,1.0e-300),
           time_normalize_anchor,        100.0*time_normalize_anchor/std::max(time_total,1.0e-300),
           time_solver,                  100.0*time_solver/std::max(time_total,1.0e-300),
           time_postprocess,             100.0*time_postprocess/std::max(time_total,1.0e-300),
           time_route_stats,             100.0*time_route_stats/std::max(time_total,1.0e-300),
           time_legacy_write,            100.0*time_legacy_write/std::max(time_total,1.0e-300),
           time_total);

    free(Ablock);
    free(A);
    return EXIT_SUCCESS;
}
