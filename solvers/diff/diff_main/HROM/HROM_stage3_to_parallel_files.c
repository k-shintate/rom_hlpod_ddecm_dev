/*
 * Convert serial Stage-3 result back to the original per-fine-subdomain
 * DDECM files used by the parallel HROM online calculation.
 *
 * IMPORTANT FIX (v6):
 *   A Stage-3 column is the adopted MPI-rank contribution, which is assembled
 *   from all fine subdomains owned by that MPI rank.  Therefore a selected
 *   physical element must be written to EVERY fine subdomain of the adopted
 *   rank that contains that element, not only to the owner fine subdomain.
 *
 * Inputs:
 *   DDECM/stage3_selected_elem.txt
 *   metagraph_parted.0/metagraph.dat.n_internal.<rank>
 *   metagraph_parted.0/metagraph.dat.id.<rank>
 *   parted.1/elem.dat.n_internal.<fine_subdomain>
 *   parted.1/elem.dat.id.<fine_subdomain>
 *   DDECM/lb_selected_elem*.txt.stage2.bak (preferred classification source)
 *
 * Outputs (overwrite current files):
 *   DDECM/lb_selected_elem_D_bc.<fine_subdomain>.txt
 *   DDECM/lb_selected_elem.<fine_subdomain>.txt
 */
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static std::string join_path(const char* directory, const std::string& relative)
{
    std::string base = directory ? directory : ".";
    if (!base.empty() && base.back() != '/') base.push_back('/');
    return base + relative;
}

static bool file_exists(const std::string& path)
{
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return false;
    fclose(fp);
    return true;
}

static void die(const char* msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}

static int read_n_internal_value(
    const char* directory,
    const std::string& relative)
{
    const std::string path = join_path(directory, relative);
    FILE* fp = fopen(path.c_str(), "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    char label[256];
    int tmp = 0;
    int value = -1;
    if (fscanf(fp, "%255s %d", label, &tmp) != 2 ||
        fscanf(fp, "%d", &value) != 1 || value < 0) {
        fclose(fp);
        fprintf(stderr, "ERROR: invalid n_internal file %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    return value;
}

static std::vector<int64_t> read_id_file(
    const char* directory,
    const std::string& relative)
{
    const std::string path = join_path(directory, relative);
    FILE* fp = fopen(path.c_str(), "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    char label[256];
    int n = -1;
    int tmp = 0;
    if (fscanf(fp, "%255s", label) != 1 ||
        fscanf(fp, "%d %d", &n, &tmp) != 2 || n < 0) {
        fclose(fp);
        fprintf(stderr, "ERROR: invalid id-file header %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    std::vector<int64_t> ids((size_t)n, -1);
    for (int i = 0; i < n; ++i) {
        long long v = -1;
        if (fscanf(fp, "%lld", &v) != 1) {
            fclose(fp);
            fprintf(stderr, "ERROR: failed to read ID %d/%d from %s\n",
                i, n, path.c_str());
            exit(EXIT_FAILURE);
        }
        ids[(size_t)i] = (int64_t)v;
    }

    fclose(fp);
    return ids;
}

struct PartitionInfo {
    std::vector<std::vector<int>> rank_fine_subdomains;
    std::vector<int> all_fine_subdomains;
    std::unordered_map<int64_t, int> adopted_rank;
    std::unordered_map<int, std::unordered_set<int64_t>> fine_membership;
};

static PartitionInfo build_partition_info(
    const char* directory,
    const int num_ranks,
    const char* metagraph_dir,
    const char* fine_partition_dir)
{
    PartitionInfo out;
    out.rank_fine_subdomains.resize((size_t)num_ranks);

    std::unordered_set<int> seen_fine;

    for (int rank = 0; rank < num_ranks; ++rank) {
        char meta_nint[512];
        char meta_id[512];
        snprintf(meta_nint, sizeof(meta_nint),
            "%s/metagraph.dat.n_internal.%d", metagraph_dir, rank);
        snprintf(meta_id, sizeof(meta_id),
            "%s/metagraph.dat.id.%d", metagraph_dir, rank);

        const int n_owned_sub =
            read_n_internal_value(directory, meta_nint);
        const std::vector<int64_t> sub_ids =
            read_id_file(directory, meta_id);

        if (n_owned_sub > (int)sub_ids.size()) {
            fprintf(stderr,
                "ERROR: rank=%d metagraph n_internal=%d exceeds id count=%zu\n",
                rank, n_owned_sub, sub_ids.size());
            exit(EXIT_FAILURE);
        }

        printf("[partition] rank=%d fine subdomains:", rank);

        for (int k = 0; k < n_owned_sub; ++k) {
            const int s = (int)sub_ids[(size_t)k];
            printf(" %d", s);
            out.rank_fine_subdomains[(size_t)rank].push_back(s);
            if (seen_fine.insert(s).second) {
                out.all_fine_subdomains.push_back(s);
            }

            char elem_nint[512];
            char elem_id[512];
            snprintf(elem_nint, sizeof(elem_nint),
                "%s/elem.dat.n_internal.%d", fine_partition_dir, s);
            snprintf(elem_id, sizeof(elem_id),
                "%s/elem.dat.id.%d", fine_partition_dir, s);

            const int n_owner_elem =
                read_n_internal_value(directory, elem_nint);
            const std::vector<int64_t> ids =
                read_id_file(directory, elem_id);

            if (n_owner_elem > (int)ids.size()) {
                fprintf(stderr,
                    "\nERROR: subdomain=%d n_internal=%d exceeds id count=%zu\n",
                    s, n_owner_elem, ids.size());
                exit(EXIT_FAILURE);
            }

            /* Full membership is needed when restoring online files. */
            std::unordered_set<int64_t> members;
            members.reserve(ids.size() * 2 + 1);
            for (int64_t gid : ids) members.insert(gid);
            out.fine_membership.emplace(s, std::move(members));

            /* First n_internal entries define the unique adopted rank. */
            for (int e = 0; e < n_owner_elem; ++e) {
                const int64_t gid = ids[(size_t)e];
                const auto ins = out.adopted_rank.emplace(gid, rank);
                if (!ins.second && ins.first->second != rank) {
                    fprintf(stderr,
                        "\nERROR: original global element %lld has multiple adopted ranks: %d and %d\n",
                        (long long)gid, ins.first->second, rank);
                    exit(EXIT_FAILURE);
                }
            }
        }
        printf("\n");
    }

    std::sort(out.all_fine_subdomains.begin(), out.all_fine_subdomains.end());

    printf("[partition] adopted elements=%zu fine subdomains=%zu\n",
        out.adopted_rank.size(), out.all_fine_subdomains.size());

    return out;
}

static std::string selection_path(
    const char* directory,
    const int s,
    const bool dbc,
    const bool prefer_stage2_backup)
{
    char rel[512];
    snprintf(rel, sizeof(rel),
        dbc ? "DDECM/lb_selected_elem_D_bc.%d.txt"
            : "DDECM/lb_selected_elem.%d.txt", s);

    const std::string current = join_path(directory, rel);
    const std::string backup = current + ".stage2.bak";

    if (prefer_stage2_backup && file_exists(backup)) {
        return backup;
    }
    return current;
}

static void read_class_file(
    const char* directory,
    const int s,
    const bool dbc,
    std::unordered_map<int64_t, int>& cls)
{
    const std::string path =
        selection_path(directory, s, dbc, true);

    FILE* fp = fopen(path.c_str(), "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    int n = -1;
    if (fscanf(fp, "%d", &n) != 1 || n < 0) {
        fclose(fp);
        fprintf(stderr, "ERROR: invalid count in %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; ++i) {
        long long gid = -1;
        double old_weight = 0.0;
        if (fscanf(fp, "%lld %lf", &gid, &old_weight) != 2) {
            fclose(fp);
            fprintf(stderr, "ERROR: invalid row in %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }
        (void)old_weight;

        const int value = dbc ? 1 : 0;
        const auto it = cls.find((int64_t)gid);
        if (it == cls.end()) {
            cls.emplace((int64_t)gid, value);
        }
        else if (it->second != value) {
            /* D_bc is conservative if inconsistent old files exist. */
            it->second = 1;
            fprintf(stderr,
                "WARNING: element %lld appears in both D_bc classes; using D_bc.\n",
                gid);
        }
    }

    fclose(fp);
}

static std::vector<std::pair<int64_t,double>> read_stage3_selected(
    const char* directory)
{
    const std::string path =
        join_path(directory, "DDECM/stage3_selected_elem.txt");

    FILE* fp = fopen(path.c_str(), "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    int n = -1;
    if (fscanf(fp, "%d", &n) != 1 || n < 0) {
        fclose(fp);
        die("invalid stage3_selected_elem header");
    }

    std::vector<std::pair<int64_t,double>> selected;
    selected.reserve((size_t)n);

    for (int i = 0; i < n; ++i) {
        long long gid = -1;
        double w = 0.0;
        if (fscanf(fp, "%lld %lf", &gid, &w) != 2) {
            fclose(fp);
            die("invalid stage3_selected_elem row");
        }
        selected.emplace_back((int64_t)gid, w);
    }

    fclose(fp);
    return selected;
}

static void copy_file_binary_if_missing(
    const std::string& src,
    const std::string& dst)
{
    if (file_exists(dst)) return;

    FILE* in = fopen(src.c_str(), "rb");
    if (!in) {
        fprintf(stderr, "ERROR: cannot open %s\n", src.c_str());
        exit(EXIT_FAILURE);
    }

    FILE* out = fopen(dst.c_str(), "wb");
    if (!out) {
        fclose(in);
        fprintf(stderr, "ERROR: cannot open %s\n", dst.c_str());
        exit(EXIT_FAILURE);
    }

    char buf[65536];
    while (true) {
        const size_t n = fread(buf, 1, sizeof(buf), in);
        if (n > 0 && fwrite(buf, 1, n, out) != n) {
            fclose(in); fclose(out);
            die("backup write failed");
        }
        if (n < sizeof(buf)) {
            if (ferror(in)) {
                fclose(in); fclose(out);
                die("backup read failed");
            }
            break;
        }
    }

    fclose(in);
    fclose(out);
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr,
            "usage: %s <directory> <num_mpi_ranks> "
            "[fine_partition_dir=parted.1] "
            "[metagraph_dir=metagraph_parted.0]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    const char* directory = argv[1];
    const int num_ranks = atoi(argv[2]);
    const char* fine_partition_dir =
        argc >= 4 ? argv[3] : "parted.1";
    const char* metagraph_dir =
        argc >= 5 ? argv[4] : "metagraph_parted.0";

    if (num_ranks <= 0) die("invalid number of MPI ranks");

    const PartitionInfo part = build_partition_info(
        directory,
        num_ranks,
        metagraph_dir,
        fine_partition_dir);

    /* Use original Stage-2 files (backup preferred) only for D_bc class. */
    std::unordered_map<int64_t, int> classification;
    for (int s : part.all_fine_subdomains) {
        read_class_file(directory, s, true, classification);
        read_class_file(directory, s, false, classification);
    }

    const auto selected = read_stage3_selected(directory);

    /* Validate before touching current online files. */
    for (const auto& p : selected) {
        const int64_t elem = p.first;
        const auto ar = part.adopted_rank.find(elem);
        if (ar == part.adopted_rank.end()) {
            fprintf(stderr,
                "ERROR: Stage-3 element %lld has no adopted MPI rank\n",
                (long long)elem);
            return EXIT_FAILURE;
        }
        if (classification.find(elem) == classification.end()) {
            fprintf(stderr,
                "ERROR: Stage-3 element %lld is absent from original Stage-2 "
                "selection files, so D_bc class cannot be recovered safely.\n",
                (long long)elem);
            return EXIT_FAILURE;
        }

        const int rank = ar->second;
        int copies_on_adopted_rank = 0;
        for (int s : part.rank_fine_subdomains[(size_t)rank]) {
            const auto mit = part.fine_membership.find(s);
            if (mit != part.fine_membership.end() &&
                mit->second.find(elem) != mit->second.end()) {
                copies_on_adopted_rank++;
            }
        }

        if (copies_on_adopted_rank == 0) {
            fprintf(stderr,
                "ERROR: Stage-3 element %lld adopted by rank %d, but it is not "
                "contained in any fine subdomain of that rank.\n",
                (long long)elem, rank);
            return EXIT_FAILURE;
        }
    }

    using EW = std::pair<int64_t,double>;
    std::unordered_map<int, std::vector<EW>> dbc_out;
    std::unordered_map<int, std::vector<EW>> nodbc_out;

    for (int s : part.all_fine_subdomains) {
        dbc_out[s] = {};
        nodbc_out[s] = {};
    }

    int distributed_copies = 0;
    int multi_subdomain_elements = 0;

    for (const auto& p : selected) {
        const int64_t elem = p.first;
        const double weight = p.second;
        const int rank = part.adopted_rank.at(elem);
        const int cls = classification.at(elem);

        int ncopy = 0;

        /*
         * Critical: reproduce the adopted rank column by activating every
         * fine-subdomain copy of this physical element on the adopted rank.
         */
        for (int s : part.rank_fine_subdomains[(size_t)rank]) {
            const auto& members = part.fine_membership.at(s);
            if (members.find(elem) == members.end()) continue;

            if (cls == 1) {
                dbc_out[s].emplace_back(elem, weight);
            }
            else {
                nodbc_out[s].emplace_back(elem, weight);
            }
            ncopy++;
            distributed_copies++;
        }

        if (ncopy > 1) multi_subdomain_elements++;

        printf("[distribute] elem=%lld weight=%.15e adopted_rank=%d copies=%d\n",
            (long long)elem, weight, rank, ncopy);
    }

    /* Write after all validation/distribution is complete. */
    for (int s : part.all_fine_subdomains) {
        auto& a = dbc_out[s];
        auto& b = nodbc_out[s];

        auto cmp = [](const EW& x, const EW& y) {
            return x.first < y.first;
        };
        std::sort(a.begin(), a.end(), cmp);
        std::sort(b.begin(), b.end(), cmp);

        char rel_a[512];
        char rel_b[512];
        snprintf(rel_a, sizeof(rel_a),
            "DDECM/lb_selected_elem_D_bc.%d.txt", s);
        snprintf(rel_b, sizeof(rel_b),
            "DDECM/lb_selected_elem.%d.txt", s);

        const std::string pa = join_path(directory, rel_a);
        const std::string pb = join_path(directory, rel_b);

        /* Preserve original Stage-2 files if they are not preserved yet. */
        copy_file_binary_if_missing(pa, pa + ".stage2.bak");
        copy_file_binary_if_missing(pb, pb + ".stage2.bak");

        FILE* fa = fopen(pa.c_str(), "w");
        FILE* fb = fopen(pb.c_str(), "w");
        if (!fa || !fb) {
            if (fa) fclose(fa);
            if (fb) fclose(fb);
            die("cannot open output selection files");
        }

        fprintf(fa, "%zu\n", a.size());
        for (const auto& q : a) {
            fprintf(fa, "%lld %.30e\n",
                (long long)q.first, q.second);
        }

        fprintf(fb, "%zu\n", b.size());
        for (const auto& q : b) {
            fprintf(fb, "%lld %.30e\n",
                (long long)q.first, q.second);
        }

        fclose(fa);
        fclose(fb);

        printf("[write] subdomain=%d selected=%zu D_bc=%zu no_D_bc=%zu\n",
            s, a.size() + b.size(), a.size(), b.size());
    }

    printf("[done] unique Stage-3 selected=%zu distributed copies=%d "
           "elements spanning >1 fine subdomain on adopted rank=%d\n",
        selected.size(), distributed_copies, multi_subdomain_elements);

    if (distributed_copies < (int)selected.size()) {
        die("distributed copy count is smaller than unique Stage-3 selection count");
    }

    return EXIT_SUCCESS;
}