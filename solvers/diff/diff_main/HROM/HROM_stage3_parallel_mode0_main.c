#include <cstdio>
#include <cstdlib>

#include <mpi.h>

#include "HROM_stage3_parallel_mode0.h"

static void usage(const char* prog)
{
    std::fprintf(stderr,
        "Usage:\n"
        "  mpirun -np <solver_ranks> %s "
        "[directory] [input_rank_files] [max_iter] "
        "[tol_outer] [tol_inner] [weight_tol] [output_relative]\n\n"
        "Defaults:\n"
        "  directory        = .\n"
        "  input_rank_files = 0   (auto-detect stage3_rank.0.bin, .1.bin, ...)\n"
        "  max_iter         = 5000\n"
        "  tol_outer        = 1e-10\n"
        "  tol_inner        = 0   (auto: 1e-14*max(rows,candidates))\n"
        "  weight_tol       = 0\n"
        "  output_relative  = (empty: disabled; legacy files are authoritative)\n",
        prog);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 8) {
        if (rank == 0) usage(argv[0]);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    const char* directory = (argc > 1) ? argv[1] : ".";
    const int input_rank_files = (argc > 2) ? std::atoi(argv[2]) : 0;
    const int max_iter = (argc > 3) ? std::atoi(argv[3]) : 5000;
    const double tol_outer = (argc > 4) ? std::atof(argv[4]) : 1.0e-10;
    const double tol_inner = (argc > 5) ? std::atof(argv[5]) : 0.0;
    const double weight_tol = (argc > 6) ? std::atof(argv[6]) : 0.0;
    const char* output_relative =
        (argc > 7) ? argv[7] : "";

    const int rc = HROM_stage3_parallel_mode0_run(
        MPI_COMM_WORLD,
        directory,
        input_rank_files,
        max_iter,
        tol_outer,
        tol_inner,
        weight_tol,
        output_relative);

    MPI_Finalize();
    return rc;
}
