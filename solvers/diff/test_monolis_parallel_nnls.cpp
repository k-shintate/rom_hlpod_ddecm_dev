#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mpi.h>

extern "C" {
#include "monolis_nnls_c.h"
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0, nproc = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    /* One distributed column per MPI rank.
       Global A is identity(nproc), b = [1,...,1]^T.
       Expected NNLS solution is x_i = 1 on every rank. */
    const int m = nproc;
    const int n_loc = 1;
    const int max_iter = 100;
    const double tol_outer = 1.0e-12;
    const double tol_inner = 1.0e-14 * (double)((m > nproc) ? m : nproc);

    std::vector<double> storage((size_t)m * n_loc, 0.0);
    std::vector<double*> A((size_t)m, nullptr);
    for (int i = 0; i < m; ++i) {
        A[(size_t)i] = &storage[(size_t)i * n_loc];
    }
    A[(size_t)rank][0] = 1.0;

    std::vector<double> b((size_t)m, 1.0);
    std::vector<double> x((size_t)n_loc, 0.0);
    double residual = -1.0;

    MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
    int monolis_comm = (int)fcomm;

    std::printf("[SMOKE-BEFORE] rank=%d/%d fcomm=%d m=%d n_loc=%d\n",
                rank, nproc, monolis_comm, m, n_loc);
    std::fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    monolis_optimize_nnls_R_with_sparse_solution(
        A.data(), b.data(), x.data(),
        m, n_loc, max_iter,
        tol_outer, tol_inner,
        &residual, monolis_comm);

    std::printf("[SMOKE-AFTER] rank=%d/%d x=%.17e residual=%.17e\n",
                rank, nproc, x[0], residual);
    std::fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
