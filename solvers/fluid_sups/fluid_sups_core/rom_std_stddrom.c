#include "rom_std_stddrom.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>


/*
 * Rank-local POD uses a rank-local LAPACK SVD.
 *
 * The validated ST-DDROM POD is independent on each spatial MPI rank.  The
 * previous implementation called the MONOLIS ScaLAPACK wrapper on MPI_COMM_SELF,
 * which unnecessarily pulled ScaLAPACK/BLACS into the final executable.
 * Keeping the SVD local removes that link-time dependency without changing the
 * POD algebra or any online reduced operator/solver routine.
 */
#ifdef __cplusplus
extern "C" {
#endif
void dgesvd_(
    char* jobu,
    char* jobvt,
    int* m,
    int* n,
    double* a,
    int* lda,
    double* s,
    double* u,
    int* ldu,
    double* vt,
    int* ldvt,
    double* work,
    int* lwork,
    int* info);
#ifdef __cplusplus
}
#endif

#define ROM_STDD_SNAPSHOT_MAGIC UINT64_C(0x53544444534E5031) /* STDDSNP1 */
#define ROM_STDD_BASIS_MAGIC    UINT64_C(0x5354444442415331) /* STDDBAS1 */
#define ROM_STDD_OPERATOR_MAGIC UINT64_C(0x535444444F505231) /* STDDOPR1 */
#define ROM_STDD_FILE_VERSION   UINT32_C(3)

#include "rom_std_stddrom/system_snapshot_io.inc"
#include "rom_std_stddrom/pod_basis_io.inc"
#include "rom_std_stddrom/operator_projection.inc"
#include "rom_std_stddrom/reduced_solver.inc"
#include "rom_std_stddrom/diagnostics.inc"
