#!/bin/bash
set -Eeuo pipefail

repo_root="${1:-$(pwd)}"
src="${repo_root}/solvers/fluid_sups/fluid_sups_core/core_NR.c"

if [ ! -f "${src}" ]; then
    echo "ERROR: not found: ${src}" >&2
    exit 2
fi

echo "============================================================"
echo "core_NR linear-algebra dependency diagnosis"
echo "source: ${src}"
echo "============================================================"

echo
echo "[1] MKL include locations"
grep -nE '^[[:space:]]*#[[:space:]]*include[[:space:]]*[<"]mkl\.h[>"]' "${src}" || true

echo
echo "[2] MKL-specific identifiers"
grep -nE '\b(MKL_[A-Za-z0-9_]*|mkl_[A-Za-z0-9_]*)\b' "${src}" || true

echo
echo "[3] LAPACKE calls"
grep -nE '\bLAPACKE_[A-Za-z0-9_]+\b' "${src}" || true

echo
echo "[4] CBLAS calls"
grep -nE '\bcblas_[A-Za-z0-9_]+\b' "${src}" || true

echo
echo "[5] likely Fortran BLAS/LAPACK calls"
grep -nE '\b[sdcz](gesv|gesvd|gesdd|getrf|getri|potrf|potrs|gemm|gemv|axpy|dot|scal)_?\b' "${src}" || true

echo
echo "[6] source around current failure"
nl -ba "${src}" | sed -n '1800,1900p'

echo
echo "============================================================"
echo "Interpretation"
echo "============================================================"
echo "A) Only '#include <mkl.h>' appears:"
echo "   The include is stale. Remove it."
echo
echo "B) cblas_* appears, but no MKL_/mkl_/LAPACKE_:"
echo "   Replace mkl.h by cblas.h and link a CBLAS provider."
echo
echo "C) LAPACKE_* appears:"
echo "   Either provide LAPACKE explicitly or convert those calls to the"
echo "   project's existing Fortran LAPACK ABI (-llapack -lblas)."
echo
echo "D) MKL_* or mkl_* appears:"
echo "   This is a real MKL-specific dependency. Keep an optional USE_MKL path"
echo "   and implement a portable fallback before requiring the full stack."
