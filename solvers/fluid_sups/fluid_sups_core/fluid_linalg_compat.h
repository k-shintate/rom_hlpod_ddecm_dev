#pragma once

/*
 * Portable BLAS / MKL compatibility layer
 *
 * IMPORTANT:
 *   The fluid sources are compiled with mpic++, even though the source
 *   extension is .c.  Therefore CBLAS declarations must have C linkage.
 *
 * Supported:
 *   Intel oneAPI / MKL: <mkl.h> or <mkl/mkl.h>
 *   macOS:              Accelerate
 *   Generic Linux:      CBLAS
 */

#include <stddef.h>
#include <stdlib.h>

/* =========================================================================
 * Intel MKL backend
 * ========================================================================= */

#if defined(FLUID_USE_MKL) && FLUID_USE_MKL

#if defined(__has_include)

#  if __has_include(<mkl.h>)
#    include <mkl.h>
#    define FLUID_HAVE_MKL 1
#  elif __has_include(<mkl/mkl.h>)
#    include <mkl/mkl.h>
#    define FLUID_HAVE_MKL 1
#  else
#    error "FLUID_USE_MKL=1 but MKL headers were not found"
#  endif

#else

#  include <mkl/mkl.h>
#  define FLUID_HAVE_MKL 1

#endif

/* =========================================================================
 * Non-MKL header backend
 * ========================================================================= */

#else

#define FLUID_HAVE_MKL 0

#if defined(__APPLE__)

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK 1
#endif

#include <Accelerate/Accelerate.h>

#else

/*
 * The project compiles .c files with mpic++, so without this wrapper a
 * CBLAS header that lacks its own __cplusplus guard declares cblas_* with
 * C++ linkage.  The library (MKL/OpenBLAS/Netlib CBLAS) exports C symbols,
 * producing linker errors such as:
 *
 *   undefined reference to
 *   cblas_ddot(int, double const*, int, double const*, int)
 *
 * Force C linkage explicitly.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <cblas.h>

#ifdef __cplusplus
}
#endif

#endif  /* __APPLE__ */

/* -------------------------------------------------------------------------
 * MKL-compatible aligned allocation for non-MKL systems
 * ------------------------------------------------------------------------- */

static inline void* fluid_aligned_malloc(
    size_t size,
    size_t alignment)
{
    void* pointer = NULL;

    if(size == 0){
        return NULL;
    }

    if(alignment < sizeof(void*)){
        alignment = sizeof(void*);
    }

    if(posix_memalign(&pointer, alignment, size) != 0){
        return NULL;
    }

    return pointer;
}

static inline void fluid_aligned_free(
    void* pointer)
{
    free(pointer);
}

#define mkl_malloc(size, alignment) \
    fluid_aligned_malloc((size), (alignment))

#define mkl_free(pointer) \
    fluid_aligned_free((pointer))

#endif  /* FLUID_USE_MKL */