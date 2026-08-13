#pragma once

/*
 * Portable linear-algebra compatibility layer for fluid_sups_core/core_NR.c.
 *
 * macOS:
 *   - Accelerate CBLAS
 *   - LP64 interface (32-bit BLAS/LAPACK integer arguments)
 *   - posix_memalign for the existing 64-byte aligned work buffers
 *
 * IMPORTANT:
 * ACCELERATE_NEW_LAPACK should preferably be supplied from the compiler
 * command line so it is defined before ANY Accelerate/vecLib header can be
 * included indirectly.  This header also defines it defensively.
 *
 * We intentionally DO NOT define ACCELERATE_LAPACK_ILP64 because the current
 * core_NR.c passes ordinary int values (nrv, msz, etc.) to CBLAS.
 */

#include <stddef.h>
#include <stdlib.h>

#if defined(FLUID_USE_MKL) && FLUID_USE_MKL

#include <mkl.h>

#else

#if defined(__APPLE__)

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK 1
#endif

#include <Accelerate/Accelerate.h>

#else
#include <cblas.h>
#endif

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

static inline void fluid_aligned_free(void* pointer)
{
    free(pointer);
}

#define mkl_malloc(size, alignment) \
    fluid_aligned_malloc((size), (alignment))

#define mkl_free(pointer) \
    fluid_aligned_free((pointer))

#endif
