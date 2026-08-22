#pragma once

/*
 * Portable linear-algebra compatibility layer.
 *
 * Linux + Intel MKL:
 *   Compile with:
 *
 *       -DFLUID_USE_MKL=1
 *
 * macOS:
 *   Accelerate CBLAS is used.
 *
 * Other platforms:
 *   Standard CBLAS is used.
 */

#include <stddef.h>
#include <stdlib.h>


/* ============================================================
 * Intel MKL
 * ============================================================ */
#if defined(FLUID_USE_MKL) && FLUID_USE_MKL

#include <mkl.h>


/* ============================================================
 * Non-MKL
 * ============================================================ */
#else


/* ------------------------------------------------------------
 * macOS Accelerate
 * ------------------------------------------------------------ */
#if defined(__APPLE__)

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK 1
#endif

#include <Accelerate/Accelerate.h>


/* ------------------------------------------------------------
 * Standard CBLAS
 * ------------------------------------------------------------ */
#else

#include <cblas.h>

#endif


/* ============================================================
 * MKL-compatible aligned allocation for non-MKL environments
 * ============================================================ */
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


#endif