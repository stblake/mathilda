/*
 * version.c — compile-time assembly of the $Version string.
 *
 * Every component version is available as a preprocessor macro at compile
 * time, so the whole descriptive string is built by concatenating adjacent
 * string literals. Optional components are guarded by the exact build defines
 * the makefile emits (USE_MPFR, USE_FLINT, NO_ECM, USE_GRAPHICS, USE_LAPACK /
 * MATHILDA_USE_ACCELERATE, NO_READLINE), so the string only ever names the
 * libraries that were actually linked. No runtime work, no generated files.
 */

#include "version.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#ifdef USE_FLINT
#include <flint/flint.h>
#endif
#ifndef NO_ECM
#include <ecm.h>
#endif
#ifdef USE_GRAPHICS
#include <raylib.h>
#endif

/* Two-level stringify so macro *values* (not names) get quoted. */
#define MVER_STR2(x) #x
#define MVER_STR(x)  MVER_STR2(x)

/* --- C compiler (leads the list, no separator) ------------------------- */
/* __VERSION__ alone is inconsistent across compilers: clang and Apple clang
 * expand it to a string that already begins with the vendor name ("clang ...",
 * "Apple LLVM ..."), but real GCC expands it to a bare version number like
 * "16.1.0" with no "GCC" prefix. Identify the compiler explicitly so the
 * $Version string always names it. Order matters: clang also defines __GNUC__,
 * so test __clang__ first. */
#if defined(__clang__)
#  define VER_COMPILER __VERSION__      /* already begins with "clang"/"Apple LLVM" */
#elif defined(__GNUC__)
#  define VER_COMPILER "GCC " MVER_STR(__GNUC__) "." \
                              MVER_STR(__GNUC_MINOR__) "." \
                              MVER_STR(__GNUC_PATCHLEVEL__)
#elif defined(__VERSION__)
#  define VER_COMPILER __VERSION__
#else
#  define VER_COMPILER "unknown compiler"
#endif

/* --- GMP (always present; exposes only integer version macros) --------- */
#define VER_GMP ", GMP " MVER_STR(__GNU_MP_VERSION) "." \
                         MVER_STR(__GNU_MP_VERSION_MINOR) "." \
                         MVER_STR(__GNU_MP_VERSION_PATCHLEVEL)

/* --- MPFR -------------------------------------------------------------- */
#ifdef USE_MPFR
#  define VER_MPFR ", MPFR " MPFR_VERSION_STRING
#else
#  define VER_MPFR ""
#endif

/* --- FLINT ------------------------------------------------------------- */
#ifdef USE_FLINT
#  define VER_FLINT ", FLINT " FLINT_VERSION
#else
#  define VER_FLINT ""
#endif

/* --- GMP-ECM (version macro if the header provides one) ---------------- */
#ifndef NO_ECM
#  ifdef ECM_VERSION
#    define VER_ECM ", ECM " ECM_VERSION
#  else
#    define VER_ECM ", ECM"
#  endif
#else
#  define VER_ECM ""
#endif

/* --- Raylib (graphics backend) ----------------------------------------- */
#ifdef USE_GRAPHICS
#  ifdef RAYLIB_VERSION
#    define VER_GRAPHICS ", Raylib " RAYLIB_VERSION
#  else
#    define VER_GRAPHICS ", Raylib"
#  endif
#else
#  define VER_GRAPHICS ""
#endif

/* --- Dense LA backend (no version macro; report the backend name) ------ */
#ifdef MATHILDA_USE_ACCELERATE
#  define VER_LAPACK ", Accelerate"
#elif defined(USE_LAPACK)
#  define VER_LAPACK ", LAPACK"
#else
#  define VER_LAPACK ""
#endif

/* --- GNU Readline (no portable compile-time version macro; report name) */
#ifndef NO_READLINE
#  define VER_READLINE ", Readline"
#else
#  define VER_READLINE ""
#endif

static const char kVersionString[] =
    "Mathilda " MATHILDA_VERSION_STRING " ("
    VER_COMPILER
    VER_GMP
    VER_MPFR
    VER_FLINT
    VER_ECM
    VER_GRAPHICS
    VER_LAPACK
    VER_READLINE
    ")";

const char* mathilda_version(void) {
    return kVersionString;
}
