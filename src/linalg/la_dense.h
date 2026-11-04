/*
 * src/linalg/la_dense.h
 *
 * Small shared helpers for the dense matrix heads added together in
 * cholesky.c, hessenberg.c, singularvaluelist.c, minors.c and
 * orthogonalize.c (and Eigensystem's argument checks).
 *
 * Each of those heads has the same three-way dispatch:
 *
 *   - machine numbers (an NDArray / packed buffer, or a boxed matrix whose
 *     inexact leaves are all <= 53 bits)  -> a raw-double LAPACK kernel;
 *   - arbitrary-precision numbers (MPFR leaves, > 53 bits)  -> a generic
 *     kernel run at that precision;
 *   - exact or symbolic entries  -> the same generic kernel over Exprs.
 *
 * The helpers below classify a matrix into one of those classes, read its
 * shape whichever representation it has, unpack it into a boxed nested List
 * when the generic kernel needs one, and build/evaluate small calls.
 *
 * Ownership: every function documents what it borrows and what it returns.
 * Nothing here frees an argument unless its name says `_free`.
 */
#ifndef MATHILDA_LINALG_LA_DENSE_H
#define MATHILDA_LINALG_LA_DENSE_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Precision class of a matrix (or vector) argument. */
typedef enum {
    LA_EXACT   = 0,   /* no inexact leaf: exact numbers and/or symbols       */
    LA_MACHINE = 1,   /* an NDArray, or inexact leaves all <= 53 bits        */
    LA_MPFR    = 2    /* some inexact leaf, and every one > 53 bits          */
} LaPrecClass;

/* Classify `e`.  For LA_MPFR, *bits receives the minimum precision over the
 * inexact leaves (the precision the answer keeps); otherwise 53. */
LaPrecClass la_prec_class(const Expr* e, long* bits);

/* Shape of a rank-2 argument: an NDArray of rank 2, or a List of equal-length
 * Lists.  Returns true with *rows / *cols set (both may be 0 for {{}}); false
 * for anything else (a vector, a ragged list, an atom). */
bool la_matrix_dims(const Expr* e, int* rows, int* cols);

/* True when `e` is a List each of whose elements is a List (including the empty
 * List {}), or an NDArray of rank >= 2 -- i.e. something that LOOKS like a
 * matrix argument, so a malformed one deserves the ::matsq / ::matrix message
 * rather than a silent decline. */
bool la_looks_like_matrix(const Expr* e);

/* A boxed nested-List copy of `e` (an NDArray is unpacked, a List is copied).
 * Caller owns the result. */
Expr* la_boxed_copy(const Expr* e);

/* Entry (i, j) of a boxed List-of-Lists matrix (borrowed). */
static inline Expr* la_at(const Expr* m, int i, int j) {
    return m->data.function.args[i]->data.function.args[j];
}

/* True when every element of the boxed matrix (or vector) `e` is a numeric
 * quantity in the NumericQ sense -- no free symbols. */
bool la_all_numeric(const Expr* e);

/* True when any leaf of `e` is a Complex number. */
bool la_has_complex(const Expr* e);

/* Build head[a] / head[a, b] / head[a, b, c].  The arguments are ADOPTED
 * (ownership passes to the new call). */
Expr* la_call1(const char* head, Expr* a);
Expr* la_call2(const char* head, Expr* a, Expr* b);
Expr* la_call3(const char* head, Expr* a, Expr* b, Expr* c);

/* Evaluate head[a] / head[a, b], adopting the arguments; caller owns result. */
Expr* la_ev1(const char* head, Expr* a);
Expr* la_ev2(const char* head, Expr* a, Expr* b);

/* An empty List. */
Expr* la_empty_list(void);

/* A rows x cols boxed matrix built from a row-major array of adopted Exprs
 * (the array itself is borrowed). */
Expr* la_matrix_from(Expr** cells, int rows, int cols);

/* The printed form of `e` (malloc'd, caller frees) -- for message text. */
char* la_str(const Expr* e);

/* Parse the standard "take" specification used by Eigenvalues and
 * SingularValueList as a positional argument: an Integer k (k >= 0 the first
 * k, k < 0 the last |k|) or UpTo[k].  On success sets *k and *upto and returns
 * true; false for anything else. */
bool la_parse_take_spec(const Expr* spec, long long* k, bool* upto);

/* ---- direct machine kernels, for Compile[] -------------------------------
 * Each takes an NDArray (borrowed) and returns the same packed result the
 * interpreter's machine path returns for it, or NULL WITHOUT a message on
 * anything it does not answer (the VM then aborts to the interpreter, which
 * reports it once). */
Expr* la_cholesky_machine(const Expr* m);            /* cholesky.c */
Expr* la_singularvaluelist_machine(const Expr* m);    /* singularvaluelist.c:
                                                       * the default-tolerance
                                                       * nonzero values (real) */
Expr* la_minors_machine(const Expr* m);               /* minors.c: Minors[m]   */
Expr* la_adjugate_machine(const Expr* m);             /* minors.c: Adjugate[m] */
Expr* la_orthogonalize_machine(const Expr* m);        /* orthogonalize.c: the
                                                       * Automatic method     */
Expr* la_hessenberg_machine(const Expr* m);           /* hessenberg.c: {p, h}
                                                       * as ONE rank-3
                                                       * 2 x n x n array      */

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_LA_DENSE_H */
