#ifndef PAD_SCHEMES_H
#define PAD_SCHEMES_H

#include "expr.h"
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* Shared 1-D padding-scheme engine for ArrayPad and ArrayReshape.
 *
 * A "fiber" is a sequence of borrowed Expr* (leaves for ArrayReshape's flat
 * element list, or whole sub-arrays for one axis of a nested ArrayPad).  The
 * engine extends a fiber by `front` elements before and `back` elements after,
 * according to a padding scheme, producing a fresh owned array of
 * (front + L + back) Expr* copies.  ArrayReshape uses front == 0. */

typedef enum {
    PS_CONSTANT,        /* a scalar c (default 0): fill with copies of c        */
    PS_CYCLIC,          /* a List {c1,c2,...}: cyclic repetition of constants   */
    PS_FIXED,           /* "Fixed": repeat the boundary element                 */
    PS_PERIODIC,        /* "Periodic": cyclic repetition of the whole fiber     */
    PS_REFLECTED,       /* "Reflected": reflect, boundary not repeated (L>=2)   */
    PS_REVERSED,        /* "Reversed": reflect, boundary repeated               */
    PS_REVERSEDNEG,     /* "ReversedNegation": negated reversals                */
    PS_REFLECTEDDIFF,   /* "ReflectedDifferences": antisym reflection (L>=2)    */
    PS_REVERSEDDIFF,    /* "ReversedDifferences": antisym reversal (L>=2)       */
    PS_EXTRAPOLATED,    /* "Extrapolated": polynomial extrapolation            */
    PS_INVALID          /* unrecognised scheme argument                        */
} PadScheme;

/* Classify a padding-scheme argument.  `s` may be NULL (the default: PS_CONSTANT
 * with value 0), a scalar (PS_CONSTANT), a List (PS_CYCLIC), or a String/Symbol
 * naming one of the schemes.  Returns PS_INVALID for anything else. */
PadScheme pad_scheme_classify(const Expr* s);

/* True for the schemes whose padding is derived from the fiber's own values
 * (everything except PS_CONSTANT and PS_CYCLIC). */
bool pad_scheme_value_dependent(PadScheme sc);

/* True for the schemes that require a padded axis of length >= 2 (the
 * antisymmetric difference schemes and Reflected).  Used to raise
 * ArrayPad::mindimsize. */
bool pad_scheme_needs_two(PadScheme sc);

/* Extend fiber v[0..L-1] by `front` before and `back` after under scheme `sc`.
 * `scheme` is the original argument (needed for the constant/cyclic block);
 * `order` is the InterpolationOrder for PS_EXTRAPOLATED (< 0 means Infinity =
 * degree L-1).  On success writes a fresh malloc'd array of
 * (front + L + back) owned Expr* to *out and its length to *out_len, and
 * returns true.  Returns false (writing nothing) when the scheme cannot be
 * applied (e.g. a value-dependent scheme on an empty/too-short fiber). */
bool pad_scheme_extend(const Expr** v, size_t L, int64_t front, int64_t back,
                       const Expr* scheme, PadScheme sc, int64_t order,
                       Expr*** out, size_t* out_len);

#endif /* PAD_SCHEMES_H */
