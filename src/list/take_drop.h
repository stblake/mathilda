#ifndef TAKE_DROP_H
#define TAKE_DROP_H

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>

Expr* builtin_take(Expr* res);
Expr* builtin_drop(Expr* res);

/* Resolve a Take/Drop-style sequence spec against a length `len` into a set of
 * 1-based indices (allocated in *out_indices, count in *out_count; the caller
 * frees *out_indices). Handles All, None, an Integer k (first/last |k|), UpTo[k]
 * (clamped), and {m}, {m,n}, {m,n,s} lists with negative-index normalization and
 * a step. Returns false — and allocates nothing — for a malformed or
 * out-of-range spec. Shared with Ordering[list, seq], which selects positions of
 * the sort permutation the same way Take selects list elements. */
bool get_seq_spec_indices(Expr* spec, int64_t len, int64_t** out_indices, size_t* out_count);

#endif /* TAKE_DROP_H */
