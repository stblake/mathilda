#ifndef MATHILDA_NDARRAY_H
#define MATHILDA_NDARRAY_H

/* ---------------------------------------------------------------------------
 * NDArray — a first-class, visibly-distinct dense machine-precision ndarray.
 *
 * Modelled on numpy's ndarray: an N-dimensional (rank 1..NDARRAY_MAX_RANK),
 * rectangular, C-order (row-major) dense array of machine-precision values.
 * Unlike Mathematica's packed arrays, which are an invisible internal
 * optimization a list may or may not have, `NDArray[...]` is always what it
 * says it is: `Head[NDArray[{{1,2},{3,4}}]]` is `NDArray`, never `List`, and
 * `ArrayQ`/`ListQ` never disagree with each other about which one a value is.
 * It is a real `EXPR_NDARRAY` node (see expr.h's `NDArrayData`), storing a
 * flat row-major `double*` buffer plus rank/dims directly in the `Expr` —
 * not a nested `List[List[...]]` tree — so builtins that recognize it can
 * operate on the buffer directly instead of paying the flatten/rebuild cost
 * every List-based tensor operation pays today (see src/linalg/dot.c).
 *
 * numpy correspondence: `Dimensions` is `shape`, `ArrayDepth` is `ndim`,
 * `Length` is `len` (= shape[0]). Rank is unbounded up to NDARRAY_MAX_RANK.
 *
 * Machine precision only. Any operation that would need to store a
 * non-double entry (a symbol, an exact/rational number, a BigInt, an MPFR
 * value) degrades to an ordinary nested List instead of forcing a lossy
 * conversion — see ndarray_from_nested_list.
 * -------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* Maximum supported rank (number of axes). */
#define NDARRAY_MAX_RANK 64

/* ---------------------------------------------------------------------------
 * dtype helpers. `NDType` itself lives in expr.h (NDArrayData needs it, and
 * expr.h cannot include ndarray.h). These operate on it. Complex dtypes store
 * elements as interleaved (re, im) pairs; the ndt_get/ndt_set pair is the
 * single choke point that widens/narrows between the stored representation and
 * a machine-precision (double re, double im) pair, so every generic consumer
 * (packing, arithmetic, sort, precision, predicates) stays dtype-agnostic.
 * -------------------------------------------------------------------------- */

/* Bytes per element: 8 (float64), 4 (float32), 16 (complex64), 8 (complex32). */
size_t ndt_elem_size(NDType dt);
/* 1 for real dtypes, 2 for complex. */
int ndt_components(NDType dt);
/* True for NDT_COMPLEX64 / NDT_COMPLEX32. */
bool ndt_is_complex(NDType dt);
/* Bytes per scalar component: 8 (float64/complex64) or 4 (float32/complex32). */
size_t ndt_comp_size(NDType dt);
/* Read element k of `buf` (dtype dt) as a (re, im) pair; *im is 0 for reals. */
void ndt_get(const void* buf, size_t k, NDType dt, double* re, double* im);
/* Write (re, im) into element k of `buf` (dtype dt), narrowing to the stored
 * component type; im is ignored for real dtypes. */
void ndt_set(void* buf, size_t k, NDType dt, double re, double im);
/* Map an option string ("float64"/"float32"/"complex64"/"complex32") to a
 * dtype. Returns true and sets *out on success; false on an unknown name. */
bool ndt_from_string(const char* s, NDType* out);
/* Canonical option-string name for a dtype (never NULL). */
const char* ndt_to_string(NDType dt);
/* Result dtype of a binary op on operands of dtypes a and b. Promotion lattice:
 * complex dominates real, 64-bit component dominates 32-bit. */
NDType ndt_promote(NDType a, NDType b);
/* Move a dtype onto the complex axis preserving component width (float32 ->
 * complex32, float64 -> complex64). */
NDType ndt_as_complex(NDType dt);
/* Build the Expr leaf for element k of NDArray `a`: expr_new_real for real
 * dtypes, Complex[re, im] for complex dtypes. Caller owns the result. */
Expr* ndarray_element_to_expr(const Expr* a, size_t k);

/* True when `e` is an EXPR_NDARRAY node. */
bool is_ndarray(const Expr* e);

/* Total element count (product of dims). */
size_t ndarray_size(const Expr* a);

/* Probe a nested List for machine-precision-packable shape: rejects jagged
 * lists and any leaf that isn't EXPR_INTEGER/EXPR_REAL. On success returns
 * a newly built NDArray[...] (caller owns it); on failure (jagged, empty,
 * non-machine-precision leaf, or not a List at all) returns NULL and the
 * caller keeps the original list unevaluated. `dtype` selects the packed
 * element type; for real dtypes leaves must be Integer/Real, for complex dtypes
 * Integer/Real/Complex are all accepted. */
Expr* ndarray_from_nested_list(const Expr* list, NDType dtype);

/* Inverse of ndarray_from_nested_list: rebuilds the equivalent nested
 * List[...] tree from an NDArray's flat buffer. Caller owns the result. */
Expr* ndarray_to_nested_list(const Expr* a);

/* NDArray[nested_list] constructor builtin: packs `res`'s sole argument, or
 * returns NULL (leave unevaluated) if it can't be packed. */
Expr* builtin_ndarray(Expr* res);

/* NDArrayQ[expr] builtin: True iff expr is an NDArray value, else False. */
Expr* builtin_ndarrayq(Expr* res);

/* DataType[ndarray] builtin: the dtype option string ("float64" etc.) of an
 * NDArray; returns NULL (stays symbolic) for any non-NDArray argument. */
Expr* builtin_datatype(Expr* res);

/* Fast C-level Dot for two NDArray operands of rank 1 or 2: contracts the
 * trailing axis of `a` with the leading axis of `b` using raw double loops
 * (no symbolic Times/Plus). Returns a new EXPR_NDARRAY (rank 1 or 2), or a
 * bare EXPR_REAL for a vector.vector contraction (rank 0). Returns NULL if
 * either rank is > 2 (caller should fall back to the generic tensor path);
 * sets *shape_error and returns NULL if the contracted dimensions disagree. */
Expr* ndarray_dot2(const Expr* a, const Expr* b, bool* shape_error);

/* Fast C-level elementwise Plus (is_plus=true) or Times (is_plus=false) over
 * n NDArray operands: loops over the flat double buffers directly, no
 * symbolic Times/Plus per element. Returns a new EXPR_NDARRAY with the same
 * shape, or NULL if any operand isn't EXPR_NDARRAY or the shapes disagree
 * (caller falls back to the generic symbolic path, which treats mismatched
 * NDArray operands as opaque non-numeric terms). Broadcasting (numpy-style
 * shape compatibility) is intentionally not yet handled here. */
Expr* ndarray_elementwise(Expr** args, size_t n, bool is_plus);

/* Fast C-level elementwise Power over two same-shape NDArray operands (a^b,
 * element by element). Promotes the result dtype to complex when either operand
 * is complex or a real base with a non-integer exponent leaves the real axis
 * (e.g. (-1.0)^0.5). Returns NULL if either operand isn't an NDArray or the
 * shapes disagree. */
Expr* ndarray_elementwise_power(const Expr* a, const Expr* b);

/* Fast C-level NDArray raised to a scalar exponent (er + ei*I): the common
 * A^2 / A^0.5 case. Same complex-promotion rule as ndarray_elementwise_power.
 * Returns NULL if `a` isn't an NDArray. */
Expr* ndarray_scalar_power(const Expr* a, double er, double ei);

/* Fast C-level scalar base (br + bi*I) raised to an NDArray exponent, elementwise
 * (numpy `2 ** arr`). Promotes to complex on a complex base/exponent or a
 * negative real base with a non-integer exponent. NULL if `exp_arr` isn't an
 * NDArray. */
Expr* ndarray_base_scalar_power(double br, double bi, const Expr* exp_arr);

/* ---------------------------------------------------------------------------
 * Element-wise scalar-function acceleration.
 *
 * Every elementary/special function is Listable, but Listable threading only
 * sees List-headed args, not NDArrays. These descriptors let a function head
 * register a machine-precision scalar kernel that the evaluator maps directly
 * over an NDArray's flat buffer (no per-element Expr / evaluator round-trip),
 * exactly as Plus/Times/Power do. A kernel writes its result through the
 * (double re, double im) choke point, so one kernel serves all four dtypes.
 * A kernel returns false for an element it cannot compute at machine precision
 * (non-finite / pole); the caller then degrades to the faithful List path.
 * -------------------------------------------------------------------------- */

/* Unary kernel: f[array]. `cplx` is always valid. `real` is an optional faster
 * path used only when the input dtype is real AND `real_closed` is true (a real
 * input is guaranteed a real output); NULL means "use cplx and take the real
 * part". When `real_closed` is false (Sqrt/Log/ArcSin/... which can leave the
 * real axis on real input) the engine computes via `cplx` and promotes the
 * result dtype to complex iff any element came out non-real. */
typedef struct {
    bool (*cplx)(double ire, double iim, double* ore, double* oim);
    bool (*real)(double in, double* out);
    bool real_closed;
    /* Projection functions (Abs/Re/Im/Arg): the result is always real even for
     * complex input. The engine writes a real-dtype array, taking `ore` from
     * `cplx` as the value (`oim` ignored). Overrides real_closed. */
    bool to_real;
} NDUnaryKernel;

/* Binary scalar-index kernel: f[scalar, array] or f[array, scalar] — one
 * operand is an NDArray, the other a broadcast numeric scalar (BesselJ[n, arr],
 * Log[b, arr], ArcTan[arr, y], ...). `real_closed` as above. */
typedef struct {
    bool (*cplx)(double are, double aim, double bre, double bim,
                 double* ore, double* oim);
    bool real_closed;
} NDBinaryKernel;

/* Map unary kernel `k` over NDArray `a`, returning a new EXPR_NDARRAY of the
 * same shape (dtype narrowed/promoted per the kernel's real-closedness). Returns
 * NULL if `a` isn't an NDArray or any element's kernel reports failure (caller
 * degrades to ndarray_delist_and_reeval). */
Expr* ndarray_map_unary(const Expr* a, const NDUnaryKernel* k);

/* Map binary kernel `k` over exactly one NDArray operand (`a0` or `a1`) and one
 * broadcast numeric scalar. Returns NULL if the operands aren't that shape or
 * any element's kernel reports failure. */
Expr* ndarray_map_binary(const Expr* a0, const Expr* a1, const NDBinaryKernel* k);

/* Faithful degrade path: rebuild `call` with every NDArray arg replaced by its
 * equivalent nested List and re-evaluate, so the result matches f[{...}] exactly
 * (Listable threads element-wise, poles yield ComplexInfinity, etc.). Borrows
 * `call` (never frees it); caller owns the returned Expr. */
Expr* ndarray_delist_and_reeval(const Expr* call);

/* When every operand in args[0..n) is an NDArray but they are not all the same
 * shape, print a one-line `NDArray::shape` warning naming the two offending
 * shapes (verb is the operation, e.g. "added"/"multiplied") and return true.
 * Returns false when the operands are not all NDArrays or the shapes are
 * uniform (nothing to warn about). Respects the arithmetic-warning mute so it
 * stays quiet inside internal computations. */
bool ndarray_warn_shape_mismatch(Expr** args, size_t n, const char* verb);

/* When at least one operand in args[0..n) is an NDArray and at least one other
 * operand is symbolic (not itself an NDArray and not numeric per
 * expr_is_numeric_like — e.g. a bare symbol, or any non-numeric expression),
 * print a one-line `NDArray::sym` warning (verb is the operation, e.g.
 * "added"/"multiplied"/"raised to a power") and return true. NDArray objects
 * are purely numeric, so such a combination can never be carried out
 * elementwise; the caller leaves the expression unevaluated. Returns false
 * (nothing to warn about) when there is no NDArray operand, or every non-array
 * operand is a numeric scalar that broadcasts. Respects the arithmetic-warning
 * mute so it stays quiet inside internal computations. */
bool ndarray_warn_symbolic(Expr** args, size_t n, const char* verb);

/* Register NDArray's builtins, attributes, and docstring. */
void ndarray_init(void);

/* Register the elementary-function NDArray element-wise kernels (ndkernels.c)
 * on their function heads (Sin, Cos, Exp, Log, ArcTan, Abs, ...). Special
 * functions register their own kernels from their own modules. */
void ndkernels_init(void);

#endif /* MATHILDA_NDARRAY_H */
