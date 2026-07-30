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
#include <stdint.h>

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
/* int64 audit tripwire -- see the comment at the definition in ndarray.c.
 * Called from the NDT_INT64 arm of ndt_get/ndt_set, i.e. exactly where a
 * value that must stay exact is being routed through `double`. A stable
 * breakpoint target: `break nd_int64_lossy_hit`. MATHILDA_PACK_DIAG=1 warns,
 * =abort aborts. `nd_int64_lossy_count` is the running total. */
void nd_int64_lossy_hit(const char* what);
extern unsigned long nd_int64_lossy_count;

/* EXACT int64 element access. The pair above routes through `double` and is
 * therefore exact only to 2^53, so every read or write of an NDT_INT64 buffer
 * that must be faithful goes through these two instead. */
int64_t ndt_get_i(const void* buf, size_t k, NDType dt);
void    ndt_set_i(void* buf, size_t k, NDType dt, int64_t v);
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
/* Build the Expr leaf for element k of a raw `dt` buffer: an exact Integer for
 * NDT_INT64 (read through ndt_get_i, since ndt_get is exact only to 2^53),
 * expr_new_real for real dtypes, Complex[re, im] for complex dtypes. Caller owns
 * the result.
 *
 * This is the ONE place that decides what head an element materialises as. Every
 * unpack path routes through it -- ndarray_element_to_expr, rebuild_level,
 * ndarray_part's two scalar leaves, ndarray_dot2's rank-0 result -- so none of
 * them can drift from the others (an Integer silently becoming a Real is the
 * exact bug class this centralisation exists to prevent). */
Expr* ndarray_buffer_element_to_expr(const void* buf, size_t k, NDType dt);

/* Build the Expr leaf for element k of NDArray `a`. Thin wrapper over
 * ndarray_buffer_element_to_expr with `a`'s buffer and dtype. */
Expr* ndarray_element_to_expr(const Expr* a, size_t k);

/* True when `e` is an EXPR_NDARRAY node -- EITHER surface, a visible
 * `NDArray[...]` or a packed List. Every existing fast path wants this, because
 * a fast path cares about the buffer, not how the value prints. */
bool is_ndarray(const Expr* e);

/* True when `e` is an EXPR_NDARRAY presenting as a packed List (present_as ==
 * NDA_HEAD_LIST): a value the user sees as an ordinary List. This is the
 * predicate for the handful of places that must tell the two surfaces apart --
 * Head, AtomQ, ListQ, printing, and the Part-assignment exactness rule. */
bool is_packed_list(const Expr* e);

/* Which of two array operands decides a binary result's presentation.
 *
 * A VISIBLE NDArray[...] dominates a packed List: if the user wrote NDArray on
 * either side of a Dot or a Plus they asked for an NDArray, and they get one
 * back. Two packed Lists give a packed List; two visible arrays give a visible
 * array. Either argument may be NULL or a non-array scalar (a broadcast
 * operand contributes nothing). Returns NULL when neither is an array, which
 * expr_new_ndarray_like treats as "raw". */
const Expr* nd_present_src2(const Expr* a, const Expr* b);

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

/* ndarray_from_nested_list, but inheriting `src`'s presentation (visible
 * NDArray[...] vs packed List). Use this for a REPACK: whenever the List being
 * packed came from materialising an existing array, because a builtin took the
 * materialise-reuse-repack detour. Packing it raw instead turns a packed list
 * into a visible NDArray at the first such builtin -- which is a user-visible
 * wrong answer with no crash to find it by. `src` may be NULL (treated as raw).
 *
 * Same relationship as expr_new_ndarray_like to expr_new_ndarray_raw. */
Expr* ndarray_from_nested_list_like(const Expr* src, const Expr* list, NDType dtype);

/* Inverse of ndarray_from_nested_list: rebuilds the equivalent nested
 * List[...] tree from an NDArray's flat buffer. Caller owns the result. */
Expr* ndarray_to_nested_list(const Expr* a);

/* Native Part[a, i1, ..., iN] over a leading run of plain-integer subscripts,
 * operating directly on the flat buffer (no full materialization). N == rank
 * yields a scalar leaf (Real or Complex[re,im]); N < rank yields a contiguous
 * rank-(rank-N) sub-NDArray. Subscripts are 1-based, negatives count from the
 * end. A non-integer subscript (Span/All/List/...) sets *degrade true and
 * returns NULL, so the caller can fall back to ndarray_delist_and_reeval; an
 * out-of-range subscript or too many subscripts returns NULL with *degrade
 * false (leave the Part unevaluated). Caller owns any returned Expr. */
Expr* ndarray_part(const Expr* a, Expr** indices, size_t nindices, bool* degrade);

/* Part[a, i1, ..., iN] = rhs, written straight into `a`'s buffer. The positions
 * come from the same per-axis selector ndarray_part gathers through, so the
 * full spec vocabulary (Integer / negative / All / Span / List of positions /
 * implicit trailing axes) names the same elements for a write as for a read.
 * `rhs` is one number broadcast to every selected position, or an NDArray with
 * a matching element count. Returns false — having written NOTHING — for an
 * unsupported spec, an out-of-range subscript, a count mismatch, or a complex
 * value into a real buffer.
 *
 * IN PLACE: only sound for an array the caller owns outright. Compiled bytecode
 * is the intended caller and checks that at emit time; the interpreter's Part
 * assignment is value semantics and rebuilds the array instead. */
bool ndarray_part_set(Expr* a, Expr* const* indices, size_t nindices,
                      const Expr* rhs);

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

    /* ---------------------------------------------------------------- *
     *  NARROWING: real in, EXACT INTEGER out.
     * ---------------------------------------------------------------- *
     * Floor, Ceiling, Round, IntegerPart, Sign and UnitStep all answer with an
     * exact Integer -- Floor[{0.25, 0.5, 1.0}] is {0, 0, 1}, not {0., 0., 1.}.
     * The two categories above cannot express that: real_closed keeps the input
     * dtype and would write 1.0 where the List gives the exact 1, and to_real
     * likewise. That is why all six sat on pack.c's NOT_AWARE list, materialising
     * rather than computing -- and why UnitStep, which has no other sensible
     * category at all, had no kernel whatsoever and cost ~500 ns/element.
     *
     * When `to_int` is set the engine writes an NDT_INT64 array:
     *   to_int_r   from a real element; false = not representable as an exact
     *              int64 (non-finite, or beyond the range), which ABANDONS the
     *              whole array so the List path answers with a bignum.
     *   to_int_i   from an exact int64 element; NULL declines, leaving int64
     *              input to materialise as before.
     *
     * These are additive: `cplx`, `real` and `real_closed` keep whatever they
     * already said, so the Compile VM and every other kernel consumer are
     * unaffected -- only ndarray_map_unary prefers this path. */
    bool (*to_int_r)(double x, int64_t* out);
    bool (*to_int_i)(int64_t x, int64_t* out);
    bool to_int;
} NDUnaryKernel;

/* Binary scalar-index kernel: f[scalar, array] or f[array, scalar] — one
 * operand is an NDArray, the other a broadcast numeric scalar (BesselJ[n, arr],
 * Log[b, arr], ArcTan[arr, y], ...). `real_closed` as above. */
typedef struct {
    bool (*cplx)(double are, double aim, double bre, double bim,
                 double* ore, double* oim);
    bool real_closed;
} NDBinaryKernel;

/* N-ary scalar kernel: f[a1, ..., an] for n >= 3 (LerchPhi, Hypergeometric1F1,
 * Hypergeometric2F1, ...).  The unary/binary descriptors above cannot express
 * these, and the Compile[] VM is the consumer that needs them — the NDArray
 * element-wise path still degrades for such heads, since "one array plus
 * broadcast scalars" is a different question from "n scalars".  Arguments arrive
 * as parallel (re, im) arrays of length n, the same choke point the other
 * kernels use. */
typedef struct {
    bool (*cplx)(const double* re, const double* im, size_t n,
                 double* ore, double* oim);
    size_t nargs;            /* exact arity this kernel accepts */
    bool real_closed;
} NDNaryKernel;

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

/* Degrade path for an exact-integer buffer whose answer is not machine
 * representable (Rational, radical, exact Complex). Returns the re-evaluated
 * call on a List, or NULL when no operand is an int64 buffer. See the
 * definition. */
Expr* ndarray_int64_delist_retry(const Expr* call);

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
