#include "sort.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "assoc.h"
#include "ndarray.h"   /* ndt_get for NDArray canonical ordering */
#include "ndstruct.h"  /* ndstruct_sort / ndstruct_ordering NDArray fast paths */
#include "pack.h"      /* pack_offer — a sorted machine list packs */
#include "common.h"    /* builtin_arg_error */
#include "take_drop.h" /* get_seq_spec_indices — Ordering[list, seq] reuses Take's spec */
#include "list/list_common.h" /* is_real_numeric / is_infinity / is_minus_infinity / is_listq */
#include "ndreduce.h"  /* ndred_ranked_min / ndred_ranked_max — NDArray order-statistic path */
#include "numeric.h"   /* numericalize / numeric_machine_spec — key for symbolic reals */
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* ====================================================================
 * Canonical comparison of expressions (Mathematica-compatible `Order`).
 *
 * `expr_compare` lives here -- rather than in expr.c -- so the canonical
 * order policy is co-located with the user-facing Sort/OrderedQ builtins
 * that surface it.  The Orderless-attribute argument-sorting pipeline in
 * eval.c calls into this same comparator via expr.h.
 *
 * Algorithm, applied uniformly to every pair:
 *   1. Numeric atoms (Integer, Real, Rational, BigInt, MPFR) come first
 *      and compare by value.  Numeric < non-numeric.
 *   2. Strings come next, compared lexicographically.
 *   3. Everything else (Symbols and Functions) is compared by polynomial
 *      degree vector:
 *        - collect every Symbol name appearing in either expression,
 *        - sort those names in REVERSE alphabetical order so the
 *          lex-last variable is the most significant dimension,
 *        - for each name v, compute the polynomial degree of each side
 *          in v as a double (see expr_poly_degree),
 *        - lex-compare the degree vectors ascending.
 *      This makes `Plus[a x^2, b x]` sort as `b x + a x^2` and gives
 *      the multivariate grevlex-with-reverse-alpha display order on
 *      `Expand[(1 + x y^2 + y x^3)^3]`.
 *   4. Equal degree vectors fall back to a structural tiebreak: bare
 *      Symbol before compound, then head/arity/args recursively (each
 *      arg re-enters this same comparator).
 * ==================================================================== */

/* True for numeric atom kinds with an extractable double value. */
static bool is_atomic_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    /* is_rational_like, not is_rational: a Rational with a bignum numerator or
     * denominator (e.g. Rational[1, 10^25]) must still sort as the number it
     * is.  With the int64-only predicate it fell out of the numeric tier and
     * was misordered by the polynomial-degree machinery below. */
    if (is_rational_like(e)) return true;
    return false;
}

/* Extract the numerator and denominator of a rational-like atom (Integer,
 * BigInt, or Rational[int-like, int-like]) as GMP integers, initialising both
 * (caller mpz_clears).  Integer/BigInt yield den = 1.  The denominator keeps
 * its stored sign; canonical rationals are positive. */
static void rat_to_mpz_parts(const Expr* e, mpz_t num, mpz_t den) {
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, num);
        mpz_init_set_ui(den, 1);
    } else { /* Rational[num, den] */
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
    }
}

/* True for a literal complex number Complex[re, im] whose real and imaginary
 * parts are both real numeric atoms (Integer / Rational / Real / BigInt /
 * MPFR). Such a node is a NUMBER, not a polynomial expression -- it must not
 * be routed through the polynomial-degree machinery, whose symbol collector
 * would otherwise mistake the wrapper heads `Complex` (and `Rational`, when a
 * component is rational) for polynomial variables. That misroute made
 * `I/2 * Pi` canonicalize as `Pi * I/2` while `I * Pi` (integer imaginary
 * part, no `Rational` head) canonicalized correctly -- an inconsistency. */
static bool is_complex_number(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Complex) return false;
    return is_atomic_numeric(e->data.function.args[0]) &&
           is_atomic_numeric(e->data.function.args[1]);
}

static double get_numeric_value(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_BIGINT) return mpz_get_d(e->data.bigint);
    if (e->type == EXPR_REAL) return e->data.real;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) return (double)n / d;
    /* Bignum-component Rational: is_rational rejected it and this used to
     * return 0, collapsing every such value to a false tie. */
    if (is_rational_like(e) && e->type == EXPR_FUNCTION) {
        mpz_t num, den;
        rat_to_mpz_parts(e, num, den);
        double v = mpz_get_d(num) / mpz_get_d(den);
        mpz_clear(num); mpz_clear(den);
        return v;
    }
    return 0;
}

/* Case-insensitive then case-sensitive lex (lowercase wins on case ties);
 * preserves the behaviour the previous expr.c implementation pinned for
 * string-vs-string ordering. */
static int string_compare_canonical(const char* s1, const char* s2) {
    const char* p1 = s1;
    const char* p2 = s2;
    while (*p1 && *p2) {
        char l1 = tolower((unsigned char)*p1);
        char l2 = tolower((unsigned char)*p2);
        if (l1 != l2) return (l1 < l2) ? -1 : 1;
        p1++; p2++;
    }
    if (*p1) return 1;
    if (*p2) return -1;

    p1 = s1; p2 = s2;
    while (*p1 && *p2) {
        if (*p1 != *p2) {
            if (islower((unsigned char)*p1) && !islower((unsigned char)*p2)) return -1;
            if (!islower((unsigned char)*p1) && islower((unsigned char)*p2)) return 1;
            return (*p1 < *p2) ? -1 : 1;
        }
        p1++; p2++;
    }
    return 0;
}

/* Return true iff `e` contains a Symbol with the interned name `sym`
 * anywhere in its tree (heads and args included). */
static bool expr_contains_symbol(const Expr* e, const char* sym) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == sym;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_contains_symbol(e->data.function.head, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_contains_symbol(e->data.function.args[i], sym)) return true;
    }
    return false;
}

/* Polynomial degree of `e` in the variable `v`, as a double:
 *   0.0    -- e does not contain v.
 *   k > 0  -- e is v (k=1), Power[v, k] for positive integer or rational
 *             k, or a Times product whose direct factors sum to k.
 *   +inf   -- e contains v in a non-polynomial position: inside a
 *             compound that isn't bare-v / Power[v, positive rational]
 *             (e.g. Sin[v], Plus[..., v, ...], Power[v, negative],
 *             Power[non-v-containing-v, _]).
 *
 * Only DIRECT factors of a Times contribute to the sum; a non-poly
 * direct factor that itself contains v (e.g. Sin[v]) saturates the
 * result at +inf.  This matches Mathematica:
 *   - `Sin[x] + x`      displays as `x + Sin[x]`     (poly-in-x before non-poly).
 *   - `Sqrt[x] + x^2`   displays as `Sqrt[x] + x^2`  (deg 1/2 < 2).
 *   - `1/x + 1`         displays as `1 + 1/x`        (negative-exp x -> +inf).
 *   - `(3 + 6 a + 3 a^2) x` sorts by deg_x = 1 regardless of the Plus
 *      coefficient's internal a-degree. */
static double expr_poly_degree(const Expr* e, const char* v) {
    if (!e) return 0.0;
    if (e->type == EXPR_SYMBOL) return (e->data.symbol.name == v) ? 1.0 : 0.0;
    if (e->type != EXPR_FUNCTION) return 0.0;
    if (e->data.function.head->type != EXPR_SYMBOL) {
        return expr_contains_symbol(e, v) ? INFINITY : 0.0;
    }
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == v) {
            if (exp->type == EXPR_INTEGER) {
                int64_t k = exp->data.integer;
                if (k < 0) return INFINITY;
                return (double)k;
            }
            int64_t num, den;
            if (is_rational((Expr*)exp, &num, &den) && den > 0) {
                if (num < 0) return INFINITY;
                return (double)num / (double)den;
            }
            return INFINITY;
        }
        return (expr_contains_symbol(base, v) || expr_contains_symbol(exp, v))
               ? INFINITY : 0.0;
    }
    if (h == SYM_Times) {
        double sum = 0.0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double d = expr_poly_degree(e->data.function.args[i], v);
            if (isinf(d)) return INFINITY;
            sum += d;
        }
        return sum;
    }
    return expr_contains_symbol(e, v) ? INFINITY : 0.0;
}

/* Growable buffer of distinct interned-symbol-name pointers.  We need
 * collection to be ORDER-INDEPENDENT across a vs b: a silent truncation
 * after `a` would leave `b`'s unique symbols unrepresented, making the
 * resulting comparator asymmetric (compare(a,b) != -compare(b,a)) and
 * causing qsort on Orderless heads with hundreds of unknowns to
 * oscillate forever between two non-converging orderings. */
typedef struct {
    const char* stack_buf[64];  /* common-case avoid-malloc path */
    const char** items;          /* points at stack_buf or heap     */
    size_t count;
    size_t cap;
    bool ok;
} SymSet;

static void symset_init(SymSet* s) {
    s->items = s->stack_buf;
    s->count = 0;
    s->cap = sizeof(s->stack_buf) / sizeof(s->stack_buf[0]);
    s->ok = true;
}

static void symset_free(SymSet* s) {
    if (s->items && s->items != s->stack_buf) free((void*)s->items);
    s->items = NULL;
}

static bool symset_grow(SymSet* s) {
    size_t new_cap = s->cap * 2;
    const char** new_buf = (const char**)malloc(sizeof(const char*) * new_cap);
    if (!new_buf) { s->ok = false; return false; }
    memcpy(new_buf, s->items, sizeof(const char*) * s->count);
    if (s->items != s->stack_buf) free((void*)s->items);
    s->items = new_buf;
    s->cap = new_cap;
    return true;
}

static void symset_add(SymSet* s, const char* sym) {
    if (!s->ok) return;
    for (size_t i = 0; i < s->count; i++) {
        if (s->items[i] == sym) return;
    }
    if (s->count >= s->cap && !symset_grow(s)) return;
    s->items[s->count++] = sym;
}

/* Append every distinct Symbol-name pointer appearing in `e` into `s`,
 * DFS-visiting heads and args.  Grows the buffer as needed; never
 * truncates. */
static void collect_symbols_in(const Expr* e, SymSet* s) {
    if (!e || !s->ok) return;
    if (e->type == EXPR_SYMBOL) {
        symset_add(s, e->data.symbol.name);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    collect_symbols_in(e->data.function.head, s);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_symbols_in(e->data.function.args[i], s);
    }
}

/* qsort comparator that sorts char* pointers in REVERSE lex order, so
 * the lex-LAST variable name (e.g. `y` in {a, b, x, y}) ends up first
 * in the sorted array and drives the most significant dimension of the
 * degree-vector compare. */
static int symbol_reverse_cmp(const void* pa, const void* pb) {
    const char* a = *(const char* const*)pa;
    const char* b = *(const char* const*)pb;
    return strcmp(b, a);
}

/* Compare two packed lists straight off their buffers. Valid ONLY when the
 * shapes match exactly and the elements materialise to the same head class,
 * because the List order is elementwise while a buffer scan would otherwise
 * order by dims first: {{2.,0.},{1.,9.}} vs {{1.,9.,9.},{2.,0.,0.}} is +1
 * elementwise (2. > 1.) but -1 by dims (2 < 3). Sets *ok false when it cannot
 * decide, and the caller materialises instead. */
static int nd_packed_compare_fast(const Expr* a, const Expr* b, bool* ok) {
    const NDArrayData* x = &a->data.ndarray;
    const NDArrayData* y = &b->data.ndarray;
    *ok = false;
    if (x->rank != y->rank) return 0;
    size_t n = 1;
    for (int i = 0; i < x->rank; i++) {
        if (x->dims[i] != y->dims[i]) return 0;
        n *= (size_t)x->dims[i];
    }
    bool xi = (x->dtype == NDT_INT64), yi = (y->dtype == NDT_INT64);
    if (xi != yi) return 0;                                    /* Integer vs Real */
    if (!xi && ndt_is_complex(x->dtype) != ndt_is_complex(y->dtype))
        return 0;                                              /* Complex vs Real */
    *ok = true;
    if (xi) {
        for (size_t i = 0; i < n; i++) {
            int64_t va = ndt_get_i(x->data, i, x->dtype);
            int64_t vb = ndt_get_i(y->data, i, y->dtype);
            if (va != vb) return (va < vb) ? -1 : 1;
        }
        return 0;
    }
    for (size_t i = 0; i < n; i++) {
        double are, aim, bre, bim;
        ndt_get(x->data, i, x->dtype, &are, &aim);
        ndt_get(y->data, i, y->dtype, &bre, &bim);
        if (are < bre) return -1;
        if (are > bre) return 1;
        if (aim < bim) return -1;
        if (aim > bim) return 1;
    }
    return 0;
}

int expr_compare(const Expr* a, const Expr* b) {
    if (a == b) return 0;
    if (!a) return -1;
    if (!b) return 1;

    /* 0. Packed lists order as the Lists they are, NOT in the NDArray tier at
     * step 5 below (which sorts after everything and treats dtype as identity).
     * Two packed lists of matching shape and element class compare straight off
     * their buffers -- the case that costs, since Orderless sorting of an
     * arithmetic head hits it. Everything else materialises: correct by
     * construction, and only paid when a packed list is genuinely being ordered
     * against structural material. */
    if (is_packed_list(a) || is_packed_list(b)) {
        if (is_packed_list(a) && is_packed_list(b)) {
            bool ok;
            int c = nd_packed_compare_fast(a, b, &ok);
            if (ok) return c;
        }
        /* ONE packed list against a NUMBER, a literal Complex or a String.
         *
         * Steps 1, 1b and 2 below settle those three tiers on membership
         * alone: whichever side is the number/complex/string wins, and the
         * other side's CONTENTS are never read. Materialising to discover that
         * is pure waste, and it is not a rare shape -- it is what `Orderless`
         * does on every call. The evaluator sorts `GCD[1234, cv]` before
         * dispatching, that comparison is packed-against-scalar, and 200000
         * Expr nodes were being built and thrown away per sort. Measured at
         * 10^6/200000 elements: GCD 38.6 ms packed against 16.7 ms visible and
         * LCM 29.4 against 8.3, purely from this -- `Mod` and `Quotient`, the
         * same kernels without ATTR_ORDERLESS, showed no gap at all.
         *
         * The stand-in is an EMPTY List rather than a hardcoded verdict. A
         * packed list IS a List, so `List[]` lands in exactly the same tier and
         * the existing code computes the answer; writing `return 1` here would
         * restate an ordering rule in a second place and let the two drift.
         * Restricted to the three tiers where the result provably cannot depend
         * on the elements -- a Symbol or a general expression on the other side
         * reaches step 3's polynomial-degree walk, which does read them, and
         * still materialises. */
        const Expr* other = is_packed_list(a) ? b : a;
        if (is_atomic_numeric(other) || is_complex_number(other) ||
            other->type == EXPR_STRING) {
            Expr* stub = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
            int c = is_packed_list(a) ? expr_compare(stub, b)
                                      : expr_compare(a, stub);
            expr_free(stub);
            return c;
        }

        Expr* ma = is_packed_list(a) ? ndarray_to_nested_list(a) : NULL;
        Expr* mb = is_packed_list(b) ? ndarray_to_nested_list(b) : NULL;
        int c = expr_compare(ma ? ma : a, mb ? mb : b);
        if (ma) expr_free(ma);
        if (mb) expr_free(mb);
        return c;
    }

    /* 1. Numeric atoms. */
    bool a_atomic = is_atomic_numeric(a);
    bool b_atomic = is_atomic_numeric(b);
    if (a_atomic && b_atomic) {
        if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
            if (a->data.integer < b->data.integer) return -1;
            if (a->data.integer > b->data.integer) return 1;
            return 0;
        }
        if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
            mpz_t ma2, mb2;
            expr_to_mpz(a, ma2);
            expr_to_mpz(b, mb2);
            int cmp = mpz_cmp(ma2, mb2);
            mpz_clear(ma2);
            mpz_clear(mb2);
            return cmp;
        }
        /* Exact rationals (bignum components included).  The double fallback
         * below reads a bignum-component Rational as 0 (get_numeric_value
         * could not decode it), collapsing distinct values to a tie and
         * misordering Sort/Min/Max/Median.  Cross-multiply exactly; both
         * denominators are positive (canonical form), so the inequality
         * direction is preserved. */
        if (is_rational_like(a) && is_rational_like(b)) {
            mpz_t an, ad, bn, bd, lhs, rhs;
            rat_to_mpz_parts(a, an, ad);
            rat_to_mpz_parts(b, bn, bd);
            mpz_init(lhs); mpz_init(rhs);
            mpz_mul(lhs, an, bd);
            mpz_mul(rhs, bn, ad);
            int c = mpz_cmp(lhs, rhs);
            mpz_clear(an); mpz_clear(ad); mpz_clear(bn); mpz_clear(bd);
            mpz_clear(lhs); mpz_clear(rhs);
            if (c != 0) return c < 0 ? -1 : 1;
            if (a->type != b->type) return (int)a->type - (int)b->type;
            return 0;
        }
        double va = get_numeric_value(a);
        double vb = get_numeric_value(b);
        if (va < vb) return -1;
        if (va > vb) return 1;
        if (a->type != b->type) return (int)a->type - (int)b->type;
        return 0;
    }
    if (a_atomic) return -1;
    if (b_atomic) return 1;

    /* 1b. Complex numbers. A literal complex sits in its own tier just after
     * the real numeric atoms and before all symbolic material (symbols, Pi,
     * Log, Power, ...), mirroring Mathematica's `real < non-real complex <
     * symbol` canonical order. Ordering among complexes is by real part then
     * imaginary part (each a real numeric atom, compared exactly by the block
     * above). This keeps `Complex[re, im]` out of the polynomial-degree
     * machinery entirely, so its ordering no longer depends on whether a
     * component happens to be an Integer or a Rational. */
    bool a_cplx = is_complex_number(a);
    bool b_cplx = is_complex_number(b);
    if (a_cplx && b_cplx) {
        int c = expr_compare(a->data.function.args[0], b->data.function.args[0]);
        if (c != 0) return c;
        return expr_compare(a->data.function.args[1], b->data.function.args[1]);
    }
    if (a_cplx) return -1;
    if (b_cplx) return 1;

    /* 2. Strings. */
    if (a->type == EXPR_STRING && b->type == EXPR_STRING) {
        return string_compare_canonical(a->data.string, b->data.string);
    }
    if (a->type == EXPR_STRING) return -1;
    if (b->type == EXPR_STRING) return 1;

    /* 3. Polynomial degree vector comparison. */
    SymSet vs;
    symset_init(&vs);
    collect_symbols_in(a, &vs);
    collect_symbols_in(b, &vs);

    qsort((void*)vs.items, vs.count, sizeof(const char*), symbol_reverse_cmp);

    int deg_cmp = 0;
    for (size_t i = 0; i < vs.count; i++) {
        double da = expr_poly_degree(a, vs.items[i]);
        double db = expr_poly_degree(b, vs.items[i]);
        if (da < db) { deg_cmp = -1; break; }
        if (da > db) { deg_cmp = 1; break; }
    }
    symset_free(&vs);
    if (deg_cmp != 0) return deg_cmp;

    /* 4. Structural tiebreak. */
    if (a->type == EXPR_SYMBOL && b->type == EXPR_SYMBOL) {
        if (a->data.symbol.name == b->data.symbol.name) return 0;
        return strcmp(a->data.symbol.name, b->data.symbol.name);
    }
    if (a->type == EXPR_SYMBOL) return -1;
    if (b->type == EXPR_SYMBOL) return 1;

    if (a->type == EXPR_FUNCTION && b->type == EXPR_FUNCTION) {
        int head_cmp = expr_compare(a->data.function.head, b->data.function.head);
        if (head_cmp != 0) return head_cmp;
        size_t ac = a->data.function.arg_count;
        size_t bc = b->data.function.arg_count;
        if (ac != bc) return (ac < bc) ? -1 : 1;
        for (size_t i = 0; i < ac; i++) {
            int c = expr_compare(a->data.function.args[i], b->data.function.args[i]);
            if (c != 0) return c;
        }
        return 0;
    }

    /* 5. NDArrays: dtype, then rank, then dims, then row-major data (re then
     * im), lexicographic — a stable total order across dtypes. */
    if (a->type == EXPR_NDARRAY && b->type == EXPR_NDARRAY) {
        NDType dta = a->data.ndarray.dtype, dtb = b->data.ndarray.dtype;
        if (dta != dtb) return (dta < dtb) ? -1 : 1;
        if (a->data.ndarray.rank != b->data.ndarray.rank) {
            return (a->data.ndarray.rank < b->data.ndarray.rank) ? -1 : 1;
        }
        size_t n = 1;
        for (int i = 0; i < a->data.ndarray.rank; i++) {
            int64_t da = a->data.ndarray.dims[i], db = b->data.ndarray.dims[i];
            if (da != db) return (da < db) ? -1 : 1;
            n *= (size_t)da;
        }
        for (size_t i = 0; i < n; i++) {
            double are, aim, bre, bim;
            ndt_get(a->data.ndarray.data, i, dta, &are, &aim);
            ndt_get(b->data.ndarray.data, i, dtb, &bre, &bim);
            if (are < bre) return -1;
            if (are > bre) return 1;
            if (aim < bim) return -1;
            if (aim > bim) return 1;
        }
        return 0;
    }
    if (a->type == EXPR_NDARRAY) return 1;
    if (b->type == EXPR_NDARRAY) return -1;

    /* 6. CompiledFunction objects: order by payload identity (stable within a
     * run; distinct objects compare unequal). */
    if (a->type == EXPR_COMPILED && b->type == EXPR_COMPILED) {
        uint64_t ia = compiled_function_identity(a->data.compiled);
        uint64_t ib = compiled_function_identity(b->data.compiled);
        if (ia != ib) return (ia < ib) ? -1 : 1;
        return 0;
    }
    if (a->type == EXPR_COMPILED) return 1;
    if (b->type == EXPR_COMPILED) return -1;

    return 0;
}

/* ==================================================================== */

static Expr* current_sort_p = NULL;

/* Compare ea against eb the way the sort family does: canonical order when p is
 * NULL, else the user ordering function p invoked as p[ea, eb]. Returns a
 * qsort-style sign — negative when ea sorts before eb, positive after, 0 equal.
 * "In order" for a custom p means it returns True or 1 (the Order/OrderedQ
 * surface convention); anything unrecognised defaults to "in order". Shared by
 * Sort's custom_compare and Ordering's index comparator so the two never drift. */
static int cmp_with_p(Expr* ea, Expr* eb, Expr* p) {
    if (p == NULL) {
        return expr_compare(ea, eb);
    }

    /* Call p[ea, eb]. */
    Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
    Expr* call = expr_new_function(expr_copy(p), args, 2);
    Expr* result = evaluate(call);
    expr_free(call);

    int cmp = -1; /* default to "in order" (ea before eb) */
    if (result->type == EXPR_INTEGER) {
        /* 1: ea before eb; 0: same; -1: ea after eb. */
        if (result->data.integer == 1) cmp = -1;
        else if (result->data.integer == 0) cmp = 0;
        else if (result->data.integer == -1) cmp = 1;
    } else if (result->type == EXPR_SYMBOL) {
        if (result->data.symbol.name == SYM_True) cmp = -1;
        else if (result->data.symbol.name == SYM_False) cmp = 1;
    }

    expr_free(result);
    return cmp;
}

static int custom_compare(const void* a, const void* b) {
    return cmp_with_p(*(Expr**)a, *(Expr**)b, current_sort_p);
}

Expr* builtin_sort(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];

    /* Sort[ndarray]: sort the flat buffer directly (see ndstruct.c). */
    if (is_ndarray(list)) return ndstruct_sort(res);

    if (list->type != EXPR_FUNCTION) {
        return expr_copy(list);
    }

    /* Sort[assoc] orders the entries by their values (keys follow). */
    if (res->data.function.arg_count == 1 && is_association(list))
        return assoc_sort_by_value(list);

    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    
    size_t count = list->data.function.arg_count;
    if (count == 0) {
        return expr_copy(list);
    }
    
    Expr** sorted_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        sorted_args[i] = expr_copy(list->data.function.args[i]);
    }
    
    // Set global context for comparison
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    qsort(sorted_args, count, sizeof(Expr*), custom_compare);
    
    current_sort_p = old_p;
    
    Expr* result = expr_new_function(expr_copy(list->data.function.head), sorted_args, count);
    free(sorted_args);

    /* A sorted machine-number list is worth packing: it is the shape most
     * likely to be reduced or indexed next. Declines when the head is not List
     * (Sort[g[3., 1.]] stays g[...]) or the elements are not uniform. */
    return pack_offer(result);
}

/* ------------------- SortBy ------------------- */

/* A payload paired with the sort key f[payload], for a key-once qsort
 * (avoids re-evaluating f during every comparison). */
typedef struct { Expr* payload; Expr* key; } SortByPair;

static int sortby_pair_cmp(const void* a, const void* b) {
    return expr_compare(((const SortByPair*)a)->key, ((const SortByPair*)b)->key);
}

/* SortBy[list, f]  — sort list elements by canonical order of f[element].
 * SortBy[assoc, f] — sort association entries by f[value] (keys follow).
 * SortBy[f]        — operator form: SortBy[f][expr] == SortBy[expr, f].
 * The key f[...] is evaluated once per element. */
Expr* builtin_sort_by(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        /* Operator form: SortBy[f] -> Function[SortBy[#, f]]. */
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner = expr_new_function(expr_new_symbol(SYM_SortBy), inner_args, 2);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc != 2) return NULL;

    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (coll->type != EXPR_FUNCTION) return expr_copy(coll);

    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;
    if (n == 0) return expr_copy(coll);

    SortByPair* pairs = malloc(sizeof(SortByPair) * n);
    /* SortBy[list, {f1, f2, ...}] sorts by f1, breaking ties with f2, ... .
     * The per-element key becomes the tuple {f1[e], ...}; expr_compare orders
     * equal-length lists lexicographically, giving exactly that behaviour. */
    bool multi = (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
                  f->data.function.head->data.symbol.name == SYM_List);

    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        /* For an association, sort by f applied to the value; else by f[elem]. */
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        pairs[i].payload = expr_copy(elem);
        if (multi) {
            size_t nc = f->data.function.arg_count;
            Expr** keyparts = malloc(sizeof(Expr*) * (nc ? nc : 1));
            for (size_t c = 0; c < nc; c++) {
                Expr* ca[1] = { expr_copy(subject) };
                Expr* call = expr_new_function(expr_copy(f->data.function.args[c]), ca, 1);
                keyparts[c] = evaluate(call);
                expr_free(call);
            }
            pairs[i].key = expr_new_function(expr_new_symbol(SYM_List), keyparts, nc);
            free(keyparts);
        } else {
            Expr* call_args[1] = { expr_copy(subject) };
            Expr* call = expr_new_function(expr_copy(f), call_args, 1);
            pairs[i].key = evaluate(call);
            expr_free(call);
        }
    }

    qsort(pairs, n, sizeof(SortByPair), sortby_pair_cmp);

    Expr** out = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) { out[i] = pairs[i].payload; expr_free(pairs[i].key); }
    free(pairs);

    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, n);
    free(out);
    return result;
}


/* ------------------- MaximalBy / MinimalBy ------------------- */

/* mode 0 = MaximalBy, mode 1 = MinimalBy.
 * MaximalBy[list, f]  -> the element(s) e maximising f[e] (all ties, in order).
 * MaximalBy[assoc, f] -> the entries whose value maximises f[value] (an assoc).
 * MaximalBy[f]        -> operator form. */
static Expr* maximal_minimal_by(Expr* res, int mode) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1) {
        const char* head = (mode == 0) ? SYM_MaximalBy : SYM_MinimalBy;
        Expr* slot_args[1] = { expr_new_integer(1) };
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
        Expr* inner_args[2] = { slot, expr_copy(res->data.function.args[0]) };
        Expr* inner = expr_new_function(expr_new_symbol(head), inner_args, 2);
        Expr* func_args[1] = { inner };
        return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
    }
    if (argc != 2) return NULL;

    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (coll->type != EXPR_FUNCTION) return NULL;

    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;
    if (n == 0) return expr_copy(coll);  /* empty in, empty out */

    Expr** keys = malloc(sizeof(Expr*) * n);   /* f-value per element (owned) */
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        Expr* call_args[1] = { expr_copy(subject) };
        Expr* call = expr_new_function(expr_copy(f), call_args, 1);
        keys[i] = evaluate(call);
        expr_free(call);
    }

    size_t best = 0;
    for (size_t i = 1; i < n; i++) {
        int c = expr_compare(keys[i], keys[best]);
        if ((mode == 0 && c > 0) || (mode == 1 && c < 0)) best = i;
    }

    Expr** out = malloc(sizeof(Expr*) * n);
    size_t nout = 0;
    for (size_t i = 0; i < n; i++)
        if (expr_compare(keys[i], keys[best]) == 0)
            out[nout++] = expr_copy(coll->data.function.args[i]);

    for (size_t i = 0; i < n; i++) expr_free(keys[i]);
    free(keys);

    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, nout);
    free(out);
    return result;
}

Expr* builtin_maximal_by(Expr* res) { return maximal_minimal_by(res, 0); }
Expr* builtin_minimal_by(Expr* res) { return maximal_minimal_by(res, 1); }

/* ------------------- ReverseSort / ReverseSortBy ------------------- */

/* Reverse the top-level elements of a freshly-built (owned) expression, in
 * place. Handles both representations of a list.
 *
 * The NDArray arm is not an optimisation, it is the correctness fix: an NDArray
 * is EXPR_NDARRAY, not EXPR_FUNCTION, so the early return silently did nothing
 * and ReverseSort[NDArray[{3.,1.,4.,1.,5.}]] answered {1.,1.,3.,4.,5.} --
 * builtin_sort's ascending result, unreversed and indistinguishable from Sort.
 * A live wrong answer on the visible NDArray surface, and it would have become
 * one on ordinary lists the moment ReverseSort was marked packed-aware.
 * ReverseSortBy shares this helper and had the same bug. */
static Expr* reverse_top_level(Expr* e) {
    if (!e) return e;

    if (is_ndarray(e)) {
        int rank = e->data.ndarray.rank;
        const int64_t* dims = e->data.ndarray.dims;
        size_t blocks = (size_t)dims[0];
        if (blocks < 2) return e;
        size_t rowelts = 1;
        for (int i = 1; i < rank; i++) rowelts *= (size_t)dims[i];
        size_t rowbytes = rowelts * ndt_elem_size(e->data.ndarray.dtype);
        char* data = (char*)e->data.ndarray.data;
        char* tmp = malloc(rowbytes);
        if (!tmp) return e;                 /* cannot reverse: leave as built */
        for (size_t i = 0; i < blocks / 2; i++) {
            char* lo = data + i * rowbytes;
            char* hi = data + (blocks - 1 - i) * rowbytes;
            memcpy(tmp, lo, rowbytes);
            memcpy(lo, hi, rowbytes);
            memcpy(hi, tmp, rowbytes);
        }
        free(tmp);
        return e;
    }

    if (e->type != EXPR_FUNCTION) return e;
    Expr** a = e->data.function.args;
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n / 2; i++) { Expr* t = a[i]; a[i] = a[n - 1 - i]; a[n - 1 - i] = t; }
    return e;
}

/* Re-head an owned call node's arguments onto `head` for delegation.
 *
 * builtin_reverse_sort used to hand its own `res` -- head ReverseSort -- straight
 * to builtin_sort. That is fine while the ascending routine only looks at the
 * arguments, and WRONG the moment it degrades: ndarray_delist_and_reeval rebuilds
 * the call from `res`'s head, so a rank-2 NDArray made it evaluate
 * ReverseSort[plainList], which is already descending, and the outer
 * reverse_top_level then flipped it back to ascending.
 *
 * Delegating under the correct head removes the whole class: whatever the
 * ascending routine does internally, including re-evaluating itself, it does as
 * Sort. Caller owns the returned node. */
static Expr* sort_call_as(const Expr* res, const char* head) {
    size_t n = res->data.function.arg_count;
    Expr** args = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) args[i] = expr_copy(res->data.function.args[i]);
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    free(args);
    return call;
}

/* ReverseSort[list] / ReverseSort[assoc] — descending order (by value for an
 * association), i.e. Reverse of Sort. */
Expr* builtin_reverse_sort(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    Expr* call = sort_call_as(res, "Sort");
    Expr* asc = builtin_sort(call);   /* borrows call, returns a new expr */
    expr_free(call);
    return reverse_top_level(asc);
}

/* ReverseSortBy[coll, f] — descending by f (of each value for an association). */
Expr* builtin_reverse_sort_by(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* call = sort_call_as(res, "SortBy");
    Expr* asc = builtin_sort_by(call);
    expr_free(call);
    return reverse_top_level(asc);
}

/* ------------------- TakeLargest / TakeSmallest (+By) ------------------- */

/* Take the `nreq` extreme elements of `coll` ranked by a key.
 *   largest = true  -> the largest, returned in descending order
 *   largest = false -> the smallest, returned in ascending order
 * The key is `f[subject]` when f != NULL, else `subject` itself, where
 * `subject` is each list element, or each value for an association. For an
 * association the result is an association of the corresponding entries. */
static Expr* take_extreme(Expr* coll, Expr* f, int64_t nreq, bool largest) {
    if (coll->type != EXPR_FUNCTION || nreq < 0) return NULL;
    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;

    SortByPair* pairs = malloc(sizeof(SortByPair) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* elem = coll->data.function.args[i];
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        pairs[i].payload = expr_copy(elem);
        if (f) {
            Expr* call_args[1] = { expr_copy(subject) };
            Expr* call = expr_new_function(expr_copy(f), call_args, 1);
            pairs[i].key = evaluate(call);
            expr_free(call);
        } else {
            pairs[i].key = expr_copy(subject);
        }
    }

    qsort(pairs, n, sizeof(SortByPair), sortby_pair_cmp);  /* ascending by key */

    size_t k = ((int64_t)n < nreq) ? n : (size_t)nreq;
    Expr** out = malloc(sizeof(Expr*) * (k ? k : 1));
    for (size_t j = 0; j < k; j++) {
        /* largest: walk from the top (descending); smallest: from the bottom. */
        size_t src = largest ? (n - 1 - j) : j;
        out[j] = expr_copy(pairs[src].payload);
    }
    for (size_t i = 0; i < n; i++) { expr_free(pairs[i].payload); expr_free(pairs[i].key); }
    free(pairs);

    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, k);
    free(out);
    return result;
}

Expr* builtin_take_largest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    /* O(n log k) through a bounded heap rather than the full O(n log n) sort
     * take_extreme does below -- and without one Expr per element to sort. */
    if (is_ndarray(res->data.function.args[0])) {
        Expr* nd = ndstruct_take_extreme(res, true);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
    Expr* nexpr = res->data.function.args[1];
    if (nexpr->type != EXPR_INTEGER) return NULL;
    return take_extreme(res->data.function.args[0], NULL, nexpr->data.integer, true);
}

Expr* builtin_take_smallest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    if (is_ndarray(res->data.function.args[0])) {   /* see builtin_take_largest */
        Expr* nd = ndstruct_take_extreme(res, false);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
    Expr* nexpr = res->data.function.args[1];
    if (nexpr->type != EXPR_INTEGER) return NULL;
    return take_extreme(res->data.function.args[0], NULL, nexpr->data.integer, false);
}

Expr* builtin_take_largest_by(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* nexpr = res->data.function.args[2];
    if (nexpr->type != EXPR_INTEGER) return NULL;
    return take_extreme(res->data.function.args[0], res->data.function.args[1], nexpr->data.integer, true);
}

Expr* builtin_take_smallest_by(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* nexpr = res->data.function.args[2];
    if (nexpr->type != EXPR_INTEGER) return NULL;
    return take_extreme(res->data.function.args[0], res->data.function.args[1], nexpr->data.integer, false);
}

/* ------------------- RankedMin / RankedMax ------------------- */

/* A machine-double sort key for one real element, for the mixed / symbolic
 * ordering path: +/-HUGE_VAL for +/-Infinity, the exact value for a real atom,
 * else the numericalized value of a symbolic real (Pi, E, Sqrt[2], Pi + E, ...).
 * Returns false when e is not a real number -- a free symbol or a non-real
 * complex -- which is the signal that RankedMin has no definite result and must
 * stay unevaluated. */
static bool ranked_numeric_key(Expr* e, double* out) {
    if (is_infinity(e)) { *out = HUGE_VAL; return true; }
    if (is_minus_infinity(e)) { *out = -HUGE_VAL; return true; }
    /* A real atom (incl. a real-valued Complex[re, 0]) reads its value directly;
     * a symbolic real is numericalized to machine precision first. The giant
     * high-precision reals stay untouched -- only their double key is taken. */
    Expr* nv = is_real_numeric(e) ? NULL : numericalize(e, numeric_machine_spec());
    Expr* v = nv ? nv : e;
    bool ok = false;
    if (is_real_numeric(v)) {
        Expr* re; Expr* im;
        *out = is_complex(v, &re, &im) ? get_numeric_value(re) : get_numeric_value(v);
        ok = true;
    }
    if (nv) expr_free(nv);
    return ok;
}

/* Context + comparator for the index quickselect. Exactly one of `keys` (the
 * approximate double-key path) and expr_compare (the exact path, keys == NULL)
 * decides order; the original index is the stable tiebreak, so the result is a
 * strict total order over distinct indices either way. cmp(i, i) == 0 so the
 * Hoare partition treats an element as equal to itself. */
typedef struct { Expr** elem; const double* keys; } RankedCtx;

static int ranked_cmp(const RankedCtx* c, size_t i, size_t j) {
    if (i == j) return 0;
    if (c->keys) {
        if (c->keys[i] < c->keys[j]) return -1;
        if (c->keys[i] > c->keys[j]) return 1;
    } else {
        int r = expr_compare(c->elem[i], c->elem[j]);
        if (r != 0) return r;
    }
    return (i < j) ? -1 : 1;   /* stable: original position breaks value ties */
}

/* Quickselect the k-th (0-based) position of idx[0..m) under ranked_cmp, so that
 * idx[k] names the k-th smallest element. O(m) average; middle pivot as in
 * nd_select_kth (ndreduce.c). Only idx[k] is guaranteed settled on return. */
static void ranked_select_idx(size_t* idx, size_t m, size_t k, const RankedCtx* c) {
    size_t lo = 0, hi = m - 1;
    while (lo < hi) {
        size_t pivot = idx[lo + (hi - lo) / 2];
        size_t i = lo, j = hi;
        while (i <= j) {                                  /* Hoare partition */
            while (ranked_cmp(c, idx[i], pivot) < 0) i++;
            while (ranked_cmp(c, idx[j], pivot) > 0) j--;
            if (i <= j) {
                size_t t = idx[i]; idx[i] = idx[j]; idx[j] = t;
                i++;
                if (j == 0) break;                        /* size_t underflow guard */
                j--;
            }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else break;                                       /* j < k < i: idx[k] settled */
    }
}

/* RankedMin[list, n] core (RankedMax delegates with is_max = true, which negates
 * n). Selects and returns the EXACT r-th ascending element, or NULL to leave the
 * call unevaluated (empty list, |n| out of range, or a non-real element -- a
 * definite result needs every element to be a real number). */
static Expr* ranked_select(Expr* list, int64_t n0, bool is_max) {
    size_t m = list->data.function.arg_count;
    if (m == 0) return NULL;

    int64_t n = is_max ? -n0 : n0;
    int64_t r = (n > 0) ? n : ((int64_t)m + n + 1);       /* ascending rank, 1-based */
    if (r < 1 || r > (int64_t)m) return NULL;             /* |n| out of range */

    Expr** elem = list->data.function.args;

    /* One scan chooses the path: exact (expr_compare, BigInt/Rational-safe) when
     * every element is a real numeric atom, else the double-key path (which also
     * validates that every element is a real number). */
    bool all_exact = true;
    for (size_t i = 0; i < m; i++)
        if (!is_real_numeric(elem[i])) { all_exact = false; break; }

    double* keys = NULL;
    if (!all_exact) {
        keys = malloc(sizeof(double) * m);
        if (!keys) return NULL;
        for (size_t i = 0; i < m; i++)
            if (!ranked_numeric_key(elem[i], &keys[i])) { free(keys); return NULL; }
    }

    size_t* idx = malloc(sizeof(size_t) * m);
    if (!idx) { free(keys); return NULL; }
    for (size_t i = 0; i < m; i++) idx[i] = i;

    RankedCtx ctx = { elem, keys };
    ranked_select_idx(idx, m, (size_t)(r - 1), &ctx);
    Expr* out = expr_copy(elem[idx[r - 1]]);

    free(idx);
    free(keys);
    return out;
}

/* RankedMin[list, n] -- the n-th smallest element (n<0: the n-th largest).
 * RankedMin[list, 1] == Min[list], RankedMin[list, -1] == Max[list]. Yields a
 * definite result iff every element is a real number; +/-Infinity rank as
 * +/-inf. An NDArray argument takes the buffer fast path (ndreduce.c). Returns
 * NULL (unevaluated) for a non-integer / out-of-range n or a symbolic element. */
Expr* builtin_ranked_min(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return builtin_arg_error("RankedMin",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 2, 2);
    Expr* list  = res->data.function.args[0];
    Expr* nexpr = res->data.function.args[1];
    if (nexpr->type != EXPR_INTEGER || nexpr->data.integer == 0) return NULL;
    if (is_ndarray(list)) return ndred_ranked_min(res);
    if (!is_listq(list)) return NULL;
    return ranked_select(list, nexpr->data.integer, false);
}

/* RankedMax[list, n] -- the n-th largest element (n<0: the n-th smallest);
 * RankedMax[list, n] == RankedMin[list, -n]. See builtin_ranked_min. */
Expr* builtin_ranked_max(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return builtin_arg_error("RankedMax",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 2, 2);
    Expr* list  = res->data.function.args[0];
    Expr* nexpr = res->data.function.args[1];
    if (nexpr->type != EXPR_INTEGER || nexpr->data.integer == 0) return NULL;
    if (is_ndarray(list)) return ndred_ranked_max(res);
    if (!is_listq(list)) return NULL;
    return ranked_select(list, nexpr->data.integer, true);
}

Expr* builtin_orderedq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];
    if (list->type != EXPR_FUNCTION) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;
    size_t count = list->data.function.arg_count;
    
    if (count <= 1) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* old_p = current_sort_p;
    current_sort_p = p;
    
    bool ordered = true;
    for (size_t i = 0; i < count - 1; i++) {
        Expr* ea = list->data.function.args[i];
        Expr* eb = list->data.function.args[i+1];
        
        int cmp = 0;
        if (current_sort_p == NULL) {
            cmp = expr_compare(ea, eb);
        } else {
            Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
            Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
            Expr* result = evaluate(call);
            expr_free(call);
            
            cmp = -1; // Default to "in order" (ea <= eb)
            if (result->type == EXPR_INTEGER) {
                if (result->data.integer == 1) cmp = -1;
                else if (result->data.integer == 0) cmp = 0;
                else if (result->data.integer == -1) cmp = 1;
            } else if (result->type == EXPR_SYMBOL) {
                if (result->data.symbol.name == SYM_True) cmp = -1;
                else if (result->data.symbol.name == SYM_False) cmp = 1;
            }
            expr_free(result);
        }
        
        if (cmp > 0) {
            ordered = false;
            break;
        }
    }
    
    current_sort_p = old_p;
    
    return expr_new_symbol(ordered ? "True" : "False");
}

/* Order[e1, e2] -- the user-facing surface of the canonical comparator that
 * every sorting routine (Sort, SortBy, OrderedQ, ...) is built on. Returns 1 if
 * e1 is before e2 in canonical order, -1 if after, 0 if identical. Compares
 * structurally, not by numerical value: Order[6, Pi] is 1 (Integer before the
 * symbol Pi) whereas Order[6, N[Pi]] is -1 (two numeric atoms, by value).
 *
 * expr_compare returns negative when e1 sorts before e2; Order's convention is
 * +1 in that case, so the sign is inverted. */
Expr* builtin_order(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) {
        return builtin_arg_error("Order",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 2, 2);
    }
    int cmp = expr_compare(res->data.function.args[0],
                           res->data.function.args[1]);
    return expr_new_integer(cmp < 0 ? 1 : (cmp > 0 ? -1 : 0));
}

/* ------------------- Ordering ------------------- */

/* File-static context for the index comparator, saved/restored around each sort
 * so a custom ordering function that itself calls Ordering nests correctly. */
static Expr** ordering_subjects = NULL;
static Expr*  ordering_p        = NULL;

/* qsort over 0-based indices: order by the subject each points at, ties broken by
 * the original index ascending. That tie-break is what turns libc's unstable
 * qsort into a STABLE argsort, so equal elements keep their input order --
 * Ordering[{2,2,1}] is {3,1,2}, and Ordering[list, 1] names the first minimum. */
static int ordering_index_compare(const void* a, const void* b) {
    int64_t ia = *(const int64_t*)a;
    int64_t ib = *(const int64_t*)b;
    int c = cmp_with_p(ordering_subjects[ia], ordering_subjects[ib], ordering_p);
    if (c != 0) return c;
    return (ia < ib) ? -1 : 1;
}

/* Ordering[list] -- the permutation of 1-based positions that sorts list, so that
 * list[[Ordering[list]]] === Sort[list]. Ordering[list, seq] == Take[Ordering[list],
 * seq] (seq an Integer n / -n, a {m,n[,s]} span, UpTo[k], or All). Ordering[list,
 * seq, p] orders by the function p as in Sort[list, p]. The result is always a
 * List of integers, regardless of list's head; for an Association the ordering is
 * of the VALUES. */
Expr* builtin_ordering(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3)
        return builtin_arg_error("Ordering", argc, 1, 3);

    Expr* list = res->data.function.args[0];

    /* Packed/visible NDArray fast path: argsort straight off the buffer. It
     * declines (rank >= 2, complex dtype, or a custom comparator) by
     * materialising and re-evaluating, which re-enters this builtin on a plain
     * List, so a NULL return can only mean "genuinely unevaluated". */
    if (is_ndarray(list)) {
        Expr* nd = ndstruct_ordering(res);
        if (nd) return nd;
    }

    if (list->type != EXPR_FUNCTION) return NULL;  /* Ordering[atom]: leave unevaluated */

    Expr* seq = (argc >= 2) ? res->data.function.args[1] : NULL;
    Expr* p   = (argc == 3) ? res->data.function.args[2] : NULL;

    size_t n = list->data.function.arg_count;
    bool assoc = is_association(list);

    /* Borrowed subject pointers -- the comparator only reads them, never frees.
     * For an Association order by each entry's VALUE (the Rule's second arg);
     * otherwise by the element itself. */
    Expr** subjects = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* elem = list->data.function.args[i];
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subjects[i] = elem->data.function.args[1];
        else
            subjects[i] = elem;
    }

    int64_t* idx = malloc(sizeof(int64_t) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) idx[i] = (int64_t)i;

    Expr** old_subjects = ordering_subjects;
    Expr*  old_p        = ordering_p;
    ordering_subjects = subjects;
    ordering_p = p;
    qsort(idx, n, sizeof(int64_t), ordering_index_compare);
    ordering_subjects = old_subjects;
    ordering_p = old_p;

    /* Which ranks of the permutation to emit: all of it, or a Take-style slice. */
    int64_t* sel = NULL;      /* 1-based ranks into the permutation */
    size_t sel_count = n;
    if (seq != NULL) {
        if (!get_seq_spec_indices(seq, (int64_t)n, &sel, &sel_count)) {
            free(subjects);
            free(idx);
            return NULL;      /* malformed / out-of-range spec: leave unevaluated */
        }
    }

    Expr** out = (sel_count > 0) ? malloc(sizeof(Expr*) * sel_count) : NULL;
    for (size_t i = 0; i < sel_count; i++) {
        int64_t rank = sel ? sel[i] : (int64_t)(i + 1);   /* 1-based rank */
        out[i] = expr_new_integer(idx[rank - 1] + 1);      /* 1-based position */
    }

    free(sel);
    free(subjects);
    free(idx);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, sel_count);
    if (out) free(out);
    return pack_offer(result);   /* a list of positions packs; often used as Part indices next */
}
