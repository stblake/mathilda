#include "sort.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "assoc.h"
#include "assoc_struct.h"   /* OrderedQ compares an association's values */
#include "ndarray.h"   /* ndt_get for NDArray canonical ordering */
#include "ndstruct.h"  /* ndstruct_sort / ndstruct_ordering NDArray fast paths */
#include "pack.h"      /* pack_offer — a sorted machine list packs */
#include "common.h"    /* builtin_arg_error */
#include "take_drop.h" /* get_seq_spec_indices — Ordering[list, seq] reuses Take's spec */
#include "list/list_common.h" /* is_real_numeric / is_infinity / is_minus_infinity / is_listq */
#include "ndreduce.h"  /* ndred_ranked_min / ndred_ranked_max — NDArray order-statistic path */
#include "numeric.h"   /* numericalize / numeric_machine_spec — key for symbolic reals */
#include "sparse/sparse.h" /* EXPR_SPARSE (SparseArray raw objects) */
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
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
 * The order is the reference's, rule by rule (numbers first by value; then
 * products, powers and sums compared as polynomial terms from the last factor;
 * atoms before compound expressions; strings before symbols; symbols by the
 * string collation; other expressions by head, length, arguments).  The rules,
 * with examples, sit with their helpers in the "Canonical order" block below,
 * just above expr_compare.
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

/* ------------------------------------------------------------------ */
/* Canonical string order (WL collation)                              */
/* ------------------------------------------------------------------ */
/*
 * Strings compare on three levels, each lexicographic with the shorter
 * sequence first; a later level is consulted only on a tie:
 *
 *   L1  the PRIMARY characters alone -- digits < letters (case-folded) <
 *       non-ASCII (by code point) < `$`; punctuation and space are ignored.
 *   L2  the per-character mask punctuation = 0 / primary = 1, so leading
 *       punctuation sorts first: "!-0" < "!0" < "0", and "" < "!" < "!@".
 *   L3  per character: a punctuation mark's rank in PUNCT_ORDER below, or a
 *       letter's case (lowercase first); the first differing position wins.
 *
 * So "a" < "A" < "a!" < "a1", "a_b" < "a-b" < "ab" < "aB" < "Ab" < "ab.",
 * and ".hid" sorts beside "hid".  The model reproduces the reference sort
 * exactly on 5,233 random printable-ASCII strings.  Equal keys mean equal
 * strings; a final byte comparison guarantees that for malformed UTF-8.
 */

/* Punctuation (and space) in collation order. */
static const char PUNCT_ORDER[] = "!@#%^&*()_-+=~`{[}]|\\<>.,;\"?'/: ";

/* Collation weights of one character. */
typedef struct {
    uint32_t prim;   /* 0 = punctuation (ignored at L1), else primary weight */
    uint32_t tert;   /* L3 weight: punctuation rank, or 1 for an uppercase letter */
} CollChar;

/* Read the character at *p (UTF-8 decoded; an invalid byte stands for
 * itself), advance *p, and fill `out`.  Returns 0 at end of string. */
static int coll_next(const unsigned char** p, CollChar* out) {
    const unsigned char* s = *p;
    unsigned c = *s;
    if (c == 0) return 0;
    if (c < 0x80) {
        *p = s + 1;
        out->tert = 0;
        if (c >= '0' && c <= '9') {
            out->prim = 1 + (c - '0');
        } else if (c >= 'a' && c <= 'z') {
            out->prim = 11 + (c - 'a');
        } else if (c >= 'A' && c <= 'Z') {
            out->prim = 11 + (c - 'A');
            out->tert = 1;
        } else if (c == '$') {
            out->prim = UINT32_MAX;
        } else {
            const char* r = strchr(PUNCT_ORDER, (int)c);
            out->prim = 0;
            /* Control characters rank after every listed mark. */
            out->tert = r ? (uint32_t)(r - PUNCT_ORDER) : (uint32_t)(sizeof PUNCT_ORDER + c);
        }
        return 1;
    }
    /* Multi-byte UTF-8 sequence -> code point; non-ASCII is primary. */
    int n = (c >= 0xF0) ? 3 : (c >= 0xE0) ? 2 : (c >= 0xC0) ? 1 : 0;
    uint32_t cp = n ? (c & (0x3Fu >> n)) : c;
    int i = 1;
    for (; i <= n && (s[i] & 0xC0) == 0x80; i++) cp = (cp << 6) | (s[i] & 0x3Fu);
    if (i <= n) { cp = c; i = 1; }      /* truncated sequence: the byte alone */
    *p = s + i;
    out->prim = 64 + cp;
    out->tert = 0;
    return 1;
}

static int string_compare_canonical(const char* s1, const char* s2) {
    const unsigned char* a;
    const unsigned char* b;
    CollChar ca, cb;
    int ha, hb;

    /* L1: primary characters only (the common, early-exit case). */
    a = (const unsigned char*)s1; b = (const unsigned char*)s2;
    for (;;) {
        while ((ha = coll_next(&a, &ca)) && ca.prim == 0) {}
        while ((hb = coll_next(&b, &cb)) && cb.prim == 0) {}
        if (!ha || !hb) {
            if (ha != hb) return ha ? 1 : -1;
            break;
        }
        if (ca.prim != cb.prim) return ca.prim < cb.prim ? -1 : 1;
    }

    /* L2: punctuation/primary mask, position by position. */
    a = (const unsigned char*)s1; b = (const unsigned char*)s2;
    for (;;) {
        ha = coll_next(&a, &ca);
        hb = coll_next(&b, &cb);
        if (!ha || !hb) {
            if (ha != hb) return ha ? 1 : -1;
            break;
        }
        int ma = ca.prim != 0, mb = cb.prim != 0;
        if (ma != mb) return ma < mb ? -1 : 1;
    }

    /* L3: punctuation rank / letter case.  The masks are equal, so each
     * position pairs like with like. */
    a = (const unsigned char*)s1; b = (const unsigned char*)s2;
    while (coll_next(&a, &ca) && coll_next(&b, &cb)) {
        if (ca.tert != cb.tert) return ca.tert < cb.tert ? -1 : 1;
    }

    int c = strcmp(s1, s2);
    return (c > 0) - (c < 0);
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

/* ------------------------------------------------------------------ */
/* Numbers                                                            */
/* ------------------------------------------------------------------ */

/* Rank of a real number's representation, for two equal values:
 * Integer < Real < Rational (1 < 1., 0.5 < 1/2). */
static int num_rank(const Expr* e) {
    if (expr_is_integer_like(e)) return 0;
    if (e->type == EXPR_REAL) return 1;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return 1;
#endif
    return 2;
}

/* Compare sa*a with sb*b, a and b real numeric atoms, sa/sb = +1 or -1.
 * Exact for two rational-likes; through doubles otherwise.  Value only. */
static int real_value_cmp(const Expr* a, int sa, const Expr* b, int sb) {
    /* Machine integers (the common case: coefficients, exponents). */
    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER &&
        a->data.integer != INT64_MIN && b->data.integer != INT64_MIN) {
        int64_t va = sa < 0 ? -a->data.integer : a->data.integer;
        int64_t vb = sb < 0 ? -b->data.integer : b->data.integer;
        return (va > vb) - (va < vb);
    }
    if (a->type == EXPR_REAL && b->type == EXPR_REAL) {
        double va = sa * a->data.real, vb = sb * b->data.real;
        return (va > vb) - (va < vb);
    }
    if (is_rational_like(a) && is_rational_like(b)) {
        mpz_t an, ad, bn, bd;
        rat_to_mpz_parts(a, an, ad);
        rat_to_mpz_parts(b, bn, bd);
        if (sa < 0) mpz_neg(an, an);
        if (sb < 0) mpz_neg(bn, bn);
        mpz_mul(an, an, bd);           /* denominators are positive */
        mpz_mul(bn, bn, ad);
        int c = mpz_cmp(an, bn);
        mpz_clear(an); mpz_clear(ad); mpz_clear(bn); mpz_clear(bd);
        return (c > 0) - (c < 0);
    }
    double va = sa * get_numeric_value(a);
    double vb = sb * get_numeric_value(b);
    return (va > vb) - (va < vb);
}

static int real_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (e->data.integer > 0) - (e->data.integer < 0);
    if (is_rational_like(e)) {
        mpz_t n, d;
        rat_to_mpz_parts(e, n, d);
        int s = mpz_sgn(n);
        mpz_clear(n); mpz_clear(d);
        return s;
    }
    double v = get_numeric_value(e);
    return (v > 0) - (v < 0);
}

static int real_cmp(const Expr* a, const Expr* b) {
    int c = real_value_cmp(a, 1, b, 1);
    if (c) return c;
    int ra = num_rank(a), rb = num_rank(b);
    return (ra > rb) - (ra < rb);
}

/* A number: a real numeric atom or a literal Complex.  Order by real part,
 * then by the magnitude of the imaginary part (a real number has none), then
 * negative imaginary part first:
 *   -1 < -1 - I < 0 < -I < I < 1 < 1 - I < 1 + I < 1 - 2 I < 1. */
static int number_cmp(const Expr* a, const Expr* b) {
    const Expr* ra = is_complex_number(a) ? a->data.function.args[0] : a;
    const Expr* rb = is_complex_number(b) ? b->data.function.args[0] : b;
    const Expr* ia = is_complex_number(a) ? a->data.function.args[1] : NULL;
    const Expr* ib = is_complex_number(b) ? b->data.function.args[1] : NULL;
    int c = real_cmp(ra, rb);
    if (c) return c;
    if (!ia || !ib) return ia ? 1 : ib ? -1 : 0;
    int sa = real_sign(ia), sb = real_sign(ib);
    c = real_value_cmp(ia, sa < 0 ? -1 : 1, ib, sb < 0 ? -1 : 1);   /* |im| */
    if (c) return c;
    if (sa != sb) return sa < sb ? -1 : 1;
    int qa = num_rank(ia), qb = num_rank(ib);
    return (qa > qb) - (qa < qb);
}

static bool is_number(const Expr* e) {
    return is_atomic_numeric(e) || is_complex_number(e);
}

/* ------------------------------------------------------------------ */
/* Products, sums, powers                                             */
/* ------------------------------------------------------------------ */

static bool has_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* View e as `coef * f1 * ... * fn` (Times) or `num + t1 + ... + tn` (Plus):
 * a leading number is split off; a non-`head` expression is one item. */
static void split_items(const Expr* const* self, const char* head,
                        const Expr** num, const Expr* const** items, size_t* n) {
    const Expr* e = *self;
    *num = NULL;
    if (has_head(e, head) && e->data.function.arg_count > 0) {
        const Expr* const* args = (const Expr* const*)e->data.function.args;
        size_t k = e->data.function.arg_count;
        if (k > 1 && is_number(args[0])) { *num = args[0]; args++; k--; }
        *items = args;
        *n = k;
    } else {
        *items = self;
        *n = 1;
    }
}

/* Compare two item lists from the LAST item backwards; the shorter list wins
 * a tie. */
static int items_cmp_rev(const Expr* const* ia, size_t na,
                         const Expr* const* ib, size_t nb) {
    size_t m = na < nb ? na : nb;
    for (size_t i = 1; i <= m; i++) {
        int c = expr_compare(ia[na - i], ib[nb - i]);
        if (c) return c;
    }
    return (na > nb) - (na < nb);
}

/* Times / Plus rule.  `unit` is the value an absent number stands for (1 for
 * a coefficient, 0 for a numeric term). */
static int assoc_cmp(const Expr* a, const Expr* b, const char* head, const Expr* unit) {
    const Expr *na, *nb;
    const Expr* const* ia;
    const Expr* const* ib;
    size_t ka, kb;
    split_items(&a, head, &na, &ia, &ka);
    split_items(&b, head, &nb, &ib, &kb);
    int c = items_cmp_rev(ia, ka, ib, kb);
    if (c) return c;
    if (na || nb) {
        c = number_cmp(na ? na : unit, nb ? nb : unit);
        if (c) return c;
    }
    bool ha = has_head(a, head), hb = has_head(b, head);
    return ha == hb ? 0 : ha ? 1 : -1;
}

/* Power rule: bases, then exponents; a non-power is base^1. */
static int power_cmp(const Expr* a, const Expr* b, const Expr* one) {
    bool pa = has_head(a, SYM_Power) && a->data.function.arg_count == 2;
    bool pb = has_head(b, SYM_Power) && b->data.function.arg_count == 2;
    const Expr* ba = pa ? a->data.function.args[0] : a;
    const Expr* bb = pb ? b->data.function.args[0] : b;
    int c = expr_compare(ba, bb);
    if (c) return c;
    c = expr_compare(pa ? a->data.function.args[1] : one,
                     pb ? b->data.function.args[1] : one);
    if (c) return c;
    return pa == pb ? 0 : pa ? 1 : -1;
}

/* Symbols by the string collation of their names without context, then by
 * full name. */
static const char* short_name(const char* s) {
    const char* p = strrchr(s, '`');
    return p ? p + 1 : s;
}

static int symbol_cmp(const Expr* a, const Expr* b) {
    const char* x = a->data.symbol.name;
    const char* y = b->data.symbol.name;
    if (x == y) return 0;
    int c = string_compare_canonical(short_name(x), short_name(y));
    if (c) return c;
    c = strcmp(x, y);
    return (c > 0) - (c < 0);
}

static int opaque_cmp(const Expr* a, const Expr* b);

static const Expr* infinity_alias(const Expr* e) {
    static Expr* k_inf = NULL;
    static Expr* k_cinf = NULL;
    if (e->type != EXPR_SYMBOL) return e;
    const char* n = e->data.symbol.name;
    if (n != SYM_Infinity && n != SYM_ComplexInfinity) return e;
    if (!k_inf) {
        Expr* one[1];
        one[0] = expr_new_integer(1);
        k_inf = expr_new_function(expr_new_symbol(SYM_DirectedInfinity), one, 1);
        k_cinf = expr_new_function(expr_new_symbol(SYM_DirectedInfinity), NULL, 0);
    }
    return n == SYM_Infinity ? k_inf : k_cinf;
}

int expr_compare(const Expr* a, const Expr* b) {
    if (a == b) return 0;
    if (!a) return -1;
    if (!b) return 1;

    /* Fast paths for the two commonest shapes, each the exact outcome of the
     * general rules below:
     *   - two machine integers: the numeric tier, by value;
     *   - two calls with the SAME symbolic head that is not an arithmetic or
     *     number head (Times, Plus, Power, Complex, Rational): no packed,
     *     sparse, infinity or number tier applies, so it is step 6 directly
     *     -- equal heads, then length, then the arguments in order.
     * Sorting a large sum compares g[x[i], x[j]] kernels millions of times;
     * each nested level used to walk the whole tier preamble. */
    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER)
        return (a->data.integer > b->data.integer) - (a->data.integer < b->data.integer);
    if (a->type == EXPR_FUNCTION && b->type == EXPR_FUNCTION) {
        const Expr* ha = a->data.function.head;
        const Expr* hb = b->data.function.head;
        if (ha && hb && ha->type == EXPR_SYMBOL && hb->type == EXPR_SYMBOL &&
            ha->data.symbol.name == hb->data.symbol.name) {
            const char* h = ha->data.symbol.name;
            if (h != SYM_Times && h != SYM_Plus && h != SYM_Power &&
                h != SYM_Complex && h != SYM_Rational) {
                size_t ac = a->data.function.arg_count;
                size_t bc = b->data.function.arg_count;
                if (ac != bc) return (ac < bc) ? -1 : 1;
                for (size_t i = 0; i < ac; i++) {
                    int c = expr_compare(a->data.function.args[i], b->data.function.args[i]);
                    if (c != 0) return c;
                }
                return 0;
            }
        }
    }

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

    /* 0b. SparseArray raw objects sort after numbers and before everything
     * else; two arrays order by dims, then size, then contents (sp_compare). */
    if (a->type == EXPR_SPARSE || b->type == EXPR_SPARSE) {
        if (a->type == EXPR_SPARSE && b->type == EXPR_SPARSE)
            return sp_compare(a->data.sparse, b->data.sparse);
        const Expr* other = a->type == EXPR_SPARSE ? b : a;
        int sign = a->type == EXPR_SPARSE ? 1 : -1;
        if (is_atomic_numeric(other) || is_complex_number(other)) return sign;
        return -sign;
    }

    /* Infinity and ComplexInfinity are symbols here; they order as the
     * DirectedInfinity[1] and DirectedInfinity[] they stand for. */
    a = infinity_alias(a);
    b = infinity_alias(b);
    if (a == b) return 0;

    /* NDArray and CompiledFunction objects sort after everything. */
    if (a->type == EXPR_NDARRAY || a->type == EXPR_COMPILED ||
        b->type == EXPR_NDARRAY || b->type == EXPR_COMPILED)
        return opaque_cmp(a, b);

    /* 1. Numbers first, by value. */
    bool a_num = is_number(a), b_num = is_number(b);
    if (a_num && b_num) return number_cmp(a, b);
    if (a_num) return -1;
    if (b_num) return 1;

    /* 2-4. Products, then powers, then sums, as polynomial terms: a power
     * is compared before a sum, so x^2 < a + x ((a + x)^1 against x^2). */
    static Expr* k_zero = NULL;
    static Expr* k_one = NULL;
    if (!k_one) { k_zero = expr_new_integer(0); k_one = expr_new_integer(1); }
    if (has_head(a, SYM_Times) || has_head(b, SYM_Times))
        return assoc_cmp(a, b, SYM_Times, k_one);
    if (has_head(a, SYM_Power) || has_head(b, SYM_Power))
        return power_cmp(a, b, k_one);
    if (has_head(a, SYM_Plus) || has_head(b, SYM_Plus))
        return assoc_cmp(a, b, SYM_Plus, k_zero);

    /* 5. Atoms before compound expressions; strings before symbols. */
    if (a->type == EXPR_STRING && b->type == EXPR_STRING)
        return string_compare_canonical(a->data.string, b->data.string);
    if (a->type == EXPR_STRING) return -1;
    if (b->type == EXPR_STRING) return 1;
    if (a->type == EXPR_SYMBOL && b->type == EXPR_SYMBOL) return symbol_cmp(a, b);
    if (a->type == EXPR_SYMBOL) return -1;
    if (b->type == EXPR_SYMBOL) return 1;

    /* 6. Head, then length, then the arguments left to right. */
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
    return 0;
}

/* NDArray / CompiledFunction ordering: after every other expression. */
static int opaque_cmp(const Expr* a, const Expr* b) {
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

/* qsort adaptor over Expr* slots, using the canonical order. */
static int sort_cmp_expr_ptrs(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}

/* Canonically sort args[0..n) in place (the Orderless canonicalisation).
 *
 * Adaptive: one pass finds the maximal ascending runs, which are then merged
 * pairwise, bottom-up -- O(n log r) comparisons for r runs instead of
 * O(n log n).  Orderless input is very often a concatenation of sorted pieces:
 * Plus of k canonical sums flattens to k sorted runs (13 runs of 10395 terms
 * in the issue #112 Dirac trace), where a full qsort re-compared every term
 * log2(135135) ~ 17 times.  A sorted input costs n - 1 comparisons; input with
 * many short runs (r > n/8, i.e. effectively unordered) goes to qsort, whose
 * constant factor is lower there.  Equal elements are expr_eq-equivalent under
 * the canonical order, so stability is immaterial. */
void expr_orderless_sort(Expr** args, size_t n) {
    if (!args || n < 2) return;
    size_t cap = 16, nr = 0;
    size_t* runs = malloc(sizeof(size_t) * cap);
    if (!runs) { qsort(args, n, sizeof(Expr*), sort_cmp_expr_ptrs); return; }
    runs[nr++] = 0;
    bool unordered = false;                 /* too many runs: plain qsort */
    for (size_t i = 1; i < n; i++) {
        if (expr_compare(args[i - 1], args[i]) > 0) {
            if (nr >= 2 && nr > n / 8) { unordered = true; break; }
            if (nr == cap) {
                size_t* r2 = realloc(runs, sizeof(size_t) * cap * 2);
                if (!r2) { unordered = true; break; }
                runs = r2; cap *= 2;
            }
            runs[nr++] = i;
        }
    }
    if (!unordered && nr == 1) { free(runs); return; }
    Expr** tmp = unordered ? NULL : malloc(sizeof(Expr*) * n);
    if (!tmp) {
        free(runs);
        qsort(args, n, sizeof(Expr*), sort_cmp_expr_ptrs);
        return;
    }
    Expr** src = args;
    Expr** dst = tmp;
    while (nr > 1) {
        size_t out = 0;
        for (size_t k = 0; k < nr; k += 2) {
            size_t lo = runs[k];
            size_t mid = (k + 1 < nr) ? runs[k + 1] : n;
            size_t hi = (k + 2 < nr) ? runs[k + 2] : n;
            size_t i = lo, j = mid, o = lo;
            while (i < mid && j < hi)
                dst[o++] = (expr_compare(src[j], src[i]) < 0) ? src[j++] : src[i++];
            while (i < mid) dst[o++] = src[i++];
            while (j < hi)  dst[o++] = src[j++];
            runs[out++] = lo;
        }
        nr = out;
        Expr** t = src; src = dst; dst = t;
    }
    if (src != args) memcpy(args, src, sizeof(Expr*) * n);
    free(tmp);
    free(runs);
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
    Expr* result = eval_callback(call);
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

/*
 * Mathematica's merge sort for a user ordering function p, over a permutation.
 *
 * Sort[list, p] is not "any sort consistent with p": a strict p such as Greater
 * says False for equal elements, so ties come out in an order fixed by the
 * algorithm -- Sort[<|a -> 2, b -> 2, c -> 1|>, Greater] is <|b -> 2, a -> 2,
 * c -> 1|>. Mathematica 15 uses a top-down merge sort (left half floor(n/2))
 * whose merge keeps the LEFT element unless p[left, right] is False or -1;
 * this reproduces its output exactly, including Sort[list, -1 &] reversing the
 * list and Order-style p (0 on ties) being stable. It is O(n log n) calls of p.
 *
 * subj[i] is what p compares for element i; idx is the permutation, sorted in
 * place, tmp scratch of the same length. With p == NULL the canonical order is
 * used (stable; `desc` flips it). `swap` calls p[right, left] instead, which is
 * how ReverseSortBy[list, f, p] orders.
 */
static void p_merge_sort(size_t* idx, size_t* tmp, size_t n, Expr** subj,
                         Expr* p, bool swap, bool desc) {
    if (n < 2) return;
    size_t m = n / 2;
    p_merge_sort(idx, tmp, m, subj, p, swap, desc);
    p_merge_sort(idx + m, tmp, n - m, subj, p, swap, desc);
    size_t i = 0, j = m, k = 0;
    while (i < m && j < n) {
        Expr* l = subj[idx[i]];
        Expr* r = subj[idx[j]];
        bool left_first;
        if (p) left_first = (swap ? cmp_with_p(r, l, p) : cmp_with_p(l, r, p)) <= 0;
        else   left_first = desc ? expr_compare(l, r) >= 0 : expr_compare(l, r) <= 0;
        tmp[k++] = left_first ? idx[i++] : idx[j++];
    }
    while (i < m) tmp[k++] = idx[i++];
    while (j < n) tmp[k++] = idx[j++];
    memcpy(idx, tmp, n * sizeof(size_t));
}

/* A fresh identity permutation of n sorted by p_merge_sort. Caller frees. */
static size_t* p_sort_perm(Expr** subj, size_t n, Expr* p, bool swap, bool desc) {
    size_t* idx = malloc(sizeof(size_t) * (n ? n : 1));
    size_t* tmp = malloc(sizeof(size_t) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) idx[i] = i;
    p_merge_sort(idx, tmp, n, subj, p, swap, desc);
    free(tmp);
    return idx;
}

/* The thing an ordering compares for element `elem` of a collection: the value
 * of an association entry, else the element itself. Borrowed. */
static Expr* sort_subject(Expr* elem, bool assoc) {
    if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
        return elem->data.function.args[1];
    return elem;
}

/* `coll` with its elements permuted by perm (entries, for an association).
 * The head is kept. Caller owns the result. */
static Expr* permute_collection(const Expr* coll, const size_t* perm) {
    size_t n = coll->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) out[i] = expr_copy(coll->data.function.args[perm[i]]);
    Expr* r = expr_new_function(expr_copy(coll->data.function.head), out, n);
    free(out);
    return r;
}

/* Sort a collection's elements (an association's entries by value) with the
 * merge sort above. Caller owns the result. */
static Expr* sort_collection_p(const Expr* coll, Expr* p, bool swap, bool desc) {
    size_t n = coll->data.function.arg_count;
    bool assoc = is_association(coll);
    Expr** subj = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) subj[i] = sort_subject(coll->data.function.args[i], assoc);
    size_t* perm = p_sort_perm(subj, n, p, swap, desc);
    Expr* r = permute_collection(coll, perm);
    free(perm);
    free(subj);
    return r;
}

Expr* builtin_sort(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) {
        return NULL;
    }
    
    Expr* list = res->data.function.args[0];

    /* Sort[ndarray]: sort the flat buffer directly (see ndstruct.c). */
    if (is_ndarray(list)) return ndstruct_sort(res);

    if (list->type != EXPR_FUNCTION) return builtin_normal_error("Sort", res, 1);

    /* Sort[assoc] orders the entries by their values (keys follow). */
    if (res->data.function.arg_count == 1 && is_association(list))
        return assoc_sort_by_value(list);

    Expr* p = (res->data.function.arg_count == 2) ? res->data.function.args[1] : NULL;

    size_t count = list->data.function.arg_count;
    if (count == 0) {
        return expr_copy(list);
    }

    /* Sort[list, p] / Sort[assoc, p]: Mathematica's merge sort, so ties under a
     * strict p land where Mathematica puts them (an association compares its
     * values; the entries move with them). */
    if (p) return pack_offer(sort_collection_p(list, p, false, false));

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
typedef struct { Expr* payload; Expr* key; size_t pos; bool multi; } SortByPair;

/* SortBy order: by the key; ties by the canonical order of the element itself
 * (for an association, of the whole entry, i.e. by key) when sorting by a
 * single f, and by input position for a list {f1, f2, ...} -- both as in
 * Mathematica 15: SortBy[{{1, "b"}, {1, "a"}}, First] puts {1, "a"} first,
 * SortBy[{{1, "b"}, {1, "a"}}, {First}] leaves the order alone. */
static int sortby_pair_cmp(const void* a, const void* b) {
    const SortByPair* pa = (const SortByPair*)a;
    const SortByPair* pb = (const SortByPair*)b;
    int c = expr_compare(pa->key, pb->key);
    if (c) return c;
    if (!pa->multi) {
        c = expr_compare(pa->payload, pb->payload);
        if (c) return c;
    }
    return (pa->pos < pb->pos) ? -1 : (pa->pos > pb->pos);
}

/* Canonical order of sort subjects, ties by the element: the pre-order
 * SortBy[coll, f, p] merges from. Stable by position last. */
static int subject_pair_cmp(const void* a, const void* b) {
    const SortByPair* pa = (const SortByPair*)a;
    const SortByPair* pb = (const SortByPair*)b;
    int c = expr_compare(pa->key, pb->key);           /* key = the subject here */
    if (c) return c;
    c = expr_compare(pa->payload, pb->payload);
    if (c) return c;
    return (pa->pos < pb->pos) ? -1 : (pa->pos > pb->pos);
}

/* SortBy[coll, f, p] (swap = false) and ReverseSortBy[coll, f, p] (swap = true).
 *
 * Mathematica first puts the elements in canonical order of their subjects
 * (values, for an association), then merge-sorts them by p applied to the
 * f-values (p[b, a] for ReverseSortBy) -- see p_merge_sort. Reproduces e.g.
 * SortBy[<|a -> {1, z}, b -> {1, y}|>, First, Greater] = <|a -> .., b -> ..|>. */
static Expr* sort_by_with_p(const Expr* coll, Expr* f, Expr* p, bool swap) {
    size_t n = coll->data.function.arg_count;
    bool assoc = is_association(coll);
    SortByPair* pre = malloc(sizeof(SortByPair) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        pre[i].payload = coll->data.function.args[i];                /* borrowed */
        pre[i].key = sort_subject(coll->data.function.args[i], assoc); /* borrowed */
        pre[i].pos = i;
        pre[i].multi = false;
    }
    qsort(pre, n, sizeof(SortByPair), subject_pair_cmp);

    Expr** fv = malloc(sizeof(Expr*) * (n ? n : 1));     /* owned f-values */
    for (size_t i = 0; i < n; i++) {
        Expr* ca[1] = { expr_copy(pre[i].key) };
        Expr* call = expr_new_function(expr_copy(f), ca, 1);
        fv[i] = eval_callback(call);
        expr_free(call);
    }
    size_t* perm = p_sort_perm(fv, n, p, swap, false);

    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) out[i] = expr_copy(pre[perm[i]].payload);
    Expr* result = expr_new_function(expr_copy(coll->data.function.head), out, n);
    for (size_t i = 0; i < n; i++) expr_free(fv[i]);
    free(out); free(perm); free(fv); free(pre);
    return result;
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
    if (argc != 2 && argc != 3) return NULL;

    /* Without this a visible NDArray fell into the atom branch below and came
     * back UNSORTED. A visible NDArray is an atom here; see ndarray_visible_delist_retry. */
    { Expr* nd_ = ndarray_visible_delist_retry(res, 0, false); if (nd_) return nd_; }

    Expr* coll = res->data.function.args[0];
    Expr* f    = res->data.function.args[1];
    if (coll->type != EXPR_FUNCTION) return builtin_normal_error("SortBy", res, 1);

    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;
    if (n == 0) return expr_copy(coll);

    /* SortBy[coll, f, p]: compare the f-values with the ordering function p. */
    if (argc == 3) return sort_by_with_p(coll, f, res->data.function.args[2], false);

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
        pairs[i].pos = i;
        pairs[i].multi = multi;
        if (multi) {
            size_t nc = f->data.function.arg_count;
            Expr** keyparts = malloc(sizeof(Expr*) * (nc ? nc : 1));
            for (size_t c = 0; c < nc; c++) {
                Expr* ca[1] = { expr_copy(subject) };
                Expr* call = expr_new_function(expr_copy(f->data.function.args[c]), ca, 1);
                keyparts[c] = eval_callback(call);
                expr_free(call);
            }
            pairs[i].key = expr_new_function(expr_new_symbol(SYM_List), keyparts, nc);
            free(keyparts);
        } else {
            Expr* call_args[1] = { expr_copy(subject) };
            Expr* call = expr_new_function(expr_copy(f), call_args, 1);
            pairs[i].key = eval_callback(call);
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
    /* A visible NDArray is an atom here; see ndarray_visible_delist_retry. */
    { Expr* nd_ = ndarray_visible_delist_retry(res, 0, false); if (nd_) return nd_; }

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
        keys[i] = eval_callback(call);
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
    /* ReverseSort[assoc]: descending by value with equal values kept in input
     * order -- not Reverse[Sort[assoc]], which would flip the ties (Mathematica
     * 15: ReverseSort[<|a -> 2, b -> 2, c -> 1|>] is <|a -> 2, b -> 2, c -> 1|>).
     * With an ordering function it is Reverse[Sort[assoc, p]], below. */
    if (res->data.function.arg_count == 1 && is_association(res->data.function.args[0]))
        return sort_collection_p(res->data.function.args[0], NULL, false, true);
    Expr* call = sort_call_as(res, "Sort");
    Expr* asc = builtin_sort(call);   /* borrows call, returns a new expr */
    expr_free(call);
    return reverse_top_level(asc);
}

/* ReverseSortBy[coll, f] — descending by f (of each value for an association). */
Expr* builtin_reverse_sort_by(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* ReverseSortBy[coll, f, p] is SortBy[coll, f, p[#2, #1] &] (Mathematica
     * 15), which differs from Reverse[SortBy[coll, f, p]] on ties. */
    if (res->data.function.arg_count == 3) {
        Expr* coll = res->data.function.args[0];
        if (coll->type != EXPR_FUNCTION || is_ndarray(coll)) return NULL;
        if (coll->data.function.arg_count == 0) return expr_copy(coll);
        return sort_by_with_p(coll, res->data.function.args[1], res->data.function.args[2], true);
    }
    if (res->data.function.arg_count != 2) return NULL;
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
/* Rank two TakeLargest/TakeSmallest keys. Numbers rank by VALUE: 3 and 3. tie
 * here although canonical order separates them (by type), so that a tie is
 * settled by position, as for any other tie. */
static int take_key_cmp(const Expr* a, const Expr* b) {
    int c = expr_compare(a, b);
    if (c && a->type != b->type && is_atomic_numeric(a) && is_atomic_numeric(b) &&
        get_numeric_value(a) == get_numeric_value(b))
        return 0;
    return c;
}

/* Ties are STABLE: the element that comes first in the collection wins it and
 * is listed first. Each comparator is a total order (position last), so the
 * result does not depend on qsort, which is not stable. */
static int take_smallest_cmp(const void* a, const void* b) {
    const SortByPair* pa = (const SortByPair*)a;
    const SortByPair* pb = (const SortByPair*)b;
    int c = take_key_cmp(pa->key, pb->key);
    if (c) return c;
    return (pa->pos < pb->pos) ? -1 : (pa->pos > pb->pos);
}

static int take_largest_cmp(const void* a, const void* b) {
    const SortByPair* pa = (const SortByPair*)a;
    const SortByPair* pb = (const SortByPair*)b;
    int c = take_key_cmp(pb->key, pa->key);                 /* descending */
    if (c) return c;
    return (pa->pos < pb->pos) ? -1 : (pa->pos > pb->pos);
}

static Expr* take_extreme(Expr* coll, Expr* f, int64_t nreq, bool largest) {
    if (coll->type != EXPR_FUNCTION || nreq < 0) return NULL;
    bool assoc = is_association(coll);
    size_t n = coll->data.function.arg_count;

    SortByPair* pairs = malloc(sizeof(SortByPair) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        pairs[i].pos = i;
        pairs[i].multi = false;
        Expr* elem = coll->data.function.args[i];
        Expr* subject = elem;
        if (assoc && elem->type == EXPR_FUNCTION && elem->data.function.arg_count == 2)
            subject = elem->data.function.args[1];
        pairs[i].payload = expr_copy(elem);
        if (f) {
            Expr* call_args[1] = { expr_copy(subject) };
            Expr* call = expr_new_function(expr_copy(f), call_args, 1);
            pairs[i].key = eval_callback(call);
            expr_free(call);
        } else {
            pairs[i].key = expr_copy(subject);
        }
    }

    /* Best first: descending for largest, ascending for smallest. */
    qsort(pairs, n, sizeof(SortByPair), largest ? take_largest_cmp : take_smallest_cmp);

    size_t k = ((int64_t)n < nreq) ? n : (size_t)nreq;
    Expr** out = malloc(sizeof(Expr*) * (k ? k : 1));
    for (size_t j = 0; j < k; j++) out[j] = expr_copy(pairs[j].payload);
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
    /* A visible NDArray is an atom here; see ndarray_visible_delist_retry. */
    { Expr* nd_ = ndarray_visible_delist_retry(res, 0, false); if (nd_) return nd_; }
    Expr* nexpr = res->data.function.args[2];
    if (nexpr->type != EXPR_INTEGER) return NULL;
    return take_extreme(res->data.function.args[0], res->data.function.args[1], nexpr->data.integer, true);
}

Expr* builtin_take_smallest_by(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    /* A visible NDArray is an atom here; see ndarray_visible_delist_retry. */
    { Expr* nd_ = ndarray_visible_delist_retry(res, 0, false); if (nd_) return nd_; }
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

    bool assoc = assoc_is_wellformed(list);
    bool ordered = true;
    for (size_t i = 0; i < count - 1; i++) {
        /* An association is ordered by its values, never its keys:
         * OrderedQ[<|a -> 2, b -> 1|>] is False (assoc_struct.h). */
        Expr* ea = struct_part(list, i, assoc);
        Expr* eb = struct_part(list, i + 1, assoc);
        
        int cmp = 0;
        if (current_sort_p == NULL) {
            cmp = expr_compare(ea, eb);
        } else {
            Expr* args[2] = { expr_copy(ea), expr_copy(eb) };
            Expr* call = expr_new_function(expr_copy(current_sort_p), args, 2);
            Expr* result = eval_callback(call);
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
    if (p) {
        /* A user ordering function: Mathematica's merge sort, so a strict p
         * places ties exactly as Sort[list, p] does (Ordering[{1, 1, 0, 1},
         * All, Greater] is {4, 2, 1, 3}). */
        size_t* perm = p_sort_perm(subjects, n, p, false, false);
        for (size_t i = 0; i < n; i++) idx[i] = (int64_t)perm[i];
        free(perm);
    } else {
        for (size_t i = 0; i < n; i++) idx[i] = (int64_t)i;
        Expr** old_subjects = ordering_subjects;
        Expr*  old_p        = ordering_p;
        ordering_subjects = subjects;
        ordering_p = NULL;
        qsort(idx, n, sizeof(int64_t), ordering_index_compare);
        ordering_subjects = old_subjects;
        ordering_p = old_p;
    }

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
