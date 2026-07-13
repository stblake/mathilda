/* flint_qqbar.c — exact algebraic-number canonicalisation for RootReduce.
 *
 * See flint_qqbar.h for the contract. The engine is FLINT's `qqbar` (exact
 * real/complex algebraic numbers). A constant algebraic-number Expr is folded
 * bottom-up into a single qqbar, whose minimal polynomial, degree, rationality,
 * quadratic-radical form and root index then produce the WL-faithful canonical
 * Expr: a rational, a quadratic radical, or Root[Function[minpoly&], k].
 */

#include "flint_qqbar.h"

#include "sym_names.h"
#include "eval.h"

#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/qqbar.h>

/* Intermediate qqbar degree above which we give up and fall back (identity /
 * the parametric engine). WL's degree-21 examples stay comfortably under. */
#define QQBAR_DEGREE_CAP 120
/* Max distinct algebraic generators collected for the NumberField method. */
#define QQBAR_MAX_ATOMS  12

/* ------------------------------------------------------------------ */
/*  Small Expr / FLINT helpers                                         */
/* ------------------------------------------------------------------ */

static int head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static Expr* expr_from_fmpz(const fmpz_t z) {
    mpz_t m; mpz_init(m); fmpz_get_mpz(m, z);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
    mpz_clear(m);
    return r;
}

static Expr* expr_from_fmpq(const fmpq_t q) {
    if (fmpz_is_one(fmpq_denref(q))) return expr_from_fmpz(fmpq_numref(q));
    Expr* args[2] = { expr_from_fmpz(fmpq_numref(q)), expr_from_fmpz(fmpq_denref(q)) };
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

/* Fill an fmpz from an integer-like Expr; 1 on success. */
static int fmpz_from_intlike(const Expr* e, fmpz_t out) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, (slong)e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

/* Fill an fmpq from an integer-like or Rational[p,q] Expr; 1 on success. */
static int fmpq_from_expr(const Expr* e, fmpq_t out) {
    fmpz_t z; fmpz_init(z);
    if (fmpz_from_intlike(e, z)) {
        fmpz_set(fmpq_numref(out), z); fmpz_one(fmpq_denref(out));
        fmpz_clear(z); return 1;
    }
    fmpz_clear(z);
    if (head_is(e, "Rational") && e->data.function.arg_count == 2) {
        fmpz_t n, d; fmpz_init(n); fmpz_init(d);
        int ok = fmpz_from_intlike(e->data.function.args[0], n) &&
                 fmpz_from_intlike(e->data.function.args[1], d);
        if (ok) fmpq_set_fmpz_frac(out, n, d);
        fmpz_clear(n); fmpz_clear(d);
        return ok;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Univariate integer poly Expr (Root body) -> fmpz_poly              */
/* ------------------------------------------------------------------ */

/* Treat `Slot[1]` (slot form) or the symbol named `var` (2-arg Function form)
 * as the polynomial variable. Only integer coefficients are accepted; a
 * rational coefficient or a foreign symbol makes the build fail (return 0),
 * which just leaves the Root object unconverted. */
static int is_poly_var(const Expr* e, const char* var) {
    if (var && e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, var) == 0;
    return head_is(e, "Slot") && e->data.function.arg_count == 1 &&
           e->data.function.args[0]->type == EXPR_INTEGER &&
           e->data.function.args[0]->data.integer == 1;
}

static int build_fmpz_poly(const Expr* e, const char* var, fmpz_poly_t out) {
    if (!e) return 0;
    fmpz_t z; fmpz_init(z);
    if (fmpz_from_intlike(e, z)) { fmpz_poly_set_fmpz(out, z); fmpz_clear(z); return 1; }
    fmpz_clear(z);
    if (is_poly_var(e, var)) {
        fmpz_poly_zero(out);                 /* out may be a reused temp */
        fmpz_poly_set_coeff_si(out, 1, 1);
        return 1;
    }
    if (e->type != EXPR_FUNCTION) return 0;

    size_t n = e->data.function.arg_count;
    if (head_is(e, "Plus")) {
        fmpz_poly_zero(out);
        fmpz_poly_t t; fmpz_poly_init(t);
        for (size_t i = 0; i < n; i++) {
            if (!build_fmpz_poly(e->data.function.args[i], var, t)) { fmpz_poly_clear(t); return 0; }
            fmpz_poly_add(out, out, t);
        }
        fmpz_poly_clear(t);
        return 1;
    }
    if (head_is(e, "Times")) {
        fmpz_poly_set_si(out, 1);
        fmpz_poly_t t; fmpz_poly_init(t);
        for (size_t i = 0; i < n; i++) {
            if (!build_fmpz_poly(e->data.function.args[i], var, t)) { fmpz_poly_clear(t); return 0; }
            fmpz_poly_mul(out, out, t);
        }
        fmpz_poly_clear(t);
        return 1;
    }
    if (head_is(e, "Power") && n == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return 0;
        fmpz_poly_t b; fmpz_poly_init(b);
        if (!build_fmpz_poly(base, var, b)) { fmpz_poly_clear(b); return 0; }
        fmpz_poly_pow(out, b, (ulong)exp->data.integer);
        fmpz_poly_clear(b);
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  WL-faithful root ordering                                          */
/* ------------------------------------------------------------------ */

/* Ordering: real roots first, ascending by real part; then non-real roots by
 * qqbar's canonical root order. Matches WL's Root indexing on the (common) case
 * of real roots; exotic all-complex orderings may differ. */
static int wl_cmp_qq(const qqbar_t a, const qqbar_t b) {
    int ra = qqbar_is_real(a), rb = qqbar_is_real(b);
    if (ra != rb) return ra ? -1 : 1;
    if (ra) return qqbar_cmp_re(a, b);
    return qqbar_cmp_root_order(a, b);
}

/* Selection-sort an index permutation of `roots` into WL order (n <= cap). */
static void wl_sort_indices(qqbar_srcptr roots, slong n, slong* idx) {
    for (slong i = 0; i < n; i++) idx[i] = i;
    for (slong i = 0; i < n; i++) {
        slong best = i;
        for (slong j = i + 1; j < n; j++)
            if (wl_cmp_qq(roots + idx[j], roots + idx[best]) < 0) best = j;
        slong t = idx[i]; idx[i] = idx[best]; idx[best] = t;
    }
}

/* 1-based WL index of x among the roots of its minimal polynomial. */
static slong wl_root_index(const qqbar_t x) {
    slong d = qqbar_degree(x);
    if (d <= 0) return 1;
    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, QQBAR_POLY(x), QQBAR_ROOTS_IRREDUCIBLE);
    slong* idx = malloc(sizeof(slong) * (size_t)d);
    wl_sort_indices(roots, d, idx);
    slong k = 1;
    for (slong j = 0; j < d; j++)
        if (qqbar_equal(roots + idx[j], x)) { k = j + 1; break; }
    free(idx);
    _qqbar_vec_clear(roots, d);
    return k;
}

/* ------------------------------------------------------------------ */
/*  Expr -> qqbar (recursive field arithmetic)                         */
/* ------------------------------------------------------------------ */

static int to_qqbar(const Expr* e, qqbar_t out);

/* Root[Function[...], k] -> the k-th root (WL order) of the body polynomial. */
static int root_object_to_qqbar(const Expr* e, qqbar_t out) {
    size_t n = e->data.function.arg_count;
    if (n < 2) return 0;
    const Expr* fn = e->data.function.args[0];
    const Expr* ke = e->data.function.args[1];
    if (!head_is(fn, "Function") || ke->type != EXPR_INTEGER) return 0;

    const char* var = NULL;
    const Expr* body;
    if (fn->data.function.arg_count == 2 &&
        fn->data.function.args[0]->type == EXPR_SYMBOL) {
        var  = fn->data.function.args[0]->data.symbol.name;   /* Function[t, body] */
        body = fn->data.function.args[1];
    } else if (fn->data.function.arg_count == 1) {
        body = fn->data.function.args[0];                /* Function[body(Slot[1])] */
    } else {
        return 0;
    }

    fmpz_poly_t P; fmpz_poly_init(P);
    if (!build_fmpz_poly(body, var, P)) { fmpz_poly_clear(P); return 0; }
    slong d = fmpz_poly_degree(P);
    slong k = (slong)ke->data.integer;
    if (d < 1 || k < 1 || k > d) { fmpz_poly_clear(P); return 0; }

    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, P, 0);
    slong* idx = malloc(sizeof(slong) * (size_t)d);
    wl_sort_indices(roots, d, idx);
    qqbar_set(out, roots + idx[k - 1]);
    free(idx);
    _qqbar_vec_clear(roots, d);
    fmpz_poly_clear(P);
    return 1;
}

static int to_qqbar(const Expr* e, qqbar_t out) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) { qqbar_set_si(out, (slong)e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT) {
        fmpz_t z; fmpz_init(z); fmpz_set_mpz(z, e->data.bigint);
        qqbar_set_fmpz(out, z); fmpz_clear(z); return 1;
    }
    if (e->type != EXPR_FUNCTION) return 0;   /* symbol / real / string */

    size_t n = e->data.function.arg_count;

    if (head_is(e, "Rational") && n == 2) {
        fmpq_t q; fmpq_init(q);
        int ok = fmpq_from_expr(e, q);
        if (ok) qqbar_set_fmpq(out, q);
        fmpq_clear(q);
        return ok;
    }
    if (head_is(e, "Complex") && n == 2) {
        qqbar_t re, im, ii; qqbar_init(re); qqbar_init(im); qqbar_init(ii);
        int ok = to_qqbar(e->data.function.args[0], re) &&
                 to_qqbar(e->data.function.args[1], im);
        if (ok) { qqbar_i(ii); qqbar_mul(im, im, ii); qqbar_add(out, re, im); }
        qqbar_clear(re); qqbar_clear(im); qqbar_clear(ii);
        return ok;
    }
    if (head_is(e, "Plus") || head_is(e, "Times")) {
        int is_plus = head_is(e, "Plus");
        qqbar_t acc, t; qqbar_init(acc); qqbar_init(t);
        if (is_plus) qqbar_set_si(acc, 0); else qqbar_set_si(acc, 1);
        int ok = 1;
        for (size_t i = 0; i < n && ok; i++) {
            if (!to_qqbar(e->data.function.args[i], t)) { ok = 0; break; }
            if (is_plus) qqbar_add(acc, acc, t); else qqbar_mul(acc, acc, t);
            if (qqbar_degree(acc) > QQBAR_DEGREE_CAP) { ok = 0; break; }
        }
        if (ok) qqbar_set(out, acc);
        qqbar_clear(acc); qqbar_clear(t);
        return ok;
    }
    if (head_is(e, "Power") && n == 2) {
        const Expr* be = e->data.function.args[0];
        const Expr* xe = e->data.function.args[1];
        qqbar_t base; qqbar_init(base);
        int ok = to_qqbar(be, base);
        if (ok) {
            if (xe->type == EXPR_INTEGER) {
                qqbar_pow_si(out, base, (slong)xe->data.integer);
            } else {
                fmpq_t ef; fmpq_init(ef);
                if (fmpq_from_expr(xe, ef)) qqbar_pow_fmpq(out, base, ef);
                else ok = 0;
                fmpq_clear(ef);
            }
        }
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) ok = 0;
        qqbar_clear(base);
        return ok;
    }
    if (head_is(e, "Sqrt") && n == 1) {
        qqbar_t base; qqbar_init(base);
        int ok = to_qqbar(e->data.function.args[0], base);
        if (ok) qqbar_sqrt(out, base);
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) ok = 0;
        qqbar_clear(base);
        return ok;
    }
    if (head_is(e, "Root")) return root_object_to_qqbar(e, out);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  qqbar -> canonical Expr (rational / quadratic radical / Root)      */
/* ------------------------------------------------------------------ */

/* Build the minimal polynomial of x as an Expr in Slot[1]. */
static Expr* minpoly_slot_expr(const qqbar_t x) {
    const fmpz_poly_struct* P = QQBAR_POLY(x);
    slong len = fmpz_poly_length(P);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
    size_t nt = 0;
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i < len; i++) {
        fmpz_poly_get_coeff_fmpz(c, P, i);
        if (fmpz_is_zero(c)) continue;
        Expr* mono;
        if (i == 0) {
            mono = NULL;
        } else {
            Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot),
                            (Expr*[]){ expr_new_integer(1) }, 1);
            mono = (i == 1) ? slot
                 : expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ slot, expr_new_integer((int64_t)i) }, 2);
        }
        Expr* ce = expr_from_fmpz(c);
        Expr* term;
        if (!mono)                          term = ce;
        else if (fmpz_is_one(c))            { term = mono; expr_free(ce); }
        else                                term = expr_new_function(expr_new_symbol(SYM_Times),
                                                        (Expr*[]){ ce, mono }, 2);
        terms[nt++] = term;
    }
    fmpz_clear(c);
    Expr* body = (nt == 0) ? expr_new_integer(0)
               : (nt == 1) ? terms[0]
               : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ body }, 1);
}

static Expr* qqbar_to_expr(const qqbar_t x) {
    if (qqbar_is_rational(x)) {
        fmpq_t q; fmpq_init(q); qqbar_get_fmpq(q, x);
        Expr* r = expr_from_fmpq(q); fmpq_clear(q);
        return r;
    }
    slong d = qqbar_degree(x);
    if (d == 2) {
        fmpz_t a, b, c, q; fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init(q);
        qqbar_get_quadratic(a, b, c, q, x, 0);        /* x = (a + b sqrt(c)) / q */
        Expr* half[2] = { expr_from_fmpz(c),
                          expr_new_function(expr_new_symbol(SYM_Rational),
                              (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) };
        Expr* sqrtc = expr_new_function(expr_new_symbol(SYM_Power), half, 2);
        Expr* bterm = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_from_fmpz(b), sqrtc }, 2);
        Expr* num = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_from_fmpz(a), bterm }, 2);
        Expr* invq = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_from_fmpz(q), expr_new_integer(-1) }, 2);
        Expr* frac = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ num, invq }, 2);
        fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(q);
        return eval_and_free(frac);
    }
    /* degree >= 3: Root[Function[minpoly&], k] (held). */
    Expr* fn = minpoly_slot_expr(x);
    slong k = wl_root_index(x);
    return expr_new_function(expr_new_symbol(SYM_Root),
               (Expr*[]){ fn, expr_new_integer((int64_t)k) }, 2);
}

/* ------------------------------------------------------------------ */
/*  Method -> "NumberField": common-number-field re-expression         */
/* ------------------------------------------------------------------ */

/* True if x is expressible in Q(alpha). */
static int in_field(const qqbar_t x, const qqbar_t alpha) {
    fmpq_poly_t f; fmpq_poly_init(f);
    int ok = qqbar_express_in_field(f, alpha, x, 100000, 0, 64);
    fmpq_poly_clear(f);
    return ok;
}

/* Collect the distinct algebraic generators (radicals, roots of unity, the
 * imaginary unit, Root objects) of `e` into `atoms`. Returns 0 on overflow /
 * conversion failure. Integers, rationals, Plus, Times and integer powers are
 * descended, not treated as atoms. */
static int collect_atoms(const Expr* e, qqbar_ptr atoms, int* na) {
    if (!e) return 1;
    if (e->type != EXPR_FUNCTION) return 1;   /* int/bigint/symbol: no atom */
    size_t n = e->data.function.arg_count;

    int is_generator = 0;
    if (head_is(e, "Power") && n == 2) {
        const Expr* xe = e->data.function.args[1];
        if (xe->type == EXPR_INTEGER)             /* integer power: descend base */
            return collect_atoms(e->data.function.args[0], atoms, na);
        is_generator = 1;                          /* fractional power: a radical */
    } else if (head_is(e, "Sqrt") || head_is(e, "Root")) {
        is_generator = 1;
    } else if (head_is(e, "Complex")) {
        is_generator = 1;                          /* the imaginary unit */
    } else if (head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!collect_atoms(e->data.function.args[i], atoms, na)) return 0;
        return 1;
    } else if (head_is(e, "Rational")) {
        return 1;
    } else {
        return 0;
    }

    if (is_generator) {
        if (*na >= QQBAR_MAX_ATOMS) return 0;
        qqbar_t v; qqbar_init(v);
        if (!to_qqbar(e, v)) { qqbar_clear(v); return 0; }
        if (qqbar_is_rational(v)) { qqbar_clear(v); return 1; }  /* not a real generator */
        for (int i = 0; i < *na; i++)
            if (qqbar_equal(atoms + i, v)) { qqbar_clear(v); return 1; }  /* dup */
        qqbar_set(atoms + *na, v);
        (*na)++;
        qqbar_clear(v);
    }
    return 1;
}

/* Evaluate the rational polynomial f at alpha via qqbar Horner. */
static void eval_poly_at(qqbar_t out, const fmpq_poly_t f, const qqbar_t alpha) {
    slong d = fmpq_poly_degree(f);
    qqbar_t acc; qqbar_init(acc); qqbar_set_si(acc, 0);
    fmpq_t c; fmpq_init(c);
    for (slong i = d; i >= 0; i--) {
        qqbar_mul(acc, acc, alpha);
        fmpq_poly_get_coeff_fmpq(c, f, i);
        qqbar_add_fmpq(acc, acc, c);
    }
    fmpq_clear(c);
    qqbar_set(out, acc);
    qqbar_clear(acc);
}

/* Re-derive `direct` through a single primitive element alpha of the field
 * generated by the atoms of `e` (the "AlgebraicNumber objects in a common
 * number field" route). Writes the reconstructed value to `out` and returns 1
 * only if it exactly matches `direct`; otherwise returns 0 and the caller uses
 * the direct value. Distinct computation, identical (verified) result. */
static int number_field_value(const Expr* e, const qqbar_t direct, qqbar_t out) {
    qqbar_ptr atoms = _qqbar_vec_init(QQBAR_MAX_ATOMS);
    int na = 0;
    int ok = collect_atoms(e, atoms, &na);
    if (ok && na == 0) { qqbar_set(out, direct); _qqbar_vec_clear(atoms, QQBAR_MAX_ATOMS); return 1; }

    qqbar_t alpha, cand, ca; qqbar_init(alpha); qqbar_init(cand); qqbar_init(ca);
    if (ok) qqbar_set(alpha, atoms + 0);
    for (int i = 1; i < na && ok; i++) {
        int found = 0;
        for (slong c = 1; c <= 8 && !found; c++) {
            qqbar_mul_si(ca, atoms + i, c);
            qqbar_add(cand, alpha, ca);
            if (in_field(atoms + i, cand) && in_field(alpha, cand)) {
                qqbar_set(alpha, cand); found = 1;
            }
        }
        if (!found) ok = 0;
    }

    if (ok) {
        fmpq_poly_t f; fmpq_poly_init(f);
        if (qqbar_express_in_field(f, alpha, direct, 100000, 0, 64)) {
            eval_poly_at(out, f, alpha);
            ok = qqbar_equal(out, direct);
        } else ok = 0;
        fmpq_poly_clear(f);
    }

    qqbar_clear(alpha); qqbar_clear(cand); qqbar_clear(ca);
    _qqbar_vec_clear(atoms, QQBAR_MAX_ATOMS);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Public API                                                         */
/* ------------------------------------------------------------------ */

int flint_qqbar_is_constant_algebraic(const Expr* e) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:  return 1;
        case EXPR_FUNCTION: break;
        default:           return 0;   /* real / string / free symbol */
    }
    size_t n = e->data.function.arg_count;
    if (head_is(e, "Rational")) return 1;
    if (head_is(e, "Root"))     return 1;   /* opaque; converter validates */
    if (head_is(e, "Complex") || head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!flint_qqbar_is_constant_algebraic(e->data.function.args[i])) return 0;
        return 1;
    }
    if (head_is(e, "Sqrt") && n == 1)
        return flint_qqbar_is_constant_algebraic(e->data.function.args[0]);
    if (head_is(e, "Power") && n == 2) {
        const Expr* xe = e->data.function.args[1];
        int exp_ok = xe->type == EXPR_INTEGER || head_is(xe, "Rational") ||
                     xe->type == EXPR_BIGINT;
        return exp_ok && flint_qqbar_is_constant_algebraic(e->data.function.args[0]);
    }
    return 0;
}

Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod method) {
    if (!flint_qqbar_is_constant_algebraic(e)) return NULL;
    qqbar_t val; qqbar_init(val);
    if (!to_qqbar(e, val)) { qqbar_clear(val); return NULL; }

    Expr* out;
    if (method == QQBAR_METHOD_NUMBERFIELD) {
        qqbar_t v2; qqbar_init(v2);
        out = qqbar_to_expr(number_field_value(e, val, v2) ? v2 : val);
        qqbar_clear(v2);
    } else {
        out = qqbar_to_expr(val);
    }
    qqbar_clear(val);
    return out;
}

int flint_qqbar_equal(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -1;
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -1;
    if (to_qqbar(a, qa) && to_qqbar(b, qb)) r = qqbar_equal(qa, qb) ? 1 : 0;
    qqbar_clear(qa); qqbar_clear(qb);
    return r;
}

int flint_qqbar_compare(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -2;
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -2;
    if (to_qqbar(a, qa) && to_qqbar(b, qb) && qqbar_is_real(qa) && qqbar_is_real(qb))
        r = qqbar_cmp_re(qa, qb);
    qqbar_clear(qa); qqbar_clear(qb);
    return r;
}

#else /* !USE_FLINT */

int   flint_qqbar_is_constant_algebraic(const Expr* e) { (void)e; return 0; }
Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod m) { (void)e; (void)m; return NULL; }
int   flint_qqbar_equal(const Expr* a, const Expr* b) { (void)a; (void)b; return -1; }
int   flint_qqbar_compare(const Expr* a, const Expr* b) { (void)a; (void)b; return -2; }

#endif /* USE_FLINT */
