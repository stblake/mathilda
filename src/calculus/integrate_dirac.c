/* integrate_dirac.c -- DiracDelta sifting under a definite integral.
 * See integrate_dirac.h for the contract and conventions. */

#include "integrate_dirac.h"
#include "expr.h"
#include "eval.h"
#include "../sym_names.h"

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

/* ---- small construction / evaluation helpers (mirrors integrate_symmetry.c) - */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn1(const char* head, Expr* a) {
    Expr* args[1] = { a }; return mk_fn(head, args, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b }; return mk_fn(head, args, 2);
}

/* Evaluate `call` (taking ownership) and return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Does e contain any DiracDelta[...] subexpression? */
static bool contains_dirac(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) {
        if (head_name_is(e, "DiracDelta")) return true;
        if (contains_dirac(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_dirac(e->data.function.args[i])) return true;
    }
    return false;
}

/* Is symbol `name` absent from e? */
static bool free_of(const Expr* e, const char* name) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) != 0;
    if (e->type == EXPR_FUNCTION) {
        if (!free_of(e->data.function.head, name)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!free_of(e->data.function.args[i], name)) return false;
    }
    return true;
}

/* Sign of an expression that may involve symbolic-but-decidable constants
 * (Pi, etc.): +1 (>0), -1 (<0), 0 (==0), or 2 (undecided).  Uses the
 * Positive/Negative/Equal predicates rather than Sign, which stays inert on a
 * symbolic constant like Pi (so Sign[Pi] would read as undecided). */
static int classify_sign(const Expr* diff) {
    Expr* z = eval_take(mk_fn2("Equal", cp(diff), expr_new_integer(0)));
    bool iszero = (z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_True);
    expr_free(z);
    if (iszero) return 0;
    Expr* p = eval_take(mk_fn1("Positive", cp(diff)));
    bool pos = (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_True);
    expr_free(p);
    if (pos) return 1;
    Expr* ng = eval_take(mk_fn1("Negative", cp(diff)));
    bool neg = (ng->type == EXPR_SYMBOL && ng->data.symbol.name == SYM_True);
    expr_free(ng);
    if (neg) return -1;
    return 2;
}

/* x0 -> value: evaluate ReplaceAll[e, x -> val]. */
static Expr* subst(const Expr* e, const Expr* x, Expr* val /*consumed*/) {
    Expr* rule = mk_fn2("Rule", cp(x), val);
    return eval_take(mk_fn2("ReplaceAll", cp(e), rule));
}

/* The captured fraction B(x0; lo, hi).  Returns a fresh Expr or NULL if the
 * interval membership of the impulse cannot be decided. */
static Expr* boundary_factor(const Expr* x0, const Expr* lo, const Expr* hi) {
    Expr* dlo_e = mk_fn2("Plus", cp(x0), mk_fn2("Times", expr_new_integer(-1), cp(lo)));
    Expr* dhi_e = mk_fn2("Plus", cp(hi), mk_fn2("Times", expr_new_integer(-1), cp(x0)));
    int dlo = classify_sign(dlo_e), dhi = classify_sign(dhi_e);
    expr_free(dlo_e); expr_free(dhi_e);

    if (dlo != 2 && dhi != 2) {
        /* Both endpoints decided (numeric limits). */
        if (dlo == 1 && dhi == 1) return expr_new_integer(1);        /* interior */
        if (dlo == -1 || dhi == -1) return expr_new_integer(0);      /* outside  */
        /* on an endpoint: half the impulse (Mathematica's convention) */
        return mk_fn2("Rational", expr_new_integer(1), expr_new_integer(2));
    }
    /* Symbolic upper limit (the convolution case): the impulse is captured for
     * hi >= x0, provided it is at or after the lower limit (dlo in {0,1}).  Use
     * the full step even at x0 == lo -- the causal Green's-function convention. */
    if (dlo == 0 || dlo == 1)
        return mk_fn1("HeavisideTheta", mk_fn2("Plus", cp(hi),
                          mk_fn2("Times", expr_new_integer(-1), cp(x0))));
    return NULL;   /* cannot decide -- leave the integral unevaluated */
}

/* Find the DiracDelta[...] nodes of e whose argument depends on x (products
 * may nest, e.g. the -1 of a negated term wraps an inner Times); returns the
 * first and counts all. A delta free of x (DiracDelta[y] under an x integral)
 * is a constant factor, not an impulse in x. */
static const Expr* find_dirac(const Expr* e, const char* xn, int* count) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* found = NULL;
    if (head_name_is(e, "DiracDelta") && !free_of(e, xn)) { (*count)++; found = e; }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* f2 = find_dirac(e->data.function.args[i], xn, count);
        if (!found) found = f2;
    }
    return found;
}

/* Does e contain an x-dependent DiracDelta? */
static bool contains_dirac_in(const Expr* e, const char* xn) {
    int c = 0;
    return find_dirac(e, xn, &c) != NULL;
}

/* DiracDelta[u1, ..., un] = DiracDelta[u1] ... DiracDelta[un]: a fresh copy of
 * e with every multivariate delta split into one-argument factors, so that
 * each variable's integral sifts its own factor. */
static Expr* split_multivariate(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return cp(e);
    size_t n = e->data.function.arg_count;
    Expr** as = malloc((n ? n : 1) * sizeof(Expr*));
    if (head_name_is(e, "DiracDelta") && n >= 2) {
        for (size_t i = 0; i < n; i++)
            as[i] = mk_fn1("DiracDelta", cp(e->data.function.args[i]));
        Expr* r = mk_fn("Times", as, n);
        free(as);
        return r;
    }
    for (size_t i = 0; i < n; i++) as[i] = split_multivariate(e->data.function.args[i]);
    Expr* r = expr_new_function(cp(e->data.function.head), as, n);
    free(as);
    return r;
}

/* Does e contain a function headed by the symbol `name`? */
static bool contains_head_named(const Expr* e, const char* name) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, name)) return true;
    if (contains_head_named(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head_named(e->data.function.args[i], name)) return true;
    return false;
}

/* The real roots of g in the closed interval [lo, hi] as a fresh List, or NULL
 * when they cannot be listed as explicit points (an undecided system, a
 * parametric family). Reduce bounds a transcendental g by the interval
 * (Sin[x] == 0 on [-1, 1] is just x == 0) where Solve alone would answer with
 * the whole periodic family. */
static Expr* roots_in(const Expr* g, const Expr* x, const Expr* lo, const Expr* hi) {
    Expr* ineq[5] = { cp(lo), mk_sym("LessEqual"), cp(x), mk_sym("LessEqual"), cp(hi) };
    Expr* sys = mk_fn2("And", mk_fn2("Equal", cp(g), expr_new_integer(0)),
                       mk_fn("Inequality", ineq, 5));
    Expr* rargs[3] = { sys, cp(x), mk_sym("Reals") };
    Expr* red = eval_take(mk_fn("Reduce", rargs, 3));
    Expr* out = NULL;
    if (red->type == EXPR_SYMBOL && strcmp(red->data.symbol.name, "False") == 0) {
        out = mk_fn("List", NULL, 0);                       /* no root */
    } else {
        bool is_or = head_name_is(red, "Or");
        size_t n = is_or ? red->data.function.arg_count : 1;
        Expr** rs = malloc(n * sizeof(Expr*));
        size_t got = 0;
        for (size_t i = 0; i < n; i++) {
            const Expr* c = is_or ? red->data.function.args[i] : red;
            if (!(head_name_is(c, "Equal") && c->data.function.arg_count == 2 &&
                  expr_eq(c->data.function.args[0], x) &&
                  free_of(c->data.function.args[1], x->data.symbol.name) &&
                  !contains_head_named(c->data.function.args[1], "C")))
                break;
            rs[got++] = cp(c->data.function.args[1]);
        }
        if (got == n) out = mk_fn("List", rs, n);
        else for (size_t j = 0; j < got; j++) expr_free(rs[j]);
        free(rs);
    }
    expr_free(red);
    return out;
}

/* Sift a single additive term that carries exactly one x-dependent
 * DiracDelta[g] factor (x-free deltas are constants and ride along in h). A
 * linear g = alpha x + beta gives h(x0)/|alpha| B(x0); a nonlinear g the
 * composition rule
 *     Integrate[h(x) DiracDelta[g(x)], {x, lo, hi}] = Sum_r h(r)/|g'(r)| B(r)
 * over the roots r of g in [lo, hi], each of which must be simple (the sign of
 * g'(r) decided nonzero) -- at a multiple root DiracDelta[g] is not defined.
 * Returns the fresh sifted value, or NULL to decline the whole recognizer. */
static Expr* sift_term(const Expr* term, const Expr* x, const Expr* lo, const Expr* hi) {
    if (x->type != EXPR_SYMBOL) return NULL;
    const char* xn = x->data.symbol.name;
    /* Locate the (unique) x-dependent DiracDelta anywhere in the term;
     * ReplaceAll below removes it regardless of how deeply the product nests. */
    int dcount = 0;
    const Expr* delta = find_dirac(term, xn, &dcount);
    if (!delta || dcount != 1 || delta->data.function.arg_count != 1) return NULL;
    const Expr* arg = delta->data.function.args[0];

    /* alpha = d(arg)/dx: free of x for a linear argument. */
    Expr* alpha = eval_take(mk_fn2("D", cp(arg), cp(x)));
    bool linear = free_of(alpha, xn);
    if (linear) {
        /* Guard against a degenerate argument: arg must equal alpha*x + beta. */
        Expr* recon = eval_take(mk_fn2("Plus",
                          mk_fn2("Times", cp(alpha), cp(x)),
                          subst(arg, x, expr_new_integer(0))));
        Expr* diff = eval_take(mk_fn1("Simplify",
                          mk_fn2("Plus", cp(arg),
                                 mk_fn2("Times", expr_new_integer(-1), recon))));
        linear = (diff->type == EXPR_INTEGER && diff->data.integer == 0);
        expr_free(diff);
        if (!linear) { expr_free(alpha); return NULL; }
    }

    /* h(x) = term with the DiracDelta factor removed. */
    Expr* h = eval_take(mk_fn2("ReplaceAll", cp(term),
                    mk_fn2("Rule", cp(delta), expr_new_integer(1))));

    if (linear) {
        /* beta = arg |_{x=0};  x0 = -beta/alpha. */
        Expr* beta = subst(arg, x, expr_new_integer(0));
        Expr* x0 = eval_take(mk_fn2("Times", expr_new_integer(-1),
                        mk_fn2("Times", beta, mk_fn2("Power", cp(alpha),
                                                     expr_new_integer(-1)))));
        Expr* B = boundary_factor(x0, lo, hi);
        if (!B) { expr_free(alpha); expr_free(x0); expr_free(h); return NULL; }
        Expr* h_at = subst(h, x, cp(x0));
        expr_free(h);
        Expr* invabs = eval_take(mk_fn2("Power", mk_fn1("Abs", cp(alpha)),
                                        expr_new_integer(-1)));
        expr_free(alpha); expr_free(x0);
        return eval_take(mk_fn2("Times", h_at, mk_fn2("Times", invabs, B)));
    }

    /* Nonlinear g: the composition rule over the roots in [lo, hi]. */
    Expr* roots = roots_in(arg, x, lo, hi);
    if (!roots) { expr_free(alpha); expr_free(h); return NULL; }
    size_t nr = roots->data.function.arg_count;
    Expr** parts = malloc((nr ? nr : 1) * sizeof(Expr*));
    size_t np = 0;
    bool ok = true;
    for (size_t i = 0; i < nr; i++) {
        const Expr* r = roots->data.function.args[i];
        Expr* d_at = subst(alpha, x, cp(r));                /* g'(r) */
        int sg = classify_sign(d_at);
        Expr* B = (sg == 1 || sg == -1) ? boundary_factor(r, lo, hi) : NULL;
        if (!B) { expr_free(d_at); ok = false; break; }   /* multiple root / undecided */
        Expr* h_at = subst(h, x, cp(r));
        Expr* inv = mk_fn2("Power", mk_fn1("Abs", d_at), expr_new_integer(-1));
        parts[np++] = mk_fn2("Times", h_at, mk_fn2("Times", inv, B));
    }
    expr_free(alpha); expr_free(h); expr_free(roots);
    if (!ok) {
        for (size_t i = 0; i < np; i++) expr_free(parts[i]);
        free(parts);
        return NULL;
    }
    Expr* val = (np == 0) ? expr_new_integer(0) : eval_take(mk_fn("Plus", parts, np));
    free(parts);
    return val;
}

Expr* integrate_dirac_try(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!contains_dirac(f)) return NULL;
    if (!x || x->type != EXPR_SYMBOL) return NULL;
    const char* xn = x->data.symbol.name;
    /* Only an impulse IN x is sifted here; an x-free delta is a constant
     * factor for the ordinary cascade. */
    Expr* fs = split_multivariate(f);
    if (!contains_dirac_in(fs, xn)) { expr_free(fs); return NULL; }

    /* Expand so a mixed integrand (K + K*delta, 1 + delta, ...) becomes a Plus of
     * additive terms we can partition into delta-bearing and delta-free. */
    Expr* fe = eval_take(mk_fn1("ExpandAll", fs));

    Expr** terms; size_t nt; Expr* one[1];
    if (head_name_is(fe, "Plus")) {
        terms = fe->data.function.args; nt = fe->data.function.arg_count;
    } else {
        one[0] = fe; terms = one; nt = 1;
    }

    Expr** sifted = malloc(nt * sizeof(Expr*)); size_t nsift = 0;
    Expr** remain = malloc(nt * sizeof(Expr*)); size_t nrem = 0;
    bool fail = false;
    for (size_t i = 0; i < nt && !fail; i++) {
        if (!contains_dirac_in(terms[i], xn)) { remain[nrem++] = cp(terms[i]); continue; }
        Expr* v = sift_term(terms[i], x, a, b);
        if (!v) { fail = true; break; }
        sifted[nsift++] = v;
    }

    if (fail || nsift == 0) {
        for (size_t i = 0; i < nsift; i++) expr_free(sifted[i]);
        for (size_t i = 0; i < nrem;  i++) expr_free(remain[i]);
        free(sifted); free(remain); expr_free(fe);
        return NULL;   /* delta-free after expansion, or undecidable */
    }

    /* Assemble: sum of sifted values + Integrate[remainder, {x, a, b}]. */
    Expr** parts = malloc((nsift + 1) * sizeof(Expr*)); size_t np = 0;
    for (size_t i = 0; i < nsift; i++) parts[np++] = sifted[i];
    if (nrem > 0) {
        Expr* rem_sum = (nrem == 1) ? remain[0]
                        : eval_take(expr_new_function(expr_new_symbol("Plus"),
                                                      remain, nrem));
        Expr* listargs[3] = { cp(x), cp(a), cp(b) };
        Expr* rng = mk_fn("List", listargs, 3);
        parts[np++] = eval_take(mk_fn2("Integrate", rem_sum, rng));
    }
    free(sifted); free(remain);

    Expr* result = (np == 1) ? parts[0]
                   : eval_take(expr_new_function(expr_new_symbol("Plus"), parts, np));
    free(parts);
    expr_free(fe);
    return result;
}

bool integrate_dirac_in(Expr* f, Expr* x) {
    if (!x || x->type != EXPR_SYMBOL || !contains_dirac(f)) return false;
    return contains_dirac_in(f, x->data.symbol.name);
}
