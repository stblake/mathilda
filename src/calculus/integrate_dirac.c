/* integrate_dirac.c -- DiracDelta sifting under a definite integral.
 * See integrate_dirac.h for the contract and conventions. */

#include "integrate_dirac.h"
#include "expr.h"
#include "eval.h"

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>

/* ---- small construction / evaluation helpers (mirrors integrate_symmetry.c) - */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
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
    bool iszero = (z->type == EXPR_SYMBOL && strcmp(z->data.symbol.name, "True") == 0);
    expr_free(z);
    if (iszero) return 0;
    Expr* p = eval_take(mk_fn1("Positive", cp(diff)));
    bool pos = (p->type == EXPR_SYMBOL && strcmp(p->data.symbol.name, "True") == 0);
    expr_free(p);
    if (pos) return 1;
    Expr* ng = eval_take(mk_fn1("Negative", cp(diff)));
    bool neg = (ng->type == EXPR_SYMBOL && strcmp(ng->data.symbol.name, "True") == 0);
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

/* Find a DiracDelta[...] node anywhere in e (products may nest, e.g. the -1 of a
 * negated term wraps an inner Times); returns the first and counts all. */
static const Expr* find_dirac(const Expr* e, int* count) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* found = NULL;
    if (head_name_is(e, "DiracDelta")) { (*count)++; found = e; }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* f2 = find_dirac(e->data.function.args[i], count);
        if (!found) found = f2;
    }
    return found;
}

/* Sift a single additive term that carries exactly one DiracDelta[linear] factor.
 * Returns the fresh sifted value, or NULL to decline the whole recognizer. */
static Expr* sift_term(const Expr* term, const Expr* x, const Expr* lo, const Expr* hi) {
    /* Locate the (unique) DiracDelta anywhere in the term; ReplaceAll below
     * removes it regardless of how deeply the product nests. */
    int dcount = 0;
    const Expr* delta = find_dirac(term, &dcount);
    if (!delta || dcount != 1 || delta->data.function.arg_count != 1) return NULL;

    const Expr* arg = delta->data.function.args[0];
    if (x->type != EXPR_SYMBOL) return NULL;
    const char* xn = x->data.symbol.name;

    /* Linear argument a*x + b: alpha = d(arg)/dx must be free of x and nonzero. */
    Expr* alpha = eval_take(mk_fn2("D", cp(arg), cp(x)));
    if (!free_of(alpha, xn)) { expr_free(alpha); return NULL; }
    Expr* alpha_at0 = alpha;                     /* alpha is constant in x */
    /* beta = arg |_{x=0};  x0 = -beta/alpha. */
    Expr* beta = subst(arg, x, expr_new_integer(0));
    Expr* x0 = eval_take(mk_fn2("Times", expr_new_integer(-1),
                    mk_fn2("Times", beta, mk_fn2("Power", cp(alpha_at0),
                                                 expr_new_integer(-1)))));
    /* Guard against a nonlinear/degenerate argument: arg must equal alpha*x + beta. */
    {
        Expr* recon = eval_take(mk_fn2("Plus",
                          mk_fn2("Times", cp(alpha_at0), cp(x)),
                          subst(arg, x, expr_new_integer(0))));
        Expr* diff = eval_take(mk_fn1("Simplify",
                          mk_fn2("Plus", cp(arg),
                                 mk_fn2("Times", expr_new_integer(-1), recon))));
        bool linear = (diff->type == EXPR_INTEGER && diff->data.integer == 0);
        expr_free(diff);
        if (!linear) { expr_free(alpha); expr_free(x0); return NULL; }
    }

    Expr* B = boundary_factor(x0, lo, hi);
    if (!B) { expr_free(alpha); expr_free(x0); return NULL; }

    /* h(x) = term with the DiracDelta factor removed; evaluate h(x0). */
    Expr* h = eval_take(mk_fn2("ReplaceAll", cp(term),
                    mk_fn2("Rule", cp(delta), expr_new_integer(1))));
    Expr* h_at = subst(h, x, cp(x0));
    expr_free(h);

    Expr* invabs = eval_take(mk_fn2("Power", mk_fn1("Abs", cp(alpha_at0)),
                                    expr_new_integer(-1)));
    expr_free(alpha); expr_free(x0);

    Expr* val = eval_take(mk_fn2("Times", h_at, mk_fn2("Times", invabs, B)));
    return val;
}

Expr* integrate_dirac_try(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!contains_dirac(f)) return NULL;
    if (!x || x->type != EXPR_SYMBOL) return NULL;

    /* Expand so a mixed integrand (K + K*delta, 1 + delta, ...) becomes a Plus of
     * additive terms we can partition into delta-bearing and delta-free. */
    Expr* fe = eval_take(mk_fn1("ExpandAll", cp(f)));

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
        if (!contains_dirac(terms[i])) { remain[nrem++] = cp(terms[i]); continue; }
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
