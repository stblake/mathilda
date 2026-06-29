/* integrate_chebychev.c
 *
 * Integration of Chebychev binomial differentials  x^p (a x^r + b)^q  (p, q, r
 * rational, a, b free of x).  See integrate_chebychev.h for the high-level
 * description and Chebychev's theorem.
 *
 * ---------------------------------------------------------------------
 * Algorithm
 * ---------------------------------------------------------------------
 * 1. Recognition (single structural scan, no expansion): split f into the
 *    constant prefactor C (factors free of x), the pure power x^p, and the
 *    binomial power (a x^r + b)^q.  Extract p, q, r as reduced rationals and
 *    a, b as expressions free of x.  Reject anything that is not exactly of
 *    this shape (trinomial base, two binomials, x inside a transcendental,
 *    non-rational exponents, ...).
 *
 * 2. Type selection (Chebychev):
 *      q in Z, N = LCM(den p, den r) >= 2   -> Type I   : x = u^N
 *      else (s = den q >= 2):
 *        (p+1)/r in Z                       -> Type II  : u^s = a x^r + b
 *        q + (p+1)/r in Z                   -> Type III : u = x^r, t^s=(a u+b)/u
 *        else  non-elementary               -> NULL (cascade falls through)
 *
 * 3. Build the rationalised integrand g in the new variable (closed forms
 *    derived in the per-type comments below), confirm it is a rational function
 *    of that variable (Numerator/Denominator are polynomials), and recurse
 *    Integrate[g, newvar].
 *
 * 4. Back-substitute.  Types I and II introduce a single principal radical and
 *    the substitution is an exact bijection on the principal-branch domain, so
 *    the antiderivative is correct by construction (mirrors LinearRadicals /
 *    LinearRatioRadicals -- no differentiate-back verification, which is
 *    branch-fragile and prohibitively expensive with symbolic a, b).  Type III
 *    is branch-sensitive: the forward substitution uses the *combined* radical
 *    t = ((a x^r + b)/x^r)^(1/s) but the antiderivative must use the *split*
 *    radical  t -> (a x^r + b)^(1/s) / x^(r/s)  so every radical in the result
 *    is the principal root of a sub-expression that appears in the original
 *    integrand.  (This is the generalisation of Sqrt[u-1]/Sqrt[u], not
 *    Sqrt[(u-1)/u], from the worked example.)  Finally multiply C back in.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take`/`replace_one`
 * consume their argument.  We never expr_free(res).  The small helper set
 * mirrors integrate_linrad.c / integrate_linratiorad.c (each integration-method
 * file keeps a private copy).
 */

#include "integrate_chebychev.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "poly.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirrors integrate_linrad.c)       */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(name), args, 2);
}

static Expr* mk_fn3(const char* name, Expr* a, Expr* b, Expr* c) {
    Expr* args[3] = { a, b, c };
    return expr_new_function(expr_new_symbol(name), args, 3);
}

/* Power[base, rn/rd], with the exponent reduced/normalised by make_rational
 * (returns an Integer when rd | rn).  Consumes base. */
static Expr* mk_pow_rat(Expr* base, int64_t rn, int64_t rd) {
    return mk_fn2("Power", base, make_rational(rn, rd));
}

/* Power[base, n] for an integer exponent.  Consumes base. */
static Expr* mk_pow_int(Expr* base, int64_t n) {
    return mk_fn2("Power", base, mk_int(n));
}

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True if `f` contains no subexpression structurally equal to `x`. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_unintegrated(e->data.function.args[i])) return true;
    }
    return false;
}

/* Cancel[Together[e]], consuming `e`. */
static Expr* cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* ReplaceAll[expr, from -> to]; consumes `expr`; returns owned (evaluated). */
static Expr* replace_one(Expr* expr, const Expr* from, const Expr* to) {
    Expr* rule = mk_fn2("Rule", expr_copy((Expr*)from), expr_copy((Expr*)to));
    Expr* call = internal_replace_all((Expr*[]){ expr, rule }, 2);
    return eval_take(call);
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2("Integrate", expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* True iff `e`, after Cancel[Together[.]], is a rational function of the single
 * variable `var` (Numerator and Denominator are polynomials in var).  Borrows
 * e and var. */
static bool is_rational_in(Expr* e, Expr* var) {
    Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(e)));
    Expr* den = eval_take(mk_fn1("Denominator", expr_copy(e)));
    Expr* vars[1] = { var };
    bool ok = num && den &&
              is_polynomial(num, vars, 1) &&
              is_polynomial(den, vars, 1);
    if (num) expr_free(num);
    if (den) expr_free(den);
    return ok;
}

/* ---------------------------------------------------------------------- */
/* Small int64 rational arithmetic (exponent bookkeeping)                 */
/* ---------------------------------------------------------------------- */

/* Reduce n/d to lowest terms with d > 0. */
static void frac_reduce(int64_t* n, int64_t* d) {
    if (*d < 0) { *n = -*n; *d = -*d; }
    int64_t g = gcd(*n < 0 ? -*n : *n, *d);
    if (g > 1) { *n /= g; *d /= g; }
}

/* (an/ad) + (bn/bd) -> on/od (reduced). */
static void frac_add(int64_t an, int64_t ad, int64_t bn, int64_t bd,
                     int64_t* on, int64_t* od) {
    *on = an * bd + bn * ad;
    *od = ad * bd;
    frac_reduce(on, od);
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define CH_MAX_DEPTH 8
static int ch_depth = 0;
static unsigned long ch_sym_counter = 0;

/* ---------------------------------------------------------------------- */
/* Recognition                                                            */
/* ---------------------------------------------------------------------- */

/* Is `g` the symbol x or Power[x, rational]?  On success fills the exponent
 * en/ed (reduced, ed > 0) and returns true.  Borrows g, x. */
static bool is_x_power(Expr* g, Expr* x, int64_t* en, int64_t* ed) {
    if (expr_eq(g, x)) { *en = 1; *ed = 1; return true; }
    if (head_is(g, SYM_Power) && g->data.function.arg_count == 2 &&
        expr_eq(g->data.function.args[0], x) &&
        is_rational(g->data.function.args[1], en, ed)) {
        frac_reduce(en, ed);
        return true;
    }
    return false;
}

/* Parse a monomial T = c * x^r into (c, r).  c is any product of factors free
 * of x (returned as an owned Expr; integer 1 when none); exactly one factor
 * must be x or Power[x, rational].  Returns false (allocating nothing) when T
 * is not such a monomial.  Borrows T, x. */
static bool parse_monomial(Expr* T, Expr* x, Expr** c_out,
                           int64_t* rn, int64_t* rd) {
    int64_t en, ed;
    if (is_x_power(T, x, &en, &ed)) {
        *c_out = mk_int(1);
        *rn = en; *rd = ed;
        return true;
    }
    if (head_is(T, SYM_Times)) {
        Expr* c = mk_int(1);
        int xcount = 0;
        int64_t lrn = 0, lrd = 1;
        for (size_t i = 0; i < T->data.function.arg_count; i++) {
            Expr* h = T->data.function.args[i];
            int64_t hn, hd;
            if (is_x_power(h, x, &hn, &hd)) {
                xcount++;
                lrn = hn; lrd = hd;
            } else if (expr_free_of(h, x)) {
                c = mk_fn2("Times", c, expr_copy(h));
            } else {
                expr_free(c);
                return false;   /* x buried in a non-monomial factor */
            }
        }
        if (xcount != 1) { expr_free(c); return false; }
        *c_out = c;
        *rn = lrn; *rd = lrd;
        return true;
    }
    return false;
}

/* Parse the binomial base = a x^r + b: a 2-term Plus with exactly one term free
 * of x (= b) and the other a monomial a x^r.  Fills owned a, b and r = rn/rd.
 * Returns false (allocating nothing) otherwise.  Borrows base, x. */
static bool parse_binom(Expr* base, Expr* x, Expr** a_out, Expr** b_out,
                        int64_t* rn, int64_t* rd) {
    if (!head_is(base, SYM_Plus) || base->data.function.arg_count != 2)
        return false;
    Expr* t0 = base->data.function.args[0];
    Expr* t1 = base->data.function.args[1];
    Expr* xterm; Expr* bterm;
    if (expr_free_of(t0, x) && !expr_free_of(t1, x))      { bterm = t0; xterm = t1; }
    else if (expr_free_of(t1, x) && !expr_free_of(t0, x)) { bterm = t1; xterm = t0; }
    else return false;   /* both free (no x) or both x-dependent (trinomial-ish) */

    Expr* a = NULL;
    if (!parse_monomial(xterm, x, &a, rn, rd)) return false;
    if (*rn == 0) { expr_free(a); return false; }   /* r == 0 */

    *a_out = a;
    *b_out = expr_copy(bterm);
    return true;
}

/* Recognise f = C * x^p * (a x^r + b)^q.  On success returns true and fills:
 *   *C  owned constant prefactor (Integer 1 when none),
 *   p = pn/pd, q = qn/qd, r = rn/rd  (reduced, denominators > 0),
 *   *a, *b owned, free of x.
 * On false nothing is allocated.  Borrows f, x. */
static bool recognise(Expr* f, Expr* x,
                      Expr** C, int64_t* pn, int64_t* pd,
                      int64_t* qn, int64_t* qd, int64_t* rn, int64_t* rd,
                      Expr** a, Expr** b) {
    /* Factor list: args of Times, or the single node f. */
    Expr** factors; size_t nfac;
    if (head_is(f, SYM_Times)) {
        factors = f->data.function.args;
        nfac    = f->data.function.arg_count;
    } else {
        factors = &f;
        nfac    = 1;
    }

    Expr* c = mk_int(1);          /* constant prefactor */
    int64_t ppn = 0, ppd = 1;     /* accumulated p */
    bool have_binom = false;
    Expr* la = NULL; Expr* lb = NULL;
    int64_t lqn = 0, lqd = 1, lrn = 0, lrd = 1;

    for (size_t i = 0; i < nfac; i++) {
        Expr* g = factors[i];

        if (expr_free_of(g, x)) {                 /* constant factor */
            c = mk_fn2("Times", c, expr_copy(g));
            continue;
        }

        int64_t en, ed;
        if (is_x_power(g, x, &en, &ed)) {         /* pure power x^p */
            int64_t nn, dd;
            frac_add(ppn, ppd, en, ed, &nn, &dd);
            ppn = nn; ppd = dd;
            continue;
        }

        if (head_is(g, SYM_Power) && g->data.function.arg_count == 2) {
            Expr* gb = g->data.function.args[0];
            Expr* ge = g->data.function.args[1];
            int64_t tqn, tqd;
            if (!expr_free_of(gb, x) && is_rational(ge, &tqn, &tqd)) {
                /* candidate binomial power */
                if (have_binom) goto fail;        /* a second binomial */
                Expr* pa = NULL; Expr* pb = NULL; int64_t prn, prd;
                if (!parse_binom(gb, x, &pa, &pb, &prn, &prd)) goto fail;
                la = pa; lb = pb;
                lqn = tqn; lqd = tqd; frac_reduce(&lqn, &lqd);
                lrn = prn; lrd = prd; frac_reduce(&lrn, &lrd);
                have_binom = true;
                continue;
            }
            goto fail;   /* x-dependent power that is not a binomial */
        }

        goto fail;       /* any other x-dependent factor */
    }

    if (!have_binom) goto fail;

    *C = c;
    *pn = ppn; *pd = ppd;
    *qn = lqn; *qd = lqd;
    *rn = lrn; *rd = lrd;
    *a = la; *b = lb;
    return true;

fail:
    expr_free(c);
    if (la) expr_free(la);
    if (lb) expr_free(lb);
    return false;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

static Expr* ch_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL)   return NULL;
    if (expr_free_of(f, x))       return NULL;   /* nothing to integrate in x */
    if (ch_depth >= CH_MAX_DEPTH) return NULL;

    Expr* C = NULL; Expr* a = NULL; Expr* b = NULL;
    int64_t pn, pd, qn, qd, rn, rd;
    if (!recognise(f, x, &C, &pn, &pd, &qn, &qd, &rn, &rd, &a, &b))
        return NULL;

    /* Fresh substitution variable. */
    char uname[96];
    snprintf(uname, sizeof(uname),
             "Integrate`ChebychevAlgebraic`u$%lu", ch_sym_counter++);
    Expr* u = expr_new_symbol(uname);

    Expr* result = NULL;
    Expr* g = NULL;   /* rationalised integrand in u */
    Expr* G = NULL;   /* its antiderivative */
    Expr* back = NULL;

    /* --- Type selection. --- */
    int64_t s = qd;   /* den(q) */

    if (qd == 1) {
        /* ---------- Type I : q in Z, x = u^N, N = LCM(den p, den r). -------- */
        int64_t N = lcm(pd, rd);
        if (N < 2) goto done;          /* p, r integral -> rational, not ours */
        int64_t P = N * pn / pd;       /* exact: pd | N */
        int64_t R = N * rn / rd;       /* exact: rd | N */
        int64_t qz = qn;               /* q as an integer */

        /* g = N u^(N-1+P) (a u^R + b)^q */
        Expr* binom = mk_fn2("Plus",
            mk_fn2("Times", expr_copy(a), mk_pow_int(expr_copy(u), R)),
            expr_copy(b));
        g = mk_fn3("Times",
            mk_int(N),
            mk_pow_int(expr_copy(u), N - 1 + P),
            mk_pow_int(binom, qz));
        g = cancel_together(g);
        if (!g || !is_rational_in(g, u)) goto done;

        ch_depth++;
        G = integrate_in(g, u);
        ch_depth--;
        if (!G) goto done;

        back = mk_pow_rat(expr_copy(x), 1, N);     /* u -> x^(1/N) */
        result = replace_one(G, u, back);          /* consumes G */
        G = NULL;
    } else if ( ((pn + pd) * rd) % (pd * rn) == 0 ) {
        /* ---------- Type II : (p+1)/r in Z, u^s = a x^r + b. --------------- */
        /* k = (p+1)/r = ((pn+pd) rd) / (pd rn) (an integer here). */
        int64_t k = ((pn + pd) * rd) / (pd * rn);

        /* g = (s/(a r)) u^(s q + s - 1) ((u^s - b)/a)^(k-1).
         *   s q = qn (since q = qn/s);  1/r = rd/rn.
         *   coeff = s * (rd/rn) * a^-1. */
        Expr* coeff = mk_fn3("Times",
            mk_int(s),
            make_rational(rd, rn),
            mk_pow_int(expr_copy(a), -1));
        Expr* usb_over_a = mk_fn2("Times",
            mk_fn2("Plus", mk_pow_int(expr_copy(u), s),
                           mk_fn2("Times", mk_int(-1), expr_copy(b))),
            mk_pow_int(expr_copy(a), -1));          /* (u^s - b)/a */
        g = mk_fn3("Times",
            coeff,
            mk_pow_int(expr_copy(u), qn + s - 1),
            mk_pow_int(usb_over_a, k - 1));
        g = cancel_together(g);
        if (!g || !is_rational_in(g, u)) goto done;

        ch_depth++;
        G = integrate_in(g, u);
        ch_depth--;
        if (!G) goto done;

        /* u -> (a x^r + b)^(1/s). */
        back = mk_pow_rat(
            mk_fn2("Plus", mk_fn2("Times", expr_copy(a),
                                  mk_pow_rat(expr_copy(x), rn, rd)),
                           expr_copy(b)), 1, s);
        result = replace_one(G, u, back);           /* consumes G */
        G = NULL;
    } else {
        /* ---------- Type III : q + (p+1)/r in Z. -------------------------- */
        /* m = q + (p+1)/r = qn/qd + ((pn+pd) rd)/(pd rn). */
        int64_t mn, md;
        frac_add(qn, qd, (pn + pd) * rd, pd * rn, &mn, &md);
        if (md != 1) goto done;        /* non-elementary */
        int64_t m = mn;

        /* g (in t) = -(s/r) b^m t^(s q + s - 1) (t^s - a)^(-(m+1)).
         *   -(s/r) = -s * (rd/rn) = Times[-1, s, rd/rn]. */
        Expr* tsa = mk_fn2("Plus", mk_pow_int(expr_copy(u), s),
                                   mk_fn2("Times", mk_int(-1), expr_copy(a)));
        g = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
            mk_int(-1),
            mk_int(s),
            make_rational(rd, rn),                 /* 1/r */
            mk_pow_int(expr_copy(b), m),           /* b^m */
            mk_pow_int(expr_copy(u), qn + s - 1),  /* t^(s q + s - 1) */
            mk_pow_int(tsa, -(m + 1))              /* (t^s - a)^-(m+1) */
        }, 6);
        g = cancel_together(g);
        if (!g || !is_rational_in(g, u)) goto done;

        ch_depth++;
        G = integrate_in(g, u);
        ch_depth--;
        if (!G) goto done;

        /* Branch-correct split back-substitution:
         *   t -> (a x^r + b)^(1/s) * x^(-r/s). */
        back = mk_fn2("Times",
            mk_pow_rat(
                mk_fn2("Plus", mk_fn2("Times", expr_copy(a),
                                      mk_pow_rat(expr_copy(x), rn, rd)),
                               expr_copy(b)), 1, s),
            mk_pow_rat(expr_copy(x), -rn, rd * s));
        result = replace_one(G, u, back);           /* consumes G */
        G = NULL;
    }

    /* Multiply the constant prefactor back in. */
    if (result) {
        result = eval_take(mk_fn2("Times", expr_copy(C), result));
    }

done:
    if (C)    expr_free(C);
    if (a)    expr_free(a);
    if (b)    expr_free(b);
    if (u)    expr_free(u);
    if (g)    expr_free(g);
    if (G)    expr_free(G);
    if (back) expr_free(back);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_chebychev_try(Expr* f, Expr* x) {
    return ch_core(f, x);
}

Expr* builtin_integrate_chebychev(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return ch_core(f, x);
}

void integrate_chebychev_init(void) {
    symtab_add_builtin("Integrate`ChebychevAlgebraic", builtin_integrate_chebychev);
    symtab_get_def("Integrate`ChebychevAlgebraic")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`ChebychevAlgebraic",
        "Integrate`ChebychevAlgebraic[f, x] integrates a Chebychev binomial\n"
        "differential x^p (a x^r + b)^q (p, q, r rational, a, b free of x). By\n"
        "Chebychev's theorem this is elementary iff one of q, (p+1)/r, q+(p+1)/r is\n"
        "an integer, giving substitutions x = u^N (Type I), u^s = a x^r + b\n"
        "(Type II), or u = x^r then t^s = (a u + b)/u (Type III) that reduce f to a\n"
        "rational function. Strict: returns unevaluated when f is not a Chebychev\n"
        "binomial or the integral is non-elementary.");
}
