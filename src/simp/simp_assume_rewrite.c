#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* Assumption-driven seed rewriters                                        */
/* ----------------------------------------------------------------------- */

/* For each direct EXPR_SYMBOL fact-target, generate context-specific
 * rewrite rules and apply them via ReplaceRepeated. The rules are
 * unconditional in pattern form: their conditional nature is captured by
 * the choice of the rule's free symbol -- e.g., we only emit
 *   Power[Power[<x>, 2], Rational[1, 2]] :> <x>
 * when <x> is the literal symbol that the assumption set says is
 * positive. So the rules are valid by construction whenever applied.
 *
 * The generated rule list is built as a string and parsed; this is
 * cheaper to maintain than constructing the AST by hand and matches the
 * style used in trigsimp.c.
 */

static bool sym_already_listed(char** list, size_t n, const char* s) {
    for (size_t i = 0; i < n; i++) if (strcmp(list[i], s) == 0) return true;
    return false;
}

/* Walk the assumption fact list and collect every EXPR_SYMBOL that the
 * context proves positive, real, integer, or even. The caller passes
 * pre-sized arrays plus the maximum count. */
static void collect_known_symbols(const AssumeCtx* ctx,
                                  char** positives, size_t* npos,
                                  char** reals,     size_t* nreal,
                                  char** integers,  size_t* nint,
                                  char** negatives, size_t* nneg,
                                  char** evens,     size_t* neven,
                                  size_t cap) {
    *npos = *nreal = *nint = *nneg = *neven = 0;
    if (!ctx) return;
    /* Mod[m, 2] == 0 hides `m` inside a function argument, so a top-level
     * "scan operands of facts" alone misses it. We additionally walk
     * Mod[s, _] for any symbol s that appears under an even-type fact. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION) continue;
        for (size_t j = 0; j < f->data.function.arg_count; j++) {
            Expr* a = f->data.function.args[j];
            if (a->type == EXPR_SYMBOL) {
                const char* nm = a->data.symbol.name;
                if (assume_known_positive(ctx, a) && *npos < cap && !sym_already_listed(positives, *npos, nm)) {
                    positives[(*npos)++] = (char*)nm;
                }
                if (assume_known_negative(ctx, a) && *nneg < cap && !sym_already_listed(negatives, *nneg, nm)) {
                    negatives[(*nneg)++] = (char*)nm;
                }
                if (assume_known_real(ctx, a) && *nreal < cap && !sym_already_listed(reals, *nreal, nm)) {
                    reals[(*nreal)++] = (char*)nm;
                }
                if (assume_known_integer(ctx, a) && *nint < cap && !sym_already_listed(integers, *nint, nm)) {
                    integers[(*nint)++] = (char*)nm;
                }
                if (assume_known_even(ctx, a) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            } else if (a->type == EXPR_FUNCTION &&
                       a->data.function.head &&
                       a->data.function.head->type == EXPR_SYMBOL &&
                       a->data.function.head->data.symbol.name == SYM_Mod &&
                       a->data.function.arg_count == 2 &&
                       a->data.function.args[0]->type == EXPR_SYMBOL) {
                Expr* sym = a->data.function.args[0];
                const char* nm = sym->data.symbol.name;
                if (assume_known_even(ctx, sym) && *neven < cap && !sym_already_listed(evens, *neven, nm)) {
                    evens[(*neven)++] = (char*)nm;
                }
            }
        }
    }
}

/* ----------------------------------------------------------------------- */
/* Structural assumption rewrites                                          */
/*                                                                         */
/* Some assumption-driven simplifications cannot be expressed as a rule    */
/* keyed on a pre-enumerated symbol: Floor[e]->e needs to know e is an     */
/* integer for an ARBITRARY subexpression e; Ceiling[x]->3 needs interval  */
/* reasoning; Mod[a,4]->1 needs modular reasoning from an Element fact;    */
/* Re[a+b I]->a needs every free symbol to be real. These are handled by a */
/* bottom-up structural walk that queries the assume_known_* API at each   */
/* node. The walk runs inside apply_assumption_rules, so Simplify (which   */
/* calls apply_assumption_rules) benefits from it too.                     */
/* ----------------------------------------------------------------------- */

static bool sr_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Numeric value of an Integer/Real/Bigint/Rational literal. */
static bool sr_num(const Expr* e, double* v) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *v = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *v = e->data.real; return true; }
    if (e->type == EXPR_BIGINT)  { *v = mpz_get_d(e->data.bigint); return true; }
    if (sr_head(e, "Rational") && e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER && q->data.integer != 0) {
            *v = (double)p->data.integer / (double)q->data.integer; return true;
        }
    }
    return false;
}

/* True iff e is a purely-imaginary numeric literal Complex[0, c]; sets *im = c. */
static bool sr_pure_imag(const Expr* e, double* im) {
    if (!sr_head(e, "Complex") || e->data.function.arg_count != 2) return false;
    double re;
    if (!sr_num(e->data.function.args[0], &re) || re != 0.0) return false;
    return sr_num(e->data.function.args[1], im);
}

/* Does ctx prove the strict order A < B by a direct inequality fact? */
static bool sr_proves_slt(const AssumeCtx* ctx, const Expr* A, const Expr* B) {
    if (!ctx) return false;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* x = f->data.function.args[0];
        const Expr* y = f->data.function.args[1];
        if (h == SYM_Less && expr_eq((Expr*)x, (Expr*)A) && expr_eq((Expr*)y, (Expr*)B))
            return true;
        if (h == SYM_Greater && expr_eq((Expr*)x, (Expr*)B) && expr_eq((Expr*)y, (Expr*)A))
            return true;
    }
    return false;
}

/* Does ctx prove the non-strict order A <= B by a direct inequality fact? A
 * strict A < B also establishes A <= B, so this accepts both. Used to gate the
 * closed-interval inverse-function collapses (e.g. ArcSin[Sin[x]] on
 * [-Pi/2, Pi/2]), where the identity holds at the endpoints too. */
static bool sr_proves_sle(const AssumeCtx* ctx, const Expr* A, const Expr* B) {
    if (!ctx) return false;
    if (sr_proves_slt(ctx, A, B)) return true;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* x = f->data.function.args[0];
        const Expr* y = f->data.function.args[1];
        if (h == SYM_LessEqual && expr_eq((Expr*)x, (Expr*)A) && expr_eq((Expr*)y, (Expr*)B))
            return true;
        if (h == SYM_GreaterEqual && expr_eq((Expr*)x, (Expr*)B) && expr_eq((Expr*)y, (Expr*)A))
            return true;
    }
    return false;
}

/* Tightest numeric bounds on the bare symbol `sym` from direct inequality
 * facts. Returns whether a lower / upper bound was found; *lo_strict indicates
 * a strict `>`. */
static void sr_bounds(const AssumeCtx* ctx, const Expr* sym,
                      bool* has_lo, double* lo, bool* lo_strict,
                      bool* has_hi, double* hi, bool* hi_strict) {
    *has_lo = *has_hi = false;
    if (!ctx) return;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* A = f->data.function.args[0];
        const Expr* B = f->data.function.args[1];
        double k;
        bool sym_left, strict; int kind; /* kind: -1 lower, +1 upper */
        if (expr_eq((Expr*)A, (Expr*)sym) && sr_num(B, &k)) { sym_left = true; }
        else if (expr_eq((Expr*)B, (Expr*)sym) && sr_num(A, &k)) { sym_left = false; }
        else continue;
        /* Normalise to "sym <op> k". */
        if (h == SYM_Less)        { kind = sym_left ? +1 : -1; strict = true;  }
        else if (h == SYM_LessEqual)   { kind = sym_left ? +1 : -1; strict = false; }
        else if (h == SYM_Greater)     { kind = sym_left ? -1 : +1; strict = true;  }
        else if (h == SYM_GreaterEqual){ kind = sym_left ? -1 : +1; strict = false; }
        else continue;
        if (kind < 0) { /* lower bound sym > k or sym >= k */
            if (!*has_lo || k > *lo || (k == *lo && strict && !*lo_strict)) {
                *has_lo = true; *lo = k; *lo_strict = strict;
            }
        } else {        /* upper bound sym < k or sym <= k */
            if (!*has_hi || k < *hi || (k == *hi && strict && !*hi_strict)) {
                *has_hi = true; *hi = k; *hi_strict = strict;
            }
        }
    }
}

/* Collect distinct bare non-constant symbols appearing in e. */
static void sr_collect_syms(const Expr* e, const Expr*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol.name)) return;
        for (size_t i = 0; i < *n; i++)
            if ((*out)[i]->data.symbol.name == e->data.symbol.name) return;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *out = realloc(*out, *cap * sizeof(Expr*)); }
        (*out)[(*n)++] = e;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        /* Descend into arguments only: a function's head (Plus, Times, Re, ...)
         * is a structural operator, not a variable whose reality we require. */
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            sr_collect_syms(e->data.function.args[i], out, n, cap);
    }
}

/* ----------------------------------------------------------------------- */
/* Symbolic-integer Pi periodicity: Sin[t + k Pi] -> (-1)^k Sin[t], etc.,    */
/* where k is a provable integer-linear combination (2n, 2n+1, ...) of       */
/* integer-assumed symbols -- the cases the coefficient-1 string rules miss. */
/* ----------------------------------------------------------------------- */

/* True iff the symbol Pi appears anywhere in e. */
static bool sr_mentions_pi(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Pi;
    if (e->type == EXPR_FUNCTION) {
        if (sr_mentions_pi(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (sr_mentions_pi(e->data.function.args[i])) return true;
    }
    return false;
}

/* If `term` is (integer)*Pi, return its integer coefficient k as a fresh Expr
 * (caller owns); else NULL. Computed as k = term/Pi (evaluated): robust to any
 * nesting/ordering of the Times. A term linear in Pi leaves a Pi-free quotient;
 * anything else (no Pi, Pi^2, Pi inside a function) keeps Pi and is rejected,
 * as is a non-integer quotient. */
static Expr* pi_integer_coeff(const Expr* term, const AssumeCtx* ctx) {
    if (!sr_mentions_pi(term)) return NULL;
    Expr* piinv = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_new_symbol("Pi"), expr_new_integer(-1) }, 2);
    Expr* k = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_copy((Expr*)term), piinv }, 2));
    if (sr_mentions_pi(k)) { expr_free(k); return NULL; }  /* not linear in Pi */
    if (prov_int(ctx, k)) return k;
    expr_free(k);
    return NULL;
}

/* Parity of a provable integer K: 0 even, 1 odd, -1 undetermined. */
static int pi_parity(const AssumeCtx* ctx, const Expr* K) {
    if (K->type == EXPR_INTEGER) return (int)(((K->data.integer % 2) + 2) % 2);
    if (sr_head(K, "Plus")) {
        int p = 0;
        for (size_t i = 0; i < K->data.function.arg_count; i++) {
            int pt = pi_parity(ctx, K->data.function.args[i]);
            if (pt < 0) return -1;
            p ^= pt;
        }
        return p;
    }
    if (prov_even(ctx, K)) return 0;
    /* K odd  <=>  K-1 even. */
    Expr* km1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy((Expr*)K), expr_new_integer(-1) }, 2));
    int km1_even = prov_even(ctx, km1);
    expr_free(km1);
    return km1_even ? 1 : -1;
}

/* ----------------------------------------------------------------------- */
/* Root-of-unity power canonicalizer under integer assumptions.            */
/*                                                                         */
/* The elegant core: for a root of unity raised to an integer-valued       */
/* exponent, reduce the exponent modulo the base's order. Everything is    */
/* routed to (-1)^Y via the identities  I = (-1)^(1/2),  -I = (-1)^(-1/2), */
/* E^(I Pi X) = (-1)^X, and (-1)^Y is reduced by a Boolean-ring mod-2       */
/* reduction of the integer-valued part of Y:                              */
/*   F2[g]/(g^2 - g):  coeffs mod 2, g^e -> g (sound: g^e == g (mod 2) for  */
/*   integer g), like monomials XOR. So n^2 -> n, n(n+1) -> 0, 3n -> n,     */
/*   2n+1 -> 1, n-m -> n+m.  A rational constant in Y (e.g. the 1/2 of I)   */
/*   stays as a literal factor (-1)^(1/2) = I.                              */
/* Sound by construction and declines (returns NULL) on anything it cannot  */
/* express as a polynomial in provably-integer generators.                 */
/* ----------------------------------------------------------------------- */

/* Complex[0, sign]: the imaginary unit (sign=1) or its negative (sign=-1). */
static bool is_imag_unit(const Expr* e, long sign) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Complex) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* re = e->data.function.args[0];
    const Expr* im = e->data.function.args[1];
    return re->type == EXPR_INTEGER && re->data.integer == 0 &&
           im->type == EXPR_INTEGER && im->data.integer == sign;
}

/* Replace Mod[t, M] by t wherever M is a nonzero even integer literal:
 * Mod[t, M] == t (mod 2) when 2 | M, so this is sound inside a mod-2
 * reduction. Recurses; caller owns the result. */
static Expr* mod2_absorb(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (sr_head(e, "Mod") && e->data.function.arg_count == 2) {
        const Expr* M = e->data.function.args[1];
        if (M->type == EXPR_INTEGER && M->data.integer != 0 &&
            (M->data.integer % 2) == 0)
            return mod2_absorb(e->data.function.args[0]);
    }
    size_t n = e->data.function.arg_count;
    Expr** a = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = mod2_absorb(e->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return out;
}

/* Parity (0/1) of an integer literal/bigint, else -1. */
static int lit_parity(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (int)(((e->data.integer % 2) + 2) % 2);
    if (e->type == EXPR_BIGINT)  return mpz_odd_p(e->data.bigint) ? 1 : 0;
    return -1;
}

#define IVP_MAX_TERMS 512
#define IVP_MAX_GENS   32

/* Decompose a single expanded monomial `term` into an integer-coefficient
 * parity and its set of (provably-integer) generators. Returns:
 *   1  success (coeff_odd and gens filled; gens are borrowed pointers)
 *   0  coefficient even -> term vanishes mod 2 (gens irrelevant)
 *  -1  not an integer-coefficient monomial in integer generators -> decline. */
static int ivp_term_decompose(const AssumeCtx* ctx, const Expr* term,
                              const Expr** gens, size_t* ngen) {
    *ngen = 0;
    int coeff_odd = 1;
    const Expr* const* facs; size_t nf; const Expr* single = term;
    if (sr_head(term, "Times")) { nf = term->data.function.arg_count;
        facs = (const Expr* const*)term->data.function.args; }
    else { nf = 1; facs = &single; }
    /* An even coefficient makes the term vanish mod 2 only once every other
     * factor is shown to be an integer: returning 0 at the coefficient let
     * (-1)^(2 n) -> 1 -- so E^(2 I n Pi) -> 1 and Arg[Exp[I Pi (2 n + 1/2)]]
     * -> Pi/2 -- for a merely real n. So scan every factor first. */
    bool even = false;
    for (size_t i = 0; i < nf; i++) {
        const Expr* f = facs[i];
        int lp = lit_parity(f);
        if (lp >= 0) { if (lp == 0) even = true; continue; }   /* integer coeff */
        /* a generator factor: g or Power[g, positive-integer] */
        const Expr* g = f;
        if (sr_head(f, "Power") && f->data.function.arg_count == 2) {
            const Expr* k = f->data.function.args[1];
            if (k->type != EXPR_INTEGER || k->data.integer <= 0) return -1;
            g = f->data.function.args[0];               /* g^k == g (mod 2) */
        }
        /* reject anything numeric-but-not-integer (Rational, Real) */
        if (g->type == EXPR_INTEGER || g->type == EXPR_BIGINT) { continue; }
        if (g->type == EXPR_REAL || sr_head(g, "Rational")) return -1;
        if (!prov_int(ctx, g)) return -1;
        if (*ngen >= IVP_MAX_GENS) return -1;
        /* multilinear: include each distinct generator once */
        bool seen = false;
        for (size_t j = 0; j < *ngen; j++) if (expr_eq((Expr*)gens[j], (Expr*)g)) { seen = true; break; }
        if (!seen) gens[(*ngen)++] = g;
    }
    (void)coeff_odd;
    return even ? 0 : 1;
}

/* Boolean-ring reduction of an integer-valued polynomial E modulo 2.
 * Returns a fresh reduced exponent (multilinear, {0,1} coefficients), or
 * NULL to decline. */
static Expr* simp_reduce_mod2(const AssumeCtx* ctx, const Expr* E) {
    Expr* pre = mod2_absorb(E);
    Expr* ex = expr_expand(pre);
    expr_free(pre);
    if (!ex) return NULL;

    const Expr* const* terms; size_t nt; const Expr* single = ex;
    if (sr_head(ex, "Plus")) { nt = ex->data.function.arg_count;
        terms = (const Expr* const*)ex->data.function.args; }
    else { nt = 1; terms = &single; }
    if (nt > IVP_MAX_TERMS) { expr_free(ex); return NULL; }

    /* Accumulate surviving monomials (canonicalised) with XOR semantics. */
    Expr** mons = (Expr**)malloc((nt ? nt : 1) * sizeof(Expr*));
    size_t nmon = 0;
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* gens[IVP_MAX_GENS]; size_t ng = 0;
        int st = ivp_term_decompose(ctx, terms[i], gens, &ng);
        if (st < 0) { ok = false; break; }
        if (st == 0) continue;                      /* even coefficient */
        /* canonical monomial = Times[gens...] (empty -> 1) */
        Expr* mon;
        if (ng == 0) mon = expr_new_integer(1);
        else if (ng == 1) mon = expr_copy((Expr*)gens[0]);
        else {
            Expr** gc = (Expr**)malloc(ng * sizeof(Expr*));
            for (size_t j = 0; j < ng; j++) gc[j] = expr_copy((Expr*)gens[j]);
            mon = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), gc, ng));
            free(gc);
        }
        /* XOR into the accumulator. */
        size_t found = nmon;
        for (size_t j = 0; j < nmon; j++) if (expr_eq(mons[j], mon)) { found = j; break; }
        if (found < nmon) { expr_free(mons[found]); mons[found] = mons[--nmon]; expr_free(mon); }
        else mons[nmon++] = mon;
    }
    expr_free(ex);
    if (!ok) { for (size_t j = 0; j < nmon; j++) expr_free(mons[j]); free(mons); return NULL; }

    Expr* R;
    if (nmon == 0) R = expr_new_integer(0);
    else if (nmon == 1) { R = mons[0]; }
    else R = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), mons, nmon));
    free(mons);
    return R;
}

/* Canonical value of (-1)^Y where Y is a rational-linear combination of
 * provably-integer generators. The non-integer rational-constant part of Y
 * is kept as a literal factor (-1)^c (e.g. (-1)^(1/2) = I); the rest is
 * reduced mod 2. Returns a fresh, evaluated Expr, or NULL to decline. */
static Expr* canon_neg1_pow(const AssumeCtx* ctx, const Expr* Y) {
    const Expr* const* terms; size_t nt; const Expr* single = Y;
    if (sr_head(Y, "Plus")) { nt = Y->data.function.arg_count;
        terms = (const Expr* const*)Y->data.function.args; }
    else { nt = 1; terms = &single; }

    Expr** ratc = (Expr**)malloc((nt ? nt : 1) * sizeof(Expr*)); size_t nrat = 0;
    Expr** intp = (Expr**)malloc((nt ? nt : 1) * sizeof(Expr*)); size_t nip = 0;
    for (size_t i = 0; i < nt; i++) {
        const Expr* t = terms[i];
        if (t->type == EXPR_REAL) {                 /* inexact: bail */
            for (size_t j = 0; j < nrat; j++) expr_free(ratc[j]);
            for (size_t j = 0; j < nip;  j++) expr_free(intp[j]);
            free(ratc); free(intp); return NULL;
        }
        /* Every numeric constant (integer or rational) becomes a literal
         * (-1)^c factor: an odd integer constant thus yields a -1 pulled OUT
         * as Times[-1, ...], the form the evaluator keeps and cancels against
         * a sibling (-1)^E (it does NOT re-absorb it into the exponent). */
        if (t->type == EXPR_INTEGER || t->type == EXPR_BIGINT || sr_head(t, "Rational"))
            ratc[nrat++] = expr_copy((Expr*)t);
        else
            intp[nip++] = expr_copy((Expr*)t);
    }

    Expr* Yint = (nip == 0) ? expr_new_integer(0)
               : (nip == 1) ? intp[0]
               : expr_new_function(expr_new_symbol(SYM_Plus), intp, nip);
    free(intp);
    Expr* R = simp_reduce_mod2(ctx, Yint);
    expr_free(Yint);
    if (!R) { for (size_t j = 0; j < nrat; j++) expr_free(ratc[j]); free(ratc); return NULL; }

    /* (-1)^R times (-1)^(sum of rational constants). */
    Expr* factors[2]; size_t nfac = 0;
    factors[nfac++] = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_new_integer(-1), R }, 2);
    if (nrat > 0) {
        Expr* c = (nrat == 1) ? ratc[0]
                : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ratc, nrat));
        factors[nfac++] = expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ expr_new_integer(-1), c }, 2);
    }
    free(ratc);
    Expr* out = (nfac == 1) ? factors[0]
              : expr_new_function(expr_new_symbol(SYM_Times), factors, nfac);
    return eval_and_free(out);
}

/* True iff a Complex[...] node appears anywhere in e (head or args). */
static bool sr_has_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sr_has_complex(e->data.function.args[i])) return true;
    return false;
}

/* If `ex` is a pure imaginary multiple of Pi, i.e. ex == I Pi X, return X
 * (I- and Pi-free), else NULL. Used for E^(I Pi X) = (-1)^X. */
static Expr* exp_ipi_coeff(const Expr* ex) {
    if (!sr_mentions_pi(ex)) return NULL;
    Expr* denom = expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Complex),
                                     (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2),
                                 expr_new_symbol(SYM_Pi) }, 2);
    Expr* dinv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ denom, expr_new_integer(-1) }, 2);
    Expr* X = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_copy((Expr*)ex), dinv }, 2));
    if (sr_mentions_pi(X) || sr_has_complex(X)) { expr_free(X); return NULL; }
    return X;
}

/* Reduce H[t + K Pi] for a trig head H, K a provable integer: Tan/Cot have
 * period Pi (no sign); Sin/Cos/Sec/Csc flip sign with the parity of K. Returns
 * a fresh Expr on success, NULL if there is no integer-Pi part or the parity
 * cannot be determined. */
static Expr* trig_pi_shift_symbolic(const AssumeCtx* ctx, const char* h,
                                    const Expr* arg) {
    bool is_tan_cot = (h == SYM_Tan || h == SYM_Cot);
    bool is_sin_cos_sec_csc = (h == SYM_Sin || h == SYM_Cos ||
                               h == SYM_Sec || h == SYM_Csc);
    if (!is_tan_cot && !is_sin_cos_sec_csc) return NULL;

    const Expr* single = arg;
    size_t nt; const Expr* const* terms;
    if (sr_head(arg, "Plus")) {
        nt = arg->data.function.arg_count;
        terms = (const Expr* const*)arg->data.function.args;
    } else { nt = 1; terms = &single; }

    Expr** kco = (Expr**)malloc(nt * sizeof(Expr*)); size_t nk = 0;
    Expr** theta = (Expr**)malloc(nt * sizeof(Expr*)); size_t nth = 0;
    for (size_t i = 0; i < nt; i++) {
        Expr* k = pi_integer_coeff(terms[i], ctx);
        if (k) kco[nk++] = k;
        else theta[nth++] = expr_copy((Expr*)terms[i]);
    }
    if (nk == 0) {
        for (size_t i = 0; i < nth; i++) expr_free(theta[i]);
        free(kco); free(theta);
        return NULL;
    }

    Expr* K = (nk == 1) ? kco[0]
            : expr_new_function(expr_new_symbol(SYM_Plus), kco, nk);
    free(kco);
    K = eval_and_free(K);
    /* K is a provable integer by construction (pi_integer_coeff gated on
     * prov_int), so H[t + K Pi] = (-1)^K H[t] is EXACT whether or not the
     * parity of K is known as a literal. Tan/Cot have period Pi (no sign).
     * Emit the symbolic (-1)^K when the parity is undetermined rather than
     * declining; the root-of-unity canonicalizer reduces it further. */
    int parity = is_tan_cot ? 0 : pi_parity(ctx, K);

    Expr* th = (nth == 0) ? expr_new_integer(0)
             : (nth == 1) ? theta[0]
             : expr_new_function(expr_new_symbol(SYM_Plus), theta, nth);
    free(theta);
    Expr* hcall = expr_new_function(expr_new_symbol(h), (Expr*[]){ th }, 1);
    if (parity == 0) { expr_free(K); return hcall; }
    if (parity == 1) { expr_free(K);
        return expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), hcall }, 2); }
    /* parity undetermined: (-1)^K H[t], K consumed. */
    Expr* sign = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_new_integer(-1), K }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
               (Expr*[]){ sign, hcall }, 2);
}

/* If `term` is (integer)*Pi*I, return its integer coefficient m as a fresh Expr
 * (caller owns); else NULL. Computed as m = term/(Pi I) (evaluated). A real
 * integer-Pi term k*Pi yields the imaginary quotient -k I, which prov_int
 * rejects, so only a genuine imaginary multiple of Pi is accepted. */
static Expr* pi_i_integer_coeff(const Expr* term, const AssumeCtx* ctx) {
    if (!sr_mentions_pi(term)) return NULL;
    Expr* deninv = eval_and_free(parse_expression("1/(Pi*I)"));  /* = -I/Pi */
    if (!deninv) return NULL;
    Expr* m = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_copy((Expr*)term), deninv }, 2));
    if (sr_mentions_pi(m)) { expr_free(m); return NULL; }  /* not linear in Pi */
    if (prov_int(ctx, m)) return m;  /* rejects the imaginary -k I of a real k*Pi */
    expr_free(m);
    return NULL;
}

/* Reduce H[t + m Pi I] for a hyperbolic head H, m a provable integer. The real
 * hyperbolic functions have an IMAGINARY period: Sinh/Cosh/Sech/Csch period
 * 2 Pi I and flip sign with the parity of m (Sinh[u + m Pi I] = (-1)^m Sinh[u]);
 * Tanh/Coth period Pi I (no sign). The circular I-bridge already peels a bare
 * Sinh[n Pi I] -> I Sin[n Pi]; this handles the shifted argument t + m Pi I with
 * a symbolic integer m that the bridge leaves alone. Returns a fresh Expr on
 * success, NULL if there is no integer-Pi-I part or the parity is undetermined. */
static Expr* hyp_pi_i_shift_symbolic(const AssumeCtx* ctx, const char* h,
                                     const Expr* arg) {
    bool is_tanh_coth = (h == SYM_Tanh || h == SYM_Coth);
    bool is_sign_flipping = (h == SYM_Sinh || h == SYM_Cosh ||
                             h == SYM_Sech || h == SYM_Csch);
    if (!is_tanh_coth && !is_sign_flipping) return NULL;

    const Expr* single = arg;
    size_t nt; const Expr* const* terms;
    if (sr_head(arg, "Plus")) {
        nt = arg->data.function.arg_count;
        terms = (const Expr* const*)arg->data.function.args;
    } else { nt = 1; terms = &single; }

    Expr** mco = (Expr**)malloc(nt * sizeof(Expr*)); size_t nm = 0;
    Expr** theta = (Expr**)malloc(nt * sizeof(Expr*)); size_t nth = 0;
    for (size_t i = 0; i < nt; i++) {
        Expr* m = pi_i_integer_coeff(terms[i], ctx);
        if (m) mco[nm++] = m;
        else theta[nth++] = expr_copy((Expr*)terms[i]);
    }
    if (nm == 0) {
        for (size_t i = 0; i < nth; i++) expr_free(theta[i]);
        free(mco); free(theta);
        return NULL;
    }

    Expr* M = (nm == 1) ? mco[0]
            : expr_new_function(expr_new_symbol(SYM_Plus), mco, nm);
    free(mco);
    M = eval_and_free(M);
    int parity = is_tanh_coth ? 0 : pi_parity(ctx, M);
    expr_free(M);
    if (parity < 0) {
        for (size_t i = 0; i < nth; i++) expr_free(theta[i]);
        free(theta);
        return NULL;
    }

    Expr* th = (nth == 0) ? expr_new_integer(0)
             : (nth == 1) ? theta[0]
             : expr_new_function(expr_new_symbol(SYM_Plus), theta, nth);
    free(theta);
    Expr* hcall = expr_new_function(expr_new_symbol(h), (Expr*[]){ th }, 1);
    if (parity == 1)
        return expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), hcall }, 2);
    return hcall;
}

/* ----------------------------------------------------------------------- */
/* Inverse-trig-of-trig adjacent-branch reduction: ArcSin[Sin[x]] -> Pi - x, */
/* ArcTan[Tan[x]] -> x - Pi, ArcCos[Cos[x]] -> x (sub-interval) and the like, */
/* when x's numeric interval lies in a single (possibly shifted) branch. The  */
/* principal-range checks above only match a literal bound fact, so they miss */
/* a sub-interval like Pi/2 < x < Pi (inside [0,Pi]) and every shifted branch.*/
/* ----------------------------------------------------------------------- */

/* Numeric value of a bound expression, numericalizing symbolic constants such
 * as Pi/2 and Pi via N[] (sr_num only reads plain literals). */
static bool sr_num_eval(const Expr* e, double* v) {
    if (sr_num(e, v)) return true;
    Expr* ne = eval_and_free(expr_new_function(expr_new_symbol("N"),
                   (Expr*[]){ expr_copy((Expr*)e) }, 1));
    if (!ne) return false;
    bool ok = sr_num(ne, v);
    if (!ok && ne->type == EXPR_REAL) { *v = ne->data.real; ok = true; }
    expr_free(ne);
    return ok;
}

/* Like sr_bounds but numericalizes symbolic bounds (Pi/2, Pi, ...). */
static void sr_bounds_eval(const AssumeCtx* ctx, const Expr* sym,
                           bool* has_lo, double* lo,
                           bool* has_hi, double* hi) {
    *has_lo = *has_hi = false;
    if (!ctx) return;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL ||
            f->data.function.arg_count != 2) continue;
        const char* h = f->data.function.head->data.symbol.name;
        const Expr* A = f->data.function.args[0];
        const Expr* B = f->data.function.args[1];
        double k; bool sym_left; int kind;
        if (expr_eq((Expr*)A, (Expr*)sym) && sr_num_eval(B, &k)) sym_left = true;
        else if (expr_eq((Expr*)B, (Expr*)sym) && sr_num_eval(A, &k)) sym_left = false;
        else continue;
        if (h == SYM_Less || h == SYM_LessEqual)
            kind = sym_left ? +1 : -1;
        else if (h == SYM_Greater || h == SYM_GreaterEqual)
            kind = sym_left ? -1 : +1;
        else continue;
        if (kind < 0) { if (!*has_lo || k > *lo) { *has_lo = true; *lo = k; } }
        else          { if (!*has_hi || k < *hi) { *has_hi = true; *hi = k; } }
    }
}

/* Build e - k Pi (k integer), unevaluated. */
static Expr* invtrig_e_minus_kpi(const Expr* e, long k) {
    if (k == 0) return expr_copy((Expr*)e);
    Expr* kpi = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer((int64_t)(-k)),
                               expr_new_symbol("Pi") }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ expr_copy((Expr*)e), kpi }, 2);
}

/* Reduce ArcF[F[inner]] with inner a bare symbol bounded to [lo, hi]. Returns a
 * fresh reduced Expr, or NULL if the interval is not contained in one branch. */
static Expr* invtrig_branch_reduce(const AssumeCtx* ctx, const char* arc,
                                   const Expr* inner) {
    if (inner->type != EXPR_SYMBOL) return NULL;
    bool hl = false, hu = false; double lo = 0, hi = 0;
    sr_bounds_eval(ctx, inner, &hl, &lo, &hu, &hi);
    if (!hl || !hu) return NULL;            /* need a bounded interval */
    const double PI = 3.14159265358979323846;
    const double T = 1e-9;
    Expr* res = NULL;

    if (arc == SYM_ArcSin) {
        /* ArcSin[Sin[x]] = (-1)^k (x - k Pi) on [k Pi - Pi/2, k Pi + Pi/2]. */
        long k = lround((lo + hi) * 0.5 / PI);
        if (lo >= k * PI - PI / 2 - T && hi <= k * PI + PI / 2 + T) {
            Expr* emk = invtrig_e_minus_kpi(inner, k);
            res = (k & 1) ? expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){ expr_new_integer(-1), emk }, 2) : emk;
        }
    } else if (arc == SYM_ArcTan) {
        /* ArcTan[Tan[x]] = x - k Pi on (k Pi - Pi/2, k Pi + Pi/2). */
        long k = lround((lo + hi) * 0.5 / PI);
        if (lo >= k * PI - PI / 2 - T && hi <= k * PI + PI / 2 + T)
            res = invtrig_e_minus_kpi(inner, k);
    } else if (arc == SYM_ArcCot) {
        /* ArcCot[Cot[x]] = x - k Pi on (k Pi, (k+1) Pi). */
        long k = (long)floor((lo + hi) * 0.5 / PI);
        if (lo >= k * PI - T && hi <= (k + 1) * PI + T)
            res = invtrig_e_minus_kpi(inner, k);
    } else if (arc == SYM_ArcCos) {
        /* ArcCos[Cos[x]] is even, period 2Pi: x - 2m Pi on [2m Pi,(2m+1)Pi];
         * 2m Pi - x on [(2m-1)Pi, 2m Pi]. */
        long seg = (long)floor((lo + hi) * 0.5 / PI);
        if ((seg & 1) == 0) {                   /* ascending half, 2m = seg */
            if (lo >= seg * PI - T && hi <= (seg + 1) * PI + T)
                res = invtrig_e_minus_kpi(inner, seg);
        } else {                                /* descending half, 2m = seg+1 */
            long twom = seg + 1;
            if (lo >= (twom - 1) * PI - T && hi <= twom * PI + T) {
                Expr* emk = invtrig_e_minus_kpi(inner, twom);   /* x - 2m Pi */
                res = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(-1), emk }, 2); /* 2m Pi - x */
            }
        }
    }
    return res ? eval_and_free(res) : NULL;
}

/* ----------------------------------------------------------------------- */
/* Reduce-verified inverse-trig multiple-angle: ArcCos[1-2x^2] -> 2 ArcSin[x], */
/* ArcCos[2x^2-1] -> Pi - 2 ArcCos[x] (= 2 ArcCos[x]), ArcSin[2x Sqrt[1-x^2]]  */
/* -> 2 ArcSin[x], etc. The candidate is applied ONLY after Reduce proves the  */
/* equality on the assumed region, so it is sound by construction regardless   */
/* of which branch the region sits in.                                         */
/* ----------------------------------------------------------------------- */

/* coeff * InvH[x]  (+ Pi if pi_off). */
static Expr* mk_scaled_inv(const char* invh, int coeff, int pi_off, const Expr* x) {
    Expr* inv = expr_new_function(expr_new_symbol(invh),
                    (Expr*[]){ expr_copy((Expr*)x) }, 1);
    Expr* term = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(coeff), inv }, 2);
    if (!pi_off) return term;
    return expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ expr_new_symbol("Pi"), term }, 2);
}

static Expr* invtrig_multiangle_reduce(const AssumeCtx* ctx, const char* arc,
                                       const Expr* u) {
    if (arc != SYM_ArcSin && arc != SYM_ArcCos) return NULL;
    /* Argument must be a polynomial/radical (Plus/Times/Power), not a trig or
     * bare symbol -- the multiple-angle arguments 1-2x^2, 2x^2-1, 2x Sqrt[...]. */
    if (!sr_head(u, "Plus") && !sr_head(u, "Times") && !sr_head(u, "Power"))
        return NULL;
    const Expr** syms = NULL; size_t ns = 0, cap = 0;
    sr_collect_syms(u, &syms, &ns, &cap);
    int one = (ns == 1);
    const Expr* x = one ? syms[0] : NULL;
    free(syms);
    if (!one) return NULL;
    bool hl = false, hu = false; double lo = 0, hi = 0;
    sr_bounds_eval(ctx, x, &hl, &lo, &hu, &hi);
    if (!hl || !hu) return NULL;             /* need a bounded region */

    static const struct { const char* h; int c; int pi; } cand[4] = {
        { "ArcSin", 2, 0 }, { "ArcSin", -2, 1 },
        { "ArcCos", 2, 0 }, { "ArcCos", -2, 1 } };
    Expr* arcu = expr_new_function(expr_new_symbol(arc),
                     (Expr*[]){ expr_copy((Expr*)u) }, 1);
    Expr* result = NULL;
    for (int i = 0; i < 4 && !result; i++) {
        Expr* c = mk_scaled_inv(cand[i].h, cand[i].c, cand[i].pi, x);
        Expr* eq = expr_new_function(expr_new_symbol("Equal"),
                       (Expr*[]){ expr_copy(arcu), expr_copy(c) }, 2);
        int v = assume_reduce_entails(ctx, eq);   /* 1 iff A |= ArcF[u] == c */
        expr_free(eq);
        if (v == 1) result = eval_and_free(c);
        else expr_free(c);
    }
    expr_free(arcu);
    return result;
}

/* Bottom-up walk applying the structural assumption rewrites. Always returns a
 * freshly owned tree; sets *changed when any node was rewritten. */
static Expr* assume_structural_rewrite(const Expr* e, const AssumeCtx* ctx, int* changed) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        args[i] = assume_structural_rewrite(e->data.function.args[i], ctx, changed);
    Expr* node = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);

    if (node->data.function.head->type != EXPR_SYMBOL) return node;
    const char* h = node->data.function.head->data.symbol.name;
    size_t nn = node->data.function.arg_count;

    if (nn == 1) {
        Expr* a0 = node->data.function.args[0];

        /* Symbolic-integer Pi periodicity: Sin[t + 2 n Pi] -> Sin[t],
         * Sin[t + (2n+1) Pi] -> -Sin[t], Tan[t + n Pi] -> Tan[t], etc. The
         * coefficient-1 string rules only match a bare Times[n, Pi]; this
         * handles any integer-linear Pi coefficient (2n, 2n+1, ...). */
        if ((h == SYM_Sin || h == SYM_Cos ||
             h == SYM_Tan || h == SYM_Cot ||
             h == SYM_Sec || h == SYM_Csc) && ctx) {
            Expr* out = trig_pi_shift_symbolic(ctx, h, a0);
            if (out) { expr_free(node); *changed = 1; return out; }
        }

        /* Imaginary-period shift for the hyperbolic family: Sinh[t + 2 n Pi I]
         * -> Sinh[t], Tanh[t + n Pi I] -> Tanh[t], etc. (period 2 Pi I for
         * Sinh/Cosh/Sech/Csch, Pi I for Tanh/Coth). */
        if ((h == SYM_Sinh || h == SYM_Cosh ||
             h == SYM_Tanh || h == SYM_Coth ||
             h == SYM_Sech || h == SYM_Csch) && ctx) {
            Expr* out = hyp_pi_i_shift_symbolic(ctx, h, a0);
            if (out) { expr_free(node); *changed = 1; return out; }
        }

        /* Integer-valued roundings of a provably-integer argument. */
        if ((h == SYM_Floor || h == SYM_Ceiling ||
             h == SYM_Round || h == SYM_IntegerPart) &&
            assume_known_integer(ctx, a0)) {
            Expr* out = expr_copy(a0); expr_free(node); *changed = 1; return out;
        }
        if (h == SYM_FractionalPart && assume_known_integer(ctx, a0)) {
            expr_free(node); *changed = 1; return expr_new_integer(0);
        }

        /* Double factorial of an even integer argument: (2q)!! = 2^q q!,
         * gated on q = arg/2 being a provable integer (so it declines for a
         * generic/odd argument). The odd case is handled elsewhere. */
        if (h == SYM_Factorial2 && ctx) {
            Expr* half = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(a0),
                           expr_new_function(expr_new_symbol(SYM_Rational),
                               (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2));
            if (prov_int(ctx, half)) {
                Expr* pow2 = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_integer(2), expr_copy(half) }, 2);
                Expr* fac = expr_new_function(expr_new_symbol("Factorial"),
                    (Expr*[]){ half }, 1);
                Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ pow2, fac }, 2));
                expr_free(node); *changed = 1; return out;
            }
            expr_free(half);
        }

        /* Integer-shift extraction: Floor[k + t] -> k + Floor[t] (same for
         * Ceiling/Round), FractionalPart[k + t] -> FractionalPart[t], pulling
         * the provably-integer part k out of a sum argument. Sound for any
         * integer k; fires only when there is both an integer part and a
         * non-integer remainder. */
        if ((h == SYM_Floor || h == SYM_Ceiling ||
             h == SYM_Round || h == SYM_FractionalPart) &&
            ctx && sr_head(a0, "Plus")) {
            size_t na = a0->data.function.arg_count;
            Expr** kp = (Expr**)malloc(na * sizeof(Expr*)); size_t nkp = 0;
            Expr** tp = (Expr**)malloc(na * sizeof(Expr*)); size_t ntp = 0;
            for (size_t i = 0; i < na; i++) {
                Expr* term = a0->data.function.args[i];
                if (prov_int(ctx, term)) kp[nkp++] = expr_copy(term);
                else                     tp[ntp++] = expr_copy(term);
            }
            if (nkp > 0 && ntp > 0) {
                Expr* K = (nkp == 1) ? kp[0]
                        : expr_new_function(expr_new_symbol(SYM_Plus), kp, nkp);
                Expr* T = (ntp == 1) ? tp[0]
                        : expr_new_function(expr_new_symbol(SYM_Plus), tp, ntp);
                free(kp); free(tp);
                Expr* out;
                if (h == SYM_FractionalPart) {
                    expr_free(K);
                    out = expr_new_function(expr_new_symbol("FractionalPart"),
                              (Expr*[]){ T }, 1);
                } else {
                    Expr* inner = expr_new_function(expr_new_symbol(h),
                                      (Expr*[]){ T }, 1);
                    out = expr_new_function(expr_new_symbol(SYM_Plus),
                              (Expr*[]){ K, inner }, 2);
                }
                out = eval_and_free(out);
                expr_free(node); *changed = 1; return out;
            }
            for (size_t i = 0; i < nkp; i++) expr_free(kp[i]);
            for (size_t i = 0; i < ntp; i++) expr_free(tp[i]);
            free(kp); free(tp);
        }

        /* Floor/Ceiling of an exact division by a positive integer:
         *   Floor[P/d]   -> (P - Mod[P,d]) / d
         *   Ceiling[P/d] -> (P + Mod[-P,d]) / d
         * for a provably-integer P and integer denominator d > 1. This is the
         * definition of floor-division; the Mod residue is reduced by the
         * Mod rules below (iterated by the Simplify confluence post-pass), so
         * e.g. Floor[n/2] + Ceiling[n/2] -> n. */
        if ((h == SYM_Floor || h == SYM_Ceiling) && ctx) {
            Expr* den = eval_and_free(expr_new_function(
                            expr_new_symbol("Denominator"),
                            (Expr*[]){ expr_copy(a0) }, 1));
            if (den->type == EXPR_INTEGER && den->data.integer > 1) {
                Expr* P = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_copy(a0), expr_copy(den) }, 2));
                if (prov_int(ctx, P)) {
                    bool is_floor = (h == SYM_Floor);
                    /* residue argument: P (Floor) or -P (Ceiling) */
                    Expr* resarg = is_floor ? expr_copy(P)
                        : eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(-1), expr_copy(P) }, 2));
                    Expr* modr = expr_new_function(expr_new_symbol("Mod"),
                                     (Expr*[]){ resarg, expr_copy(den) }, 2);
                    /* P +/- Mod, all over d */
                    Expr* signed_mod = is_floor
                        ? expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(-1), modr }, 2)
                        : modr;
                    Expr* numer = expr_new_function(expr_new_symbol(SYM_Plus),
                                      (Expr*[]){ expr_copy(P), signed_mod }, 2);
                    Expr* out = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ numer, expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
                    /* Distribute the 1/d so the P/d term collects against an
                     * enclosing Plus (Floor[n/2] + Ceiling[n/2] - n -> 0). */
                    out = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                              (Expr*[]){ out }, 1));
                    expr_free(P); expr_free(den);
                    expr_free(node); *changed = 1; return out;
                }
                expr_free(P);
            }
            expr_free(den);
        }

        /* Ceiling[x] / Floor[x] pinned to a single integer by numeric bounds. */
        if ((h == SYM_Ceiling || h == SYM_Floor) &&
            a0->type == EXPR_SYMBOL) {
            bool hl, hu, ls = false, us = false; double lo = 0, hi = 0;
            sr_bounds(ctx, a0, &hl, &lo, &ls, &hu, &hi, &us);
            if (h == SYM_Ceiling && hu && hl) {
                long c = (long)ceil(hi);
                if ((double)c >= hi &&                       /* s <= c */
                    (lo > (double)(c - 1) || (lo == (double)(c - 1) && ls))) { /* s > c-1 */
                    expr_free(node); *changed = 1; return expr_new_integer(c);
                }
            }
            if (h == SYM_Floor && hl && hu) {
                long fl = (long)floor(lo);
                if ((double)fl <= lo &&                      /* s >= f */
                    (hi < (double)(fl + 1) || (hi == (double)(fl + 1) && us))) { /* s < f+1 */
                    expr_free(node); *changed = 1; return expr_new_integer(fl);
                }
            }
        }

        /* FractionalPart[a] from a Mod[a,1] fact plus the sign of a. */
        if (h == SYM_FractionalPart && a0->type == EXPR_SYMBOL && ctx) {
            for (size_t i = 0; i < ctx->count; i++) {
                const Expr* f = ctx->facts[i];
                if (!fact_is_function(f, "Equal", 2)) continue;
                const Expr* L = f->data.function.args[0];
                const Expr* R = f->data.function.args[1];
                const Expr* modv = NULL; const Expr* rv = NULL;
                if (sr_head(L, "Mod") && L->data.function.arg_count == 2) { modv = L; rv = R; }
                else if (sr_head(R, "Mod") && R->data.function.arg_count == 2) { modv = R; rv = L; }
                if (!modv) continue;
                if (!expr_eq((Expr*)modv->data.function.args[0], a0)) continue;
                double one; const Expr* m1 = modv->data.function.args[1];
                if (!(sr_num(m1, &one) && one == 1.0)) continue;
                double r;
                if (!sr_num(rv, &r) || !(r > 0.0 && r < 1.0)) continue;
                Expr* out = NULL;
                if (assume_known_negative(ctx, a0))
                    out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                              (Expr*[]){ expr_copy((Expr*)rv), expr_new_integer(-1) }, 2));
                else if (assume_known_nonneg(ctx, a0))
                    out = expr_copy((Expr*)rv);
                if (out) { expr_free(node); *changed = 1; return out; }
            }
        }

        /* Arg of a quantity whose sign/axis the assumptions pin, which
         * ComplexExpand (sign-agnostic) cannot settle: Arg[positive] = 0,
         * Arg[negative] = Pi, and Arg[(i*c) * positive] = +-Pi/2 (the sign from
         * c).  prov_pos already descends Times/Power, so a^2, 2a, 1/a^2, ... are
         * recognised as positive under a > 0.  Fixes the Arg[...]-laden keyhole
         * forms (e.g. Log[x]/(x^2+a^2), x/Sinh[a x]). */
        if (h == SYM_Arg && ctx) {
            if (assume_known_positive(ctx, a0)) {
                expr_free(node); *changed = 1; return expr_new_integer(0);
            }
            if (assume_known_negative(ctx, a0)) {
                expr_free(node); *changed = 1;
                return eval_and_free(parse_expression("Pi"));
            }
            if (sr_head(a0, "Times")) {
                size_t n = a0->data.function.arg_count;
                double imv = 0; bool got = false;
                Expr** rp = malloc(n * sizeof(Expr*)); size_t nr = 0;
                for (size_t i = 0; i < n; i++) {
                    Expr* fe = a0->data.function.args[i];
                    if (!got && sr_pure_imag(fe, &imv)) { got = true; continue; }
                    rp[nr++] = expr_copy(fe);
                }
                if (got && imv != 0.0) {
                    Expr* r = (nr == 0) ? expr_new_integer(1)
                            : (nr == 1) ? rp[0]
                            : expr_new_function(expr_new_symbol(SYM_Times), rp, nr);
                    free(rp);
                    bool pos = assume_known_positive(ctx, r);
                    expr_free(r);
                    if (pos) {
                        expr_free(node); *changed = 1;
                        return eval_and_free(parse_expression(imv > 0.0 ? "Pi/2" : "-Pi/2"));
                    }
                } else {
                    for (size_t i = 0; i < nr; i++) expr_free(rp[i]);
                    free(rp);
                }
            }
        }

        /* Schwarz reflection: Conjugate[H[z]] -> H[Conjugate[z]] for H entire or
         * meromorphic with real Taylor coefficients (trig, hyperbolic, Exp), when
         * every free symbol of z is real. ComplexExpand of a pole-bearing head
         * (Tan, ...) yields a double-angle rational form that does NOT structurally
         * match the single-angle form the round loop gives the sibling H[x - I y],
         * so the difference would not cancel; the reflection produces the matching
         * single-angle form directly. Runs BEFORE the generic ComplexExpand block
         * below so it wins for this safe family. */
        if (h == SYM_Conjugate && a0->type == EXPR_FUNCTION &&
            a0->data.function.head &&
            a0->data.function.head->type == EXPR_SYMBOL &&
            a0->data.function.arg_count == 1) {
            const char* fh = a0->data.function.head->data.symbol.name;
            bool reflect_safe =
                fh == SYM_Sin  || fh == SYM_Cos  ||
                fh == SYM_Tan  || fh == SYM_Cot  ||
                fh == SYM_Sec  || fh == SYM_Csc  ||
                fh == SYM_Sinh || fh == SYM_Cosh ||
                fh == SYM_Tanh || fh == SYM_Coth ||
                fh == SYM_Sech || fh == SYM_Csch ||
                fh == SYM_Exp;
            if (reflect_safe) {
                const Expr* z = a0->data.function.args[0];
                const Expr** syms = NULL; size_t ns = 0, cap = 0;
                sr_collect_syms(z, &syms, &ns, &cap);
                bool all_real = (ns > 0);
                for (size_t i = 0; i < ns; i++)
                    if (!prov_re(ctx, syms[i])) { all_real = false; break; }
                free(syms);
                if (all_real) {
                    /* Conjugate[z] via ComplexExpand (z real symbols => conj). */
                    Expr* cz = eval_and_free(expr_new_function(
                        expr_new_symbol("ComplexExpand"),
                        (Expr*[]){ expr_new_function(expr_new_symbol("Conjugate"),
                            (Expr*[]){ expr_copy((Expr*)z) }, 1) }, 1));
                    Expr* out = eval_and_free(expr_new_function(
                        expr_copy(a0->data.function.head),
                        (Expr*[]){ cz }, 1));
                    expr_free(node); *changed = 1; return out;
                }
            }
        }

        /* Re/Im/Conjugate/Arg/Abs of an expression whose every free symbol is
         * real: ComplexExpand assumes reality, so it yields the refined form
         * (e.g. Abs[a + b I] -> Sqrt[a^2 + b^2]). When ComplexExpand cannot
         * reduce (unknown sign, as in Abs[x]) it returns the input unchanged,
         * which the expr_eq guard below discards. */
        /* Re / Im of a provably-real quantity, or of a sum with provably-real
         * terms (Re and Im are additive).  Runs before ComplexExpand, which is
         * sign-agnostic: it writes Re[Log[b]] as Log[b^2]/2 even under b > 0,
         * where the assumptions make Log[b] itself real.  This is what turns the
         * principal value Re[-I Pi - Log[a] + Log[b]] (a, b > 0) into
         * Log[b] - Log[a] instead of (Log[b^2] - Log[a^2])/2. */
        if ((h == SYM_Re || h == SYM_Im) && ctx) {
            bool is_re = h == SYM_Re;
            if (prov_re(ctx, a0)) {
                Expr* out = is_re ? expr_copy(a0) : expr_new_integer(0);
                expr_free(node); *changed = 1; return out;
            }
            if (sr_head(a0, "Plus")) {
                size_t n = a0->data.function.arg_count;
                Expr** re = malloc(n * sizeof(Expr*)); size_t nre = 0;
                Expr** rest = malloc(n * sizeof(Expr*)); size_t nrest = 0;
                for (size_t i = 0; i < n; i++) {
                    Expr* t = a0->data.function.args[i];
                    if (prov_re(ctx, t)) re[nre++] = expr_copy(t);
                    else rest[nrest++] = expr_copy(t);
                }
                if (nre > 0 && nrest > 0) {
                    Expr* rsum = expr_new_function(expr_new_symbol(SYM_Plus), rest, nrest);
                    Expr* part0 = eval_and_free(expr_new_function(
                        expr_new_symbol(is_re ? SYM_Re : SYM_Im), (Expr*[]){ rsum }, 1));
                    /* The remainder gets the rest of this rewrite (the
                     * ComplexExpand path below for real-symbol remainders:
                     * Re[a + b I] -> a + Re[I b] -> a). */
                    int sub_changed = 0;
                    Expr* part = assume_structural_rewrite(part0, ctx, &sub_changed);
                    expr_free(part0);
                    Expr* out;
                    if (is_re) {
                        re = realloc(re, (nre + 1) * sizeof(Expr*));
                        re[nre++] = part;
                        out = expr_new_function(expr_new_symbol(SYM_Plus), re, nre);
                    } else {
                        for (size_t i = 0; i < nre; i++) expr_free(re[i]);
                        out = part;
                    }
                    free(re); free(rest);
                    out = eval_and_free(out);
                    expr_free(node); *changed = 1; return out;
                }
                for (size_t i = 0; i < nre; i++) expr_free(re[i]);
                for (size_t i = 0; i < nrest; i++) expr_free(rest[i]);
                free(re); free(rest);
            }
        }
        if ((h == SYM_Re || h == SYM_Im ||
             h == SYM_Conjugate || h == SYM_Arg ||
             h == SYM_Abs)) {
            const Expr** syms = NULL; size_t ns = 0, cap = 0;
            sr_collect_syms(a0, &syms, &ns, &cap);
            bool all_real = (ns > 0);
            for (size_t i = 0; i < ns; i++)
                if (!prov_re(ctx, syms[i])) { all_real = false; break; }
            free(syms);
            if (all_real) {
                Expr* ce = eval_and_free(expr_new_function(expr_new_symbol("ComplexExpand"),
                               (Expr*[]){ expr_copy(node) }, 1));
                if (ce && !expr_eq(ce, node)) { expr_free(node); *changed = 1; return ce; }
                if (ce) expr_free(ce);
            }
        }

        /* ArcTan[Tan[e]] -> e when Re[e] (or e itself) lies in (-Pi/2, Pi/2). */
        if (h == SYM_ArcTan && sr_head(a0, "Tan") &&
            a0->data.function.arg_count == 1) {
            const Expr* inner = a0->data.function.args[0];
            Expr* pihalf = eval_and_free(parse_expression("Pi/2"));
            Expr* neghalf = eval_and_free(parse_expression("-Pi/2"));
            Expr* u_re = eval_and_free(expr_new_function(expr_new_symbol(SYM_Re),
                             (Expr*[]){ expr_copy((Expr*)inner) }, 1));
            bool ok = (sr_proves_slt(ctx, neghalf, u_re) && sr_proves_slt(ctx, u_re, pihalf)) ||
                      (sr_proves_slt(ctx, neghalf, inner) && sr_proves_slt(ctx, inner, pihalf));
            expr_free(pihalf); expr_free(neghalf); expr_free(u_re);
            if (ok) { Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out; }
        }

        /* Inverse-of-forward for the hyperbolic family, mirroring ArcTan[Tan]
         * above.  Sinh and Tanh are bijective on the reals (ArcSinh/ArcTanh are
         * their real inverses), so ArcSinh[Sinh[e]] = ArcTanh[Tanh[e]] = e for
         * real e.  Cosh is even with ArcCosh[Cosh[e]] = |e|, so = e for e >= 0,
         * -e for e <= 0.  Each gate is an EXISTING prover; no new machinery. */
        /* Sinh/Tanh and their reciprocals Csch/Coth are odd bijections of the
         * reals (Csch/Coth on R\{0}), so Arc[.][.[e]] = e for real e. */
        if ((h == SYM_ArcSinh && sr_head(a0, "Sinh") &&
             a0->data.function.arg_count == 1) ||
            (h == SYM_ArcTanh && sr_head(a0, "Tanh") &&
             a0->data.function.arg_count == 1) ||
            (h == SYM_ArcCsch && sr_head(a0, "Csch") &&
             a0->data.function.arg_count == 1) ||
            (h == SYM_ArcCoth && sr_head(a0, "Coth") &&
             a0->data.function.arg_count == 1)) {
            const Expr* inner = a0->data.function.args[0];
            if (prov_re(ctx, inner)) {
                Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out;
            }
        }
        /* Cosh is even with ArcCosh range [0,inf); Sech is even with ArcSech
         * range [0,inf).  So Arc[.][.[e]] = |e| = e for e >= 0, -e for e <= 0. */
        if ((h == SYM_ArcCosh && sr_head(a0, "Cosh") &&
             a0->data.function.arg_count == 1) ||
            (h == SYM_ArcSech && sr_head(a0, "Sech") &&
             a0->data.function.arg_count == 1)) {
            const Expr* inner = a0->data.function.args[0];
            if (assume_known_nonneg(ctx, inner)) {
                Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out;
            }
            if (assume_known_nonpos(ctx, inner)) {
                Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)inner) }, 2));
                expr_free(node); *changed = 1; return out;
            }
        }

        /* Inverse-of-forward for the periodic circular family, mirroring
         * ArcTan[Tan] above. Each principal inverse is a bijection from exactly
         * one strip: ArcSin from Re in [-Pi/2, Pi/2], ArcCos from [0, Pi],
         * ArcCot from (0, Pi). So ArcSin[Sin[e]] = e, ArcCos[Cos[e]] = e,
         * ArcCot[Cot[e]] = e exactly when Re[e] (or e) lies in that strip. The
         * trig ranges are closed (identity holds at the endpoints) -> sr_proves_sle;
         * the ArcCot range is open (Cot singular at 0, Pi) -> sr_proves_slt. */
        if ((h == SYM_ArcSin && sr_head(a0, "Sin") &&
             a0->data.function.arg_count == 1) ||
            (h == SYM_ArcCos && sr_head(a0, "Cos") &&
             a0->data.function.arg_count == 1)) {
            const Expr* inner = a0->data.function.args[0];
            bool is_sin = (h == SYM_ArcSin);
            Expr* lo = eval_and_free(parse_expression(is_sin ? "-Pi/2" : "0"));
            Expr* hi = eval_and_free(parse_expression(is_sin ? "Pi/2" : "Pi"));
            Expr* u_re = eval_and_free(expr_new_function(expr_new_symbol(SYM_Re),
                             (Expr*[]){ expr_copy((Expr*)inner) }, 1));
            bool ok = (sr_proves_sle(ctx, lo, u_re) && sr_proves_sle(ctx, u_re, hi)) ||
                      (sr_proves_sle(ctx, lo, inner) && sr_proves_sle(ctx, inner, hi));
            expr_free(lo); expr_free(hi); expr_free(u_re);
            if (ok) { Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out; }
        }
        if (h == SYM_ArcCot && sr_head(a0, "Cot") &&
            a0->data.function.arg_count == 1) {
            const Expr* inner = a0->data.function.args[0];
            Expr* lo = eval_and_free(parse_expression("0"));
            Expr* hi = eval_and_free(parse_expression("Pi"));
            Expr* u_re = eval_and_free(expr_new_function(expr_new_symbol(SYM_Re),
                             (Expr*[]){ expr_copy((Expr*)inner) }, 1));
            bool ok = (sr_proves_slt(ctx, lo, u_re) && sr_proves_slt(ctx, u_re, hi)) ||
                      (sr_proves_slt(ctx, lo, inner) && sr_proves_slt(ctx, inner, hi));
            expr_free(lo); expr_free(hi); expr_free(u_re);
            if (ok) { Expr* out = expr_copy((Expr*)inner); expr_free(node); *changed = 1; return out; }
        }

        /* Adjacent-branch / sub-interval reduction for ArcF[F[x]] when x's
         * numeric interval lies in a single (possibly shifted) branch. Runs
         * after the principal-range checks above, which only match a literal
         * bound fact (so e.g. ArcCos[Cos[x]] on Pi/2 < x < Pi, a sub-interval of
         * [0,Pi], and ArcSin[Sin[x]] -> Pi - x on the k=1 branch fall here). */
        if (((h == SYM_ArcSin && sr_head(a0, "Sin")) ||
             (h == SYM_ArcCos && sr_head(a0, "Cos")) ||
             (h == SYM_ArcTan && sr_head(a0, "Tan")) ||
             (h == SYM_ArcCot && sr_head(a0, "Cot"))) &&
            a0->data.function.arg_count == 1 && ctx) {
            Expr* out = invtrig_branch_reduce(ctx, h, a0->data.function.args[0]);
            if (out) { expr_free(node); *changed = 1; return out; }
        }

        /* Reduce-verified inverse-trig multiple-angle: ArcCos[1-2x^2]->2ArcSin[x],
         * ArcCos[2x^2-1]->Pi-2ArcSin[x] (= 2ArcCos[x]), ArcSin[2x Sqrt[1-x^2]]->
         * 2ArcSin[x]. Applied only on a Reduce proof of the equality over the
         * region, so sound on whichever branch the region occupies. */
        if ((h == SYM_ArcSin || h == SYM_ArcCos) && ctx) {
            Expr* out = invtrig_multiangle_reduce(ctx, h, a0);
            if (out) { expr_free(node); *changed = 1; return out; }
        }

        /* ArcCosh[u] -> Log[u + Sqrt[u^2 - 1]] when u >= 1. This is the principal
         * -value identity, exact on [1, inf) but UNSOUND for u <= -1 (the combined
         * radical Sqrt[u^2-1] sits on the wrong branch there), so the u >= 1 gate
         * is essential. It exposes exactly the combined-radical Log form a user
         * writes, so ArcCosh[u] - Log[u + Sqrt[u^2-1]] cancels. (ArcCosh[Cosh[e]]
         * is handled sign-wise above; this fires for any other argument.) */
        if (h == SYM_ArcCosh && !sr_head(a0, "Cosh")) {
            const Expr* u = a0;
            Expr* um1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_copy((Expr*)u), expr_new_integer(-1) }, 2));
            bool ge1 = assume_known_nonneg(ctx, um1) || (assume_reduce_nonneg(ctx, um1) == 1);
            expr_free(um1);
            if (ge1) {
                Expr* u2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy((Expr*)u), expr_new_integer(2) }, 2));
                Expr* u2m1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                               (Expr*[]){ u2, expr_new_integer(-1) }, 2));
                Expr* rad = eval_and_free(expr_new_function(expr_new_symbol(SYM_Sqrt),
                               (Expr*[]){ u2m1 }, 1));
                Expr* sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                               (Expr*[]){ expr_copy((Expr*)u), rad }, 2));
                Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Log),
                               (Expr*[]){ sum }, 1));
                expr_free(node); *changed = 1; return out;
            }
        }
    }

    /* Power[Power[E, w], r] -> Power[E, r w] when Im[w] lies in (-Pi, Pi], the
     * strip on which the principal Log recovers w, so (E^w)^r = E^(r w) for ANY
     * exponent r (Sqrt included). This is the sound, branch-correct counterpart
     * of PowerExpand's unconditional strip-ignoring collapse: e.g.
     * Sqrt[E^(2 I x)] -> E^(I x) under -Pi/2 < x < Pi/2. The structural guard
     * (base is a literal E-power) keeps the Reduce entailment call rare. */
    if (h == SYM_Power && nn == 2 && ctx) {
        const Expr* base = node->data.function.args[0];
        const Expr* r = node->data.function.args[1];
        if (base->type == EXPR_FUNCTION && base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL &&
            base->data.function.head->data.symbol.name == SYM_Power &&
            base->data.function.arg_count == 2 &&
            base->data.function.args[0]->type == EXPR_SYMBOL &&
            base->data.function.args[0]->data.symbol.name == SYM_E) {
            const Expr* w = base->data.function.args[1];
            Expr* imw = expr_new_function(expr_new_symbol(SYM_Im),
                            (Expr*[]){ expr_copy((Expr*)w) }, 1);
            Expr* negpi = eval_and_free(parse_expression("-Pi"));
            Expr* pi = eval_and_free(parse_expression("Pi"));
            Expr* lo = expr_new_function(expr_new_symbol(SYM_Less),
                           (Expr*[]){ negpi, expr_copy(imw) }, 2);
            Expr* hi = expr_new_function(expr_new_symbol(SYM_LessEqual),
                           (Expr*[]){ imw, pi }, 2);
            Expr* pred = expr_new_function(expr_new_symbol(SYM_And),
                             (Expr*[]){ lo, hi }, 2);
            bool instrip = (assume_reduce_entails(ctx, pred) == 1);
            expr_free(pred);
            if (instrip) {
                Expr* rw = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                               (Expr*[]){ expr_copy((Expr*)r), expr_copy((Expr*)w) }, 2));
                Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_new_symbol(SYM_E), rw }, 2));
                expr_free(node); *changed = 1; return out;
            }
        }
    }

    /* Power[Power[b, p], q] -> Power[b, p q] when b >= 0 AND the inner exponent p
     * is provably REAL. For a nonnegative real base, b^p is a nonnegative real,
     * so the outer power meets no branch cut and (b^p)^q = b^(p q). The p-real
     * requirement is essential and was the bug in the former string rules, which
     * folded for ANY inner exponent: for complex p, Im[p Log b] can leave
     * (-Pi, Pi], the principal Log of b^p no longer equals p Log b, and
     * (b^p)^q != b^(p q) (e.g. a > 0 alone does NOT license (a^p)^q -> a^(p q)).
     * The outer-integer case (sound for any base/inner exponent) is handled
     * separately in simp_power.c. */
    if (h == SYM_Power && nn == 2 && ctx) {
        const Expr* base = node->data.function.args[0];
        const Expr* q = node->data.function.args[1];
        if (base->type == EXPR_FUNCTION && base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL &&
            base->data.function.head->data.symbol.name == SYM_Power &&
            base->data.function.arg_count == 2) {
            const Expr* b = base->data.function.args[0];
            const Expr* p = base->data.function.args[1];
            if (prov_nn(ctx, b) && prov_re(ctx, p)) {
                Expr* pq = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                               (Expr*[]){ expr_copy((Expr*)p), expr_copy((Expr*)q) }, 2));
                Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy((Expr*)b), pq }, 2));
                expr_free(node); *changed = 1; return out;
            }
        }
    }

    /* Root-of-unity power: (-1)^E, I^E, (-I)^E, E^(I Pi X) under integer
     * assumptions, reduced by canon_neg1_pow. All route to (-1)^Y. */
    if (h == SYM_Power && nn == 2 && ctx) {
        const Expr* base = node->data.function.args[0];
        const Expr* ex   = node->data.function.args[1];
        Expr* Y = NULL;
        if (base->type == EXPR_INTEGER && base->data.integer == -1 &&
            prov_int(ctx, ex)) {
            Y = expr_copy((Expr*)ex);                    /* (-1)^ex */
        } else if (is_imag_unit(base, 1) && prov_int(ctx, ex)) {
            /* I^ex = (-1)^(ex/2) */
            Y = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy((Expr*)ex),
                               expr_new_function(expr_new_symbol(SYM_Rational),
                                   (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2));
        } else if (is_imag_unit(base, -1) && prov_int(ctx, ex)) {
            /* (-I)^ex = (-1)^(-ex/2) */
            Y = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy((Expr*)ex),
                               expr_new_function(expr_new_symbol(SYM_Rational),
                                   (Expr*[]){ expr_new_integer(-1), expr_new_integer(2) }, 2) }, 2));
        } else if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E) {
            Y = exp_ipi_coeff(ex);                       /* E^(I Pi X) = (-1)^X */
        }
        if (Y) {
            Expr* out = canon_neg1_pow(ctx, Y);
            expr_free(Y);
            if (out && !expr_eq(out, node)) { expr_free(node); *changed = 1; return out; }
            if (out) expr_free(out);
        }
    }

    /* Mod of a provably-integer argument under integer assumptions:
     *   Mod[c, m] -> 0         when m | c  (c/m is a provable integer),
     *   Mod[P, 2] -> Mod[R, 2] where R is the Boolean-ring reduction of P
     *                          (parity-preserving; e.g. Mod[n^2,2] -> Mod[n,2],
     *                           Mod[-n,2] -> Mod[n,2]). */
    if (h == SYM_Mod && nn == 2 && ctx) {
        Expr* c = node->data.function.args[0];
        Expr* m = node->data.function.args[1];
        if (prov_int(ctx, c) && prov_int(ctx, m)) {
            /* divisibility: c/m a provable integer => Mod == 0 */
            Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_copy(c),
                                     expr_new_function(expr_new_symbol(SYM_Power),
                                         (Expr*[]){ expr_copy(m), expr_new_integer(-1) }, 2) }, 2));
            bool divides = prov_int(ctx, q);
            expr_free(q);
            if (divides) { expr_free(node); *changed = 1; return expr_new_integer(0); }
            /* parity-preserving reduction of the Mod-2 argument */
            if (m->type == EXPR_INTEGER && m->data.integer == 2) {
                Expr* R = simp_reduce_mod2(ctx, c);
                if (R && !expr_eq(R, c)) {
                    Expr* out = expr_new_function(expr_new_symbol("Mod"),
                                    (Expr*[]){ R, expr_new_integer(2) }, 2);
                    expr_free(node); *changed = 1; return eval_and_free(out);
                }
                if (R) expr_free(R);
            }
        }
    }

    /* Mod[a, m] -> r from an Element[(a + c)/m, Integers] fact. */
    if (h == SYM_Mod && nn == 2 && ctx) {
        Expr* a = node->data.function.args[0];
        Expr* m = node->data.function.args[1];
        if (m->type == EXPR_INTEGER && m->data.integer > 0) {
            long mm = m->data.integer;
            for (size_t i = 0; i < ctx->count; i++) {
                const Expr* f = ctx->facts[i];
                if (!fact_is_function(f, "Element", 2)) continue;
                const Expr* E = f->data.function.args[0];
                const Expr* dom = f->data.function.args[1];
                if (dom->type != EXPR_SYMBOL || dom->data.symbol.name != SYM_Integers) continue;
                /* d = m*E - a; valid iff d is a constant integer (a cancels). */
                Expr* d = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(mm), expr_copy((Expr*)E) }, 2),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), expr_copy(a) }, 2) }, 2));
                if (d->type == EXPR_INTEGER) {
                    long dv = d->data.integer;
                    long r = ((-dv) % mm + mm) % mm;
                    expr_free(d); expr_free(node); *changed = 1; return expr_new_integer(r);
                }
                expr_free(d);
            }
        }
    }

    /* LegendreP[n, -1] -> (-1)^n for integer n (P_n(1)=1 is unconditional and
     * handled in the evaluator; at x=-1 the closed form needs n integer). */
    if (h == SYM_LegendreP && nn == 2 && ctx) {
        Expr* n_e = node->data.function.args[0];
        Expr* x_e = node->data.function.args[1];
        if (x_e->type == EXPR_INTEGER && x_e->data.integer == -1 &&
            prov_int(ctx, n_e)) {
            Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_new_integer(-1), expr_copy(n_e) }, 2));
            expr_free(node); *changed = 1; return out;
        }
    }

    /* Symbolic GCD/LCM of provable integers:
     *   GCD[a,a] = LCM[a,a] = Abs[a],  GCD[a,1] = 1,  LCM[a,1] = Abs[a],
     *   GCD[a,b] = 1 when |a-b| = 1 (consecutive integers are coprime). */
    if ((h == SYM_GCD || h == SYM_LCM) && nn == 2 && ctx) {
        Expr* a = node->data.function.args[0];
        Expr* b = node->data.function.args[1];
        if (prov_int(ctx, a) && prov_int(ctx, b)) {
            bool isgcd = (h == SYM_GCD);
            bool a1 = (a->type == EXPR_INTEGER && a->data.integer == 1);
            bool b1 = (b->type == EXPR_INTEGER && b->data.integer == 1);
            Expr* out = NULL;
            if (expr_eq(a, b)) {
                out = eval_and_free(expr_new_function(expr_new_symbol("Abs"),
                          (Expr*[]){ expr_copy(a) }, 1));
            } else if (a1 || b1) {
                out = isgcd ? expr_new_integer(1)
                    : eval_and_free(expr_new_function(expr_new_symbol("Abs"),
                          (Expr*[]){ expr_copy(a1 ? b : a) }, 1));
            } else if (isgcd) {
                Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy(b), expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(a) }, 2) }, 2));
                if (diff->type == EXPR_INTEGER &&
                    (diff->data.integer == 1 || diff->data.integer == -1))
                    out = expr_new_integer(1);
                expr_free(diff);
            }
            if (out) { expr_free(node); *changed = 1; return out; }
        }
    }

    return node;
}

/* Produce a rewritten expression by applying assumption-derived rules via
 * ReplaceRepeated. Returns a newly owned expression, or NULL if no rules
 * were generated. The input is not consumed. */
/* True if `name` (an interned symbol pointer) occurs anywhere in `e`. */
static bool ar_expr_has_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (ar_expr_has_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (ar_expr_has_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* Compact a symbol bucket down to those symbols that actually appear in the
 * target expression. A per-symbol rewrite rule for a symbol absent from `input`
 * can never fire, so dropping it loses nothing -- and it keeps the O(n^2)
 * pairwise positive rules and the 8192-byte rule buffer from overflowing when
 * many symbols are assumed but few are used (e.g. Abs[a1] under 20 positivity
 * facts, which otherwise dropped ALL rules on overflow). */
static void ar_filter_to_input(char** arr, size_t* n, const Expr* input) {
    size_t w = 0;
    for (size_t i = 0; i < *n; i++)
        if (ar_expr_has_symbol(input, arr[i])) arr[w++] = arr[i];
    *n = w;
}

Expr* apply_assumption_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!ctx) return NULL;

    /* Conservative caps for the per-symbol rule synthesis. */
    enum { MAX_SYM = 16 };
    char* positives[MAX_SYM]; size_t npos;
    char* reals    [MAX_SYM]; size_t nreal;
    char* integers [MAX_SYM]; size_t nint;
    char* negatives[MAX_SYM]; size_t nneg;
    char* evens    [MAX_SYM]; size_t neven;
    collect_known_symbols(ctx, positives, &npos, reals, &nreal,
                          integers, &nint, negatives, &nneg,
                          evens, &neven, MAX_SYM);

    /* Three further buckets not gathered by collect_known_symbols:
     *   nonnegs -- x >= 0  (for (x^m)^r -> x^(m r) and Sqrt[x^2] -> x at x = 0)
     *   nonpos  -- x <= 0  (for Sqrt[x^2] -> -x and Abs[x] -> -x at x = 0)
     *   abslt1  -- -1 < b < 1  (for (a^b)^c -> a^(b c)). */
    char* nonnegs[MAX_SYM]; size_t nnn = 0;
    char* nonpos [MAX_SYM]; size_t nnp = 0;
    char* abslt1 [MAX_SYM]; size_t nab = 0;
    if (ctx) {
        Expr* neg1 = expr_new_integer(-1);
        Expr* pos1 = expr_new_integer(1);
        for (size_t i = 0; i < ctx->count; i++) {
            const Expr* f = ctx->facts[i];
            if (f->type != EXPR_FUNCTION) continue;
            for (size_t j = 0; j < f->data.function.arg_count; j++) {
                Expr* a = f->data.function.args[j];
                if (a->type != EXPR_SYMBOL) continue;
                const char* nm = a->data.symbol.name;
                if (assume_known_nonneg(ctx, a) && nnn < MAX_SYM &&
                    !sym_already_listed(nonnegs, nnn, nm)) nonnegs[nnn++] = (char*)nm;
                if (assume_known_nonpos(ctx, a) && nnp < MAX_SYM &&
                    !sym_already_listed(nonpos, nnp, nm)) nonpos[nnp++] = (char*)nm;
                if (assume_known_gt(ctx, a, neg1) && assume_known_lt(ctx, a, pos1) &&
                    nab < MAX_SYM && !sym_already_listed(abslt1, nab, nm)) abslt1[nab++] = (char*)nm;
            }
        }
        expr_free(neg1);
        expr_free(pos1);
    }

    /* Keep only symbols the target actually mentions: rules for the rest can
     * never match, and pruning them bounds the rule-string size (see the helper). */
    ar_filter_to_input(positives, &npos,  input);
    ar_filter_to_input(reals,     &nreal, input);
    ar_filter_to_input(integers,  &nint,  input);
    ar_filter_to_input(negatives, &nneg,  input);
    ar_filter_to_input(evens,     &neven, input);
    ar_filter_to_input(nonnegs,   &nnn,   input);
    ar_filter_to_input(nonpos,    &nnp,   input);
    ar_filter_to_input(abslt1,    &nab,   input);

    /* Build a single rule list "{r1, r2, ...}" as a string, then parse.  The
     * user's symbols are pasted into the text, so every pattern variable lives
     * in the AssumeRules` context: a bare r_ would capture a user symbol r in
     * the replacement (Power[Power[r, 2], r_] :> Power[r, 2 r] turned
     * Refine[1/Sqrt[r^2], r > 0] into -2). */
    char buf[32768];
    size_t off = 0;
    int wrote_any = 0;

    #define EMIT(...) do { \
        int w = snprintf(buf + off, sizeof(buf) - off, __VA_ARGS__); \
        if (w < 0 || (size_t)w >= sizeof(buf) - off) goto overflow; \
        off += (size_t)w; \
    } while (0)
    #define SEP() do { if (wrote_any) EMIT(", "); wrote_any = 1; } while (0)

    EMIT("{");

    for (size_t i = 0; i < npos; i++) {
        const char* x = positives[i];
        /* Sqrt[x^2] forms */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, -1], Rational[1, 2]] :> Power[%s, Rational[-1, 2]]", x, x);
        SEP(); EMIT("Power[Power[%s, -2], Rational[1, 2]] :> Power[%s, -1]", x, x);
        /* General (x^2)^r -> x^(2r) for x > 0 (valid for any rational r, since
         * x^2 > 0): covers 1/Sqrt[x^2] = (x^2)^(-1/2) -> 1/x and the
         * (x^2)^(3/2) -> x^3 that arise from Laurent/fractional radical powers. */
        SEP(); EMIT("Power[Power[%s, 2], AssumeRules`r_] :> Power[%s, 2 AssumeRules`r]", x, x);
        /* Sqrt[x^2 * rest] -> x * Sqrt[rest] for x > 0; lets multi-factor
         * radicals like Sqrt[x^2 y^2] reduce one symbol at a time. */
        SEP(); EMIT("Power[Times[Power[%s, 2], AssumeRules`rest___], Rational[1, 2]] :> %s Power[Times[AssumeRules`rest], Rational[1, 2]]", x, x);
        /* Abs[x] -> x  for x > 0 */
        SEP(); EMIT("Abs[%s] :> %s", x, x);
        /* Sign[x] -> 1  for x > 0 */
        SEP(); EMIT("Sign[%s] :> 1", x);
        /* Arg[x] -> 0  for x > 0 (positive real) */
        SEP(); EMIT("Arg[%s] :> 0", x);
        /* Conjugate[x] -> x  for x > 0 (x is real) */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Log[x^p] -> p Log[x]  for x > 0 (any real p; v1 accepts symbolic p too) */
        SEP(); EMIT("Log[Power[%s, AssumeRules`p_]] :> AssumeRules`p Log[%s]", x, x);
        /* Log[x rest] -> Log[x] + Log[rest]  for x > 0: a positive real factor
         * never shifts the branch (arg[x rest] = arg[rest]), so it splits off
         * exactly. Peels one positive factor at a time; Log[x] (bare) does not
         * re-match this Times pattern, so it terminates. */
        SEP(); EMIT("Log[Times[%s, AssumeRules`rest___]] :> Log[%s] + Log[Times[AssumeRules`rest]]", x, x);
        /* Inverse-trig sum identity: ArcTan[x] + ArcTan[1/x] -> Pi/2  for x > 0.
         * Mathilda's matcher does NOT perform orderless-Plus subset matching
         * out of the box (unlike Mathematica), so the rule must explicitly
         * absorb the trailing terms via `+ rest___` and re-emit them. The
         * BlankNullSequence pattern matches 0 or more remaining terms, so
         * the bare two-term sum reduces too. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + AssumeRules`rest___ :> Pi/2 + AssumeRules`rest", x, x);
        /* ArcCosh double-angle reduction: ArcCosh[2x^2 - 1] -> 2 ArcCosh[x]
         * for x > 0. The identity holds for any complex x with the
         * principal branch (verified for x in [0,1] via the i*ArcCos bridge
         * and for x >= 1 directly), so the x > 0 condition is the weakest
         * sufficient assumption. The Plus arg list is canonical
         * (Plus[-1, Times[2, x^2]]). */
        SEP(); EMIT("ArcCosh[Plus[-1, Times[2, Power[%s, 2]]]] :> 2 ArcCosh[%s]", x, x);
        /* NOTE: (x^p)^q -> x^(p q) for x > 0 is handled structurally in
         * assume_structural_rewrite, which can require the inner exponent p to
         * be provably REAL (a string rule here cannot see the assumptions and
         * would fold for complex p too -- an unsound result). */
        /* x^p y^p -> (x y)^p for two positive bases x, y (any common
         * exponent p). Emitted once per unordered pair. */
        for (size_t j = i + 1; j < npos; j++) {
            const char* y = positives[j];
            SEP(); EMIT("Times[Power[%s, AssumeRules`p_], Power[%s, AssumeRules`p_], AssumeRules`rest___] :> "
                        "Power[Times[%s, %s], AssumeRules`p] AssumeRules`rest", x, y, x, y);
        }
    }

    for (size_t i = 0; i < nneg; i++) {
        const char* x = negatives[i];
        /* Abs[x] -> -x  for x < 0 */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        /* Sign[x] -> -1  for x < 0 */
        SEP(); EMIT("Sign[%s] :> -1", x);
        /* Arg[x] -> Pi  for x < 0 (negative real, principal branch) */
        SEP(); EMIT("Arg[%s] :> Pi", x);
        /* Conjugate[x] -> x  for x < 0 (x is real) */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Sqrt[x^2] -> -x  for x < 0  (Power[Power[x,2], 1/2]) */
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        /* General (x^2)^r -> (-x)^(2r) for x < 0 (x^2 = (-x)^2, -x > 0):
         * covers 1/Sqrt[x^2] = (x^2)^(-1/2) -> 1/(-x) and (x^2)^(3/2) -> -x^3. */
        SEP(); EMIT("Power[Power[%s, 2], AssumeRules`r_] :> Power[Times[-1, %s], 2 AssumeRules`r]", x, x);
        /* Sqrt[x^2 * rest] -> -x * Sqrt[rest] for x < 0. */
        SEP(); EMIT("Power[Times[Power[%s, 2], AssumeRules`rest___], Rational[1, 2]] :> -%s Power[Times[AssumeRules`rest], Rational[1, 2]]", x, x);
        /* Mirror of the x > 0 inverse-trig sum: ArcTan[x] + ArcTan[1/x]
         * -> -Pi/2 for x < 0. Same `+ rest___` trick as above to handle
         * the embedded-in-larger-sum case. */
        SEP(); EMIT("ArcTan[%s] + ArcTan[Power[%s, -1]] + AssumeRules`rest___ :> -Pi/2 + AssumeRules`rest", x, x);
        /* Log[x] -> I Pi + Log[-x]  for x < 0 (principal branch). */
        SEP(); EMIT("Log[%s] :> I Pi + Log[-%s]", x, x);
        /* Half-integer power of a negative real, principal branch:
         *   x^(p/2) = (-1)^(p/2) (-x)^(p/2) = I^p (-x)^(p/2)  for x < 0.
         * Covers Sqrt[x] -> I Sqrt[-x], 1/Sqrt[x] -> -I/Sqrt[-x],
         * x^(3/2) -> -I (-x)^(3/2). The RHS base Times[-1,x] cannot re-match
         * Power[x, _] and Power[I, p] folds to a concrete Complex, so it
         * terminates. The multiplicative analogue of the Log rule above. */
        SEP(); EMIT("Power[%s, Rational[AssumeRules`p_, 2]] :> Power[I, AssumeRules`p] Power[Times[-1, %s], Rational[AssumeRules`p, 2]]", x, x);
    }

    /* Non-positive (x <= 0, but not proven strictly negative): Abs[x] -> -x,
     * Sqrt[x^2] -> -x, (x^2)^r -> (-x)^(2r) all hold at x <= 0 including x = 0
     * (where -x = 0 = |x|). Sign/Arg are NOT constant on x <= 0 (0 at x = 0),
     * so they are deliberately omitted. Skip symbols already strictly negative. */
    for (size_t i = 0; i < nnp; i++) {
        const char* x = nonpos[i];
        if (sym_already_listed(negatives, nneg, x)) continue;
        if (sym_already_listed(positives, npos, x)) continue;   /* 0 only: handled by nonneg too */
        SEP(); EMIT("Abs[%s] :> -%s", x, x);
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> -%s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], AssumeRules`r_] :> Power[Times[-1, %s], 2 AssumeRules`r]", x, x);
        SEP(); EMIT("Power[Times[Power[%s, 2], AssumeRules`rest___], Rational[1, 2]] :> -%s Power[Times[AssumeRules`rest], Rational[1, 2]]", x, x);
    }

    /* For real-but-unknown-sign, Sqrt[x^2] -> Abs[x]. Skip symbols already
     * proven positive, negative, or nonpositive (their stronger rule above wins). */
    for (size_t i = 0; i < nreal; i++) {
        const char* x = reals[i];
        if (sym_already_listed(positives, npos, x)) continue;
        if (sym_already_listed(negatives, nneg, x)) continue;
        if (sym_already_listed(nonpos, nnp, x)) continue;
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> Abs[%s]", x, x);
        /* Sqrt[x^2 * rest] -> Abs[x] * Sqrt[rest] for real x. */
        SEP(); EMIT("Power[Times[Power[%s, 2], AssumeRules`rest___], Rational[1, 2]] :> Abs[%s] Power[Times[AssumeRules`rest], Rational[1, 2]]", x, x);
        /* Conjugate[x] -> x  for real x. */
        SEP(); EMIT("Conjugate[%s] :> %s", x, x);
        /* Log[x^2] -> 2 Log[Abs[x]] for real x (|x| keeps the argument real
         * and positive across the sign of x). */
        SEP(); EMIT("Log[Power[%s, 2]] :> 2 Log[Abs[%s]]", x, x);
        /* Log[E^x] -> x for real x (E^x is a positive real). */
        SEP(); EMIT("Log[Power[E, %s]] :> %s", x, x);
    }

    /* Sin[n Pi] -> 0, Cos[n Pi] -> (-1)^n, Tan[n Pi] -> 0 for integer n.
     * Plus: (-1)^(even_int * n) -> 1 and ((-1)^n)^even_int -> 1, so the
     * Cos rule can collapse all the way (the standalone Cos result
     * Cos[k Pi]^4 -> Power[-1, 4 k], for instance). */
    for (size_t i = 0; i < nint; i++) {
        const char* n = integers[i];
        SEP(); EMIT("Sin[%s Pi] :> 0", n);
        SEP(); EMIT("Sin[Pi %s] :> 0", n);
        SEP(); EMIT("Cos[%s Pi] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Cos[Pi %s] :> Power[-1, %s]", n, n);
        SEP(); EMIT("Tan[%s Pi] :> 0", n);
        SEP(); EMIT("Tan[Pi %s] :> 0", n);
        SEP(); EMIT("Power[-1, Times[AssumeRules`m_Integer /; EvenQ[AssumeRules`m], %s]] :> 1", n);
        SEP(); EMIT("Power[Power[-1, %s], AssumeRules`m_Integer /; EvenQ[AssumeRules`m]] :> 1", n);
        /* (-1)^(k + n) -> (-1)^k (-1)^n for integer n and any k: an integer
         * additive shift in a (-1) exponent splits out so the sign factor can
         * fold and cancel (e.g. (-1)^n + (-1)^(n+1) -> 0). Sound for integer k
         * since (-1)^(k+r) = (-1)^k (-1)^r exactly. */
        SEP(); EMIT("Power[-1, Plus[AssumeRules`k_Integer, %s]] :> Power[-1, AssumeRules`k] Power[-1, %s]", n, n);
        /* Shift identities: Cos[x + n Pi] -> (-1)^n Cos[x], likewise Sin, and
         * Tan[x + n Pi] -> Tan[x] (period Pi), for integer n. The rest___
         * absorbs the remaining Plus terms (Plus[] -> 0 gives the bare
         * Cos[n Pi] case too). */
        SEP(); EMIT("Cos[Plus[Times[%s, Pi], AssumeRules`rest___]] :> Power[-1, %s] Cos[Plus[AssumeRules`rest]]", n, n);
        SEP(); EMIT("Sin[Plus[Times[%s, Pi], AssumeRules`rest___]] :> Power[-1, %s] Sin[Plus[AssumeRules`rest]]", n, n);
        SEP(); EMIT("Tan[Plus[Times[%s, Pi], AssumeRules`rest___]] :> Tan[Plus[AssumeRules`rest]]", n);
        /* Exp periodicity for integer n: e^(2 I Pi n + rest) -> e^rest, and the
         * odd-multiple e^(I Pi n + rest) -> (-1)^n e^rest. The evaluator keeps
         * Exp[sum] fused and does not know n is an integer, so this peel must be
         * assumption-gated here. */
        SEP(); EMIT("Power[E, Plus[Times[2 I, Pi, %s], AssumeRules`rest___]] :> Power[E, Plus[AssumeRules`rest]]", n);
        SEP(); EMIT("Power[E, Plus[Times[I, Pi, %s], AssumeRules`rest___]] :> Power[-1, %s] Power[E, Plus[AssumeRules`rest]]", n, n);
    }

    /* Even-exponent identities: (-1)^m = 1 when m is even, and the
     * lifted forms ((-1)^k)^m and (-1)^(k m) when additionally k is an
     * integer (so k m is also even). The pair-wise integer/even rules
     * cover the common Cos[k Pi]^m -> 1 path -- the existing Cos rule
     * above rewrites Cos[k Pi] to Power[-1, k], so the Power surface form
     * we land on is Power[Power[-1, k], m]. */
    for (size_t i = 0; i < neven; i++) {
        const char* m = evens[i];
        SEP(); EMIT("Power[-1, %s] :> 1", m);
        /* Literal-integer multiplier (handles concrete numbers without
         * needing them in the integer set). */
        SEP(); EMIT("Power[-1, AssumeRules`k_Integer %s] :> 1", m);
        SEP(); EMIT("Power[Power[-1, AssumeRules`k_Integer], %s] :> 1", m);
        for (size_t j = 0; j < nint; j++) {
            const char* k = integers[j];
            SEP(); EMIT("Power[-1, %s %s] :> 1", k, m);
            SEP(); EMIT("Power[Power[-1, %s], %s] :> 1", k, m);
        }
    }

    /* Non-negative (but not proven strictly positive) symbols: (x^m)^r ->
     * x^(m r), Sqrt[x^2] -> x, Abs[x] -> x all hold at x >= 0 including x = 0.
     * Strictly-positive symbols already received the stronger rule set above. */
    for (size_t i = 0; i < nnn; i++) {
        const char* x = nonnegs[i];
        if (sym_already_listed(positives, npos, x)) continue;
        if (sym_already_listed(negatives, nneg, x)) continue;
        SEP(); EMIT("Abs[%s] :> %s", x, x);
        SEP(); EMIT("Power[Power[%s, 2], Rational[1, 2]] :> %s", x, x);
        /* (x^p)^q -> x^(p q) handled structurally (needs p real); see above. */
        SEP(); EMIT("Power[Times[Power[%s, 2], AssumeRules`rest___], Rational[1, 2]] :> %s Power[Times[AssumeRules`rest], Rational[1, 2]]", x, x);
    }

    /* Symbols with -1 < b < 1: (a^b)^c -> a^(b c) on the principal branch,
     * for any base a. */
    for (size_t i = 0; i < nab; i++) {
        const char* b = abslt1[i];
        SEP(); EMIT("Power[Power[AssumeRules`base_, %s], AssumeRules`c_] :> Power[AssumeRules`base, %s AssumeRules`c]", b, b);
    }

    /* Equal[u, v] facts -> two-way substitution rules. We use immediate
     * Rule (->) so the pattern uses exact structural matching. */
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        /* We can't easily re-emit arbitrary Expr* into our string buffer;
         * instead, build these rules as Expr* and merge them in below. */
        (void)f; /* handled in the Expr* merge step below */
    }

    EMIT("}");

    if (!wrote_any) {
        /* No string-built rules. We may still have Equal substitutions. */
    }

    Expr* string_rules = wrote_any ? parse_expression(buf) : NULL;

    /* Now build Equal-substitution rules.
     *
     * Two complementary rules per equation:
     *   1. The direct rule heavier(lhs,rhs) -> lighter (catches cases
     *      where the equation's LHS appears verbatim as a subterm).
     *   2. ONE monomial-isolation rule when diff = lhs - rhs is a Plus
     *      with >= 3 terms: pick the heaviest non-numeric term t and
     *      emit t -> -(other terms). Polynomial relations like
     *      a^2 + b^2 == 1 then rewrite occurrences of a^2 even when
     *      the full "a^2 + b^2" sum is not present in the input.
     *
     * Emitting only one monomial rule (instead of one per term) avoids
     * the bidirectional cycle a^2 -> 1-b^2 ; b^2 -> 1-a^2 that
     * ReplaceRepeated would chase up to its 65536 iteration cap. */
    Expr** eq_diffs = (Expr**)calloc(ctx->count, sizeof(Expr*));
    size_t eq_count = 0;
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* sub_args[2] = { expr_copy(lhs),
                              expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_new_integer(-1), expr_copy(rhs) }, 2) };
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), sub_args, 2);
        Expr* diff = eval_and_free(sum);
        eq_diffs[i] = diff;
        /* Always emit the direct heavier->lighter rule. */
        eq_count++;
        /* Plus extra monomial-isolation rule for polynomial relations. */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol.name == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            for (size_t j = 0; j < diff->data.function.arg_count; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                eq_count++;
                break; /* one monomial rule per equation */
            }
        }
    }

    /* Apply the synthesized string + Equal-substitution rules (if any) to
     * obtain base_result; the structural pass below runs regardless. */
    Expr* base_result = NULL;
    if (string_rules || eq_count > 0) {
    size_t string_len = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        string_len = string_rules->data.function.arg_count;
    }
    size_t total = string_len + eq_count;
    Expr** all = (Expr**)calloc(total, sizeof(Expr*));
    size_t fill = 0;
    if (string_rules && string_rules->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < string_len; i++) {
            all[fill++] = expr_copy(string_rules->data.function.args[i]);
        }
    }
    for (size_t i = 0; i < ctx->count; i++) {
        const Expr* f = ctx->facts[i];
        if (!fact_is_function(f, "Equal", 2)) continue;
        Expr* lhs = f->data.function.args[0];
        Expr* rhs = f->data.function.args[1];
        Expr* diff = eq_diffs[i];

        /* Direct heavier->lighter rule. */
        Expr *src, *dst;
        if (simp_default_complexity(lhs) >= simp_default_complexity(rhs)) {
            src = lhs; dst = rhs;
        } else {
            src = rhs; dst = lhs;
        }
        Expr* direct[2] = { expr_copy(src), expr_copy(dst) };
        all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), direct, 2);

        /* Polynomial-relation monomial-isolation rule (one per fact). */
        if (diff->type == EXPR_FUNCTION &&
            diff->data.function.head &&
            diff->data.function.head->type == EXPR_SYMBOL &&
            diff->data.function.head->data.symbol.name == SYM_Plus &&
            diff->data.function.arg_count >= 3) {
            size_t n = diff->data.function.arg_count;
            /* Pick the first non-numeric term, breaking ties by canonical
             * (Plus-Orderless) order which is already applied by the
             * evaluator. */
            size_t pick = (size_t)-1;
            size_t pick_score = 0;
            for (size_t j = 0; j < n; j++) {
                Expr* term = diff->data.function.args[j];
                if (term->type == EXPR_INTEGER || term->type == EXPR_BIGINT ||
                    term->type == EXPR_REAL) continue;
                size_t s = simp_default_complexity(term);
                if (pick == (size_t)-1 || s > pick_score) {
                    pick = j;
                    pick_score = s;
                }
            }
            if (pick != (size_t)-1) {
                Expr* term = diff->data.function.args[pick];
                Expr** other_args = (Expr**)calloc(n - 1, sizeof(Expr*));
                size_t oi = 0;
                for (size_t k = 0; k < n; k++) {
                    if (k == pick) continue;
                    other_args[oi++] = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1),
                                   expr_copy(diff->data.function.args[k]) }, 2);
                }
                Expr* iso_rhs;
                if (n - 1 == 1) {
                    iso_rhs = other_args[0];
                    free(other_args);
                } else {
                    iso_rhs = expr_new_function(expr_new_symbol(SYM_Plus), other_args, n - 1);
                    free(other_args);
                }
                Expr* iso[2] = { expr_copy(term), iso_rhs };
                all[fill++] = expr_new_function(expr_new_symbol(SYM_Rule), iso, 2);
            }
        }
    }
    for (size_t i = 0; i < ctx->count; i++) if (eq_diffs[i]) expr_free(eq_diffs[i]);
    free(eq_diffs);
    if (string_rules) expr_free(string_rules);

    Expr* rules_list = expr_new_function(expr_new_symbol(SYM_List), all, fill);
    free(all);

    Expr* call_args[2] = { expr_copy((Expr*)input), rules_list };
    Expr* call = expr_new_function(expr_new_symbol(SYM_ReplaceRepeated), call_args, 2);
    Expr* out = evaluate(call);
    expr_free(call);
    base_result = out;
    } else {
        free(eq_diffs);   /* eq_count == 0: no entries were populated */
    }

    /* Structural assumption rewrites (Floor/Ceiling/Round/IntegerPart/
     * FractionalPart/Mod, Re/Im/Conjugate/Arg, ArcTan[Tan]) run whether or not
     * any string/Equal rule was generated, so a lone structural fact still
     * fires and Simplify inherits the same rewrites. */
    {
        const Expr* walk_in = base_result ? base_result : input;
        int schanged = 0;
        Expr* walked = assume_structural_rewrite(walk_in, ctx, &schanged);
        if (schanged && walked) {
            Expr* ev = eval_and_free(walked);
            if (base_result) expr_free(base_result);
            return ev;
        }
        if (walked) expr_free(walked);
    }
    return base_result;

overflow:
    /* Buffer was too small; bail out, no rules applied. */
    return NULL;

    #undef EMIT
    #undef SEP
}

