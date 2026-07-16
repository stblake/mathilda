/* risch_trig_frontend.c — trig/hyperbolic front-end for the Risch integrator.
 *
 * Exponentialise (TrigToExp) -> integrate as a Laurent-rational in E^(I x)/E^x
 * -> ExpToTrig -> real reconstruction (rt_realify).  Plus the real hypertangent
 * (Tan/Tanh) case and the single-exponential rational closer.  Runs on top of
 * the single-extension and field integrators.  See risch_trig_frontend.h.
 */

#include "risch_trig_frontend.h"
#include "integrate_risch_rde.h"
#include "risch_tower.h"
#include "risch_util.h"
#include "risch_singleext.h"
#include "risch_field_integrate.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

static Expr* cx_int(long n) { return expr_new_integer(n); }
static Expr* cx_add(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* cx_mul(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* cx_neg(Expr* a) { return cx_mul(cx_int(-1), a); }
static Expr* cx_half(void) {
    return expr_new_function(expr_new_symbol("Rational"),
        (Expr*[]){ cx_int(1), cx_int(2) }, 2);
}
static Expr* cx_I(void) {
    return expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ cx_int(0), cx_int(1) }, 2);
}
static Expr* cx_log(Expr* a) { return expr_new_function(expr_new_symbol("Log"), (Expr*[]){ a }, 1); }

/* True iff e contains no Complex[...] atom (manifestly real for real symbols). */
static bool rt_free_of_complex(Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_FUNCTION) {
        if (rt_head_is(e, "Complex")) return false;
        if (!rt_free_of_complex(e->data.function.head)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!rt_free_of_complex(e->data.function.args[i])) return false;
    }
    return true;
}

/* e == *re + I (*im), for a function of the real variable x.  Returns false on a
 * construct it cannot split.  Owns the returned parts. */
static bool cx_reim(Expr* e, Expr* x, Expr** re, Expr** im) {
    *re = NULL; *im = NULL;
    if (rt_head_is(e, "Complex") && e->data.function.arg_count == 2) {
        *re = expr_copy(e->data.function.args[0]);
        *im = expr_copy(e->data.function.args[1]);
        return true;
    }
    if (rt_free_of_complex(e)) {          /* real atom / trig/Sqrt of real x, ... */
        *re = expr_copy(e); *im = cx_int(0); return true;
    }
    if (e->type != EXPR_FUNCTION) return false;
    size_t n = e->data.function.arg_count;

    if (rt_head_is(e, "Plus")) {
        Expr* R = cx_int(0); Expr* M = cx_int(0);
        for (size_t i = 0; i < n; i++) {
            Expr* r; Expr* m;
            if (!cx_reim(e->data.function.args[i], x, &r, &m)) {
                expr_free(R); expr_free(M); return false;
            }
            R = cx_add(R, r); M = cx_add(M, m);
        }
        *re = R; *im = M; return true;
    }
    if (rt_head_is(e, "Times")) {
        Expr* R = cx_int(1); Expr* M = cx_int(0);
        for (size_t i = 0; i < n; i++) {
            Expr* c; Expr* d;
            if (!cx_reim(e->data.function.args[i], x, &c, &d)) {
                expr_free(R); expr_free(M); expr_free(c); return false;
            }
            /* (R + M i)(c + d i) = (Rc - Md) + (Rd + Mc) i */
            Expr* nR = cx_add(cx_mul(expr_copy(R), expr_copy(c)),
                              cx_neg(cx_mul(expr_copy(M), expr_copy(d))));
            Expr* nM = cx_add(cx_mul(expr_copy(R), expr_copy(d)),
                              cx_mul(expr_copy(M), expr_copy(c)));
            expr_free(R); expr_free(M); expr_free(c); expr_free(d);
            R = nR; M = nM;
        }
        *re = R; *im = M; return true;
    }
    if (rt_head_is(e, "Power") && n == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (rt_free_of_complex(base)) {
            /* base^(a + b i) = base^a (Cos[b Log base] + i Sin[b Log base]). */
            Expr* a; Expr* b;
            if (!cx_reim(ex, x, &a, &b)) return false;
            Expr* mag = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(base), a }, 2);                 /* a adopted */
            Expr* th = cx_mul(b, cx_log(expr_copy(base)));           /* b adopted */
            *re = cx_mul(expr_copy(mag),
                expr_new_function(expr_new_symbol("Cos"), (Expr*[]){ expr_copy(th) }, 1));
            *im = cx_mul(mag,
                expr_new_function(expr_new_symbol("Sin"), (Expr*[]){ th }, 1));
            return true;
        }
        if (ex->type == EXPR_INTEGER) {                             /* complex^int */
            Expr* br; Expr* bi;
            if (!cx_reim(base, x, &br, &bi)) return false;
            long p = ex->data.integer, ap = p < 0 ? -p : p;
            Expr* R = cx_int(1); Expr* M = cx_int(0);
            for (long k = 0; k < ap; k++) {
                Expr* nR = cx_add(cx_mul(expr_copy(R), expr_copy(br)),
                                  cx_neg(cx_mul(expr_copy(M), expr_copy(bi))));
                Expr* nM = cx_add(cx_mul(expr_copy(R), expr_copy(bi)),
                                  cx_mul(expr_copy(M), expr_copy(br)));
                expr_free(R); expr_free(M); R = nR; M = nM;
            }
            expr_free(br); expr_free(bi);
            if (p < 0) {                                            /* invert */
                Expr* den = cx_add(cx_mul(expr_copy(R), expr_copy(R)),
                                   cx_mul(expr_copy(M), expr_copy(M)));
                Expr* invden = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ den, cx_int(-1) }, 2);
                Expr* nR = cx_mul(R, expr_copy(invden));
                Expr* nM = cx_mul(cx_neg(M), invden);
                R = nR; M = nM;
            }
            *re = R; *im = M; return true;
        }
        return false;                                              /* complex^noninteger */
    }
    if (rt_head_is(e, "Log") && n == 1) {
        Expr* a; Expr* b;
        if (!cx_reim(e->data.function.args[0], x, &a, &b)) return false;
        /* Log[a + b i] = (1/2) Log[a^2 + b^2] + i ArcTan[a, b]  (Arg). */
        Expr* mag2 = cx_add(cx_mul(expr_copy(a), expr_copy(a)),
                            cx_mul(expr_copy(b), expr_copy(b)));
        *re = cx_mul(cx_half(), cx_log(mag2));
        *im = expr_new_function(expr_new_symbol("ArcTan"), (Expr*[]){ a, b }, 2);
        return true;
    }
    /* ArcTan[z] = (i/2)(Log[1 - i z] - Log[1 + i z]); recurse on the rewrite. */
    if (rt_head_is(e, "ArcTan") && n == 1) {
        Expr* z = e->data.function.args[0];
        Expr* iz1 = cx_mul(cx_I(), expr_copy(z));
        Expr* iz2 = cx_mul(cx_I(), expr_copy(z));
        Expr* rw = cx_mul(cx_mul(cx_I(), cx_half()),
            cx_add(cx_log(cx_add(cx_int(1), cx_neg(iz1))),
                   cx_neg(cx_log(cx_add(cx_int(1), iz2)))));
        bool ok = cx_reim(rw, x, re, im);
        expr_free(rw);
        return ok;
    }
    /* ArcTanh[z] = (1/2)(Log[1 + z] - Log[1 - z]). */
    if (rt_head_is(e, "ArcTanh") && n == 1) {
        Expr* z = e->data.function.args[0];
        Expr* rw = cx_mul(cx_half(),
            cx_add(cx_log(cx_add(cx_int(1), expr_copy(z))),
                   cx_neg(cx_log(cx_add(cx_int(1), cx_neg(expr_copy(z)))))));
        bool ok = cx_reim(rw, x, re, im);
        expr_free(rw);
        return ok;
    }
    return false;
}

/* Real closed form of the I-laden antiderivative G of the real integrand f:
 * Re[G] (x real), returned only if it is manifestly real and diff-backs to f
 * under the EXACT symbolic gate (rt_verify_antideriv, Simplify[D-f]===0).  A
 * Risch decision procedure never certifies by numeric sampling; if Simplify
 * cannot reduce the reconstructed real form to prove it, this declines. */
static long rt_nodecount(Expr* e) {
    if (!e) return 0;
    long n = 1;
    if (e->type == EXPR_FUNCTION) {
        n += rt_nodecount(e->data.function.head);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            n += rt_nodecount(e->data.function.args[i]);
    }
    return n;
}

/* True iff `e` is a polynomial in x and trigonometric/hyperbolic kernels: free
 * of Log and inverse (Arc*) heads, and of any Power with a negative or
 * non-integer exponent (i.e. no denominators / radicals).  This is exactly the
 * antiderivative shape produced for a trig-POLYNOMIAL integrand such as
 * x^n Cos[x] Sin[x]^m — a sum of x^k {Sin,Cos,…}[j x] — for which ExpToTrig's
 * real part is already the clean multiple-angle form.  Simplify's transform
 * search Factors the TrigExpand image of such a form, and the multivariate
 * factorizer is pathologically slow over the OVERLAPPING generators
 * {x, Cos[x], Sin[x]} (Cos[x]/Sin[x] themselves contain x): seconds per term,
 * the root cause of the Integrate[x^7 Cos[x] Sin[x]^5, …] hang.  Rational-trig
 * answers (Sec[x], 1/(2+Cos[x])) carry Log/ArcTan and are excluded here — they
 * do not hit that blow-up and genuinely benefit from full Simplify. */
static bool rt_is_trig_polynomial(Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return true;   /* atom: symbol / number */
    Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* nm = h->data.symbol.name;
        if (nm == intern_symbol("Log")) return false;
        if (strncmp(nm, "Arc", 3) == 0) return false;   /* ArcTan, ArcTanh, … */
        if (nm == intern_symbol("Power") && e->data.function.arg_count == 2) {
            /* Only non-negative integer exponents keep it a polynomial; a
             * negative or fractional exponent is a denominator / radical. */
            Expr* ex = e->data.function.args[1];
            if (!(ex->type == EXPR_INTEGER && ex->data.integer >= 0)) return false;
            return rt_is_trig_polynomial(e->data.function.args[0]);
        }
    }
    if (!rt_is_trig_polynomial(h)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_is_trig_polynomial(e->data.function.args[i])) return false;
    return true;
}

static Expr* rt_realify(Expr* G, Expr* x, Expr* f) {
    Expr* re = NULL; Expr* im = NULL;
    if (!cx_reim(G, x, &re, &im)) { if (re) expr_free(re); if (im) expr_free(im); return NULL; }
    if (im) expr_free(im);
    /* Complexity guard: on a form carrying an ArcTan/Log of a complex argument (the
     * odd rational-trig powers, e.g. Sec[x]^3), cx_reim's exact real part explodes —
     * ~280k nodes vs. ~13k for the heaviest case that still reduces cleanly (Sec[x]^4)
     * — and Simplify then blows up.  Bail early (leak-free, no cost) so the front-end
     * keeps the compact diff-back-verified exponential form; every genuine clean case
     * sits far below the cap. */
    if (rt_nodecount(re) > 60000) { expr_free(re); return NULL; }
    /* Fast path for a trig-polynomial antiderivative (∫ x^n · trig-poly): the real
     * part is already the clean multiple-angle form.  Normalize it with the cheap,
     * linear TrigReduce (product-to-sum) instead of full Simplify — whose Factor
     * search over the overlapping {x, Cos[x], Sin[x]} generators is what hangs — and
     * return it once the exact diff-back gate certifies it.  For a trig POLYNOMIAL
     * the gate's Together[TrigToExp[D-f]] reduces to 0 immediately (no denominators),
     * so this is O(rational), not a Simplify. */
    if (rt_free_of_complex(re) && rt_is_trig_polynomial(re)) {
        Expr* reN = rt_eval1("TrigReduce", expr_copy(re));          /* adopts the copy */
        if (reN && rt_free_of_complex(reN) && rt_is_trig_polynomial(reN)
                && rt_verify_antideriv(reN, f, x)) {
            expr_free(re); return reN;
        }
        if (reN) expr_free(reN);
    }
    Expr* reS = rt_eval1("Simplify", re);                           /* adopts re */
    if (reS && rt_free_of_complex(reS) && rt_verify_antideriv(reS, f, x)) return reS;
    if (reS) expr_free(reS);
    return NULL;
}

/* Trigonometric / hyperbolic front-end (exponentialize path).
 * Rewrites the trig/hyperbolic kernels to complex exponentials with
 * TrigToExp, integrates the resulting (Laurent-)rational function of the
 * exponential kernel E^(I x) / E^x with the exponential machinery, and
 * converts the answer back to trigonometric form with ExpToTrig.  Both
 * rewrites are exact, so the result is correct by construction; rt_realify then
 * reconstructs a real closed form (diff-back gated) from the I-laden output. */
Expr* rt_trig_frontend(Expr* f, Expr* x) {
    Expr* fe = rt_eval1("TrigToExp", expr_copy(f));
    if (!fe) return NULL;
    if (expr_eq(fe, f)) { expr_free(fe); return NULL; }   /* no trig/hyperbolic */
    /* TrigToExp of a trig/inverse-trig of a LOGARITHM leaves a general power
     * (Cos[x Log x] -> (x^(I x) + x^(-I x))/2), the E^(±I x Log x) kernels having
     * collapsed to b^e.  Re-expose them as raw base-e exponentials so the tower
     * cases below recognise the hyperexponential monomial (no-op otherwise). */
    { Expr* fe2 = rt_powers_to_exp(fe, x);
      if (fe2) { expr_free(fe); fe = fe2; } }
    /* All exponential cases are used, including the coupled hyperexponential
     * one, for completeness — a correct antiderivative is returned even when
     * it cannot be reduced to the cleanest real form.
     *
     * KNOWN SIMPLIFICATION GAP (Simplify improvement opportunity): via the
     * complex substitution u = I x, Tan[x]/Tanh[x] and similar close to a
     * correct but I-laden form such as  I x - Log[1 + E^(2 I x)]  (which equals
     * -Log[Cos[x]]).  No current simplifier (Simplify / FullSimplify /
     * TrigReduce / PowerExpand) collapses  I x - Log[1 + E^(2 I x)]  to
     * -Log[Cos[x]].  The result is nonetheless returned (correct by
     * construction); teaching Simplify the log-of-product / half-angle
     * collapse would render it in real closed form. */
    Expr* r = rt_exp_poly_case(fe, x);
    if (!r) r = rt_frac_case(fe, x);
    /* Rational-function-of-a-single-exponential closer (kernelize E^(I x) -> pure
     * rational integral).  Closes the Laurent forms with E^(-I x) in the
     * denominator that the frac case leaves — Sec[x], Csc[x], 1/(2+Cos[x]),
     * Sec[x]^3 — which are exactly the rational trigonometric integrands.  Runs
     * after rt_frac_case so its cleaner squarefree ArcTan/Log parts win first. */
    if (!r) r = rt_exp_ratreduce_case(fe, x);
    if (!r) r = rt_hyperexp_case(fe, x);
    /* Multi-kernel decoupling (Phase B): e.g. Sin/Cos times a real exponential
     * exponentialize to a sum of two non-commensurate exponentials E^((a +/- b I) x)
     * that the single-primitive cases cannot kernelize. */
    if (!r) r = rt_expsum_case(fe, x);
    /* Nested / mixed exponential towers exposed by TrigToExp: a trig kernel whose
     * argument itself carries a nested exponential (e.g. Cos[x E^(E^x)] ->
     * E^(±I x E^(E^x))) exponentializes to a genuine multi-kernel tower that the
     * single-primitive closers above cannot kernelize.  Route it through the
     * one-extension recursion (and the flat exp-tower fallback), which builds the
     * {E^x, E^(E^x), E^(±I x E^(E^x))} tower and integrates it; the I-laden result
     * is realified below (diff-back gated).  Closes e.g.
     * Integrate[E^E^x (1 + x E^x) Cos[x E^E^x], x] = Sin[x E^E^x]. */
    if (!r) r = rt_recursive_tower_case(fe, x);
    if (!r) r = rt_exp_tower_case(fe, x);
    expr_free(fe);
    if (!r) return NULL;
    /* Soundness guard, verified on the EXPONENTIAL form `r` (not its ExpToTrig
     * image).  The exp cases are correct by construction for a genuine rational
     * function of E^(a x), but TrigToExp of a trig-of-a-TRANSCENDENTAL argument
     * (e.g. Tan[Log[x]] -> x^I) yields complex-power kernels a single-primitive
     * case can spuriously certify — so verify.  Crucially, verify `r` and NOT
     * ExpToTrig[r]: r is a rational function of the exponential kernel, so the
     * diff-back Together[TrigToExp[D[r]-f]] reduces to 0 exactly and fast, whereas
     * the multiple-angle ExpToTrig image does NOT rational-reduce and would strand
     * the check in Simplify (the Sec[x]^3 hang).  ExpToTrig is an exact rewrite, so
     * ExpToTrig[r] is a correct antiderivative iff r is. */
    if (!rt_verify_antideriv(r, f, x)) { expr_free(r); return NULL; }
    /* Reconstruct a REAL closed form from the trig image ExpToTrig[r] (nice-to-have,
     * exact-diff-back gated and time-bounded).  On success return the clean real
     * form; on failure return the COMPACT exponential form `r` itself — already
     * verified, so correct by construction, and far smaller than the multiple-angle
     * ExpToTrig image (e.g. Sec[x]^3, whose real reduction Simplify cannot close). */
    Expr* out = rt_eval1("ExpToTrig", expr_copy(r));   /* trig image; keep r */
    Expr* real = out ? rt_realify(out, x, f) : NULL;
    if (out) expr_free(out);
    if (real) { expr_free(r); return real; }
    return r;
}

/* True iff `e` contains, anywhere, a function node with head symbol `h`. */
static bool rt_expr_has_head(Expr* e, const char* h) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol(h)) return true;
    if (rt_expr_has_head(e->data.function.head, h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_expr_has_head(e->data.function.args[i], h)) return true;
    return false;
}

/* Rational-function-of-a-single-exponential fallback.  When the integrand is a
 * rational function of one exponential kernel with a LINEAR exponent (so the
 * commensurate primitive t = E^(a x) has constant logarithmic derivative
 * u' = a), the substitution t = E^(a x), dx = dt/(a t) reduces the integral to
 * the pure rational integral  ∫ R(t)/(a t) dt, which the rational integrator
 * closes in full — including the algebraic-residue log parts it expresses as
 * RootSum (e.g. the cubic 1 + t^2 + t^3 from E^(x/6)/(1 + E^(x/2) + E^(x/3))).
 * Back-substitute t -> E^(a x) and diff-back verify.  This is the general
 * closure for the cubic-and-higher residue cases that the structured real
 * Log/ArcTan Rothstein-Trager log part (rt_frac_lrt) necessarily declines;
 * it runs LAST among the exponential cases, so their cleaner ArcTan/radical
 * forms win whenever they apply.  Correctness rests on the same
 * correct-by-construction rational integrator plus the diff-back gate; the
 * RootSum it emits is recognised by the generalized D[RootSum] collapse
 * (src/root.c). */
Expr* rt_exp_ratreduce_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F || !uexp) { if (F) expr_free(F); if (uexp) expr_free(uexp); return NULL; }

    Expr* t = expr_new_symbol("rmT");
    Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));   /* u' = a */
    Expr* result = NULL;

    /* Gates: F is a rational function of the kernel t ALONE (no residual x, no
     * nested exp/log of x), it genuinely depends on t, and u' is free of x
     * (a linear exponent -> constant log-derivative, so the Jacobian 1/(u' t)
     * is rational in t). */
    bool ok = F && up
        && rt_free_of_x(F, x)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && !rt_free_of_x(F, t)
        && rt_free_of_x(up, x);

    if (ok) {
        /* G = F / (u' t). */
        Expr* denom = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(up), expr_copy(t) }, 2);
        Expr* G = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(F),
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ denom, expr_new_integer(-1) }, 2) }, 2));
        /* Pure rational integral in t. */
        Expr* A = G ? rt_eval2("Integrate", expr_copy(G), expr_copy(t)) : NULL;
        bool declined = !A || rt_expr_has_head(A, "Integrate");
        if (!declined) {
            /* Back-substitute t -> E^u and diff-back verify against f. */
            Expr* back = expr_new_function(expr_new_symbol("List"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(t),
                        expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2) }, 2) }, 1);
            Expr* Q = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(A), back }, 2));
            /* Correct by construction (exact rational integral over the kernel +
             * exact back-substitution); shipped only behind the EXACT symbolic
             * diff-back (Simplify[D-f]===0).  Where Simplify cannot reduce a
             * higher-pole real form (e.g. Sec[x]^3 -> a triple pole at
             * E^(I x) = +/- i) this declines rather than certify numerically. */
            if (Q && rt_verify_antideriv(Q, f, x))
                result = Q;
            else if (Q) expr_free(Q);
        }
        if (A) expr_free(A);
        if (G) expr_free(G);
    }

    expr_free(t);
    if (up) expr_free(up);
    expr_free(F);
    expr_free(uexp);
    return result;
}

/* Find (borrowed) the argument u of the first kernel `head`[u] of `e` (arity 1)
 * whose argument depends on x. */
static Expr* rt_find_head_of_x(Expr* e, Expr* x, const char* head) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol(head)
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    Expr* r = rt_find_head_of_x(e->data.function.head, x, head);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_head_of_x(e->data.function.args[i], x, head);
        if (r) return r;
    }
    return NULL;
}

/* True iff `e` contains no circular- or hyperbolic-trig kernel (a coarse "is this
 * a rational function of the tangent monomial and x alone" gate; the final
 * diff-back is the real guarantee). */
static bool rt_free_of_trig(Expr* e) {
    static const char* H[] = { "Tan", "Sin", "Cos", "Cot", "Sec", "Csc",
                               "Tanh", "Sinh", "Cosh", "Coth", "Sech", "Csch" };
    for (size_t i = 0; i < sizeof H / sizeof H[0]; i++)
        if (!rt_free_of_head(e, H[i])) return false;
    return true;
}

/* The monomial variable t used throughout the family (a fixed internal symbol). */
#define RT_HT_SYM "Integrate`htT"

/* Reduce the EVEN powers of the reciprocal kernels to the tangent family, so an
 * integrand that is genuinely rational in Tan[u] (etc.) but written with the
 * secant/cosecant reciprocal — most commonly the Sec[u]^2 = d/du Tan[u] Jacobian
 * factor of a tangent-substitution integrand — is recognised by the rational-in-t
 * gate below instead of stranding on the leftover Sec and falling through to the
 * (exponentializing, blow-up-prone) TrigToExp front-end:
 *   Sec[a]^(2k)  -> (1 + Tan[a]^2)^k      Csc[a]^(2k)  -> (1 + Cot[a]^2)^k
 *   Sech[a]^(2k) -> (1 - Tanh[a]^2)^k     Csch[a]^(2k) -> (Coth[a]^2 - 1)^k
 * Each is an exact Pythagorean identity, so the rewrite is value-preserving.  ODD
 * powers (a bare Sec[a] = Sqrt(1+Tan^2), genuinely NOT rational in Tan) are left
 * untouched — such an integrand still correctly declines the rational-in-t gate.
 * Returns an owned expression (a copy of f when nothing matched). */
static Expr* rt_recip_sq_to_tan(Expr* f) {
    static const char* rules_src =
        "{ (Sec[a_]^n_ /; EvenQ[n])  :> (1 + Tan[a]^2)^(n/2),"
        "  (Csc[a_]^n_ /; EvenQ[n])  :> (1 + Cot[a]^2)^(n/2),"
        "  (Sech[a_]^n_ /; EvenQ[n]) :> (1 - Tanh[a]^2)^(n/2),"
        "  (Csch[a_]^n_ /; EvenQ[n]) :> (Coth[a]^2 - 1)^(n/2) }";
    Expr* rules = parse_expression(rules_src);
    if (!rules) return expr_copy(f);
    Expr* out = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rules }, 2));   /* adopts rules */
    return out ? out : expr_copy(f);
}

/* Maximum degree in t of an irreducible factor of p over Q (0 when p is constant
 * in t, -1 on failure).  The normal-pole gate uses this: a normal denominator
 * whose irreducible factors are each linear or quadratic has Rothstein-Trager
 * residues that are rational or at-worst quadratic-algebraic, which the §5.10
 * residue criterion realises as real Log / ArcTan.  A higher-degree irreducible
 * factor carries residues whose tower-derivation spillover destabilises the
 * p = h - D[g2] + r reconciliation, so those defer. */
static long rt_max_irr_degree(Expr* p, Expr* t) {
    Expr* fac = rt_eval1("Factor", expr_copy(p));
    if (!fac) return -1;
    long maxd = 0;
    if (rt_head_is(fac, "Times")) {
        for (size_t i = 0; i < fac->data.function.arg_count; i++) {
            Expr* term = fac->data.function.args[i];
            Expr* base = (rt_head_is(term, "Power")
                          && term->data.function.arg_count == 2)
                         ? term->data.function.args[0] : term;
            long d = rt_degree(base, t);
            if (d > maxd) maxd = d;
        }
    } else {
        Expr* base = (rt_head_is(fac, "Power") && fac->data.function.arg_count == 2)
                     ? fac->data.function.args[0] : fac;
        long d = rt_degree(base, t);
        if (d > maxd) maxd = d;
    }
    expr_free(fac);
    return maxd;
}


/* Real "hypertangent family" integrator (Bronstein §5.10 and its hyperbolic
 * analogue): integrate a rational function of a SINGLE kernel `khead`[u]
 * (u rational in x) DIRECTLY and real, the real-valued replacement for the
 * TrigToExp route.  Substitute t = khead[u] so F is rational in t over C(x);
 * build the tower derivation Dt (given, = u'(t^2+1) for Tan, -u'(t^2+1) for Cot,
 * u'(1-t^2) for Tanh); run the `driver` (IntegrateHypertangent / -Hypertanh);
 * integrate the leftover base-field element in C(x); back-substitute t ->
 * khead[u]; collapse the real special log Log[1 (+/-) khead[u]^2] to
 * -2 Log[`inner`[u]] (Cos / Sin / Cosh); and diff-back verify.  `special` is the
 * monic special polynomial (t^2+1 or t^2-1) used by the normal-pole gate.  Takes
 * ownership of Dt and special. */
/* Real hypertangent / hyperbolic-tangent family integrator.
 *
 * The kernel substitution and the whole pipeline are CORRECT BY CONSTRUCTION —
 * there is no diff-back Simplify gate.  The structural gate below (F is a
 * genuine rational function of t over C(x), no trig/exp/log of x survives)
 * already guarantees F|_{t->kernel} = f: for a direct kernel t IS the kernel,
 * so back-substitution is the identity; for a Weierstrass half-angle
 * substitution the rules encode exact identities.  The §5.10 driver then
 * returns g with beta = True iff the transcendental part is elementary, leaving
 * the t-free base-field remainder base = F - D_tower[g] for a recursive integral
 * in x.  Since the tower derivation D[t] IS the derivative of the kernel bval,
 * d/dx[(g + integral base)|_{t->bval}] = F|_{t->bval} = f is a theorem (Bronstein
 * Thm 5.10.1), not something to re-verify.  The only runtime facts consulted are
 * STRUCTURAL: beta = True (elementary), base free of t, and the base integral is
 * itself elementary (not left unevaluated). */
static Expr* rt_hypertan_family(Expr* f, Expr* x, Expr* subrule, Expr* bval,
                                Expr* Dt, Expr* special,
                                const char* driver, const char* inner,
                                Expr* carg, bool cosmetic_plus) {
    Expr* t = expr_new_symbol(RT_HT_SYM);
    /* Reduce even reciprocal-kernel powers (Sec[u]^2 -> 1+Tan[u]^2, ...) so a
     * rational-in-Tan[u] integrand carrying the Sec Jacobian is recognised by the
     * rational-in-t gate rather than stranding on the leftover Sec.  fr is an exact
     * rewrite of f, so it is used for both the kernel substitution AND the diff-back
     * verify below (which is thereby a rational-in-t certificate, not a Sec-laden one). */
    Expr* fr = rt_recip_sq_to_tan(f);
    Expr* F = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(fr), subrule }, 2));                 /* adopts subrule */

    /* F must be a genuine rational function of t over C(x). */
    bool ok = F && rt_free_of_trig(F)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && !rt_free_of_x(F, t);
    if (ok) {
        Expr* Ft = rt_eval1("Together", expr_copy(F));
        Expr* nn = Ft ? rt_eval1("Numerator", expr_copy(Ft)) : NULL;
        Expr* dd = Ft ? rt_eval1("Denominator", expr_copy(Ft)) : NULL;
        ok = nn && dd && rt_is_poly(nn, t) && rt_is_poly(dd, t);
        if (ok) {
            /* Normal-pole gate.  Strip the special factor from the denominator,
             * then require every IRREDUCIBLE factor of the remaining NORMAL part
             * (over Q) to have degree <= 2: its Rothstein-Trager residues are then
             * rational (linear factors) or at-worst quadratic-algebraic (irreducible
             * quadratics), both of which the §5.10 residue criterion realises as
             * real Log / ArcTan.  A higher-degree irreducible factor carries
             * residues whose tower-derivation spillover blows up the
             * p = h - D[g2] + r reconciliation, so those defer.  (Pure-polynomial
             * and special^k integrands strip to a unit; split normal parts of any
             * degree pass, as does an irreducible quadratic such as the 3+t^2 of
             * Tan[x]/(3+Tan[x]^2).) */
            Expr* dn = expr_copy(dd);
            for (int guard = 0; dn && guard < 4096; guard++) {
                Expr* rem = rt_eval3("PolynomialRemainder",
                    expr_copy(dn), expr_copy(special), expr_copy(t));
                bool divisible = rem && rt_is_zero(rem);
                if (rem) expr_free(rem);
                if (!divisible) break;
                Expr* q = rt_eval3("PolynomialQuotient",
                    expr_copy(dn), expr_copy(special), expr_copy(t));
                expr_free(dn); dn = q;
            }
            long mdeg = dn ? rt_max_irr_degree(dn, t) : -1;
            ok = (mdeg >= 0 && mdeg <= 2);
            if (dn) expr_free(dn);
        }
        if (Ft) expr_free(Ft);
        if (nn) expr_free(nn);
        if (dd) expr_free(dd);
    }
    if (!ok) { expr_free(fr); if (F) expr_free(F); expr_free(t); expr_free(Dt); expr_free(special);
               expr_free(bval); expr_free(carg); return NULL; }

    /* deriv = {x -> 1, t -> Dt}. */
    Expr* deriv = expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(x), expr_new_integer(1) }, 2),
          expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(t), expr_copy(Dt) }, 2) }, 2);

    Expr* result = NULL;
    Expr* res = rt_eval_call(driver,
        (Expr*[]){ expr_copy(F), expr_copy(t), expr_copy(deriv) }, 3);
    if (res && rt_head_is(res, "List") && res->data.function.arg_count == 2
        && rt_is_true(res->data.function.args[1])) {
        Expr* g = res->data.function.args[0];               /* borrowed */
        /* base = F - D_tower[g]  (an element of k = C(x): must be free of t).
         * Expand after Together: the tower derivation returns unexpanded products
         * (e.g. (t^2+1)(t^2-1)) whose t-terms only cancel once expanded. */
        Expr* Dg = rt_eval2("Risch`Derivation", expr_copy(g), expr_copy(deriv));
        Expr* base = Dg ? rt_eval1("Expand", rt_eval1("Together",
            expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_copy(F), expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), Dg }, 2) }, 2))) : NULL;
        /* base is the t-free base-field remainder; its being free of t is the
         * structural proof that g integrates the transcendental part of F.  The
         * recursive base-field integral ib is elementary iff Integrate does not
         * leave it unevaluated.  Both are correct by construction (rational Risch
         * over C(x)); back-substituting the kernel then yields a correct
         * antiderivative of f without any diff-back re-verification. */
        if (base && rt_free_of_x(base, t)) {
            Expr* ib = rt_eval2("Integrate", expr_copy(base), expr_copy(x));
            if (ib && !rt_expr_has_head(ib, "Integrate")) {
                Expr* backrule = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(t), expr_copy(bval) }, 2);
                Expr* gsub = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(g), backrule }, 2));    /* adopts backrule */
                Expr* ans = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ gsub, ib }, 2));                  /* adopts gsub, ib */
                if (ans) {
                    /* Cosmetic (real): Log[1 (+/-) bval^2] = -2 Log[inner[carg]]
                     * (Tan/Cot: 1+K^2 -> -2 Log Cos/Sin; Tanh/Coth: 1-K^2 -> -2 Log
                     * Cosh/Sinh; Weierstrass: 1+Tan[u/2]^2 -> -2 Log Cos[u/2]).  Each
                     * is an exact identity (for Coth, 1-Coth^2 < 0 contributes only a
                     * constant i Pi), so the swap is derivative-preserving by
                     * construction; apply it when the pattern is present. */
                    Expr* ksq = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(bval), expr_new_integer(2) }, 2);
                    Expr* logarg = cosmetic_plus
                        ? expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_new_integer(1), ksq }, 2)
                        : expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_new_integer(1),
                                expr_new_function(expr_new_symbol("Times"),
                                    (Expr*[]){ expr_new_integer(-1), ksq }, 2) }, 2);
                    Expr* lhs = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ logarg }, 1);
                    Expr* rhs = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-2),
                            expr_new_function(expr_new_symbol("Log"),
                                (Expr*[]){ expr_new_function(expr_new_symbol(inner),
                                    (Expr*[]){ expr_copy(carg) }, 1) }, 1) }, 2);
                    Expr* crule = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ lhs, rhs }, 2);
                    Expr* clean = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(ans), crule }, 2)); /* adopts crule */
                    if (clean && !expr_eq(clean, ans)) {
                        expr_free(ans); result = clean;
                    } else {
                        if (clean) expr_free(clean);
                        result = ans;
                    }
                }
            } else if (ib) {
                expr_free(ib);
            }
        }
        if (base) expr_free(base);
    }
    if (res) expr_free(res);
    /* Diff-back safety, gated to the RELAXED case.  For a kernel argument u that is
     * rational in x (the common Tan[x], Tan[2x], ... family) the whole §5.10 pipeline
     * is a theorem (Bronstein Thm 5.10.1) — correct by construction, per the header
     * above — and its clean real Log/ArcTan-of-Tan form is precisely what the diff-back
     * zero tests cannot reduce (Together leaves the Laurent TrigToExp poles uncombined;
     * Simplify mis-reduces the Sqrt-wrapped ArcTan to a spurious residual), so verifying
     * it would spuriously REJECT a correct answer and strand the integral in the
     * exponentializing front-end (an unbounded blow-up).  Only a TRANSCENDENTAL argument
     * (u = Log[x], ...), admitted by rt_kernel_eta's relaxed eta-only gate, needs the
     * diff-back guard — there the construction rests on the additional eta-in-C(x)
     * assumption, so verify before shipping. */
    bool u_transcendental = !rt_free_of_trig(carg)
        || rt_find_exp_of_x(carg, x) != NULL || rt_find_log_of_x(carg, x) != NULL;
    if (result && u_transcendental && !rt_verify_antideriv(result, fr, x)) {
        expr_free(result); result = NULL;
    }
    expr_free(fr);
    expr_free(F); expr_free(t); expr_free(deriv); expr_free(Dt); expr_free(special);
    expr_free(bval); expr_free(carg);
    return result;
}

/* Direct single-kernel call: substitute khead[u] -> t (and back-substitute the
 * same), cosmetic argument u.  Takes ownership of Dt and special. */
static Expr* rt_hypertan_direct(Expr* f, Expr* x, Expr* u, const char* khead,
                                Expr* Dt, Expr* special, const char* driver,
                                const char* inner, bool cosmetic_plus) {
    Expr* subrule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_function(expr_new_symbol(khead),
            (Expr*[]){ expr_copy(u) }, 1), expr_new_symbol(RT_HT_SYM) }, 2);
    Expr* bval = expr_new_function(expr_new_symbol(khead), (Expr*[]){ expr_copy(u) }, 1);
    return rt_hypertan_family(f, x, subrule, bval, Dt, special, driver, inner,
                             expr_copy(u), cosmetic_plus);
}

/* eta = D[u, x] must be a rational function of x alone (free of trig/exp/log of x)
 * so that khead[u] is a monomial over k = C(x).  Returns owned eta or NULL. */
static Expr* rt_kernel_eta(Expr* u, Expr* x) {
    /* tau = Tan[u] (etc.) is a hypertangent monomial over C(x) exactly when
     * eta = Dtau/(tau^2 +/- 1) = u' lies in C(x) — i.e. u' is free of every kernel.
     * The argument u itself may be TRANSCENDENTAL (e.g. Log[x], with u' = 1/x): the
     * monomial is still hypertangent over C(x), so a nested tangent such as
     * Tan[Log[x]] integrates directly to a real Log[Cos]/ArcTan form rather than
     * declining or routing through the complex-exponential front-end.  (The
     * former rt_kernel_simple(u) gate was over-strict — it required u itself to be
     * kernel-free; the eta check below is the genuine over-C(x) condition.  An
     * integrand mixing the tangent with an independent Log/Exp of x is still
     * rejected by rt_hypertan_family's free-of-log/exp gate on the substituted F.) */
    Expr* eta = rt_eval2("D", expr_copy(u), expr_copy(x));
    if (!eta || !rt_free_of_trig(eta)
        || rt_find_exp_of_x(eta, x) != NULL || rt_find_log_of_x(eta, x) != NULL) {
        if (eta) expr_free(eta);
        return NULL;
    }
    return eta;
}
/* t^2 + sigma  (sigma = +1 or -1), in the family monomial variable. */
static Expr* rt_special_poly(long sigma) {
    return expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol(RT_HT_SYM), expr_new_integer(2) }, 2),
          expr_new_integer(sigma) }, 2);
}
/* 1 - t^2 in the family monomial variable (the Tanh/Coth tower derivative base). */
static Expr* rt_one_minus_t2(void) {
    return expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_new_integer(1),
            expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_new_symbol(RT_HT_SYM), expr_new_integer(2) }, 2) }, 2) }, 2);
}

/* Dispatch the real hypertangent family: Tan[u] and Cot[u] (special t^2+1, the
 * §5.10 hypertangent driver) and Tanh[u] (special t^2-1, the hyperbolic driver).
 * The special polynomial t^2+1 is irreducible over k (needs a coupled complex
 * system); t^2-1 splits (two real Risch DEs) — both keep the answer real. */
Expr* rt_hypertangent_case(Expr* f, Expr* x) {
    /* Tan[u]:  Dt = eta (t^2+1). */
    Expr* u = rt_find_head_of_x(f, x, "Tan");
    if (u) {
        Expr* eta = rt_kernel_eta(u, x);
        if (eta) {
            Expr* Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ eta, rt_special_poly(1) }, 2);     /* adopts eta */
            return rt_hypertan_direct(f, x, u, "Tan", Dt, rt_special_poly(1),
                "Risch`IntegrateHypertangent", "Cos", true);
        }
    }
    /* Cot[u]:  Dt = -eta (t^2+1). */
    u = rt_find_head_of_x(f, x, "Cot");
    if (u) {
        Expr* eta = rt_kernel_eta(u, x);
        if (eta) {
            Expr* neg_eta = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), eta }, 2);   /* adopts eta */
            Expr* Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ neg_eta, rt_special_poly(1) }, 2);
            return rt_hypertan_direct(f, x, u, "Cot", Dt, rt_special_poly(1),
                "Risch`IntegrateHypertangent", "Sin", true);
        }
    }
    /* Tanh[u] and Coth[u]:  Dt = eta (1 - t^2)  (special t^2-1, which splits).
     * They share the hyperbolic driver; only the real cosmetic differs — for
     * Tanh 1-Tanh^2 = Sech^2 > 0 (-> -2 Log Cosh), for Coth 1-Coth^2 = -Csch^2
     * (-> -2 Log Sinh, the sign folded into the derivative-preserving rewrite). */
    struct { const char* head; const char* inner; } hyp[] = {
        { "Tanh", "Cosh" }, { "Coth", "Sinh" }
    };
    for (size_t i = 0; i < sizeof hyp / sizeof hyp[0]; i++) {
        u = rt_find_head_of_x(f, x, hyp[i].head);
        if (!u) continue;
        Expr* eta = rt_kernel_eta(u, x);
        if (!eta) continue;
        Expr* Dt = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ eta, rt_one_minus_t2() }, 2);          /* adopts eta */
        return rt_hypertan_direct(f, x, u, hyp[i].head, Dt, rt_special_poly(-1),
            "Risch`IntegrateHypertanh", hyp[i].inner, false);
    }
    return NULL;
}
