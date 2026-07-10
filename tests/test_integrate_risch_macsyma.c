/* test_integrate_risch_macsyma.c
 *
 * Tests for the Maxima-ported Risch integrator (Integrate`RischMacsyma,
 * Method -> "RischMacsyma").  Correctness is asserted by the universal
 * predicate Simplify[D[Integrate[f, x], x] - f] === 0 rather than by
 * fixed output strings, so the tests survive surface-form changes.
 *
 * The suite grows with the implementation phases:
 *   Phase 1 — rational case + dispatch plumbing + verification gate.
 *   Phase 2 — logarithmic case.
 *   Phase 3 — exponential case.
 *   Phase 4 — polynomial-part recursion / hardening.
 *   Phase 5 — special functions (gated) + trig/hyperbolic front-end.
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert that the evaluated `input` is an unevaluated call with head
 * `head` (used for strict no-fallback / out-of-scope cases). */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* Assert Simplify[D[Integrate`RischMacsyma[f, x], x] - (f)] === 0, i.e.
 * the explicit package head produces a correct antiderivative. */
static void assert_rm_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate`RischMacsyma[%s, x], x] - (%s)]", f, f);
    assert_eval_eq(buf, "0", 0);
}

/* Assert the same via the Method-option surface form. */
static void assert_rm_method_diff_zero(const char* f) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Simplify[D[Integrate[%s, x, Method -> \"RischMacsyma\"], x] - (%s)]",
        f, f);
    assert_eval_eq(buf, "0", 0);
}

/* Fast numeric verification for special-function antiderivatives, whose
 * exact Simplify[D - f] is prohibitively slow (Erf/Gaussian).  Confirms
 * closure (non-unevaluated) and correctness at two interior points. */
static void assert_rm_num(const char* f) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
        "Abs[N[(D[Integrate`RischMacsyma[%s, x], x] - (%s)) /. x -> 13/10]] +"
        " Abs[N[(D[Integrate`RischMacsyma[%s, x], x] - (%s)) /. x -> 17/10]]"
        " < 1/100000", f, f, f, f);
    assert_eval_eq(buf, "True", 0);
}

/* ---------------- Phase 1: rational case ---------------- */

static void test_rational_case(void) {
    /* Log-producing residues, polynomial parts, repeated poles, and the
     * pure-power case (which must NOT emit a spurious Log). */
    assert_rm_diff_zero("1/x");
    assert_rm_diff_zero("1/(x^2 - 1)");
    assert_rm_diff_zero("1/(x^2 + 1)");
    assert_rm_diff_zero("x/(x^2 + 1)");
    assert_rm_diff_zero("(x^3 + 1)/(x^2 + 1)");
    assert_rm_diff_zero("1/(x^3 - 1)");
    assert_rm_diff_zero("1/(x^2 + 1)^2");
    assert_rm_diff_zero("1/x^2");
    assert_rm_diff_zero("1/(x^4 + 1)");
    assert_rm_diff_zero("(2 x + 3)/(x^2 + 3 x + 2)");
    assert_rm_diff_zero("5");
    assert_rm_diff_zero("x^3");

    /* RischMacsyma agrees with BronsteinRational on the rational case. */
    assert_eval_eq(
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischMacsyma\"]"
        " - Integrate[1/(x^2 - 1), x, Method -> \"BronsteinRational\"]]",
        "0", 0);
}

/* ---------------- Phase 2: recursive logarithmic case ---------------- */
/* Genuine recursive Risch (primitive polynomial coefficient matching),
 * NOT the parallel-Risch pmint engine. */
static void test_logarithmic_case(void) {
    assert_rm_diff_zero("Log[x]");            /* x Log[x] - x            */
    assert_rm_diff_zero("Log[x]^2");
    assert_rm_diff_zero("Log[x]^3");
    assert_rm_diff_zero("x Log[x]");
    assert_rm_diff_zero("(x^2 + 1) Log[x]");
    assert_rm_diff_zero("x^2 Log[x]");
    assert_rm_diff_zero("Log[x^2 + 1]");      /* -> x Log[..] - 2x + 2 ArcTan[x] */
    assert_rm_diff_zero("Log[2 x + 3]");      /* -> (x + 3/2) Log[2x+3] - x      */
    /* The limited-integration oracle folds the would-be new logarithm
     * back into the tower, closing this as Log[x]^2/2. */
    assert_rm_diff_zero("Log[x]/x");
    assert_rm_diff_zero("Log[x]^2/x");        /* -> Log[x]^3/3                   */
}

/* ---------------- Phase 4: fractional (Rothstein-Trager) log-part ------- */
/* Proper rational functions of a single transcendental theta with a
 * squarefree denominator integrate to sum_i c_i Log(g_i) with constant
 * residues, found by an exact SolveAlways over t and x. */
static void test_fractional_case(void) {
    assert_rm_diff_zero("1/(x (1 + Log[x]))");        /* Log[1 + Log[x]]     */
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])");       /* Log[1 + E^x]        */
    assert_rm_diff_zero("Exp[x]/(Exp[x] - 1)");       /* Log[-1 + E^x]       */
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[x]))"); /* two residues        */
    assert_rm_diff_zero("1/(x (Log[x]^2 - 1))");      /* two residues        */
    assert_rm_diff_zero("Exp[x]/(2 + 3 Exp[x])");
    assert_rm_diff_zero("(2 Exp[x])/(1 + Exp[x])");
    /* (repeated poles are handled by the Hermite case, see test_hermite_case) */
}

/* ---------------- Phase 7: Hermite reduction (repeated poles) ---------- */
/* A proper rational function of theta = Log[u] with a REPEATED pole reduces
 * to Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j), Hden = gcd(D, dD/dtheta),
 * solved by SolveAlways over {t, x}. */
static void test_hermite_case(void) {
    assert_rm_diff_zero("1/(x (1 + Log[x])^2)");   /* -1/(1 + Log[x])      */
    assert_rm_diff_zero("1/(x Log[x]^2)");         /* -1/Log[x]           */
    assert_rm_diff_zero("1/(x (1 + Log[x])^3)");   /* -1/(2 (1+Log[x])^2) */
    assert_rm_diff_zero("(1 + Log[x])/(x Log[x]^2)"); /* Hermite + log part */
    assert_rm_diff_zero("1/(x (2 + Log[x])^2)");
    assert_rm_diff_zero("Log[x]/(x (1 + Log[x])^2)");
    /* Exponential-kernel Hermite (D coprime to theta = E^u). */
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^2");     /* -1/(1 + E^x)        */
    assert_rm_diff_zero("Exp[x]/(1 - Exp[x])^2");
    assert_rm_diff_zero("Exp[x]/(2 + Exp[x])^2");
}

/* ---------------- Phase 6: coupled hyperexponential case --------------- */
/* A rational function of E^u whose integral mixes a Laurent-polynomial part
 * with a log part (D(g)/g improper in theta), solved by the unified ansatz
 * Q = sum_i w_i(x) theta^i + sum_j c_j Log(g_j) via SolveAlways over {t, x}. */
static void test_hyperexponential_case(void) {
    assert_rm_diff_zero("1/(1 + Exp[x])");       /* x - Log[1 + E^x]        */
    assert_rm_diff_zero("1/(1 - Exp[x])");       /* x - Log[-1 + E^x]       */
    assert_rm_diff_zero("1/(2 + Exp[x])");
    assert_rm_diff_zero("Exp[2 x]/(1 + Exp[x])");  /* improper: poly + log part */
    assert_rm_diff_zero("(1 + Exp[x])/(2 + Exp[x])");

    /* Phase A: exponential Hermite coupled with the Laurent (theta = 0 pole)
     * part.  A REPEATED transcendental pole for theta = E^u, possibly together
     * with a theta = 0 pole, needs the unified ansatz
     *   Q = sum_i w_i(x) theta^i + H(theta)/Hden(theta) + sum_j c_j Log(g_j).
     * These decline under the pure Hermite (which cannot supply the x/Laurent
     * term) and the squarefree hyperexp path. */
    assert_rm_diff_zero("1/(1 + Exp[x])^2");        /* x + 1/(1+E^x) - Log[1+E^x] */
    assert_rm_diff_zero("1/(Exp[x] (1 + Exp[x])^2)"); /* theta=0 pole + repeated  */
    assert_rm_diff_zero("1/(2 + Exp[x])^2");
    assert_rm_diff_zero("1/(1 + Exp[x])^3");        /* higher repeated pole     */
    assert_rm_diff_zero("Exp[x]/(1 + Exp[x])^3");
}

/* ---------------- Phase 2: special functions (Maxima erfarg/dilog/Ei/li) --- */
/* These close integrals the whole elementary cascade leaves open, using
 * special functions Mathilda already provides.  On by default. */
static void test_special_functions(void) {
    assert_rm_num("Exp[-x^2]");              /* -> (Sqrt[Pi]/2) Erf[x] */
    assert_rm_num("Exp[x^2]");               /* -> Erfi / Erf[I x]     */
    assert_rm_num("Exp[-x^2/2]");
    assert_rm_num("Exp[-3 x^2 + 2 x + 1]");  /* general Gaussian       */
    assert_rm_num("Exp[x]/x");               /* -> ExpIntegralEi[x]    */
    assert_rm_num("Exp[2 x]/(x - 1)");       /* -> shifted Ei          */
    assert_rm_num("1/Log[x]");               /* -> LogIntegral[x]      */
    assert_rm_num("3/Log[x]");
    assert_rm_num("Log[1 - x]/x");           /* -> -PolyLog[2, x]      */
    assert_rm_num("Log[1 + x]/x");           /* -> -PolyLog[2, -x]     */
}

/* ---------------- Phase 3: exponential-polynomial case (Risch DE) ------- */
/* theta = E^u, theta' = u' theta; powers of theta decouple into
 * q_i' + i u' q_i = p_i solved exactly by a polynomial ansatz. */
static void test_exponential_case(void) {
    assert_rm_diff_zero("x Exp[x]");            /* (x - 1) E^x            */
    assert_rm_diff_zero("x^2 Exp[x]");          /* (x^2 - 2x + 2) E^x     */
    assert_rm_diff_zero("(x^2 + 1) Exp[x]");
    assert_rm_diff_zero("x^3 Exp[x]");
    assert_rm_diff_zero("Exp[2 x]");            /* (1/2) E^(2 x)          */
    assert_rm_diff_zero("x Exp[2 x]");
    assert_rm_diff_zero("x Exp[x^2]");          /* (1/2) E^(x^2)          */
    assert_rm_diff_zero("x Exp[-x^2]");         /* -(1/2) E^(-x^2)        */
    assert_rm_diff_zero("x Exp[3 x + 1]");
    assert_rm_diff_zero("x + x Exp[x]");        /* i=0 and i=1 terms      */
    /* E^(-x^2) itself has no polynomial RDE solution (N < 0): the exp case
     * declines and the Erf recognizer takes over.  Confirm it still closes
     * (to an Erf form) rather than being wrongly forced through the DE. */
    assert_rm_num("Exp[-x^2]");
    /* E^x/x^2 needs ExpIntegralEi with a non-linear denominator, which no
     * implemented case handles: it stays unevaluated (never garbage). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Exp[x]/x^2, x]", "Integrate`RischMacsyma");
}

/* ---------------- Phase 5: trig / hyperbolic front-end ----------------- */
/* TrigToExp -> exponential (Laurent) machinery -> ExpToTrig.  Cases that
 * reduce to a Laurent polynomial in the exponential kernel close; those that
 * reduce to a coupled fractional exponential (Tan, Tanh) decline. */
static void test_trig_frontend(void) {
    assert_rm_num("Sin[x]");            /* -Cos[x]                */
    assert_rm_num("Cos[x]");            /* Sin[x]                 */
    assert_rm_num("Sinh[x]");           /* Cosh[x]                */
    assert_rm_num("Cosh[x]");           /* Sinh[x]                */
    assert_rm_num("Sin[x]^2");          /* x/2 - Sin[2x]/4        */
    assert_rm_num("Cos[x]^2");
    assert_rm_num("Sin[x] Cos[x]");     /* -Cos[2x]/4             */
    assert_rm_num("Sin[2 x]");
    assert_rm_num("Cosh[x]^2");
    /* Also reachable in exponential form directly. */
    assert_rm_diff_zero("(Exp[x] + Exp[-x])/2");   /* Sinh[x]     */
    assert_rm_diff_zero("Exp[x] + Exp[-x] + Exp[2 x]");
    /* Tan / Tanh close via the coupled hyperexponential case: the result is
     * correct (diff-back is 0) though, through the complex substitution, it is
     * left in an I-laden form such as I x - Log[1 + E^(2 I x)] (= -Log[Cos[x]])
     * that no current simplifier reduces to real closed form.  Completeness is
     * preferred over declining; the residual form is a Simplify improvement
     * opportunity (see rm_trig_frontend). */
    assert_rm_diff_zero("Tan[x]");
    assert_rm_diff_zero("Tanh[x]");
}

/* ---------------- Phase B (first increment): multi-kernel exponential sums --- */
/* Integrands that exponentialize to a sum of NON-commensurate exponentials
 * E^((a +/- b I) x) — e.g. a real exponential times a trig factor — are two
 * independent transcendental extensions.  The distinct exponentials decouple
 * (d/dx never mixes them), so each term p_k(x) E^(W_k) integrates via its own
 * Risch DE q_k' + W_k' q_k = p_k (rm_expsum_case).  The single-primitive
 * exponential cases cannot kernelize them.  Correct by construction (each q_k
 * is SolveAlways-certified); through the complex exponentials the ExpToTrig
 * output is left in an I-laden Cosh/Sinh form (a Simplify opportunity), but the
 * diff-back is exactly 0. */
static void test_multikernel_case(void) {
    assert_rm_diff_zero("Exp[x] Sin[x]");        /* (E^x/2)(Sin[x] - Cos[x]) */
    assert_rm_diff_zero("Exp[x] Cos[x]");
    assert_rm_diff_zero("x Exp[x] Sin[x]");      /* polynomial coefficient   */
    assert_rm_diff_zero("Exp[2 x] Sin[3 x]");    /* non-unit real/imag parts */
    assert_rm_diff_zero("Exp[2 x] Cos[3 x]");
    assert_rm_diff_zero("x^2 Exp[x] Cos[x]");
    /* Directly as a sum of independent exponentials (no trig front-end). */
    assert_rm_diff_zero("Exp[x] + Exp[2 x] Sin[x]");
    /* A non-elementary term (E^(x^2)) inside the sum makes the whole case
     * decline cleanly rather than emit a partial/garbage answer. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Exp[x] Sin[x] + Exp[x^2], x]",
        "Integrate`RischMacsyma");
}

/* ---------------- Phase B (second increment): nested logarithmic towers ---- */
/* A rational function of a chain of NESTED logarithms t_1 = Log[u_1](x),
 * t_2 = Log[u_2](x,t_1), ... is integrated by generalizing the single-kernel
 * derivation D = d/dx + Dt d/dt to the whole tower D = d/dx + sum_i Dt_i d/dt_i
 * and solving one unified ansatz over all tower variables (rm_log_tower_case).
 * Correct by construction (SolveAlways certifies the polynomial identity in the
 * algebraically independent transcendentals). */
static void test_log_tower_case(void) {
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]])");   /* Log[Log[Log[x]]]   */
    assert_rm_diff_zero("Log[Log[x]]/(x Log[x])");     /* Log[Log[x]]^2/2    */
    assert_rm_diff_zero("Log[Log[x]]/x");              /* Log[x](Log[Log[x]]-1) */
    assert_rm_diff_zero("1/(x Log[x] Log[Log[x]] Log[Log[Log[x]]])"); /* depth 3 */
    /* Independent (non-nested) logs are also a valid triangular tower. */
    assert_rm_diff_zero("Log[x] + Log[x + 1]");
    /* CRITICAL: a residual NON-rational kernel of an inner argument
     * (Sin[Log[x]] -> Sin[t_1]) must make the case DECLINE, never certify a
     * wrong closed form.  Same for a genuinely non-elementary nested integrand
     * (Log[x] Log[Log[x]] needs ExpIntegralEi). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Sin[Log[x]] Log[Log[x]], x]",
        "Integrate`RischMacsyma");
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Log[x] Log[Log[x]], x]",
        "Integrate`RischMacsyma");
}

/* ---------------- Phase B (third increment): nested exponential towers ----- */
/* Dual of the log tower for a chain of nested exponentials t_i = E^(u_i), with
 * the exponential tower derivation Dt_i = (D[u_i,x]|_{->t}) t_i and a Laurent
 * ansatz (the inner exp kernels are invertible, so their coefficient exponents
 * are Laurent).  These are bounded-ansatz searches, so the result is diff-back
 * verified (rm_verify_antideriv) — a non-elementary integrand declines rather
 * than emitting a spurious closed form. */
static void test_exp_tower_case(void) {
    assert_rm_diff_zero("E^x E^(E^x)");            /* E^(E^x)                */
    assert_rm_diff_zero("E^(2 E^x) E^x");          /* (1/2) E^(2 E^x)        */
    /* Non-elementary nested exponentials decline cleanly (never a wrong form). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[E^(E^x), x]", "Integrate`RischMacsyma");
    /* REGRESSION (pre-existing single-kernel bug, fixed by rm_kernel_simple):
     * E^(E^x)/(1+E^(E^x)) is NON-elementary, but rm_frac_case used to substitute
     * only the outer kernel and, with the inner E^x left in the derivation as a
     * free SolveAlways parameter, certified the WRONG answer Log[1+E^(E^x)]/E^x.
     * It must now decline. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[E^(E^x)/(1 + E^(E^x)), x]", "Integrate`RischMacsyma");
}

/* ---------------- Genuine one-extension recursion (mixed / rational) ---- */
static void test_recursive_tower_case(void) {
    /* MIXED exp/log tower, independent kernels — the flat single-kind tower
     * cases cannot model a chain with both an E^x and a Log[x]. */
    assert_rm_diff_zero("E^x/x + E^x Log[x]");        /* E^x Log[x]            */
    /* Nested mixed tower (a logarithm of an exponential expression). */
    assert_rm_diff_zero("Log[1 + E^x] + x E^x/(1 + E^x)");   /* x Log[1+E^x]   */
    /* RATIONAL lower-field coefficient (answer's t_n coefficient is 1/x, not a
     * polynomial) — impossible for the flat ansatz, natural for the recursion. */
    assert_rm_diff_zero("1/(x^2 Log[x]) - Log[Log[x]]/x^2");  /* Log[Log[x]]/x */
    /* PROPER rational part in the top kernel: tower Hermite (repeated pole) and
     * the polynomial-part + proper-part split with new-log fold-back. */
    assert_rm_diff_zero("1/(x Log[x] (1 + Log[Log[x]])^2)");  /* -1/(1+Log[Log[x]]) */
    assert_rm_diff_zero("1/(x (1 + Log[x])) + E^x");          /* E^x + Log[1+Log[x]] */
    assert_rm_diff_zero(
        "Log[Log[x]]/(x Log[x] (1 + Log[Log[x]]))");   /* Log[Log[x]] - Log[1+Log[Log[x]]] */
    /* A proper fraction whose Rothstein-Trager residue is non-constant is
     * non-elementary and must decline (never certify a wrong constant residue). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[1/(Log[x] (1 + Log[Log[x]])), x]", "Integrate`RischMacsyma");
    /* EXPONENTIAL top level with a proper part (t_n-pole): the Laurent and log
     * parts couple, so a unified coupled hyperexponential ansatz over the tower
     * derivation is needed (t_n = E^(Log[x]^2), a nested exp over a log). */
    assert_rm_diff_zero(
        "(2 Log[x]/x) E^(Log[x]^2)/(1 + E^(Log[x]^2))");   /* Log[1+E^(Log[x]^2)] */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[1/(1 + E^(Log[x]^2)), x]", "Integrate`RischMacsyma");
    /* General field Risch DE with a RATIONAL lower-field solution: q = 1/Log[x]
     * (monomial denominator) and q = 1/(1+Log[x]) (NON-monomial denominator) —
     * both beyond the coupled ansatz's polynomial coefficients, solved by the
     * q = h/Denominator[p] field-RDE ansatz. */
    assert_rm_diff_zero(
        "(2/x - 1/(x Log[x]^2)) E^(Log[x]^2)");            /* E^(Log[x]^2)/Log[x] */
    assert_rm_diff_zero(
        "(2 Log[x]/(x (1 + Log[x])) - 1/(x (1 + Log[x])^2)) E^(Log[x]^2)");
                                                           /* E^(Log[x]^2)/(1+Log[x]) */
    /* Field RDEs with no elementary solution decline (never a wrong answer). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[E^(Log[x]^2), x]", "Integrate`RischMacsyma");
    assert_head_unevaluated(
        "Integrate`RischMacsyma[E^(Log[x]^2)/(1 + Log[x]), x]", "Integrate`RischMacsyma");
    /* Non-elementary mixed integrands decline cleanly (never a wrong form):
     * E^x Log[x] needs ExpIntegralEi; Log[x] Log[Log[x]] is non-elementary. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[E^x Log[x], x]", "Integrate`RischMacsyma");
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Log[x] Log[Log[x]], x]", "Integrate`RischMacsyma");
    /* A non-rational kernel of an inner argument must not certify a wrong form. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Sin[Log[x]] Log[Log[x]], x]", "Integrate`RischMacsyma");
}

/* ---------------- Automatic cascade now closes special functions ------- */
static void test_cascade_default(void) {
    /* The cascade insertion adds capability by default (WL-faithful):
     * these previously returned unevaluated.  Verified numerically (exact
     * Simplify on Erf/Ei/PolyLog derivatives is too slow for a unit test). */
    assert_eval_eq(
        "Abs[N[(D[Integrate[Exp[-x^2], x], x] - Exp[-x^2]) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[Exp[x]/x, x], x] - Exp[x]/x) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Log[x], x], x] - 1/Log[x]) /. x -> 13/10]]"
        " < 1/100000", "True", 0);
    assert_eval_eq(
        "Abs[N[(D[Integrate[Log[1 - x]/x, x], x] - Log[1 - x]/x) /. x -> 3/10]]"
        " < 1/100000", "True", 0);
    /* Regression: an ordinary rational integrand is unchanged. */
    assert_eval_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]", 0);
}

/* ---------------- Dispatch plumbing ---------------- */

static void test_method_plumbing(void) {
    /* Both surface forms route to the routine and agree. */
    assert_rm_diff_zero("1/(x^2 + 1)");
    assert_rm_method_diff_zero("1/(x^2 + 1)");
    /* The two forms produce the same antiderivative. */
    assert_eval_eq(
        "Simplify[Integrate[1/(x^2 - 1), x, Method -> \"RischMacsyma\"]"
        " - Integrate`RischMacsyma[1/(x^2 - 1), x]]",
        "0", 0);
}

/* ---------------- Strict / out-of-scope behaviour ---------------- */

static void test_strict_unevaluated(void) {
    /* Genuinely non-elementary integrands bubble the explicit package head
     * back unevaluated (never garbage). */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Sin[x^2], x]", "Integrate`RischMacsyma");  /* Fresnel */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[Sin[x]/x, x]", "Integrate`RischMacsyma");  /* SinIntegral */
    /* A non-symbol integration variable is rejected. */
    assert_head_unevaluated(
        "Integrate`RischMacsyma[1/x, x + 1]", "Integrate`RischMacsyma");
}

void test_integrate_risch_macsyma(void) {
    symtab_init();
    core_init();

    TEST(test_rational_case);
    TEST(test_logarithmic_case);
    TEST(test_exponential_case);
    TEST(test_fractional_case);
    TEST(test_hermite_case);
    TEST(test_hyperexponential_case);
    TEST(test_trig_frontend);
    TEST(test_multikernel_case);
    TEST(test_log_tower_case);
    TEST(test_exp_tower_case);
    TEST(test_recursive_tower_case);
    TEST(test_special_functions);
    TEST(test_cascade_default);
    TEST(test_method_plumbing);
    TEST(test_strict_unevaluated);

    printf("All Integrate RischMacsyma tests passed!\n");
}

int main(void) {
    test_integrate_risch_macsyma();
    return 0;
}
