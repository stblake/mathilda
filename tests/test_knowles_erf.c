/* test_knowles_erf.c — Knowles' error-function integration of transcendental
 * Liouvillian functions (KNOWLES_DESIGN.md, K2).
 *
 * The positive battery pins the paper examples (Part II §4, in Mathilda's classical
 * Erf convention: erf(u) = (Sqrt[Pi]/2) Erf[u], so Knowles' erf(erf(x)) is
 * (Pi/4) Erf[Erf[x]] here) plus systematically harder generated cases (deeper
 * nesting, mixed li/erf towers, varied constants).  Each asserts BOTH that the
 * result contains the special-function head AND that it diff-backs exactly to the
 * integrand.  The decision battery asserts the engine DECLINES soundly on
 * integrands with no erf-elementary antiderivative.
 *
 * These I/O pairs double as the worked-example corpus for a future
 * "Integration in terms of error functions" tutorial (KNOWLES_DESIGN.md §4).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

static int failures = 0;

static Expr* eval_str(const char* s) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static bool eval_is(const char* s, const char* expected) {
    Expr* r = eval_str(s);
    char* got = expr_to_string_fullform(r);
    bool ok = strcmp(got, expected) == 0;
    if (!ok) printf("  FAIL [%s] -> %s (expected %s)\n", s, got, expected);
    free(got); expr_free(r);
    return ok;
}

/* Positive pin: Integrate closes to a form containing `head` and diff-backs to 0. */
static void pos(const char* f, const char* head) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{r=Integrate[%s, x]}, Head[r]=!=Integrate && !FreeQ[r,%s]]", f, head);
    if (!eval_is(buf, "True")) { printf("    (no closed %s form for %s)\n", head, f); failures++; return; }
    snprintf(buf, sizeof(buf), "Simplify[D[Integrate[%s, x], x] - (%s)]", f, f);
    if (!eval_is(buf, "0")) { printf("    (diff-back nonzero for %s)\n", f); failures++; }
    else printf("  ok   INT %s   (-> %s)\n", f, head);
}

/* Exact-form pin: Integrate[f] prints exactly `expect` (FullForm). */
static void exact(const char* f, const char* expect) {
    char buf[1600];
    snprintf(buf, sizeof(buf), "Integrate[%s, x]", f);
    if (eval_is(buf, expect)) printf("  ok   INT %s  ==  %s\n", f, expect);
    else failures++;
}

/* Decision pin: Integrate declines (stays unevaluated). */
static void declines(const char* f) {
    char buf[1600];
    snprintf(buf, sizeof(buf), "Head[Integrate[%s, x]]", f);
    if (eval_is(buf, "Integrate")) printf("  ok   INT %s   declines (sound)\n", f);
    else { printf("  FAIL [%s]: expected decline\n", f); failures++; }
}

int main(void) {
    core_init();
    printf("Running test: knowles_erf (K2 erf-Liouvillian)\n");

    /* --- Part II §4 paper pins (Mathilda Erf convention) --------------------- */
    /* Ex 4.1: INT exp(-x^2 - erf^2 x) = erf(erf x)  ->  (Pi/4) Erf[Erf[x]]. */
    pos("E^(-x^2 - Erf[x]^2)", "Erf");
    /* Ex 4.3: elementary over the tower (no erf term) -> -(Sqrt[Pi]/4) E^(-Erf[x]^2). */
    pos("Erf[x] E^(-x^2 - Erf[x]^2)", "Erf");
    /* Ex 4.4: INT [2 e^-x^2 erf x - 3 e^-1/x^2 / x^2] = erf^2 x + 3 erf(1/x). */
    pos("2 E^(-x^2) Erf[x] - 3 E^(-1/x^2)/x^2", "Erf");

    /* --- Exact-form pins (flagship results) --------------------------------- */
    exact("E^(-x^2 - Erf[x]^2)", "Times[Rational[1, 4], Pi, Erf[Erf[x]]]");
    exact("Erf[x] E^(-x^2 - Erf[x]^2)",
          "Times[Rational[-1, 4], Power[Pi, Rational[1, 2]], Power[E, Times[-1, Power[Erf[x], 2]]]]");

    /* --- Generated harder cases (nesting / mixed / Erfi) --------------------- */
    /* Triple nesting: INT exp(-x^2 - erf^2 x - erf^2(erf x)) = (const) erf(erf(erf x)). */
    pos("E^(-x^2 - Erf[x]^2 - Erf[Erf[x]]^2)", "Erf");
    /* Scaled argument: INT e^{-4 x^2} = (Sqrt[Pi]/4) Erf[2 x]. */
    pos("E^(-4 x^2)", "Erf");
    /* Reciprocal case INT e^{-1/x^2}/x^2 = -(Sqrt[Pi]/2) Erf[1/x].  Regression guard
     * for the cherry_ei complex-r gate + the Q(i) Together fast path: this integrand
     * used to HANG the whole Integrate cascade (a Gaussian-rational Together blow-up
     * inside the ei engine's degenerate erf ansatz). */
    pos("E^(-1/x^2)/x^2", "Erf");

    /* --- Radical / quasiquadratic case (Part I Ex 8.1; 1986 Ex 3.2) --------- */
    /* INT exp(1/2 loglog x - 1/log x)/(x log^2 x) = -Sqrt[Pi] Erf[1/Sqrt[Log x]]
     * (Knowles: 2 erf(1/sqrt log x)).  The algebraic factor E^(1/2 loglog x) =
     * Sqrt[Log x] hides a half-integer power of Log[x], so the erf argument is a
     * RADICAL 1/Sqrt[Log x] — the engine solves it in s = Sqrt[Log x] (Log x -> s^2). */
    pos("E^(1/2 Log[Log[x]] - 1/Log[x])/(x Log[x]^2)", "Erf");
    /* Already-reduced spelling (explicit half-integer power). */
    pos("E^(-1/Log[x])/(x Log[x]^(3/2))", "Erf");
    exact("E^(-1/Log[x])/(x Log[x]^(3/2))",
          "Times[-1, Power[Pi, Rational[1, 2]], Erf[Power[Log[x], Rational[-1, 2]]]]");
    /* Erfi radical dual: +1/log x = u^2 -> Erfi[1/Sqrt[Log x]]. */
    pos("E^(1/Log[x])/(x Log[x]^(3/2))", "Erfi");
    exact("E^(1/Log[x])/(x Log[x]^(3/2))",
          "Times[-1, Power[Pi, Rational[1, 2]], Erfi[Power[Log[x], Rational[-1, 2]]]]");

    /* --- Decision battery (sound declines) ---------------------------------- */
    /* Ex 4.2: no erf-elementary antiderivative. */
    declines("x E^(-x^2 - Erf[x]^2)");
    /* A non-integrable variant. */
    declines("x^2 E^(-x^2 - Erf[x]^2)");

    if (failures) { printf("knowles_erf: %d FAIL\n", failures); return 1; }
    printf("knowles_erf: all passed\n");
    return 0;
}
