/* test_trager_integrate.c -- end to end: src/internal/trager/Pipeline.m and
 * the C glue in src/calculus/integrate.c (try_trager /
 * integrate_trager_common): Integrate[f, x, Method -> "Trager" |
 * "TragerKauers"] and the qualified symbols Integrate`Trager /
 * Integrate`TragerKauers.
 *
 * Pipeline functions covered: integrateComponent, integrateCurve,
 * rationalIntegral, componentExpression, reduceRadicalPowers, fieldValue,
 * tragerLogExpressions, TragerIntegrateCurve, TragerIntegrate,
 * tragerIntegrateBody, spotCheck, numericZeroQ, antiderivativeQ, presentable.
 *
 * Correctness of every antiderivative is asserted by differentiation
 * (TR_ANTIDERIV: PossibleZeroQ[D[F, x] - f]).  For the answers on which
 * Mathilda's PossibleZeroQ gives a FALSE NEGATIVE (nested radicals such as
 * (1 + x^4)^(1/4) inside ArcTan/Log; all verified zero to 40+ digits), the
 * check is numeric instead (TI_ANTIDERIV_NUM: |D[F, x] - f| < 10^-30 at four
 * rational points with 60-digit arithmetic).
 *
 * Non-elementary verdicts are taken from theory, never from the package:
 * Chebyshev's theorem on binomial differentials x^m (a + b x^n)^p (elementary
 * iff p, (m+1)/n or (m+1)/n + p is an integer); [B] Ex. 5 (x/Sqrt[1-x^3]);
 * elliptic/hyperelliptic differentials of the first kind (1/Sqrt[cubic],
 * 1/Sqrt[x^5 - x], ...); residue divisors of infinite order (Nagell-Lutz).
 * Each was also cross-checked against FriCAS 1.3.6 (unevaluated there).
 */

/* dup/dup2/fileno (stderr capture for Integrate::nonelem) are POSIX. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "trager_test_utils.h"
#include "message.h"
#include <unistd.h>

/* Numeric antiderivative check (PossibleZeroQ false negatives). */
static void ti_antideriv_num(const char* f, const char* method, const char* file, int line) {
    char buf[8192];
    snprintf(buf, sizeof buf,
             "Module[{trF = Integrate[%s, x, Method -> \"%s\"], trD},"
             " trD = D[trF, x] - (%s);"
             " FreeQ[trF, Integrate] && AllTrue[{37/23, 113/71, 29/7, 3/2},"
             "   TrueQ[Abs[N[trD /. x -> #, 60]] < 10^-30] &]]",
             f, method, f);
    tr_eq(buf, "True", file, line);
}
#define TI_ANTIDERIV_NUM(f, m) ti_antideriv_num((f), (m), __FILE__, __LINE__)

/* Integrate stays unevaluated (and keeps its Method option). */
static void ti_uneval(const char* f, const char* method, const char* file, int line) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "MatchQ[Quiet[Integrate[%s, x, Method -> \"%s\"]], HoldPattern[Integrate[_, x, Method -> \"%s\"]]]",
             f, method, method);
    tr_eq(buf, "True", file, line);
}
#define TI_UNEVAL(f, m) ti_uneval((f), (m), __FILE__, __LINE__)

/* Integrate`Trager / TragerKauers returns {tag, reason} with this tag. */
static void ti_verdict(const char* sym, const char* f, const char* tag, const char* file, int line) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "MatchQ[Quiet[%s[%s, x]], {\"%s\", _String}]", sym, f, tag);
    tr_eq(buf, "True", file, line);
}
#define TI_NONELEM(f) ti_verdict("Integrate`Trager", (f), "not elementary", __FILE__, __LINE__)
#define TI_NONELEM_K(f) ti_verdict("Integrate`TragerKauers", (f), "not elementary", __FILE__, __LINE__)
#define TI_FAILED(f) ti_verdict("Integrate`Trager", (f), "failed", __FILE__, __LINE__)
#define TI_FAILED_K(f) ti_verdict("Integrate`TragerKauers", (f), "failed", __FILE__, __LINE__)

/* Evaluate `input` with stderr redirected to a temp file and report whether
 * `needle` was printed (copied from tests/test_parallelmixedtower.c). */
static int eval_stderr_contains(const char* input, const char* needle) {
    fflush(stderr);
    int saved = dup(fileno(stderr));
    if (saved == -1) { fprintf(stderr, "FAIL: dup\n"); exit(1); }
    FILE* cap = tmpfile();
    if (!cap) { fprintf(stderr, "FAIL: tmpfile\n"); exit(1); }
    dup2(fileno(cap), fileno(stderr));

    Expr* parsed = parse_expression(input);
    Expr* evaluated = parsed ? evaluate(parsed) : NULL;
    if (parsed) expr_free(parsed);
    if (evaluated) expr_free(evaluated);

    fflush(stderr);
    dup2(saved, fileno(stderr));
    close(saved);
    if (!parsed) { fprintf(stderr, "FAIL (parse): %s\n", input); exit(1); }

    rewind(cap);
    char buf[16384];
    size_t n = fread(buf, 1, sizeof(buf) - 1, cap);
    buf[n] = '\0';
    fclose(cap);
    return strstr(buf, needle) != NULL;
}
static void ti_stderr(const char* input, const char* needle, int want, const char* file, int line) {
    int got = eval_stderr_contains(input, needle);
    tr_assertions++;
    if (got != want) {
        fprintf(stderr, "FAIL %s:%d\n  input: %s\n  expected stderr %s \"%s\"\n",
                file, line, input, want ? "to contain" : "NOT to contain", needle);
        exit(1);
    }
}
#define TI_PRINTS(in, needle) ti_stderr((in), (needle), 1, __FILE__, __LINE__)
#define TI_SILENT(in, needle) ti_stderr((in), (needle), 0, __FILE__, __LINE__)

/* ================================================================ internals */

static void test_integrate_component(void) {
    /* [B] Ex. 6: Int dx/(x Sqrt[1-x^3]) = 1/3 Log[(y-1)/(y+1)] (component c Y,
     * c = 1/(x (1-x^3))) */
    TR_TRUE("Module[{r0 = integrateComponent[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Trager\"], F0},"
            " F0 = r0[\"Algebraic\"] + Total[#[[1]] Log[#[[2]]] & /@ r0[\"Logs\"]] /. r0[\"Y\"] -> Sqrt[1 - x^3];"
            " r0[\"Algebraic\"] === 0 && Length[r0[\"Logs\"]] == 1 && r0[\"Logs\"][[1, 1]] === 1/3 &&"
            " TrueQ[PossibleZeroQ[D[F0, x] - 1/(x Sqrt[1 - x^3])]]]");
    TR_TRUE("Sort[Keys[integrateComponent[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Trager\"]]] === {\"Algebraic\", \"Logs\", \"Y\"}");
    /* [B] Ex. 3/4: algebraic part -y/(4 (x^8+1)), log part 1/4 Log[x^4/(1+y)] */
    TR_TRUE("Module[{r0 = integrateComponent[(2 x^8 + 1)/(x^17 + 2 x^9 + x), x, 2, x^8 + 1, \"Trager\"], Y0},"
            " Y0 = r0[\"Y\"];"
            " TrueQ[PossibleZeroQ[r0[\"Algebraic\"] + Y0/(4 (x^8 + 1))]] && Length[r0[\"Logs\"]] == 1 && r0[\"Logs\"][[1, 1]] === 1/4]");
    /* Kauers engine on the same component */
    TR_TRUE("Module[{r0 = integrateComponent[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Kauers\"], F0},"
            " F0 = r0[\"Algebraic\"] + Total[#[[1]] Log[#[[2]]] & /@ r0[\"Logs\"]] /. r0[\"Y\"] -> Sqrt[1 - x^3];"
            " Length[r0[\"Logs\"]] == 2 && TrueQ[PossibleZeroQ[D[F0, x] - 1/(x Sqrt[1 - x^3])]]]");
    /* a shifted component (pole at infinity before the shift): Sqrt[x^2+1] */
    TR_TRUE("Module[{r0 = integrateComponent[1, x, 2, x^2 + 1, \"Trager\"], F0},"
            " F0 = r0[\"Algebraic\"] + Total[#[[1]] Log[#[[2]]] & /@ r0[\"Logs\"]] /. r0[\"Y\"] -> Sqrt[x^2 + 1];"
            " FreeQ[r0, z] && TrueQ[PossibleZeroQ[D[F0, x] - Sqrt[x^2 + 1]]]]");
    /* purely algebraic: x^2/Sqrt[x^3+1] = x^2 Y/(x^3+1), integral 2/3 Y */
    TR_TRUE("Module[{r0 = integrateComponent[x^2/(x^3 + 1), x, 2, x^3 + 1, \"Trager\"]},"
            " r0[\"Logs\"] === {} && Together[r0[\"Algebraic\"] - 2/3 r0[\"Y\"]] === 0]");
    /* a RootSum term from Kauers on a degree-2 residue factor */
    TR_TRUE("Module[{r0 = integrateComponent[1/(x^2 - 2), x, 2, x^2 + 1, \"Kauers\"]},"
            " MemberQ[Length /@ r0[\"Logs\"], 3]]");
    /* certificates */
    TR_TRUE("MatchQ[catchTrager[integrateComponent[x/(1 - x^3), x, 2, 1 - x^3, \"Trager\"]],"
            " tragerFail[\"NonElementary\", \"the Hermite remainder has a pole at infinity\", _]]");   /* [B] Ex. 5 */
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1/(x^3 + 1), x, 2, x^3 + 1, \"Trager\"]],"
            " tragerFail[\"NonElementary\", s_String /; StringContainsQ[s, \"no residues\"], _]]");
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1/((x - 2) (x^3 + 1)), x, 2, x^3 + 1, \"Trager\"]],"
            " tragerFail[\"NonElementary\", s_String /; StringContainsQ[s, \"holomorphic\"], _]]");
    /* Kauers declines instead of certifying */
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1/((x - 2) (x^3 + 1)), x, 2, x^3 + 1, \"Kauers\"]],"
            " tragerFail[\"Decline\", _, _]]");
    /* invariants */
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1, x, 2, x^2, \"Trager\"]], tragerFail[\"Internal\", \"CurveNotReduced\", _]]");
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1, x, 1, x + 1, \"Trager\"]], tragerFail[\"Internal\", \"CurveNotReduced\", _]]");
    TR_TRUE("MatchQ[catchTrager[integrateComponent[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Foo\"]], tragerFail[\"Internal\", \"UnknownMethod\", _]]");
}

static void test_integrate_curve_rational(void) {
    TR_EQ("integrateCurve[1/(x^2 + 1), x, 1, 1, \"Trager\"][\"Algebraic\"]", "ArcTan[x]");
    TR_EQ("integrateCurve[1/(x^2 + 1), x, 1, 1, \"Trager\"][\"Logs\"]", "{}");
    TR_TRUE("Module[{r0 = integrateCurve[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Trager\"]},"
            " r0[\"Logs\"][[1, 1]] === 1/3]");
    TR_EQ("rationalIntegral[0, x]", "0");
    TR_EQ("rationalIntegral[x - x, x]", "0");
    TR_EQ("rationalIntegral[3 x^2, x]", "x^3");
    TR_TRUE("TrueQ[PossibleZeroQ[D[rationalIntegral[1/(x^2 - 1), x], x] - 1/(x^2 - 1)]]");
    TR_TRUE("TrueQ[PossibleZeroQ[D[rationalIntegral[(x^3 + 1)/(x^2 + x + 1)^2, x], x] - (x^3 + 1)/(x^2 + x + 1)^2]]");
}

static void test_assembly_helpers(void) {
    /* fieldValue */
    TR_EQ("fieldValue[fieldQ[], 3/2]", "3/2");
    TR_EQ("fieldValue[fieldNF[tiT, tiT^2 - 2, Sqrt[2]], 1 + tiT]", "1 + Sqrt[2]");
    TR_EQ("fieldValue[fieldNF[tiT, tiT^2 - 2, Sqrt[2]], tiT^3]", "2 Sqrt[2]");   /* reduced mod mu first */
    TR_EQ("fieldValue[fieldNF[tiT, tiT^2 + 1, I], 2 tiT - 1]", "-1 + 2*I");
    /* componentExpression: generator substituted, tidy applied to logands */
    TR_EQ("componentExpression[<|\"Algebraic\" -> 2 tiY, \"Logs\" -> {{1/2, tiY + 1}}, \"Y\" -> tiY|>, Sqrt[x], Identity]",
          "2 Sqrt[x] + 1/2 Log[1 + Sqrt[x]]");
    TR_EQ("componentExpression[<|\"Algebraic\" -> 0, \"Logs\" -> {{1, RootSum[#^3 - 2 &, Log[tiY - #] &], \"RootSum\"}}, \"Y\" -> tiY|>, x^(1/3), Identity]",
          "RootSum[#1^3 - 2 &, Log[x^(1/3) - #1] &]");
    TR_TRUE("Module[{e = componentExpression[<|\"Algebraic\" -> 0, \"Logs\" -> {{1, 2 tiY + 4}}, \"Y\" -> tiY|>, Sqrt[x], dropNumericContent]},"
            " MatchQ[e, Log[_]] && FreeQ[Together[e[[1]]/(2 + Sqrt[x])], x]]");
    TR_EQ("componentExpression[<|\"Algebraic\" -> tiY^2, \"Logs\" -> {}, \"Y\" -> tiY|>, Sqrt[x], Identity]", "x");
    /* reduceRadicalPowers: yh^2 = x */
    TR_TRUE("Module[{e = reduceRadicalPowers[tiy^5/(tiy^3 + 1), tiy, 2, x]},"
            " Exponent[Numerator[e], tiy] < 2 && Exponent[Denominator[e], tiy] < 2 &&"
            " Together[e - x^2 tiy/(x tiy + 1)] === 0]");
    TR_TRUE("Module[{e = reduceRadicalPowers[tiy^7 + x tiy^4, tiy, 3, x^2 + 1]},"
            " Exponent[e, tiy] < 3 && Expand[e - ((x^2 + 1)^2 tiy + x (x^2 + 1) tiy)] === 0]");
    TR_EQ("reduceRadicalPowers[x/2, tiy, 2, x]", "1/2 x");
    /* tragerLogExpressions: Ex. 6, expressions in x and Y */
    TR_TRUE("Module[{B0 = curveBasis[2, 1 - x^3, x], lg, r0},"
            " r0 = integrateComponent[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Trager\"];"
            " lg = tragerLogExpressions[fieldQ[], B0, {{1/3, {(2 (-2 + x^3))/x^3, 4/x^3}}}, tiY];"
            " lg[[1, 1]] === 1/3 && TrueQ[PossibleZeroQ[D[lg[[1, 1]] Log[lg[[1, 2]]] /. tiY -> Sqrt[1 - x^3], x] - 1/(x Sqrt[1 - x^3])]]]");
    TR_TRUE("Module[{lg = tragerLogExpressions[fieldNF[tiT, tiT^2 - 6, Sqrt[6]], curveBasis[2, x^2 + 1, x],"
            " {{tiT/4, {tiT x, 1}}}, tiY]}, lg === {{Sqrt[6]/4, Sqrt[6] x + tiY}}]");
}

static void test_checks(void) {
    /* numericZeroQ */
    TR_TRUE("numericZeroQ[0, x]");
    TR_TRUE("numericZeroQ[10^-40, x]");
    TR_FALSE("numericZeroQ[10^-20, x]");
    TR_FALSE("numericZeroQ[Sqrt[x^2] - x, x]");                    /* nonzero at negative points */
    TR_TRUE("numericZeroQ[((1 + x^4)^(1/4))^4 - 1 - x^4, x]");
    TR_TRUE("numericZeroQ[Sqrt[2 + Sqrt[3]] - (Sqrt[6] + Sqrt[2])/2, x]");   /* denesting identity */
    TR_TRUE("numericZeroQ[If[x > 1, 0, Indeterminate], x]");       /* exactly 4 numeric points */
    TR_FALSE("numericZeroQ[If[x > 2, 0, Indeterminate], x]");      /* only 2 */
    TR_FALSE("numericZeroQ[Sin[x] - x, x]");
    /* antiderivativeQ */
    TR_TRUE("antiderivativeQ[1/Sqrt[x^2 + 1], Log[x + Sqrt[x^2 + 1]], x]");
    TR_TRUE("antiderivativeQ[1/Sqrt[x^2 + 1], ArcSinh[x] + 7, x]");
    TR_FALSE("antiderivativeQ[1/Sqrt[x^2 + 1], Log[x + Sqrt[x^2 + 1]] + x/10^20, x]");
    TR_FALSE("antiderivativeQ[1/Sqrt[x^2 + 1], Log[x - Sqrt[x^2 + 1]], x]");
    /* the (1 + x^4)^(-1/4) answer: numericZeroQ rescues PossibleZeroQ */
    TR_TRUE("Module[{f0 = (x^4 + 1)^(-1/4), F0}, F0 = TragerIntegrate[f0, x, \"Trager\"];"
            " !ListQ[F0] && antiderivativeQ[f0, F0, x]]");
    /* spotCheck: True or an Internal failure, never False */
    TR_TRUE("spotCheck[x/Sqrt[x^2 + 1], Sqrt[x^2 + 1], x]");
    TR_TRUE("MatchQ[catchTrager[spotCheck[x/Sqrt[x^2 + 1], 2 Sqrt[x^2 + 1], x]],"
            " tragerFail[\"Internal\", \"AntiderivativeCheckFailed\", _]]");
    /* presentable: a verified real form replaces the complex one */
    TR_TRUE("Module[{G0 = presentable[1/(x^2 + 1), I/2 Log[(x + I)/(x - I)], x]},"
            " FreeQ[G0, Complex] && TrueQ[PossibleZeroQ[D[G0, x] - 1/(x^2 + 1)]]]");
    TR_EQ("presentable[1/x, Log[2 x], x]", "Log[x]");
    TR_EQ("presentable[2 x, x^2, x]", "x^2");
    /* the real form is kept only when it re-verifies: Re(d/dx I Log[x + 2 I])
     * = 2/(x^2 + 4) is not 1/(x^2 + 1), so F comes back unchanged ... */
    TR_TRUE("Module[{F0 = I Log[x + 2 I]}, presentable[1/(x^2 + 1), F0, x] === F0]");
    /* ... while I Log[x + I] has real part -ArcTan[1/x], an antiderivative */
    TR_TRUE("Module[{G0 = presentable[1/(x^2 + 1), I Log[x + I], x]}, FreeQ[G0, Complex] && TrueQ[PossibleZeroQ[D[G0, x] - 1/(x^2 + 1)]]]");
}

static void test_public_entry_points(void) {
    TR_EQ("TragerIntegrate[1/Sqrt[x^2 + 1], x]", "Log[x + Sqrt[1 + x^2]]");   /* default method */
    TR_EQ("TragerIntegrate[1/Sqrt[x^2 + 1], x, \"Trager\"]", "Log[x + Sqrt[1 + x^2]]");
    TR_TRUE("TrueQ[PossibleZeroQ[D[TragerIntegrate[1/Sqrt[x^2 + 1], x, \"Kauers\"], x] - 1/Sqrt[x^2 + 1]]]");
    TR_TRUE("TragerIntegrate[1/Sqrt[x^2 + 1], x, \"Foo\"] === {\"failed\", \"unknown method Foo\"}");
    TR_TRUE("TragerIntegrate[1/(x^2 + 1), x] === {\"failed\", \"no radical in the integrand (a rational function)\"}");
    TR_TRUE("TragerIntegrate[1/Sqrt[x^3 + 1], x] === {\"not elementary\", \"the Hermite remainder is a nonzero holomorphic differential (no residues)\"}");
    TR_TRUE("TragerIntegrate[x/Sqrt[1 - x^3], x] === {\"not elementary\", \"the Hermite remainder has a pole at infinity\"}");
    TR_TRUE("MatchQ[TragerIntegrate[(x^4 + 1)^(-1/4), x, \"Kauers\"], {\"failed\", _String}]");
    TR_TRUE("MatchQ[TragerIntegrate[Sqrt[x] + Sqrt[x + 1], x], {\"failed\", _String}]");
    /* tragerIntegrateBody throws instead */
    TR_TRUE("MatchQ[catchTrager[tragerIntegrateBody[1/Sqrt[x^3 + 1], x, \"Trager\"]], tragerFail[\"NonElementary\", _, _]]");
    TR_TRUE("MatchQ[catchTrager[tragerIntegrateBody[1/Sqrt[x^2 + 1], x, \"Bar\"]], tragerFail[\"Decline\", \"unknown method Bar\", _]]");
    TR_EQ("catchTrager[tragerIntegrateBody[x/Sqrt[x^2 + 1], x, \"Trager\"]]", "Sqrt[1 + x^2]");
    /* a rational part plus a radical part */
    TR_TRUE("Module[{f0 = 1/(x^2 + 1) + 1/Sqrt[x^2 + 1], F0}, F0 = TragerIntegrate[f0, x];"
            " TrueQ[PossibleZeroQ[D[F0, x] - f0]]]");
    /* TragerIntegrateCurve */
    TR_TRUE("Module[{r0 = TragerIntegrateCurve[1/(x (1 - x^3)), x, 2, 1 - x^3]}, AssociationQ[r0] && r0[\"Logs\"][[1, 1]] === 1/3]");
    TR_TRUE("Module[{r0 = TragerIntegrateCurve[1/(x (1 - x^3)), x, 2, 1 - x^3, \"Kauers\"]}, AssociationQ[r0] && Length[r0[\"Logs\"]] == 2]");
    TR_TRUE("TragerIntegrateCurve[x/(1 - x^3), x, 2, 1 - x^3] === {\"not elementary\", \"the Hermite remainder has a pole at infinity\"}");
    TR_TRUE("TragerIntegrateCurve[1, x, 2, x^2] === {\"failed\", \"internal: CurveNotReduced\"}");
    TR_TRUE("MatchQ[TragerIntegrateCurve[1/(x^4 + 1), x, 4, (x^4 + 1)^3, \"Kauers\"], {\"failed\", _String}]");
}

/* ================================================================ end to end */

static const char* const trager_ok[] = {
    /* quadratic radicands, genus 0 */
    "1/Sqrt[x^2 + 1]", "1/Sqrt[x^2 - 1]", "1/Sqrt[1 - x^2]", "x/Sqrt[x^2 + 1]", "x/Sqrt[1 - x^2]",
    "1/(x Sqrt[x^2 + 1])", "1/(x Sqrt[x^2 - 1])", "1/(x Sqrt[1 - x^2])", "Sqrt[x^2 + 1]/(x^2 - 2)",
    "Sqrt[x^2 - 4 x - 1]/(x^2 + 1)", "Sqrt[x^2 + 1]", "Sqrt[1 - x^2]", "x^3 Sqrt[x^2 + 1]",
    "Sqrt[x^2 + 2 x + 5]", "1/Sqrt[x^2 + 2 x + 5]", "1/(x^2 Sqrt[x^2 + 1])", "1/((x + 1) Sqrt[x^2 + 1])",
    "1/((x^2 + 1) Sqrt[x^2 + 2])", "x/((x^2 + 3) Sqrt[x^2 + 1])", "(x^2 + x) Sqrt[x^2 + 1]",
    "1/(x^2 + 1)^(3/2)", "Sqrt[x]", "1/Sqrt[x]", "Sqrt[2 x + 3]/x", "x^2/Sqrt[3 - 2 x - x^2]",
    "1/((2 x + 1) Sqrt[x^2 + x + 1])", "1/(x^2 + 1) + 1/Sqrt[x^2 + 1]",
    /* monomial radicands, several powers of one radicand */
    "x^(-1/3)", "x^(-2/3)", "x^(1/3)", "x^(3/4)", "x^(-3/4)", "(x^3)^(-1/4)", "(x^2)^(1/3)",
    "x^(1/2) + x^(1/3)", "x^(1/2) + x^(1/4)", "1/(1 + x^(1/3))", "x^(1/3)/(1 + x)",
    /* rescale: rational radicands, square factors */
    "Sqrt[(x + 1)/(x - 1)]", "1/Sqrt[(x + 1)/(x - 1)]", "Sqrt[(1 - x)/(1 + x)]",
    "Sqrt[x^2 (x^2 + 1)]", "((x^2 + 1)^2)^(1/4)",
    /* purely algebraic parts, polynomial times radical, higher genus */
    "x^2/Sqrt[x^3 + 1]", "x^2 Sqrt[x^3 + 1]", "(3 x^2 + 1)/Sqrt[x^3 + x]", "x^2/(x^3 + 1)^(1/3)",
    "x (x^2 + 1)^(1/3)", "x^3/(x^4 + 1)^(1/4)", "x^4/Sqrt[x^5 + 1]", "(5 x^4 + 1)/(x^5 + x + 3)^(2/3)",
    /* Bronstein Ex. 3, 4, 6 and relatives (positive genus, torsion > 1) */
    "(2 x^8 + 1) Sqrt[x^8 + 1]/(x^17 + 2 x^9 + x)", "Sqrt[x^8 + 1]/(x (x^8 + 1))", "1/(x Sqrt[1 - x^3])",
    "1/(x Sqrt[x^3 + 1])", "1/(x Sqrt[x^5 + 1])", "1/(x Sqrt[x^6 + 1])",
    /* Chebyshev/Goursat-type pseudo-elliptic integrals */
    "(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1])", "(x^2 + 1)/((x^2 - 1) Sqrt[x^4 + 1])",
    "(1 + x^2)/((1 - x^2) Sqrt[1 + x^4])", "x/((x^2 - 1) Sqrt[x^4 + 1])", "(x^2 - 1)/(x Sqrt[x^4 + 1])",
    "x/Sqrt[x^4 + 1]", "(x^2 + 1)/((x^2 - 1) Sqrt[x^4 + 3 x^2 + 1])", "(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 3 x^2 + 1])",
    "(x^2 + 1)/((x^2 - 1) Sqrt[x^4 - x^2 + 1])", "x/((x^2 + 2) Sqrt[x^4 + 1])",
    "(2 x + 1)/((x^2 + x) Sqrt[(x^2 + x)^2 + 1])", "x^2/((x^3 + 2) Sqrt[x^3 + 1])",
    "x/Sqrt[x^4 + 10 x^2 - 96 x - 71]",
    /* higher-index radicals with logarithmic parts */
    "1/(x (x^3 + 1)^(1/3))", "1/(x (x^6 + 1)^(1/3))", "1/(x (x^5 + 1)^(1/3))", "3 x^2/((x^3 - 1) (x^3 + 1)^(1/3))",
    /* exact derivatives D[v + c Log[u]] with u in Q(x, y): elementary by construction */
    "D[(x + 1) (x^2 - x + 2)^(2/3)/(x - 3) + 2/3 Log[x + (x^2 - x + 2)^(1/3)], x]",
    "D[(x^2 + 1) Sqrt[x^3 - 2 x + 2]/(x - 1) + 1/2 Log[x - 2 + x Sqrt[x^3 - 2 x + 2]], x]",
    NULL
};

/* Same answers, but Mathilda's PossibleZeroQ cannot see that D[F] - f is 0
 * (checked numerically; see the file header). */
static const char* const trager_ok_numeric[] = {
    "(x^4 + 1)^(-1/4)", "1/(x^3 + 1)^(1/3)", "x/(x^3 + 1)^(2/3)", "1/(x^3 - 1)^(1/3)",
    "1/(x^4 - 1)^(1/4)", "1/(x^3 + x)^(1/3)", "1/(x (x^4 + 2)^(1/4))", "1/(x (x^4 + 1)^(1/4))",
    NULL
};

static const char* const kauers_ok[] = {
    "1/Sqrt[x^2 + 1]", "1/Sqrt[x^2 - 1]", "1/Sqrt[1 - x^2]", "x/Sqrt[x^2 + 1]", "1/(x Sqrt[x^2 + 1])",
    "1/(x Sqrt[x^2 - 1])", "1/(x Sqrt[1 - x^2])", "Sqrt[x^2 + 1]/(x^2 - 2)", "Sqrt[x^2 - 4 x - 1]/(x^2 + 1)",
    "x^(-1/3)", "x^(1/2) + x^(1/3)", "x^(3/4)", "1/(x Sqrt[1 - x^3])", "1/(x Sqrt[x^8 + 1])",
    "(2 x^8 + 1) Sqrt[x^8 + 1]/(x^17 + 2 x^9 + x)", "Sqrt[(x + 1)/(x - 1)]", "Sqrt[(1 - x)/(1 + x)]",
    "Sqrt[x^2 (x^2 + 1)]", "((x^2 + 1)^2)^(1/4)", "x^2/Sqrt[x^3 + 1]", "x^2 Sqrt[x^3 + 1]",
    "(3 x^2 + 1)/Sqrt[x^3 + x]", "Sqrt[x^2 + 1]", "Sqrt[1 - x^2]", "Sqrt[x^2 + 2 x + 5]",
    "1/((x + 1) Sqrt[x^2 + 1])", "1/((x^2 + 1) Sqrt[x^2 + 2])", "x/((x^2 + 3) Sqrt[x^2 + 1])",
    "x/(x^3 + 1)^(2/3)", "x/Sqrt[x^4 + 1]", "(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1])",
    "(1 + x^2)/((1 - x^2) Sqrt[1 + x^4])", "x/((x^2 - 1) Sqrt[x^4 + 1])", "(x^2 - 1)/(x Sqrt[x^4 + 1])",
    "(x^2 + 1)/((x^2 - 1) Sqrt[x^4 + 3 x^2 + 1])", "1/(x Sqrt[x^6 + 1])", "1/(x Sqrt[x^5 + 1])",
    "x^2/((x^3 + 2) Sqrt[x^3 + 1])", "(2 x + 1)/((x^2 + x) Sqrt[(x^2 + x)^2 + 1])", "x/((x^2 + 2) Sqrt[x^4 + 1])",
    NULL
};

/* Not elementary (proof in the comment), certified by both engines. */
static const char* const nonelem[] = {
    "1/Sqrt[x^3 + 1]",                /* first kind, genus 1 */
    "x/Sqrt[1 - x^3]",                /* [B] Ex. 5: second kind */
    "1/Sqrt[x^4 + 1]",                /* first kind, genus 1 */
    "1/Sqrt[x^5 - x]",                /* first kind, genus 2 */
    "1/Sqrt[x^3 - x]",
    "Sqrt[x^3 + 1]",
    "Sqrt[1 - x^4]",
    "x^2/Sqrt[x^4 + 1]",              /* Chebyshev: m=2,n=4,p=-1/2: 3/4, 1/4 */
    "x/(x^4 + 1)^(1/4)",              /* m=1,n=4,p=-1/4: 1/2, 1/4 */
    "x/(1 - x^4)^(1/4)",
    "1/(x^2 (x^4 + 1)^(1/4))",        /* m=-2: -1/4, -1/2 */
    "1/(x^4 - 1)^(3/4)",              /* m=0,p=-3/4: 1/4, -1/2 */
    "1/(x^4 + 1)^(1/3)",              /* 1/4, -1/12 */
    "1/(x^2 + 1)^(1/3)",              /* 1/2, 1/6 */
    "(x^3 + 1)^(1/4)",                /* 1/3, 7/12 */
    "1/(x^5 + 1)^(1/2)",
    "1/((x + 1) Sqrt[x^3 + 1])",      /* = 1/3 (x+1)/y dx - d(2y/(3(x+1))): first+second kind */
    "Sqrt[x^3 + 1]/(x - 2)",
    NULL
};

/* Not elementary, residue divisor of infinite order or torsion with a
 * nonzero holomorphic remainder: certified by Trager; Kauers only declines. */
static const char* const nonelem_trager_only[] = {
    "5/((x - 3) Sqrt[x^3 - 2])",      /* (3, 5) on y^2 = x^3 - 2: Nagell-Lutz, infinite order */
    "1/((x - 3) Sqrt[x^3 - 2])",
    "1/((x - 1) Sqrt[x^3 + 8])",      /* (1, 3): 2P has x = -7/4 */
    "1/((x - 1) Sqrt[x^3 + 1])",      /* (1, Sqrt[2]): infinite order */
    "1/((2 x - 1) Sqrt[x^3 + 1])",
    "1/((x - 2) Sqrt[x^3 + 1])",      /* (2, 3) has order 6, but the holomorphic part is nonzero (FriCAS agrees) */
    "1/((x^2 + 1) Sqrt[x^4 + 4])",
    "1/((x^2 - 2) Sqrt[x^4 + 1])",
    NULL
};

/* Out of scope: declined silently (no certificate, no message). */
static const char* const declines[] = {
    "Sqrt[x] + Sqrt[x + 1]",          /* compositum */
    "Exp[x]/Sqrt[x]",
    "Sqrt[2] x^(1/2)",                /* algebraic constant */
    "Sin[x] Sqrt[x]",
    "Sqrt[x^2 + a]",                  /* parameter */
    "x^Pi",
    "x^2",                            /* rational: no radical */
    "1/(x^2 + 1)",
    "Log[x] Sqrt[x]",
    "Sqrt[x] + (x + 1)^(1/3)",
    NULL
};

/* Kauers' heuristic does not close (elementary, but no principal power of
 * the residue ideal found): {"failed", ...} and unevaluated. */
static const char* const kauers_fails[] = {
    "(x^4 + 1)^(-1/4)", "1/(x^3 + 1)^(1/3)", "1/(x (x^3 + 1)^(1/3))", "1/(1 + x^(1/3))",
    "x/Sqrt[x^4 + 10 x^2 - 96 x - 71]",
    NULL
};

static void test_end_to_end(void) {
    int i;
    for (i = 0; trager_ok[i]; i++) TR_ANTIDERIV(trager_ok[i], "Trager");
    for (i = 0; trager_ok_numeric[i]; i++) TI_ANTIDERIV_NUM(trager_ok_numeric[i], "Trager");
    for (i = 0; kauers_ok[i]; i++) TR_ANTIDERIV(kauers_ok[i], "TragerKauers");
    for (i = 0; nonelem[i]; i++) {
        TI_NONELEM(nonelem[i]);
        TI_NONELEM_K(nonelem[i]);
        TI_UNEVAL(nonelem[i], "Trager");
        TI_UNEVAL(nonelem[i], "TragerKauers");
    }
    for (i = 0; nonelem_trager_only[i]; i++) {
        TI_NONELEM(nonelem_trager_only[i]);
        TI_UNEVAL(nonelem_trager_only[i], "Trager");
        TI_FAILED_K(nonelem_trager_only[i]);
        TI_UNEVAL(nonelem_trager_only[i], "TragerKauers");
    }
    for (i = 0; declines[i]; i++) {
        TI_FAILED(declines[i]);
        TI_FAILED_K(declines[i]);
        TI_UNEVAL(declines[i], "Trager");
    }
    for (i = 0; kauers_fails[i]; i++) {
        TI_FAILED_K(kauers_fails[i]);
        TI_UNEVAL(kauers_fails[i], "TragerKauers");
    }
    /* specific answers, derived by hand */
    TR_EQ("Integrate[1/Sqrt[x^2 + 1], x, Method -> \"Trager\"]", "Log[x + Sqrt[1 + x^2]]");
    TR_EQ("Integrate[x/Sqrt[x^2 + 1], x, Method -> \"Trager\"]", "Sqrt[1 + x^2]");
    TR_EQ("Integrate[x^(-2/3), x, Method -> \"Trager\"]", "3 x^(1/3)");
    TR_EQ("Integrate[x^(1/2) + x^(1/3), x, Method -> \"TragerKauers\"]", "3/4 x^(4/3) + 2/3 x^(3/2)");
    /* [B] Ex. 6: exactly one logarithm with coefficient 1/3 */
    TR_TRUE("Module[{F0 = Integrate[1/(x Sqrt[1 - x^3]), x, Method -> \"Trager\"]}, MatchQ[F0, Rational[1, 3] Log[_]]]");
    /* [B] Ex. 4: 1/4 Log[...] */
    TR_TRUE("Module[{F0 = Integrate[1/(x Sqrt[x^8 + 1]), x, Method -> \"Trager\"]}, MatchQ[F0, Rational[1, 4] Log[_]]]");
    /* answers are real (ArcTan, no I) for real integrands where possible */
    TR_TRUE("FreeQ[Integrate[1/Sqrt[1 - x^2], x, Method -> \"Trager\"], Complex]");
    TR_TRUE("FreeQ[Integrate[(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1]), x, Method -> \"Trager\"], Complex]");
    TR_TRUE("FreeQ[Integrate[(x^4 + 1)^(-1/4), x, Method -> \"Trager\"], Complex]");
    /* the variable need not be x */
    TR_EQ("Integrate[1/Sqrt[t^2 + 1], t, Method -> \"Trager\"]", "Log[t + Sqrt[1 + t^2]]");
    TR_TRUE("TrueQ[PossibleZeroQ[D[Integrate[1/(s Sqrt[1 - s^3]), s, Method -> \"TragerKauers\"], s] - 1/(s Sqrt[1 - s^3])]]");
    /* the two engines agree up to a constant on Ex. 6 */
    TR_TRUE("TrueQ[PossibleZeroQ[D[Integrate[1/(x Sqrt[1 - x^3]), x, Method -> \"Trager\"] -"
            " Integrate[1/(x Sqrt[1 - x^3]), x, Method -> \"TragerKauers\"], x]]]");
}

static void test_qualified_symbols(void) {
    TR_TRUE("Integrate`Trager[1/Sqrt[x^3 + 1], x] === {\"not elementary\", \"the Hermite remainder is a nonzero holomorphic differential (no residues)\"}");
    TR_TRUE("Integrate`TragerKauers[1/Sqrt[x^3 + 1], x] === {\"not elementary\", \"the Hermite remainder is a nonzero holomorphic differential (no residues)\"}");
    TR_TRUE("Integrate`TragerKauers[(x^4 + 1)^(-1/4), x] === {\"failed\", \"Kauers' heuristic: no principal power of the residue ideal\"}");
    TR_TRUE("Integrate`Trager[Sqrt[x] + Sqrt[x + 1], x] === {\"failed\", \"more than one radicand (a radical compositum)\"}");
    TR_TRUE("Integrate`Trager[x^2, x] === {\"failed\", \"no radical in the integrand (a rational function)\"}");
    TR_TRUE("MatchQ[Integrate`Trager[5/((x - 3) Sqrt[x^3 - 2]), x], {\"not elementary\", s_String /; StringContainsQ[s, \"infinite order\"]}]");
    TR_EQ("Integrate`Trager[1/Sqrt[x^2 + 1], x]", "Log[x + Sqrt[1 + x^2]]");
    TR_EQ("Integrate`TragerKauers[x/Sqrt[x^2 + 1], x]", "Sqrt[1 + x^2]");
    /* wrong arity: unevaluated */
    TR_TRUE("Head[Integrate`Trager[1/Sqrt[x^2 + 1]]] === Integrate`Trager");
    TR_TRUE("Head[Integrate`TragerKauers[1/Sqrt[x^2 + 1], x, y]] === Integrate`TragerKauers");
    TR_TRUE("MemberQ[Attributes[Integrate`Trager], Protected] && MemberQ[Attributes[Integrate`TragerKauers], Protected]");
}

static void test_definite_and_options(void) {
    /* definite integrals pass through the indefinite method:
     * Int_0^1 dx/Sqrt[x^2+1] = ArcSinh[1] = Log[1 + Sqrt[2]] */
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[1/Sqrt[x^2 + 1], {x, 0, 1}, Method -> \"Trager\"] - Log[1 + Sqrt[2]]]]");
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[1/Sqrt[x^2 + 1], {x, 0, 1}, Method -> \"TragerKauers\"] - Log[1 + Sqrt[2]]]]");
    /* Int_0^1 x/Sqrt[x^2+1] = Sqrt[2] - 1 */
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[x/Sqrt[x^2 + 1], {x, 0, 1}, Method -> \"Trager\"] - (Sqrt[2] - 1)]]");
    /* Int_0^1 Sqrt[1 - x^2] = Pi/4 */
    TR_TRUE("TrueQ[Abs[N[Integrate[Sqrt[1 - x^2], {x, 0, 1}, Method -> \"Trager\"] - Pi/4, 30]] < 10^-25]");
    /* Int_1^8 x^(-2/3) = 3 (2 - 1) = 3 */
    TR_EQ("Integrate[x^(-2/3), {x, 1, 8}, Method -> \"Trager\"]", "3");
    /* a non-elementary definite integral stays unevaluated */
    TR_TRUE("MatchQ[Quiet[Integrate[1/Sqrt[x^3 + 1], {x, 0, 1}, Method -> \"Trager\"]], HoldPattern[Integrate[_, {x, 0, 1}, ___]]]");
    /* an unknown method is still an error, and unevaluated */
    TR_TRUE("MatchQ[Quiet[Integrate[1/Sqrt[x^2 + 1], x, Method -> \"Tragerr\"]], HoldPattern[Integrate[_, x, Method -> \"Tragerr\"]]]");
    TR_EQ("Check[Integrate[1/Sqrt[x^2 + 1], x, Method -> \"Bogus\"], tiCaught]", "tiCaught");
}

/* Messages: run with the suite-wide suppression lifted so they print. */
static void test_messages(void) {
    mth_msg_suppress_pop();
    TI_PRINTS("Integrate[1/Sqrt[x^3 + 1], x, Method -> \"Trager\"]", "Integrate::nonelem");
    TI_PRINTS("Integrate[x/Sqrt[1 - x^3], x, Method -> \"TragerKauers\"]", "Integrate::nonelem");
    TI_PRINTS("Integrate`Trager[1/Sqrt[x^4 + 1], x]", "Integrate::nonelem");
    TI_PRINTS("Integrate[5/((x - 3) Sqrt[x^3 - 2]), x, Method -> \"Trager\"]", "Integrate::nonelem");
    TI_PRINTS("Integrate[1/Sqrt[x^3 + 1], {x, 0, 1}, Method -> \"Trager\"]", "Integrate::nonelem");
    /* Quiet silences it */
    TI_SILENT("Quiet[Integrate[1/Sqrt[x^3 + 1], x, Method -> \"Trager\"]]", "nonelem");
    TI_SILENT("Quiet[Integrate[1/Sqrt[x^3 + 1], x, Method -> \"Trager\"], Integrate::nonelem]", "nonelem");
    /* declines are silent, and so is Kauers' "failed" */
    TI_SILENT("Integrate[Sqrt[x] + Sqrt[x + 1], x, Method -> \"Trager\"]", "Integrate::");
    TI_SILENT("Integrate[Exp[x]/Sqrt[x], x, Method -> \"Trager\"]", "Integrate::");
    TI_SILENT("Integrate[(x^4 + 1)^(-1/4), x, Method -> \"TragerKauers\"]", "Integrate::");
    TI_SILENT("Integrate[5/((x - 3) Sqrt[x^3 - 2]), x, Method -> \"TragerKauers\"]", "Integrate::");
    TI_SILENT("Integrate[1/Sqrt[x^2 + 1], x, Method -> \"Trager\"]", "Integrate::");
    /* the internal stages' own messages never leak */
    TI_SILENT("Integrate[(x^4 + 1)^(-1/4), x, Method -> \"Trager\"]", "::");
    TI_SILENT("Integrate[Sqrt[x^2 - 4 x - 1]/(x^2 + 1), x, Method -> \"Trager\"]", "::");
    /* an unknown method prints Integrate::method */
    TI_PRINTS("Integrate[1/Sqrt[x^2 + 7], x, Method -> \"Bogus7\"]", "Integrate::method");
    mth_msg_suppress_push();
    /* Check catches the certificate, not a decline (a Check entered inside a
     * quiet region still sees the messages fired in it) */
    TR_EQ("Check[Integrate[1/Sqrt[x^3 + 1], x, Method -> \"Trager\"], tiCaught]", "tiCaught");
    TR_EQ("Check[Integrate[x/Sqrt[1 - x^3], x, Method -> \"TragerKauers\"], tiCaught]", "tiCaught");
    TR_TRUE("Check[Integrate[Sqrt[x] + Sqrt[x + 1], x, Method -> \"Trager\"], tiCaught] =!= tiCaught");
    TR_TRUE("Check[Integrate[(x^4 + 1)^(-1/4), x, Method -> \"TragerKauers\"], tiCaught] =!= tiCaught");
    TR_EQ("Check[Integrate[x/Sqrt[x^2 + 1], x, Method -> \"Trager\"], tiCaught]", "Sqrt[1 + x^2]");
    TR_EQ("Quiet[Check[Integrate[1/Sqrt[x^3 + 1], x, Method -> \"Trager\"], tiCaught]]", "tiCaught");
}

/* ================================================================ known failures
 * Genuine bugs, with the CORRECT expectation; run only when TRAGER_KNOWN_FAIL
 * is set in the environment. */
static void known_failures_trager_integrate(void) {
    /* BUG 2 (> 60 s).  Int dx/((x - 2) (x^3 + 1)^(1/3)): a genus-1 curve with
     * residues in a degree-6 field; no verdict after 400 s (torsion order 21
     * found mod p in seconds, then the principality test of 21 delta over the
     * sextic residue field does not finish). */
    TR_TRUE("Module[{r0 = TimeConstrained[Integrate`Trager[1/((x - 2) (x^3 + 1)^(1/3)), x], 60, tiTimeout]},"
            " r0 =!= tiTimeout && !MatchQ[r0, {\"failed\", _}]]");
    /* BUG 3 (> 60 s on elementary integrands of index 5 and 6).  Each is
     * D[Log[u]] for u in Q(x, y), so elementary by construction; the first is
     * on the GENUS-0 curve y^6 = x + 1 (rational by t = (x+1)^(1/6)).  The
     * time goes into the residue-ideal arithmetic / torsion search for m = 6
     * (7 of 40 random index-5/6 exact derivatives also time out). */
    TR_TRUE("Module[{f0 = D[Log[1 + x (x + 1)^(1/6)], x], r0},"
            " r0 = TimeConstrained[Integrate`Trager[f0, x], 60, tiTimeout];"
            " r0 =!= tiTimeout && !ListQ[r0] && TrueQ[PossibleZeroQ[D[r0, x] - f0]]]");
    TR_TRUE("Module[{f0 = D[Log[1 + (x^3 + 1)^(1/6)], x], r0},"
            " r0 = TimeConstrained[Integrate`Trager[f0, x], 60, tiTimeout];"
            " r0 =!= tiTimeout && !ListQ[r0]]");
    TR_TRUE("Module[{f0 = D[Log[x + (x^4 + 1)^(1/5)], x], r0},"
            " r0 = TimeConstrained[Integrate`Trager[f0, x], 60, tiTimeout];"
            " r0 =!= tiTimeout && !ListQ[r0]]");
}

/* Formerly known failures, now fixed.
 * - Definite integrals: the real form chooses the ArcTan orientation without a
 *   jump inside the domain, and Newton-Leibniz now cross-checks ArcTan of a
 *   quotient (integrate_newton_leibniz.c nl_has_arctan_of_quotient).
 * - Residue fields that ToNumberField could not build: the per-orbit
 *   construction works over Q[Z]/(q) (FriCAS stress cases random-027/030/053). */
static void test_fixed_known_failures(void) {
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[1/Sqrt[1 - x^2], {x, -1/2, 1/2}, Method -> \"Trager\"] - Pi/3]]");
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[Sqrt[1 - x^2], {x, -1, 1}, Method -> \"Trager\"] - Pi/2]]");
    TR_TRUE("TrueQ[PossibleZeroQ[Integrate[1/Sqrt[1 - x^2], {x, -1/2, 1/2}, Method -> \"TragerKauers\"] - Pi/3]]");
    TR_EQ("Integrate[1/Sqrt[1 - x^2], x, Method -> \"Trager\"]", "ArcTan[x/Sqrt[1 - x^2]]");
    TR_ANTIDERIV("(-2)/(-x^2 + 2*x - 3)*(3*x + 2)^(2/3)", "Trager");
    TR_ANTIDERIV("(-3*x^2 + 3*x - 3)/(2*x^2 - x - 1)*(x + 3)^(1/4)", "Trager");
    TR_ANTIDERIV("(-3*x + 2)/(2*x - 2)*(3*x)^(-4/5)", "Trager");
    TR_ANTIDERIV("((2*x + 1)/(x + 3))^(1/6)", "Trager");
}

/* The slowdowns closed at v0.6084 (FLINT` polynomial-matrix kernels, the
 * torsion screen, the algebraic residue-field tower, no search caps).  Each
 * verdict comes with a time budget well above its measured time (0.2 - 5 s),
 * so a regression to the former 60 - 240 s shows. */
static void tr_timed_verdict(const char* f, const char* verdict, int seconds, const char* file, int line) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{tt, rr}, {tt, rr} = AbsoluteTiming[TimeConstrained[TragerIntegrate[%s, x], %d, $Aborted]];"
             " ListQ[rr] && rr[[1]] === \"%s\" && tt < %d]", f, seconds, verdict, seconds);
    tr_eq(buf, "True", file, line);
}
#define TR_TIMED_VERDICT(f, v, s) tr_timed_verdict((f), (v), (s), __FILE__, __LINE__)

static void test_former_slowdowns(void) {
    /* genus 1, cubic residue field, predicted order 21: the class is refuted at
     * a third good prime (torsionScreen) instead of 21 delta being formed over
     * the residue field (was > 240 s) */
    TR_TIMED_VERDICT("1/((x - 2) (x^3 + 1)^(1/3))", "not elementary", 20);
    TR_TIMED_VERDICT("x/((x - 2) (x^3 + 1)^(2/3))", "not elementary", 20);
    /* ToNumberField gave up on the octic residue polynomial: the field of
     * degree 16 is built as a tower (was a decline) */
    TR_TIMED_VERDICT("1/((x^2 + 2) (x^3 + 1)^(1/4))", "not elementary", 40);
    /* GF(p) torsion searches at primes 150 - 430 (were 60 - 110 s or timeouts) */
    TR_TIMED_VERDICT("x^3/(x - 2) (x^2 + 1)^(1/6)", "not elementary", 30);
    TR_TIMED_VERDICT("1/(x + 1) (x^3 - 2 x)^(1/4)", "not elementary", 30);
    TR_TIMED_VERDICT("(x^2 - 1)/((x^2 + 1) (x^4 + 1)^(1/4))", "not elementary", 30);
    TR_TIMED_VERDICT("(2 x - 1)/(x - 2) (x^2 + 1)^(-3/5)", "not elementary", 30);
    /* elementary, formerly 1.4 - 4.2 s: verified antiderivatives */
    TR_ANTIDERIV("1/(x^4 + 1)^(1/4)", "Trager");
    TR_ANTIDERIV("1/x*(x^2 + 1)^(-1/6)", "Trager");
    TR_ANTIDERIV("x^(-1)*(2*x^4 - 1)^(2/3)", "Trager");
    TR_ANTIDERIV("x^3*(x^3 - 2*x)^(-2/3)", "Trager");
    TR_ANTIDERIV("(x^2 + 2)/((x^2 + 1)*(x^3 + x)^(1/3))", "Trager");
}

int main(void) {
    tr_init();
    test_integrate_component();
    test_integrate_curve_rational();
    test_assembly_helpers();
    test_checks();
    test_public_entry_points();
    test_end_to_end();
    test_qualified_symbols();
    test_definite_and_options();
    test_messages();
    test_fixed_known_failures();
    test_former_slowdowns();
    if (getenv("TRAGER_KNOWN_FAIL")) known_failures_trager_integrate();
    printf("trager_integrate: %d assertions passed\n", tr_assertions);
    return 0;
}
