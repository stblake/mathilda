/* test_trager_torsion_stress.c -- nightly stress for the Trager method's hard
 * cases (torsion orders over large residue fields, high genus, residue fields
 * built as towers).  Nightly by name ("stress"): the whole run takes minutes.
 *
 *   1. Every FriCAS-stress case (tests/data/trager_stress_corpus.m) that timed
 *      out or took more than 10 s at v0.6083, with its verdict and a time
 *      budget of several times its v0.6084 time: a certificate stays a
 *      certificate, an antiderivative is verified, and a regression to the
 *      former 60-240 s fails the budget.
 *   2. Seeded random integrands d/dx Log[p + q g^(1/m)] (m = 2..6): elementary
 *      by construction, so the method must return a verified antiderivative --
 *      a certificate or a decline would be a wrong answer or a lost one.  A
 *      correctness check: the budget (1500 s) is generous because for m = 5, 6
 *      the answers reach 10^5 leaves and take 350 - 700 s (v0.6086), almost all
 *      of it in the generator over the residue field in characteristic 0.
 *   3. Residue fields with no degree limit: the degree-36 splitting field of a
 *      sextic of the shape Z^6 + a Z^3 + b, built by the algebraic tower.
 */
#include "trager_test_utils.h"

/* f must get `verdict` ("not elementary" / "failed" prefix) within `budget` s. */
static void st_verdict(const char* f, const char* verdict, int budget, const char* file, int line) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{tt, rr}, {tt, rr} = AbsoluteTiming[TimeConstrained[TragerIntegrate[%s, x], %d, $Aborted]];"
             " If[ListQ[rr] && rr[[1]] === \"%s\" && tt < %d, True, {tt, rr}]]",
             f, budget, verdict, budget);
    printf("  %-70s ", f); fflush(stdout);
    tr_eq(buf, "True", file, line);
    printf("ok\n");
}
#define ST_NONELEM(f, b) st_verdict((f), "not elementary", (b), __FILE__, __LINE__)

/* f must integrate within `budget` s to a verified antiderivative. */
static void st_elem(const char* f, int budget, const char* file, int line) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{tt, rr, d}, {tt, rr} = AbsoluteTiming[TimeConstrained[TragerIntegrate[%s, x], %d, $Aborted]];"
             " If[!ListQ[rr] && rr =!= $Aborted && tt < %d &&"
             "   (d = D[rr, x] - (%s); TrueQ[PossibleZeroQ[d]] || numericZeroQ[d, x]), True, {tt, rr}]]",
             f, budget, budget, f);
    printf("  %-70s ", f); fflush(stdout);
    tr_eq(buf, "True", file, line);
    printf("ok\n");
}
#define ST_ELEM(f, b) st_elem((f), (b), __FILE__, __LINE__)

/* ---------------------------------------- 1. the former slow cases */
static void test_corpus_hard_cases(void) {
    /* certificates: residue divisors of infinite order (orders mod p, q
     * disagree, or the torsion screen refutes the only possible order) */
    ST_NONELEM("x^3/(x - 2)*(x^2 + 1)^(1/6)", 60);                         /* ratrad-002 */
    ST_NONELEM("1/(x + 1)*(x^3 - 2*x)^(1/4)", 60);                         /* ratrad-007 */
    ST_NONELEM("(2*x - 1)/(x + 1)*(x^2 - x + 1)^(1/5)", 60);               /* ratrad-010 */
    ST_NONELEM("x/(x - 2)*(x^3 + 1)^(-2/3)", 60);                          /* ratrad-012 */
    ST_NONELEM("(x^2 + 1)/x*(x^4 + 3)^(1/6)", 60);                         /* ratrad-016 */
    ST_NONELEM("1/(x^2 + 2)*(x^3 + 1)^(-1/4)", 60);                        /* ratrad-019 */
    ST_NONELEM("x/(x + 1)*(x^2 + 1)^(1/5)", 60);                           /* ratrad-031 */
    ST_NONELEM("(2*x - 1)/(x + 1)*(x^2 + 1)^(2/5)", 60);                   /* ratrad-032 */
    ST_NONELEM("1/(x + 1)*(x^3 - 2*x)^(-1/4)", 60);                        /* ratrad-037 */
    ST_NONELEM("(2*x - 1)/(x - 2)*(x^2 + 1)^(-3/5)", 60);                  /* ratrad-049 */
    ST_NONELEM("1/((x^2 + 1)*(x^4 + 1)^(1/4))", 60);                       /* goursat-014 */
    ST_NONELEM("(x^2 - 1)/((x^2 + 1)*(x^4 + 1)^(1/4))", 60);               /* goursat-015 */
    ST_NONELEM("1/((x^2 - 2)*(x^4 - 4*x^2 + 1)^(1/4))", 60);               /* goursat-021 */
    ST_NONELEM("(3*x + 2)/(-2*x + 3)*(x^4 + 2*x^2 - x + 3)^(-3/4)", 60);   /* random-003 */
    ST_NONELEM("(3*x - 3)/(2*x + 1)*(x^2 + x)^(-5/6)", 60);                /* random-012 */
    ST_NONELEM("(-x^2 - 3*x - 3)/(2*x + 1)*(x^4 + 3*x^3 + 3*x^2 - 3*x - 2)^(-5/6)", 120); /* random-037 */
    ST_NONELEM("(-1)/(-3*x + 1)*(-2*x^4 - 2*x^3 + x^2 - x - 3)^(-3/4)", 60);/* random-038 */
    ST_NONELEM("1/((x - 2)*(x^3 + 1)^(1/3))", 60);                         /* order-21 screen */
    ST_NONELEM("(x^2 + 1)/x*(x^4 + 2)^(1/6)", 60);                         /* torsion-009: genus 7, Weil bound > 2^63 */
    ST_NONELEM("(-3*x^2 - 3*x)/(-2*x^2 + 2*x - 1)*(-2*x^4 + 3*x^3 - 3*x^2 + 3*x + 2)^(-5/6)", 120); /* random-010: trace functional over a degree-12 field */
    /* elementary, formerly 1.3 - 4.2 s */
    ST_ELEM("1/(x^4 + 1)^(1/4)", 60);                                      /* binom-032 */
    ST_ELEM("x^(-1)*(2*x^4 - 1)^(2/3)", 60);                               /* cheb-020 */
    ST_ELEM("1/x*(x^2 + 1)^(-1/6)", 60);                                   /* ratrad-045 */
    ST_ELEM("((2*x + 1)/(x + 3))^(1/6)", 60);                              /* ratradicand-009 */
    ST_ELEM("(x^2 + 2)/((x^2 + 1)*(x^3 + x)^(1/3))", 60);                  /* goursat-020 */
}

/* ------------------------------------- 2. seeded elementary integrands */
static void test_seeded_logs(void) {
    TR_TRUE("SeedRandom[20261030];"
            "stRp[d_] := Sum[RandomInteger[{-3, 3}] x^k, {k, 0, d}];"
            "stBad = {};"
            "Do[Module[{m = RandomChoice[{2, 3, 4, 5, 6}], g, p, q, f, tt, rr, d},"
            "   g = stRp[RandomInteger[{2, 4}]]; If[Exponent[g, x] < 1, g = x^2 + 1];"
            "   p = stRp[RandomInteger[{1, 3}]]; q = stRp[RandomInteger[{0, 2}]]; If[q === 0, q = 1];"
            "   f = Together[D[Log[p + q g^(1/m)], x]];"
            "   {tt, rr} = AbsoluteTiming[TimeConstrained[TragerIntegrate[f, x], 1500, $Aborted]];"
            "   Print[\"  seeded \", k, \" m=\", m, \" \", Round[tt, 0.01], \" s\"];"
            "   If[ListQ[rr] || rr === $Aborted ||"
            "      !(d = D[rr, x] - f; TrueQ[PossibleZeroQ[d]] || numericZeroQ[d, x]),"
            "     AppendTo[stBad, {k, m, g, p, q, tt, If[ListQ[rr], rr, Head[rr]]}]]],"
            " {k, 40}];"
            "If[stBad === {}, True, stBad]");
}

/* ------------------------------------------- 3. large residue fields */
static void test_large_residue_field(void) {
    /* the sextic of ratrad-048's remainder: splitting field of degree 36 */
    TR_TRUE("Module[{R = Expand[(2187 - 832 Z^3 + 512 Z^6)/512], r},"
            " r = residueField[R, Z, t];"
            " fieldDegree[r[\"F\"]] === 36 && Length[r[\"roots\"]] === 6 &&"
            " AllTrue[r[\"roots\"], kc[r[\"F\"], R /. Z -> #] === 0 &] && DuplicateFreeQ[r[\"roots\"]]]");
}

int main(void) {
    alarm(0);                              /* nightly: no 120 s unit watchdog */
    setvbuf(stdout, NULL, _IOLBF, 0);
    tr_init();
    test_corpus_hard_cases();
    test_seeded_logs();
    test_large_residue_field();
    printf("trager_torsion_stress: %d assertions passed\n", tr_assertions);
    return 0;
}
