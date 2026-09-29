/*
 * bench_field_gcd.c
 * -----------------
 * Stress and scaling harness for flint_field_gcd -- the native multivariate GCD
 * over a number field K = Q(theta) in src/poly/flint_bridge.c.
 *
 * WHY THIS IS NOT A PLAIN TIMING BENCH.  The characteristic failure of a modular
 * GCD is not a wrong answer and not a slow one: it is a DECLINE.  When the engine
 * cannot certify a candidate within its prime budget it returns NULL, and
 * poly_gcd_internal's post-check then answers 1 -- a valid common divisor, so
 * nothing downstream complains, and the real gcd is simply lost.  A harness that
 * measured only time would have reported the pre-v0.232 engine as healthy while
 * it silently gave up on every gcd whose coefficients exceeded ~831 bits and on
 * every degree-6 radical field.  So every case here is SELF-CERTIFYING: the
 * operands are built as d*u and d*v with u, v coprime by construction, and the
 * answer must come back an associate of the planted d.  A decline is a failure,
 * not a slow row.
 *
 * The engine is called DIRECTLY rather than through PolynomialGCD.  That isolates
 * it from evaluator overhead, keeps the evaluator's result memo from making the
 * repeat trials vacuous, and distinguishes a decline (NULL) from a genuine gcd of
 * 1, which the surface spelling cannot.
 *
 * Verification is NUMERIC, by the ratio of the answer to the planted factor at
 * two rational points.  An exact check would have to divide polynomials over K,
 * which is both the thing under test and -- through PolynomialQuotient, which
 * cannot divide by a Root-coefficient polynomial -- a source of false negatives.
 * The ratio being the same nonzero constant at independent points says associate
 * without relying on either.
 *
 * Two gates, so this fails CI rather than merely reporting:
 *   1. every case certifies (correctness + no decline);
 *   2. the term-count doubling ratio stays under SCALE_MAX, which catches a
 *      return to O(L^2) image construction machine-independently.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "flint_bridge.h"

#define N_TRIALS   5

/* t(2n)/t(n): ~2 is linear, ~4 is quadratic.  Measured 2.33-2.44 across runs, and
 * a return to a keyed insert per term would read ~3.8, so the limit sits between
 * the two with enough headroom that scheduling noise on a loaded CI machine
 * cannot trip it.  A gate that fails at random stops being read. */
#define SCALE_MAX  3.2

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1.0e9 + (double)ts.tv_nsec;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static Expr* eval_str(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) { fprintf(stderr, "parse failed: %s\n", src); exit(2); }
    Expr* r = evaluate(p);
    expr_free(p);
    return r;
}

static void eval_discard(const char* src) { expr_free(eval_str(src)); }

static int eval_is_true(const char* src) {
    Expr* r = eval_str(src);
    char* s = expr_to_string(r);
    int t = (s && strcmp(s, "True") == 0);
    free(s);
    expr_free(r);
    return t;
}

/* Two independent rational points; unused variables in the rule list are inert,
 * so one list serves every arity. */
#define PT1 "{x -> 3, y -> 5, z -> 7, w -> 11, t -> 13}"
#define PT2 "{x -> 4, y -> 9, z -> 2, w -> 3,  t -> 5}"

/* Is `got` an associate of `want` -- i.e. got/want a nonzero constant?  Checked
 * at two points, at 30 digits, with a generous tolerance: the question is
 * constant-or-not, and the two ratios are either bit-comparable or differ by
 * orders of magnitude. */
static int is_associate(const char* got, const char* want) {
    /* Both operands appear TWICE in the template below, so the buffer must hold
     * two copies of each.  Sizing it for one was a heap overflow that only fired
     * on the large-coefficient rows, where a single coefficient is ~900 digits. */
    size_t need = 2 * (strlen(got) + strlen(want)) + 512;
    char* buf = malloc(need);
    if (!buf) return 0;
    snprintf(buf, need,
        "Module[{r1, r2}, "
        "r1 = N[((%s)/(%s)) /. %s, 30]; "
        "r2 = N[((%s)/(%s)) /. %s, 30]; "
        "NumericQ[r1] && NumericQ[r2] && Abs[r1] > 10^-15 && "
        "Abs[r1 - r2] < 10^-15 Max[1, Abs[r1]]]",
        got, want, PT1, got, want, PT2);
    int ok = eval_is_true(buf);
    free(buf);
    return ok;
}

/* One stress case: operands d*u and d*v, so the gcd must be an associate of d. */
typedef struct {
    const char* label;
    const char* d;      /* planted common factor  */
    const char* u;      /* cofactor of the first operand  */
    const char* v;      /* cofactor of the second         */
} Case;

/* Coverage rationale, case by case:
 *   - degrees 2,3,4,5,6,8,12: 6 is the one that regressed (a compositum whose
 *     membership relation needs >64 bits, denied escalation), and 12 is the
 *     largest field the qqbar layer builds cheaply;
 *   - Q(sqrt2, sqrt3) is NON-CYCLIC, Galois group (Z/2)^2 with no element of
 *     order 4, so NO prime leaves M irreducible -- it is the case that makes
 *     splitting the residue ring mandatory rather than an optimisation;
 *   - the nested radical and the Root spelling exercise the two surface forms
 *     that reach the engine by different routes;
 *   - Q(i) is the one field the classical path also handles (get_int_content
 *     understands Gaussian integers), so it must not regress;
 *   - the coefficient rows are the old ceiling: ~831 bits was the cliff, so 300
 *     bits must pass, 1000 must pass, and 3000 exercises well past it;
 *   - algebraic cofactors matter because Extension -> Automatic used to answer 1
 *     as soon as the COFACTORS carried algebraic constants, not the gcd. */
static Case CASES[] = {
  { "n=2   Q(sqrt2)",            "x^3 + Sqrt[2] x y + y^2 + 1",  "x^2 + Sqrt[2] y + 3",   "x^2 + 2 Sqrt[2] y - 1" },
  { "n=3   Q(2^1/3)",            "x^3 + 2^(1/3) x y + y^2 + 1",  "x^2 + 2^(1/3) y + 3",   "x^2 + 2 2^(1/3) y - 1" },
  { "n=4   Q(2^1/4)",            "x^3 + 2^(1/4) x y + y^2 + 1",  "x^2 + 2^(1/4) y + 3",   "x^2 + 2 2^(1/4) y - 1" },
  { "n=5   Q(2^1/5)",            "x^3 + 2^(1/5) x y + y^2 + 1",  "x^2 + 2^(1/5) y + 3",   "x^2 + 2 2^(1/5) y - 1" },
  { "n=6   Q(2^1/6)",            "x^3 + 2^(1/6) x y + y^2 + 1",  "x^2 + 2^(1/6) y + 3",   "x^2 + 2 2^(1/6) y - 1" },
  { "n=6   Q(sqrt2, 2^1/3)",     "x^3 + (Sqrt[2] + 2^(1/3)) x y + y^2 + 1", "x^2 + y + 3", "x^2 + 2 y - 1" },
  { "n=8   Q(2^1/8)",            "x^3 + 2^(1/8) x y + y^2 + 1",  "x^2 + 2^(1/8) y + 3",   "x^2 + 2 2^(1/8) y - 1" },
  { "n=12  Q(2^1/12)",           "x^3 + 2^(1/12) x y + y^2 + 1", "x^2 + y + 3",           "x^2 + 2 y - 1" },
  { "n=4   Q(sqrt2,sqrt3) split","x^3 + Sqrt[2] x y + y^2 + 1",  "x^2 + Sqrt[3] y + 3",   "x^2 + Sqrt[2] y - 1" },
  { "nested radical",            "x^2 + Sqrt[1 + Sqrt[2]] y + 1","x + y + 1",             "x - y + 2" },
  { "Root generator",            "x^2 + Root[-2 + #1^3 &, 1] y + 1", "x + y + 1",          "x - y + 2" },
  { "AlgebraicNumber spelling",  "x^2 + AlgebraicNumber[Sqrt[2], {0, 1}] y + 1", "x + y + 1", "x - y + 2" },
  { "Q(i) gaussian",             "x^2 + I y + 1",                "x + y + 1",             "x - y + 2" },
  { "3 variables",               "x^2 + Sqrt[2] y z + z + 1",    "x + y + z",             "x - y + 2 z" },
  { "4 variables",               "x + Sqrt[2] y + z + w + 1",    "x + y + z + w",         "x - y + 2 z - w" },
  { "5 variables",               "x + Sqrt[2] y + z + w + t",    "x + y + z + w + t",     "x - y + 2 z - w + t" },
  { "algebraic cofactors",       "x^3 + Sqrt[2] x y + y^2 + 1",  "x^2 + Sqrt[2] y + 3",   "x^2 + 2 Sqrt[2] y - 1" },
  { "one operand divides other", "x^3 + Sqrt[2] x y + y^2 + 1",  "1",                     "x + 2" },
  { "coeff ~300 bits",           "x^3 + 10^90 Sqrt[2] x y + y^2 + 1", "x^2 + Sqrt[2] y + 3", "x^2 + 2 Sqrt[2] y - 1" },
  { "coeff ~1000 bits",          "x^3 + 10^301 Sqrt[2] x y + y^2 + 1","x^2 + Sqrt[2] y + 3", "x^2 + 2 Sqrt[2] y - 1" },
  { "coeff ~3000 bits",          "x^3 + 10^903 Sqrt[2] x y + y^2 + 1","x^2 + Sqrt[2] y + 3", "x^2 + 2 Sqrt[2] y - 1" },
  { "high degree",               "x^8 + Sqrt[2] x^4 y^3 + y^6 + 1", "x^3 + Sqrt[2] y^2 + 3", "x^3 + 2 y^2 - 1" },
};
#define N_CASES ((int)(sizeof(CASES) / sizeof(CASES[0])))

/* Build Expand[d*u] and Expand[d*v], time flint_field_gcd over them, and check
 * the answer.  Returns 0 on success, 1 on a decline, 2 on a wrong answer. */
static int run_case(const Case* c, double* out_us, long* out_terms) {
    /* Sized from the strings rather than fixed: a truncated assignment would
     * silently change the problem being measured. */
    size_t need = strlen(c->d) + strlen(c->u) + strlen(c->v) + 256;
    char* buf = malloc(need);
    if (!buf) { fprintf(stderr, "out of memory\n"); exit(2); }
    snprintf(buf, need, "fgA = Expand[(%s) (%s)];", c->d, c->u);
    eval_discard(buf);
    snprintf(buf, need, "fgB = Expand[(%s) (%s)];", c->d, c->v);
    eval_discard(buf);
    free(buf);

    Expr* A = eval_str("fgA");
    Expr* B = eval_str("fgB");

    {   /* operand size, for the report */
        Expr* n = eval_str("If[Head[fgA] === Plus, Length[fgA], 1]");
        *out_terms = (n->type == EXPR_INTEGER) ? (long)n->data.integer : 0;
        expr_free(n);
    }

    /* One untimed warm-up: the first call in a process builds the shared prime
     * pool and the field cache, which belong to neither this case nor the
     * steady-state cost being measured. */
    Expr* warm = flint_field_gcd(A, B);
    int declined = (warm == NULL);
    char* got = NULL;
    if (warm) got = expr_to_string(warm);
    if (warm) expr_free(warm);

    double s[N_TRIALS];
    for (int i = 0; i < N_TRIALS; i++) {
        double t0 = now_ns();
        Expr* g = flint_field_gcd(A, B);
        double t1 = now_ns();
        if (g) expr_free(g);
        s[i] = (t1 - t0) / 1000.0;
    }
    qsort(s, N_TRIALS, sizeof(double), cmp_double);
    *out_us = s[N_TRIALS / 2];

    expr_free(A);
    expr_free(B);

    if (declined) return 1;
    int ok = is_associate(got, c->d);
    free(got);
    return ok ? 0 : 2;
}

/* Term-count scaling: the same problem shape at roughly n and 2n terms.  Uses
 * (x + y + sqrt2)^D, whose expansion has (D+1)(D+2)/2 terms, so D and
 * D*sqrt(2) roughly double it. */
static double scaling_us(int D, long* terms) {
    char buf[512];
    snprintf(buf, sizeof buf,
             "fgD = Expand[(x + y + Sqrt[2])^%d]; "
             "fgA = Expand[fgD (x^2 + Sqrt[2] y + 3)]; "
             "fgB = Expand[fgD (x^2 + 2 Sqrt[2] y - 1)];", D);
    eval_discard(buf);
    Expr* A = eval_str("fgA");
    Expr* B = eval_str("fgB");
    { Expr* n = eval_str("Length[fgA]");
      *terms = (n->type == EXPR_INTEGER) ? (long)n->data.integer : 0;
      expr_free(n); }

    Expr* w = flint_field_gcd(A, B);
    if (w) expr_free(w); else { fprintf(stderr, "scaling probe D=%d DECLINED\n", D); }

    double s[N_TRIALS];
    for (int i = 0; i < N_TRIALS; i++) {
        double t0 = now_ns();
        Expr* g = flint_field_gcd(A, B);
        double t1 = now_ns();
        if (g) expr_free(g);
        s[i] = (t1 - t0) / 1000.0;
    }
    qsort(s, N_TRIALS, sizeof(double), cmp_double);
    expr_free(A); expr_free(B);
    return s[N_TRIALS / 2];
}

int main(void) {
    symtab_init();
    core_init();

    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping field-gcd bench.\n");
        return 0;
    }

    printf("flint_field_gcd -- stress and scaling\n");
    printf("%-30s %10s %12s  %s\n", "case", "terms", "median(us)", "result");
    printf("--------------------------------------------------------------------\n");

    int declines = 0, wrong = 0;
    for (int i = 0; i < N_CASES; i++) {
        double us = 0; long terms = 0;
        int rc = run_case(&CASES[i], &us, &terms);
        const char* verdict = (rc == 0) ? "ok" : (rc == 1) ? "DECLINED" : "WRONG";
        if (rc == 1) declines++;
        if (rc == 2) wrong++;
        printf("%-30s %10ld %12.1f  %s\n", CASES[i].label, terms, us, verdict);
    }

    printf("\nterm-count scaling\n");
    long t1 = 0, t2 = 0;
    double u1 = scaling_us(22, &t1);
    double u2 = scaling_us(32, &t2);
    double ratio = (u1 > 0) ? u2 / u1 : 0.0;
    printf("  %ld terms -> %.1f us\n", t1, u1);
    printf("  %ld terms -> %.1f us\n", t2, u2);
    printf("  ratio %.2f for %.2fx the terms (limit %.2f)\n",
           ratio, t1 ? (double)t2 / (double)t1 : 0.0, SCALE_MAX);

    int fail = 0;
    if (declines || wrong) {
        printf("\nFAIL: %d declined, %d wrong\n", declines, wrong);
        fail = 1;
    }
    if (ratio > SCALE_MAX) {
        printf("\nFAIL: scaling ratio %.2f exceeds %.2f -- image construction may "
               "have returned to a keyed insert per term\n", ratio, SCALE_MAX);
        fail = 1;
    }
    if (!fail) printf("\nOK: %d cases certified, scaling within limit\n", N_CASES);
    return fail;
}
