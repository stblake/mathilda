
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_trig_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[0]", "0"},
        {"Sin[Pi/2]", "1"},
        {"Sin[Pi/6]", "Rational[1, 2]"},
        {"Sin[Pi/4]", "Power[2, Rational[-1, 2]]"},
        {"Sin[Pi/3]", "Times[Rational[1, 2], Power[3, Rational[1, 2]]]"},
        {"Cos[0]", "1"},
        {"Cos[Pi/2]", "0"},
        {"Cos[Pi/3]", "Rational[1, 2]"},
        {"Tan[Pi/4]", "1"},
        {"Cot[Pi/4]", "1"},
        {"Sec[0]", "1"},
        {"Csc[Pi/2]", "1"},
        // Exact values for d=12
        {"Sin[Pi/12]", "Times[Rational[1, 4], Plus[Times[-1, Power[2, Rational[1, 2]]], Power[6, Rational[1, 2]]]]"},
        {"Cos[Pi/12]", "Times[Rational[1, 4], Plus[Power[2, Rational[1, 2]], Power[6, Rational[1, 2]]]]"},
        // Exact values for d=10, 5
        {"Sin[Pi/10]", "Times[Rational[1, 4], Plus[-1, Power[5, Rational[1, 2]]]]"},
        {"Cos[Pi/5]", "Times[Rational[1, 4], Plus[1, Power[5, Rational[1, 2]]]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_inverse() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"ArcSin[0]", "0"},
        {"ArcSin[1]", "Times[Rational[1, 2], Pi]"},
        {"ArcSin[1/2]", "Times[Rational[1, 6], Pi]"},
        {"ArcCos[0]", "Times[Rational[1, 2], Pi]"},
        {"ArcCos[1]", "0"},
        {"ArcTan[0]", "0"},
        {"ArcTan[1]", "Times[Rational[1, 4], Pi]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_forward_of_inverse() {
    /* f[f_inv[x]] = x for each direct trig / inverse-trig pair. Uses
     * InputForm-style strings since the simplification strips the inverse
     * and leaves the inner argument untouched. */
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[ArcSin[x]]", "x"},
        {"Cos[ArcCos[y]]", "y"},
        {"Tan[ArcTan[a]]", "a"},
        {"Cot[ArcCot[b]]", "b"},
        {"Sec[ArcSec[c]]", "c"},
        {"Csc[ArcCsc[d]]", "d"},
        /* Composite arguments survive as-is */
        {"Sin[ArcSin[x^2 + 1]]", "1 + x^2"},
        {"Cos[ArcCos[2 y]]", "2 y"},
        /* ArcTan[x, y] two-arg form must NOT be stripped by the
         * single-argument rule: symbolic Tan[ArcTan[x, y]] stays.  With exact
         * x > 0, ArcTan[3, 4] is ArcTan[4/3] (ArcTan#2, as the reference), and
         * then Tan of it is 4/3. */
        {"Tan[ArcTan[x, y]]", "Tan[ArcTan[x, y]]"},
        {"Tan[ArcTan[3, 4]]", "4/3"},
        /* Opposite direction is NOT folded: ArcSin[Sin[x]] stays put */
        {"ArcSin[Sin[x]]", "ArcSin[Sin[x]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0,
                   "Forward-of-inverse %s: expected %s, got %s",
                   cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

/* Phase 4: Sin / Cos / Tan / Cot / Sec / Csc on Complex[MPFR, MPFR] must
 * preserve MPFR precision rather than coerce to double via csin/ccos/etc.
 *
 * The asserted prefixes are the first ~40 digits of the closed-form
 * value computed at 50-digit precision; precision is verified via a
 * Precision[...] round-trip (which itself produces an MPFR with the
 * expected magnitude, ~50.27 decimal digits for a 50-digit-input
 * result). */
void test_trig_mpfr_complex(void) {
    struct {
        const char* input;
        const char* expected_prefix;
        size_t prefix_len;
    } cases[] = {
        /* sin(1+i) = sin(1)cosh(1) + i cos(1)sinh(1)
         *          ~ 1.29846 + 0.63496 i */
        {"Sin[Complex[N[1, 50], N[1, 50]]]",
         "1.29845758141597729482604236580781562031343656163522", 40},
        /* cos(1+i) = cos(1)cosh(1) - i sin(1)sinh(1) */
        {"Cos[Complex[N[1, 50], N[1, 50]]]",
         "0.83373002513114904888388539433509447980987478520962", 40},
        /* tan(1+i) ~ 0.27175 + 1.08392 i */
        {"Tan[Complex[N[1, 50], N[1, 50]]]",
         "0.27175258531951171652884372249858892070946411146178", 40},
        /* cot(1+i) = 1/tan(1+i) ~ 0.21762 - 0.86801 i */
        {"Cot[Complex[N[1, 50], N[1, 50]]]",
         "0.21762156185440268136513424360523807352075436916785", 40},
        /* sec(1+i) = 1/cos(1+i) ~ 0.49834 + 0.59108 i */
        {"Sec[Complex[N[1, 50], N[1, 50]]]",
         "0.49833703055518678521380589177216953443287793247109", 40},
        /* csc(1+i) = 1/sin(1+i) ~ 0.62152 - 0.30393 i */
        {"Csc[Complex[N[1, 50], N[1, 50]]]",
         "0.62151801717042842123490780585592014816751214181073", 40},
        {NULL, NULL, 0}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* r = evaluate(e);
        char* s = expr_to_string(r);
        ASSERT_MSG(strncmp(s, cases[i].expected_prefix, cases[i].prefix_len) == 0,
                   "%s: expected leading %s, got %s",
                   cases[i].input, cases[i].expected_prefix, s);
        free(s);
        expr_free(e);
        expr_free(r);
    }

    /* Spot-check precision is preserved for Sin (representative). */
    Expr* ep = parse_expression(
        "Precision[Sin[Complex[N[1, 50], N[1, 50]]]]");
    Expr* rp = evaluate(ep);
    char* sp = expr_to_string(rp);
    ASSERT_MSG(strncmp(sp, "50.", 3) == 0,
               "Precision[Sin[1+i]] (50 digits): expected 50.*, got %s", sp);
    free(sp); expr_free(ep); expr_free(rp);

    /* Real MPFR input still takes the real path (no I in result). */
    Expr* eR = parse_expression("Sin[N[1, 50]]");
    Expr* rR = evaluate(eR);
    char* sR = expr_to_string(rR);
    ASSERT_MSG(strchr(sR, 'I') == NULL,
               "Sin[N[1, 50]]: expected pure real, got %s", sR);
    ASSERT_MSG(strncmp(sR, "0.8414709848078965066525023216302989996225",
                       42) == 0,
               "Sin[N[1, 50]] (50 digits): expected leading 0.84147..., got %s",
               sR);
    free(sR); expr_free(eR); expr_free(rR);
}

/* Phase 6: inverse-trig builtins on Complex[MPFR, MPFR] and on MPFR
 * reals outside the real domain must produce MPFR-precision results
 * via the log-form complex identities. Pre-Phase-6 these silently
 * returned NaN (real path) or were stuck in the symbolic form. */
void test_arc_trig_mpfr_complex(void) {
    struct {
        const char* input;
        const char* expected_prefix;
        size_t prefix_len;
    } cases[] = {
        /* ArcSin[2] = Pi/2 - i acosh(2) ~ 1.57080 - 1.31696 i — real
         * input outside [-1, 1], complex result at MPFR precision. */
        {"ArcSin[N[2, 50]]",
         "1.57079632679489661923132169163975144209858469968755", 40},
        /* ArcCos[2] = i acosh(2) ~ 1.31696 i, the real part an exact 0 */
        {"ArcCos[N[2, 50]]",
         "1.3169578969248167086250463473079684440269819714675*I", 30},
        /* ArcTan[Complex[1, 1]] ~ 1.01722 + 0.40236 i */
        {"ArcTan[Complex[N[1, 50], N[1, 50]]]",
         "1.01722196789785136772278896155048292206356087698684", 40},
        /* ArcCot[Complex[1, 1]] = ArcTan[1/(1+i)] ~ 0.55357 - 0.40236 i */
        {"ArcCot[Complex[N[1, 50], N[1, 50]]]",
         "0.55357435889704525150853273008926852003502382270071", 40},
        /* ArcSec[Complex[1, 1]] ~ 1.11852 + 0.53064 i */
        {"ArcSec[Complex[N[1, 50], N[1, 50]]]",
         "1.11851787964370593716766329380877208138303741920675", 40},
        /* ArcCsc[Complex[1, 1]] ~ 0.45228 - 0.53064 i */
        {"ArcCsc[Complex[N[1, 50], N[1, 50]]]",
         "0.45227844715119068206365839783097936071554728048080", 40},
        {NULL, NULL, 0}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* r = evaluate(e);
        char* s = expr_to_string(r);
        ASSERT_MSG(strncmp(s, cases[i].expected_prefix, cases[i].prefix_len) == 0,
                   "%s: expected leading %s, got %s",
                   cases[i].input, cases[i].expected_prefix, s);
        free(s);
        expr_free(e);
        expr_free(r);
    }

    /* Precision[] round-trip on a representative complex case. */
    Expr* ep = parse_expression(
        "Precision[ArcSin[Complex[N[1, 50], N[1, 50]]]]");
    Expr* rp = evaluate(ep);
    char* sp = expr_to_string(rp);
    ASSERT_MSG(strncmp(sp, "50.", 3) == 0,
               "Precision[ArcSin[1+i]] (50 digits): expected 50.*, got %s",
               sp);
    free(sp); expr_free(ep); expr_free(rp);

    /* In-domain ArcSin still takes the real path. */
    Expr* eR = parse_expression("ArcSin[N[0.5, 50]]");
    Expr* rR = evaluate(eR);
    char* sR = expr_to_string(rR);
    ASSERT_MSG(strchr(sR, 'I') == NULL,
               "ArcSin[N[0.5, 50]]: expected pure real, got %s", sR);
    free(sR); expr_free(eR); expr_free(rR);
}

/* Capture stderr while `input` is parsed + evaluated. Returns the collected
 * stderr text (caller frees) and writes the printed result into
 * *out_result_str (also heap-allocated). Mirrors test_logexp.c. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_trig_stderr.log";
    int saved_stderr = test_stderr_redirect(path);
    if (saved_stderr < 0) {
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
    expr_free(p);
    expr_free(e);
    test_stderr_restore(saved_stderr);

    FILE* f = fopen(path, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    remove(path);
    return buf;
}

/* Elementary functions called with the wrong number of arguments must emit
 * a Mathematica-style argument-count diagnostic (`argx` for the fixed-arity
 * unary functions, `argt` for the 1-or-2-argument ArcTan) and leave the call
 * unevaluated. */
void test_elementary_argx(void) {
    struct { const char* input; const char* tag; const char* expect_phrase; size_t argc; } cases[] = {
        {"Sin[]",          "Sin::argx",    "1 argument is expected",      0},
        {"Sin[1, 2, 3]",   "Sin::argx",    "1 argument is expected",      3},
        {"Cos[1, 2]",      "Cos::argx",    "1 argument is expected",      2},
        {"Sec[]",          "Sec::argx",    "1 argument is expected",      0},
        {"ArcCsc[1, 2]",   "ArcCsc::argx", "1 argument is expected",      2},
        {"ArcTan[]",       "ArcTan::argt", "1 or 2 arguments are expected", 0},
        {"ArcTan[1, 2, 3]","ArcTan::argt", "1 or 2 arguments are expected", 3},
        {NULL, NULL, NULL, 0}
    };
    for (int i = 0; cases[i].input != NULL; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(cases[i].input, &result);
        ASSERT(result != NULL);
        ASSERT_MSG(strcmp(result, cases[i].input) == 0,
                   "%s: expected call unchanged, got %s", cases[i].input, result);
        ASSERT_MSG(err && strstr(err, cases[i].tag) != NULL,
                   "%s: expected %s diagnostic, got: %s",
                   cases[i].input, cases[i].tag, err ? err : "(null)");
        char needle[64];
        snprintf(needle, sizeof needle, "called with %zu argument%s",
                 cases[i].argc, cases[i].argc == 1 ? "" : "s");
        ASSERT_MSG(strstr(err, needle) != NULL,
                   "%s: missing arg-count phrase '%s' in: %s",
                   cases[i].input, needle, err);
        ASSERT_MSG(strstr(err, cases[i].expect_phrase) != NULL,
                   "%s: missing expectation phrase '%s' in: %s",
                   cases[i].input, cases[i].expect_phrase, err);
        free(result);
        free(err);
    }

    /* Valid arities must NOT emit a diagnostic. */
    struct { const char* input; } ok[] = {
        {"Sin[1]"}, {"ArcTan[3, 2]"}, {NULL}
    };
    for (int i = 0; ok[i].input != NULL; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(ok[i].input, &result);
        ASSERT_MSG(!(err && strstr(err, "::arg")),
                   "%s: unexpected arg diagnostic: %s", ok[i].input, err ? err : "");
        free(result);
        free(err);
    }
}

/* Degree (= Pi/180) must fold inside the trig heads exactly like the
 * corresponding rational multiple of Pi. See extract_pi_multiplier in
 * src/trig.c. */
void test_trig_degree() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[30 Degree]", "Rational[1, 2]"},               /* = Sin[Pi/6]  */
        {"Cos[45 Degree]", "Power[2, Rational[-1, 2]]"},     /* = Cos[Pi/4]  */
        {"Tan[60 Degree]", "Power[3, Rational[1, 2]]"},      /* = Tan[Pi/3]  */
        {"Cot[45 Degree]", "1"},
        {"Sec[60 Degree]", "2"},
        {"Csc[30 Degree]", "2"},
        {"Sin[90 Degree]", "1"},
        {"Cos[180 Degree]", "-1"},
        {"Sin[360 Degree]", "0"},
        {"Sin[-30 Degree]", "Rational[-1, 2]"},              /* odd fold      */
        {"Cos[120 Degree]", "Rational[-1, 2]"},
        {"Sin[210 Degree]", "Rational[-1, 2]"},              /* n >= d branch */
        {"Sin[15 Degree]",                                   /* d=12 case     */
         "Times[Rational[1, 4], Plus[Times[-1, Power[2, Rational[1, 2]]], Power[6, Rational[1, 2]]]]"},
        /* No closed form: 3/2 Degree = Pi/120 stays symbolic (not folded). */
        {"Sin[3/2 Degree]", "Sin[Times[Rational[3, 2], Degree]]"},
        /* Degree alone (= Pi/180) has no nice value: stays symbolic. */
        {"Sin[Degree]", "Sin[Degree]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing degree: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Degree %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

/* Every exact table value of the six forward heads at k Pi/d (d in the
 * tables' denominators, k over two periods, poles excluded) agrees with the
 * MPFR value of the same head at the numeric angle.  Csc[k Pi/5] for
 * k = +-1 mod 5 used to return Csc[2 Pi/5]'s magnitude. */
void test_trig_exact_table_numeric() {
    static const char* heads[] = { "Sin", "Cos", "Tan", "Cot", "Sec", "Csc" };
    static const int dens[] = { 3, 4, 5, 6, 10, 12 };
    char buf[512];
    for (size_t h = 0; h < sizeof heads / sizeof heads[0]; h++)
        for (size_t di = 0; di < sizeof dens / sizeof dens[0]; di++) {
            int d = dens[di];
            for (int k = -2 * d; k <= 2 * d; k++) {
                if (k % d == 0) continue;                       /* Sin/Csc/Tan/Cot poles */
                if (2 * k % d == 0 && (heads[h][0] == 'T' || heads[h][1] == 'e')) continue; /* Tan/Sec */
                snprintf(buf, sizeof buf, "Chop[N[%s[%d Pi/%d], 30] - %s[N[%d Pi/%d, 40]], 10^-25]",
                         heads[h], k, d, heads[h], k, d);
                Expr* e = parse_expression(buf);
                Expr* res = evaluate(e);
                char* s = expr_to_string(res);
                ASSERT_MSG(strcmp(s, "0") == 0, "%s[%d Pi/%d]: exact value off by %s", heads[h], k, d, s);
                free(s);
                expr_free(res);
                expr_free(e);
            }
        }
}

/* Every f[g[u]] of the six circular functions and their six inverses
 * evaluates to its principal-branch algebraic form, and agrees with
 * direct numerical evaluation at real and complex points (18 of the 36
 * used to stay unevaluated, e.g. Csc[ArcSin[u]], Sec[ArcTan[u]]). */
void test_trig_of_inverse_all_pairs() {
    assert_eval_eq("fs = {Sin, Cos, Tan, Cot, Sec, Csc}; "
                   "gs = {ArcSin, ArcCos, ArcTan, ArcCot, ArcSec, ArcCsc}; "
                   "FreeQ[Table[f[g[u]], {f, fs}, {g, gs}], "
                   "(Sin | Cos | Tan | Cot | Sec | Csc)[(ArcSin | ArcCos | ArcTan | ArcCot | ArcSec | ArcCsc)[_]]]",
                   "True", 0);
    assert_eval_eq("Max[Flatten[Table[Abs[N[(f[g[v]] /. v -> p) - f[N[g[p], 40]], 30]], "
                   "{f, fs}, {g, gs}, {p, {3/10, -7/10, 2, -3, 1/2 + I, -2 - 3 I/2}}]]] < 10^-25",
                   "True", 0);
    assert_eval_eq("{Csc[ArcSin[u]], Sec[ArcTan[u]], Cos[ArcSec[u]], Sin[ArcCsc[u]]}",
                   "{1/u, Sqrt[1 + u^2], 1/u, 1/u}", 0);
}

/* Arbitrary-precision ArcTan / ArcCot on the cuts agree with the exact value and
 * the machine kernel (v0.6702): ArcTan[0 - 2 I] read +Pi/2 - 0.549 I at MPFR
 * precision, its 1 - I z = -1 - 0 I carrying a negated zero into Log. */
static void test_arctan_mpfr_lower_cut(void) {
    static const char* cases[] = {
        "Abs[ArcTan[SetPrecision[-2 I, 40]] - N[ArcTan[-2 I], 40]] < 10^-30",
        "Abs[ArcTan[SetPrecision[-3 I, 40]] + ArcTan[SetPrecision[3 I, 40]]] < 10^-30",
        "Abs[ArcTan[SetPrecision[-2 I, 40]] - ArcTan[-2. I]] < 10^-14",
        "Abs[ArcCot[SetPrecision[I/2, 40]] - N[ArcCot[I/2], 40]] < 10^-30",
        "Abs[ArcCot[SetPrecision[-I/2, 40]] - N[ArcCot[-I/2], 40]] < 10^-30",
        "Abs[ArcTan[SetPrecision[2 I, 40]] - N[ArcTan[2 I], 40]] < 10^-30",
        NULL
    };
    for (int i = 0; cases[i]; i++) {
        Expr* e = parse_expression(cases[i]);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, "True") == 0, "%s: got %s", cases[i], s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

static void assert_all_true(const char* const* cases) {
    for (int i = 0; cases[i]; i++) {
        Expr* e = parse_expression(cases[i]);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, "True") == 0, "%s: got %s", cases[i], s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

/* The machine kernels and Compile[] on the imaginary-axis cuts of ArcTan,
 * ArcCot, ArcSinh, ArcCsch: C99 catan / casinh read the side from the sign
 * of the zero real part, so ArcTan[{-2. I}] was +Pi/2 - 0.549 I where the
 * scalar, the exact value and MPFR give -Pi/2 - 0.549 I. */
static void test_inverse_imag_cut_all_surfaces(void) {
    static const char* cases[] = {
        "And @@ Flatten[Table[With[{s = f[z]}, Abs[f[Developer`ToPackedArray[{z, 1. + I}]][[1]] - s] < 10^-15 && "
        "Abs[Normal[f[NDArray[{z, 1. + I}]]][[1]] - s] < 10^-15 && "
        "Abs[With[{g = f}, Compile[{{w, _Complex}}, g[w]]][z] - s] < 10^-15 && "
        "Abs[With[{g = f}, Compile[{{v, _Complex, 1}}, g[v]]][{z, 1. + I}][[1]] - s] < 10^-15 && "
        "Abs[s - N[f[SetPrecision[z, Infinity]]]] < 10^-14], "
        "{f, {ArcTan, ArcCot, ArcSinh, ArcCsch}}, "
        "{z, {2. I, -2. I, 0.5 I, -0.5 I, 3. I, -3. I, Complex[-0., 2.], Complex[-0., -2.], 10.^4 I, -10.^4 I}}]]",
        "Abs[ArcTan[Developer`ToPackedArray[{-2. I, 1.}]][[1]] - (-Pi/2 - ArcTanh[1/2] I)] < 10^-15",
        "Abs[ArcSinh[NDArray[{-2. I, 1.}]][[1]] - (-ArcCosh[2] - Pi/2 I)] < 10^-15",
        NULL
    };
    assert_all_true(cases);
}

/* MPFR complex inverses through cancellation-free forms: the log forms lost
 * every digit on half the plane (ArcSin[SetPrecision[-10^30, 40]] was
 * Overflow[]), near zero (ArcTan of 10^-30 I: six digits of 30), next to
 * the branch points (ArcCot at I + 10^-20: 23 of 30) and in Log near |z| = 1.
 * Each value is checked against the exact value at the same (rounded) input,
 * so only the algorithm's error is measured. */
static void test_mpfr_complex_inverse_accuracy(void) {
    static const char* cases[] = {
        "rel[f_, p_, n_] := With[{q = SetPrecision[p, n]}, Abs[f[q] - N[f[SetPrecision[q, Infinity]], n + 20]]/"
        "Abs[N[f[SetPrecision[q, Infinity]], n + 20]]]; True",
        "rel[ArcSin, -10^30, 40] < 10^-36",
        "rel[ArcSin, -10^30, 60] < 10^-56",
        "rel[ArcSin, 10^30 I, 60] < 10^-56",
        "rel[ArcCos, -10^30 + 2 10^30 I, 60] < 10^-56",
        "rel[ArcSinh, -10^30 I, 60] < 10^-56",
        "rel[ArcTan, 10^-30 I, 30] < 10^-26",
        "rel[ArcTan, 10^-30 (1 + I), 30] < 10^-26",
        "rel[ArcTanh, 10^-30 (2 - I), 30] < 10^-26",
        "rel[ArcTan, 10^-20 + I, 30] < 10^-26",
        "rel[ArcCot, 10^-20 + I, 30] < 10^-26",
        "rel[ArcCoth, 1 + 10^-20 I, 30] < 10^-26",
        "rel[ArcTanh, -1 + 10^-25 + 10^-25 I, 60] < 10^-56",
        "rel[ArcSin, 1 + 10^-25 + 10^-25 I, 30] < 10^-26",
        "rel[ArcSec, 1 + 10^-25, 60] < 10^-56",
        "rel[ArcCsc, -1 - 10^-25, 60] < 10^-56",
        "rel[ArcCoth, 1 + 10^-25, 60] < 10^-56",
        "rel[ArcSech, 1 - 10^-25, 60] < 10^-56",
        "rel[Log, 1 + 10^-20 I, 30] < 10^-26",
        "rel[ArcSech, -10^100 + I, 30] < 10^-26",   /* Smith's division put it on the cut */
        "rel[ArcSec, -1, 30] < 10^-26",
        NULL
    };
    assert_all_true(cases);
}

/* MPFR forward functions at a huge imaginary part: cosh^2 b overflowed and
 * inf/inf (or 0 inf) came back Indeterminate. */
static void test_mpfr_complex_forward_huge(void) {
    static const char* cases[] = {
        "Abs[Tan[SetPrecision[10^30 I, 30]] - I] < 10^-29",
        "Abs[Cot[SetPrecision[-10^9 I, 30]] - I] < 10^-29",
        "Abs[Tanh[SetPrecision[10^30 + 10^30 I, 30]] - 1] < 10^-29",
        "Abs[Coth[SetPrecision[-10^100 + I, 30]] + 1] < 10^-29",
        "Abs[Sec[SetPrecision[10^30 I, 30]]] < 10^-29",
        "Abs[Csch[SetPrecision[10^9 + I, 30]]] < 10^-29",
        "Sin[SetPrecision[10^30 I, 30]] =!= Indeterminate",
        "Abs[Tan[SetPrecision[1/3 + 10^4 I, 30]] - N[Tan[1/3 + 10^4 I], 30]] < 10^-28",
        /* E^(complex) in N: exp(y log e) multiplied the rounding of log e by |y| */
        "With[{d = SetPrecision[1.*^100, Infinity]}, Abs[N[Exp[I d], 30] - N[Cos[d] + I Sin[d], 30]] < 10^-28]",
        NULL
    };
    assert_all_true(cases);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_trig_exact_table_numeric);

    TEST(test_trig_forward);
    TEST(test_trig_inverse);
    TEST(test_trig_forward_of_inverse);
    TEST(test_trig_of_inverse_all_pairs);
    TEST(test_trig_mpfr_complex);
    TEST(test_arc_trig_mpfr_complex);
    TEST(test_arctan_mpfr_lower_cut);
    TEST(test_inverse_imag_cut_all_surfaces);
    TEST(test_mpfr_complex_inverse_accuracy);
    TEST(test_mpfr_complex_forward_huge);
    TEST(test_elementary_argx);
    TEST(test_trig_degree);

    return 0;
}
