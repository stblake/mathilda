/* test_options.c — Options, SetOptions, OptionValue, and OptionsPattern. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

static int failures = 0;

/* Evaluate `input`; assert its printed form equals `expected`. Uses exit-on-
 * fail rather than libc assert so it still fires under -DNDEBUG. */
static void chk(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL (parse): %s\n", input); failures++; return; }
    Expr* val = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(val);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(val);
}

/* Evaluate and discard (for definitions / side effects). */
static void run(const char* input) {
    Expr* parsed = parse_expression(input);
    if (!parsed) { fprintf(stderr, "FAIL (parse): %s\n", input); failures++; return; }
    Expr* val = evaluate(parsed);
    expr_free(parsed);
    expr_free(val);
}

int main(void) {
    symtab_init();
    core_init();

    /* ---- Options[sym]: empty, then assigned ---- */
    chk("Options[fa]", "{}");
    chk("Options[fa] = {a -> 1, b -> 2}", "{a -> 1, b -> 2}");
    chk("Options[fa]", "{a -> 1, b -> 2}");

    /* ---- registered builtin defaults (comprehensive sweep) ---- */
    chk("Options[LinearSolve]", "{Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}");
    chk("Options[NSeries]", "{Radius -> 1.0, WorkingPrecision -> MachinePrecision}");
    /* Symbolic builtins now carry their honored options (was empty {}). */
    chk("Options[Integrate]", "{Method -> Automatic}");
    chk("Options[Limit]", "{Direction -> Automatic, Assumptions -> Automatic}");
    chk("Options[D]", "{NonConstants -> {}}");
    chk("Options[GroebnerBasis]",
        "{MonomialOrder -> Lexicographic, CoefficientDomain -> Rationals, "
        "Method -> Automatic, Sort -> False, Modulus -> 0}");
    chk("Options[Eigenvalues]", "{Cubics -> True, Quartics -> True}");
    chk("Options[FactorInteger]", "{GaussianIntegers -> False}");
    chk("Options[SingularValueDecomposition]",
        "{Tolerance -> Automatic, TargetStructure -> \"Dense\"}");
    /* Heads default differs across the structural family. */
    chk("Options[Cases]", "{Heads -> False}");
    chk("Options[Position]", "{Heads -> True}");
    /* OptionValue resolves a registered default with no explicit options. */
    chk("OptionValue[Cases, Heads]", "False");
    chk("OptionValue[Eigenvalues, Cubics]", "True");
    /* SetOptions edits a registered default in place. */
    chk("SetOptions[Integrate, Method -> \"RischNorman\"]", "{Method -> \"RischNorman\"}");
    chk("Options[Integrate]", "{Method -> \"RischNorman\"}");

    /* ---- Options[obj, name] / Options[obj, {names}] ---- */
    chk("Options[fa, a]", "{a -> 1}");
    chk("Options[fa, {a, b}]", "{a -> 1, b -> 2}");
    chk("Options[fa, zzz]", "{}");

    /* ---- Options[expr]: explicit options of a compound expression ---- */
    chk("Options[Graphics[Circle[], Axes -> True, PlotRange -> All]]",
        "{Axes -> True, PlotRange -> All}");

    /* ---- SetOptions: update existing, return new list, leave unknown ---- */
    chk("SetOptions[fa, a -> 9]", "{a -> 9, b -> 2}");
    chk("Options[fa]", "{a -> 9, b -> 2}");
    /* Unknown option: emits SetOptions::optnf and stays unevaluated. */
    chk("SetOptions[fa, c -> 2]", "SetOptions[fa, c -> 2]");
    chk("Options[fa]", "{a -> 9, b -> 2}");       /* unchanged after the error */

    /* ---- AppendTo / PrependTo grow the option list, then SetOptions works ---- */
    chk("AppendTo[Options[fa], c -> 3]", "{a -> 9, b -> 2, c -> 3}");
    chk("SetOptions[fa, c -> 4]", "{a -> 9, b -> 2, c -> 4}");
    chk("PrependTo[Options[fa], z -> 0]", "{z -> 0, a -> 9, b -> 2, c -> 4}");

    /* ---- OptionValue explicit forms ---- */
    chk("OptionValue[fb, {a -> 5}, a]", "5");
    chk("OptionValue[fb, {a -> 5}, a, Hold]", "Hold[5]");
    chk("OptionValue[fb, {a -> 5}, \"a\"]", "5");  /* string/symbol name interchangeable */
    run("Options[fb] = {a -> 1, b -> 2}");
    chk("OptionValue[fb, a]", "1");                /* falls back to Options[fb] defaults */
    chk("OptionValue[fb, {a -> 7}, b]", "2");      /* explicit miss -> default */

    /* ---- End-to-end OptionsPattern + OptionValue ---- */
    run("Options[g] = {a -> 1, b -> 2}");
    run("g[OptionsPattern[]] := {OptionValue[a], OptionValue[b]}");
    chk("g[]", "{1, 2}");
    chk("g[a -> 17]", "{17, 2}");
    chk("g[b -> 18]", "{1, 18}");
    chk("g[a -> 17, b -> 18]", "{17, 18}");
    chk("{g[], g[a -> 17], g[b -> 18], g[a -> 17, b -> 18]}",
        "{{1, 2}, {17, 2}, {1, 18}, {17, 18}}");
    /* A non-option argument does not match the OptionsPattern rule. */
    chk("g[x]", "g[x]");
    /* Options passed as a list are flattened by OptionsPattern. */
    chk("g[{a -> 3, b -> 4}]", "{3, 4}");

    /* ---- OptionsPattern after a fixed argument ---- */
    run("Options[h] = {opt -> 0}");
    run("h[x_, OptionsPattern[]] := {x, OptionValue[opt]}");
    chk("h[5]", "{5, 0}");
    chk("h[5, opt -> 9]", "{5, 9}");

    if (failures) {
        fprintf(stderr, "\n%d test(s) FAILED\n", failures);
        return 1;
    }
    printf("\nAll options tests passed.\n");
    return 0;
}
