/* Compile[] arbitrary precision (WorkingPrecision -> n MPFR reals/complex,
 * "BigIntegers" -> True GMP integers).  The compiled managed path must reproduce
 * the interpreter's FIXED-PRECISION numeric evaluation of the same body on the
 * same argument to full precision (not exact-then-round), stay exact for bignum,
 * and leave the machine path (WorkingPrecision -> MachinePrecision / no option)
 * unchanged.  See docs/design/compile-arbitrary-precision.md.
 *
 * Soft asserts: prints FAIL and keeps going. Run: ./compile_arbprec_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "numeric.h"                 /* numeric_digits_to_bits */
#include "compile/autocompile.h"     /* autocompile_new_prec / autocompiled_eval_mpfr */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int failures = 0;

static void run(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) { printf("FAIL: could not parse setup: %s\n", src); failures++; return; }
    expr_free(evaluate(p));
    expr_free(p);
}

/* Evaluate `src` to a freshly owned Expr (caller frees). */
static Expr* eval_str(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) return NULL;
    Expr* v = evaluate(p);
    expr_free(p);
    return v;
}

/* Evaluate `src`; require the result to be the symbol True. */
static void check_true(const char* src, const char* label) {
    Expr* parsed = parse_expression(src);
    if (!parsed) { printf("FAIL: %s (parse)\n", label); failures++; return; }
    Expr* v = evaluate(parsed);
    expr_free(parsed);
    bool ok = v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "True") == 0;
    if (!ok) {
        char* s = v ? expr_to_string(v) : NULL;
        printf("FAIL: %s -> %s\n", label, s ? s : "(null)");
        free(s);
        failures++;
    } else {
        printf("ok:   %s\n", label);
    }
    expr_free(v);
}

int main(void) {
    symtab_init();
    core_init();

    /* Faithful equivalence: compiled == interpreter of the SAME body, SAME
     * precision, SAME numeric argument. */
    run("fR = Compile[{{x, _Real}}, x*x + x, WorkingPrecision -> 50]");
    run("xr = N[1/3, 50]");
    check_true("fR[xr] == (xr*xr + xr)", "mpfr real: x*x+x @ 50");

    run("fK = Compile[{{x, _Real}}, Sin[x]*Cos[x] - Sqrt[x] + x^3, WorkingPrecision -> 45]");
    run("yk = N[7/5, 45]");
    check_true("fK[yk] == (Sin[yk]*Cos[yk] - Sqrt[yk] + yk^3)", "mpfr real: mixed kernels @ 45");

    run("fD = Compile[{{x, _Real}}, 1/x + x/3, WorkingPrecision -> 40]");
    run("xd = N[2, 40]");
    check_true("fD[xd] == (1/xd + xd/3)", "mpfr real: divide @ 40");

    run("fP = Compile[{{x, _Real}}, Pi*x + E, WorkingPrecision -> 40]");
    run("xp = N[2, 40]");
    check_true("fP[xp] == (Pi*xp + E)", "mpfr real: named constants @ 40");

    /* Precision of the result is the requested working precision. */
    check_true("Precision[fP[xp]] > 39", "mpfr real: result precision >= requested");

    /* Bignum: exact past int64. */
    run("bI = Compile[{{n, _Integer}}, n^20, \"BigIntegers\" -> True]");
    check_true("bI[99] == 99^20", "bignum: n^20 exact");
    run("bS = Compile[{{n, _Integer}}, n*n*n - 2*n, \"BigIntegers\" -> True]");
    check_true("bS[1000000000000] == (1000000000000^3 - 2*1000000000000)", "bignum: polynomial exact");

    /* MPFR complex. */
    run("cE = Compile[{{z, _Complex}}, Exp[z] + z^2, WorkingPrecision -> 45]");
    run("z0 = N[1/2 + I/3, 45]");
    check_true("cE[z0] == (Exp[z0] + z0^2)", "mpfr complex: exp+z^2 @ 45");
    run("cM = Compile[{{z, _Complex}}, z*z + z, WorkingPrecision -> 40]");
    run("z1 = N[1 + I, 40]");
    check_true("cM[z1] == (z1*z1 + z1)", "mpfr complex: z^2+z @ 40");

    /* Interpreter fallback on an unsupported head, still correct at precision. */
    run("gF = Compile[{{x, _Real}}, Gamma[x], WorkingPrecision -> 30]");
    run("xg = N[7/2, 30]");
    check_true("gF[xg] == N[Gamma[7/2], 30]", "fallback: Gamma bails to interpreter");

    /* Machine path unchanged: no option, and WorkingPrecision -> MachinePrecision. */
    run("mA = Compile[{{x, _Real}}, x*x + x]");
    check_true("mA[0.5] == 0.75", "machine: value unchanged");
    check_true("Head[mA[0.5]] === Real", "machine: result Head is Real");
    run("mB = Compile[{{x, _Real}}, x^2, WorkingPrecision -> MachinePrecision]");
    check_true("Head[mB[3.0]] === Real", "machine: WorkingPrecision->MachinePrecision stays machine");

    /* CompileDiagnostics: verify the managed subset actually LOWERS (Compiled ->
     * True) rather than silently falling back — a value check cannot see this,
     * since a bailed body still answers correctly through the interpreter.  This
     * is the regression guard against the managed subset rotting. */
    check_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Sin[x]*Cos[x] + x^3, WorkingPrecision -> 40]) === True",
               "diag: managed real subset lowers");
    check_true("(\"ResultType\" /. CompileDiagnostics[{{x, _Real}}, x + 1, WorkingPrecision -> 40]) === \"MPFRReal\"",
               "diag: managed result type is MPFRReal");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Exp[z] + z^2, WorkingPrecision -> 40]) === True",
               "diag: managed complex subset lowers");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{n, _Integer}}, n^4 - 2 n, \"BigIntegers\" -> True]) === True",
               "diag: bignum subset lowers");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Gamma[x], WorkingPrecision -> 40]) === False",
               "diag: unsupported head reports the bail");

    /* Phase 3 substrate: precision-aware auto-compilation (the entry a sampler
     * running at high WorkingPrecision uses).  Compile Sin[x]+x^2 in MPFR, feed
     * one MPFR sample point, and require the result to equal the interpreter's. */
    {
        long bits = numeric_digits_to_bits(40.0);
        Expr* body = parse_expression("Sin[x] + x^2");
        Expr* xsym = parse_expression("x");
        const Expr* vars[1] = { xsym };
        AutoCompiled* ac = autocompile_new_prec(body, vars, 1, bits);
        if (!ac) { printf("FAIL: autocompile: new_prec returned NULL\n"); failures++; }
        else {
            Expr* xval = eval_str("N[7/5, 40]");           /* an EXPR_MPFR sample point */
            const Expr* xs[1] = { xval };
            Expr* got = autocompiled_eval_mpfr(ac, xs);    /* owned result */
            Expr* ref = eval_str("Sin[N[7/5, 40]] + N[7/5, 40]^2");
            bool ok = got && ref;
            if (ok) {
                Expr* eqargs[2] = { got, ref };            /* consumes got, ref */
                Expr* eq = expr_new_function(expr_new_symbol("Equal"), eqargs, 2);
                Expr* r = evaluate(eq);
                ok = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
                expr_free(eq);
                expr_free(r);
            } else { expr_free(got); expr_free(ref); }
            if (ok) printf("ok:   autocompile: MPFR sampler substrate matches interpreter\n");
            else { printf("FAIL: autocompile: MPFR sampler substrate mismatch\n"); failures++; }
            expr_free(xval);
            autocompiled_free(ac);
        }
        /* prec_bits <= 0 must decline (the machine path is autocompile_new). */
        AutoCompiled* none = autocompile_new_prec(body, vars, 1, 0);
        if (none) { printf("FAIL: autocompile: new_prec(0) should be NULL\n"); failures++; autocompiled_free(none); }
        else printf("ok:   autocompile: new_prec(prec<=0) declines\n");
        expr_free(body);
        expr_free(xsym);
    }

    if (failures == 0) printf("\nAll Compile[] arbitrary-precision tests passed.\n");
    else               printf("\n%d Compile[] arbitrary-precision test(s) FAILED.\n", failures);
    symtab_clear();
    return 0;   /* soft-assert harness: grep FAIL */
}
