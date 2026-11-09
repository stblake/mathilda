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

    /* Contour-shaped complex bodies: integer powers by binary powering, shared
     * subtrees lowered once (value numbering), E^(k I t) as (E^(I t))^k, a
     * bignum Rational coefficient, and the reciprocal trig functions -- each
     * against the interpreter on the same body at the same precision. */
    run("cH = Compile[{{t, _Real}}, (1 + E^(I t))^200 (1 + E^(I t) + E^(2 I t))^70 E^(-3 I t) "
        "+ E^(-I t) + (2 + 3^40/7^30 E^(5 I t))^-9, WorkingPrecision -> 60]");
    run("t0 = N[3/7, 60]");
    check_true("Abs[cH[t0]/((1 + E^(I t0))^200 (1 + E^(I t0) + E^(2 I t0))^70 E^(-3 I t0) "
               "+ E^(-I t0) + (2 + 3^40/7^30 E^(5 I t0))^-9) - 1] < 10^-50",
               "mpfr complex: powers, shared harmonics, bignum rational @ 60");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{t, _Complex}}, Evaluate[(1 + E^(I t))^200 "
               "E^(-3 I t) + 3^40/7^30], WorkingPrecision -> 40]) === True",
               "diag: harmonics and a bignum Rational lower");
    run("cT = Compile[{{z, _Complex}}, Csc[z]^3 + Tan[z] + Sec[z] - Cot[z], WorkingPrecision -> 40]");
    run("z2 = N[7/10 + 2 I/5, 40]");
    check_true("Abs[cT[z2]/(Csc[z2]^3 + Tan[z2] + Sec[z2] - Cot[z2]) - 1] < 10^-35",
               "mpfr complex: Csc/Tan/Sec/Cot @ 40");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Csc[z]^3 + Tan[z], "
               "WorkingPrecision -> 40]) === True", "diag: complex reciprocal trig lowers");
    /* real rational powers by q-th roots; a negative base is interpreted (its
     * principal power is complex), not answered with the real root */
    run("cR = Compile[{{x, _Real}}, x^(1/3) + x^(-2/5) + x^(7/4), WorkingPrecision -> 50]");
    run("x4 = N[17/7, 50]");
    check_true("Abs[cR[x4]/(x4^(1/3) + x4^(-2/5) + x4^(7/4)) - 1] < 10^-45", "mpfr real: x^(p/q) @ 50");
    run("cQ = Compile[{{x, _Real}}, x^(1/3), WorkingPrecision -> 40]");
    check_true("Abs[cQ[N[-8, 40]] - (1 + Sqrt[3] I)] < 10^-35", "mpfr real: (-8)^(1/3) stays principal");
    /* complex Sqrt (algebraic form): every quadrant, the principal branch */
    run("cS = Compile[{{z, _Complex}}, Sqrt[z], WorkingPrecision -> 40]");
    check_true("Abs[cS[N[3 - 4 I, 40]] - (2 - I)] < 10^-35 && Abs[cS[N[-3 + 4 I, 40]] - (1 + 2 I)] < 10^-35"
               " && Abs[cS[N[-3 - 4 I, 40]] - (1 - 2 I)] < 10^-35 && Abs[cS[N[5/4 + 3 I, 40]] - (3/2 + I)] < 10^-35",
               "mpfr complex: Sqrt in all quadrants @ 40");
    run("rT = Compile[{{x, _Real}}, Csc[x] + Sec[x] + Cot[x], WorkingPrecision -> 40]");
    run("x3 = N[6/5, 40]");
    check_true("Abs[rT[x3]/(Csc[x3] + Sec[x3] + Cot[x3]) - 1] < 10^-35", "mpfr real: Csc/Sec/Cot @ 40");

    /* the reciprocal / inverse hyperbolic functions and With (lowered as the
     * substituted body, the shared binding evaluated once) */
    run("hY = Compile[{{x, _Real}}, Sech[x] + Csch[x]^2 + Coth[x] + ArcSinh[x] + ArcTanh[x/3] "
        "+ ArcCosh[x + 1], WorkingPrecision -> 45]");
    run("xh = N[6/5, 45]");
    check_true("Abs[hY[xh]/(Sech[xh] + Csch[xh]^2 + Coth[xh] + ArcSinh[xh] + ArcTanh[xh/3] "
               "+ ArcCosh[xh + 1]) - 1] < 10^-40", "mpfr real: Sech/Csch/Coth/ArcSinh/ArcTanh/ArcCosh @ 45");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Sech[x] + ArcTanh[x], "
               "WorkingPrecision -> 40]) === True", "diag: real hyperbolic reciprocals lower");
    run("wW = Compile[{{j, _Real}, {k, _Real}}, With[{r = Pi Sqrt[j^2 + k^2], q = j}, "
        "(Exp[-r/4] - Exp[-7 r/4])/(1 - Exp[-2 r]) + q], WorkingPrecision -> 50]");
    run("j0 = N[3, 50]; k0 = N[5, 50]; r0 = Pi Sqrt[j0^2 + k0^2]");
    check_true("Abs[wW[j0, k0]/((Exp[-r0/4] - Exp[-7 r0/4])/(1 - Exp[-2 r0]) + j0) - 1] < 10^-45",
               "mpfr real: With binding @ 50");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{j, _Real}}, With[{r = Sqrt[j]}, r + r^2], "
               "WorkingPrecision -> 40]) === True", "diag: With lowers on the managed path");
    /* a With local shadowing an argument is the local */
    run("wS = Compile[{{x, _Real}}, With[{x = 2}, x + 1], WorkingPrecision -> 30]");
    check_true("wS[N[7, 30]] == 3", "mpfr real: With local shadows the argument");
    /* ArcCosh below 1 is complex: the managed real kernel declines, the
     * interpreter answers */
    run("aC = Compile[{{x, _Real}}, ArcCosh[x], WorkingPrecision -> 30]");
    check_true("Abs[aC[N[1/2, 30]] - ArcCosh[N[1/2, 30]]] < 10^-25", "mpfr real: ArcCosh[1/2] complex");

    /* A special function (Gamma: compiled through MG_SPECIAL since v0.6482), still
     * the interpreter's value at precision. */
    run("gF = Compile[{{x, _Real}}, Gamma[x], WorkingPrecision -> 30]");
    run("xg = N[7/2, 30]");
    check_true("gF[xg] == N[Gamma[7/2], 30]", "special function: Gamma at precision");

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
    /* Gamma lowers since the special-function registry (MG_SPECIAL); an
     * integer-valued head such as Floor is still outside the managed subset */
    check_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Gamma[x], WorkingPrecision -> 40]) === True",
               "diag: special function (Gamma) lowers");
    check_true("(\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Floor[x], WorkingPrecision -> 40]) === False",
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
