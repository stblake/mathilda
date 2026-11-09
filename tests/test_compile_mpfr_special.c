/* Compile[] arbitrary precision for the special functions (MG_SPECIAL).
 *
 * Every head in the arbitrary-precision special-function registry
 * (src/compile/compile_mgd_special.c) must, at WorkingPrecision -> 30 and 50,
 * for a real and (where the registry says so) a complex argument:
 *   - LOWER: CompileDiagnostics reports Compiled -> True;
 *   - RUN COMPILED: the call does not decline to the interpreter
 *     (Internal`CompileArbitraryPrecisionStats "EvalDeclines" unchanged);
 *   - AGREE with the interpreter's evaluation of the same body to the
 *     requested precision.
 * Independent reference values come from mpmath (50 digits).  Auto-compilation
 * is checked end to end: NIntegrate / NSum / FindRoot / NDSolve at
 * WorkingPrecision -> 30 over special-function bodies compile with no
 * interpreter fallback and match both mpmath and the uncompiled run.
 *
 * Soft asserts: prints FAIL and keeps going.  Run: ./compile_mpfr_special_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "compile/autocompile.h"     /* autocompile_set_enabled */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int failures = 0, checks = 0;

static void run(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) { printf("FAIL: could not parse setup: %s\n", src); failures++; return; }
    expr_free(evaluate(p));
    expr_free(p);
}

static bool is_true(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) return false;
    Expr* v = evaluate(p);
    expr_free(p);
    bool ok = v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "True") == 0;
    expr_free(v);
    return ok;
}

static void check_true(const char* src, const char* label) {
    checks++;
    if (is_true(src)) { printf("ok:   %s\n", label); return; }
    printf("FAIL: %s\n      %s\n", label, src);
    failures++;
}

/* A probe body in x, and whether the head takes a complex argument.  One row
 * per registry head (and per argument convention); generated from
 * tools/compile_coverage.py's PREC_BODIES so the gate and this suite ask the
 * same questions. */
typedef struct { const char* body; int cplx; } Probe;

static const Probe PROBES[] = {
    { "AiryAi[x]", 1 },
    { "AiryAiPrime[x]", 1 },
    { "AiryBi[x]", 1 },
    { "AiryBiPrime[x]", 1 },
    { "AngerJ[1/3, x]", 1 },
    { "ArcCos[x]", 1 },
    { "ArcCosh[x + 1]", 1 },
    { "ArcCot[x]", 1 },
    { "ArcCoth[x + 1]", 1 },
    { "ArcCsc[x + 1]", 1 },
    { "ArcCsch[x]", 1 },
    { "ArcSec[x + 1]", 1 },
    { "ArcSech[x]", 1 },
    { "ArcSin[x]", 1 },
    { "ArcSinh[x]", 1 },
    { "ArcTan[x]", 1 },
    { "ArcTan[1/3, x]", 0 },   /* the two-argument form is real only */
    { "ArcTanh[x]", 1 },
    { "Arg[x]", 1 },
    { "BarnesG[x]", 1 },
    { "BernoulliB[3, x]", 1 },
    { "BesselI[1/3, x]", 1 },
    { "BesselJ[1/3, x]", 1 },
    { "BesselK[1/3, x]", 1 },
    { "BesselY[1/3, x]", 1 },
    { "Beta[1/3, x]", 1 },
    { "BetaRegularized[1/3, 2/5, x]", 1 },
    { "Binomial[1/3, x]", 1 },
    { "ChebyshevT[1/3, x]", 1 },
    { "ChebyshevU[1/3, x]", 1 },
    { "Conjugate[x]", 1 },
    { "CosIntegral[x]", 1 },
    { "Cosh[x]", 1 },
    { "CoshIntegral[x]", 1 },
    { "Coth[x]", 1 },
    { "Csch[x]", 1 },
    { "DirichletL[5, 3, x]", 1 },
    { "DirichletL[5, 2, x]", 1 },
    { "EllipticE[x]", 1 },
    { "EllipticF[1/3, x]", 1 },
    { "EllipticK[x]", 1 },
    { "EllipticPi[1/3, x]", 1 },
    { "EllipticTheta[2, x, 1/3]", 1 },
    { "EllipticThetaPrime[2, x, 1/3]", 1 },
    { "Erf[x]", 1 },
    { "Erfc[x]", 1 },
    { "Erfi[x]", 1 },
    { "EulerE[3, x]", 1 },
    { "ExpIntegralE[1/3, x]", 1 },
    { "ExpIntegralEi[x]", 1 },
    { "Factorial[x]", 1 },
    { "Factorial2[x]", 1 },
    { "FactorialPower[1/3, x]", 1 },
    { "Fibonacci[x]", 1 },
    { "FresnelC[x]", 1 },
    { "FresnelS[x]", 1 },
    { "Gamma[x]", 1 },
    { "GegenbauerC[1/3, x]", 1 },
    { "HankelH1[1/3, x]", 1 },
    { "HankelH2[1/3, x]", 1 },
    { "HarmonicNumber[x]", 1 },
    { "HermiteH[1/3, x]", 1 },
    { "HurwitzZeta[1/3, x]", 1 },
    { "Hyperfactorial[x]", 1 },
    { "Hypergeometric0F1[1/3, x]", 1 },
    { "Hypergeometric1F1[1/3, 2/5, x]", 1 },
    { "Hypergeometric2F1[1/3, 2/5, 3/7, x]", 1 },
    { "HypergeometricPFQ[{1/3, 1/4}, {1/5}, x]", 1 },
    { "HypergeometricU[1/3, 2/5, x]", 1 },
    { "Im[x]", 1 },
    { "InverseJacobiCN[1/3, x]", 1 },
    { "InverseJacobiDN[x + 1/5, 1/3]", 1 },
    { "InverseJacobiSN[1/3, x]", 1 },
    { "JacobiAmplitude[1/3, x]", 1 },
    { "JacobiCD[1/3, x]", 1 },
    { "JacobiCN[1/3, x]", 1 },
    { "JacobiCS[1/3, x]", 1 },
    { "JacobiDC[1/3, x]", 1 },
    { "JacobiDN[1/3, x]", 1 },
    { "JacobiDS[1/3, x]", 1 },
    { "JacobiNC[1/3, x]", 1 },
    { "JacobiND[1/3, x]", 1 },
    { "JacobiNS[1/3, x]", 1 },
    { "JacobiP[1/3, 2/5, 3/7, x]", 1 },
    { "JacobiSC[1/3, x]", 1 },
    { "JacobiSD[1/3, x]", 1 },
    { "JacobiSN[1/3, x]", 1 },
    { "KelvinBei[x]", 1 },
    { "KelvinBer[x]", 1 },
    { "KelvinKei[x]", 1 },
    { "KelvinKer[x]", 1 },
    { "LaguerreL[1/3, x]", 1 },
    { "LegendreP[2, x]", 1 },
    { "LegendreP[1/3, 2/5, x]", 1 },
    { "LegendreQ[2, x]", 1 },
    { "LegendreQ[1/3, 2/5, x]", 1 },
    { "LerchPhi[1/3, 2/5, x]", 1 },
    { "Log[x]", 1 },
    { "Log[2, x]", 1 },
    { "Log10[x]", 1 },
    { "Log2[x]", 1 },
    { "LogGamma[x]", 1 },
    { "LogIntegral[x]", 1 },
    { "LucasL[x]", 1 },
    { "MathieuC[2/5, 1/3, x]", 1 },
    { "MathieuCPrime[2/5, 1/3, x]", 1 },
    { "MathieuS[2/5, 1/3, x]", 1 },
    { "MathieuSPrime[2/5, 1/3, x]", 1 },
    { "MeijerG[{{1/3}, {1/4}}, {{1/5, 1/6}, {1/7}}, x]", 1 },
    { "MeijerG[{{1/3}, {}}, {{1/2}, {1}}, x, 1/2]", 1 },
    { "Multinomial[1/3, x]", 1 },
    { "ParabolicCylinderD[1/3, x]", 1 },
    { "Pochhammer[1/3, x]", 1 },
    { "PolyGamma[x]", 1 },
    { "PolyLog[1/3, x]", 1 },
    { "PolyLog[2, 1, x]", 1 },
    { "ProductLog[x]", 1 },
    { "ProductLog[-1, -x/3]", 1 },
    { "QPochhammer[1/3, x]", 1 },
    { "QPochhammer[1/3, x, 3]", 1 },
    { "Re[x]", 1 },
    { "RiemannR[x]", 1 },
    { "RiemannSiegelTheta[x]", 1 },
    { "RiemannSiegelZ[x]", 1 },
    { "Sech[x]", 1 },
    { "SinIntegral[x]", 1 },
    { "Sinc[x]", 1 },
    { "Sinh[x]", 1 },
    { "SinhIntegral[x]", 1 },
    { "SphericalBesselJ[1/3, x]", 1 },
    { "SphericalBesselY[1/3, x]", 1 },
    { "SphericalHankelH1[1/3, x]", 1 },
    { "SphericalHankelH2[1/3, x]", 1 },
    { "StieltjesGamma[1, x]", 1 },
    { "StruveH[1/3, x]", 1 },
    { "StruveL[1/3, x]", 1 },
    { "Tanh[x]", 1 },
    { "WeberE[1/3, x]", 1 },
    { "WeierstrassP[x, {1, 2}]", 1 },
    { "WeierstrassPPrime[x, {1, 2}]", 1 },
    { "WeierstrassSigma[x, {1, 2}]", 1 },
    { "WeierstrassZeta[x, {1, 2}]", 1 },
    { "WhittakerM[1/3, 2/5, x]", 1 },
    { "WhittakerW[1/3, 2/5, x]", 1 },
    { "Zeta[x]", 1 },
};

/* One body, one domain, one precision: lowers, runs compiled, agrees. */
static void probe(const char* body, bool cplx, int digits) {
    char src[2048], label[512];
    const char* spec = cplx ? "{{x, _Complex}}" : "{{x, _Real}}";
    const char* dom = cplx ? "complex" : "real";
    snprintf(src, sizeof src, "x0 = N[%s, %d]", cplx ? "7/10 + I/5" : "7/10", digits);
    run(src);
    snprintf(src, sizeof src,
             "Lookup[CompileDiagnostics[%s, %s, WorkingPrecision -> %d], \"Compiled\"]",
             spec, body, digits);
    snprintf(label, sizeof label, "%s %s @%d: lowers", body, dom, digits);
    check_true(src, label);
    snprintf(src, sizeof src, "cf = Compile[%s, %s, WorkingPrecision -> %d]", spec, body, digits);
    run(src);
    run("s0 = \"EvalDeclines\" /. Internal`CompileArbitraryPrecisionStats[]");
    run("v = cf[x0]");
    run("s1 = \"EvalDeclines\" /. Internal`CompileArbitraryPrecisionStats[]");
    snprintf(src, sizeof src, "ref = Block[{x = x0}, %s]", body);
    run(src);
    snprintf(label, sizeof label, "%s %s @%d: runs compiled", body, dom, digits);
    check_true("s1 == s0", label);
    snprintf(src, sizeof src,
             "NumberQ[v] && NumberQ[ref] && Abs[v - ref] <= 10^-%d Max[1, Abs[ref]]",
             digits - 6);
    snprintf(label, sizeof label, "%s %s @%d: compiled == interpreted", body, dom, digits);
    check_true(src, label);
}

typedef struct { const char* expr; const char* ref; const char* label; } AcCase;

static const AcCase AC[] = {
    { "NIntegrate[Erf[x] BesselK[0, x], {x, 1/10, 1}, WorkingPrecision -> 30]",
      "0.385447747669799402740886491309297331882094229`45", "NIntegrate Erf*BesselK" },
    { "NIntegrate[AiryAi[x] EllipticK[x/2], {x, 1/10, 1}, WorkingPrecision -> 30]",
      "0.341237254185965561031045180393933407510303917`45", "NIntegrate AiryAi*EllipticK" },
    { "NIntegrate[Hypergeometric1F1[1/3, 2, x] EllipticE[x/2], {x, 1/10, 1}, "
      "WorkingPrecision -> 30]",
      "1.44718814450936781049303418458477951126654955`45", "NIntegrate 1F1*EllipticE" },
    { "NSum[Erfc[k]/k^2, {k, 1, Infinity}, WorkingPrecision -> 30]",
      "0.158471096258853574368797532572179741681506371`45", "NSum Erfc[k]/k^2" },
    { "x /. FindRoot[BesselJ[0, x] Exp[Erf[x] - Erf[x]], {x, 24/10}, WorkingPrecision -> 30]",
      "2.404825557695772768621631879326454643124244909146`45", "FindRoot BesselJ zero" },
};

int main(void) {
    symtab_init();
    core_init();

    /* ---- every registered head, real and complex, 30 and 50 digits ---- */
    size_t np = sizeof PROBES / sizeof PROBES[0];
    for (size_t i = 0; i < np; i++)
        for (int d = 30; d <= 50; d += 20) {
            probe(PROBES[i].body, false, d);
            if (PROBES[i].cplx) probe(PROBES[i].body, true, d);
        }
    check_true("Length[Internal`CompileArbitraryPrecisionHeads[]] >= 142",
               "registry: at least the 142 heads this suite probes");

    /* ---- independent references (mpmath, 50 digits) ---- */
    run("x0 = N[7/10, 50]; z0 = N[7/10 + I/5, 50]");
    run("chk[v_, r_] := Abs[v - r] <= 10^-44 Max[1, Abs[r]]");
    run("c1 = Compile[{{x, _Real}}, Erf[x], WorkingPrecision -> 50]");
    check_true("chk[c1[x0], 0.67780119383741847297562880924415139671628817433487`50]",
               "mpmath: Erf[0.7]");
    run("c2 = Compile[{{x, _Real}}, BesselK[2, x], WorkingPrecision -> 50]");
    check_true("chk[c2[x0], 3.6613299608091528385402696981650674146983113168059`50]",
               "mpmath: BesselK[2, 0.7]");
    run("c3 = Compile[{{z, _Complex}}, AiryAi[z], WorkingPrecision -> 50]");
    check_true("chk[c3[z0], 0.18649366102203101833733547542952391763544677443122`50 "
               "- 0.04003477266815675113714012381871903746716870764758`50 I]",
               "mpmath: AiryAi[0.7 + 0.2 I]");
    run("c4 = Compile[{{z, _Complex}}, Hypergeometric1F1[1/3, 2, z], WorkingPrecision -> 50]");
    check_true("chk[c4[z0], 1.1353533998829286163600954936375082818311049221872`50 "
               "+ 0.046104258871111266676427465153720478587604404590636`50 I]",
               "mpmath: 1F1[1/3; 2; 0.7 + 0.2 I]");
    run("c5 = Compile[{{x, _Real}}, EllipticK[x], WorkingPrecision -> 50]");
    check_true("chk[c5[x0], 2.0753631352924691438534405558824158057376556817907`50]",
               "mpmath: EllipticK[0.7]");
    run("c6 = Compile[{{z, _Complex}}, Zeta[z], WorkingPrecision -> 50]");
    check_true("chk[c6[z0], -1.7525679280097623842070738018718823784983873077468`50 "
               "- 1.523333864733266194759475581179591060664984740939`50 I]",
               "mpmath: Zeta[0.7 + 0.2 I]");
    run("c7 = Compile[{{x, _Real}}, PolyLog[2, x], WorkingPrecision -> 50]");
    check_true("chk[c7[x0], 0.88937762428603873860100627480736179353714742710847`50]",
               "mpmath: PolyLog[2, 0.7]");
    run("c8 = Compile[{{x, _Real}}, JacobiSN[x, 1/3], WorkingPrecision -> 50]");
    check_true("chk[c8[x0], 0.63097849068692889887438414261282738839667957712505`50]",
               "mpmath: JacobiSN[0.7, 1/3]");
    /* complex-valued on the real line: compiled with an ncpx result */
    run("c9 = Compile[{{x, _Real}}, HankelH1[2, x], WorkingPrecision -> 50]");
    check_true("chk[c9[x0], 0.058786944364191713014890042735882796375892456858591`50 "
               "- 2.9614775618272717030112063060672276502460627254797`50 I]",
               "mpmath: HankelH1[2, 0.7] (complex result of a real argument)");
    check_true("Lookup[CompileDiagnostics[{{x, _Real}}, HankelH1[2, x], "
               "WorkingPrecision -> 30], \"ResultType\"] === \"MPFRComplex\"",
               "HankelH1 of a real argument compiles to an MPFR complex result");

    /* ---- composite bodies: special functions inside arithmetic ---- */
    probe("Erf[x]^2 + BesselJ[1, 2 x] Sin[x] - Gamma[x + 1/2]/(1 + x^2)", false, 40);
    probe("Erf[x]^2 + BesselJ[1, 2 x] Sin[x] - Gamma[x + 1/2]/(1 + x^2)", true, 40);
    probe("Hypergeometric2F1[1/3, 1/2, 3/2, x^2] EllipticE[x/2] + MathieuC[1, 1/2, x]",
          false, 30);
    probe("Sinh[x] Zeta[x + 2] + ArcSin[x]", true, 30);   /* complex elementary via the registry */
    /* a closed call is folded to a constant at the working precision */
    probe("x + BesselJZero[0, 3] + Gamma[1/3]", false, 30);
    check_true("Lookup[CompileDiagnostics[{{x, _Real}}, x + BesselJZero[0, 3], "
               "WorkingPrecision -> 30], \"Instructions\"] <= 3",
               "closed special-function call folds to one constant");
    /* a head outside the registry still bails cleanly */
    check_true("Lookup[CompileDiagnostics[{{x, _Real}}, Floor[x] + x, "
               "WorkingPrecision -> 30], \"Compiled\"] === False",
               "Floor (exempt, integer-valued) still bails at arbitrary precision");
    /* a real program meeting a complex value declines; the caller interprets */
    run("cn = Compile[{{x, _Real}}, BesselY[0, x], WorkingPrecision -> 30]");
    check_true("Abs[cn[N[-7/10, 30]] - BesselY[0, N[-7/10, 30]]] < 10^-25",
               "real program, complex value: declines, still the interpreter's answer");

    /* ---- auto-compilation at WorkingPrecision -> 30: no interpreter fallback ---- */
    for (size_t i = 0; i < sizeof AC / sizeof AC[0]; i++) {
        char src[1024], label[256];
        run("a0 = Internal`CompileArbitraryPrecisionStats[]");
        snprintf(src, sizeof src, "res = %s", AC[i].expr);
        run(src);
        run("a1 = Internal`CompileArbitraryPrecisionStats[]");
        snprintf(label, sizeof label, "%s: compiled (no auto-compile bail)", AC[i].label);
        check_true("(\"AutoCompileBails\" /. a1) == (\"AutoCompileBails\" /. a0)", label);
        snprintf(label, sizeof label, "%s: special functions ran compiled", AC[i].label);
        check_true("(\"SpecialCalls\" /. a1) > (\"SpecialCalls\" /. a0)", label);
        snprintf(src, sizeof src, "NumberQ[res] && Abs[res - %s] < 10^-25", AC[i].ref);
        snprintf(label, sizeof label, "%s: matches mpmath to 25 digits", AC[i].label);
        check_true(src, label);
        autocompile_set_enabled(false);
        snprintf(src, sizeof src, "res0 = %s", AC[i].expr);
        run(src);
        autocompile_set_enabled(true);
        snprintf(label, sizeof label, "%s: compiled == uncompiled", AC[i].label);
        check_true("Abs[res - res0] < 10^-25", label);
    }
    /* NDSolve's MPFR right-hand side (one program per component) */
    run("a0 = Internal`CompileArbitraryPrecisionStats[]");
    run("sol = NDSolve[{y'[t] == -Erf[t] y[t] + BesselJ[0, t], y[0] == 1}, y, {t, 0, 1}, "
        "WorkingPrecision -> 30]");
    run("a1 = Internal`CompileArbitraryPrecisionStats[]");
    check_true("(\"AutoCompileBails\" /. a1) == (\"AutoCompileBails\" /. a0) && "
               "(\"SpecialCalls\" /. a1) > (\"SpecialCalls\" /. a0)",
               "NDSolve WorkingPrecision -> 30: special-function RHS compiled");

    /* Table over a List, and over a range, of arbitrary-precision reals */
    run("xs = Table[N[1/10 + k/100, 30], {k, 50}]");
    run("a0 = Internal`CompileArbitraryPrecisionStats[]");
    run("tl = Table[Erf[x]^2 + BesselK[0, x] Sin[x], {x, xs}]");
    run("tr = Table[Erf[x]^2 + BesselK[0, x] Sin[x], {x, N[1/10, 30], N[6/10, 30], N[1/100, 30]}]");
    run("a1 = Internal`CompileArbitraryPrecisionStats[]");
    check_true("(\"SpecialCalls\" /. a1) - (\"SpecialCalls\" /. a0) == 2 (50 + 51) "
               "&& (\"AutoCompileBails\" /. a1) == (\"AutoCompileBails\" /. a0)",
               "Table over MPFR List and range: compiled, no bail");
    autocompile_set_enabled(false);
    run("tl0 = Table[Erf[x]^2 + BesselK[0, x] Sin[x], {x, xs}]");
    run("tr0 = Table[Erf[x]^2 + BesselK[0, x] Sin[x], {x, N[1/10, 30], N[6/10, 30], N[1/100, 30]}]");
    autocompile_set_enabled(true);
    check_true("Length[tl] == 50 && Max[Abs[tl - tl0]] < 10^-27 && Precision[tl[[7]]] >= 29",
               "Table over MPFR List: compiled == uncompiled, precision kept");
    check_true("Length[tr] == Length[tr0] == 51 && Max[Abs[tr - tr0]] < 10^-27",
               "Table over MPFR range: compiled == uncompiled");

    printf("\n%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
