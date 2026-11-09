/* Mathilda — Compile[]: arbitrary-precision special functions.  See
 * compile_mgd_special.h and docs/design/compile-arbitrary-precision.md §19.
 *
 * WHY ONE GENERIC OPCODE AND NOT ONE PER HEAD.  Every head below already has an
 * arbitrary-precision evaluator (arb / MPFR / an expression template) behind its
 * builtin, and that evaluator is what the interpreter runs at a WorkingPrecision
 * sample point.  Re-binding each to a raw C signature would duplicate ~130
 * dispatchers and — worse — let the compiled value drift from the interpreted
 * one whenever a head's evaluator is improved.  Instead OP_MG_SPECIAL boxes the
 * operand containers into numeric Exprs at the program precision and evaluates
 * f[args] through the ordinary evaluator: the SAME routine, so compiled ==
 * uncompiled by construction.  The boxing is a few allocations against an arb
 * evaluation that costs microseconds at 30+ digits, and what the compiled
 * program removes — the tree-walk of the REST of the body (arithmetic,
 * elementary functions, every other subexpression) on every sample — is kept.
 *
 * The table, not the opcode, decides what lowers: a head is lowered only if it
 * is listed here, at an arity it lists, with complex operands only where its
 * evaluator takes them.  Anything else bails at emit time exactly as before. */
#include "compile_mgd_special.h"
#include "compile_internal.h"      /* mgd_mpfr_from_expr */
#include "../eval.h"
#include "../sym_intern.h"
#include "../arithmetic.h"        /* is_complex, is_rational */
#include "../symtab.h"
#include "../attr.h"
#include "../sym_names.h"
#include "autocompile.h"         /* autocompile_prec_stats */
#ifdef USE_MPFR
#include "../numeric_complex.h"   /* ncpx, numeric_mpfr_make_complex */
#endif
#include <stdlib.h>
#include <string.h>

#define R_  MGS_REAL
#define RC_ (MGS_REAL | MGS_CPLX)
#define RCO_ (MGS_REAL | MGS_CPLX | MGS_COUT)

/* name, min arity, max arity, flags.  Complex support was measured per head
 * (N[f[reals..., 7/10 + I/5], 30] returns an arbitrary-precision number);
 * the real-only rows are the heads whose evaluator leaves a complex argument
 * unevaluated.  Kept alphabetical so the audit's diff reads cleanly. */
static const MgdSpecialEntry g_table[] = {
    { "AiryAi",                 1, 1, RC_ },
    { "AiryAiPrime",            1, 1, RC_ },
    { "AiryBi",                 1, 1, RC_ },
    { "AiryBiPrime",            1, 1, RC_ },
    { "AngerJ",                 2, 3, RC_ },
    { "ArcCos",                 1, 1, RC_ },   /* complex: elementary rows   */
    { "ArcCosh",                1, 1, RC_ },   /* are only reached where the */
    { "ArcCot",                 1, 1, RC_ },   /* mpfr / ncpx opcode set has */
    { "ArcCoth",                1, 1, RC_ },   /* no direct lowering         */
    { "ArcCsc",                 1, 1, RC_ },
    { "ArcCsch",                1, 1, RC_ },
    { "ArcSec",                 1, 1, RC_ },
    { "ArcSech",                1, 1, RC_ },
    { "ArcSin",                 1, 1, RC_ },
    { "ArcSinh",                1, 1, RC_ },
    { "ArcTan",                 1, 1, RC_ },
    { "ArcTan",                 2, 2, R_  },   /* ArcTan[x, y]: a real pair only */
    { "ArcTanh",                1, 1, RC_ },
    { "Arg",                    1, 1, RC_ },   /* complex Arg is native   */
    { "BarnesG",                1, 1, RC_ },
    { "BernoulliB",             2, 2, RC_ },
    { "BesselI",                2, 2, RC_ },
    { "BesselJ",                2, 2, RC_ },
    { "BesselK",                2, 2, RC_ },
    { "BesselY",                2, 2, RC_ },
    { "Beta",                   2, 4, RC_ },
    { "BetaRegularized",        3, 4, RC_ },
    { "Binomial",               2, 2, RC_ },
    { "ChebyshevT",             2, 2, RC_ },
    { "ChebyshevU",             2, 2, RC_ },
    { "Conjugate",              1, 1, RC_ },
    { "CosIntegral",            1, 1, RC_ },
    { "Cosh",                   1, 1, RC_ },
    { "CoshIntegral",           1, 1, RC_ },
    { "Coth",                   1, 1, RC_ },
    { "Csch",                   1, 1, RC_ },
    { "CubeRoot",               1, 1, R_  },
    { "DirichletL",             3, 3, RCO_ },  /* complex characters */
    { "EllipticE",              1, 2, RC_ },
    { "EllipticF",              2, 2, RC_ },
    { "EllipticK",              1, 1, RC_ },
    { "EllipticPi",             2, 3, RC_ },
    { "EllipticTheta",          3, 3, RC_ },
    { "EllipticThetaPrime",     3, 3, RC_ },
    { "Erf",                    1, 2, RC_ },
    { "Erfc",                   1, 1, RC_ },
    { "Erfi",                   1, 1, RC_ },
    { "EulerE",                 2, 2, RC_ },
    { "ExpIntegralE",           2, 2, RC_ },
    { "ExpIntegralEi",          1, 1, RC_ },
    { "Factorial",              1, 1, RC_ },
    { "Factorial2",             1, 1, RC_ },
    { "FactorialPower",         2, 3, RC_ },
    { "Fibonacci",              1, 2, RC_ },
    { "FresnelC",               1, 1, RC_ },
    { "FresnelS",               1, 1, RC_ },
    { "Gamma",                  1, 3, RC_ },
    { "GegenbauerC",            2, 3, RC_ },
    { "HankelH1",               2, 2, RCO_ },
    { "HankelH2",               2, 2, RCO_ },
    { "HarmonicNumber",         1, 2, RC_ },
    { "HermiteH",               2, 2, RC_ },
    { "HurwitzZeta",            2, 2, RC_ },
    { "Hyperfactorial",         1, 1, RC_ },
    { "Hypergeometric0F1",      2, 2, RC_ },
    { "Hypergeometric1F1",      3, 3, RC_ },
    { "Hypergeometric2F1",      4, 4, RC_ },
    { "HypergeometricPFQ",      3, 3, RC_ },   /* literal parameter Lists */
    { "HypergeometricU",        3, 3, RC_ },
    { "Im",                     1, 1, RC_ },
    { "InverseErf",             1, 1, R_  },
    { "InverseErfc",            1, 1, R_  },
    { "InverseJacobiCN",        2, 2, RC_ },
    { "InverseJacobiDN",        2, 2, RC_ },
    { "InverseJacobiSN",        2, 2, RC_ },
    { "JacobiAmplitude",        2, 2, RC_ },
    { "JacobiCD",               2, 2, RC_ },
    { "JacobiCN",               2, 2, RC_ },
    { "JacobiCS",               2, 2, RC_ },
    { "JacobiDC",               2, 2, RC_ },
    { "JacobiDN",               2, 2, RC_ },
    { "JacobiDS",               2, 2, RC_ },
    { "JacobiNC",               2, 2, RC_ },
    { "JacobiND",               2, 2, RC_ },
    { "JacobiNS",               2, 2, RC_ },
    { "JacobiP",                4, 4, RC_ },
    { "JacobiSC",               2, 2, RC_ },
    { "JacobiSD",               2, 2, RC_ },
    { "JacobiSN",               2, 2, RC_ },
    { "KelvinBei",              1, 2, RC_ },
    { "KelvinBer",              1, 2, RC_ },
    { "KelvinKei",              1, 2, RC_ },
    { "KelvinKer",              1, 2, RC_ },
    { "LaguerreL",              2, 3, RC_ },
    { "LegendreP",              2, 4, RC_ },   /* [n, m, type, x]: literal type */
    { "LegendreQ",              2, 4, RC_ },
    { "LerchPhi",               3, 3, RC_ },
    { "Log",                    1, 2, RC_ },   /* Log[b, x]; Log[x] complex     */
    { "Log10",                  1, 1, RC_ },
    { "Log2",                   1, 1, RC_ },
    { "LogGamma",               1, 1, RC_ },
    { "LogIntegral",            1, 1, RC_ },
    { "LucasL",                 1, 2, RC_ },
    { "MathieuC",               3, 3, RC_ },
    { "MathieuCPrime",          3, 3, RC_ },
    { "MathieuCharacteristicA", 2, 2, R_  },
    { "MathieuCharacteristicB", 2, 2, R_  },
    { "MathieuS",               3, 3, RC_ },
    { "MathieuSPrime",          3, 3, RC_ },
    { "MeijerG",                3, 4, RC_ },   /* literal parameter Lists */
    { "Multinomial",            1, 16, RC_ },
    { "ParabolicCylinderD",     2, 2, RC_ },
    { "Pochhammer",             2, 2, RC_ },
    { "PolyGamma",              1, 2, RC_ },
    { "PolyLog",                2, 3, RC_ },   /* PolyLog[n, p, z]: literal n, p */
    { "ProductLog",             1, 2, RC_ },   /* ProductLog[k, z]: literal k    */
    { "QPochhammer",            2, 3, RC_ },
    { "Re",                     1, 1, RC_ },
    { "RiemannR",               1, 1, RC_ },
    { "RiemannSiegelTheta",     1, 1, RC_ },
    { "RiemannSiegelZ",         1, 1, RC_ },
    { "Sech",                   1, 1, RC_ },
    { "Sinc",                   1, 1, RC_ },
    { "Sinh",                   1, 1, RC_ },
    { "SinIntegral",            1, 1, RC_ },
    { "SinhIntegral",           1, 1, RC_ },
    { "SphericalBesselJ",       2, 2, RC_ },
    { "SphericalBesselY",       2, 2, RC_ },
    { "SphericalHankelH1",      2, 2, RCO_ },
    { "SphericalHankelH2",      2, 2, RCO_ },
    { "StieltjesGamma",         1, 2, RC_ },   /* literal index n               */
    { "StruveH",                2, 2, RC_ },
    { "StruveL",                2, 2, RC_ },
    { "Surd",                   2, 2, R_  },
    { "Tanh",                   1, 1, RC_ },
    { "WeberE",                 2, 3, RC_ },
    { "WeierstrassP",           2, 2, RC_ },   /* literal invariants {g2, g3}   */
    { "WeierstrassPPrime",      2, 2, RC_ },
    { "WeierstrassSigma",       2, 2, RC_ },
    { "WeierstrassZeta",        2, 2, RC_ },
    { "WhittakerM",             3, 3, RC_ },
    { "WhittakerW",             3, 3, RC_ },
    { "Zeta",                   1, 2, RC_ },
};
#define NTABLE (sizeof g_table / sizeof g_table[0])

/* Interned names, filled on first lookup (compile time, single-threaded). */
static const char* g_interned[NTABLE];
static bool g_interned_ready = false;

static void ensure_interned(void) {
    if (g_interned_ready) return;
    for (size_t i = 0; i < NTABLE; i++) g_interned[i] = intern_symbol(g_table[i].name);
    g_interned_ready = true;
}

const MgdSpecialEntry* mgd_special_lookup(const char* h, size_t nargs) {
    if (!h) return NULL;
    ensure_interned();
    for (size_t i = 0; i < NTABLE; i++)
        if (g_interned[i] == h && nargs >= g_table[i].min_args
            && nargs <= g_table[i].max_args)
            return &g_table[i];       /* a head may have one row per arity range */
    return NULL;
}

const MgdSpecialEntry* mgd_special_table(size_t* n) {
    if (n) *n = NTABLE;
    return g_table;
}

void mgd_special_call_free(MgdSpecialCall* sc) {
    if (!sc) return;
    if (sc->lit)
        for (size_t i = 0; i < sc->nargs; i++) expr_free(sc->lit[i]);
    free(sc->lit);
    free(sc->cplx);
    free(sc);
}

static long g_calls = 0, g_declines = 0;

void mgd_special_stats(long* calls, long* declines) {
    if (calls) *calls = g_calls;
    if (declines) *declines = g_declines;
}

#ifdef USE_MPFR
/* An exact or arbitrary-precision real part: an EXPR_MPFR, or an exact
 * integer / rational (the evaluator answers BesselJ[1, 0.`30] with 0).  A
 * machine Real is refused: it means precision was lost, and the interpreter's
 * downstream arithmetic would then run at machine precision too, which the
 * fixed-precision program cannot mirror. */
static bool part_to_mpfr(mpfr_ptr d, const Expr* e) {
    if (e->type == EXPR_REAL) return false;
    return mgd_mpfr_from_expr(d, e);   /* MPFR, Integer, BigInt, Rational */
}

bool mgd_special_invoke(const MgdSpecialCall* sc, void* const* ops,
                        void* dst, bool dst_cplx) {
    g_calls++;
    Expr** av = malloc((sc->nargs ? sc->nargs : 1) * sizeof *av);
    if (!av) { g_declines++; return false; }
    size_t j = 0;
    for (size_t i = 0; i < sc->nargs; i++) {
        if (sc->lit[i]) { av[i] = expr_copy(sc->lit[i]); continue; }
        if (sc->cplx[i]) {
            const ncpx* z = (const ncpx*)ops[j++];
            av[i] = numeric_mpfr_make_complex(z->re, z->im);
        } else {
            av[i] = expr_new_mpfr_copy((mpfr_srcptr)ops[j++]);
        }
    }
    Expr* call = expr_new_function(expr_new_symbol(sc->head), av, sc->nargs);
    free(av);
    /* Through the ordinary evaluator, not the builtin pointer: the evaluator's
     * dispatch around the builtin is part of what N[f[...], n] does — the
     * exact-parameter retry (numeric_mixed_args_call), DownValues, the numeric
     * call memo, and the dependency bookkeeping that keeps an enclosing fixed
     * point honest.  Calling the builtin directly saved ~8% of a lone call and
     * would have skipped all of it. */
    Expr* r = evaluate(call);
    expr_free(call);
    bool ok = false;
    if (r) {
        Expr *re, *im;
        if (dst_cplx) {
            ncpx* d = (ncpx*)dst;
            if (is_complex(r, &re, &im))
                ok = part_to_mpfr(d->re, re) && part_to_mpfr(d->im, im);
            else if ((ok = part_to_mpfr(d->re, r)))
                mpfr_set_ui(d->im, 0, MPFR_RNDN);
            ok = ok && mpfr_number_p(d->re) && mpfr_number_p(d->im);
        } else {
            mpfr_ptr d = (mpfr_ptr)dst;
            /* a genuinely complex value where a real was compiled declines,
             * and the caller interprets the point (the _z program, or the
             * interpreter, then supplies it) */
            if (is_complex(r, &re, &im))
                ok = (im->type == EXPR_INTEGER && im->data.integer == 0)
                  && part_to_mpfr(d, re);
            else
                ok = part_to_mpfr(d, r);
            ok = ok && mpfr_number_p(d);
        }
        expr_free(r);
    }
    if (!ok) g_declines++;
    return ok;
}
#endif /* USE_MPFR */

/* ---- introspection (tests and tools/compile_coverage.py --precision) ---- */

static Expr* mk_rule(const char* k, Expr* v) {
    Expr* a[2] = { expr_new_string(k), v };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* Internal`CompileArbitraryPrecisionHeads[] — the registry as
 * {{"head", minArgs, maxArgs, complexQ}, ...}. */
static Expr* builtin_ap_heads(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    Expr** rows = malloc(NTABLE * sizeof *rows);
    if (!rows) return NULL;
    for (size_t i = 0; i < NTABLE; i++) {
        Expr* f[4] = {
            expr_new_string(g_table[i].name),
            expr_new_integer(g_table[i].min_args),
            expr_new_integer(g_table[i].max_args),
            expr_new_symbol((g_table[i].flags & MGS_CPLX) ? SYM_True : SYM_False),
        };
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), f, 4);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, NTABLE);
    free(rows);
    return out;
}

/* Internal`CompileArbitraryPrecisionStats[] — counters since startup. */
static Expr* builtin_ap_stats(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    long bails = 0, declines = 0;
    autocompile_prec_stats(&bails, &declines);
    Expr* r[5] = {
        mk_rule("SpecialCalls", expr_new_integer(g_calls)),
        mk_rule("SpecialDeclines", expr_new_integer(g_declines)),
        mk_rule("EvalDeclines", expr_new_integer(compiled_managed_eval_failures())),
        mk_rule("AutoCompileBails", expr_new_integer(bails)),
        mk_rule("AutoCompileDeclines", expr_new_integer(declines)),
    };
    return expr_new_function(expr_new_symbol(SYM_List), r, 5);
}

void compile_mgd_special_init(void) {
    symtab_add_builtin("Internal`CompileArbitraryPrecisionHeads", builtin_ap_heads);
    symtab_get_def("Internal`CompileArbitraryPrecisionHeads")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Internal`CompileArbitraryPrecisionHeads",
        "Internal`CompileArbitraryPrecisionHeads[] lists the special functions that "
        "Compile[..., WorkingPrecision -> n] and arbitrary-precision auto-compilation "
        "lower to their own arbitrary-precision evaluator, as {head, minArgs, "
        "maxArgs, complexQ} rows.");
    symtab_add_builtin("Internal`CompileArbitraryPrecisionStats", builtin_ap_stats);
    symtab_get_def("Internal`CompileArbitraryPrecisionStats")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Internal`CompileArbitraryPrecisionStats",
        "Internal`CompileArbitraryPrecisionStats[] gives counters since startup: "
        "compiled special-function calls run and declined, arbitrary-precision "
        "compiled calls that declined to the interpreter, arbitrary-precision "
        "auto-compilations that fell back to the interpreter, and auto-compiled "
        "sample points that declined to the interpreter.");
}
