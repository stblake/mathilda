#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/*
 * simp.c -- Simplify, Assuming, $Assumptions, and AssumeCtx.
 *
 * This translation unit is the module entry point (simp_init). The actual
 * implementation has been split into a family of sibling .c files (see
 * simp_internal.h for the cross-module surface):
 *
 *   simp_complexity.c       Default complexity measure + SimplifyCount
 *   simp_assume.c           AssumeCtx + domain queries
 *   simp_util.c             Helpers, $SimplifyDebug tracing, CandSet, scoring
 *   simp_assume_rewrite.c   Assumption-driven seed rewriters
 *   simp_trig_roundtrip.c   Trig/exp roundtrip composite
 *   simp_denest.c           Radical canon + sqrt/multi-extension denesting
 *   simp_cuberoot.c         Cube-root denesting
 *   simp_rationalize.c      Denominator rationalisation
 *   simp_algebraic.c        Algebraic-extension reduction
 *   simp_canon.c            Common-factor lift + sign canon
 *   simp_trig_pi.c          Trig at rational Pi + Pythag transforms
 *   simp_power.c            PrimeRebase, PowerOneify, PowerDistribute, RadicalCanon
 *   simp_tan_add.c          TanAddition
 *   simp_logexp_abs.c       Log/Power rewriter + Abs + Sqrt[e^2]
 *   simp_search.c           Heuristic search, shape classifier, pipelines, additive splitter
 *   simp_factorial.c        simp_factorial
 *   simp_bottomup.c         SimpMemo + bottom-up driver
 *   simp_builtins.c         read_$Assumptions, Element, rebalancing, Simplify, Assuming
 *
 * Simplify implements a small heuristic search over the existing battery
 * of algebraic transforms. The default complexity measure is
 * LeafCount(expr) + decimal-digit count of integer leaves; this matches
 * Mathematica's default and stops e.g. "100 Log[2]" from being rewritten
 * to "Log[2^100]". A user-supplied ComplexityFunction option overrides
 * this. See simp.h for the AssumeCtx contract.
 *
 * Assuming desugars to Block[{$Assumptions = $Assumptions && a}, body],
 * which reuses Block's existing scope-restoration code path. Nested
 * Assuming calls compose because each Block reads the current
 * $Assumptions OwnValue before extending it.
 */

/* ----------------------------------------------------------------------- */
/* simp_init                                                               */
/* ----------------------------------------------------------------------- */

void simp_init(void) {
    /* $Assumptions defaults to True. */
    Expr* dollar_pat = expr_new_symbol("$Assumptions");
    Expr* dollar_val = expr_new_symbol("True");
    symtab_add_own_value("$Assumptions", dollar_pat, dollar_val);
    expr_free(dollar_pat);
    expr_free(dollar_val);

    /* $SimplifyDebug defaults to False. When set to True, Simplify emits
     * one stderr line per transform invocation in the form
     *   /<TransformName>/: <input> -> <output> [<elapsed> ms]
     * to help diagnose hangs and runaway candidate explosion. */
    Expr* dbg_pat = expr_new_symbol("$SimplifyDebug");
    Expr* dbg_val = expr_new_symbol("False");
    symtab_add_own_value("$SimplifyDebug", dbg_pat, dbg_val);
    expr_free(dbg_pat);
    expr_free(dbg_val);
    symtab_set_docstring("$SimplifyDebug",
        "$SimplifyDebug\n\tWhen set to True, Simplify prints one stderr line per\n"
        "\ttransform invocation: /Name/: <input> -> <output> [<ms> ms].\n"
        "\tDefaults to False. Useful for diagnosing slow Simplify calls.");

    symtab_add_builtin("Simplify", builtin_simplify);
    symtab_get_def("Simplify")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("SimplifyCount", builtin_simplify_count);
    symtab_get_def("SimplifyCount")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_set_docstring("SimplifyCount",
        "SimplifyCount[expr]\n\tThe complexity measure used by Simplify when no\n"
        "\tComplexityFunction option (or ComplexityFunction -> Automatic) is\n"
        "\tgiven. Counts subexpressions; integers contribute their decimal\n"
        "\tdigit count plus a constant for the sign. Real numbers contribute\n"
        "\t2 (NumberQ but not Integer/Rational).");

    symtab_add_builtin("Assuming", builtin_assuming);
    symtab_get_def("Assuming")->attributes |= (ATTR_HOLDREST | ATTR_PROTECTED);

    symtab_add_builtin("Element", builtin_element);
    symtab_get_def("Element")->attributes |= ATTR_PROTECTED;
}
