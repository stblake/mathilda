/* TimeConstrained uses sigaction / sigsetjmp / siglongjmp / setitimer,
 * none of which are in C99.  These feature-test macros expose them under
 * glibc's -std=c99; Darwin exposes them implicitly.  Must come before
 * the first #include. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "core.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "arithmetic.h"
#include "numbertheory.h"
#include "gamma.h"
#include "beta.h"
#include "erf.h"
#include "erfc.h"
#include "erfi.h"
#include "expintegralei.h"
#include "logintegral.h"
#include "sinintegral.h"
#include "cosintegral.h"
#include "sinhintegral.h"
#include "coshintegral.h"
#include "fresnel.h"
#include "sinc.h"
#include "inverf.h"
#include "inverfc.h"
#include "loggamma.h"
#include "polygamma.h"
#include "harmonicnumber.h"
#include "pochhammer.h"
#include "eulergamma.h"
#include "zeta.h"
#include "hurwitzzeta.h"
#include "stieltjesgamma.h"
#include "bernoullib.h"
#include "eulere.h"
#include "polylog.h"
#include "lerchphi.h"
#include "productlog.h"
#include "legendre.h"
#include "airyai.h"
#include "airybi.h"
#include "bessel.h"
#include "comparisons.h"
#include "boolean.h"
#include "names.h"
#include "list.h"
#include "assoc.h"
#include "ndarray.h"
#include "replace.h"
#include "patterns.h"
#include "cond.h"
#include "iter.h"
#include "complex.h"
#include "trig.h"
#include "hyperbolic.h"
#include "logexp.h"
#include "piecewise.h"
#include "int.h"
#include "real.h"
#include "attr.h"
#include "purefunc.h"
#include "modular.h"
#include "sort.h"
#include "stats.h"
#include "partitions.h"
#include "fit.h"
#include "info.h"
#include "expand.h"
#include "expand_power.h"
#include "poly.h"
#include "rat.h"
#include "facint.h"
#include "fibonacci.h"
#include "lucas.h"
#include "facpoly.h"
#include "flint_bridge.h"
#include "solve.h"
#include "findroot.h"
#include "datetime.h"
#include "linalg.h"
#include "readwrite.h"
#include "loadmodule.h"
#include "files.h"
#include "part.h"
#include "plus.h"
#include "times.h"
#include "power.h"
#include "funcprog.h"
#include "match.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>
#include <setjmp.h>
#include <sys/time.h>
#include <time.h>
#include <float.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include "rat.h"
#include "parfrac.h"
#include "contfrac.h"
#include "random.h"
#include "picostrings.h"
#include "series.h"
#include "deriv.h"
#include "limit.h"
#include "context.h"
#include "numeric.h"
#include "precision.h"
#include "rationalize.h"
#include "numeric.h"
#include "linsolve.h"
#include "common.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "repl_hooks.h"
#include "version.h"

/*
 * register_system_constant:
 * Helper for system_constants_init below — bind `name` to `value` as an
 * OwnValue, mark the symbol Protected, and free both helper Exprs. The
 * value is now owned by the symbol table.
 */
static void register_system_constant(const char* name, Expr* value) {
    Expr* sym = expr_new_symbol(name);
    symtab_add_own_value(name, sym, value);
    expr_free(sym);
    expr_free(value);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
}

/*
 * system_constants_init:
 * Registers the read-only floating-point system constants:
 *   $MachinePrecision, $MachineEpsilon, $MinMachineNumber, $MaxMachineNumber,
 *   $MaxNumber, $MinNumber
 *
 * Values are derived from the platform's <float.h> macros (and MPFR's
 * current exponent range, when USE_MPFR is enabled) rather than hard-coded,
 * so they reflect the actual representation Mathilda was compiled against.
 */
static void system_constants_init(void) {
    /* Machine-precision quantities — exact IEEE 754 limits of the local
     * double type. */
    register_system_constant("$MachinePrecision", expr_new_real(NUMERIC_MACHINE_PRECISION_DIGITS));
    /* $MaxExtraPrecision bounds the auto-precision a numeric search (e.g.
     * FindIntegerNullVector) may add beyond $MachinePrecision.  Unlike the
     * other constants it is user-settable, so it is NOT marked Protected. */
    {
        Expr* sym = expr_new_symbol(SYM_DollarMaxExtraPrecision);
        symtab_add_own_value("$MaxExtraPrecision", sym, expr_new_real(50.0));
        expr_free(sym);
    }
    register_system_constant("$MachineEpsilon",   expr_new_real(DBL_EPSILON));
    register_system_constant("$MinMachineNumber", expr_new_real(DBL_MIN));
    register_system_constant("$MaxMachineNumber", expr_new_real(DBL_MAX));

    /* $MaxNumber / $MinNumber describe the arbitrary-precision range. With
     * MPFR available we walk one step in from the current emax/emin bounds
     * at machine precision (~53 bits). Without MPFR there is no arb-prec
     * representation, so the values collapse onto the machine extrema. */
#ifdef USE_MPFR
    {
        mpfr_prec_t bits = (mpfr_prec_t)DBL_MANT_DIG;
        mpfr_t maxv, minv;
        mpfr_init2(maxv, bits);
        mpfr_set_inf(maxv, +1);
        mpfr_nextbelow(maxv);  /* largest finite value at this precision */
        register_system_constant("$MaxNumber", expr_new_mpfr_move(maxv));

        mpfr_init2(minv, bits);
        mpfr_set_zero(minv, +1);
        mpfr_nextabove(minv);  /* smallest positive finite value */
        register_system_constant("$MinNumber", expr_new_mpfr_move(minv));
    }
#else
    register_system_constant("$MaxNumber", expr_new_real(DBL_MAX));
    register_system_constant("$MinNumber", expr_new_real(DBL_MIN));
#endif

    /* Release identity. $VersionNumber is the single source of truth (a Real);
     * $Version is the descriptive string assembled at compile time in
     * version.c, listing Mathilda's version and every library it links. Both
     * are read-only (Protected via register_system_constant). */
    register_system_constant("$Version", expr_new_string(mathilda_version()));
    register_system_constant("$VersionNumber", expr_new_real(MATHILDA_VERSION_NUMBER));
    symtab_set_docstring("$Version",
        "$Version\n\tgives a string describing the version of Mathilda, "
        "including the versions of the libraries it was built against.");
    symtab_set_docstring("$VersionNumber",
        "$VersionNumber\n\tgives the Mathilda version number as a real number.");
}

void core_init(void) {
    /* Route GMP/MPFR allocation through the guarded wrappers before anything
     * else can allocate a bignum, so TimeConstrained's async abort can never
     * unwind out of a held libmalloc lock (see tc_install_alloc_guard). */
    tc_install_alloc_guard();
    /* Cache canonical pointers for well-known symbol names. Doing this
     * before any other init step guarantees every later expr_new_symbol
     * with one of these names returns the same pointer the SYM_*
     * constants hold, so eval-time pointer comparisons are valid. */
    sym_names_init();
    /* Context system must come first: the parser calls context_resolve_name
     * on every identifier, including those produced while loading init.m. */
    context_init();
    /* Seed user-visible system variables ($RecursionLimit, ...) early so
     * any subsequent init step that triggers evaluation can already see
     * them. */
    eval_init();
    system_constants_init();
    repl_hooks_init();
    parfrac_init();
    contfrac_init();
    modular_init();
    symtab_add_builtin("AtomQ", builtin_atomq);
    symtab_add_builtin("Identity", builtin_identity);
    symtab_get_def("Identity")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Composition", builtin_composition);
    symtab_get_def("Composition")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_add_builtin("ComposeList", builtin_compose_list);
    symtab_get_def("ComposeList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Length", builtin_length);
    symtab_add_builtin("Dimensions", builtin_dimensions);
    symtab_add_builtin("Clear", builtin_clear);
    symtab_add_builtin("Unset", builtin_unset);
    symtab_add_builtin("ClearAll", builtin_clear_all);
    symtab_add_builtin("Remove", builtin_remove);
    symtab_add_builtin("Protect", builtin_protect);
    symtab_add_builtin("Unprotect", builtin_unprotect);
    symtab_add_builtin("Part", builtin_part);
    symtab_add_builtin("Extract", builtin_extract);
    symtab_add_builtin("Head", builtin_head);
    symtab_add_builtin("First", builtin_first);
    symtab_add_builtin("Last", builtin_last);
    symtab_add_builtin("Most", builtin_most);
    symtab_add_builtin("Rest", builtin_rest);

    symtab_get_def("Clear")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Unset")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("ClearAll")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Remove")->attributes |= ATTR_HOLDALL | ATTR_LOCKED | ATTR_PROTECTED;
    symtab_get_def("Protect")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Unprotect")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_get_def("Part")->attributes |= ATTR_NHOLDREST | ATTR_PROTECTED;
    symtab_get_def("Extract")->attributes |= ATTR_NHOLDREST | ATTR_PROTECTED;
    symtab_add_builtin("Insert", builtin_insert);
    symtab_add_builtin("Delete", builtin_delete);
    symtab_add_builtin("Append", builtin_append);
    symtab_add_builtin("Prepend", builtin_prepend);
    symtab_add_builtin("AppendTo", builtin_append_to);
    symtab_add_builtin("PrependTo", builtin_prepend_to);
    symtab_add_builtin("Increment", builtin_increment);
    symtab_add_builtin("Decrement", builtin_decrement);
    symtab_add_builtin("PreIncrement", builtin_preincrement);
    symtab_add_builtin("PreDecrement", builtin_predecrement);
    symtab_add_builtin("AddTo", builtin_addto);
    symtab_add_builtin("SubtractFrom", builtin_subtractfrom);
    symtab_get_def("Increment")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("Decrement")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("PreIncrement")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("PreDecrement")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("AddTo")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_get_def("SubtractFrom")->attributes |= ATTR_HOLDFIRST | ATTR_PROTECTED;
    symtab_add_builtin("TimeConstrained", builtin_time_constrained);
    symtab_get_def("TimeConstrained")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_add_builtin("Chop", builtin_chop);
    symtab_get_def("Chop")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Clip", builtin_clip);
    symtab_get_def("Clip")->attributes |= ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_add_builtin("Attributes", builtin_attributes);
    symtab_add_builtin("SetAttributes", builtin_set_attributes);
    symtab_add_builtin("OwnValues", builtin_own_values);
    symtab_add_builtin("DownValues", builtin_down_values);
    symtab_add_builtin("Out", builtin_out);
    symtab_add_builtin("Plus", builtin_plus);
    symtab_add_builtin("Times", builtin_times);
    symtab_add_builtin("Divide", builtin_divide);
    symtab_add_builtin("Subtract", builtin_subtract);
    symtab_add_builtin("Complex", builtin_complex);
    symtab_add_builtin("Rational", builtin_rational);
    symtab_add_builtin("Power", builtin_power);
    symtab_add_builtin("Sqrt", builtin_sqrt);
    symtab_add_builtin("Apply", builtin_apply);
    symtab_add_builtin("Map", builtin_map);
    symtab_add_builtin("MapIndexed", builtin_mapindexed);
    symtab_get_def("MapIndexed")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MapIndexed",
        "MapIndexed[f, list]\n\tGives {f[e1, {1}], f[e2, {2}], ...}. Over an\n"
        "\tassociation, f[value, {Key[k]}] keeping keys.");
    symtab_add_builtin("MapAll", builtin_map_all);
    symtab_add_builtin("MapAt", builtin_map_at);
    symtab_get_def("MapAt")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MapAt",
        "MapAt[f, expr, n]\n"
        "\tapplies f to the element at position n in expr. Negative n counts from the end.\n"
        "MapAt[f, expr, {i, j, ...}]\n"
        "\tapplies f to the part of expr at position {i, j, ...}.\n"
        "MapAt[f, expr, {{i1, j1, ...}, {i2, j2, ...}, ...}]\n"
        "\tapplies f to the parts of expr at each of the listed positions.\n"
        "\n"
        "Positions may contain All or Span specifications. MapAt[f, expr, 0]\n"
        "applies f to the head of expr. Repeated positions apply f repeatedly.");
    symtab_add_builtin("Nest", builtin_nest);
    symtab_get_def("Nest")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Nest",
        "Nest[f, expr, n]\n"
        "\tgives an expression with f applied n times to expr.\n"
        "\n"
        "n must be a non-negative integer. Nest[f, expr, 0] returns expr. The\n"
        "function f may be a symbol or a pure function. Each iteration evaluates\n"
        "f applied to the current value before proceeding.");
    symtab_add_builtin("NestList", builtin_nestlist);
    symtab_get_def("NestList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NestList",
        "NestList[f, expr, n]\n"
        "\tgives a list of the results of applying f to expr 0 through n times.\n"
        "\n"
        "The result is a list of length n+1 whose first element is expr and\n"
        "whose (k+1)-th element is f applied k times to expr. n must be a\n"
        "non-negative integer. f may be a symbol or a pure function; each\n"
        "intermediate application is evaluated before the next one.");
    symtab_add_builtin("NestWhile", builtin_nestwhile);
    symtab_get_def("NestWhile")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NestWhile",
        "NestWhile[f, expr, test]\n"
        "\tstarts with expr and repeatedly applies f while test still yields True.\n"
        "NestWhile[f, expr, test, m]\n"
        "\tsupplies the most recent m results as arguments to test.\n"
        "NestWhile[f, expr, test, All]\n"
        "\tsupplies all results so far as arguments to test.\n"
        "NestWhile[f, expr, test, {mmin, mmax}]\n"
        "\tdelays testing until at least mmin results exist, then passes up to mmax.\n"
        "NestWhile[f, expr, test, m, max]\n"
        "\tapplies f at most max times.\n"
        "NestWhile[f, expr, test, m, max, n]\n"
        "\tapplies f an additional n times after the loop terminates.\n"
        "NestWhile[f, expr, test, m, max, -n]\n"
        "\treturns the result found when f had been applied n fewer times.\n"
        "\n"
        "If test[expr] does not yield True initially, NestWhile returns expr.\n"
        "NestWhile[f, expr, UnsameQ, 2] is equivalent to FixedPoint[f, expr].");
    symtab_add_builtin("NestWhileList", builtin_nestwhilelist);
    symtab_get_def("NestWhileList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NestWhileList",
        "NestWhileList[f, expr, test]\n"
        "\tgenerates the list {expr, f[expr], f[f[expr]], ...} continuing while\n"
        "\ttest applied to the most recent result yields True.\n"
        "NestWhileList[f, expr, test, m]\n"
        "\tsupplies the most recent m results as arguments to test.\n"
        "NestWhileList[f, expr, test, All]\n"
        "\tsupplies all results so far as arguments to test.\n"
        "NestWhileList[f, expr, test, {mmin, mmax}]\n"
        "\tdelays testing until at least mmin results exist, then passes up to mmax.\n"
        "NestWhileList[f, expr, test, m, max]\n"
        "\tapplies f at most max times.\n"
        "NestWhileList[f, expr, test, m, max, n]\n"
        "\tappends n additional applications of f to the list.\n"
        "NestWhileList[f, expr, test, m, max, -n]\n"
        "\tdrops the last n elements from the list.\n"
        "\n"
        "NestWhileList[f, expr, UnsameQ, 2] is equivalent to FixedPointList[f, expr].\n"
        "NestWhileList[f, expr, test, All] is equivalent to\n"
        "NestWhileList[f, expr, test, {1, Infinity}].");
    symtab_add_builtin("FixedPointList", builtin_fixedpointlist);
    symtab_get_def("FixedPointList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FixedPointList",
        "FixedPointList[f, expr]\n"
        "\tgenerates the list {expr, f[expr], f[f[expr]], ...} of successive\n"
        "\tapplications of f, stopping when two consecutive results are SameQ.\n"
        "\tThe last two elements of the result are always the same.\n"
        "FixedPointList[f, expr, n]\n"
        "\tstops after at most n applications of f. If n is reached before\n"
        "\tconvergence, the last two elements may not be equal.\n"
        "FixedPointList[f, expr, SameTest -> s]\n"
        "FixedPointList[f, expr, n, SameTest -> s]\n"
        "\tuses the binary predicate s instead of SameQ to test successive pairs.\n"
        "\n"
        "FixedPointList[f, expr] is equivalent to\n"
        "NestWhileList[f, expr, UnsameQ, 2].");
    symtab_add_builtin("FixedPoint", builtin_fixedpoint);
    symtab_get_def("FixedPoint")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FixedPoint",
        "FixedPoint[f, expr]\n"
        "\tstarts with expr and applies f repeatedly until the result no longer\n"
        "\tchanges, returning the final value.\n"
        "FixedPoint[f, expr, n]\n"
        "\tstops after at most n applications of f, returning the last value\n"
        "\tobtained even if a fixed point has not been reached.\n"
        "FixedPoint[f, expr, SameTest -> s]\n"
        "FixedPoint[f, expr, n, SameTest -> s]\n"
        "\tuses the binary predicate s instead of SameQ to test successive pairs.\n"
        "\n"
        "FixedPoint[f, expr] gives the last element of FixedPointList[f, expr].\n"
        "Throw can be used inside f to exit early.");
    symtab_add_builtin("Fold", builtin_fold);
    symtab_get_def("Fold")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Fold",
        "Fold[f, x, list]\n"
        "\tgives the last element of FoldList[f, x, list]:\n"
        "\tf[...f[f[f[x, list[[1]]], list[[2]]], list[[3]]]..., list[[n]]].\n"
        "Fold[f, list]\n"
        "\tis equivalent to Fold[f, First[list], Rest[list]].\n"
        "\n"
        "The head of list need not be List. Fold[f, x, {}] returns x, and\n"
        "Fold[f, {a}] returns a. Fold[f, {}] remains unevaluated. f may be a\n"
        "symbol or a pure function; each intermediate application is evaluated\n"
        "before the next one.");
    symtab_add_builtin("FoldList", builtin_foldlist);
    symtab_get_def("FoldList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FoldList",
        "FoldList[f, x, list]\n"
        "\tgives {x, f[x, list[[1]]], f[f[x, list[[1]]], list[[2]]], ...}.\n"
        "FoldList[f, list]\n"
        "\tgives {list[[1]], f[list[[1]], list[[2]]], ...}.\n"
        "\n"
        "For a length-n list, FoldList generates a list of length n+1. The\n"
        "head of list is preserved in the output:\n"
        "\tFoldList[f, x, p[a, b]] -> p[x, f[x, a], f[f[x, a], b]].\n"
        "FoldList[f, {}] returns an empty list {}. f may be a symbol or a\n"
        "pure function; each intermediate application is evaluated before the\n"
        "next one.");
    symtab_add_builtin("Through", builtin_through);
    symtab_add_builtin("Thread", builtin_thread);
    symtab_get_def("Thread")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Distribute", builtin_distribute);
    symtab_get_def("Distribute")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Inner", builtin_inner);
    symtab_get_def("Inner")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Outer", builtin_outer);
    symtab_get_def("Outer")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Tuples", builtin_tuples);
    symtab_get_def("Tuples")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Permutations", builtin_permutations);
    symtab_get_def("Permutations")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Select", builtin_select);
    symtab_add_builtin("TakeWhile", builtin_takewhile);
    symtab_get_def("TakeWhile")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TakeWhile",
        "TakeWhile[list, crit]\n\tGives the longest leading run of elements e for\n"
        "\twhich crit[e] is True. Over an association, tests the values and keeps\n"
        "\tthe matching leading entries (keys preserved).");
    symtab_add_builtin("LengthWhile", builtin_lengthwhile);
    symtab_get_def("LengthWhile")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LengthWhile",
        "LengthWhile[list, crit]\n\tGives the length of the longest leading run of\n"
        "\telements e for which crit[e] is True. Over an association, tests values.");
    symtab_add_builtin("Scan", builtin_scan);
    symtab_get_def("Scan")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Scan",
        "Scan[f, expr]\n\tApplies f to each element of expr for its side effects\n"
        "\tand returns Null. Over an association, applies f to each value.");
    symtab_add_builtin("SelectFirst", builtin_select_first);
    symtab_get_def("SelectFirst")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SelectFirst",
        "SelectFirst[list, pred]\n\tGives the first element e of list for which\n"
        "\tpred[e] is True, or Missing[\"NotFound\"]. SelectFirst[list, pred, default]\n"
        "\tuses default. Over an association, tests values and returns the first match.");
    symtab_add_builtin("AllTrue", builtin_all_true);
    symtab_get_def("AllTrue")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AllTrue",
        "AllTrue[list, test]\n\tGives True if test[e] is True for every element e\n"
        "\t(True for an empty list). Over an association, tests the values.");
    symtab_add_builtin("AnyTrue", builtin_any_true);
    symtab_get_def("AnyTrue")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("AnyTrue",
        "AnyTrue[list, test]\n\tGives True if test[e] is True for some element e\n"
        "\t(False for an empty list). Over an association, tests the values.");
    symtab_add_builtin("NoneTrue", builtin_none_true);
    symtab_get_def("NoneTrue")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NoneTrue",
        "NoneTrue[list, test]\n\tGives True if test[e] is True for no element e\n"
        "\t(True for an empty list). Over an association, tests the values.");
    symtab_add_builtin("FreeQ", builtin_freeq);
    symtab_add_builtin("Sort", builtin_sort);
    symtab_add_builtin("SortBy", builtin_sort_by);
    symtab_get_def("SortBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SortBy",
        "SortBy[list, f]\n\tSorts the elements of list by the canonical order of\n"
        "\tf applied to each element.\n"
        "SortBy[assoc, f]\n\tSorts an association by f applied to each value.\n"
        "SortBy[f]\n\tOperator form: SortBy[f][expr] is SortBy[expr, f].");
    symtab_add_builtin("MaximalBy", builtin_maximal_by);
    symtab_get_def("MaximalBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MaximalBy",
        "MaximalBy[list, f]\n\tGives the element(s) of list for which f is maximal\n"
        "\t(all ties, in order). Over an association, gives the entries whose\n"
        "\tvalue maximises f. MaximalBy[f] is the operator form.");
    symtab_add_builtin("MinimalBy", builtin_minimal_by);
    symtab_get_def("MinimalBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MinimalBy",
        "MinimalBy[list, f]\n\tGives the element(s) of list for which f is minimal\n"
        "\t(all ties, in order). Over an association, gives the entries whose\n"
        "\tvalue minimises f. MinimalBy[f] is the operator form.");
    symtab_add_builtin("TakeLargest", builtin_take_largest);
    symtab_get_def("TakeLargest")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TakeLargest",
        "TakeLargest[list, n]\n\tGives the n largest elements of list, in\n"
        "\tdescending order. Over an association, gives the n entries with the\n"
        "\tlargest values (as an association).");
    symtab_add_builtin("TakeSmallest", builtin_take_smallest);
    symtab_get_def("TakeSmallest")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TakeSmallest",
        "TakeSmallest[list, n]\n\tGives the n smallest elements of list, in\n"
        "\tascending order. Over an association, gives the n entries with the\n"
        "\tsmallest values (as an association).");
    symtab_add_builtin("TakeLargestBy", builtin_take_largest_by);
    symtab_get_def("TakeLargestBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TakeLargestBy",
        "TakeLargestBy[list, f, n]\n\tGives the n elements of list for which f is\n"
        "\tlargest, in descending order of f. Over an association, ranks by f of\n"
        "\teach value.");
    symtab_add_builtin("TakeSmallestBy", builtin_take_smallest_by);
    symtab_get_def("TakeSmallestBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("TakeSmallestBy",
        "TakeSmallestBy[list, f, n]\n\tGives the n elements of list for which f is\n"
        "\tsmallest, in ascending order of f. Over an association, ranks by f of\n"
        "\teach value.");
    symtab_add_builtin("ReverseSort", builtin_reverse_sort);
    symtab_get_def("ReverseSort")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReverseSort",
        "ReverseSort[list]\n\tSorts into descending order (Reverse of Sort).\n"
        "\tOver an association, sorts the entries by value, descending.");
    symtab_add_builtin("ReverseSortBy", builtin_reverse_sort_by);
    symtab_get_def("ReverseSortBy")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ReverseSortBy",
        "ReverseSortBy[list, f]\n\tSorts by f in descending order. Over an\n"
        "\tassociation, sorts by f of each value, descending.");
    symtab_add_builtin("OrderedQ", builtin_orderedq);
    symtab_get_def("OrderedQ")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialQ", builtin_polynomialq);
    symtab_add_builtin("Variables", builtin_variables);
    symtab_add_builtin("Level", builtin_level);
    symtab_add_builtin("Depth", builtin_depth);
    symtab_add_builtin("LeafCount", builtin_leafcount);
    symtab_get_def("LeafCount")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ByteCount", builtin_bytecount);
    symtab_get_def("ByteCount")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatchQ", builtin_matchq);
    symtab_add_builtin("CompoundExpression", builtin_compoundexpression);
    symtab_add_builtin("NumberQ", builtin_numberq);
    symtab_add_builtin("NumericQ", builtin_numericq);
    symtab_add_builtin("Positive", builtin_positive);
    symtab_add_builtin("Negative", builtin_negative);
    symtab_add_builtin("NonNegative", builtin_nonnegative);
    symtab_add_builtin("NonPositive", builtin_nonpositive);
    symtab_add_builtin("IntegerQ", builtin_integerq);
    symtab_add_builtin("EvenQ", builtin_evenq);
    symtab_add_builtin("OddQ", builtin_oddq);
    symtab_add_builtin("Mod", builtin_mod);
    symtab_add_builtin("Quotient", builtin_quotient);
    symtab_add_builtin("QuotientRemainder", builtin_quotientremainder);
    /* Number-theory builtins (GCD, LCM, ExtendedGCD, PowerMod, Factorial,
       Factorial2, FactorialPower, Binomial, PrimitiveRoot, PrimitiveRootList,
       MultiplicativeOrder) are registered by numbertheory_init() below. */
    symtab_add_builtin("Print", builtin_print);
    symtab_add_builtin("FullForm", builtin_fullform);
    symtab_add_builtin("InputForm", builtin_inputform);
    symtab_add_builtin("TeXForm", builtin_texform);
    symtab_add_builtin("Information", builtin_information);
    symtab_add_builtin("Evaluate", builtin_evaluate);
    symtab_get_def("Evaluate")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReleaseHold", builtin_releasehold);
    symtab_get_def("ReleaseHold")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ToString", builtin_tostring);
    symtab_get_def("ToString")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ToExpression", builtin_toexpression);
    symtab_get_def("ToExpression")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Symbol", builtin_symbol);
    symtab_get_def("Symbol")->attributes |= ATTR_PROTECTED;

    symtab_get_def("AtomQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NumberQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NumericQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Positive")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_get_def("Negative")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_get_def("NonNegative")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_get_def("NonPositive")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
    symtab_get_def("IntegerQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("EvenQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("OddQ")->attributes |= ATTR_PROTECTED;
    /* HoldAll (matching real Mathematica): ?x/Information[x] must inspect
     * the symbol x itself, not its evaluated value -- otherwise a symbol
     * with an eager OwnValue (e.g. the named color constants in
     * src/graphics/graphics_init.c) would have its value, not its name,
     * reach builtin_information, which only recognizes a bare symbol. */
    symtab_get_def("Information")->attributes |= (ATTR_HOLDALL | ATTR_PROTECTED);

    symtab_get_def("Mod")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Quotient")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("QuotientRemainder")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    /* Number-theory builtin attributes + docstrings live in numbertheory_init(). */
    symtab_get_def("Print")->attributes |= ATTR_PROTECTED;
    symtab_get_def("FullForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("InputForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("TeXForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("HoldForm")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    facint_init();
    numbertheory_init();
    gamma_init();
    beta_init();
    erf_init();
    erfc_init();
    erfi_init();
    expintegralei_init();
    logintegral_init();
    sinintegral_init();
    cosintegral_init();
    sinhintegral_init();
    coshintegral_init();
    fresnel_init();
    sinc_init();
    inverf_init();
    inverfc_init();
    loggamma_init();
    polygamma_init();
    harmonicnumber_init();
    pochhammer_init();
    eulergamma_init();
    zeta_init();
    hurwitzzeta_init();
    stieltjesgamma_init();
    bernoullib_init();
    eulere_init();
    polylog_init();
    lerchphi_init();
    productlog_init();
    legendre_init();
    airyai_init();
    airybi_init();
    bessel_init();
    void besseljzero_init(void);    besseljzero_init();
    fibonacci_init();
    lucas_init();
    void hyperfactorial_init(void); hyperfactorial_init();
    void barnesg_init(void);        barnesg_init();
    void qpochhammer_init(void);    qpochhammer_init();

    Expr* zero = expr_new_integer(0);
    Expr* one = expr_new_integer(1);
    Expr* sym_I = expr_new_symbol(SYM_I);
    Expr* val_I = make_complex(zero, one);
    symtab_add_own_value("I", sym_I, val_I);
    expr_free(sym_I);
    expr_free(val_I);
    
    comparisons_init();
    boolean_init();
    names_init();
    list_init();
    assoc_init();
    ndarray_init();
    replace_init();
    patterns_init();
    cond_init();
    iter_init();
    complex_init();
    trig_init();
    hyperbolic_init();
    void trigsimp_init(void);
    trigsimp_init();
    void simp_init(void);
    simp_init();
    logexp_init();
    piecewise_init();
    void interp_init(void);
    interp_init();
    int_init();
    real_init();
    attr_init();
    purefunc_init();
    stats_init();
    partitions_init();
    poly_init();
    facpoly_init();
    flint_bridge_init();
    void flint_mat_bridge_init(void);
    flint_mat_bridge_init();
    void flint_num_bridge_init(void);
    flint_num_bridge_init();
    void blas_bridge_init(void);
    blas_bridge_init();
    void lapack_bridge_init(void);
    lapack_bridge_init();
    void squarefreeq_init(void);
    squarefreeq_init();
    void irrpolyq_init(void);
    irrpolyq_init();
    void minpoly_init(void);
    minpoly_init();
    rat_init();
    void rootreduce_init(void);
    rootreduce_init();
    expand_init();
    expand_power_init();
    solve_init();
    findroot_init();
    void findmin_init(void);
    findmin_init();
    void nresidue_init(void);
    nresidue_init();
    void nd_init(void);
    nd_init();
    void nseries_init(void);
    nseries_init();
    void nlimit_init(void);
    nlimit_init();
    void nsum_init(void);
    nsum_init();
    void nprod_init(void);
    nprod_init();
    void nintegrate_init(void);
    nintegrate_init();
    void nroots_init(void);
    nroots_init();
    void nsolve_init(void);
    nsolve_init();
    info_init();
    datetime_init();
    linalg_init();
    void matsol_init(void);
    matsol_init();
    void matinv_init(void);
    matinv_init();
    void matlstsq_init(void);
    matlstsq_init();
    void mateigen_init(void);
    mateigen_init();
    void matnull_init(void);
    matnull_init();
    void matrank_init(void);
    matrank_init();
    void qrdecomp_init(void);
    qrdecomp_init();
    void ludecomp_init(void);
    ludecomp_init();
    void svdecomp_init(void);
    svdecomp_init();
    fit_init();
    readwrite_init();
    loadmodule_init();
    files_init();
    random_init();
    strings_init();
    regex_init();
    series_init();
    deriv_init();
    limit_init();
    void residue_init(void);
    residue_init();
    numeric_init();
    precision_init();
    rationalize_init();
    void root_init(void);
    root_init();
    void radicals_init(void);
    radicals_init();
    void groebner_init(void);
    groebner_init();
    void eliminate_init(void);
    eliminate_init();
    void solvealways_init(void);
    solvealways_init();
    void integrate_init(void);
    integrate_init();
    void sum_init(void);
    sum_init();
    void product_init(void);
    product_init();
    void hypergeopfq_init(void);
    hypergeopfq_init();
    void zero_test_init(void);
    zero_test_init();
    void graphics_init(void);
    graphics_init();
    void graph_init(void);
    graph_init();

    /* Options/SetOptions/OptionValue + the default-options registry. Runs last
     * so every option-name symbol used by the registry is already interned. */
    void options_builtin_init(void);
    options_builtin_init();

    /* Flag every symbol interned so far as a System symbol. At this point in
     * startup the interner holds exactly the kernel's built-in names (cached
     * SYM_* pointers, registered builtins, option names with docstrings, ...).
     * The context resolver consults this flag so that a bare System name --
     * including pure structural heads and constants (List, Rule, Integer,
     * Automatic, Infinity, None, Heads, ...) that have no SymbolDef until used
     * -- resolves to System` rather than being qualified into the current
     * context. Without it, a package referencing such a name in its prologue
     * (e.g. a `_List` pattern head or an Automatic option default) would bind a
     * stray private symbol. This only sets a flag; it creates no SymbolDefs and
     * changes no evaluation behaviour. Must run before the first package load. */
    intern_mark_all_system();

    /* Eager-load the FullSimplify package (driver + manifest only; the heavy
     * per-function identity libraries load lazily on demand). Done last so
     * every builtin it relies on -- Simplify above all -- is already
     * registered. Resolved working-directory-independently so the REPL and the
     * CMake test binaries both pick it up.
     *
     * FullSimplify.m is a proper package: BeginPackage["FullSimplify`"]
     * exports FullSimplify and RegisterTransforms (with their usage messages,
     * which become the symbols' docstrings), implements everything in the
     * FullSimplify`Private` context, and runs Protect[FullSimplify] itself.
     * So no C-side docstring or attribute wiring is needed here -- the package
     * is self-documenting and self-protecting, and EndPackage[] leaves
     * FullSimplify` on $ContextPath so bare FullSimplify[...] resolves. */
    mathilda_load_module("simp/FullSimplify.m");
}

Expr* builtin_compoundexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 0) return expr_new_symbol(SYM_Null);

    Expr* last_val = NULL;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (last_val) expr_free(last_val);
        last_val = evaluate(res->data.function.args[i]);
        if (last_val->type == EXPR_FUNCTION && last_val->data.function.head->type == EXPR_SYMBOL) {
            const char* hname = last_val->data.function.head->data.symbol.name;
            if (hname == SYM_Return || hname == SYM_Break || 
                hname == SYM_Continue || hname == SYM_Throw || 
                hname == SYM_Abort || hname == SYM_Quit) {
                break;
            }
        }
    }
    return last_val;
}
Expr* builtin_clear(Expr* res) {
    if (res->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < res->data.function.arg_count; i++) {
            Expr* arg = res->data.function.args[i];
            if (arg->type == EXPR_SYMBOL) {
                symtab_clear_symbol(arg->data.symbol.name);
            }
        }
        return expr_new_symbol(SYM_Null);
    }
    return NULL;
}

/* Unset[lhs] / `lhs =.`: remove the single rule whose left-hand side is
 * `lhs` (up to renaming of bound pattern variables). A bare symbol clears
 * its OwnValue; a function form clears the matching DownValue on the head
 * symbol. Unset carries HoldFirst, so `lhs` arrives unevaluated. Always
 * returns Null, matching Mathematica -- whether or not a rule was found. */
Expr* builtin_unset(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }
    Expr* lhs = res->data.function.args[0];

    const char* symbol_name = NULL;
    bool own_value = false;

    if (lhs->type == EXPR_SYMBOL) {
        symbol_name = lhs->data.symbol.name;
        own_value = true;
    } else if (lhs->type == EXPR_FUNCTION &&
               lhs->data.function.head &&
               lhs->data.function.head->type == EXPR_SYMBOL) {
        const char* head = lhs->data.function.head->data.symbol.name;
        /* `f[x_] /; cond =.` parses to Unset[Condition[f[x_], cond]]; the
         * stored rule is keyed by the inner head, so resolve through it. */
        if (head == SYM_Condition && lhs->data.function.arg_count == 2) {
            Expr* inner = lhs->data.function.args[0];
            if (inner->type == EXPR_SYMBOL) {
                symbol_name = inner->data.symbol.name;
                own_value = true;
            } else if (inner->type == EXPR_FUNCTION &&
                       inner->data.function.head &&
                       inner->data.function.head->type == EXPR_SYMBOL) {
                symbol_name = inner->data.function.head->data.symbol.name;
            } else {
                return NULL;
            }
        } else {
            symbol_name = head;
        }
    } else {
        /* e.g. Unset[5] -- not an assignable left-hand side. */
        return NULL;
    }

    /* Protected/Locked symbols (every builtin among them) cannot be Unset,
     * mirroring Set's wrsym guard. */
    uint32_t attrs = get_attributes(symbol_name);
    if (attrs & (ATTR_PROTECTED | ATTR_LOCKED)) {
        fprintf(stderr, "Unset::wrsym: Symbol %s is Protected.\n", symbol_name);
        return expr_new_symbol(SYM_Null);
    }

    symtab_remove_matching_rule(symbol_name, lhs, own_value);
    return expr_new_symbol(SYM_Null);
}

/* ============================================================
 * ClearAll / Remove / Protect / Unprotect
 *
 * The four symbol-management builtins share an argument shape: each
 * argument is a symbol, a string naming a symbol, or a List of such
 * specs (e.g. ClearAll[{a, b}]). All four carry HoldAll so the symbols
 * arrive at the builtin unevaluated.
 * ============================================================ */

/* Extract a symbol name from a spec element, or NULL if it is neither a
 * bare symbol nor a string. */
static const char* core_symbol_name_of(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name;
    if (e->type == EXPR_STRING) return e->data.string;
    return NULL;
}

/* ClearAll[s]: drop all OwnValues/DownValues, attributes and the usage
 * message (docstring) for s -- unless s is Protected or Locked, which
 * ClearAll never touches. The builtin C function pointer (if any) is
 * left intact; the Protected guard already shields every builtin. */
static void core_clear_all_one(const char* name) {
    if (!name) return;
    uint32_t attrs = get_attributes(name);
    if (attrs & (ATTR_PROTECTED | ATTR_LOCKED)) return;

    symtab_clear_symbol(name);          /* values */
    SymbolDef* def = symtab_get_def(name);
    if (def->attributes != 0) {         /* attributes */
        def->attributes = 0;
        eval_clock_bump();
    }
    if (def->docstring) {               /* usage / messages */
        free(def->docstring);
        def->docstring = NULL;
    }
}

/* Remove[s]: delete the symbol's definition entirely. The Protected /
 * Locked guard is what keeps Remove from ever deleting a builtin. */
static void core_remove_one(const char* name) {
    if (!name) return;
    uint32_t attrs = get_attributes(name);
    if (attrs & (ATTR_PROTECTED | ATTR_LOCKED)) return;
    symtab_remove_symbol(name);
}

/* Protect[s]: set the Protected attribute. Returns true iff the bit was
 * newly set (so the caller can report the changed name, as WL does).
 * Locked symbols are left untouched. */
static bool core_protect_one(const char* name) {
    if (!name) return false;
    SymbolDef* def = symtab_get_def(name);
    if (def->attributes & ATTR_LOCKED) return false;
    if (def->attributes & ATTR_PROTECTED) return false;
    def->attributes |= ATTR_PROTECTED;
    eval_clock_bump();
    return true;
}

/* Unprotect[s]: clear the Protected attribute. Returns true iff the bit
 * was actually cleared. Locked symbols are left untouched. */
static bool core_unprotect_one(const char* name) {
    if (!name) return false;
    SymbolDef* def = symtab_get_def(name);
    if (def->attributes & ATTR_LOCKED) return false;
    if (!(def->attributes & ATTR_PROTECTED)) return false;
    def->attributes &= ~ATTR_PROTECTED;
    eval_clock_bump();
    return true;
}

/* Apply `action` to every symbol spec in `res`'s argument list, where a
 * spec is a symbol/string or a flat List of them. Used by ClearAll and
 * Remove, which both return Null. */
static Expr* core_apply_symbol_action(Expr* res, void (*action)(const char*)) {
    if (res->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (head_is(arg, SYM_List)) {
            for (size_t j = 0; j < arg->data.function.arg_count; j++) {
                action(core_symbol_name_of(arg->data.function.args[j]));
            }
        } else {
            action(core_symbol_name_of(arg));
        }
    }
    return expr_new_symbol(SYM_Null);
}

Expr* builtin_clear_all(Expr* res) {
    return core_apply_symbol_action(res, core_clear_all_one);
}

Expr* builtin_remove(Expr* res) {
    return core_apply_symbol_action(res, core_remove_one);
}

/* Shared driver for Protect / Unprotect. Walks the spec arguments,
 * applies the change, and returns a List of the names (as strings) whose
 * Protected state actually changed -- matching WL's return value. */
static Expr* core_protect_unprotect(Expr* res, bool protecting) {
    if (res->type != EXPR_FUNCTION) return NULL;

    Expr** changed = NULL;
    size_t count = 0, cap = 0;

    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        size_t n = head_is(arg, SYM_List) ? arg->data.function.arg_count : 1;
        for (size_t j = 0; j < n; j++) {
            const char* name = core_symbol_name_of(
                head_is(arg, SYM_List) ? arg->data.function.args[j] : arg);
            bool did = protecting ? core_protect_one(name)
                                  : core_unprotect_one(name);
            if (did) {
                if (count == cap) {
                    cap = cap ? cap * 2 : 4;
                    changed = realloc(changed, cap * sizeof(Expr*));
                }
                changed[count++] = expr_new_string(name);
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), changed, count);
    free(changed);
    return list;
}

Expr* builtin_protect(Expr* res) {
    return core_protect_unprotect(res, true);
}

Expr* builtin_unprotect(Expr* res) {
    return core_protect_unprotect(res, false);
}

Expr* builtin_length(Expr* res) {
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];
        int64_t len = 0;
        if (arg->type == EXPR_FUNCTION) {
            len = (int64_t)arg->data.function.arg_count;
        } else if (arg->type == EXPR_NDARRAY) {
            /* numpy len(): the leading-axis length (rank >= 1 always). */
            len = arg->data.ndarray.dims[0];
        }
        return expr_new_integer(len);
    }
    return NULL;
}

#define DIMENSIONS_MAX_DEPTH 64

/*
 * get_dimensions:
 * Recursively determines the dimensions of a rectangular nested structure
 * whose every level shares the head named by head_name. Stops as soon as
 * the structure becomes ragged (sub-arrays differ in shape) or once
 * max_depth levels have been recorded. Returns the number of dimensions
 * written into dims (which must hold at least max_depth slots).
 */
static int get_dimensions(Expr* e, int64_t* dims, int max_depth, const char* head_name) {
    if (max_depth <= 0) return 0;
    if (!head_is(e, intern_symbol(head_name))) {
        return 0;
    }

    dims[0] = (int64_t)e->data.function.arg_count;
    if (dims[0] == 0 || max_depth == 1) return 1;

    int64_t sub_dims[DIMENSIONS_MAX_DEPTH];
    int sub_depth = get_dimensions(e->data.function.args[0], sub_dims, max_depth - 1, head_name);

    if (sub_depth == 0) return 1;

    for (size_t i = 1; i < e->data.function.arg_count; i++) {
        int64_t cur_dims[DIMENSIONS_MAX_DEPTH];
        int cur_depth = get_dimensions(e->data.function.args[i], cur_dims, max_depth - 1, head_name);

        if (cur_depth != sub_depth) return 1;
        for (int j = 0; j < sub_depth; j++) {
            if (cur_dims[j] != sub_dims[j]) return 1;
        }
    }

    for (int i = 0; i < sub_depth; i++) {
        dims[i + 1] = sub_dims[i];
    }
    return sub_depth + 1;
}

Expr* builtin_dimensions(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    /* Determine the level cap from the optional second argument. */
    int max_depth = DIMENSIONS_MAX_DEPTH;
    if (argc == 2) {
        Expr* n = res->data.function.args[1];
        int64_t nv;
        if (n->type == EXPR_INTEGER) {
            nv = n->data.integer;
        } else if (head_is(n, SYM_DirectedInfinity) &&
                   n->data.function.arg_count == 1 &&
                   n->data.function.args[0]->type == EXPR_INTEGER &&
                   n->data.function.args[0]->data.integer == 1) {
            /* Dimensions[expr, Infinity] -> use the full depth cap. */
            nv = DIMENSIONS_MAX_DEPTH;
        } else {
            return NULL;
        }
        if (nv < 0) return NULL;
        if (nv > DIMENSIONS_MAX_DEPTH) nv = DIMENSIONS_MAX_DEPTH;
        max_depth = (int)nv;
    }

    Expr* arg = res->data.function.args[0];
    int depth = 0;
    int64_t dims[DIMENSIONS_MAX_DEPTH];
    if (arg->type == EXPR_NDARRAY) {
        /* O(1): rank/dims are already stored directly, no probing needed. */
        depth = arg->data.ndarray.rank;
        if (depth > max_depth) depth = max_depth;
        for (int i = 0; i < depth; i++) dims[i] = arg->data.ndarray.dims[i];
    } else if (max_depth > 0 && arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL) {
        depth = get_dimensions(arg, dims, max_depth, arg->data.function.head->data.symbol.name);
    }

    Expr** dim_args = (depth > 0) ? malloc(sizeof(Expr*) * (size_t)depth) : NULL;
    for (int i = 0; i < depth; i++) {
        dim_args[i] = expr_new_integer(dims[i]);
    }

    Expr* ret = expr_new_function(expr_new_symbol(SYM_List), dim_args, (size_t)depth);
    free(dim_args);
    return ret;
}

/*
 * is_hold_head:
 * Returns true if the symbol name is one of the standard hold wrappers
 * that ReleaseHold should strip.
 */
static bool is_hold_head(const char* name) {
    return (strcmp(name, "Hold") == 0 ||
            strcmp(name, "HoldForm") == 0 ||
            strcmp(name, "HoldPattern") == 0 ||
            strcmp(name, "HoldComplete") == 0);
}

/*
 * release_hold_recursive:
 * Traverses the expression tree and replaces any Hold/HoldForm/HoldPattern/
 * HoldComplete wrapper with its contents. Does NOT recurse into the
 * replacement contents (only removes one layer).
 * Returns a new expression (caller owns it).
 */
static Expr* release_hold_recursive(Expr* e) {
    if (!e) return NULL;

    /* Check if e itself is a hold wrapper */
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        is_hold_head(e->data.function.head->data.symbol.name)) {
        /* Strip the wrapper: for single arg, return the arg as-is (copy).
         * For multiple args, wrap in Sequence. */
        if (e->data.function.arg_count == 1) {
            return expr_copy(e->data.function.args[0]);
        } else {
            /* Multiple args: wrap in Sequence */
            Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                args[i] = expr_copy(e->data.function.args[i]);
            }
            Expr* seq = expr_new_function(expr_new_symbol(SYM_Sequence), args, e->data.function.arg_count);
            free(args);
            return seq;
        }
    }

    /* Not a hold wrapper: recurse into arguments (but not head) */
    if (e->type == EXPR_FUNCTION) {
        bool changed = false;
        Expr** new_args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            new_args[i] = release_hold_recursive(e->data.function.args[i]);
            if (!expr_eq(new_args[i], e->data.function.args[i])) {
                changed = true;
            }
        }
        if (changed) {
            Expr* result = expr_new_function(expr_copy(e->data.function.head), new_args, e->data.function.arg_count);
            free(new_args);
            return result;
        }
        /* No changes: free copies and return original copy */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            expr_free(new_args[i]);
        }
        free(new_args);
        return expr_copy(e);
    }

    /* Atomic expression: return copy */
    return expr_copy(e);
}

/*
 * Chop: replace approximate-real numbers smaller in magnitude than delta
 * by the exact integer 0.
 *
 * Only EXPR_REAL (and EXPR_MPFR when USE_MPFR is enabled) values are
 * candidates for chopping; integers, bigints, rationals, and symbolic
 * subexpressions pass through untouched.
 *
 * "Machine complex" Complex[re, im] -- where both components are
 * machine reals -- gets the Mathematica special case: chopping the
 * real part yields the machine zero 0.0 (preserving the machine-complex
 * shape) rather than the exact integer 0, but chopping the imaginary
 * part drops the whole Complex wrapper so the result is just the
 * machine real `re`. Non-machine Complex[re, im] (e.g. Complex[1, 1.e-12])
 * is recursed into so the outer evaluator's Complex[r, 0]->r simplification
 * still fires.
 */
static bool chop_is_machine_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static double chop_to_double(const Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    return 0.0;
}

static Expr* chop_recursive(Expr* e, double delta) {
    if (!e) return NULL;

    if (chop_is_machine_real(e)) {
        if (fabs(chop_to_double(e)) < delta) {
            return expr_new_integer(0);
        }
        return expr_copy(e);
    }

    Expr *re = NULL, *im = NULL;
    if (is_complex(e, &re, &im)) {
        if (chop_is_machine_real(re) && chop_is_machine_real(im)) {
            double rv = chop_to_double(re);
            double iv = chop_to_double(im);
            bool re_small = fabs(rv) < delta;
            bool im_small = fabs(iv) < delta;

            if (re_small && im_small) {
                return expr_new_integer(0);
            }
            if (im_small) {
                /* Drop the small imaginary part: machine real survives. */
                return expr_copy(re);
            }
            if (re_small) {
                /* Preserve the machine-complex shape with a machine 0.0
                 * for the real part. Constructing Complex[0., im]
                 * directly avoids any re-evaluation that would collapse
                 * the wrapper. */
                Expr* zero = expr_new_real(0.0);
                Expr* im_copy = expr_copy(im);
                return make_complex(zero, im_copy);
            }
            /* Neither component chops: return the Complex unchanged. */
            return expr_copy(e);
        }
        /* Non-machine Complex: recurse into each part and let
         * builtin_complex auto-simplify Complex[r, 0] -> r on the next
         * evaluator pass. */
        Expr* new_re = chop_recursive(re, delta);
        Expr* new_im = chop_recursive(im, delta);
        return make_complex(new_re, new_im);
    }

    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        Expr** new_args = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = chop_recursive(e->data.function.args[i], delta);
        }
        Expr* result = expr_new_function(expr_copy(e->data.function.head), new_args, n);
        if (new_args) free(new_args);
        return result;
    }

    /* Atom (Integer, BigInt, Symbol, String, Rational components): no change. */
    return expr_copy(e);
}

/* Convert a delta argument to a positive double tolerance. Accepts
 * Integer, BigInt, Real, MPFR (when enabled), and Rational[n, d].
 * Returns -1.0 if the argument cannot be coerced. */
static double chop_extract_delta(Expr* d) {
    if (!d) return -1.0;
    if (d->type == EXPR_REAL) return fabs(d->data.real);
    if (d->type == EXPR_INTEGER) return fabs((double)d->data.integer);
    if (d->type == EXPR_BIGINT) return fabs(mpz_get_d(d->data.bigint));
#ifdef USE_MPFR
    if (d->type == EXPR_MPFR) return fabs(mpfr_get_d(d->data.mpfr, MPFR_RNDN));
#endif
    int64_t n, den;
    if (is_rational(d, &n, &den) && den != 0) {
        return fabs((double)n / (double)den);
    }
    return -1.0;
}

/*
 * builtin_chop:
 * Chop[expr]        replaces approximate real numbers in expr that are
 *                   close to zero by the exact integer 0 (default
 *                   tolerance 10^-10).
 * Chop[expr, delta] uses |delta| as the threshold.
 */
Expr* builtin_chop(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;

    double delta = 1.0e-10;
    if (argc == 2) {
        double d = chop_extract_delta(res->data.function.args[1]);
        if (d < 0.0) return NULL;
        delta = d;
    }

    return chop_recursive(res->data.function.args[0], delta);
}

/*
 * Clip: clamp a numeric value to an interval, optionally substituting
 * named replacement values when the input falls outside.
 *
 * Forms supported:
 *   Clip[x]                              -- clamp to [-1, 1]
 *   Clip[x, {min, max}]                  -- clamp to [min, max]
 *   Clip[x, {min, max}, {v_min, v_max}]  -- return v_min when x < min,
 *                                            v_max when x > max,
 *                                            x otherwise
 *
 * The interval endpoints and replacements may be symbolic; Clip only
 * needs to decide which side of `min` / `max` `x` is on, not the
 * algebraic identity of the bounds themselves. The decision is made by
 * numericalizing `x`, `min`, and `max` to machine doubles via the
 * existing `numericalize` machinery. Symbolic constants such as `Pi`
 * are therefore handled (Pi -> ~3.14 -> Clip[Pi] = 1), and the
 * original symbolic `x` is returned when it lies inside the interval
 * (rather than the numeric approximation).
 *
 * Special cases:
 *   - Infinity / -Infinity are handled before numericalization.
 *   - A complex (non-real) input emits Clip::ncompl and the call stays
 *     unevaluated, matching Mathematica.
 *   - When `x`, `min`, or `max` cannot be reduced to a number, the
 *     call stays unevaluated so user-supplied DownValues / pattern
 *     rules can still take over.
 */
static bool clip_to_double_value(Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        return true;
    }
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d) && d != 0) {
        *out = (double)n / (double)d;
        return true;
    }
    /* Symbolic forms (Pi, E, Sqrt[2], ...) -- ask numericalize. */
    Expr* approx = numericalize(e, numeric_machine_spec());
    if (!approx) return false;
    bool ok = false;
    if (approx->type == EXPR_INTEGER) {
        *out = (double)approx->data.integer; ok = true;
    } else if (approx->type == EXPR_REAL) {
        *out = approx->data.real; ok = !isnan(*out) && !isinf(*out);
    } else if (approx->type == EXPR_BIGINT) {
        *out = mpz_get_d(approx->data.bigint); ok = true;
#ifdef USE_MPFR
    } else if (approx->type == EXPR_MPFR) {
        *out = mpfr_get_d(approx->data.mpfr, MPFR_RNDN);
        ok = !isnan(*out) && !isinf(*out);
#endif
    } else if (is_rational(approx, &n, &d) && d != 0) {
        *out = (double)n / (double)d; ok = true;
    }
    expr_free(approx);
    return ok;
}

/* True if `e` is a Complex[re, im] with an imaginary part that is
 * provably non-zero. Returns false for plain real numbers and for the
 * (rare) Complex shape whose imaginary part is structurally zero. */
static bool clip_has_imaginary_part(Expr* e) {
    Expr *re = NULL, *im = NULL;
    if (!is_complex(e, &re, &im)) return false;
    if (im->type == EXPR_INTEGER && im->data.integer == 0) return false;
    if (im->type == EXPR_REAL && im->data.real == 0.0) return false;
    return true;
}

/* +1 for Infinity, -1 for -Infinity, 0 otherwise. Recognises the
 * bare symbol `Infinity` as well as the standard Times[c, Infinity]
 * negative-infinity form. */
static int clip_classify_infinity(Expr* e) {
    if (is_infinity_sym(e)) return 1;
    if (is_neg_infinity_form(e)) return -1;
    return 0;
}

Expr* builtin_clip(Expr* res) {
    static uint64_t clip_last_warn = 0;

    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* x = res->data.function.args[0];

    /* Manual first-argument threading: Clip[{x1, x2, ...}, ...]
     *   -> {Clip[x1, ...], Clip[x2, ...], ...}.
     * This is the Listable-style behaviour the spec promises for the
     * first argument, implemented in-builtin so the {min, max} and
     * {vmin, vmax} configuration lists are NOT threaded over. */
    if (x->type == EXPR_FUNCTION
        && x->data.function.head->type == EXPR_SYMBOL
        && x->data.function.head->data.symbol.name == SYM_List) {
        size_t n = x->data.function.arg_count;
        Expr** out = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** call_args = malloc(sizeof(Expr*) * argc);
            call_args[0] = expr_copy(x->data.function.args[i]);
            for (size_t k = 1; k < argc; k++) {
                call_args[k] = expr_copy(res->data.function.args[k]);
            }
            Expr* call = expr_new_function(expr_new_symbol(SYM_Clip), call_args, argc);
            free(call_args);
            out[i] = eval_and_free(call);
        }
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), out, n);
        if (out) free(out);
        return list;
    }

    /* Default interval [-1, 1] and identity replacements. */
    Expr* min_expr = NULL;
    Expr* max_expr = NULL;
    Expr* vmin_expr = NULL;
    Expr* vmax_expr = NULL;
    bool min_owned = false, max_owned = false;

    if (argc >= 2) {
        Expr* iv = res->data.function.args[1];
        if (iv->type != EXPR_FUNCTION
            || iv->data.function.head->type != EXPR_SYMBOL
            || iv->data.function.head->data.symbol.name != SYM_List
            || iv->data.function.arg_count != 2) {
            return NULL;
        }
        min_expr = iv->data.function.args[0];
        max_expr = iv->data.function.args[1];
    } else {
        min_expr = expr_new_integer(-1); min_owned = true;
        max_expr = expr_new_integer(1);  max_owned = true;
    }

    if (argc == 3) {
        Expr* rv = res->data.function.args[2];
        if (rv->type != EXPR_FUNCTION
            || rv->data.function.head->type != EXPR_SYMBOL
            || rv->data.function.head->data.symbol.name != SYM_List
            || rv->data.function.arg_count != 2) {
            if (min_owned) expr_free(min_expr);
            if (max_owned) expr_free(max_expr);
            return NULL;
        }
        vmin_expr = rv->data.function.args[0];
        vmax_expr = rv->data.function.args[1];
    } else {
        vmin_expr = min_expr;
        vmax_expr = max_expr;
    }

    Expr* result = NULL;

    /* Reject complex-valued x. Matches Mathematica's Clip::ncompl. */
    if (clip_has_imaginary_part(x)) {
        matsol_warn_once(&clip_last_warn, res,
                         "Clip::ncompl: Symbolic or noncomplex numerical "
                         "arguments are expected.\n");
        goto cleanup;
    }

    /* Infinity is handled before numericalize because Infinity itself
     * is not a machine real. */
    int inf = clip_classify_infinity(x);
    if (inf > 0) {
        result = expr_copy(vmax_expr);
        goto cleanup;
    }
    if (inf < 0) {
        result = expr_copy(vmin_expr);
        goto cleanup;
    }

    /* Disallow degenerate complex bounds the same way. */
    if (clip_has_imaginary_part(min_expr) || clip_has_imaginary_part(max_expr)) {
        goto cleanup;
    }

    double dx, dmin, dmax;
    if (!clip_to_double_value(x, &dx)
        || !clip_to_double_value(min_expr, &dmin)
        || !clip_to_double_value(max_expr, &dmax)) {
        /* Can't decide -- leave unevaluated. */
        goto cleanup;
    }

    if (dx < dmin) {
        result = expr_copy(vmin_expr);
    } else if (dx > dmax) {
        result = expr_copy(vmax_expr);
    } else {
        result = expr_copy(x);
    }

cleanup:
    if (min_owned) expr_free(min_expr);
    if (max_owned) expr_free(max_expr);
    return result;
}

/*
 * builtin_releasehold:
 * ReleaseHold[expr] removes Hold, HoldForm, HoldPattern, and HoldComplete
 * wrappers from expr. Only removes one layer; does not strip nested holds.
 */
Expr* builtin_releasehold(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    return release_hold_recursive(arg);
}

/*
 * builtin_evaluate:
 * Evaluate[expr] causes expr to be evaluated even if it appears as the
 * argument of a function whose attributes specify that it should be held
 * unevaluated. When Evaluate appears outside a held context, it acts as
 * identity since args are already evaluated by the standard pipeline.
 */
Expr* builtin_evaluate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return expr_copy(res->data.function.args[0]);
}

Expr* builtin_append(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (expr->type != EXPR_FUNCTION) return NULL;
    
    size_t new_count = expr->data.function.arg_count + 1;
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    for (size_t i = 0; i < expr->data.function.arg_count; i++) {
        new_args[i] = expr_copy(expr->data.function.args[i]);
    }
    new_args[new_count - 1] = expr_copy(elem);
    Expr* final_res = expr_new_function(expr_copy(expr->data.function.head), new_args, new_count);
    free(new_args);
    return final_res;
}

Expr* builtin_prepend(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (expr->type != EXPR_FUNCTION) return NULL;
    
    size_t new_count = expr->data.function.arg_count + 1;
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    new_args[0] = expr_copy(elem);
    for (size_t i = 0; i < expr->data.function.arg_count; i++) {
        new_args[i + 1] = expr_copy(expr->data.function.args[i]);
    }
    Expr* final_res = expr_new_function(expr_copy(expr->data.function.head), new_args, new_count);
    free(new_args);
    return final_res;
}

/* Assign new_val back to the target `sym`. For a plain symbol this is a direct
 * OwnValue update; for any other assignable form (e.g. Options[f], a Part) it
 * routes through Set so the appropriate handler in apply_assignment runs. In
 * the symbol case it adopts new_val and returns it; otherwise it consumes
 * new_val (handed to Set) and returns the assigned value. */
static Expr* inplace_assign_back(Expr* sym, Expr* new_val) {
    if (sym->type == EXPR_SYMBOL) {
        symtab_add_own_value(sym->data.symbol.name, sym, new_val);
        return new_val;
    }
    Expr** set_args = malloc(sizeof(Expr*) * 2);
    set_args[0] = expr_copy(sym);
    set_args[1] = new_val;                 /* adopt */
    Expr* set_expr = expr_new_function(expr_new_symbol(SYM_Set), set_args, 2);
    free(set_args);
    Expr* assigned = evaluate(set_expr);
    expr_free(set_expr);
    return assigned;
}

Expr* builtin_append_to(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* sym = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (sym->type != EXPR_SYMBOL && sym->type != EXPR_FUNCTION) return NULL;

    Expr* current_val = evaluate(sym);
    if (!current_val || current_val->type != EXPR_FUNCTION) {
        if (current_val) expr_free(current_val);
        return NULL;
    }

    size_t new_count = current_val->data.function.arg_count + 1;
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    for (size_t i = 0; i < current_val->data.function.arg_count; i++) {
        new_args[i] = expr_copy(current_val->data.function.args[i]);
    }
    new_args[new_count - 1] = expr_copy(elem);
    Expr* new_val = expr_new_function(expr_copy(current_val->data.function.head), new_args, new_count);
    free(new_args);

    Expr* result = inplace_assign_back(sym, new_val);

    expr_free(current_val);
    return result;
}

Expr* builtin_prepend_to(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* sym = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (sym->type != EXPR_SYMBOL && sym->type != EXPR_FUNCTION) return NULL;

    Expr* current_val = evaluate(sym);
    if (!current_val || current_val->type != EXPR_FUNCTION) {
        if (current_val) expr_free(current_val);
        return NULL;
    }

    size_t new_count = current_val->data.function.arg_count + 1;
    Expr** new_args = malloc(sizeof(Expr*) * new_count);
    new_args[0] = expr_copy(elem);
    for (size_t i = 0; i < current_val->data.function.arg_count; i++) {
        new_args[i + 1] = expr_copy(current_val->data.function.args[i]);
    }
    Expr* new_val = expr_new_function(expr_copy(current_val->data.function.head), new_args, new_count);
    free(new_args);

    Expr* result = inplace_assign_back(sym, new_val);

    expr_free(current_val);
    return result;
}

static Expr* rules_to_list(Rule* r) {
    size_t count = 0;
    Rule* curr = r;
    while (curr) {
        count++;
        curr = curr->next;
    }
    
    Expr** rule_exprs = malloc(sizeof(Expr*) * count);
    curr = r;
    for (size_t i = 0; i < count; i++) {
        Expr** rule_args = malloc(sizeof(Expr*) * 2);
        rule_args[0] = expr_copy(curr->pattern);
        rule_args[1] = expr_copy(curr->replacement);
        rule_exprs[i] = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
        free(rule_args);
        curr = curr->next;
    }
    
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), rule_exprs, count);
    if (count > 0) free(rule_exprs); else { free(rule_exprs); } 
    /* IMPORTANT: Rules in DownValues/OwnValues should be returned unevaluated */
    return list;
}

Expr* builtin_own_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_SYMBOL) return NULL;
    
    Rule* r = symtab_get_own_values(arg->data.symbol.name);
    return rules_to_list(r);
}

Expr* builtin_down_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_SYMBOL) return NULL;
    
    Rule* r = symtab_get_down_values(arg->data.symbol.name);
    return rules_to_list(r);
}

Expr* builtin_out(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_INTEGER) return NULL;
    
    int64_t n = arg->data.integer;
    if (n < 0) {
        Expr* line_sym = expr_new_symbol(SYM_DollarLine);
        Expr* line_expr = evaluate(line_sym);
        expr_free(line_sym);
        if (line_expr->type == EXPR_INTEGER) {
            int64_t current_line = line_expr->data.integer;
            n = current_line + n; 
        }
        expr_free(line_expr);
    }
    
    if (n <= 0) return NULL;
    
    Expr* out_head = expr_new_symbol(SYM_Out);
    Expr* out_arg = expr_new_integer(n);
    Expr* out_call = expr_new_function(out_head, &out_arg, 1);
    
    Expr* val = apply_down_values(out_call);
    expr_free(out_call);
    
    return val;
}

Expr* builtin_identity(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }
    return expr_copy(res->data.function.args[0]);
}

/* Composition[f1, f2, ...] represents the symbolic composition f1 . f2 . ...
 * Composition has the attributes Flat and OneIdentity. The evaluator handles
 * the application Composition[f1, ..., fn][args...] -> f1[f2[...[fn[args...]]]]
 * directly (see eval.c). This builtin only performs the algebraic
 * simplifications: Composition[]  -> Identity, Composition[f] -> f, drop
 * Identity arguments, and cancel adjacent f / InverseFunction[f] pairs. */
Expr* builtin_composition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;

    if (n == 0) return expr_new_symbol(SYM_Identity);
    if (n == 1) return expr_copy(res->data.function.args[0]);

    Expr** args = res->data.function.args;
    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;
    size_t out_count = 0;
    bool changed = false;

    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Identity) {
            changed = true;
            continue;
        }
        out[out_count++] = a; /* shallow ref into res; copied below */
    }

    /* Cancel adjacent f / InverseFunction[f] pairs (in either order). */
    bool pair_changed;
    do {
        pair_changed = false;
        for (size_t i = 0; i + 1 < out_count; i++) {
            Expr* a = out[i];
            Expr* b = out[i + 1];
            bool cancel = false;
            if (a->type == EXPR_FUNCTION &&
                a->data.function.head->type == EXPR_SYMBOL &&
                a->data.function.head->data.symbol.name == SYM_InverseFunction &&
                a->data.function.arg_count == 1 &&
                expr_eq(a->data.function.args[0], b)) {
                cancel = true;
            }
            if (!cancel &&
                b->type == EXPR_FUNCTION &&
                b->data.function.head->type == EXPR_SYMBOL &&
                b->data.function.head->data.symbol.name == SYM_InverseFunction &&
                b->data.function.arg_count == 1 &&
                expr_eq(b->data.function.args[0], a)) {
                cancel = true;
            }
            if (cancel) {
                for (size_t j = i; j + 2 < out_count; j++) out[j] = out[j + 2];
                out_count -= 2;
                pair_changed = true;
                changed = true;
                break;
            }
        }
    } while (pair_changed);

    if (!changed) {
        free(out);
        return NULL;
    }

    if (out_count == 0) {
        free(out);
        return expr_new_symbol(SYM_Identity);
    }
    if (out_count == 1) {
        Expr* result = expr_copy(out[0]);
        free(out);
        return result;
    }

    Expr** copies = malloc(sizeof(Expr*) * out_count);
    if (!copies) {
        free(out);
        return NULL;
    }
    for (size_t i = 0; i < out_count; i++) copies[i] = expr_copy(out[i]);
    Expr* result = expr_new_function(expr_new_symbol(SYM_Composition), copies, out_count);
    free(copies);
    free(out);
    return result;
}

/* ComposeList[{f1, f2, ..., fn}, x] generates a list of length n+1 of the
 * form {x, f1[x], f2[f1[x]], ..., fn[...f2[f1[x]]...]}. The list of
 * functions can be any expression with head List; the symbolic
 * applications are constructed and the outer evaluator reduces them
 * to fixed point. */
Expr* builtin_compose_list(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) {
        return NULL;
    }
    Expr* fns = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    if (fns->type != EXPR_FUNCTION ||
        fns->data.function.head->type != EXPR_SYMBOL ||
        fns->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }

    size_t n = fns->data.function.arg_count;
    size_t out_count = n + 1;
    Expr** out = malloc(sizeof(Expr*) * out_count);
    if (!out) return NULL;

    out[0] = expr_copy(x);
    for (size_t i = 0; i < n; i++) {
        Expr* fi = fns->data.function.args[i];
        Expr* arg = expr_copy(out[i]);
        out[i + 1] = expr_new_function(expr_copy(fi), &arg, 1);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, out_count);
    free(out);
    return result;
}

Expr* builtin_atomq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_FUNCTION) {
        if (arg->data.function.head->type == EXPR_SYMBOL) {
            const char* head_name = arg->data.function.head->data.symbol.name;
            if (head_name == SYM_Complex || head_name == SYM_Rational) {
                return expr_new_symbol(SYM_True);
            }
        }
        return expr_new_symbol(SYM_False);
    }

    return expr_new_symbol(SYM_True);
}

Expr* builtin_numberq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || arg->type == EXPR_BIGINT) {
        return expr_new_symbol(SYM_True);
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        return expr_new_symbol(SYM_True);
    }
#endif

    if (arg->type == EXPR_FUNCTION) {
        if (arg->data.function.head->type == EXPR_SYMBOL) {
            const char* head_name = arg->data.function.head->data.symbol.name;
            if (head_name == SYM_Complex || head_name == SYM_Rational) {
                return expr_new_symbol(SYM_True);
            }
        }
    }

    return expr_new_symbol(SYM_False);
}

static bool is_numeric_quantity(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol.name;
        if (name == SYM_Pi || name == SYM_E || name == SYM_I ||
            name == SYM_Infinity || name == SYM_ComplexInfinity ||
            name == SYM_EulerGamma || name == SYM_GoldenRatio ||
            name == SYM_Catalan || name == SYM_Degree ||
            name == SYM_GoldenAngle || name == SYM_Glaisher ||
            name == SYM_Khinchin) {
            return true;
        }
        return false;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_SYMBOL) {
            const char* head_name = e->data.function.head->data.symbol.name;
            if (head_name == SYM_Complex || head_name == SYM_Rational) return true;
            
            SymbolDef* def = symtab_get_def(head_name);
            if (def && (def->attributes & ATTR_NUMERICFUNCTION)) {
                for (size_t i = 0; i < e->data.function.arg_count; i++) {
                    if (!is_numeric_quantity(e->data.function.args[i])) return false;
                }
                return true;
            }
        }
        return false;
    }
    return false;
}

Expr* builtin_numericq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL; 
    }

    Expr* arg = res->data.function.args[0];
    if (is_numeric_quantity(arg)) {
        return expr_new_symbol(SYM_True);
    }
    return expr_new_symbol(SYM_False);
}

/* Classify a numeric-quantity expression by reality and sign. Returns:
 *   1  -> arg is real; *sign set to -1, 0, or +1
 *   0  -> arg numericalizes to a genuinely complex (non-real) value
 *  -1  -> arg could not be numericalized to a concrete number
 * Exact integer/rational/real (and MPFR) inputs are read directly; everything
 * else is numericalized to machine precision, mirroring is_numeric_real() in
 * complex.c. The caller is expected to have already established that arg is a
 * numeric quantity via is_numeric_quantity(). */
static int numeric_real_sign(Expr* arg, int* sign) {
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL ||
        arg->type == EXPR_BIGINT || is_rational(arg, &n, &d)) {
        *sign = expr_numeric_sign(arg);
        return 1;
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        *sign = mpfr_zero_p(arg->data.mpfr) ? 0 : mpfr_sgn(arg->data.mpfr);
        return 1;
    }
#endif
    Expr* approx = numericalize(arg, numeric_machine_spec());
    if (!approx) return -1;
    int kind;
    if (approx->type == EXPR_INTEGER || approx->type == EXPR_REAL ||
        approx->type == EXPR_BIGINT || is_rational(approx, &n, &d)) {
        *sign = expr_numeric_sign(approx);
        kind = 1;
    }
#ifdef USE_MPFR
    else if (approx->type == EXPR_MPFR) {
        *sign = mpfr_zero_p(approx->data.mpfr) ? 0 : mpfr_sgn(approx->data.mpfr);
        kind = 1;
    }
#endif
    else {
        kind = 0;  /* Complex (or otherwise non-real) numeric value. */
    }
    expr_free(approx);
    return kind;
}

/* Positive[x] gives True when x is a real positive number, False when x is a
 * manifestly non-positive numeric quantity (negative, zero, or a non-real
 * complex value). Non-numeric arguments are left unevaluated so symbolic
 * expressions flow through the evaluator unchanged. Mirrors Wolfram's
 * Positive. */
Expr* builtin_positive(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        size_t n = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
        fprintf(stderr,
                "Positive::argx: Positive called with %zu argument%s; "
                "1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    /* Only decide for numeric quantities; symbolic arguments stay unevaluated. */
    if (!is_numeric_quantity(arg)) return NULL;

    int sign;
    int kind = numeric_real_sign(arg, &sign);
    if (kind < 0) return NULL;                       /* numeric but un-numericalizable */
    if (kind == 0) return expr_new_symbol(SYM_False);  /* non-real complex */
    return expr_new_symbol(sign > 0 ? "True" : "False");
}

Expr* builtin_negative(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        size_t n = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
        fprintf(stderr,
                "Negative::argx: Negative called with %zu argument%s; "
                "1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    /* Only decide for numeric quantities; symbolic arguments stay unevaluated. */
    if (!is_numeric_quantity(arg)) return NULL;

    int sign;
    int kind = numeric_real_sign(arg, &sign);
    if (kind < 0) return NULL;                       /* numeric but un-numericalizable */
    if (kind == 0) return expr_new_symbol(SYM_False);  /* non-real complex */
    return expr_new_symbol(sign < 0 ? "True" : "False");
}

/* NonNegative[x] gives True when x is a real number that is positive or zero,
 * False when x is a manifestly negative real or a non-real complex value.
 * Non-numeric arguments are left unevaluated. Mirrors Wolfram's NonNegative. */
Expr* builtin_nonnegative(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        size_t n = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
        fprintf(stderr,
                "NonNegative::argx: NonNegative called with %zu argument%s; "
                "1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    /* Only decide for numeric quantities; symbolic arguments stay unevaluated. */
    if (!is_numeric_quantity(arg)) return NULL;

    int sign;
    int kind = numeric_real_sign(arg, &sign);
    if (kind < 0) return NULL;                       /* numeric but un-numericalizable */
    if (kind == 0) return expr_new_symbol(SYM_False);  /* non-real complex */
    return expr_new_symbol(sign >= 0 ? "True" : "False");
}

/* NonPositive[x] gives True when x is a real number that is negative or zero,
 * False when x is a manifestly positive real or a non-real complex value.
 * Non-numeric arguments are left unevaluated. Mirrors Wolfram's NonPositive. */
Expr* builtin_nonpositive(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        size_t n = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
        fprintf(stderr,
                "NonPositive::argx: NonPositive called with %zu argument%s; "
                "1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    /* Only decide for numeric quantities; symbolic arguments stay unevaluated. */
    if (!is_numeric_quantity(arg)) return NULL;

    int sign;
    int kind = numeric_real_sign(arg, &sign);
    if (kind < 0) return NULL;                       /* numeric but un-numericalizable */
    if (kind == 0) return expr_new_symbol(SYM_False);  /* non-real complex */
    return expr_new_symbol(sign <= 0 ? "True" : "False");
}

Expr* builtin_integerq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (expr_is_integer_like(arg)) {
        return expr_new_symbol(SYM_True);
    }

    return expr_new_symbol(SYM_False);
}

Expr* builtin_information(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    const char* sym_name = NULL;
    if (arg->type == EXPR_SYMBOL) sym_name = arg->data.symbol.name;
    else if (arg->type == EXPR_STRING) sym_name = arg->data.string;
    
    if (!sym_name) return NULL;
    
    const char* doc = symtab_get_docstring(sym_name);
    if (!doc) {
        char buf[256];
        /* Show the short (context-shortened) name in the diagnostic so that
         * symbols reached via the context path read naturally. */
        snprintf(buf, sizeof(buf), "No information available for symbol \"%s\".",
                 context_display_name(sym_name));
        return expr_new_string(buf);
    }
    return expr_new_string(doc);
}

Expr* builtin_evenq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer % 2 == 0) {
            return expr_new_symbol(SYM_True);
        }
    }

    if (arg->type == EXPR_BIGINT) {
        return mpz_even_p(arg->data.bigint) ? expr_new_symbol(SYM_True) : expr_new_symbol(SYM_False);
    }

    return expr_new_symbol(SYM_False);
}

Expr* builtin_oddq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer % 2 != 0) {
            return expr_new_symbol(SYM_True);
        }
    }

    if (arg->type == EXPR_BIGINT) {
        return mpz_odd_p(arg->data.bigint) ? expr_new_symbol(SYM_True) : expr_new_symbol(SYM_False);
    }

    return expr_new_symbol(SYM_False);
}

/* Lift Integer, BigInt, or Rational[Integer|BigInt, Integer|BigInt] into an
 * initialized + canonicalized mpq_t. Returns false (and leaves out cleared)
 * for any other shape, or if the rational has a zero denominator. Used by
 * the BigInt-aware Mod and Quotient paths so Rational[10^50, 3] no longer
 * falls through to a lossy double. */
static bool mod_quot_expr_to_mpq(const Expr* e, mpq_t out) {
    mpq_init(out);
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1UL);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpq_set_z(out, e->data.bigint);
        return true;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den;
        mpz_init(num);
        mpz_init(den);
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        if (mpz_sgn(den) == 0) {
            mpz_clears(num, den, NULL);
            mpq_clear(out);
            return false;
        }
        mpq_set_num(out, num);
        mpq_set_den(out, den);
        mpq_canonicalize(out);
        mpz_clears(num, den, NULL);
        return true;
    }
    mpq_clear(out);
    return false;
}

/* Build an Expr from a canonical mpq_t: returns a normalized Integer (or
 * BigInt) when the denominator is 1, else Rational[num, den] with both
 * components normalized. */
static Expr* mod_quot_mpq_to_expr(const mpq_t q) {
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) {
        return expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_numref(q)));
    }
    Expr* num = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_numref(q)));
    Expr* den = expr_bigint_normalize(expr_new_bigint_from_mpz(mpq_denref(q)));
    Expr* args[2] = { num, den };
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

Expr* builtin_mod(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 2 && res->data.function.arg_count != 3)) {
        return NULL;
    }

    Expr* m_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];

    if (res->data.function.arg_count == 2) {
#ifdef USE_MPFR
        bool m_is_mpfr = (m_expr->type == EXPR_MPFR);
        bool n_is_mpfr = (n_expr->type == EXPR_MPFR);
#else
        bool m_is_mpfr = false;
        bool n_is_mpfr = false;
#endif
        bool m_is_num = (m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || m_expr->type == EXPR_BIGINT || m_is_mpfr || is_rational_like(m_expr));
        bool n_is_num = (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_REAL || n_expr->type == EXPR_BIGINT || n_is_mpfr || is_rational_like(n_expr));
        if (!m_is_num || !n_is_num) return NULL;

#ifdef USE_MPFR
        /* MPFR-precision Mod: m - n * floor(m / n), computed at the
         * maximum of the input precisions. */
        if (m_is_mpfr || n_is_mpfr) {
            mpfr_prec_t prec = 53;
            if (m_is_mpfr) prec = mpfr_get_prec(m_expr->data.mpfr);
            if (n_is_mpfr) {
                mpfr_prec_t np = mpfr_get_prec(n_expr->data.mpfr);
                if (np > prec) prec = np;
            }
            mpfr_t m, n, q, r;
            mpfr_init2(m, prec);
            mpfr_init2(n, prec);
            mpfr_init2(q, prec);
            mpfr_init2(r, prec);
            if (m_is_mpfr) mpfr_set(m, m_expr->data.mpfr, MPFR_RNDN);
            else if (m_expr->type == EXPR_INTEGER) mpfr_set_si(m, (long)m_expr->data.integer, MPFR_RNDN);
            else if (m_expr->type == EXPR_BIGINT) mpfr_set_z(m, m_expr->data.bigint, MPFR_RNDN);
            else /* EXPR_REAL */ mpfr_set_d(m, m_expr->data.real, MPFR_RNDN);
            if (n_is_mpfr) mpfr_set(n, n_expr->data.mpfr, MPFR_RNDN);
            else if (n_expr->type == EXPR_INTEGER) mpfr_set_si(n, (long)n_expr->data.integer, MPFR_RNDN);
            else if (n_expr->type == EXPR_BIGINT) mpfr_set_z(n, n_expr->data.bigint, MPFR_RNDN);
            else /* EXPR_REAL */ mpfr_set_d(n, n_expr->data.real, MPFR_RNDN);
            if (mpfr_zero_p(n)) {
                mpfr_clears(m, n, q, r, (mpfr_ptr)0);
                return NULL;
            }
            mpfr_div(q, m, n, MPFR_RNDN);
            mpfr_floor(q, q);
            mpfr_mul(r, q, n, MPFR_RNDN);
            mpfr_sub(r, m, r, MPFR_RNDN);
            Expr* out = expr_new_mpfr_bits(prec);
            mpfr_set(out->data.mpfr, r, MPFR_RNDN);
            mpfr_clears(m, n, q, r, (mpfr_ptr)0);
            return out;
        }
#endif

        if ((m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_BIGINT) &&
            (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_BIGINT)) {

            mpz_t m, n, r;
            expr_to_mpz(m_expr, m);
            expr_to_mpz(n_expr, n);
            if (mpz_cmp_ui(n, 0) == 0) {
                mpz_clears(m, n, NULL);
                return NULL;
            }

            mpz_init(r);
            mpz_fdiv_r(r, m, n);

            Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clears(m, n, r, NULL);
            return out;
        } else if (is_rational_like(m_expr) && is_rational_like(n_expr)) {
            /* Rational-aware Mod: at least one of m, n is a Rational (with
             * BigInt-capable components). Computes m - n * Floor[m/n] in
             * mpq. Catches Mod[Rational[10^50, 3], 7] and Mod[100/3, 7],
             * which the integer-only path skipped and the double fallback
             * either rejected or answered with garbage. */
            mpq_t m, n;
            if (!mod_quot_expr_to_mpq(m_expr, m)) return NULL;
            if (!mod_quot_expr_to_mpq(n_expr, n)) { mpq_clear(m); return NULL; }
            if (mpq_sgn(n) == 0) { mpq_clear(m); mpq_clear(n); return NULL; }

            mpq_t qmpq, prod, rmpq;
            mpz_t qint;
            mpq_inits(qmpq, prod, rmpq, NULL);
            mpz_init(qint);
            mpq_div(qmpq, m, n);
            mpz_fdiv_q(qint, mpq_numref(qmpq), mpq_denref(qmpq));
            mpq_set_z(prod, qint);
            mpq_mul(prod, prod, n);
            mpq_sub(rmpq, m, prod);
            mpq_canonicalize(rmpq);
            Expr* out = mod_quot_mpq_to_expr(rmpq);
            mpz_clear(qint);
            mpq_clears(m, n, qmpq, prod, rmpq, NULL);
            return out;
        } else if (m_expr->type == EXPR_BIGINT || n_expr->type == EXPR_BIGINT || m_expr->type == EXPR_INTEGER || n_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || n_expr->type == EXPR_REAL) {
            double m_val = (m_expr->type == EXPR_REAL) ? m_expr->data.real : (m_expr->type == EXPR_INTEGER) ? (double)m_expr->data.integer : (m_expr->type == EXPR_BIGINT) ? mpz_get_d(m_expr->data.bigint) : 0.0;
            double n_val = (n_expr->type == EXPR_REAL) ? n_expr->data.real : (n_expr->type == EXPR_INTEGER) ? (double)n_expr->data.integer : (n_expr->type == EXPR_BIGINT) ? mpz_get_d(n_expr->data.bigint) : 0.0;
            if (n_val == 0.0) return NULL;
            double result = m_val - n_val * floor(m_val / n_val);
            return expr_new_real(result);
        }
    } else {
        Expr* d_expr = res->data.function.args[2];
        bool m_is_num = (m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || m_expr->type == EXPR_BIGINT || is_rational_like(m_expr));
        bool n_is_num = (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_REAL || n_expr->type == EXPR_BIGINT || is_rational_like(n_expr));
        bool d_is_num = (d_expr->type == EXPR_INTEGER || d_expr->type == EXPR_REAL || d_expr->type == EXPR_BIGINT || is_rational_like(d_expr));
        if (!m_is_num || !n_is_num || !d_is_num) return NULL;

        if ((m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_BIGINT) && 
            (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_BIGINT) && 
            (d_expr->type == EXPR_INTEGER || d_expr->type == EXPR_BIGINT)) {
            
            mpz_t m, n, d, m_minus_d, r;
            expr_to_mpz(m_expr, m);
            expr_to_mpz(n_expr, n);
            expr_to_mpz(d_expr, d);
            
            if (mpz_cmp_ui(n, 0) == 0) {
                mpz_clears(m, n, d, NULL);
                return NULL;
            }
            
            mpz_inits(m_minus_d, r, NULL);
            mpz_sub(m_minus_d, m, d);
            mpz_fdiv_r(r, m_minus_d, n);
            mpz_add(r, r, d);

            Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clears(m, n, d, m_minus_d, r, NULL);
            return out;
        } else if (is_rational_like(m_expr) && is_rational_like(n_expr) && is_rational_like(d_expr)) {
            /* Rational-aware 3-arg Mod: result = d + ((m-d) mod n). */
            mpq_t m, n, d;
            if (!mod_quot_expr_to_mpq(m_expr, m)) return NULL;
            if (!mod_quot_expr_to_mpq(n_expr, n)) { mpq_clear(m); return NULL; }
            if (!mod_quot_expr_to_mpq(d_expr, d)) { mpq_clear(m); mpq_clear(n); return NULL; }
            if (mpq_sgn(n) == 0) { mpq_clears(m, n, d, NULL); return NULL; }

            mpq_t diff, qmpq, prod, rmpq;
            mpz_t qint;
            mpq_inits(diff, qmpq, prod, rmpq, NULL);
            mpz_init(qint);
            mpq_sub(diff, m, d);
            mpq_div(qmpq, diff, n);
            mpz_fdiv_q(qint, mpq_numref(qmpq), mpq_denref(qmpq));
            mpq_set_z(prod, qint);
            mpq_mul(prod, prod, n);
            mpq_sub(rmpq, diff, prod);
            mpq_add(rmpq, rmpq, d);
            mpq_canonicalize(rmpq);
            Expr* out = mod_quot_mpq_to_expr(rmpq);
            mpz_clear(qint);
            mpq_clears(m, n, d, diff, qmpq, prod, rmpq, NULL);
            return out;
        } else if (m_expr->type == EXPR_BIGINT || n_expr->type == EXPR_BIGINT || d_expr->type == EXPR_BIGINT || m_expr->type == EXPR_INTEGER || n_expr->type == EXPR_INTEGER || d_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || n_expr->type == EXPR_REAL || d_expr->type == EXPR_REAL) {
            double m_val = (m_expr->type == EXPR_REAL) ? m_expr->data.real : (m_expr->type == EXPR_INTEGER) ? (double)m_expr->data.integer : (m_expr->type == EXPR_BIGINT) ? mpz_get_d(m_expr->data.bigint) : 0.0;
            double n_val = (n_expr->type == EXPR_REAL) ? n_expr->data.real : (n_expr->type == EXPR_INTEGER) ? (double)n_expr->data.integer : (n_expr->type == EXPR_BIGINT) ? mpz_get_d(n_expr->data.bigint) : 0.0;
            double d_val = (d_expr->type == EXPR_REAL) ? d_expr->data.real : (d_expr->type == EXPR_INTEGER) ? (double)d_expr->data.integer : (d_expr->type == EXPR_BIGINT) ? mpz_get_d(d_expr->data.bigint) : 0.0;
            if (n_val == 0.0) return NULL;
            double m_minus_d = m_val - d_val;
            double mod_val = m_minus_d - n_val * floor(m_minus_d / n_val);
            double result = d_val + mod_val;
            return expr_new_real(result);
        }
    }
    return NULL;
}

typedef struct { double re; double im; } Cplx;

static bool get_numeric_as_complex(Expr* e, Cplx* out) {
    if (e->type == EXPR_INTEGER) {
        *out = (Cplx){ .re = (double)e->data.integer, .im = 0.0 };
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        *out = (Cplx){ .re = mpz_get_d(e->data.bigint), .im = 0.0 };
        return true;
    }
    if (e->type == EXPR_REAL) {
        *out = (Cplx){ .re = e->data.real, .im = 0.0 };
        return true;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *out = (Cplx){ .re = mpfr_get_d(e->data.mpfr, MPFR_RNDN), .im = 0.0 };
        return true;
    }
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        *out = (Cplx){ .re = (double)n / d, .im = 0.0 };
        return true;
    }
    Expr *re_expr, *im_expr;
    if (is_complex(e, &re_expr, &im_expr)) {
        Cplx re_c, im_c;
        if (!get_numeric_as_complex(re_expr, &re_c) || !get_numeric_as_complex(im_expr, &im_c)) {
            return false; 
        }
        *out = (Cplx){ .re = re_c.re, .im = im_c.re };
        return true;
    }
    return false; 
}

Expr* builtin_quotient(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 2 && res->data.function.arg_count != 3)) {
        return NULL;
    }

    Expr* m_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];
    Expr* d_expr = (res->data.function.arg_count == 3) ? res->data.function.args[2] : NULL;

    bool is_m_complex = is_complex(m_expr, NULL, NULL);
    bool is_n_complex = is_complex(n_expr, NULL, NULL);
    bool is_d_complex = d_expr ? is_complex(d_expr, NULL, NULL) : false;

    if (is_m_complex || is_n_complex || is_d_complex) {
        Cplx m, n, d = {0.0, 0.0};
        if (!get_numeric_as_complex(m_expr, &m) || !get_numeric_as_complex(n_expr, &n)) return NULL;
        if (d_expr && !get_numeric_as_complex(d_expr, &d)) return NULL;

        Cplx m_minus_d = { m.re - d.re, m.im - d.im };
        double n_norm_sq = n.re * n.re + n.im * n.im;
        if (n_norm_sq == 0.0) return NULL;

        Cplx z = {
            .re = (m_minus_d.re * n.re + m_minus_d.im * n.im) / n_norm_sq,
            .im = (m_minus_d.im * n.re - m_minus_d.re * n.im) / n_norm_sq
        };

        Cplx result_cplx = { round(z.re), round(z.im) };

        if (result_cplx.im == 0.0) {
            return expr_new_integer((int64_t)result_cplx.re);
        } else {
            Expr* re_part = expr_new_integer((int64_t)result_cplx.re);
            Expr* im_part = expr_new_integer((int64_t)result_cplx.im);
            Expr* result = make_complex(re_part, im_part);
            return result;
        }
    }

    /* Exact-integer fast path. Covers Integer and BigInt without ever
     * collapsing through a lossy double, which previously silently
     * corrupted Quotient[10^50, 7] to INT64_MIN (and any large-int64
     * input where the double quotient overflowed int64). */
    bool m_int = (m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_BIGINT);
    bool n_int = (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_BIGINT);
    bool d_int = d_expr ? (d_expr->type == EXPR_INTEGER || d_expr->type == EXPR_BIGINT) : true;
    if (m_int && n_int && d_int) {
        mpz_t m, n, d, q;
        expr_to_mpz(m_expr, m);
        expr_to_mpz(n_expr, n);
        mpz_init(d);
        if (d_expr) expr_to_mpz(d_expr, d);
        if (mpz_sgn(n) == 0) {
            mpz_clears(m, n, d, NULL);
            return NULL;
        }
        mpz_init(q);
        mpz_sub(q, m, d);     /* q = m - d */
        mpz_fdiv_q(q, q, n);  /* q = floor((m-d)/n) */
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(q));
        mpz_clears(m, n, d, q, NULL);
        return out;
    }

    /* Rational-aware quotient: any of m, n, d is a Rational with BigInt-
     * capable components. Computes Floor[(m-d)/n] in mpq; result is an
     * integer (denominator drops via floor). Catches
     * Quotient[Rational[10^50, 3], 7] which the integer-only path
     * skipped and the double fallback corrupted via mpz_get_d. */
    if (is_rational_like(m_expr) && is_rational_like(n_expr) &&
        (!d_expr || is_rational_like(d_expr))) {
        mpq_t m, n, d;
        if (!mod_quot_expr_to_mpq(m_expr, m)) return NULL;
        if (!mod_quot_expr_to_mpq(n_expr, n)) { mpq_clear(m); return NULL; }
        if (d_expr) {
            if (!mod_quot_expr_to_mpq(d_expr, d)) { mpq_clear(m); mpq_clear(n); return NULL; }
        } else {
            mpq_init(d);  /* d defaults to 0 */
        }
        if (mpq_sgn(n) == 0) { mpq_clears(m, n, d, NULL); return NULL; }

        mpq_t diff, qmpq;
        mpq_inits(diff, qmpq, NULL);
        mpq_sub(diff, m, d);
        mpq_div(qmpq, diff, n);
        mpz_t qint;
        mpz_init(qint);
        mpz_fdiv_q(qint, mpq_numref(qmpq), mpq_denref(qmpq));
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(qint));
        mpz_clear(qint);
        mpq_clears(m, n, d, diff, qmpq, NULL);
        return out;
    }

    /* Fallback for any input that contains a Real / Rational / etc. The
     * previous double-only path is retained here. */
    double m_val, n_val, d_val = 0.0;
    Cplx temp;
    if (!get_numeric_as_complex(m_expr, &temp)) return NULL;
    m_val = temp.re;
    if (!get_numeric_as_complex(n_expr, &temp)) return NULL;
    n_val = temp.re;
    if (d_expr) {
        if (!get_numeric_as_complex(d_expr, &temp)) return NULL;
        d_val = temp.re;
    }

    if (n_val == 0.0) return NULL;

    int64_t result = (int64_t)floor((m_val - d_val) / n_val);
    return expr_new_integer(result);
}

static int64_t get_expr_depth(Expr* e, bool heads) {
    /* An NDArray is the flat-storage equivalent of a rank-deep nested List,
     * so its Depth matches: a rank-r array has depth r + 1 (the atom level). */
    if (e->type == EXPR_NDARRAY) return (int64_t)e->data.ndarray.rank + 1;
    if (e->type != EXPR_FUNCTION) return 1;

    // Rational and Complex are considered atomic in terms of Depth/Length in Mathematica
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) return 1;
    }

    int64_t max_d = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int64_t d = get_expr_depth(e->data.function.args[i], heads);
        if (d > max_d) max_d = d;
    }
    if (heads) {
        int64_t d_head = get_expr_depth(e->data.function.head, heads);
        if (d_head > max_d) max_d = d_head;
    }
    return max_d + 1;
}

Expr* builtin_depth(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    Expr* e = res->data.function.args[0];
    bool heads = false;
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule) {
            if (opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
            }
        }
    }
    return expr_new_integer(get_expr_depth(e, heads));
}

int64_t leaf_count_internal(Expr* e, bool heads) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    
    int64_t count = 0;
    if (heads) {
        count += leaf_count_internal(e->data.function.head, heads);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        count += leaf_count_internal(e->data.function.args[i], heads);
    }
    
    if (!heads && count == 0 && e->data.function.arg_count == 0) return 0;
    return count;
}

Expr* builtin_leafcount(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    
    Expr* expr = res->data.function.args[0];
    bool heads = true;
    
    if (res->data.function.arg_count == 2) {
        Expr* opt = res->data.function.args[1];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->type == EXPR_SYMBOL && opt->data.function.head->data.symbol.name == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol.name == SYM_False) {
                    heads = false;
                }
            }
        }
    }
    
    int64_t count = leaf_count_internal(expr, heads);
    return expr_new_integer(count);
}

static int64_t byte_count_internal(Expr* e) {
    if (!e) return 0;
    int64_t total = sizeof(Expr);
    switch (e->type) {
        case EXPR_SYMBOL:
            if (e->data.symbol.name) total += strlen(e->data.symbol.name) + 1;
            break;
        case EXPR_STRING:
            if (e->data.string) total += strlen(e->data.string) + 1;
            break;
        case EXPR_FUNCTION:
            if (e->data.function.head) total += byte_count_internal(e->data.function.head);
            if (e->data.function.arg_count > 0 && e->data.function.args) {
                total += sizeof(Expr*) * e->data.function.arg_count;
                for (size_t i = 0; i < e->data.function.arg_count; i++) {
                    total += byte_count_internal(e->data.function.args[i]);
                }
            }
            break;
        case EXPR_BIGINT:
            /* GMP significand: mpz_size limbs of mp_limb_t each. */
            total += (int64_t)(mpz_size(e->data.bigint) * sizeof(mp_limb_t));
            break;
        case EXPR_NDARRAY:
            /* Two owned allocations beyond the Expr node: the dims[] array
             * (rank int64_t entries) and the flat row-major data buffer
             * (element count * bytes-per-element for the dtype). The buffer
             * is the dominant cost and the whole point of the object. */
            total += (int64_t)((size_t)e->data.ndarray.rank * sizeof(int64_t));
            total += (int64_t)(ndarray_size(e) * ndt_elem_size(e->data.ndarray.dtype));
            break;
#ifdef USE_MPFR
        case EXPR_MPFR:
            /* Significand storage scales with the value's precision in bits. */
            total += (int64_t)mpfr_custom_get_size(mpfr_get_prec(e->data.mpfr));
            break;
#endif
        default:
            break;
    }
    return total;
}

Expr* builtin_bytecount(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return expr_new_integer(byte_count_internal(res->data.function.args[0]));
}

static void level_rec(Expr* e, int64_t current_level, int64_t min_l, int64_t max_l, bool heads, Expr*** results, size_t* count, size_t* cap) {
    if (max_l >= 0 && current_level > max_l) return;

    int64_t d = get_expr_depth(e, heads);

    /* Rational[n,d] and Complex[re,im] are atomic per Mathematica's AtomQ
     * semantics: Level must not descend into their parts and must not visit
     * their head, regardless of the Heads option. */
    bool atomic = false;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Rational || h == SYM_Complex) atomic = true;
    }

    if (e->type == EXPR_FUNCTION && !atomic) {
        if (heads) level_rec(e->data.function.head, current_level + 1, min_l, max_l, heads, results, count, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            level_rec(e->data.function.args[i], current_level + 1, min_l, max_l, heads, results, count, cap);
        }
    }

    bool match = true;
    if (min_l >= 0) {
        if (current_level < min_l) match = false;
    } else {
        if (d > -min_l) match = false;
    }

    if (max_l >= 0) {
        if (current_level > max_l) match = false;
    } else {
        if (d < -max_l) match = false;
    }

    if (match) {
        if (*count == *cap) {
            *cap *= 2;
            *results = realloc(*results, sizeof(Expr*) * (*cap));
        }
        (*results)[(*count)++] = expr_copy(e);
    }
}

Expr* builtin_level(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* e = res->data.function.args[0];
    Expr* ls = res->data.function.args[1];

    int64_t min_l = 1, max_l = 1;
    if (ls->type == EXPR_INTEGER) {
        max_l = ls->data.integer;
    } else if (ls->type == EXPR_SYMBOL && ls->data.symbol.name == SYM_Infinity) {
        max_l = 1000000;
    } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol.name == SYM_List) {
        if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
            min_l = max_l = ls->data.function.args[0]->data.integer;
        } else if (ls->data.function.arg_count == 2) {
            if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
            if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
            else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol.name == SYM_Infinity) max_l = 1000000;
        }
    }

    bool heads = false;
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol.name == SYM_Rule) {
            if (opt->data.function.args[0]->data.symbol.name == SYM_Heads) {
                if (opt->data.function.args[1]->data.symbol.name == SYM_True) heads = true;
            }
        }
    }

    size_t count = 0, cap = 16;
    Expr** results = malloc(sizeof(Expr*) * cap);
    level_rec(e, 0, min_l, max_l, heads, &results, &count, &cap);

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, count);
    free(results);
    return list;
}

/*
 * Increment / Decrement / PreIncrement / PreDecrement / AddTo / SubtractFrom
 *
 * All six operators follow the same pattern:
 *   1. Resolve the ultimate symbol that holds the mutable value. For a plain
 *      symbol argument this is the symbol itself; for a Part[sym, ...]
 *      argument it is the underlying sym.
 *   2. If that symbol has no OwnValue, emit an X::rvalue message and leave
 *      the expression unevaluated (returning NULL from the builtin).
 *   3. Evaluate the l-value to obtain the current value.
 *   4. Combine it with the delta (Plus for Increment/AddTo, Plus[v, Times[-1,dx]]
 *      for Decrement/SubtractFrom) and re-evaluate so list threading, symbolic
 *      simplification etc. all happen naturally.
 *   5. Write the new value back via Set, which already knows how to update a
 *      plain symbol or an indexed Part target.
 *   6. Return the old value (Increment/Decrement) or the new value
 *      (PreIncrement/PreDecrement, AddTo, SubtractFrom).
 */

/* Return the name of the symbol ultimately holding the mutable value for an
 * l-value used with Increment-family operators. Supports plain symbols and
 * Part[sym, ...] lvalues. Returns NULL for anything else. */
static const char* lvalue_symbol_name(Expr* lhs) {
    if (!lhs) return NULL;
    if (lhs->type == EXPR_SYMBOL) return lhs->data.symbol.name;
    if (lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.arg_count >= 1 &&
        lhs->data.function.head->data.symbol.name == SYM_Part) {
        return lvalue_symbol_name(lhs->data.function.args[0]);
    }
    return NULL;
}

/* Shared worker. dx is the amount to add (caller owns; we make our own copies).
 * If negate is true, dx is subtracted instead. If pre is true, the new value
 * is returned; otherwise the old value is returned. op_name is used only to
 * compose the X::rvalue diagnostic. */
static Expr* increment_core(Expr* lhs, Expr* dx, bool negate, bool pre, const char* op_name) {
    const char* sym = lvalue_symbol_name(lhs);
    if (!sym || symtab_get_own_values(sym) == NULL) {
        fprintf(stderr,
                "%s::rvalue: %s is not a variable with a value, so its value cannot be changed.\n",
                op_name, sym ? sym : "argument");
        return NULL;
    }

    Expr* old_val = evaluate(lhs);

    /* Build Plus[old_val, dx_or_-dx] and evaluate. */
    Expr* delta_term;
    if (negate) {
        Expr** neg_args = malloc(sizeof(Expr*) * 2);
        neg_args[0] = expr_new_integer(-1);
        neg_args[1] = expr_copy(dx);
        delta_term = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        free(neg_args);
    } else {
        delta_term = expr_copy(dx);
    }
    Expr** plus_args = malloc(sizeof(Expr*) * 2);
    plus_args[0] = expr_copy(old_val);
    plus_args[1] = delta_term;
    Expr* plus_expr = expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
    free(plus_args);
    Expr* new_val = evaluate(plus_expr);
    expr_free(plus_expr);

    /* Write the new value back via Set. Set has HoldFirst so the lhs shape
     * (e.g. Part[list, 2]) is preserved for apply_assignment to handle. */
    Expr** set_args = malloc(sizeof(Expr*) * 2);
    set_args[0] = expr_copy(lhs);
    set_args[1] = expr_copy(new_val);
    Expr* set_call = expr_new_function(expr_new_symbol(SYM_Set), set_args, 2);
    free(set_args);
    Expr* set_result = evaluate(set_call);
    expr_free(set_call);
    if (set_result) expr_free(set_result);

    if (pre) {
        expr_free(old_val);
        return new_val;
    } else {
        expr_free(new_val);
        return old_val;
    }
}

Expr* builtin_increment(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* one = expr_new_integer(1);
    Expr* out = increment_core(res->data.function.args[0], one, false, false, "Increment");
    expr_free(one);
    return out;
}

Expr* builtin_decrement(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* one = expr_new_integer(1);
    Expr* out = increment_core(res->data.function.args[0], one, true, false, "Decrement");
    expr_free(one);
    return out;
}

Expr* builtin_preincrement(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* one = expr_new_integer(1);
    Expr* out = increment_core(res->data.function.args[0], one, false, true, "PreIncrement");
    expr_free(one);
    return out;
}

Expr* builtin_predecrement(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* one = expr_new_integer(1);
    Expr* out = increment_core(res->data.function.args[0], one, true, true, "PreDecrement");
    expr_free(one);
    return out;
}

Expr* builtin_addto(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return increment_core(res->data.function.args[0], res->data.function.args[1], false, true, "AddTo");
}

Expr* builtin_subtractfrom(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return increment_core(res->data.function.args[0], res->data.function.args[1], true, true, "SubtractFrom");
}

Expr* builtin_quotientremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) {
        return NULL;
    }

    Expr* quot_expr = expr_new_function(expr_new_symbol(SYM_Quotient), res->data.function.args, 2);
    Expr* mod_expr = expr_new_function(expr_new_symbol(SYM_Mod), res->data.function.args, 2);

    Expr* quotient = builtin_quotient(quot_expr);
    Expr* remainder = builtin_mod(mod_expr);

    quot_expr->data.function.arg_count = 0;
    mod_expr->data.function.arg_count = 0;
    expr_free(quot_expr);
    expr_free(mod_expr);

    if (!quotient || !remainder) {
        if (quotient) expr_free(quotient);
        if (remainder) expr_free(remainder);
        return NULL;
    }

    Expr** results = malloc(sizeof(Expr*) * 2);
    results[0] = quotient;
    results[1] = remainder;

    Expr* final_res = expr_new_function(expr_new_symbol(SYM_List), results, 2);
    free(results);
    return final_res;
}

/* ToString[expr] and ToString[expr, form].
 *
 * Returns a String containing the printed representation of expr.  The
 * optional second argument selects a form: InputForm (the default),
 * FullForm, or TeXForm.  StandardForm and OutputForm are accepted as
 * aliases for InputForm because Mathilda does not distinguish them when
 * rendering to a flat string.
 *
 * Returns NULL (leaving the call unevaluated) when the form is not one
 * we recognise, so the user sees ToString[..., UnknownForm] in the
 * output rather than a misleading silent fallback.
 */
Expr* builtin_tostring(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* expr = res->data.function.args[0];
    const char* form_name = SYM_InputForm;

    if (argc == 2) {
        Expr* form_arg = res->data.function.args[1];
        if (form_arg->type != EXPR_SYMBOL) return NULL;
        form_name = form_arg->data.symbol.name;
    }

    char* str = NULL;
    if (form_name == SYM_FullForm) {
        str = expr_to_string_fullform(expr);
    } else if (form_name == SYM_TeXForm) {
        Expr* inner[1];
        inner[0] = expr_copy(expr);
        Expr* wrapper = expr_new_function(expr_new_symbol(SYM_TeXForm), inner, 1);
        str = expr_to_string(wrapper);
        expr_free(wrapper);
    } else if (form_name == SYM_InputForm
               || strcmp(form_name, "StandardForm") == 0
               || strcmp(form_name, "OutputForm") == 0) {
        str = expr_to_string(expr);
    } else {
        return NULL;
    }

    if (!str) return NULL;
    Expr* out = expr_new_string(str);
    free(str);
    return out;
}

/* ToExpression[input], ToExpression[input, form], ToExpression[input, form, h].
 *
 * Parses input (a string) and returns the resulting expression for the
 * evaluator to further reduce.  The optional `form` argument may be
 * InputForm, FullForm, or StandardForm — all currently route through the
 * same Mathilda parser because the parser is form-agnostic.  When the
 * three-argument form is used, the head h is wrapped around the parsed
 * expression before it is handed back; common use cases include
 * Hold (to delay evaluation) and HoldComplete.
 *
 * On a parse failure we return the symbol $Failed, matching Mathematica.
 * A non-string input (or unrecognised form) leaves the call unevaluated
 * by returning NULL so the user sees the call shape they wrote.
 */
Expr* builtin_toexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    Expr* input = res->data.function.args[0];
    if (input->type != EXPR_STRING) return NULL;

    if (argc >= 2) {
        Expr* form_arg = res->data.function.args[1];
        if (form_arg->type != EXPR_SYMBOL) return NULL;
        const char* form_name = form_arg->data.symbol.name;
        if (form_name != SYM_InputForm
            && form_name != SYM_FullForm
            && strcmp(form_name, "StandardForm") != 0) {
            return NULL;
        }
    }

    Expr* parsed = parse_expression(input->data.string);
    if (!parsed) {
        return expr_new_symbol(SYM_DollarFailed);
    }

    if (argc == 3) {
        Expr* h = res->data.function.args[2];
        Expr* wrap_args[1];
        wrap_args[0] = parsed;
        Expr* wrapped = expr_new_function(expr_copy(h), wrap_args, 1);
        return wrapped;
    }

    return parsed;
}

/* A single Mathematica symbol-name segment (the part between two backticks,
 * or the only part of an unqualified name) must start with a letter or '$'
 * and continue with letters, digits, or '$'. Empty segments are invalid. */
static bool symbol_segment_is_valid(const char* start, size_t len) {
    if (len == 0) return false;
    unsigned char c0 = (unsigned char)start[0];
    if (!(isalpha(c0) || c0 == '$')) return false;
    for (size_t i = 1; i < len; i++) {
        unsigned char c = (unsigned char)start[i];
        if (!(isalnum(c) || c == '$')) return false;
    }
    return true;
}

/* Validate a string for use as a Symbol[] name argument. Accepts:
 *   - a bare segment  e.g. "x", "Foo$1"
 *   - an absolutely-qualified name  e.g. "a`x", "MyPkg`Private`x"
 *   - a relatively-qualified name with a leading backtick  e.g. "`x"
 * Rejects empty input and any name with an empty or syntactically invalid
 * segment (must start with letter or '$', then letters / digits / '$'). */
static bool symbol_name_is_valid(const char* name) {
    if (!name || !*name) return false;
    const char* p = name;
    if (*p == '`') p++;            /* leading backtick: relative-context form */
    while (*p) {
        const char* seg = p;
        while (*p && *p != '`') p++;
        if (!symbol_segment_is_valid(seg, (size_t)(p - seg))) return false;
        if (*p == '`') p++;
    }
    return true;
}

Expr* builtin_symbol(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    const char* name = arg->data.string;
    if (!symbol_name_is_valid(name)) {
        fprintf(stderr,
                "Symbol::symname: The string \"%s\" cannot be used for a symbol "
                "name. A symbol name must start with a letter followed by "
                "letters and numbers.\n",
                name);
        return NULL;
    }

    char* resolved = context_resolve_name(name);
    Expr* out = expr_new_symbol(resolved ? resolved : name);
    free(resolved);
    return out;
}

/* --------------------------- TimeConstrained ----------------------------
 *
 * TimeConstrained[expr, t]            evaluates expr, aborting after t s
 * TimeConstrained[expr, t, failexpr]  ... and returns failexpr on abort
 * TimeConstrained[expr, Infinity]     imposes no time constraint
 *
 * Without failexpr, an aborted call returns the symbol `$Aborted`.
 *
 * Mechanism (two cooperating layers):
 *
 *   1. SIGPROF + ITIMER_PROF: install a one-shot SIGPROF handler and
 *      an ITIMER_PROF timer; run evaluate() inside sigsetjmp.  If the
 *      CPU-time budget is exhausted, the handler siglongjmp's back
 *      here and we abort.  ITIMER_PROF counts only CPU time of this
 *      process, matching the Mathematica semantic ("CPU time spent
 *      inside the main Mathilda kernel process").  SIGPROF is used
 *      (not SIGALRM) so the intrischnorman pmint timeout -- which is
 *      wall-clock based on SIGALRM/ITIMER_REAL -- and the test-harness
 *      alarm(60) keep working.
 *
 *   2. Wall-clock cooperative deadline: capture
 *      clock_gettime(CLOCK_MONOTONIC) + budget on entry; the evaluator
 *      calls tc_check_deadline() once per rewrite step in its
 *      fixed-point loop.  If the wall-clock deadline has passed, the
 *      check siglongjmp's exactly as the signal handler would.
 *
 * Why both: on real Linux and macOS, ITIMER_PROF is precise and the
 * cooperative check is a cheap no-op.  On hosts that mishandle
 * ITIMER_PROF -- notably WSL 1, whose syscall-translation layer
 * under-counts CPU time and delivers SIGPROF ~15x late -- the
 * cooperative wall-clock backstop catches the deadline at the next
 * rewrite step.  The only case that escapes both layers is a single
 * long-running C builtin (e.g. FactorInteger on a huge composite),
 * which cannot be aborted cooperatively and must wait for the late
 * SIGPROF on broken hosts.
 *
 * Old handler / timer state is saved and restored so calls compose
 * under nesting and don't clobber other timer users in the address
 * space.
 *
 * Memory: the longjmp unwind cannot run destructors, so any Expr nodes
 * the in-flight evaluator allocated leak.  This is the documented
 * Mathematica behaviour ("may give different results on different
 * occasions"); the leak is bounded by the size of the abandoned
 * computation and the next session-level GC is `Quit[]`.
 * ------------------------------------------------------------------- */

static sigjmp_buf  tc_jmp_env;
static volatile sig_atomic_t tc_timed_out = 0;

/* Async-signal-safety guard for the timeout siglongjmp.
 *
 * TimeConstrained aborts an over-budget computation with a siglongjmp taken
 * from the asynchronous SIGPROF handler.  SIGPROF can land on ANY instruction,
 * including deep inside malloc/realloc/free while libmalloc holds an internal
 * zone lock.  Unwinding past that point with siglongjmp never releases the
 * lock, so the next allocation deadlocks -- on macOS libmalloc actively detects
 * the still-owned lock and traps with _os_unfair_lock_recursive_abort (SIGILL).
 * This was reproducible (~1 in 8 runs) on integrands whose algebraic-extension
 * Together/Cancel churns giant GMP integers: the SIGPROF budget expires while a
 * multi-megabyte __gmpz_mul is inside large_malloc.
 *
 * The large allocations that leak the lock originate exclusively from GMP (and
 * MPFR, which shares GMP's allocators), so route GMP through a guarded
 * allocator: while a real malloc/realloc/free is in flight tc_gmp_alloc_busy is
 * set, and the handler DEFERS the jump instead of taking it mid-lock.  The
 * allocator then takes the deferred jump from tc_alloc_safepoint(), once the
 * real call has returned and released libmalloc's lock -- a safe point.  Code
 * outside a GMP allocation keeps the immediate async jump. */
static volatile sig_atomic_t tc_gmp_alloc_busy = 0;

/* Cooperative-backstop state.  tc_deadline_active is read on every
 * rewrite step by tc_check_deadline; when zero, the check returns
 * immediately without even reading the clock.  Both fields are updated
 * only on the main thread between sigsetjmp installs, so no
 * sig_atomic_t / memory-barrier dance is required. */
static struct timespec tc_deadline;
static int             tc_deadline_active = 0;

static void tc_sigprof_handler(int sig) {
    (void)sig;
    tc_timed_out = 1;
    /* If the signal landed inside a GMP allocation, do NOT unwind here: that
     * would abandon libmalloc's held zone lock.  Record the timeout and let
     * tc_alloc_safepoint() take the jump the instant the allocator returns. */
    if (tc_gmp_alloc_busy) return;
    siglongjmp(tc_jmp_env, 1);
}

/* Take the deferred timeout jump from a safe point (no allocator lock held).
 * Called by the guarded GMP allocators right after the real malloc/realloc/free
 * returns.  tc_timed_out is 1 only inside an active, timed-out TimeConstrained
 * body -- it is reset on entry and saved/restored around nesting -- so outside
 * that window this is a single-branch no-op that never touches tc_jmp_env. */
static void tc_alloc_safepoint(void) {
    if (tc_timed_out) siglongjmp(tc_jmp_env, 1);
}

/* GMP memory functions: the system allocator wrapped in the busy guard so an
 * interrupting SIGPROF is deferred past libmalloc's critical section (see the
 * tc_gmp_alloc_busy comment above).  Match __gmp_default_allocate's contract:
 * GMP never passes a NULL pointer to realloc/free and requires allocate to
 * abort rather than return NULL on exhaustion. */
static void* tc_gmp_allocate(size_t size) {
    tc_gmp_alloc_busy = 1;
    void* p = malloc(size);
    tc_gmp_alloc_busy = 0;
    if (!p) { fputs("Mathilda: GMP memory exhausted\n", stderr); abort(); }
    tc_alloc_safepoint();   /* deferred timeout jump, if any (leaks p -- bounded) */
    return p;
}

static void* tc_gmp_reallocate(void* ptr, size_t old_size, size_t new_size) {
    (void)old_size;
    tc_gmp_alloc_busy = 1;
    void* p = realloc(ptr, new_size);
    tc_gmp_alloc_busy = 0;
    if (!p && new_size) { fputs("Mathilda: GMP memory exhausted\n", stderr); abort(); }
    tc_alloc_safepoint();
    return p;
}

static void tc_gmp_deallocate(void* ptr, size_t size) {
    (void)size;
    tc_gmp_alloc_busy = 1;
    free(ptr);
    tc_gmp_alloc_busy = 0;
    tc_alloc_safepoint();
}

/* Install the guarded GMP allocators.  Called once from core_init before any
 * GMP/MPFR use.  The wrappers delegate to the same malloc/realloc/free the GMP
 * default used, so memory allocated before this call is still freed correctly. */
void tc_install_alloc_guard(void) {
    mp_set_memory_functions(tc_gmp_allocate, tc_gmp_reallocate, tc_gmp_deallocate);
}

void tc_check_deadline(void) {
    if (!tc_deadline_active) return;
    struct timespec now;
    /* clock_gettime can in principle fail (EINVAL on unsupported clock);
     * treat any failure as "no time elapsed" so we don't spuriously
     * abort. CLOCK_MONOTONIC is mandatory in POSIX.1-2001 and supported
     * on every host we target (Linux/glibc, macOS 10.12+, WSL 1+2). */
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return;
    if (now.tv_sec < tc_deadline.tv_sec) return;
    if (now.tv_sec == tc_deadline.tv_sec &&
        now.tv_nsec < tc_deadline.tv_nsec) return;
    tc_timed_out = 1;
    siglongjmp(tc_jmp_env, 1);
}

/* Portable "do not inline" hint.  GCC/Clang both define __GNUC__. */
#if defined(__GNUC__)
#  define TC_NOINLINE __attribute__((noinline))
#else
#  define TC_NOINLINE
#endif

/* Evaluate `body` under the already-armed SIGPROF timer / cooperative
 * deadline.  The sigsetjmp lives HERE, isolated from builtin_time_constrained,
 * for two reasons:
 *   1. Correctness — `result` is the only automatic whose value must survive
 *      a siglongjmp, and it is qualified `volatile` so its post-jump value is
 *      well-defined (C99 §7.13.2.1p3).
 *   2. Keeping the setjmp out of the caller stops GCC's conservative
 *      -Wclobbered from flagging every non-volatile local (and inlined helper
 *      temporary) in builtin_time_constrained.
 * On the timeout path evaluate() never returns, so the partially-built tree
 * rooted at `body` is leaked — an unavoidable consequence of unwinding with
 * siglongjmp, identical to the behaviour before this refactor.  Marked
 * TC_NOINLINE so the isolation cannot be undone by the optimiser. */
static TC_NOINLINE Expr* tc_run_guarded(Expr* body) {
    Expr* volatile result = NULL;
    if (sigsetjmp(tc_jmp_env, 1) == 0) {
        result = evaluate(body);
    }
    return result;
}

/* Returns:
 *   +1  success, *out_seconds is a finite, possibly zero, possibly
 *       negative numeric time budget.
 *    0  the argument is the symbol Infinity (or DirectedInfinity[1]):
 *       caller should impose no timeout.
 *   -1  the argument is not a recognisable time -- caller should leave
 *       the call unevaluated. */
static int tc_parse_time(Expr* arg, double* out_seconds) {
    if (!arg) return -1;
    if (arg->type == EXPR_INTEGER) {
        *out_seconds = (double)arg->data.integer;
        return 1;
    }
    if (arg->type == EXPR_REAL) {
        *out_seconds = arg->data.real;
        return 1;
    }
    if (arg->type == EXPR_BIGINT) {
        *out_seconds = mpz_get_d(arg->data.bigint);
        return 1;
    }
    if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_Infinity) {
        return 0;
    }
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
        const char* h = arg->data.function.head->data.symbol.name;
        /* DirectedInfinity[1] is the canonical representation of +Infinity. */
        if (h == SYM_DirectedInfinity && arg->data.function.arg_count == 1) {
            Expr* sign = arg->data.function.args[0];
            if (sign->type == EXPR_INTEGER && sign->data.integer > 0) return 0;
        }
        /* Rational[p, q] with integer numerator/denominator. */
        if (h == SYM_Rational && arg->data.function.arg_count == 2) {
            Expr* n = arg->data.function.args[0];
            Expr* d = arg->data.function.args[1];
            double dn = 0.0, dd = 0.0;
            int ok = 1;
            if      (n->type == EXPR_INTEGER) dn = (double)n->data.integer;
            else if (n->type == EXPR_BIGINT)  dn = mpz_get_d(n->data.bigint);
            else ok = 0;
            if      (d->type == EXPR_INTEGER) dd = (double)d->data.integer;
            else if (d->type == EXPR_BIGINT)  dd = mpz_get_d(d->data.bigint);
            else ok = 0;
            if (ok && dd != 0.0) {
                *out_seconds = dn / dd;
                return 1;
            }
        }
    }
    return -1;
}

Expr* builtin_time_constrained(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* expr_arg = res->data.function.args[0];
    Expr* time_arg = res->data.function.args[1];
    Expr* fail_arg = (argc == 3) ? res->data.function.args[2] : NULL;

    /* HoldAll: every argument was kept literal.  Evaluate only the time
     * budget here so the caller can write e.g. TimeConstrained[..., 2*x]
     * with a globally-set x.  The body and failexpr stay unevaluated
     * until they are actually needed -- failexpr in particular MUST NOT
     * be evaluated unless the timeout fires (Mathematica semantic). */
    Expr* time_eval = evaluate(expr_copy(time_arg));
    double seconds = 0.0;
    int kind = tc_parse_time(time_eval, &seconds);
    expr_free(time_eval);

    if (kind < 0) {
        /* Unparseable time budget: leave the whole call unevaluated. */
        return NULL;
    }

    /* Infinity -> evaluate the body with no constraint.  The caller in
     * evaluate_step frees `res` once we return non-NULL; the builtin
     * convention is "don't double-free." */
    if (kind == 0) {
        return eval_and_free(expr_copy(expr_arg));
    }

    /* Non-positive (zero, negative, NaN) budgets abort immediately. */
    if (!(seconds > 0.0)) {
        if (fail_arg) {
            return eval_and_free(expr_copy(fail_arg));
        }
        return expr_new_symbol(SYM_DollarAborted);
    }

    /* Translate seconds -> (it_value.tv_sec, it_value.tv_usec).  Round
     * sub-microsecond positive budgets up to one microsecond so that
     * setitimer doesn't see (0, 0), which would disarm rather than
     * arm.  Cap the value below INT_MAX seconds (~68 years) to stay
     * inside time_t / suseconds_t range on every platform we target. */
    if (seconds > 2147483.0) seconds = 2147483.0;
    long secs  = (long)seconds;
    long usecs = (long)((seconds - (double)secs) * 1e6);
    if (usecs < 0) usecs = 0;
    if (usecs >= 1000000) { secs += usecs / 1000000; usecs %= 1000000; }
    if (secs == 0 && usecs == 0) usecs = 1;

    /* Save previous SIGPROF handler, ITIMER_PROF state, and tc_jmp_env so
     * nested TimeConstrained calls compose and any external profiler is
     * left exactly as we found it. */
    struct sigaction sa, old_sa;
    struct itimerval old_it;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = tc_sigprof_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGPROF, &sa, &old_sa);

    struct itimerval new_it;
    memset(&new_it, 0, sizeof(new_it));
    new_it.it_value.tv_sec  = (time_t)secs;
    new_it.it_value.tv_usec = (suseconds_t)usecs;
    setitimer(ITIMER_PROF, &new_it, &old_it);

    sigjmp_buf  saved_jmp_env;
    memcpy(&saved_jmp_env, &tc_jmp_env, sizeof(saved_jmp_env));
    int saved_depth = eval_get_recursion_depth();
    volatile sig_atomic_t saved_timed_out = tc_timed_out;
    tc_timed_out = 0;

    /* Arm the cooperative wall-clock deadline as a portability backstop.
     * On hosts where ITIMER_PROF / SIGPROF are reliable, the signal
     * fires first and this state is just read-and-discarded; on hosts
     * where they aren't (WSL 1), tc_check_deadline catches the timeout
     * at the next rewrite step.  Nested calls save/restore so the outer
     * deadline resumes correctly when the inner one returns. */
    struct timespec saved_deadline      = tc_deadline;
    int             saved_deadline_active = tc_deadline_active;
    if (clock_gettime(CLOCK_MONOTONIC, &tc_deadline) == 0) {
        time_t add_secs  = (time_t)secs;
        long   add_nsecs = usecs * 1000L;
        tc_deadline.tv_sec  += add_secs;
        tc_deadline.tv_nsec += add_nsecs;
        if (tc_deadline.tv_nsec >= 1000000000L) {
            tc_deadline.tv_sec  += tc_deadline.tv_nsec / 1000000000L;
            tc_deadline.tv_nsec %= 1000000000L;
        }
        tc_deadline_active = 1;
    } else {
        /* If the monotonic clock is unavailable we silently fall back
         * to signal-only enforcement; this is no worse than before
         * this change. */
        tc_deadline_active = 0;
    }

    Expr* result = tc_run_guarded(expr_copy(expr_arg));

    /* Restore the cooperative-deadline state first, so that any
     * tc_check_deadline call racing during the teardown below sees the
     * outer-call deadline (or no deadline if there was no outer call). */
    tc_deadline        = saved_deadline;
    tc_deadline_active = saved_deadline_active;

    /* Disarm OUR timer first, then reinstall the prior timer and handler.
     * Order matters: leaving our timer armed while we swap the handler
     * back could deliver SIGPROF into someone else's handler. */
    struct itimerval disarm;
    memset(&disarm, 0, sizeof(disarm));
    setitimer(ITIMER_PROF, &disarm, NULL);
    setitimer(ITIMER_PROF, &old_it, NULL);
    sigaction(SIGPROF, &old_sa, NULL);

    /* Restore the outer-call jmp_buf so a parent TimeConstrained can
     * still receive its own SIGPROF correctly. */
    memcpy(&tc_jmp_env, &saved_jmp_env, sizeof(saved_jmp_env));

    bool aborted_here = (tc_timed_out != 0);
    tc_timed_out = saved_timed_out;

    /* If we longjmp'd out of evaluate(), the matching depth decrements
     * never ran -- restore the counter so the rest of the session is
     * not falsely "$RecursionLimit exceeded". */
    if (aborted_here) {
        eval_reset_recursion_depth(saved_depth);
    }

    if (aborted_here) {
        if (result) expr_free(result);   /* defensive; should be NULL */
        if (fail_arg) {
            /* failexpr is evaluated AFTER the timeout context is torn
             * down, with a fresh CPU budget.  Mathematica's semantics
             * say "TimeConstrained evaluates failexpr only if the
             * evaluation is aborted." */
            return eval_and_free(expr_copy(fail_arg));
        }
        return expr_new_symbol(SYM_DollarAborted);
    }

    /* `res` is owned by the caller (evaluate_step) and will be freed
     * after we return.  Returning the evaluated body directly does not
     * double-free anything because `result` is a fresh tree built from
     * expr_copy(expr_arg). */
    return result;
}

