#include "core.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "arithmetic.h"
#include "comparisons.h"
#include "boolean.h"
#include "list.h"
#include "replace.h"
#include "patterns.h"
#include "cond.h"
#include "iter.h"
#include "complex.h"
#include "trig.h"
#include "hyperbolic.h"
#include "logexp.h"
#include "piecewise.h"
#include "attr.h"
#include "purefunc.h"
#include "modular.h"
#include "sort.h"
#include "stats.h"
#include "info.h"
#include "expand.h"
#include "poly.h"
#include "rat.h"
#include "facint.h"
#include "facpoly.h"
#include "datetime.h"
#include "linalg.h"
#include "load.h"
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
#include <gmp.h>
#include "rat.h"
#include "parfrac.h"
#include "random.h"
#include "picostrings.h"
#include "series.h"
#include "deriv.h"
#include "limit.h"
#include "context.h"
#include "numeric.h"
#include "precision.h"
#include "rationalize.h"
#include "sym_names.h"
#include "repl_hooks.h"

void core_init(void) {
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
    repl_hooks_init();
    parfrac_init();
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
    symtab_add_builtin("Part", builtin_part);
    symtab_add_builtin("Extract", builtin_extract);
    symtab_add_builtin("Head", builtin_head);
    symtab_add_builtin("First", builtin_first);
    symtab_add_builtin("Last", builtin_last);
    symtab_add_builtin("Most", builtin_most);
    symtab_add_builtin("Rest", builtin_rest);

    symtab_get_def("Clear")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
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
    symtab_add_builtin("FreeQ", builtin_freeq);
    symtab_add_builtin("Sort", builtin_sort);
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
    symtab_add_builtin("IntegerQ", builtin_integerq);
    symtab_add_builtin("EvenQ", builtin_evenq);
    symtab_add_builtin("OddQ", builtin_oddq);
    symtab_add_builtin("Mod", builtin_mod);
    symtab_add_builtin("Quotient", builtin_quotient);
    symtab_add_builtin("QuotientRemainder", builtin_quotientremainder);
    symtab_add_builtin("GCD", builtin_gcd);
    symtab_add_builtin("LCM", builtin_lcm);
    symtab_add_builtin("PowerMod", builtin_powermod);
    symtab_add_builtin("Factorial", builtin_factorial);
    symtab_add_builtin("Factorial2", builtin_factorial2);
    symtab_add_builtin("Binomial", builtin_binomial);
    symtab_add_builtin("Print", builtin_print);
    symtab_add_builtin("FullForm", builtin_fullform);
    symtab_add_builtin("InputForm", builtin_inputform);
    symtab_add_builtin("TeXForm", builtin_texform);
    symtab_add_builtin("Information", builtin_information);
    symtab_add_builtin("Evaluate", builtin_evaluate);
    symtab_get_def("Evaluate")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ReleaseHold", builtin_releasehold);
    symtab_get_def("ReleaseHold")->attributes |= ATTR_PROTECTED;

    symtab_get_def("AtomQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NumberQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NumericQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("IntegerQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("EvenQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("OddQ")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Information")->attributes |= ATTR_PROTECTED;

    symtab_get_def("Mod")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Quotient")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("QuotientRemainder")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("GCD")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("LCM")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE | ATTR_FLAT | ATTR_ORDERLESS | ATTR_ONEIDENTITY);
    symtab_get_def("PowerMod")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_get_def("Factorial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Factorial2")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Binomial")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Print")->attributes |= ATTR_PROTECTED;
    symtab_get_def("FullForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("InputForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("TeXForm")->attributes |= ATTR_PROTECTED;
    symtab_get_def("HoldForm")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    facint_init();

    Expr* zero = expr_new_integer(0);
    Expr* one = expr_new_integer(1);
    Expr* sym_I = expr_new_symbol("I");
    Expr* val_I = make_complex(zero, one);
    symtab_add_own_value("I", sym_I, val_I);
    expr_free(sym_I);
    expr_free(val_I);
    
    comparisons_init();
    boolean_init();
    list_init();
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
    attr_init();
    purefunc_init();
    stats_init();
    poly_init();
    facpoly_init();
    rat_init();
    expand_init();
    info_init();
    datetime_init();
    linalg_init();
    load_init();
    random_init();
    strings_init();
    series_init();
    deriv_init();
    limit_init();
    numeric_init();
    precision_init();
    rationalize_init();
}

Expr* builtin_compoundexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 0) return expr_new_symbol("Null");

    Expr* last_val = NULL;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (last_val) expr_free(last_val);
        last_val = evaluate(res->data.function.args[i]);
        if (last_val->type == EXPR_FUNCTION && last_val->data.function.head->type == EXPR_SYMBOL) {
            const char* hname = last_val->data.function.head->data.symbol;
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
                symtab_clear_symbol(arg->data.symbol);
            }
        }
        return expr_new_symbol("Null");
    }
    return NULL;
}

Expr* builtin_length(Expr* res) {
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];
        int64_t len = 0;
        if (arg->type == EXPR_FUNCTION) {
            len = (int64_t)arg->data.function.arg_count;
        }
        return expr_new_integer(len);
    }
    return NULL;
}

static int get_dimensions(Expr* e, int64_t* dims, const char* head_name) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        strcmp(e->data.function.head->data.symbol, head_name) != 0) {
        return 0;
    }
    
    dims[0] = (int64_t)e->data.function.arg_count;
    if (dims[0] == 0) return 1;
    
    int64_t sub_dims[64];
    int sub_depth = get_dimensions(e->data.function.args[0], sub_dims, head_name);
    
    if (sub_depth == 0) return 1;
    
    for (size_t i = 1; i < e->data.function.arg_count; i++) {
        int64_t cur_dims[64];
        int cur_depth = get_dimensions(e->data.function.args[i], cur_dims, head_name);
        
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
    if (res->type == EXPR_FUNCTION && res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];
        if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
            int64_t dims[64];
            int depth = get_dimensions(arg, dims, arg->data.function.head->data.symbol);
            
            Expr** dim_args = malloc(sizeof(Expr*) * depth);
            for (int i = 0; i < depth; i++) {
                dim_args[i] = expr_new_integer(dims[i]);
            }
            
            Expr* ret = expr_new_function(expr_new_symbol("List"), dim_args, depth);
            free(dim_args);
            return ret;
        } else {
            return expr_new_function(expr_new_symbol("List"), NULL, 0);
        }
    }
    return NULL;
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
        is_hold_head(e->data.function.head->data.symbol)) {
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
            Expr* seq = expr_new_function(expr_new_symbol("Sequence"), args, e->data.function.arg_count);
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

Expr* builtin_append_to(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* sym = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (sym->type != EXPR_SYMBOL) return NULL;
    
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
    
    symtab_add_own_value(sym->data.symbol, sym, new_val);
    
    expr_free(current_val);
    return new_val;
}

Expr* builtin_prepend_to(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* sym = res->data.function.args[0];
    Expr* elem = res->data.function.args[1];
    if (sym->type != EXPR_SYMBOL) return NULL;
    
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
    
    symtab_add_own_value(sym->data.symbol, sym, new_val);
    
    expr_free(current_val);
    return new_val;
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
        rule_exprs[i] = expr_new_function(expr_new_symbol("Rule"), rule_args, 2);
        free(rule_args);
        curr = curr->next;
    }
    
    Expr* list = expr_new_function(expr_new_symbol("List"), rule_exprs, count);
    if (count > 0) free(rule_exprs); else { free(rule_exprs); } 
    /* IMPORTANT: Rules in DownValues/OwnValues should be returned unevaluated */
    return list;
}

Expr* builtin_own_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_SYMBOL) return NULL;
    
    Rule* r = symtab_get_own_values(arg->data.symbol);
    return rules_to_list(r);
}

Expr* builtin_down_values(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_SYMBOL) return NULL;
    
    Rule* r = symtab_get_down_values(arg->data.symbol);
    return rules_to_list(r);
}

Expr* builtin_out(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_INTEGER) return NULL;
    
    int64_t n = arg->data.integer;
    if (n < 0) {
        Expr* line_sym = expr_new_symbol("$Line");
        Expr* line_expr = evaluate(line_sym);
        expr_free(line_sym);
        if (line_expr->type == EXPR_INTEGER) {
            int64_t current_line = line_expr->data.integer;
            n = current_line + n; 
        }
        expr_free(line_expr);
    }
    
    if (n <= 0) return NULL;
    
    Expr* out_head = expr_new_symbol("Out");
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

    if (n == 0) return expr_new_symbol("Identity");
    if (n == 1) return expr_copy(res->data.function.args[0]);

    Expr** args = res->data.function.args;
    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;
    size_t out_count = 0;
    bool changed = false;

    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol == SYM_Identity) {
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
                a->data.function.head->data.symbol == SYM_InverseFunction &&
                a->data.function.arg_count == 1 &&
                expr_eq(a->data.function.args[0], b)) {
                cancel = true;
            }
            if (!cancel &&
                b->type == EXPR_FUNCTION &&
                b->data.function.head->type == EXPR_SYMBOL &&
                b->data.function.head->data.symbol == SYM_InverseFunction &&
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
        return expr_new_symbol("Identity");
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
    Expr* result = expr_new_function(expr_new_symbol("Composition"), copies, out_count);
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
        fns->data.function.head->data.symbol != SYM_List) {
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

    Expr* result = expr_new_function(expr_new_symbol("List"), out, out_count);
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
            const char* head_name = arg->data.function.head->data.symbol;
            if (head_name == SYM_Complex || head_name == SYM_Rational) {
                return expr_new_symbol("True");
            }
        }
        return expr_new_symbol("False");
    }

    return expr_new_symbol("True");
}

Expr* builtin_numberq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL; 
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || arg->type == EXPR_BIGINT) {
        return expr_new_symbol("True");
    }

    if (arg->type == EXPR_FUNCTION) {
        if (arg->data.function.head->type == EXPR_SYMBOL) {
            const char* head_name = arg->data.function.head->data.symbol;
            if (head_name == SYM_Complex || head_name == SYM_Rational) {
                return expr_new_symbol("True");
            }
        }
    }

    return expr_new_symbol("False");
}

static bool is_numeric_quantity(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol;
        if (name == SYM_Pi || name == SYM_E || name == SYM_I ||
            name == SYM_Infinity || name == SYM_ComplexInfinity ||
            name == SYM_EulerGamma || name == SYM_GoldenRatio ||
            name == SYM_Catalan || name == SYM_Degree) {
            return true;
        }
        return false;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_SYMBOL) {
            const char* head_name = e->data.function.head->data.symbol;
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
        return expr_new_symbol("True");
    }
    return expr_new_symbol("False");
}

Expr* builtin_integerq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (expr_is_integer_like(arg)) {
        return expr_new_symbol("True");
    }

    return expr_new_symbol("False");
}

Expr* builtin_information(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    const char* sym_name = NULL;
    if (arg->type == EXPR_SYMBOL) sym_name = arg->data.symbol;
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
            return expr_new_symbol("True");
        }
    }

    if (arg->type == EXPR_BIGINT) {
        return mpz_even_p(arg->data.bigint) ? expr_new_symbol("True") : expr_new_symbol("False");
    }

    return expr_new_symbol("False");
}

Expr* builtin_oddq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) {
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer % 2 != 0) {
            return expr_new_symbol("True");
        }
    }

    if (arg->type == EXPR_BIGINT) {
        return mpz_odd_p(arg->data.bigint) ? expr_new_symbol("True") : expr_new_symbol("False");
    }

    return expr_new_symbol("False");
}

Expr* builtin_mod(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 2 && res->data.function.arg_count != 3)) {
        return NULL;
    }

    Expr* m_expr = res->data.function.args[0];
    Expr* n_expr = res->data.function.args[1];

    if (res->data.function.arg_count == 2) {
        bool m_is_num = (m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || m_expr->type == EXPR_BIGINT);
        bool n_is_num = (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_REAL || n_expr->type == EXPR_BIGINT);
        if (!m_is_num || !n_is_num) return NULL;

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
        } else if (m_expr->type == EXPR_BIGINT || n_expr->type == EXPR_BIGINT || m_expr->type == EXPR_INTEGER || n_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || n_expr->type == EXPR_REAL) {
            double m_val = (m_expr->type == EXPR_REAL) ? m_expr->data.real : (m_expr->type == EXPR_INTEGER) ? (double)m_expr->data.integer : (m_expr->type == EXPR_BIGINT) ? mpz_get_d(m_expr->data.bigint) : 0.0;
            double n_val = (n_expr->type == EXPR_REAL) ? n_expr->data.real : (n_expr->type == EXPR_INTEGER) ? (double)n_expr->data.integer : (n_expr->type == EXPR_BIGINT) ? mpz_get_d(n_expr->data.bigint) : 0.0;
            if (n_val == 0.0) return NULL;
            double result = m_val - n_val * floor(m_val / n_val);
            return expr_new_real(result);
        }
    } else { 
        Expr* d_expr = res->data.function.args[2];
        bool m_is_num = (m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_REAL || m_expr->type == EXPR_BIGINT);
        bool n_is_num = (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_REAL || n_expr->type == EXPR_BIGINT);
        bool d_is_num = (d_expr->type == EXPR_INTEGER || d_expr->type == EXPR_REAL || d_expr->type == EXPR_BIGINT);
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
    } else { 
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
}

static int64_t get_expr_depth(Expr* e, bool heads) {
    if (e->type != EXPR_FUNCTION) return 1;

    // Rational and Complex are considered atomic in terms of Depth/Length in Mathematica
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
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
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule) {
            if (opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
            }
        }
    }
    return expr_new_integer(get_expr_depth(e, heads));
}

static int64_t leaf_count_internal(Expr* e, bool heads) {
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
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->type == EXPR_SYMBOL && opt->data.function.head->data.symbol == SYM_Rule && opt->data.function.arg_count == 2) {
            if (opt->data.function.args[0]->type == EXPR_SYMBOL && opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->type == EXPR_SYMBOL && opt->data.function.args[1]->data.symbol == SYM_False) {
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
            if (e->data.symbol) total += strlen(e->data.symbol) + 1;
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

    if (e->type == EXPR_FUNCTION) {
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
    } else if (ls->type == EXPR_SYMBOL && ls->data.symbol == SYM_Infinity) {
        max_l = 1000000;
    } else if (ls->type == EXPR_FUNCTION && ls->data.function.head->data.symbol == SYM_List) {
        if (ls->data.function.arg_count == 1 && ls->data.function.args[0]->type == EXPR_INTEGER) {
            min_l = max_l = ls->data.function.args[0]->data.integer;
        } else if (ls->data.function.arg_count == 2) {
            if (ls->data.function.args[0]->type == EXPR_INTEGER) min_l = ls->data.function.args[0]->data.integer;
            if (ls->data.function.args[1]->type == EXPR_INTEGER) max_l = ls->data.function.args[1]->data.integer;
            else if (ls->data.function.args[1]->type == EXPR_SYMBOL && ls->data.function.args[1]->data.symbol == SYM_Infinity) max_l = 1000000;
        }
    }

    bool heads = false;
    for (size_t i = 2; i < res->data.function.arg_count; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION && opt->data.function.head->data.symbol == SYM_Rule) {
            if (opt->data.function.args[0]->data.symbol == SYM_Heads) {
                if (opt->data.function.args[1]->data.symbol == SYM_True) heads = true;
            }
        }
    }

    size_t count = 0, cap = 16;
    Expr** results = malloc(sizeof(Expr*) * cap);
    level_rec(e, 0, min_l, max_l, heads, &results, &count, &cap);

    Expr* list = expr_new_function(expr_new_symbol("List"), results, count);
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
    if (lhs->type == EXPR_SYMBOL) return lhs->data.symbol;
    if (lhs->type == EXPR_FUNCTION &&
        lhs->data.function.head->type == EXPR_SYMBOL &&
        lhs->data.function.arg_count >= 1 &&
        lhs->data.function.head->data.symbol == SYM_Part) {
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
        delta_term = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
        free(neg_args);
    } else {
        delta_term = expr_copy(dx);
    }
    Expr** plus_args = malloc(sizeof(Expr*) * 2);
    plus_args[0] = expr_copy(old_val);
    plus_args[1] = delta_term;
    Expr* plus_expr = expr_new_function(expr_new_symbol("Plus"), plus_args, 2);
    free(plus_args);
    Expr* new_val = evaluate(plus_expr);
    expr_free(plus_expr);

    /* Write the new value back via Set. Set has HoldFirst so the lhs shape
     * (e.g. Part[list, 2]) is preserved for apply_assignment to handle. */
    Expr** set_args = malloc(sizeof(Expr*) * 2);
    set_args[0] = expr_copy(lhs);
    set_args[1] = expr_copy(new_val);
    Expr* set_call = expr_new_function(expr_new_symbol("Set"), set_args, 2);
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

    Expr* quot_expr = expr_new_function(expr_new_symbol("Quotient"), res->data.function.args, 2);
    Expr* mod_expr = expr_new_function(expr_new_symbol("Mod"), res->data.function.args, 2);

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
    
    Expr* final_res = expr_new_function(expr_new_symbol("List"), results, 2);
    free(results);
    return final_res;
}

