/* Mathilda — declared argument counts and option validation for builtins.
 *
 * A builtin declares how many POSITIONAL arguments it takes, and whether it
 * takes options, once at registration:
 *
 *     symtab_set_arity("Length", 1, 1);
 *     symtab_set_arity("Catch", 1, 3);
 *     symtab_set_arity_opts("Map", 1, 3, 2);   // + options, names checked from 2
 *
 * (most declarations are rows of the table in src/arity_table.c) and the
 * evaluator checks every call on the head against it, once, after the head's
 * UpValues and DownValues have declined and before its C code runs
 * (evaluate_step, step 5-pre).  A call outside the declared range prints the
 * reference's argument-count message and comes to rest unevaluated; the C
 * code never sees it.
 *
 * The messages (General::argx, argr, argrx, argt, argtu, argb, argbu, argm,
 * argmu; see builtin_arg_error in common.h):
 *
 *     Length::argx: Length called with 0 arguments; 1 argument is expected.
 *     Catch::argb: Catch called with 4 arguments; between 1 and 3 arguments
 *                  are expected.
 *
 * Options.  A head declared with symtab_set_arity_opts, or one with a
 * non-empty Options[head], takes options: trailing OPTION ARGUMENTS -- a Rule,
 * a RuleDelayed, or a non-empty List of option arguments -- are not counted.
 * A trailing rule in one of the first `min` slots is still positional
 * (Limit[f, x -> 0]) unless every name in it is one of Options[head]
 * (Plot[f, PlotRange -> All] has ONE positional argument: Plot::argr).  A
 * call with too many positional arguments prints General::nonopt
 *
 *     Map::nonopt: Options expected (instead of 2) beyond position 3 in
 *                  Map[f, {a, b}, 1, 2]. An option must be a rule or a list
 *                  of rules.
 *
 * and, when `validate_from` >= 0, an option at argument index (0-based)
 * >= validate_from whose name is not among Options[head] prints
 * General::optx ("Unknown option Foo in call.").  `validate_from` is the first
 * index at which a rule can only be an option: Replace[expr, rules, ...] takes
 * a rule positionally at index 1, so it validates from 2.  Pass -1 when
 * Options[head] does not list every option the head reads.
 *
 * Cost.  The no-defect path is one byte load (an undeclared head), then two
 * compares, then -- for a head with options -- a head-pointer test of the
 * last argument (arity_call_suspect below).  The rest is out of line.
 */
#ifndef MATHILDA_ARITY_H
#define MATHILDA_ARITY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"
#include "symtab.h"
#include "sym_names.h"

#define ARITY_F_SET   0x01u   /* an arity is declared */
#define ARITY_F_OPTS  0x02u   /* trailing options are not positional */
#define ARITY_F_OPTX  0x04u   /* option names are checked against Options[head] */
#define ARITY_F_BETWEEN 0x08u /* a choice of two counts is worded argb, not argt */

#define ARITY_UNBOUNDED 0xFFu     /* arity_max: any number of arguments */

/* Declare the positional argument count [min, max] of `name`; max < 0 means
 * unbounded. */
void symtab_set_arity(const char* name, int min, int max);

/* As symtab_set_arity, for a head that takes options.  validate_from: see
 * the file comment (-1: do not check option names). */
void symtab_set_arity_opts(const char* name, int min, int max, int validate_from);

/* Word a declared choice of two counts as a range (argb: "between 1 and 2
 * arguments are expected"), for the heads the reference words so. */
void symtab_set_arity_between(const char* name);

/* Drop any declaration (for a head whose C code must see every call). */
void symtab_clear_arity(const char* name);

/* True when `e` is an option argument: a two-argument Rule or RuleDelayed,
 * or a non-empty List of option arguments. */
bool arity_is_option_arg(const Expr* e);

/* Out-of-line check of `call` against `d`'s declaration: prints the
 * message and returns true when the call is rejected.  A call it accepts
 * for a head whose option names it checked has its option lists spliced
 * in place (f[x, {a -> 1}] becomes f[x, a -> 1]). */
bool arity_check_slow(const SymbolDef* d, Expr* call);

/* The register-time table (src/arity_table.c): declarations for the builtins
 * of every subsystem, applied once all of them are registered. */
void arity_table_init(void);

/* Could `e` be an option argument?  One level deep; arity_is_option_arg
 * decides. */
static inline bool arity_maybe_option(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    if (hn == SYM_Rule || hn == SYM_RuleDelayed) return true;
    if (hn != SYM_List || e->data.function.arg_count == 0) return false;
    const Expr* f = e->data.function.args[0];
    if (f->type != EXPR_FUNCTION || f->data.function.head->type != EXPR_SYMBOL) return false;
    hn = f->data.function.head->data.symbol.name;
    return hn == SYM_Rule || hn == SYM_RuleDelayed || hn == SYM_List;
}

/* The inline fast test: could `call` (argument count n, arguments `args`) be
 * outside `d`'s declaration?  False is final; true sends the call to
 * arity_check_slow, which decides. */
static inline bool arity_call_suspect(const SymbolDef* d, size_t n, Expr* const* args) {
    if (n < (size_t)(d->arity_fast & 15u) || n > (size_t)(d->arity_fast >> 4)) return true;
    if (n == 0 || (!(d->arity_flags & ARITY_F_OPTS) && !d->default_options)) return false;
    return arity_maybe_option(args[n - 1]);
}

#endif /* MATHILDA_ARITY_H */
