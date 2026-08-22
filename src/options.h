/* options.h — shared option-extraction helpers for builtins.
 *
 * Mathematica-style option arguments are written as
 *   Foo[a, b, OptionName -> value, OtherOption -> other]
 * with options appearing as trailing Rule[]/RuleDelayed[] arguments.
 * The functions in this header walk a function call's argument list and
 * extract specific options without mutating the caller's Expr.
 *
 * Currently exposes:
 *   extract_extension_option — used by PolynomialGCD, PolynomialLCM,
 *     Together, Cancel, Apart, and the BronsteinRational pipeline.
 */
#ifndef MATHILDA_OPTIONS_H
#define MATHILDA_OPTIONS_H

#include "expr.h"
#include <stddef.h>
#include <stdbool.h>

/* Extract a trailing `Rule[Extension, alpha]` (or `RuleDelayed[...]`)
 * argument from `res` if present.
 *
 * Returns the alpha expression (a borrowed pointer into `res`; the
 * caller must NOT free it) when an `Extension` option with a non-trivial
 * value is found.
 *
 * Returns NULL when:
 *   - `res` is NULL or not a function call;
 *   - no trailing `Extension -> _` option is present; or
 *   - the option's value is the symbol `None` or `Automatic` (these are
 *     treated as "no extension"; the `Automatic` case is currently
 *     identical to `None` — true auto-detection of algebraic numbers
 *     in the input is deferred).
 *
 * On entry, `*new_argc` is overwritten with the number of arguments
 * remaining after stripping consumed option arguments — i.e. only
 * `res->data.function.args[0..*new_argc - 1]` should be treated as
 * non-option (polynomial / variable) arguments.
 *
 * Multiple trailing `Extension -> _` rules are all consumed; the
 * effective alpha is the LAST one (rightmost), matching Mathematica's
 * "last setting wins" option semantics.
 */
const Expr* extract_extension_option(const Expr* res, size_t* new_argc);
/* ---- generic named options -------------------------------------------------
 *
 * The Extension helper above reads ONE option for one family. This reads any set of them, so a
 * builtin can accept `Foo[a, b, Name -> v, Other -> w]` without growing a positional tail. It exists
 * because ImageCorners had reached five positional arguments, four of them settings.
 *
 * Options must be TRAILING, as in Mathematica: the first Rule/RuleDelayed argument ends the
 * positional list. `*new_argc` receives how many leading arguments are positional, so the caller
 * parses only those.
 *
 * Values are BORROWED pointers into `res` -- do not free them. An entry whose option does not appear
 * is filled from the symbol's registered defaults (`Options[head]`, via symtab_set_options), so those
 * defaults are the single source of truth and SetOptions works with no further code. An entry still
 * unset is left as the caller initialised it.
 *
 * Returns false on an unknown option name, a non-Rule trailing argument after options have started,
 * or a rule whose left side is not a symbol. Mathematica warns and continues in that case; declining
 * is the more conservative reading and matches this tree's rule of refusing rather than guessing.
 * Last setting wins among duplicates, matching Mathematica.
 */
typedef struct {
    /** Option name, e.g. "MaxFeatures". */
    const char* name;
    /** Set to a borrowed pointer to the option's value, or left alone if absent everywhere. */
    const Expr** out;
    /** Optional. Set to true only when THE CALL supplied this option, false when the value came from
     *  the registered defaults.
     *
     *  This distinction is load-bearing for any head that also accepts the setting positionally. The
     *  first version of this reader omitted it, and the consequence was immediate and silent:
     *  `CornerFilter[img, 2, "Harris"]` computed the MinimumEigenvalue response, because the default
     *  Method filled `out` and the caller could not tell that apart from an explicit option, so the
     *  default overrode the positional argument. Caught by asserting the option form and the
     *  positional form agree -- which is exactly the kind of equivalence worth asserting. */
    bool* given;
} OptEntry;

bool options_extract(const Expr* res, const char* head_name,
                     const OptEntry* entries, size_t n_entries, size_t* new_argc);


/* Variant of `extract_extension_option` that distinguishes
 * `Extension -> Automatic` (explicit auto-detect request) from
 * `Extension -> None` and the absence of any `Extension` option.
 *
 * On return:
 *   - the function's return value is the explicit α (borrowed) when
 *     `Extension -> α` was given with α neither `None` nor `Automatic`,
 *     and NULL otherwise.
 *   - `*automatic_out` (if non-NULL) is set to `true` iff
 *     `Extension -> Automatic` appeared in the trailing options.
 *   - `*new_argc` is overwritten exactly as in `extract_extension_option`.
 *
 * Callers wishing to support `Extension -> Automatic` should call this
 * variant and, when `*automatic_out == true && returned alpha == NULL`,
 * run their own auto-detection (typically `extension_autodetect`).
 *
 * `Extension -> None` overrides any earlier `Extension -> Automatic`
 * (rightmost-wins) so `Foo[poly, Extension -> Automatic, Extension -> None]`
 * leaves `*automatic_out == false`.
 */
const Expr* extract_extension_option_full(const Expr* res, size_t* new_argc,
                                          bool* automatic_out);

/* ----------------------------------------------------------------------
 * Options / SetOptions / OptionValue builtins (src/options_builtin.c)
 * -------------------------------------------------------------------- */

/* Options[sym] / Options[expr] / Options[obj, name] / Options[obj, {names}].
 * For a symbol, returns a copy of its registered default options (or {}); for
 * a compound expression, returns the option rules explicitly present in its
 * arguments. The two-argument forms select named settings. */
Expr* builtin_options(Expr* res);

/* SetOptions[s, name->val, ...] updates s's default options in place (bypassing
 * Protected, refusing Locked) and returns the new Options[s]. Unknown option
 * names emit SetOptions::optnf and leave the call unevaluated. */
Expr* builtin_setoptions(Expr* res);

/* OptionValue[name] / [f,name] / [f,opts,name] / [f,opts,name,Hold] — resolve
 * an option value from an explicit rule list and/or the defaults derived from
 * f. The bare and two-argument forms only resolve inside a rule whose LHS used
 * OptionsPattern (see optionvalue_inject_context). */
Expr* builtin_optionvalue(Expr* res);

/* Register the three builtins, their attributes, and the comprehensive table
 * of default options for option-accepting builtins. Call after every module
 * _init so all option-name symbols are interned. */
void options_builtin_init(void);

/* Functional rewrite used by apply_down_values when a fired DownValue's LHS
 * carried an OptionsPattern. Returns a NEW tree (caller frees the old one) in
 * which each context-dependent OptionValue node is rewritten into the explicit
 * OptionValue[head_sym, opts, name] form so builtin_optionvalue can resolve it.
 * `opts` is the List of options matched by the OptionsPattern. The input tree
 * is never mutated (it may share refcounted nodes with the stored rule). */
Expr* optionvalue_inject_context(const Expr* e, const char* head_sym,
                                 const Expr* opts);

#endif /* MATHILDA_OPTIONS_H */
