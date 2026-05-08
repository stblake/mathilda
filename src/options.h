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
 *     Together, Cancel, Apart, and the IntegrateRational pipeline.
 */
#ifndef PICOCAS_OPTIONS_H
#define PICOCAS_OPTIONS_H

#include "expr.h"
#include <stddef.h>

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

#endif /* PICOCAS_OPTIONS_H */
