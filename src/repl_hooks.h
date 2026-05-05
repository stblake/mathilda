#ifndef REPL_HOOKS_H
#define REPL_HOOKS_H

#include "expr.h"

/*
 * Mathematica-style REPL session hooks: $PreRead, $Pre, $Post, $PrePrint,
 * $Epilog. Each is a global symbol whose OwnValue (if assigned by the
 * user) is consulted at a specific phase of the REPL loop and applied to
 * that phase's payload.
 *
 *   $PreRead   -- applied to the raw input string before parsing.
 *   $Pre       -- applied to the parsed Expr before evaluation.
 *   $Post      -- applied to the evaluated result after evaluation.
 *   $PrePrint  -- applied to the result expression just before printing.
 *   $Epilog    -- evaluated once when the session terminates.
 *
 * If the corresponding symbol has no OwnValue assigned, the hook is a
 * pass-through. The hooks are exposed as helper functions so they are
 * unit-testable independently of the readline-driven REPL loop.
 */

/*
 * Initialise the hook symbols in the global symbol table. Adds
 * docstrings for ?$PreRead etc. and ensures the symbols exist; does
 * NOT install any default OwnValue (so the hooks default to no-op).
 *
 * Must be called once after symtab_init(), normally from core_init().
 */
void repl_hooks_init(void);

/*
 * Apply $PreRead to an input string.
 *
 * Returns a freshly malloc'd C string the caller must free(). If
 * $PreRead is unset, the result is a strdup of `input`. If $PreRead is
 * set, the function builds $PreRead[input_string], evaluates it, and
 * returns a copy of the resulting String. If the hook returns a
 * non-string value, a $PreRead::strret diagnostic is printed and the
 * original input is returned instead.
 */
char* repl_apply_pre_read(const char* input);

/*
 * Apply $Pre to a parsed expression.
 *
 * Takes ownership of `expr` (consumes one ref); returns a fresh Expr
 * the caller must expr_free(). If $Pre is unset, returns `expr`
 * unchanged (still owned by the caller). If $Pre is set, the function
 * builds $Pre[expr], evaluates it, frees the wrapper call, and returns
 * the result.
 */
Expr* repl_apply_pre(Expr* expr);

/*
 * Apply $Post to an evaluated result. Same ownership semantics as
 * repl_apply_pre.
 */
Expr* repl_apply_post(Expr* expr);

/*
 * Apply $PrePrint to a result expression about to be printed. Same
 * ownership semantics as repl_apply_pre. The return value is what the
 * REPL ultimately renders, but the unmodified result is what gets
 * stored in Out[n].
 */
Expr* repl_apply_pre_print(Expr* expr);

/*
 * If $Epilog has an OwnValue assigned, evaluate the symbol $Epilog
 * once and discard the result. Called from the REPL just before
 * shutdown (Quit[] or EOF). Safe to call when no $Epilog is set; in
 * that case the function does nothing.
 */
void repl_apply_epilog(void);

#endif /* REPL_HOOKS_H */
