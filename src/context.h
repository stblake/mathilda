#ifndef CONTEXT_H
#define CONTEXT_H

#include "expr.h"
#include "symtab.h"
#include <stdbool.h>
#include <stddef.h>

/* Context system for Mathilda.
 *
 * A context is a namespace identified by a string ending in a backtick, e.g.
 * "Global`" or "MyPackage`Private`". Every symbol conceptually lives in a
 * context, addressed as "context`name".
 *
 * Implementation notes:
 * - $Context holds the current context string; new unqualified symbols are
 *   created in it.
 * - $ContextPath is an ordered list of contexts consulted when resolving a
 *   bare identifier that has not yet been seen.
 * - For backward compatibility with existing code that stores builtins under
 *   bare names (e.g. "Plus"), bare names that already resolve to an existing
 *   builtin are left bare; they are treated as implicitly belonging to
 *   "System`".
 */

/* Initialize/shut down the context subsystem. context_init() must be called
 * from core_init() before any .m file is loaded. */
void context_init(void);
void context_shutdown(void);

/* Current $Context value (read-only borrow of its OwnValue, or "Global`" when
 * that is not a context string; valid until $Context is next assigned). */
const char* context_current(void);

/* True for $Context and $ContextPath, whose Set/SetDelayed go through
 * context_assign so the context state is never left invalid. */
bool context_is_state_var(const char* name);

/* Assign `value` (borrowed, already evaluated) to $Context / $ContextPath.
 * An invalid value emits $Context::cxset / $ContextPath::cxlist and leaves the
 * variable unchanged; returns whether the assignment happened. */
bool context_assign(const char* name, Expr* value);

/* Home context of a symbol stored under the bare name `bare` (def may be
 * NULL): "System`" for a builtin, or for a System-flagged name spelled like a
 * System name (leading capital or `$`); "Global`" otherwise. The single rule
 * behind Context[] and Names, so the two cannot disagree. The spelling guard
 * keeps scratch names the kernel interns (x, t, ...) out of System. */
const char* context_home(const char* bare, const SymbolDef* def);

/* Path access (read-only): the valid context strings of $ContextPath, in
 * order. Pointers borrow the OwnValue (valid until it is next assigned). */
size_t context_path_count(void);
const char* context_path_at(size_t i);

/* Resolve a parsed identifier to its canonical symbol-table name.
 *
 * Rules applied in order:
 *  1. Names starting with '$' are returned unchanged (system variables).
 *  2. A leading backtick is treated as relative to the current context:
 *     "`Private`f"  in context "foo`"  becomes "foo`Private`f".
 *  3. Names containing a backtick elsewhere are absolute. "System`X" is
 *     canonicalized to bare "X" (since builtins are stored under bare names).
 *  4. Otherwise (bare identifier):
 *     a. Search each prefix in { $Context } U $ContextPath. If prefix+name
 *        already exists in the symbol table, return that full name.
 *     b. If a bare symbol with a builtin function is already registered,
 *        return bare (implicit "System`").
 *     c. If $Context == "Global`", return bare (legacy behavior).
 *     d. Otherwise, qualify with $Context.
 *
 * Returned string is malloc'd; caller must free(). */
char* context_resolve_name(const char* raw);

/* Display helper: return a pointer to the short name for printing if the
 * symbol's context prefix is currently visible (equal to $Context or on
 * $ContextPath). Otherwise return `name` unchanged. The pointer is either
 * `name` or points inside it; never NULL; never allocated. */
const char* context_display_name(const char* name);

/* Built-in operations (also exposed for use outside the Context[] builtin). */
bool context_begin(const char* ctx);               /* Begin["ctx`"] */
bool context_begin_package(const char* ctx,
                           const char** needs,
                           size_t needs_count);    /* BeginPackage */
/* On success, *out_ended_ctx is set to a newly malloc'd string containing
 * the just-closed context (caller frees). */
bool context_end(char** out_ended_ctx);
bool context_end_package(char** out_ended_ctx);

/* Depth of the Begin / BeginPackage frame stack, and a forced unwind back to a
 * saved depth (each popped frame restores its saved $Context / $ContextPath).
 * For non-local exits only: a TimeConstrained timeout that siglongjmps out of
 * a package being read (a lazily loaded src/internal package) never reaches its
 * End[] / EndPackage[], and without the unwind the session is left in the
 * package's private context -- every later bare symbol then resolves to a fresh
 * Private` symbol, so user variables appear to vanish. */
size_t context_stack_depth(void);
void   context_unwind_to(size_t depth);
/* The current $Context and $ContextPath values, captured and later put back --
 * for the same non-local exits: internal packages also ASSIGN
 * $ContextPath = {"System`"} before their Begin, which the frame stack alone
 * cannot undo.  The capture is two reference-count bumps (no allocation: it
 * runs on every TimeConstrained body and Integrate probe); context_restore
 * consumes both, context_release drops them on a normal exit.  Either may be
 * NULL. */
void context_capture(Expr** ctx, Expr** path);
void context_restore(Expr* ctx, Expr* path);
void context_release(Expr* ctx, Expr* path);

/* True when `ctx` is listed in $Packages (the contexts of loaded packages;
 * BeginPackage prepends to it).  Needs loads only contexts not yet listed. */
bool context_packages_contains(const char* ctx);

/* Produce Expr* copies of current state, for use as $Context/$ContextPath
 * OwnValues. Caller owns the returned expressions. */
Expr* context_as_string(void);
Expr* context_path_as_list(void);

#endif /* CONTEXT_H */
