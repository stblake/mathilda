#ifndef MODULAR_H
#define MODULAR_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

Expr* builtin_module(Expr* res);
Expr* builtin_block(Expr* res);
Expr* builtin_with(Expr* res);
Expr* builtin_unique(Expr* res);

/* Capture-avoiding pattern substitution support, used by replace_bindings.
 * expr_is_binding_scope: is `e` a scoping construct (Module/Block/With/Function/
 * Table) that actually binds names? scoping_capture_avoid: if substituting a
 * value carrying any name in danger[] into `e`'s body could capture one of e's
 * bound locals, return a fresh copy with the colliding locals alpha-renamed,
 * else NULL. See docs/notes/MATHILDA_DIVERGENCES.md A11. */
bool expr_is_binding_scope(Expr* e);
Expr* scoping_capture_avoid(Expr* e, const char** danger, size_t ndanger);

/* Block's dynamic-scope frames, exposed for the TimeConstrained unwind only.
 * Block restores its locals' OwnValues, DownValues and attributes on the way
 * out, but a timeout siglongjmps past every C frame in between, so
 * tc_run_guarded records the depth before the body and drains back to it after
 * a jump -- otherwise a Block-installed temporary rule (a rewrite hook, say)
 * would stay installed for the rest of the session. Same pattern as the
 * async-defer count and the message-suppression depth. */
/* C-side Block: localise the symbols named by specs[0..n) (a Symbol, or
 * Set[sym, value] with the value evaluated) exactly as Block[{specs}, ...]
 * does, and restore them. Every enter is paired with one leave; the frame
 * also unwinds with mth_block_depth_unwind like a Block's. Used by Dialog. */
void* mth_block_enter(Expr* const* specs, size_t n);
void  mth_block_leave(void* frame);
int  mth_block_depth_save(void);
void mth_block_depth_unwind(int depth);
/* Free the Block rule lists parked by a TimeConstrained unwind (see modular.c).
 * Call only at a top-level statement boundary, with no evaluation in flight. */
void mth_block_graveyard_flush(void);

void modular_init(void);

#endif // MODULAR_H
