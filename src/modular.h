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
void* mth_inherited_block_enter(Expr* const* specs, size_t n);
Expr* builtin_inheritedblock(Expr* res);
void  mth_block_leave(void* frame);
/* Depth of the dynamic-scope stacks (Block frames + binding frames, below). */
int  mth_block_depth_save(void);
void mth_block_depth_unwind(int depth);
/* Binding frames: unwind-safe temporary localisation of ONE variable by a C
 * builtin that binds it itself (an iterator, an integration variable).
 * mth_bind_push detaches the variable's current OwnValues -- a Symbol, or the
 * DownValue of an indexed iterator v[k] -- and returns a token (-1 for a NULL
 * or non-variable `var`, a no-op for mth_bind_pop).  Bind values in between
 * as before (symtab_add_own_value / iter_bind); mth_bind_pop frees what was
 * bound and puts the old value back, returning whether there was one.  It
 * does NOT bump the eval clock: the caller does, as it did before.  The frame
 * lives on the stack the TimeConstrained unwind drains (with Block's frames,
 * newest first), so a timeout can no longer leave the loop value installed.
 * No allocation per push once the stack has grown; pops may be out of order.
 *   MTH_BIND_ATTRS  also save the attributes and restore them on pop;
 *   MTH_BIND_DOWN   also detach and restore the DownValues (symbols only). */
#define MTH_BIND_ATTRS 0x1u
#define MTH_BIND_DOWN  0x2u
int  mth_bind_push(Expr* var, unsigned flags);
int  mth_bind_push_name(const char* name, unsigned flags);
bool mth_bind_pop(int token);

/* Free the Block rule lists parked by a TimeConstrained unwind (see modular.c).
 * Call only at a top-level statement boundary, with no evaluation in flight. */
void mth_block_graveyard_flush(void);

void modular_init(void);

#endif // MODULAR_H
