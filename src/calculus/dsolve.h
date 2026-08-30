/*
 * dsolve.h — symbolic differential-equation solver (DSolve).
 *
 * DSolve is the symbolic analogue of Integrate: a cascade polyalgorithm whose
 * top-level dispatcher (builtin_dsolve) tries a sequence of methods, each of
 * which is ALSO a REPL-callable builtin `DSolve`Method[]`.  The dispatcher and
 * the shared problem substrate live here / in dsolve_common.{c,h}; each method
 * is one file src/calculus/dsolve_<method>.c following the same three-function
 * contract as the Integrate methods.
 *
 * Phase 1 (this file's scope) covers ordinary differential equations; PDEs are
 * Phase 2.  See DSOLVE_PLAN.md for the full method catalog.
 */
#ifndef MATHILDA_DSOLVE_H
#define MATHILDA_DSOLVE_H

#include "../expr.h"

/* Top-level DSolve[eqn, u|u[x]|{...}, x|{x,...}|{x,xmin,xmax}, opts] builtin. */
Expr* builtin_dsolve(Expr* res);

/* Nesting depth of the DSolve method cascade (a DSolve invoked recursively by
 * another method — e.g. a PDE reducing to characteristic ODEs — runs at
 * depth > 0).  Mirrors g_integrate_depth: only the outermost call arms the
 * per-command fail-memo. */
extern int g_dsolve_depth;

/* Register DSolve and fan out to every method's init.  Called from core_init. */
void dsolve_init(void);

#endif /* MATHILDA_DSOLVE_H */
