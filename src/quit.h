#ifndef QUIT_H
#define QUIT_H

#include "expr.h"

/* Quit[] / Quit[n] and its synonym Exit -- end the Mathilda session.
 *
 * Quit is a control-flow sentinel: Quit[...] stays unevaluated so that
 * CompoundExpression, the iteration constructs and FixedPoint stop at it
 * exactly as they stop at Abort[] or Throw[...].  Evaluating it additionally
 * records a process-wide quit request, which each top-level driver (the
 * interactive loop, `-file` script mode, and Get[]'s file loop) polls after
 * every evaluation.  Recording a flag rather than calling exit() keeps
 * $Epilog, stdout flushing and leak accounting on the normal shutdown path. */
Expr* builtin_quit(Expr* res);
Expr* builtin_exit(Expr* res);

/* Nonzero once Quit[] / Quit[n] has been evaluated. */
int mth_quit_requested(void);

/* Process exit status requested by Quit[n] (0 for Quit[] or no request). */
int mth_quit_status(void);

/* Forget a quit request (the REPL ignores Quit[] evaluated by a mod while
 * it starts up). */
void mth_quit_cancel(void);

void quit_init(void);

#endif /* QUIT_H */
