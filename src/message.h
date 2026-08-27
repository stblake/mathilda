/*
 * message.h
 *
 * Global message-suppression depth.  While the depth is greater than zero,
 * evaluator diagnostics (Power::infy, Infinity::indet, Solve::nsdim, the
 * FindMinimum / NMinimize warnings, ...) are silenced at their emission sites.
 *
 * The motivating client is FindInstance, which runs speculative internal probes
 * -- Reduce / Solve / NMinimize / FindMinimum over candidate reformulations of
 * the user's statement -- and must do so quietly, exactly as Mathematica
 * evaluates its own internals under an implicit Quiet.  A probe that divides by a
 * sampled zero or hands NMinimize an unsupported constraint shape is expected and
 * harmless; leaking dozens of Power::infy lines to the user is not.
 *
 * The counter is a plain process-global depth: Mathilda's evaluator is
 * single-threaded, and push / pop nest.  Always pair a push with a pop (a single
 * wrapper around the whole quiet region is the safe discipline -- a missed pop
 * would silence every later message).
 */
#ifndef MATHILDA_MESSAGE_H
#define MATHILDA_MESSAGE_H

/* Enter a quiet region (increment the suppression depth). */
void mth_msg_suppress_push(void);
/* Leave a quiet region (decrement; saturates at zero). */
void mth_msg_suppress_pop(void);
/* Non-zero iff messages are currently suppressed. */
int  mth_msg_suppressed(void);

#endif /* MATHILDA_MESSAGE_H */
