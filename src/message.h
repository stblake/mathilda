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
/* Save / restore the whole suppression depth (used by TimeConstrained to undo a
 * Quiet[] region abandoned by a timeout siglongjmp). */
int  mth_msg_suppress_depth_save(void);
void mth_msg_suppress_depth_load(int d);

/*
 * Narrow scope for the `Solve::ifun` advisory ("Inverse functions are being
 * used by Solve, ... use Reduce for complete solution information").  Distinct
 * from the general suppression above: Reduce IS the complete-solution path, so
 * that one message is self-contradictory when a Solve running *inside* a Reduce
 * emits it -- yet Reduce must still surface genuine Solve diagnostics
 * (Solve::svars / ::nsdim / ::nongen).  builtin_reduce brackets its whole
 * evaluation with a push/pop; only the ifun emission site consults it.
 * Depth-counted and single-threaded, exactly like the general flag above; pair
 * every push with a pop.
 */
void mth_msg_ifun_suppress_push(void);
void mth_msg_ifun_suppress_pop(void);
int  mth_msg_ifun_suppressed(void);

/*
 * Message-fired counter.  Distinct from the suppression depth above: a message
 * still "fires" (and bumps this counter) even while suppressed, so that
 * Check[expr, failexpr] -- typically wrapped as Quiet[Check[expr, failexpr]] --
 * can detect that a diagnostic was generated during expr's evaluation without
 * that diagnostic reaching the user.  Emission sites call mth_msg_note_fired()
 * at the point the message is generated (before the suppression check);
 * builtin_check snapshots mth_msg_fired_count() around the evaluation of its
 * first argument.  Coverage is incremental -- a site that does not yet note is
 * simply invisible to Check, never a wrong answer.
 */
void          mth_msg_note_fired(void);
unsigned long mth_msg_fired_count(void);

/* Registers the Quiet / Check / Message builtins.  Called from core_init(). */
void message_init(void);

#endif /* MATHILDA_MESSAGE_H */
