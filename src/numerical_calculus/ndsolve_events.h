/* Mathilda — NDSolve events: WhenEvent[event, action].
 *
 * An event is a condition on the solution -- a zero crossing  lhs == rhs
 * (either direction), an inequality that BECOMES true (a > b, a <= b, ...), or
 * any other Boolean condition that becomes True -- together with actions run
 * at the precisely located event time:
 *
 *   "StopIntegration"     the solution's domain ends at the event;
 *   "RemoveEvent"         the event is disarmed for the rest of the run;
 *   u'[t] -> expr, ...    state reset (all right-hand sides are evaluated at
 *                         the pre-event state); integration restarts from the
 *                         event with the new state, and the solution carries a
 *                         discontinuity there (a repeated abscissa);
 *   any other expression  evaluated at the event with t and the state bound
 *                         (Sow[t], a = u[t], Print[...], ...);
 *   {action, ...}         several of the above.
 *
 * The module is integrator-agnostic.  After every ACCEPTED step (ta, Ya) ->
 * (tb, Yb) a loop calls nd_events_locate.  When an event condition changes
 * across the step the event time is found by a safeguarded Illinois
 * (modified regula falsi) iteration on the event function along the step,
 * where the state at an interior time s is supplied by the loop's own
 * `advance` callback: a one-step re-integration from (ta, Ya) for the
 * one-step methods (a Runge-Kutta step of length s - ta, whose local error is
 * below the accepted one), the step's Taylor polynomial for the Taylor
 * integrator, and the method's interpolating polynomial for the multistep
 * methods.  The event time is therefore located to the integrator's own
 * accuracy, not to that of the cubic dense output.  The loop then truncates
 * the step at the event and calls nd_events_fire for the actions.
 *
 * Both a machine (double) and an arbitrary-precision (MPFR) interface are
 * provided; every integration loop (one-step, BDF, Adams; MPFR DOPRI5/RK4,
 * BDF and Taylor) is wired to one of them. */
#ifndef NDSOLVE_EVENTS_H
#define NDSOLVE_EVENTS_H

#include "ndsolve_common.h"

/* Is `e` a WhenEvent[...] item of an NDSolve equation list? */
bool nd_is_when_event(const Expr* e);

/* Compile the WhenEvent items into P->events.  `lits`/`subs` (nsub pairs) map
 * the state forms u[t], u'[t], ... (and the top derivatives) to the reduced
 * state expressions; `fnames`/`base`/`maxorder` (nfun entries) locate the
 * reduced index of a reset target u^(m)[t].  Returns false after a message
 * when an item is malformed (the call is then left unevaluated). */
bool nd_events_build(NdProblem* P, Expr** items, size_t nitems,
                     Expr** lits, Expr** subs, size_t nsub,
                     const char** fnames, const size_t* base,
                     const int* maxorder, size_t nfun);
void nd_events_free(struct NdEvents* E);

/* Results of nd_events_fire. */
#define ND_EVF_STOP  1u   /* "StopIntegration": end the integration here   */
#define ND_EVF_RESET 2u   /* the state was reset: restart from Ynew         */

/* ---- machine precision ---- */

/* State at time s along the current step, from the loop's own method. */
typedef bool (*NdEvAdvance)(void* ctx, double s, double* Ys);

/* Arm every event at the start (t0, Y0) of an integration direction. */
void nd_events_begin(NdProblem* P, double t0, const double* Y0);

/* Check the accepted step (ta, Ya) -> (tb, Yb).  Returns the index of the
 * first event that fires inside it, writing the located time to *te and the
 * state there to Ye (the LAST advance() call is made at *te, so a loop's FSAL
 * stage cache corresponds to it), or -1 (no event; the stored values advance
 * to tb). */
int nd_events_locate(NdProblem* P, double ta, const double* Ya,
                     double tb, const double* Yb,
                     NdEvAdvance adv, void* ctx, double* te, double* Ye);

/* Run event k's actions at (te, Ye).  Ynew receives the post-event state
 * (a copy of Ye when nothing is reset).  Returns ND_EVF_* flags. */
unsigned nd_events_fire(NdProblem* P, int k, double te, const double* Ye, double* Ynew);

#ifdef USE_MPFR
#include <mpfr.h>
typedef bool (*NdEvAdvanceMp)(void* ctx, const mpfr_t s, mpfr_t* Ys);
void nd_events_begin_mp(NdProblem* P, const mpfr_t t0, mpfr_t* Y0, long bits);
int  nd_events_locate_mp(NdProblem* P, const mpfr_t ta, mpfr_t* Ya,
                         const mpfr_t tb, mpfr_t* Yb,
                         NdEvAdvanceMp adv, void* ctx, mpfr_t te, mpfr_t* Ye,
                         long bits, long ubits);
unsigned nd_events_fire_mp(NdProblem* P, int k, const mpfr_t te, mpfr_t* Ye,
                           mpfr_t* Ynew, long bits);
#endif

#endif /* NDSOLVE_EVENTS_H */
