/*
 * message.c -- global message-suppression depth.  See message.h for rationale.
 */
#include "message.h"

/* Suppression nesting depth.  Zero => messages print; > 0 => suppressed. */
static int g_msg_suppress_depth = 0;

void mth_msg_suppress_push(void) { g_msg_suppress_depth++; }

void mth_msg_suppress_pop(void) {
    if (g_msg_suppress_depth > 0) g_msg_suppress_depth--;
}

int mth_msg_suppressed(void) { return g_msg_suppress_depth > 0; }

/* Depth for the `Solve::ifun` advisory scope (see message.h).  Zero => the
 * message may print; > 0 => a complete-solution caller (Reduce) is on the
 * stack, so the "use Reduce" advice is moot and suppressed. */
static int g_ifun_suppress_depth = 0;

void mth_msg_ifun_suppress_push(void) { g_ifun_suppress_depth++; }

void mth_msg_ifun_suppress_pop(void) {
    if (g_ifun_suppress_depth > 0) g_ifun_suppress_depth--;
}

int mth_msg_ifun_suppressed(void) { return g_ifun_suppress_depth > 0; }
