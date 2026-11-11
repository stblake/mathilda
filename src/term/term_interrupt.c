/*
 * term_interrupt.c -- see term_interrupt.h.
 */
#include "term_interrupt.h"

#include "funcprog.h"   /* abort_protect_depth */

volatile sig_atomic_t g_mth_interrupt_pending = 0;

bool mth_interrupt_pending(void) {
    return g_mth_interrupt_pending != 0;
}

void mth_interrupt_request(void) {
    g_mth_interrupt_pending = 1;
}

void mth_interrupt_clear(void) {
    g_mth_interrupt_pending = 0;
}

/* The depth is read only when the flag is set, so the common (no interrupt)
 * case costs one load. */
bool mth_interrupt_should_stop(void) {
    return g_mth_interrupt_pending != 0 && abort_protect_depth() == 0;
}
