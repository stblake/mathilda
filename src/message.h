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

#include <stdarg.h>

/* printf-style format checking for the emit funnel below. GUARDED: CLAUDE.md
 * bans an unguarded __attribute__ (breaks strict C99 on a non-GCC compiler and
 * trips make check-c99); it is a no-op everywhere but GCC/Clang. */
#if defined(__GNUC__)
#  define MTH_PRINTF_FMT(a, b) __attribute__((format(printf, a, b)))
#else
#  define MTH_PRINTF_FMT(a, b)
#endif

/* Enter a quiet region that silences every message, as Quiet[expr] does
 * (push a Quiet frame; see the frame model in message.c).  A message silenced
 * here is also hidden from any Check[] entered outside the region. */
void mth_msg_suppress_push(void);
/* Leave a quiet region (pop a frame; saturates at zero). */
void mth_msg_suppress_pop(void);
/* Non-zero iff a quiet region silencing every message is active. */
int  mth_msg_suppressed(void);
/* Save / restore the Quiet and Check frame depths as one value (used by
 * TimeConstrained to undo regions abandoned by a timeout siglongjmp). */
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
/* Save/restore the ifun depth across a TimeConstrained unwind (see the Quiet
 * pair above); builtin_reduce's push/pop can be skipped by a timeout siglongjmp. */
int  mth_msg_ifun_suppress_depth_save(void);
void mth_msg_ifun_suppress_depth_load(int d);

/*
 * Message-fired counter.  Distinct from the suppression depth above: a message
 * still "fires" (and bumps this counter) even while suppressed, so internal
 * code can snapshot mth_msg_fired_count() around a call to learn whether it
 * reported anything.  An emission site that does not route through the funnel
 * calls mth_msg_note_fired(): that counts it and marks the enclosing Check[]
 * frames that can see an anonymous diagnostic (those inside the innermost
 * Quiet that silences everything).  Check itself does not use the counter: it
 * keeps a frame of its own, so a message silenced by a Quiet INSIDE the Check
 * does not count (Check[Quiet[1/0]; 5, 7] is 5), while one silenced by a Quiet
 * OUTSIDE it does (Quiet[Check[1/0; 5, 7]] is 7).
 */
void          mth_msg_note_fired(void);
/* Note head::tag as fired without printing it: for a site that has already
 * printed it once (a de-duplicated repeat) or skips the formatting work under
 * a Quiet that covers everything.  Quiet / Check frames see it by name. */
void          mth_msg_note_named(const char* head, const char* tag);
unsigned long mth_msg_fired_count(void);
/* Roll the counter back to an earlier snapshot.  For an internal PROBE region
 * whose diagnostics are sampling noise (a speculative sub-solve whose result the
 * caller may discard): bracket it with mth_msg_suppress_push/pop and restore the
 * pre-region count, so the noise is neither printed nor seen by Check[] -- the
 * region-wide analogue of mth_message_gated's extra_mute.  Only ever lowers. */
void          mth_msg_fired_restore(unsigned long n);

/* ------------------------------------------------------------- The funnel */
/*
 * mth_message is the ONE place the two-step Quiet/Check contract lives: it
 * notes the firing (so every enclosing Check[] entered inside the Quiet that
 * silences it -- all of them when none does -- sees it) and then prints
 * "Head::tag: <formatted>\n" to stderr UNLESS a Quiet[] covering it is active.  Every
 * user-facing diagnostic routes through it (or through a subsystem helper that
 * delegates to it); a raw fprintf(stderr, "Head::tag: ...") bypasses both
 * Quiet[] and Check[] and is a bug.  `head` and `tag` are separate arguments,
 * so no "::"-bearing literal ever reaches an fprintf -- which is what lets
 * make check-messages tell a routed site from a bypassing one.
 */
void mth_message(const char* head, const char* tag, const char* fmt, ...)
    MTH_PRINTF_FMT(3, 4);

/*
 * As mth_message, but with an EXTRA subsystem-local *probe* mute -- pass the RAW
 * local flag (g_arith_warnings_muted for Power/Plus/Times, g_fm_quiet for
 * FindMinimum), NOT `flag || suppressed`.  When that flag is set an internal
 * probe is poking at divergent forms and the diagnostic is sampling noise, not a
 * user-facing event, so it is suppressed ENTIRELY: neither printed NOR noted, so
 * an enclosing Check[] does not catch the noise.  (Quiet, by contrast, silences
 * the print and hides the message only from the Checks outside the Quiet.)
 */
void mth_message_gated(int extra_mute, const char* head, const char* tag,
                       const char* fmt, ...) MTH_PRINTF_FMT(4, 5);

/* va_list core, exported so subsystem helpers (fit_warn, fm_warn, ...) delegate
 * to the funnel rather than re-implement note+guard.  extra_mute as above. */
void mth_message_v(int extra_mute, const char* head, const char* tag,
                   const char* fmt, va_list ap);

/* A continuation line under a diagnostic already emitted this call: prints an
 * extra line with NO "Head::tag:" prefix, honours the same mute, and does NOT
 * note again (the diagnostic already fired).  For the few messages a single
 * fmt cannot express (built across a loop). */
void mth_message_cont(int extra_mute, const char* fmt, ...) MTH_PRINTF_FMT(2, 3);

/* Wrap the printed "Head::tag" of every message in head_on ... head_off (the
 * interactive REPL's colour). The strings are borrowed and must be static;
 * NULL or "" (the default) leaves the output unchanged. */
void mth_message_set_style(const char* head_on, const char* head_off);

/* Registers the Quiet / Check / Message builtins.  Called from core_init(). */
void message_init(void);

/* Off[head::tag] / Off[General::tag] / Off[head]: non-zero iff head::tag is
 * switched off.  The funnel drops such a message entirely. */
int mth_msg_is_off(const char* head, const char* tag);

/* Reset $MessageList to {} -- called by the front ends at each new input line. */
void mth_message_list_reset(void);

/* Emit a preformatted "Head::tag: text" line through the funnel (with the same
 * extra_mute meaning as mth_message_gated).  For arith_warn's fixed strings. */
void mth_message_line(int extra_mute, const char* line);
/* printf-style mth_message_line: formats "Head::tag: text" (any length) and
 * emits it through the funnel by name.  For legacy helpers whose format
 * strings carry the "Head::tag: " prefix. */
void mth_message_linev(int extra_mute, const char* fmt, va_list ap);

#endif /* MATHILDA_MESSAGE_H */
