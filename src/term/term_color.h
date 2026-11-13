/*
 * term_color.h -- colour for the interactive REPL.
 *
 * Pure C99 decision logic plus one cached, once-per-process query of the
 * environment. Colour is decided once, in main, for the interactive loop
 * only; pipe mode and -file never call term_color_init, so they stay plain
 * by construction. Rules, checked in order:
 *   not the interactive loop               -> off
 *   NO_COLOR set and non-empty             -> off
 *   CLICOLOR_FORCE set and not "0"         -> on, for both streams
 *   TERM unset, empty or "dumb"            -> off
 *   otherwise stdout colour = isatty(1), stderr colour = isatty(2).
 * Windows: always off.
 */
#ifndef MATHILDA_TERM_COLOR_H
#define MATHILDA_TERM_COLOR_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    const char* no_color;
    const char* clicolor_force;
    const char* term;
    bool interactive;
    bool stdout_tty;
    bool stderr_tty;
} TermColorEnv;

typedef struct { bool out, err; } TermColorDecision;

/* The rule table above, as a pure function of its inputs. */
TermColorDecision term_color_decide(const TermColorEnv* env);

/* Read getenv / isatty once and cache the decision. Before (or without) a
 * call, everything reports false. */
void term_color_init(bool interactive);
bool term_color_stdout(void);
bool term_color_stderr(void);
bool term_tty_stdout(void);          /* isatty(1), cached by term_color_init */

typedef enum {
    TERM_STYLE_IN,          /* In[n]:=          */
    TERM_STYLE_OUT,         /* Out[n]=          */
    TERM_STYLE_MSG_HEAD,    /* Head::tag        */
    TERM_STYLE_DIM,         /* the "took" note  */
    TERM_STYLE_COUNT
} TermStyle;

/* The SGR sequence for `s`, or "" when !enabled. */
const char* term_style_on(TermStyle s, bool enabled);
/* "\033[0m", or "" when !enabled. */
const char* term_style_off(bool enabled);

/* `text` styled into buf (snprintf semantics: returns the length it needed).
 * for_readline wraps every SGR in \001...\002 so readline does not count it
 * as prompt width. Leading spaces (Dialog depth) stay outside the style. */
int term_style_wrap(char* buf, size_t cap, TermStyle s, const char* text,
                    bool enabled, bool for_readline);

/* "1.23 s", "1 min 15 s", "1 h 2 min". */
int term_format_elapsed(char* buf, size_t cap, double seconds);

/* ---- tones for the live status line (REPL PR 3) --------------------------
 * Appended rather than added to TermStyle so PR 2 / PR 4 edits of the enum
 * above do not collide. Same contract: "" when !enabled. */
typedef enum {
    TERM_TONE_SPIN,     /* spinner glyph : cyan   */
    TERM_TONE_LABEL,    /* evaluated head: bold   */
    TERM_TONE_DIM,      /* times, memory : dim    */
    TERM_TONE_BAR,      /* progress      : cyan   */
    TERM_TONE_WARN,     /* memory >= 60% : yellow */
    TERM_TONE_ALERT,    /* memory >= 85% : red    */
    TERM_TONE_COUNT
} TermTone;
const char* term_tone_on(TermTone t, bool enabled);

#endif /* MATHILDA_TERM_COLOR_H */
