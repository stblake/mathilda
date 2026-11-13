/*
 * term_color.c -- see term_color.h.
 */
/* isatty()/fileno() are POSIX, hidden by glibc under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "term_color.h"
#include "term_theme.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <unistd.h>
#endif

/* The SGR of every style and tone comes from the current theme
 * (term_theme.h); before a theme is chosen these are the original colours. */
static const char STYLE_RESET[] = "\033[0m";

TermColorDecision term_color_decide(const TermColorEnv* e) {
    TermColorDecision d = { false, false };
#ifdef _WIN32
    (void)e;
    return d;
#else
    if (!e || !e->interactive) return d;
    if (e->no_color && *e->no_color) return d;
    if (e->clicolor_force && *e->clicolor_force && strcmp(e->clicolor_force, "0") != 0) {
        d.out = d.err = true;
        return d;
    }
    if (!e->term || !*e->term || strcmp(e->term, "dumb") == 0) return d;
    d.out = e->stdout_tty;
    d.err = e->stderr_tty;
    return d;
#endif
}

static TermColorDecision g_color = { false, false };
static bool g_tty_out = false;

void term_color_init(bool interactive) {
    TermColorEnv e;
    e.no_color = getenv("NO_COLOR");
    e.clicolor_force = getenv("CLICOLOR_FORCE");
    e.term = getenv("TERM");
    e.interactive = interactive;
#ifndef _WIN32
    e.stdout_tty = isatty(fileno(stdout)) != 0;
    e.stderr_tty = isatty(fileno(stderr)) != 0;
#else
    e.stdout_tty = e.stderr_tty = false;
#endif
    g_tty_out = interactive && e.stdout_tty;
    g_color = term_color_decide(&e);
}

bool term_color_stdout(void) { return g_color.out; }
bool term_color_stderr(void) { return g_color.err; }
bool term_tty_stdout(void)   { return g_tty_out; }

const char* term_style_on(TermStyle s, bool enabled) {
    if (!enabled || (int)s < 0 || (int)s >= (int)TERM_STYLE_COUNT) return "";
    return term_theme_sgr((TermRole)(TERM_ROLE_STYLE_BASE + (int)s));
}

const char* term_style_off(bool enabled) {
    return enabled ? STYLE_RESET : "";
}

int term_style_wrap(char* buf, size_t cap, TermStyle s, const char* text,
                    bool enabled, bool for_readline) {
    if (!enabled) return snprintf(buf, cap, "%s", text);
    size_t lead = 0;
    while (text[lead] == ' ') lead++;
    const char* o1 = for_readline ? "\001" : "";
    const char* o2 = for_readline ? "\002" : "";
    return snprintf(buf, cap, "%.*s%s%s%s%s%s%s%s", (int)lead, text,
                    o1, term_style_on(s, true), o2, text + lead,
                    o1, STYLE_RESET, o2);
}

int term_format_elapsed(char* buf, size_t cap, double seconds) {
    if (seconds < 0) seconds = 0;
    if (seconds < 60.0) return snprintf(buf, cap, "%.2f s", seconds);
    long t = (long)seconds;
    if (t < 3600) return snprintf(buf, cap, "%ld min %ld s", t / 60, t % 60);
    return snprintf(buf, cap, "%ld h %ld min", t / 3600, (t % 3600) / 60);
}

/* ---- tones for the live status line (REPL PR 3) -------------------------- */

const char* term_tone_on(TermTone t, bool enabled) {
    if (!enabled || (int)t < 0 || (int)t >= (int)TERM_TONE_COUNT) return "";
    return term_theme_sgr((TermRole)(TERM_ROLE_TONE_BASE + (int)t));
}
