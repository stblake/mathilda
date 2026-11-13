/*
 * term_theme.h -- colour themes for the interactive REPL.
 *
 * Every colour the session draws comes from here: the In/Out prompts and
 * message heads (term_color), the live status line's tones (term_spinner /
 * term_status), and the editor's syntax, bracket and popup styles
 * (term_highlight). A theme gives each of those roles a set of attributes and
 * a colour in two forms: a 24-bit RGB value and a 16-colour SGR code. The
 * role's SGR sequence is built once per theme change for the terminal's
 * colour depth:
 *
 *   24-bit   38;2;R;G;B / 48;2;R;G;B
 *   256      38;5;N / 48;5;N, the nearest xterm-256 cube or grey entry
 *   16       the theme's 16-colour code
 *   none     ""  (term_color decides whether colour is on at all: NO_COLOR,
 *                 TERM=dumb and non-terminals always win)
 *
 * The depth comes from MATHILDA_COLOR_DEPTH (0, 16, 256, 24), else COLORTERM
 * (truecolor / 24bit), else TERM (*-direct: 24, *256color*: 256, set: 16).
 *
 * Shipped themes: "dark" (the default), "light", "high-contrast", "mono"
 * (attributes only) and "ansi" (the terminal's own 16-colour palette; also
 * what "dark" is below 24-bit colour, so a light 256-colour terminal still
 * gets colours that suit its background). The choice is read from
 * MATHILDA_THEME, else the "theme" key of the config file (term_config.h),
 * and saved there by /theme.
 *
 * Pure C99. Not thread-safe to change: the REPL changes the theme only at a
 * prompt, while the status thread is idle.
 */
#ifndef MATHILDA_TERM_THEME_H
#define MATHILDA_TERM_THEME_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    /* term_color's TermStyle, in its order */
    TERM_ROLE_IN, TERM_ROLE_OUT, TERM_ROLE_MSG_HEAD, TERM_ROLE_DIM,
    /* term_color's TermTone, in its order */
    TERM_ROLE_SPIN, TERM_ROLE_LABEL, TERM_ROLE_TONE_DIM, TERM_ROLE_BAR,
    TERM_ROLE_WARN, TERM_ROLE_ALERT,
    /* term_highlight's lexer and editor styles */
    TERM_ROLE_BUILTIN, TERM_ROLE_USER, TERM_ROLE_NUMBER, TERM_ROLE_STRING,
    TERM_ROLE_COMMENT, TERM_ROLE_PATTERN, TERM_ROLE_BRACKET,
    TERM_ROLE_MATCH, TERM_ROLE_BAD,
    TERM_ROLE_POPUP, TERM_ROLE_POPUP_SEL, TERM_ROLE_POPUP_DOC, TERM_ROLE_HL_DIM,
    TERM_ROLE_COUNT
} TermRole;

/* The first role of each group, so term_color can index by its own enums. */
#define TERM_ROLE_STYLE_BASE TERM_ROLE_IN
#define TERM_ROLE_TONE_BASE  TERM_ROLE_SPIN

typedef enum { TERM_DEPTH_NONE = 0, TERM_DEPTH_16 = 16, TERM_DEPTH_256 = 256,
               TERM_DEPTH_24BIT = 24 } TermDepth;

/* The colour depth from the environment (rules above). Pure. */
TermDepth   term_theme_depth_decide(const char* env_depth, const char* env_colorterm,
                                    const char* env_term);

/* The shipped themes, in menu order. */
size_t      term_theme_count(void);
const char* term_theme_name(size_t i);
const char* term_theme_summary(size_t i);    /* one line for the picker     */
/* The index of the theme called `name` (case-insensitive), or -1. */
int         term_theme_find(const char* name);

/* Build every role's SGR for theme `i` at `depth` (an unknown index selects
 * the default). Until the first call the roles hold "dark" at 16 colours,
 * which are the REPL's original colours. */
void        term_theme_apply(size_t i, TermDepth depth);
size_t      term_theme_current(void);
TermDepth   term_theme_depth(void);

/* The SGR sequence that selects `role` from a reset state; "" for a plain
 * role or an index out of range. The pointer stays valid for the process
 * and its contents follow term_theme_apply. */
const char* term_theme_sgr(TermRole role);

/* For tests: the SGR theme `i` gives `role` at `depth`, into buf. */
int         term_theme_sgr_for(size_t i, TermRole role, TermDepth depth, char* buf, size_t cap);
/* For tests: the xterm-256 index nearest to 0xRRGGBB. */
int         term_theme_rgb_to_256(unsigned long rgb);

/* The session's start-up choice: MATHILDA_THEME, else the config file's
 * "theme" key, else "light" when COLORFGBG says the background is light,
 * else "dark"; the depth from the environment. Applies it. */
void        term_theme_init(const char* config_path);

/* COLORFGBG ("fg;bg", set by rxvt, Konsole and some others) names a light
 * background (7 or 15). Pure. */
bool        term_theme_light_background(const char* colorfgbg);

#endif /* MATHILDA_TERM_THEME_H */
