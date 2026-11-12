/*
 * term_theme.c -- see term_theme.h.
 */
#include "term_theme.h"
#include "term_config.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { A_BOLD = 1, A_DIM = 2, A_ITALIC = 4, A_UNDER = 8, A_REV = 16 };

/* One role: attributes, a 16-colour fg/bg SGR code (0 = none) and a 24-bit
 * fg/bg (-1 = none: the 16-colour code is used at every depth). */
typedef struct {
    unsigned char attr;
    unsigned char fg16, bg16;
    long fg, bg;
} RoleSpec;

#define S(a, f16, b16, f, b) { (a), (f16), (b16), (f), (b) }
#define N (-1L)

typedef struct {
    const char* name;
    const char* summary;
    RoleSpec role[TERM_ROLE_COUNT];
} Theme;

/* 16-colour codes of "dark" (and all of "ansi") are the REPL's original
 * colours, so a 16-colour terminal looks exactly as it did before themes. */
static const Theme THEMES[] = {
    { "dark", "the default, for dark backgrounds", {
        S(A_BOLD, 34, 0, 0x61AFEF, N),      /* In            */
        S(A_BOLD, 35, 0, 0xC678DD, N),      /* Out           */
        S(A_BOLD, 31, 0, 0xE06C75, N),      /* Head::tag     */
        S(A_DIM,   0, 0, N, N),             /* took note     */
        S(0,      36, 0, 0x56B6C2, N),      /* spinner glyph */
        S(A_BOLD,  0, 0, N, N),             /* status label  */
        S(A_DIM,   0, 0, N, N),             /* status text   */
        S(0,      36, 0, 0x56B6C2, N),      /* progress bar  */
        S(0,      33, 0, 0xE5C07B, N),      /* memory warn   */
        S(0,      31, 0, 0xE06C75, N),      /* memory alert  */
        S(0,      36, 0, 0x56B6C2, N),      /* builtin       */
        S(0,      34, 0, 0x61AFEF, N),      /* user symbol   */
        S(0,      35, 0, 0xD19A66, N),      /* number        */
        S(0,      32, 0, 0x98C379, N),      /* string        */
        S(0,      90, 0, 0x7F848E, N),      /* comment       */
        S(0,      33, 0, 0xE5C07B, N),      /* pattern       */
        S(0,       0, 0, N, N),             /* bracket       */
        S(A_BOLD | A_UNDER, 0, 0, N, N),    /* bracket match */
        S(A_BOLD, 31, 0, 0xE06C75, N),      /* stray closer  */
        S(0,      97, 100, 0xD7DAE0, 0x2C313A),          /* popup     */
        S(A_BOLD, 97, 44, 0xFFFFFF, 0x3D6FB4),           /* selection */
        S(0,      37, 100, 0x9DA5B4, 0x2C313A),          /* popup doc */
        S(A_DIM,   0, 0, N, N),             /* editor dim    */
    } },
    { "light", "for light backgrounds", {
        S(A_BOLD, 34, 0, 0x0550AE, N),
        S(A_BOLD, 35, 0, 0x8250DF, N),
        S(A_BOLD, 31, 0, 0xCF222E, N),
        S(0,      90, 0, 0x6E7781, N),
        S(0,      34, 0, 0x0969DA, N),
        S(A_BOLD,  0, 0, N, N),
        S(0,      90, 0, 0x6E7781, N),
        S(0,      34, 0, 0x0969DA, N),
        S(0,      33, 0, 0x9A6700, N),
        S(0,      31, 0, 0xCF222E, N),
        S(0,      34, 0, 0x0550AE, N),
        S(0,      35, 0, 0x8250DF, N),
        S(0,      31, 0, 0x953800, N),
        S(0,      32, 0, 0x116329, N),
        S(0,      90, 0, 0x6E7781, N),
        S(0,      33, 0, 0x9A6700, N),
        S(0,       0, 0, N, N),
        S(A_BOLD | A_UNDER, 0, 0, N, N),
        S(A_BOLD, 31, 0, 0xCF222E, N),
        S(0,      30, 47, 0x24292F, 0xEAEEF2),
        S(A_BOLD, 97, 44, 0xFFFFFF, 0x0969DA),
        S(0,      90, 47, 0x57606A, 0xEAEEF2),
        S(0,      90, 0, 0x6E7781, N),
    } },
    { "high-contrast", "bright, bold and never dim", {
        S(A_BOLD, 94, 0, 0x5FAFFF, N),
        S(A_BOLD, 95, 0, 0xFF87FF, N),
        S(A_BOLD, 91, 0, 0xFF5F5F, N),
        S(0,      37, 0, 0xD0D0D0, N),
        S(A_BOLD, 96, 0, 0x00FFFF, N),
        S(A_BOLD, 97, 0, 0xFFFFFF, N),
        S(0,      37, 0, 0xD0D0D0, N),
        S(A_BOLD, 96, 0, 0x00FFFF, N),
        S(A_BOLD, 93, 0, 0xFFFF00, N),
        S(A_BOLD, 91, 0, 0xFF5F5F, N),
        S(0,      96, 0, 0x00FFFF, N),
        S(0,      94, 0, 0x5FAFFF, N),
        S(0,      95, 0, 0xFF87FF, N),
        S(0,      92, 0, 0x5FFF5F, N),
        S(0,      37, 0, 0xD0D0D0, N),
        S(0,      93, 0, 0xFFFF00, N),
        S(A_BOLD,  0, 0, N, N),
        S(A_BOLD | A_UNDER | A_REV, 0, 0, N, N),
        S(A_BOLD | A_UNDER, 91, 0, 0xFF5F5F, N),
        S(0,      97, 40, 0xFFFFFF, 0x000000),
        S(A_BOLD, 30, 103, 0x000000, 0xFFFF00),
        S(0,      97, 40, 0xFFFFFF, 0x000000),
        S(0,      37, 0, 0xD0D0D0, N),
    } },
    { "mono", "no colour: bold, dim, underline and reverse only", {
        S(A_BOLD, 0, 0, N, N),
        S(A_BOLD, 0, 0, N, N),
        S(A_BOLD, 0, 0, N, N),
        S(A_DIM,  0, 0, N, N),
        S(0,      0, 0, N, N),
        S(A_BOLD, 0, 0, N, N),
        S(A_DIM,  0, 0, N, N),
        S(0,      0, 0, N, N),
        S(A_BOLD, 0, 0, N, N),
        S(A_BOLD | A_UNDER, 0, 0, N, N),
        S(0,      0, 0, N, N),
        S(0,      0, 0, N, N),
        S(0,      0, 0, N, N),
        S(0,      0, 0, N, N),
        S(A_DIM,  0, 0, N, N),
        S(0,      0, 0, N, N),
        S(0,      0, 0, N, N),
        S(A_BOLD | A_UNDER, 0, 0, N, N),
        S(A_REV,  0, 0, N, N),
        S(0,      0, 0, N, N),
        S(A_REV,  0, 0, N, N),
        S(A_DIM,  0, 0, N, N),
        S(A_DIM,  0, 0, N, N),
    } },
    { "ansi", "the terminal's own 16-colour palette", {
        S(A_BOLD, 34, 0, N, N),
        S(A_BOLD, 35, 0, N, N),
        S(A_BOLD, 31, 0, N, N),
        S(A_DIM,   0, 0, N, N),
        S(0,      36, 0, N, N),
        S(A_BOLD,  0, 0, N, N),
        S(A_DIM,   0, 0, N, N),
        S(0,      36, 0, N, N),
        S(0,      33, 0, N, N),
        S(0,      31, 0, N, N),
        S(0,      36, 0, N, N),
        S(0,      34, 0, N, N),
        S(0,      35, 0, N, N),
        S(0,      32, 0, N, N),
        S(0,      90, 0, N, N),
        S(0,      33, 0, N, N),
        S(0,       0, 0, N, N),
        S(A_BOLD | A_UNDER, 0, 0, N, N),
        S(A_BOLD, 31, 0, N, N),
        S(0,      97, 100, N, N),
        S(A_BOLD, 97, 44, N, N),
        S(0,      37, 100, N, N),
        S(A_DIM,   0, 0, N, N),
    } },
};

#undef S
#undef N

#define N_THEMES (sizeof(THEMES) / sizeof(THEMES[0]))

size_t term_theme_count(void) { return N_THEMES; }
const char* term_theme_name(size_t i) { return i < N_THEMES ? THEMES[i].name : NULL; }
const char* term_theme_summary(size_t i) { return i < N_THEMES ? THEMES[i].summary : NULL; }

int term_theme_find(const char* name) {
    if (!name) return -1;
    for (size_t i = 0; i < N_THEMES; i++) {
        const char* a = THEMES[i].name;
        const char* b = name;
        while (*a && *b && tolower((unsigned char)*a) == tolower((unsigned char)*b)) { a++; b++; }
        if (!*a && !*b) return (int)i;
    }
    return -1;
}

/* ---- depth --------------------------------------------------------------- */

static bool has(const char* s, const char* needle) { return s && strstr(s, needle) != NULL; }

TermDepth term_theme_depth_decide(const char* env_depth, const char* env_colorterm,
                                  const char* env_term) {
    if (env_depth && *env_depth) {
        if (strcmp(env_depth, "24") == 0 || strcmp(env_depth, "24bit") == 0 ||
            strcmp(env_depth, "truecolor") == 0) return TERM_DEPTH_24BIT;
        if (strcmp(env_depth, "256") == 0) return TERM_DEPTH_256;
        if (strcmp(env_depth, "16") == 0 || strcmp(env_depth, "8") == 0) return TERM_DEPTH_16;
        if (strcmp(env_depth, "0") == 0 || strcmp(env_depth, "none") == 0) return TERM_DEPTH_NONE;
    }
    if (!env_term || !*env_term || strcmp(env_term, "dumb") == 0) return TERM_DEPTH_NONE;
    if (env_colorterm && (strcmp(env_colorterm, "truecolor") == 0 ||
                          strcmp(env_colorterm, "24bit") == 0)) return TERM_DEPTH_24BIT;
    if (has(env_term, "-direct")) return TERM_DEPTH_24BIT;
    if (has(env_term, "256")) return TERM_DEPTH_256;
    return TERM_DEPTH_16;
}

/* ---- 256-colour approximation ------------------------------------------- */

static const int CUBE[6] = { 0, 95, 135, 175, 215, 255 };

static int cube_index(int v) {
    int best = 0;
    for (int i = 1; i < 6; i++)
        if (abs(CUBE[i] - v) < abs(CUBE[best] - v)) best = i;
    return best;
}

static long dist2(int r1, int g1, int b1, int r2, int g2, int b2) {
    long dr = r1 - r2, dg = g1 - g2, db = b1 - b2;
    return dr * dr * 2 + dg * dg * 4 + db * db * 3;   /* weighted for the eye */
}

int term_theme_rgb_to_256(unsigned long rgb) {
    int r = (int)((rgb >> 16) & 0xFF), g = (int)((rgb >> 8) & 0xFF), b = (int)(rgb & 0xFF);
    int ri = cube_index(r), gi = cube_index(g), bi = cube_index(b);
    long dc = dist2(r, g, b, CUBE[ri], CUBE[gi], CUBE[bi]);
    /* the grey ramp 232..255: 8, 18, ..., 238 */
    int avg = (r + g + b) / 3;
    int gs = avg < 8 ? 0 : (avg > 238 ? 23 : (avg - 8 + 5) / 10);
    if (gs > 23) gs = 23;
    int gv = 8 + 10 * gs;
    long dg = dist2(r, g, b, gv, gv, gv);
    return dg < dc ? 232 + gs : 16 + 36 * ri + 6 * gi + bi;
}

/* ---- building SGR -------------------------------------------------------- */

/* Append ";piece" (or "piece" first) to the parameter list. */
static void param(char* list, size_t cap, const char* piece) {
    size_t n = strlen(list);
    snprintf(list + n, cap - n, "%s%s", n ? ";" : "", piece);
}

static int build(const RoleSpec* r, TermDepth depth, char* buf, size_t cap) {
    if (!cap) return 0;
    buf[0] = '\0';
    if (depth == TERM_DEPTH_NONE) return 0;
    char list[64] = "";
    if (r->attr & A_BOLD)   param(list, sizeof list, "1");
    if (r->attr & A_DIM)    param(list, sizeof list, "2");
    if (r->attr & A_ITALIC) param(list, sizeof list, "3");
    if (r->attr & A_UNDER)  param(list, sizeof list, "4");
    if (r->attr & A_REV)    param(list, sizeof list, "7");
    for (int k = 0; k < 2; k++) {
        long rgb = k == 0 ? r->fg : r->bg;
        int c16 = k == 0 ? r->fg16 : r->bg16;
        char p[24] = "";
        if (rgb >= 0 && depth == TERM_DEPTH_24BIT)
            snprintf(p, sizeof p, "%d;2;%ld;%ld;%ld", k ? 48 : 38,
                     (rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF);
        else if (rgb >= 0 && depth == TERM_DEPTH_256)
            snprintf(p, sizeof p, "%d;5;%d", k ? 48 : 38, term_theme_rgb_to_256((unsigned long)rgb));
        else if (c16)
            snprintf(p, sizeof p, "%d", c16);
        if (p[0]) param(list, sizeof list, p);
    }
    if (!list[0]) return 0;                    /* a plain role */
    return snprintf(buf, cap, "\033[%sm", list);
}

int term_theme_sgr_for(size_t i, TermRole role, TermDepth depth, char* buf, size_t cap) {
    if (i >= N_THEMES) i = 0;
    if ((int)role < 0 || role >= TERM_ROLE_COUNT) { if (cap) buf[0] = '\0'; return 0; }
    return build(&THEMES[i].role[role], depth, buf, cap);
}

/* ---- the current theme --------------------------------------------------- */

#define SGR_CAP 48
static char      g_sgr[TERM_ROLE_COUNT][SGR_CAP];
static bool      g_applied = false;
static size_t    g_theme = 0;
static TermDepth g_depth = TERM_DEPTH_16;

void term_theme_apply(size_t i, TermDepth depth) {
    if (i >= N_THEMES) i = 0;
    for (int r = 0; r < TERM_ROLE_COUNT; r++)
        build(&THEMES[i].role[r], depth, g_sgr[r], SGR_CAP);
    g_theme = i;
    g_depth = depth;
    g_applied = true;
}

size_t    term_theme_current(void) { return g_theme; }
TermDepth term_theme_depth(void)   { return g_depth; }

const char* term_theme_sgr(TermRole role) {
    if (!g_applied) term_theme_apply(0, TERM_DEPTH_16);
    if ((int)role < 0 || role >= TERM_ROLE_COUNT) return "";
    return g_sgr[role];
}

void term_theme_init(const char* config_path) {
    TermDepth d = term_theme_depth_decide(getenv("MATHILDA_COLOR_DEPTH"), getenv("COLORTERM"),
                                          getenv("TERM"));
    /* NO_COLOR and friends are term_color's call; a terminal with colour on
     * but an unrecognised depth still gets the 16 basic colours. */
    if (d == TERM_DEPTH_NONE) d = TERM_DEPTH_16;
    int t = term_theme_find(getenv("MATHILDA_THEME"));
    if (t < 0) {
        char name[64];
        if (term_config_get(config_path, "theme", name, sizeof name)) t = term_theme_find(name);
    }
    if (t < 0 && term_theme_light_background(getenv("COLORFGBG"))) t = term_theme_find("light");
    term_theme_apply(t < 0 ? 0 : (size_t)t, d);
}

bool term_theme_light_background(const char* colorfgbg) {
    if (!colorfgbg || !*colorfgbg) return false;
    const char* bg = strrchr(colorfgbg, ';');
    bg = bg ? bg + 1 : colorfgbg;
    if (*bg < '0' || *bg > '9') return false;
    int v = atoi(bg);
    return v == 7 || v == 15;
}
