/*
 * term_spinner.c -- see term_spinner.h.
 */
#include "term_spinner.h"
#include "term_color.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* ---- small formatters ---------------------------------------------------- */

int term_fmt_bytes(char* buf, size_t cap, uint64_t b) {
    const double k = 1024.0;
    double v = (double)b;
    if (v < k * k)           return snprintf(buf, cap, "%.0f KB", v / k);
    if (v < k * k * k)       return snprintf(buf, cap, "%.0f MB", v / (k * k));
    if (v < 10.0 * k * k * k) return snprintf(buf, cap, "%.2f GB", v / (k * k * k));
    return snprintf(buf, cap, "%.1f GB", v / (k * k * k));
}

int term_fmt_secs(char* buf, size_t cap, double s) {
    if (!(s >= 0)) s = 0;
    if (s < 60.0) return snprintf(buf, cap, "%.1fs", s);
    long t = (long)s;
    if (t < 3600) return snprintf(buf, cap, "%ldm%02lds", t / 60, t % 60);
    return snprintf(buf, cap, "%ldh%02ldm", t / 3600, (t % 3600) / 60);
}

int term_fmt_count(char* buf, size_t cap, long n) {
    if (n < 0) n = 0;
    if (n < 1000000L) return snprintf(buf, cap, "%ld", n);
    /* 1.25M, 22.4M, 1G: two decimals, trailing zeros dropped */
    char t[48];
    double v = n < 1000000000L ? (double)n / 1e6 : (double)n / 1e9;
    snprintf(t, sizeof t, "%.2f", floor(v * 100.0) / 100.0);
    size_t l = strlen(t);
    while (l > 0 && t[l - 1] == '0') t[--l] = '\0';
    if (l > 0 && t[l - 1] == '.') t[--l] = '\0';
    return snprintf(buf, cap, "%s%s", t, n < 1000000000L ? "M" : "G");
}

int term_fmt_bar(char* buf, size_t cap, double frac, int cells, bool unicode) {
    if (cells < 1) cells = 1;
    if (!(frac >= 0)) frac = 0;
    if (frac > 1) frac = 1;
    size_t o = 0;
    char tmp[256];
    if (cap == 0) return 0;
#define PUT(str) do { const char* p_ = (str); size_t l_ = strlen(p_); \
        if (o + l_ < sizeof tmp) { memcpy(tmp + o, p_, l_); o += l_; } } while (0)
    if (unicode) {
        /* eighths: U+258F..U+2589, from one eighth to seven eighths */
        static const char* const part[8] = { "", "\xe2\x96\x8f", "\xe2\x96\x8e", "\xe2\x96\x8d",
                                             "\xe2\x96\x8c", "\xe2\x96\x8b", "\xe2\x96\x8a", "\xe2\x96\x89" };
        int eighths = (int)floor(frac * cells * 8 + 1e-9);
        int full = eighths / 8, rem = eighths % 8;
        PUT("\xe2\x96\x95");                       /* U+2595 right one eighth: left edge */
        for (int i = 0; i < full; i++) PUT("\xe2\x96\x88");   /* full block */
        int used = full;
        if (rem && used < cells) { PUT(part[rem]); used++; }
        for (; used < cells; used++) PUT(" ");
        PUT("\xe2\x96\x8f");                       /* U+258F left one eighth: right edge */
    } else {
        int k = (int)floor(frac * cells + 1e-9);
        PUT("[");
        for (int i = 0; i < cells; i++) PUT(i < k ? "#" : ".");
        PUT("]");
    }
#undef PUT
    tmp[o] = '\0';
    return snprintf(buf, cap, "%s", tmp);
}

int term_display_width(const char* s) {
    int w = 0;
    for (const unsigned char* p = (const unsigned char*)s; *p; p++) {
        if (*p == 0x1b) {                       /* CSI ... final byte */
            if (p[1] == '[') {
                p += 2;
                while (*p && !(*p >= 0x40 && *p <= 0x7e)) p++;
                if (!*p) break;
            }
            continue;
        }
        if ((*p & 0xc0) != 0x80) w++;           /* not a continuation byte */
    }
    return w;
}

/* ---- fraction and ETA ------------------------------------------------------ */

double term_spinner_fraction(const TermSpinnerLoop* o, const TermSpinnerLoop* in, int depth) {
    if (!o || !o->head || o->total <= 0) return -1;
    double f = (double)o->done;
    if (depth >= 2 && in && in->total > 0) {
        double g = (double)in->done / (double)in->total;
        if (g > 1) g = 1;
        if (g > 0) f += g;
    }
    f /= (double)o->total;
    if (f < 0) f = 0;
    if (f > 1) f = 1;
    return f;
}

double term_eta_update(TermEta* e, unsigned long gen, double frac, double now) {
    if (frac < 0) { e->have = false; return -1; }
    if (!e->have || e->gen != gen || frac < e->f0) {
        e->have = true; e->gen = gen; e->t0 = now; e->f0 = frac;
        return -1;
    }
    double dt = now - e->t0, df = frac - e->f0;
    if (dt < 1.0 || df <= 0) return -1;
    double eta = dt / df * (1.0 - frac);
    return eta < 0 ? 0 : eta;
}

/* ---- the live line --------------------------------------------------------- */

enum { SEG_TEXT = 160 };

typedef struct {
    char text[SEG_TEXT];   /* plain text, no escapes */
    const char* sgr;       /* style for the whole segment ("" = none) */
    int prio;              /* drop order: lowest goes first; 99 = never */
    bool sep_before;       /* " · " before it */
    bool on;
} Seg;

static void seg(Seg* s, int prio, bool sep, const char* sgr, const char* fmt_text) {
    snprintf(s->text, sizeof s->text, "%s", fmt_text);
    s->sgr = sgr ? sgr : "";
    s->prio = prio;
    s->sep_before = sep;
    s->on = fmt_text[0] != '\0';
}

/* Cut s to at most `cells` display cells (UTF-8 safe). */
static void cut_cells(char* s, int cells) {
    int w = 0;
    unsigned char* p = (unsigned char*)s;
    for (; *p; p++) {
        if ((*p & 0xc0) != 0x80) {
            if (w == cells) break;
            w++;
        }
    }
    *p = '\0';
}

static double pct_of(uint64_t v, uint64_t budget) {
    return budget ? 100.0 * (double)v / (double)budget : -1;
}

int term_spinner_render(const TermSpinnerSnap* s, int cols, TermSpinnerStyle st,
                        char* buf, size_t cap) {
    if (!buf || cap == 0) return 0;
    buf[0] = '\0';
    if (!s) return 0;
    int room = cols - 1;
    if (room < 1) return 0;

    const char* sep = st.unicode ? " \xc2\xb7 " : " | ";
    const int sepw = 3;
    const char* dim = term_tone_on(TERM_TONE_DIM, st.color);
    Seg g[10];
    int n = 0;
    char t[SEG_TEXT], a[48], b[48];

    /* 0: glyph + label (never dropped; the label may be cut) */
    static const char* const BRAILLE[10] = {
        "\xe2\xa0\x8b", "\xe2\xa0\x99", "\xe2\xa0\xb9", "\xe2\xa0\xb8", "\xe2\xa0\xbc",
        "\xe2\xa0\xb4", "\xe2\xa0\xa6", "\xe2\xa0\xa7", "\xe2\xa0\x87", "\xe2\xa0\x8f" };
    static const char ASCII[4] = { '|', '/', '-', '\\' };
    char glyph[8];
    if (st.unicode) snprintf(glyph, sizeof glyph, "%s", BRAILLE[s->frame % 10]);
    else            snprintf(glyph, sizeof glyph, "%c", ASCII[s->frame % 4]);
    const char* label = (s->label && *s->label) ? s->label : "Evaluating";

    /* 2: progress: "Table 3400/10000 ▕███▌    ▏ 34%" */
    t[0] = '\0';
    if (s->outer.head) {
        size_t o = 0;
        term_fmt_count(a, sizeof a, s->outer.done);
        if (s->outer.total > 0) {
            term_fmt_count(b, sizeof b, s->outer.total);
            o += (size_t)snprintf(t + o, sizeof t - o, "%s %s/%s", s->outer.head, a, b);
        } else {
            o += (size_t)snprintf(t + o, sizeof t - o, "%s %s", s->outer.head, a);
        }
        if (s->inner.head && o < sizeof t) {
            term_fmt_count(a, sizeof a, s->inner.done);
            const char* arrow = st.unicode ? " \xe2\x96\xb8 " : " > ";
            if (s->inner.total > 0) {
                term_fmt_count(b, sizeof b, s->inner.total);
                o += (size_t)snprintf(t + o, sizeof t - o, "%s%s %s/%s", arrow, s->inner.head, a, b);
            } else {
                o += (size_t)snprintf(t + o, sizeof t - o, "%s%s %s", arrow, s->inner.head, a);
            }
        }
        if (s->frac >= 0 && o < sizeof t) {
            char bar[96];
            term_fmt_bar(bar, sizeof bar, s->frac, 10, st.unicode);
            snprintf(t + o, sizeof t - o, " %s %2d%%", bar, (int)floor(s->frac * 100 + 1e-9));
        }
    }
    /* A loop that is the input itself ("Do[...]") names it once: the
     * progress text takes the label's place and is never dropped. */
    char lab[SEG_TEXT];
    const char* lab_tone = term_tone_on(TERM_TONE_LABEL, st.color);
    if (t[0] && s->outer.head && strcmp(label, s->outer.head) == 0) {
        snprintf(lab, sizeof lab, "%s", t);
        lab_tone = term_tone_on(TERM_TONE_BAR, st.color);
        t[0] = '\0';
    } else {
        snprintf(lab, sizeof lab, "%s", label);
    }
    seg(&g[n++], 60, true, term_tone_on(TERM_TONE_BAR, st.color), t);

    /* Monitor[expr, mon]: what the user asked to watch, kept until only the
     * elapsed time and the label are left; cut to half the line first. */
    t[0] = '\0';
    if (s->monitor && *s->monitor) {
        snprintf(t, sizeof t, "%s", s->monitor);
        int keep = room / 2 > 16 ? room / 2 : 16;
        if (term_display_width(t) > keep) cut_cells(t, keep);
    }
    seg(&g[n++], 80, true, term_tone_on(TERM_TONE_LABEL, st.color), t);

    /* the builtin running now, when it adds something */
    const char* live = s->live;
    if (live && (!*live || strcmp(live, label) == 0 ||
                 (s->outer.head && strcmp(live, s->outer.head) == 0) ||
                 (s->inner.head && strcmp(live, s->inner.head) == 0)))
        live = NULL;
    seg(&g[n++], 10, true, dim, live ? live : "");

    /* 3: eta */
    t[0] = '\0';
    if (s->outer.head && s->eta >= 0) {
        term_fmt_secs(a, sizeof a, s->eta);
        snprintf(t, sizeof t, "eta %s", a);
    }
    seg(&g[n++], 40, true, dim, t);

    /* 4: elapsed (kept longest after the label) */
    term_fmt_secs(a, sizeof a, s->elapsed);
    seg(&g[n++], 90, true, dim, a);

    /* 5: cpu */
    t[0] = '\0';
    if (s->cpu >= 0) { term_fmt_secs(a, sizeof a, s->cpu); snprintf(t, sizeof t, "cpu %s", a); }
    seg(&g[n++], 20, true, dim, t);

    /* 6: mem (share of the budget, coloured under pressure) */
    t[0] = '\0';
    const char* memtone = dim;
    if (s->mem) {
        term_fmt_bytes(a, sizeof a, s->mem);
        double p = pct_of(s->mem, s->budget);
        if (p >= 0 && p < 0.1) snprintf(t, sizeof t, "mem %s (<0.1%%)", a);
        else if (p >= 0)       snprintf(t, sizeof t, "mem %s (%.*f%%)", a, p < 10 ? 1 : 0, p);
        else        snprintf(t, sizeof t, "mem %s", a);
        if (p >= 85)      memtone = term_tone_on(TERM_TONE_ALERT, st.color);
        else if (p >= 60) memtone = term_tone_on(TERM_TONE_WARN, st.color);
    }
    seg(&g[n++], 50, true, memtone, t);

    /* 7: peak */
    t[0] = '\0';
    const char* peaktone = dim;
    if (s->peak && s->peak > s->mem) {
        term_fmt_bytes(a, sizeof a, s->peak);
        snprintf(t, sizeof t, "peak %s", a);
        double p = pct_of(s->peak, s->budget);
        if (p >= 85)      peaktone = term_tone_on(TERM_TONE_ALERT, st.color);
        else if (p >= 60) peaktone = term_tone_on(TERM_TONE_WARN, st.color);
    }
    seg(&g[n++], 15, true, peaktone, t);

    /* 8: worker threads */
    t[0] = '\0';
    if (s->workers > 0) snprintf(t, sizeof t, "%d thread%s", s->workers, s->workers == 1 ? "" : "s");
    seg(&g[n++], 30, true, dim, t);

    /* Fit: the head part is glyph + " " + label. */
    int head_w = term_display_width(glyph) + 1 + term_display_width(lab);
    for (;;) {
        int w = head_w;
        for (int i = 0; i < n; i++)
            if (g[i].on) w += (g[i].sep_before ? sepw : 0) + term_display_width(g[i].text);
        if (w <= room) break;
        int lo = -1;
        for (int i = 0; i < n; i++)
            if (g[i].on && (lo < 0 || g[i].prio < g[lo].prio)) lo = i;
        if (lo < 0) break;
        g[lo].on = false;
    }
    int glyph_w = term_display_width(glyph);
    if (glyph_w + 1 + term_display_width(lab) > room) {
        int keep = room - glyph_w - 1;
        if (keep < 0) keep = 0;
        cut_cells(lab, keep);
    }

    /* Assemble. */
    size_t o = 0;
    const char* off = term_style_off(st.color);
#define APP(...) do { \
        if (o < cap) { \
            int k_ = snprintf(buf + o, cap - o, __VA_ARGS__); \
            if (k_ > 0) o += (size_t)k_; \
            if (o >= cap) o = cap - 1; \
        } } while (0)
    int width = 0;
    if (glyph_w <= room) {
        APP("%s%s%s", term_tone_on(TERM_TONE_SPIN, st.color), glyph, off);
        width += glyph_w;
        if (width + 1 + term_display_width(lab) <= room && lab[0]) {
            APP(" %s%s%s", lab_tone, lab, off);
            width += 1 + term_display_width(lab);
        }
    }
    for (int i = 0; i < n; i++) {
        if (!g[i].on) continue;
        int w = term_display_width(g[i].text) + (g[i].sep_before ? sepw : 0);
        if (width + w > room) continue;
        if (g[i].sep_before) APP("%s%s%s", dim, sep, *dim ? off : "");
        APP("%s%s%s", g[i].sgr, g[i].text, *g[i].sgr ? off : "");
        width += w;
    }
#undef APP
    return width;
}

/* ---- the took note and the status line ------------------------------------- */

int term_took_line(char* buf, size_t cap, double elapsed, double cpu, uint64_t peak,
                   TermSpinnerStyle st) {
    char el[64], c[96] = "", pk[96] = "", tmp[64];
    const char* sep = st.unicode ? " \xc2\xb7 " : " | ";
    term_format_elapsed(el, sizeof el, elapsed);
    if (cpu >= 0) { term_format_elapsed(tmp, sizeof tmp, cpu); snprintf(c, sizeof c, "%scpu %s", sep, tmp); }
    if (peak)     { term_fmt_bytes(tmp, sizeof tmp, peak);     snprintf(pk, sizeof pk, "%speak %s", sep, tmp); }
    return snprintf(buf, cap, "%stook %s%s%s%s\n", term_style_on(TERM_STYLE_DIM, st.color),
                    el, c, pk, term_style_off(st.color));
}

int term_status_line(char* buf, size_t cap, double last, uint64_t mem,
                     const char* context, TermSpinnerStyle st) {
    char el[64], m[96] = "", tmp[64];
    const char* sep = st.unicode ? " \xc2\xb7 " : " | ";
    term_format_elapsed(el, sizeof el, last);
    if (mem) { term_fmt_bytes(tmp, sizeof tmp, mem); snprintf(m, sizeof m, "%smem %s", sep, tmp); }
    return snprintf(buf, cap, "%slast %s%s%s%s%s\n", term_style_on(TERM_STYLE_DIM, st.color),
                    el, m, context && *context ? sep : "", context ? context : "",
                    term_style_off(st.color));
}
