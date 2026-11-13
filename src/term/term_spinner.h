/*
 * term_spinner.h -- formatting of the REPL's live status line, the "took"
 * note and the between-evaluations status line.
 *
 * Pure C99: every function here is a function of its arguments (no clock, no
 * terminal, no globals), so the whole layout is unit-tested without a tty.
 * term_status.c samples the clocks, memory and the progress stack and hands a
 * snapshot to term_spinner_render ten times a second.
 *
 * Layout, left to right (segments separated by " · ", or " | " in ASCII):
 *
 *   <glyph> <label> [<live head>] · <progress> [· <monitor>] [· eta 6s] · 3.2s
 *     · cpu 3.1s · mem 412 MB (2.5%) · peak 1.2 GB · 4 threads
 *
 * When the line does not fit in `cols - 1` cells (the last column is never
 * written, so the terminal cannot auto-wrap), whole segments are dropped in
 * this order: live head, peak, cpu, threads, eta, mem, progress, monitor,
 * elapsed; then the label is cut. The monitor text is first cut to half the
 * line. With colour, the glyph and bar are cyan, the label bold,
 * the rest dim, and memory turns yellow at 60 % of the budget, red at 85 %.
 * Without colour no escape sequence is produced at all. Without Unicode the
 * glyph is |/-\ and the bar [#####.....].
 */
#ifndef MATHILDA_TERM_SPINNER_H
#define MATHILDA_TERM_SPINNER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    bool color;     /* SGR allowed                         */
    bool unicode;   /* braille glyphs, block bar, "·"      */
} TermSpinnerStyle;

typedef struct {
    const char* head;   /* NULL: no loop */
    long done, total;   /* total 0 = unknown */
} TermSpinnerLoop;

typedef struct {
    unsigned frame;         /* advances once per redraw            */
    const char* label;      /* the input's head ("Integrate"); may be NULL */
    const char* live;       /* the builtin running now; may be NULL */
    double elapsed;         /* seconds                             */
    double cpu;             /* seconds; < 0 unknown                */
    TermSpinnerLoop outer;  /* outermost loop                      */
    TermSpinnerLoop inner;  /* the loop directly inside it         */
    double frac;            /* overall fraction 0..1; < 0 unknown  */
    double eta;             /* seconds; < 0 unknown                */
    uint64_t mem;           /* resident bytes; 0 unknown           */
    uint64_t peak;          /* peak resident bytes; 0 unknown      */
    uint64_t budget;        /* RAM or memory limit; 0 unknown      */
    int workers;            /* worker threads running              */
    const char* monitor;    /* Monitor[expr, mon]'s text; may be NULL */
} TermSpinnerSnap;

/* Render the line into buf (NUL-terminated, truncated to cap). Returns the
 * number of terminal cells it occupies (<= cols - 1). */
int term_spinner_render(const TermSpinnerSnap* s, int cols, TermSpinnerStyle st,
                        char* buf, size_t cap);

/* Fraction 0..1 of a loop nest, or -1 when the outer total is unknown. The
 * inner loop's share is added only when its total is known. */
double term_spinner_fraction(const TermSpinnerLoop* outer, const TermSpinnerLoop* inner,
                             int depth);

/* ETA tracker: reset whenever `gen` (the outer loop's identity) changes; an
 * estimate is given after one second of observed progress. */
typedef struct { unsigned long gen; double t0, f0; bool have; } TermEta;
double term_eta_update(TermEta* e, unsigned long gen, double frac, double now);

/* The pieces, exposed for tests and for the took / status lines. */
int term_fmt_bytes(char* buf, size_t cap, uint64_t bytes);       /* "412 MB", "1.4 GB" */
int term_fmt_secs(char* buf, size_t cap, double s);              /* "3.2s", "1m05s"    */
int term_fmt_count(char* buf, size_t cap, long n);               /* "3400", "1.25M"    */
int term_fmt_bar(char* buf, size_t cap, double frac, int cells, bool unicode);
int term_display_width(const char* s);   /* UTF-8 code points, escapes skipped */

/* "took 3.21 s · cpu 3.18 s · peak 412 MB" (cpu < 0 or peak 0 left out), styled
 * dim when color. Ends with "\n". */
int term_took_line(char* buf, size_t cap, double elapsed, double cpu, uint64_t peak,
                   TermSpinnerStyle st);

/* "last 3.21 s · mem 1.4 GB · Global`", dim when color. Ends with "\n". */
int term_status_line(char* buf, size_t cap, double last, uint64_t mem,
                     const char* context, TermSpinnerStyle st);

#endif /* MATHILDA_TERM_SPINNER_H */
