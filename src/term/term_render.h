/*
 * term_render.h -- the line editor's screen model: a cell grid diffed against
 * the previous frame, emitted as one buffered write.
 *
 * Pure C99, no terminal calls: output goes to a write callback, so the editor
 * runs unchanged against an in-memory sink in the unit tests.
 *
 * The editor owns a rectangular "area" of the screen whose first row is the
 * prompt row. Rows are addressed relative to that top row: moving down past
 * the area's last row emits "\r\n" (which scrolls the screen when the area is
 * at the bottom), moving within it emits CUU/CUD. Relative positions stay valid
 * under scrolling because the area is never taller than the screen. The
 * absolute row of the area's top is tracked when known (a cursor position
 * report at the start of a read), else assumed to be the last screen row; it
 * only decides where a popup goes.
 *
 * term_screen_flush compares each frame row with the previous one and redraws
 * a row only from its first changed cell; unchanged rows cost nothing.
 */
#ifndef MATHILDA_TERM_RENDER_H
#define MATHILDA_TERM_RENDER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    char    ch[4];
    uint8_t len;     /* bytes in ch; 0 for the right half of a wide character */
    uint8_t width;   /* 1 or 2; 0 for the right half                          */
    uint8_t style;   /* a TermHlStyle                                         */
    uint8_t pad;     /* 1: put by the caller (content), 0: untouched blank     */
} TermCell;

typedef void (*TermWriteFn)(void* ctx, const char* s, size_t n);
typedef const char* (*TermSgrFn)(int style);

/* Allocation hooks, so the editor can count allocations (latency gate). */
typedef struct {
    void* (*realloc_fn)(void* p, size_t n);
    void  (*free_fn)(void* p);
} TermAlloc;

typedef struct {
    int       cols, rows;         /* terminal size                               */
    TermCell* grid;               /* the frame being built, rows x cols          */
    TermCell* prev;               /* what is on screen                           */
    int       prev_rows;          /* rows of prev that are valid (0: redraw all) */
    int       frame_rows;         /* rows of grid in use for this frame          */
    int       area_rows;          /* rows the area occupies on screen            */
    int       cur_row, cur_col;   /* terminal cursor, relative to the area top;
                                     cur_col == cols: pending wrap / unknown     */
    int       top_abs;            /* absolute screen row of the area top         */
    bool      top_known;
    bool      color;
    int       cur_style;          /* SGR state on the terminal                   */
    TermSgrFn sgr;
    char*     out;                /* pending output                              */
    size_t    out_len, out_cap;
    TermWriteFn write;
    void*     write_ctx;
    size_t    writes;             /* sink writes so far (the one-write gate)     */
    bool      oom;
    TermAlloc alloc;
} TermScreen;

bool term_screen_init(TermScreen* s, int cols, int rows, bool color, TermSgrFn sgr,
                      TermWriteFn w, void* ctx, const TermAlloc* a);
void term_screen_free(TermScreen* s);
/* New terminal size: the grids are reallocated and the whole area will be
 * redrawn by the next flush. The caller has repositioned the terminal. */
bool term_screen_resize(TermScreen* s, int cols, int rows);

/* Start a frame of `rows` rows (clamped to the screen), all blank. */
void term_screen_begin_frame(TermScreen* s, int rows);
/* Grow the current frame to `rows` rows (clamped), the new rows blank. */
void term_screen_extend_frame(TermScreen* s, int rows);
/* Put one character: `width` 1 or 2 (a wide character that does not fit at
 * the row's end is not drawn). Out-of-range positions are ignored. */
void term_screen_put(TermScreen* s, int row, int col, const char* bytes, size_t len,
                     int width, int style);
/* Put an ASCII string from (row, col), clipped at the row end; returns the
 * column after it. */
int  term_screen_puts(TermScreen* s, int row, int col, const char* str, int style);
/* Restyle cells [c0, c1) of a row. */
void term_screen_style_range(TermScreen* s, int row, int c0, int c1, int style);

/* Diff the frame against the screen, append the drawing to the output buffer,
 * then move the terminal cursor to (crow, ccol). Nothing is written yet. */
void term_screen_render(TermScreen* s, int crow, int ccol);
/* Append raw bytes / move the cursor / append SGR reset. */
void term_screen_raw(TermScreen* s, const char* bytes, size_t n);
void term_screen_move(TermScreen* s, int row, int col);
void term_screen_reset_style(TermScreen* s);
/* Write the output buffer with one call of the sink (no call when empty). */
void term_screen_commit(TermScreen* s);

/* Forget what is on screen: the next frame starts a fresh area at the
 * terminal cursor's current line (after Ctrl-L, a resume, an accepted line). */
void term_screen_reset_area(TermScreen* s, int top_abs, bool top_known);
/* Force every row of the next frame to be drawn in full. */
void term_screen_invalidate(TermScreen* s);

/* Rows available below / above area row `row` without scrolling. */
int  term_screen_rows_below(const TermScreen* s, int row);

#endif /* MATHILDA_TERM_RENDER_H */
