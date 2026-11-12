/*
 * term_render.c -- see term_render.h.
 */
#include "term_render.h"

#include <stdio.h>
#include <string.h>

static const TermCell BLANK = { {' ', 0, 0, 0}, 1, 1, 0, 0 };

static bool cell_eq(const TermCell* a, const TermCell* b) {
    return a->len == b->len && a->width == b->width && a->style == b->style
        && a->pad == b->pad && memcmp(a->ch, b->ch, 4) == 0;
}
/* Never drawn: a cell no put() touched (pad 0). A space the editor put there
 * (typed text, the prompt) is content and is drawn. */
static bool cell_blank(const TermCell* c) {
    return c->pad == 0 && c->len == 1 && c->ch[0] == ' ' && c->style == 0;
}

static void emit(TermScreen* s, const char* b, size_t n) {
    if (s->out_len + n > s->out_cap) {
        size_t nc = s->out_cap ? s->out_cap : 4096;
        while (nc < s->out_len + n) nc *= 2;
        char* no = s->alloc.realloc_fn(s->out, nc);
        if (!no) { s->oom = true; return; }
        s->out = no;
        s->out_cap = nc;
    }
    memcpy(s->out + s->out_len, b, n);
    s->out_len += n;
}

static void emits(TermScreen* s, const char* b) { emit(s, b, strlen(b)); }

static void emit_csi(TermScreen* s, int n, char final) {
    char b[24];
    int k = snprintf(b, sizeof(b), "\033[%d%c", n, final);
    emit(s, b, (size_t)k);
}

bool term_screen_init(TermScreen* s, int cols, int rows, bool color, TermSgrFn sgr,
                      TermWriteFn w, void* ctx, const TermAlloc* a) {
    memset(s, 0, sizeof(*s));
    s->alloc = *a;
    s->color = color;
    s->sgr = sgr;
    s->write = w;
    s->write_ctx = ctx;
    s->top_abs = rows > 0 ? rows - 1 : 0;
    if (!term_screen_resize(s, cols, rows)) { term_screen_free(s); return false; }
    s->out_cap = 8192;
    s->out = a->realloc_fn(NULL, s->out_cap);
    if (!s->out) { term_screen_free(s); return false; }
    s->cur_row = 0;
    s->cur_col = 0;
    s->area_rows = 1;
    return true;
}

void term_screen_free(TermScreen* s) {
    if (!s->alloc.free_fn) return;
    s->alloc.free_fn(s->grid);
    s->alloc.free_fn(s->prev);
    s->alloc.free_fn(s->out);
    s->grid = s->prev = NULL;
    s->out = NULL;
}

bool term_screen_resize(TermScreen* s, int cols, int rows) {
    if (cols < 1) cols = 1;
    if (rows < 1) rows = 1;
    size_t n = (size_t)cols * (size_t)rows;
    TermCell* g = s->alloc.realloc_fn(s->grid, n * sizeof(TermCell));
    if (!g) return false;
    s->grid = g;
    TermCell* p = s->alloc.realloc_fn(s->prev, n * sizeof(TermCell));
    if (!p) return false;
    s->prev = p;
    s->cols = cols;
    s->rows = rows;
    s->prev_rows = 0;
    s->frame_rows = 0;
    if (s->area_rows > rows) s->area_rows = rows;
    if (s->cur_row >= rows) s->cur_row = rows - 1;
    if (s->top_abs >= rows) s->top_abs = rows - 1;
    s->cur_col = cols;            /* unknown: the next move starts with \r */
    return true;
}

void term_screen_begin_frame(TermScreen* s, int rows) {
    if (rows < 1) rows = 1;
    if (rows > s->rows) rows = s->rows;
    s->frame_rows = rows;
    size_t n = (size_t)rows * (size_t)s->cols;
    for (size_t i = 0; i < n; i++) s->grid[i] = BLANK;
}

void term_screen_extend_frame(TermScreen* s, int rows) {
    if (rows > s->rows) rows = s->rows;
    if (rows <= s->frame_rows) return;
    for (size_t i = (size_t)s->frame_rows * (size_t)s->cols; i < (size_t)rows * (size_t)s->cols; i++)
        s->grid[i] = BLANK;
    s->frame_rows = rows;
}

void term_screen_put(TermScreen* s, int row, int col, const char* bytes, size_t len,
                     int width, int style) {
    if (row < 0 || row >= s->frame_rows || col < 0 || col >= s->cols || len == 0 || len > 4)
        return;
    if (width == 2 && col + 1 >= s->cols) return;
    TermCell* c = &s->grid[(size_t)row * (size_t)s->cols + (size_t)col];
    memset(c->ch, 0, 4);
    memcpy(c->ch, bytes, len);
    c->len = (uint8_t)len;
    c->width = (uint8_t)(width == 2 ? 2 : 1);
    c->style = (uint8_t)style;
    c->pad = 1;
    if (width == 2) {
        TermCell* r = c + 1;
        memset(r->ch, 0, 4);
        r->len = 0;
        r->width = 0;
        r->style = (uint8_t)style;
        r->pad = 1;
    }
}

int term_screen_puts(TermScreen* s, int row, int col, const char* str, int style) {
    for (; *str && col < s->cols; str++, col++) term_screen_put(s, row, col, str, 1, 1, style);
    return col;
}

void term_screen_style_range(TermScreen* s, int row, int c0, int c1, int style) {
    if (row < 0 || row >= s->frame_rows) return;
    if (c0 < 0) c0 = 0;
    if (c1 > s->cols) c1 = s->cols;
    TermCell* r = &s->grid[(size_t)row * (size_t)s->cols];
    for (int c = c0; c < c1; c++) r[c].style = (uint8_t)style;
}

void term_screen_raw(TermScreen* s, const char* bytes, size_t n) { emit(s, bytes, n); }

void term_screen_reset_style(TermScreen* s) {
    if (s->cur_style != 0) {
        if (s->color) emits(s, "\033[0m");
        s->cur_style = 0;
    }
}

static void set_style(TermScreen* s, int st) {
    if (st == s->cur_style) return;
    if (s->color) {
        emits(s, "\033[0m");
        const char* q = s->sgr ? s->sgr(st) : "";
        if (q && *q) emits(s, q);
    }
    s->cur_style = st;
}

/* A "\r\n" from area row `row`: when that row is the screen's last one the
 * terminal scrolls, which moves the area's top up by one. */
static void newline_from(TermScreen* s, int row) {
    emits(s, "\r\n");
    if (s->top_abs + row >= s->rows - 1) {
        if (s->top_abs > 0) s->top_abs--;
    }
}

void term_screen_move(TermScreen* s, int row, int col) {
    if (row < 0) row = 0;
    if (row >= s->rows) row = s->rows - 1;
    if (row > s->cur_row) {
        int existing = s->area_rows - 1 - s->cur_row;   /* rows CUD can reach */
        int d = row - s->cur_row;
        int down = d < existing ? d : (existing > 0 ? existing : 0);
        if (down > 0) {
            emit_csi(s, down, 'B');
            s->cur_row += down;
        }
        while (s->cur_row < row) {
            term_screen_reset_style(s);
            newline_from(s, s->cur_row);
            s->cur_row++;
            s->cur_col = 0;
            if (s->cur_row >= s->area_rows) s->area_rows = s->cur_row + 1;
        }
    } else if (row < s->cur_row) {
        emit_csi(s, s->cur_row - row, 'A');
        s->cur_row = row;
    }
    if (col != s->cur_col) {
        if (col == 0 || s->cur_col >= s->cols) {
            emits(s, "\r");
            if (col > 0) emit_csi(s, col, 'C');
        } else if (col > s->cur_col) {
            emit_csi(s, col - s->cur_col, 'C');
        } else {
            emit_csi(s, s->cur_col - col, 'D');
        }
        s->cur_col = col;
    }
}

/* Index past the last non-blank cell of a row. */
static int row_extent(const TermCell* r, int cols) {
    int e = cols;
    while (e > 0 && cell_blank(&r[e - 1])) e--;
    return e;
}

void term_screen_render(TermScreen* s, int crow, int ccol) {
    int cols = s->cols;
    int nrows = s->frame_rows;
    if (s->area_rows > nrows) nrows = s->area_rows < s->rows ? s->area_rows : s->rows;
    for (int r = s->frame_rows; r < nrows; r++)          /* rows below the frame: blank */
        for (int c = 0; c < cols; c++) s->grid[(size_t)r * (size_t)cols + (size_t)c] = BLANK;

    for (int r = 0; r < nrows; r++) {
        TermCell* g = &s->grid[(size_t)r * (size_t)cols];
        TermCell* p = &s->prev[(size_t)r * (size_t)cols];
        bool known = r < s->prev_rows && r < s->area_rows;
        int first = 0;
        int old_ext;
        if (known) {
            while (first < cols && cell_eq(&g[first], &p[first])) first++;
            if (first == cols) continue;                  /* unchanged row */
            while (first > 0 && g[first].len == 0) first--; /* not mid wide char */
            old_ext = row_extent(p, cols);
        } else {
            old_ext = cols;                               /* unknown: clear the rest */
        }
        int new_ext = row_extent(g, cols);
        if (!known && r >= s->area_rows && new_ext == 0) {
            /* A new, empty row: just make it part of the area. */
            term_screen_move(s, r, 0);
            continue;
        }
        if (first > new_ext && first >= old_ext) continue;
        term_screen_move(s, r, first < new_ext ? first : new_ext);
        for (int c = first; c < new_ext; c++) {
            const TermCell* cell = &g[c];
            if (cell->len == 0) continue;                 /* right half of a wide char */
            set_style(s, cell->style);
            emit(s, cell->ch, cell->len);
            s->cur_col += cell->width;
        }
        if (s->cur_col > cols) s->cur_col = cols;
        if (old_ext > new_ext || !known) {
            term_screen_reset_style(s);
            if (s->cur_col < cols) emits(s, "\033[K");
        }
        term_screen_reset_style(s);
    }
    memcpy(s->prev, s->grid, (size_t)nrows * (size_t)cols * sizeof(TermCell));
    s->prev_rows = nrows;
    if (s->area_rows < nrows) s->area_rows = nrows;
    term_screen_reset_style(s);
    term_screen_move(s, crow, ccol);
}

void term_screen_commit(TermScreen* s) {
    if (s->out_len == 0) return;
    if (s->write) s->write(s->write_ctx, s->out, s->out_len);
    s->writes++;
    s->out_len = 0;
}

void term_screen_reset_area(TermScreen* s, int top_abs, bool top_known) {
    s->cur_row = 0;
    s->cur_col = 0;
    s->area_rows = 1;
    s->prev_rows = 0;
    s->cur_style = 0;
    s->top_abs = top_abs < 0 ? s->rows - 1 : (top_abs >= s->rows ? s->rows - 1 : top_abs);
    s->top_known = top_known && top_abs >= 0;
}

void term_screen_invalidate(TermScreen* s) { s->prev_rows = 0; }

int term_screen_rows_below(const TermScreen* s, int row) {
    int r = s->rows - 1 - (s->top_abs + row);
    return r < 0 ? 0 : r;
}
