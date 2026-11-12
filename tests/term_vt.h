/*
 * term_vt.h -- a minimal VT100 screen for the editor tests.
 *
 * It interprets exactly what the editor emits (printable UTF-8, CR, LF with
 * scrolling, CSI A/B/C/D/H/J/K, SGR and the private modes, which it records
 * or ignores) into a rows x cols character grid, so a test can assert what
 * the user would see and where the cursor is. Header-only, test use only.
 */
#ifndef MATHILDA_TERM_VT_H
#define MATHILDA_TERM_VT_H

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "term/term_width.h"

#define VT_MAX_ROWS 64
#define VT_MAX_COLS 200

typedef struct {
    int  rows, cols;
    char cell[VT_MAX_ROWS][VT_MAX_COLS][5];   /* UTF-8 of each cell, "" = right half */
    int  row, col;
    bool wrap_pending;
    int  scrolled;                            /* lines scrolled off the top          */
    size_t bytes, writes;
    char esc[64];
    size_t esc_len;
    char sgr_log[8192];                       /* every SGR sequence, concatenated   */
    size_t sgr_len;
} Vt;

static inline void vt_init(Vt* v, int cols, int rows) {
    memset(v, 0, sizeof(*v));
    v->cols = cols;
    v->rows = rows;
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) strcpy(v->cell[r][c], " ");
}

static inline void vt_scroll(Vt* v) {
    memmove(v->cell[0], v->cell[1], sizeof(v->cell[0]) * (size_t)(v->rows - 1));
    for (int c = 0; c < v->cols; c++) strcpy(v->cell[v->rows - 1][c], " ");
    v->scrolled++;
}

static inline void vt_clear_eol(Vt* v, int r, int c) {
    for (; c < v->cols; c++) strcpy(v->cell[r][c], " ");
}

static inline void vt_csi(Vt* v, const char* p, size_t n, char f) {
    int a = 0;
    bool priv = n > 0 && p[0] == '?';
    for (size_t i = priv ? 1 : 0; i < n && p[i] >= '0' && p[i] <= '9'; i++) a = a * 10 + (p[i] - '0');
    int k = a ? a : 1;
    if (priv) return;
    switch (f) {
    case 'A': v->row -= k; if (v->row < 0) v->row = 0; v->wrap_pending = false; break;
    case 'B': v->row += k; if (v->row >= v->rows) v->row = v->rows - 1; v->wrap_pending = false; break;
    case 'C': v->col += k; if (v->col >= v->cols) v->col = v->cols - 1; v->wrap_pending = false; break;
    case 'D': v->col -= k; if (v->col < 0) v->col = 0; v->wrap_pending = false; break;
    case 'H': v->row = 0; v->col = 0; v->wrap_pending = false; break;
    case 'K': vt_clear_eol(v, v->row, v->col); break;
    case 'J':
        if (a == 2) { for (int r = 0; r < v->rows; r++) vt_clear_eol(v, r, 0); }
        else { vt_clear_eol(v, v->row, v->col); for (int r = v->row + 1; r < v->rows; r++) vt_clear_eol(v, r, 0); }
        break;
    case 'm':
        if (v->sgr_len + n + 4 < sizeof(v->sgr_log)) {
            v->sgr_log[v->sgr_len++] = '[';
            memcpy(v->sgr_log + v->sgr_len, p, n);
            v->sgr_len += n;
            v->sgr_log[v->sgr_len++] = 'm';
            v->sgr_log[v->sgr_len] = '\0';
        }
        break;
    default: break;                      /* n (position report): no reply */
    }
}

static inline void vt_put(Vt* v, const char* s, size_t n, int w) {
    if (v->wrap_pending) {
        v->col = 0;
        if (v->row == v->rows - 1) vt_scroll(v); else v->row++;
        v->wrap_pending = false;
    }
    if (w == 2 && v->col == v->cols - 1) { strcpy(v->cell[v->row][v->col], " "); v->col = 0;
        if (v->row == v->rows - 1) vt_scroll(v); else v->row++; }
    memcpy(v->cell[v->row][v->col], s, n);
    v->cell[v->row][v->col][n] = '\0';
    if (w == 2) v->cell[v->row][v->col + 1][0] = '\0';
    v->col += w;
    if (v->col >= v->cols) { v->col = v->cols - 1; v->wrap_pending = true; }
}

static inline void vt_feed(Vt* v, const char* s, size_t n) {
    v->bytes += n;
    v->writes++;
    size_t i = 0;
    while (i < n) {
        unsigned char c = (unsigned char)s[i];
        if (v->esc_len) {
            if (v->esc_len < sizeof(v->esc) - 1) v->esc[v->esc_len++] = (char)c;
            i++;
            if (v->esc_len == 2 && c != '[') { v->esc_len = 0; continue; }
            if (v->esc_len > 2 && c >= 0x40 && c <= 0x7E) {
                vt_csi(v, v->esc + 2, v->esc_len - 3, (char)c);
                v->esc_len = 0;
            }
            continue;
        }
        if (c == 0x1b) { v->esc[0] = 0x1b; v->esc_len = 1; i++; continue; }
        if (c == '\r') { v->col = 0; v->wrap_pending = false; i++; continue; }
        if (c == '\n') {
            if (v->row == v->rows - 1) vt_scroll(v); else v->row++;
            v->wrap_pending = false;
            i++;
            continue;
        }
        if (c < 0x20) { i++; continue; }
        uint32_t cp;
        size_t j = term_utf8_next(s, n, i, &cp);
        int w = term_cp_width(cp);
        if (w <= 0) { i = j; continue; }
        vt_put(v, s + i, j - i, w);
        i = j;
    }
}

/* Row r as a string, trailing spaces removed. */
static inline const char* vt_row(const Vt* v, int r) {
    static char buf[VT_MAX_COLS * 4 + 1];
    size_t k = 0;
    for (int c = 0; c < v->cols; c++) {
        const char* s = v->cell[r][c];
        size_t n = strlen(s);
        memcpy(buf + k, s, n);
        k += n;
    }
    while (k > 0 && buf[k - 1] == ' ') k--;
    buf[k] = '\0';
    return buf;
}

static inline void vt_dump(const Vt* v) {
    for (int r = 0; r < v->rows; r++) fprintf(stderr, "%2d|%s|\n", r, vt_row(v, r));
    fprintf(stderr, "cursor %d,%d\n", v->row, v->col);
}

static inline void vt_write_cb(void* ctx, const char* s, size_t n) { vt_feed((Vt*)ctx, s, n); }

#endif /* MATHILDA_TERM_VT_H */
