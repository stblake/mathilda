/*
 * term_editor.c -- see term_editor.h.
 *
 * Layout of the state:
 *   text/sty    the buffer and one style id per byte, NUL-terminated text;
 *   lstart      byte offset of each logical line;
 *   lstate      the lexer state at the start of each line (term_highlight);
 *   lrows       visual rows each line occupies at the current width (0 =
 *               not yet laid out);
 *   inq         input bytes not yet turned into keys (a partial escape
 *               sequence, or typeahead after an accepted line);
 *   ops/arena   the undo log.
 * Every buffer grows by doubling and is kept across reads, so steady-state
 * editing allocates nothing.
 */
#include "term_editor.h"
#include "term_history.h"
#include "term_lex.h"
#include "term_width.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- counted allocation --------------------------------------------------- */

static size_t g_ed_allocs = 0;

static void* ed_realloc(void* p, size_t n) {
    g_ed_allocs++;
    return realloc(p, n);
}
static void ed_free(void* p) { free(p); }
static const TermAlloc ED_ALLOC = { ed_realloc, ed_free };

size_t term_editor_alloc_count(void) { return g_ed_allocs; }

typedef struct { char* p; size_t len, cap; } Bytes;

static bool bytes_reserve(Bytes* b, size_t need) {
    if (need <= b->cap) return true;
    size_t nc = b->cap ? b->cap : 256;
    while (nc < need) nc *= 2;
    char* np = ed_realloc(b->p, nc);
    if (!np) return false;
    b->p = np;
    b->cap = nc;
    return true;
}
static bool bytes_append(Bytes* b, const char* s, size_t n) {
    if (!bytes_reserve(b, b->len + n + 1)) return false;
    memcpy(b->p + b->len, s, n);
    b->len += n;
    b->p[b->len] = '\0';
    return true;
}
static bool bytes_set(Bytes* b, const char* s, size_t n) {
    b->len = 0;
    return bytes_append(b, s, n);
}

/* ---- keys ------------------------------------------------------------------ */

enum {
    K_NONE = 0, K_CHAR, K_ENTER, K_FORCE_ENTER, K_TAB, K_BTAB, K_BS, K_DEL,
    K_UP, K_DOWN, K_LEFT, K_RIGHT, K_HOME, K_END, K_PGUP, K_PGDN,
    K_WLEFT, K_WRIGHT, K_ESC, K_PASTE, K_ALT_BS, K_ALT, K_CTRL
};

typedef struct {
    int    code;
    char   b[4];     /* K_CHAR: the UTF-8 bytes; K_ALT / K_CTRL: the letter */
    size_t n;
} Key;

/* ---- undo ----------------------------------------------------------------- */

enum { EK_NONE = 0, EK_INSERT, EK_DEL_BACK, EK_DEL_FWD, EK_OTHER };

typedef struct {
    size_t   pos;
    size_t   del_off, del_len;   /* the text removed, in the arena  */
    size_t   ins_off, ins_len;   /* the text inserted, in the arena */
    size_t   cur_before;
    unsigned group;
} UndoOp;

#define UNDO_MAX_OPS   2048
#define UNDO_MAX_BYTES (1u << 20)
#define BR_STACK       256
#define BR_MAX_BAD     64

/* ---- the editor ------------------------------------------------------------ */

struct TermEditor {
    TermEditorConfig cfg;
    TermScreen   scr;

    char*        text;
    uint8_t*     sty;
    size_t       len, cap;
    size_t       cur;

    size_t*      lstart;
    TermHlState* lstate;
    uint32_t*    lrows;
    size_t       nlines, lcap;
    int          lay_g0, lay_gn, lay_cols;   /* the geometry lrows were computed for */

    char         prompt[160];
    int          prompt_w;

    Bytes        inq;
    bool         in_paste;
    bool         paste_cr;

    bool         hbrowse;
    size_t       hidx;
    Bytes        draft;

    int          goal;
    bool         goal_ok, goal_keep;

    Bytes        kill;
    int          last_kill, this_kill;       /* 0, +1 forward, -1 backward */

    UndoOp*      ops;
    size_t       nops, opcap;
    Bytes        arena;
    unsigned     group;
    int          last_kind;
    bool         edited;                    /* an edit happened during this key */
    bool         undoing;

    int          vtop;
    bool         final_frame;

    bool         br_dirty;
    size_t       br_stack[BR_STACK];
    size_t       bad[BR_MAX_BAD];
    size_t       nbad;

    TermEdStatus status;
};

/* ---- lines and relexing ---------------------------------------------------- */

static size_t line_of(const TermEditor* ed, size_t pos) {
    size_t lo = 0, hi = ed->nlines;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (ed->lstart[mid] <= pos) lo = mid; else hi = mid;
    }
    return lo;
}

static size_t line_end(const TermEditor* ed, size_t i) {
    return i + 1 < ed->nlines ? ed->lstart[i + 1] - 1 : ed->len;
}

static bool lines_reserve(TermEditor* ed, size_t n) {
    if (n <= ed->lcap) return true;
    size_t nc = ed->lcap ? ed->lcap : 64;
    while (nc < n) nc *= 2;
    size_t* a = ed_realloc(ed->lstart, nc * sizeof(size_t));
    if (!a) return false;
    ed->lstart = a;
    TermHlState* b = ed_realloc(ed->lstate, nc * sizeof(TermHlState));
    if (!b) return false;
    ed->lstate = b;
    uint32_t* c = ed_realloc(ed->lrows, nc * sizeof(uint32_t));
    if (!c) return false;
    ed->lrows = c;
    ed->lcap = nc;
    return true;
}

static void relex(TermEditor* ed, size_t from, size_t must_to) {
    TermHlState st = from == 0 ? term_hl_state_initial() : ed->lstate[from];
    for (size_t i = from; i < ed->nlines; i++) {
        if (i > must_to && term_hl_state_eq(ed->lstate[i], st)) break;
        ed->lstate[i] = st;
        size_t a = ed->lstart[i], b = line_end(ed, i);
        st = term_hl_line(ed->text + a, b - a, st, ed->sty + a,
                          ed->cfg.classify, ed->cfg.classify_ctx);
        if (b < ed->len) ed->sty[b] = TERM_HL_PLAIN;
    }
}

static size_t count_nl(const char* s, size_t n) {
    size_t k = 0;
    for (size_t i = 0; i < n; i++) if (s[i] == '\n') k++;
    return k;
}

/* ---- undo log -------------------------------------------------------------- */

static void undo_compact(TermEditor* ed) {
    /* Drop the oldest half, at a group boundary. */
    size_t k = ed->nops / 2;
    while (k < ed->nops && k > 0 && ed->ops[k].group == ed->ops[k - 1].group) k++;
    if (k == 0 || k >= ed->nops) { ed->nops = 0; ed->arena.len = 0; return; }
    size_t base = ed->ops[k].del_off;
    memmove(ed->arena.p, ed->arena.p + base, ed->arena.len - base);
    ed->arena.len -= base;
    memmove(ed->ops, ed->ops + k, (ed->nops - k) * sizeof(UndoOp));
    ed->nops -= k;
    for (size_t i = 0; i < ed->nops; i++) {
        ed->ops[i].del_off -= base;
        ed->ops[i].ins_off -= base;
    }
}

static void undo_record(TermEditor* ed, size_t pos, size_t dlen, const char* ins,
                        size_t ilen, int kind) {
    UndoOp* last = ed->nops ? &ed->ops[ed->nops - 1] : NULL;
    bool same = last && kind != EK_OTHER && kind == ed->last_kind;
    if (same && kind == EK_INSERT && dlen == 0 && last->pos + last->ins_len == pos
        && last->ins_off + last->ins_len == ed->arena.len
        && !(ilen == 1 && ins[0] == '\n')) {
        if (bytes_append(&ed->arena, ins, ilen)) last->ins_len += ilen;
        return;
    }
    if (same && kind == EK_DEL_BACK && !(last->pos == pos + dlen)) same = false;
    if (same && kind == EK_DEL_FWD && !(last->pos == pos)) same = false;
    if (same && kind == EK_INSERT && !(last->pos + last->ins_len == pos)) same = false;
    if (!same) ed->group++;
    ed->last_kind = kind;
    if (ed->nops >= UNDO_MAX_OPS || ed->arena.len > UNDO_MAX_BYTES) undo_compact(ed);
    if (ed->nops == ed->opcap) {
        size_t nc = ed->opcap ? ed->opcap * 2 : 256;
        UndoOp* no = ed_realloc(ed->ops, nc * sizeof(UndoOp));
        if (!no) return;
        ed->ops = no;
        ed->opcap = nc;
    }
    UndoOp op;
    op.pos = pos;
    op.cur_before = ed->cur;
    op.group = ed->group;
    op.del_off = ed->arena.len;
    if (!bytes_append(&ed->arena, ed->text + pos, dlen)) return;
    op.del_len = dlen;
    op.ins_off = ed->arena.len;
    if (!bytes_append(&ed->arena, ins, ilen)) return;
    op.ins_len = ilen;
    ed->ops[ed->nops++] = op;
}

/* ---- the one buffer mutation ----------------------------------------------- */

static bool ed_replace(TermEditor* ed, size_t pos, size_t dlen, const char* ins,
                       size_t ilen, int kind) {
    if (pos > ed->len) pos = ed->len;
    if (dlen > ed->len - pos) dlen = ed->len - pos;
    if (dlen == 0 && ilen == 0) return true;
    size_t nlen = ed->len - dlen + ilen;
    if (nlen + 1 > ed->cap) {
        size_t nc = ed->cap ? ed->cap : 1024;
        while (nc < nlen + 1) nc *= 2;
        char* nt = ed_realloc(ed->text, nc);
        if (!nt) return false;
        ed->text = nt;
        uint8_t* ns = ed_realloc(ed->sty, nc);
        if (!ns) return false;
        ed->sty = ns;
        ed->cap = nc;
    }
    size_t dn = count_nl(ed->text + pos, dlen);
    size_t ni = count_nl(ins, ilen);
    if (!lines_reserve(ed, ed->nlines - dn + ni)) return false;
    if (!ed->undoing) undo_record(ed, pos, dlen, ins, ilen, kind);
    ed->edited = true;

    size_t L = line_of(ed, pos);
    memmove(ed->text + pos + ilen, ed->text + pos + dlen, ed->len - pos - dlen + 1);
    memmove(ed->sty + pos + ilen, ed->sty + pos + dlen, ed->len - pos - dlen);
    memcpy(ed->text + pos, ins, ilen);
    ed->len = nlen;

    if (dn != ni) {
        size_t tail = ed->nlines - (L + 1 + dn);      /* lines after the edited range */
        memmove(ed->lstart + L + 1 + ni, ed->lstart + L + 1 + dn, tail * sizeof(size_t));
        memmove(ed->lstate + L + 1 + ni, ed->lstate + L + 1 + dn, tail * sizeof(TermHlState));
        memmove(ed->lrows + L + 1 + ni, ed->lrows + L + 1 + dn, tail * sizeof(uint32_t));
        ed->nlines = ed->nlines - dn + ni;
    }
    for (size_t i = L + 1 + ni; i < ed->nlines; i++) ed->lstart[i] = ed->lstart[i] - dlen + ilen;
    size_t at = pos;
    for (size_t k = 1; k <= ni; k++) {
        const char* nl = memchr(ed->text + at, '\n', pos + ilen - at);
        at = (size_t)(nl - ed->text) + 1;
        ed->lstart[L + k] = at;
        ed->lstate[L + k].invalid = 1;
    }
    for (size_t k = 0; k <= ni; k++) ed->lrows[L + k] = 0;
    relex(ed, L, L + ni);
    ed->br_dirty = true;
    return true;
}

static void insert(TermEditor* ed, const char* s, size_t n, int kind) {
    if (ed_replace(ed, ed->cur, 0, s, n, kind)) ed->cur += n;
}

static void set_buffer(TermEditor* ed, const char* s, size_t n) {
    if (ed_replace(ed, 0, ed->len, s, n, EK_OTHER)) ed->cur = ed->len;
}

/* ---- cursor helpers ---------------------------------------------------------- */

static size_t prev_char(const TermEditor* ed, size_t i) {
    return i == 0 ? 0 : term_utf8_prev(ed->text, i);
}
static size_t next_char(const TermEditor* ed, size_t i) {
    uint32_t cp;
    return i >= ed->len ? ed->len : term_utf8_next(ed->text, ed->len, i, &cp);
}
static bool word_byte(char c) {
    return term_lex_is_word_char(c) || (unsigned char)c >= 0x80;
}
static size_t word_right(const TermEditor* ed, size_t i) {
    while (i < ed->len && !word_byte(ed->text[i])) i++;
    while (i < ed->len && word_byte(ed->text[i])) i++;
    return i;
}
static size_t word_left(const TermEditor* ed, size_t i) {
    while (i > 0 && !word_byte(ed->text[i - 1])) i--;
    while (i > 0 && word_byte(ed->text[i - 1])) i--;
    return i;
}

static int char_cols(uint32_t cp) {
    int w = term_cp_width(cp);
    return w < 0 ? 2 : w;
}

/* Display column of byte offset `pos` within its logical line. */
static int col_in_line(const TermEditor* ed, size_t pos) {
    size_t a = ed->lstart[line_of(ed, pos)];
    return (int)term_text_width(ed->text + a, pos - a);
}

/* Byte offset in line `l` nearest to display column `goal` (not past it). */
static size_t pos_at_col(const TermEditor* ed, size_t l, int goal) {
    size_t i = ed->lstart[l], e = line_end(ed, l);
    int col = 0;
    while (i < e) {
        uint32_t cp;
        size_t j = term_utf8_next(ed->text, e, i, &cp);
        int w = char_cols(cp);
        if (col + w > goal) break;
        col += w;
        i = j;
    }
    return i;
}

/* ---- kill ring and history --------------------------------------------------- */

static void kill_range(TermEditor* ed, size_t a, size_t b, int dir) {
    if (b <= a) return;
    if (ed->last_kill != 0) {
        if (dir < 0) {                                   /* prepend */
            if (bytes_reserve(&ed->kill, ed->kill.len + (b - a) + 1)) {
                memmove(ed->kill.p + (b - a), ed->kill.p, ed->kill.len + 1);
                memcpy(ed->kill.p, ed->text + a, b - a);
                ed->kill.len += b - a;
            }
        } else {
            bytes_append(&ed->kill, ed->text + a, b - a);
        }
    } else {
        bytes_set(&ed->kill, ed->text + a, b - a);
    }
    ed->this_kill = dir;
    if (ed_replace(ed, a, b - a, "", 0, EK_OTHER)) ed->cur = a;
}

static size_t hist_count(const TermEditor* ed) {
    return ed->cfg.history ? term_history_count(ed->cfg.history) : 0;
}

static void history_show(TermEditor* ed, size_t idx, bool cursor_first_line) {
    const char* s;
    size_t n;
    if (idx >= hist_count(ed)) { s = ed->draft.p ? ed->draft.p : ""; n = ed->draft.len; }
    else { s = term_history_get(ed->cfg.history, idx); n = strlen(s); }
    set_buffer(ed, s, n);
    ed->last_kind = EK_NONE;
    ed->cur = cursor_first_line ? line_end(ed, 0) : ed->len;
}

static void history_up(TermEditor* ed) {
    size_t n = hist_count(ed);
    if (n == 0) return;
    if (!ed->hbrowse) {
        bytes_set(&ed->draft, ed->text, ed->len);
        ed->hbrowse = true;
        ed->hidx = n;
    }
    if (ed->hidx == 0) return;
    ed->hidx--;
    history_show(ed, ed->hidx, true);
}

static void history_down(TermEditor* ed) {
    if (!ed->hbrowse) return;
    size_t n = hist_count(ed);
    if (ed->hidx >= n) return;
    ed->hidx++;
    history_show(ed, ed->hidx, false);
    if (ed->hidx >= n) ed->hbrowse = false;
}

static void move_vertical(TermEditor* ed, int dir) {
    size_t l = line_of(ed, ed->cur);
    if (dir < 0 && l == 0) { history_up(ed); return; }
    if (dir > 0 && l + 1 >= ed->nlines) { history_down(ed); return; }
    if (!ed->goal_ok) { ed->goal = col_in_line(ed, ed->cur); ed->goal_ok = true; }
    ed->cur = pos_at_col(ed, dir < 0 ? l - 1 : l + 1, ed->goal);
    ed->goal_keep = true;
}

/* ---- undo -------------------------------------------------------------------- */

static void undo(TermEditor* ed) {
    if (ed->nops == 0) return;
    unsigned g = ed->ops[ed->nops - 1].group;
    size_t k = ed->nops;
    while (k > 0 && ed->ops[k - 1].group == g) k--;
    ed->undoing = true;
    for (size_t j = ed->nops; j > k; j--) {
        const UndoOp* op = &ed->ops[j - 1];
        ed_replace(ed, op->pos, op->ins_len, ed->arena.p + op->del_off, op->del_len, EK_OTHER);
    }
    ed->undoing = false;
    ed->cur = ed->ops[k].cur_before;
    if (ed->cur > ed->len) ed->cur = ed->len;
    ed->arena.len = ed->ops[k].del_off;
    ed->nops = k;
    ed->last_kind = EK_NONE;
}

/* ---- brackets ---------------------------------------------------------------- */

static bool is_open(char c)  { return c == '(' || c == '[' || c == '{'; }
static bool pair_ok(char o, char c) {
    return (o == '(' && c == ')') || (o == '[' && c == ']') || (o == '{' && c == '}');
}

static void brackets_bad(TermEditor* ed) {
    size_t depth = 0;
    ed->nbad = 0;
    for (size_t i = 0; i < ed->len; i++) {
        if (ed->sty[i] != TERM_HL_BRACKET) continue;
        char c = ed->text[i];
        if (is_open(c)) {
            if (depth < BR_STACK) ed->br_stack[depth] = i;
            depth++;
            continue;
        }
        bool bad;
        if (depth == 0) bad = true;
        else {
            depth--;
            bad = depth < BR_STACK && !pair_ok(ed->text[ed->br_stack[depth]], c);
        }
        if (bad && ed->nbad < BR_MAX_BAD) ed->bad[ed->nbad++] = i;
    }
    ed->br_dirty = false;
}

/* The bracket at or just before the cursor and its partner. */
static void brackets_match(const TermEditor* ed, size_t* ma, size_t* mb) {
    *ma = *mb = SIZE_MAX;
    size_t p;
    if (ed->cur < ed->len && ed->sty[ed->cur] == TERM_HL_BRACKET) p = ed->cur;
    else if (ed->cur > 0 && ed->sty[ed->cur - 1] == TERM_HL_BRACKET) p = ed->cur - 1;
    else return;
    char c = ed->text[p];
    size_t depth = 0;
    if (is_open(c)) {
        for (size_t j = p + 1; j < ed->len; j++) {
            if (ed->sty[j] != TERM_HL_BRACKET) continue;
            if (is_open(ed->text[j])) depth++;
            else if (depth == 0) { if (pair_ok(c, ed->text[j])) { *ma = p; *mb = j; } return; }
            else depth--;
        }
    } else {
        for (size_t j = p; j > 0; j--) {
            size_t q = j - 1;
            if (ed->sty[q] != TERM_HL_BRACKET) continue;
            if (!is_open(ed->text[q])) depth++;
            else if (depth == 0) { if (pair_ok(ed->text[q], c)) { *ma = q; *mb = p; } return; }
            else depth--;
        }
    }
}

/* ---- layout and rendering ---------------------------------------------------- */

/* Gutter widths: g0 for the very first row (the prompt), gn for every other
 * row (continuation indent, dropped on a narrow terminal). */
static void gutters(const TermEditor* ed, int* g0, int* gn) {
    int cols = ed->scr.cols;
    *g0 = ed->prompt_w < cols - 1 ? ed->prompt_w : (cols > 1 ? cols - 1 : 0);
    *gn = (cols - ed->prompt_w >= 20) ? ed->prompt_w : 0;
}

/* Walk line l's characters with wrapping. Returns the number of visual rows;
 * if `upto` lies in the line, crow and ccol receive the cursor position there. */
static uint32_t walk_line(const TermEditor* ed, size_t l, int g0, int gn, size_t upto,
                          int* crow, int* ccol) {
    int cols = ed->scr.cols;
    int avail = (l == 0 ? cols - g0 : cols - gn);
    int availn = cols - gn;
    if (avail < 1) avail = 1;
    if (availn < 1) availn = 1;
    size_t i = ed->lstart[l], e = line_end(ed, l);
    int row = 0, col = 0;
    while (i < e) {
        uint32_t cp;
        size_t j = term_utf8_next(ed->text, e, i, &cp);
        int w = char_cols(cp);
        if (col + w > avail) { row++; col = 0; avail = availn; }
        if (i == upto && crow) { *crow = row; *ccol = col; }
        col += w;
        i = j;
    }
    if (col >= avail) { row++; col = 0; avail = availn; }
    if (upto == e && crow) { *crow = row; *ccol = col; }
    return (uint32_t)row + 1;
}

static int style_at(const TermEditor* ed, size_t i, size_t ma, size_t mb) {
    int st = ed->sty[i];
    if (st != TERM_HL_BRACKET || ed->final_frame) return st;
    if (i == ma || i == mb) return TERM_HL_BRACKET_MATCH;
    for (size_t k = 0; k < ed->nbad; k++) if (ed->bad[k] == i) return TERM_HL_BRACKET_BAD;
    return st;
}

/* Draw the visual rows of line l whose global visual index lies in
 * [vtop, vtop + nrows), the line's first visual row being `vbase`. */
static void draw_line(TermEditor* ed, size_t l, int vbase, int g0, int gn, size_t ma, size_t mb) {
    TermScreen* s = &ed->scr;
    int cols = s->cols;
    int nrows = s->frame_rows;
    int avail = (l == 0 ? cols - g0 : cols - gn);
    int availn = cols - gn;
    if (avail < 1) avail = 1;
    if (availn < 1) availn = 1;
    int gut = (l == 0 ? g0 : gn);
    int row = 0, col = 0;
    int srow = vbase - ed->vtop;
    if (l == 0 && srow >= 0 && srow < nrows) {
        for (int k = 0; k < g0; k++)
            term_screen_put(s, srow, k, &ed->prompt[k], 1, 1, TERM_HL_PROMPT);
    }
    size_t i = ed->lstart[l], e = line_end(ed, l);
    while (i < e) {
        uint32_t cp;
        size_t j = term_utf8_next(ed->text, e, i, &cp);
        int w = term_cp_width(cp);
        int cw = w < 0 ? 2 : w;
        if (col + cw > avail) { row++; col = 0; avail = availn; gut = gn; }
        int r = vbase + row - ed->vtop;
        if (r >= nrows) return;
        if (r >= 0) {
            int st = style_at(ed, i, ma, mb);
            if (w < 0) {
                char hat[2] = { '^', (char)((cp & 0x7F) ^ 0x40) };
                if (cp >= 0x80) hat[1] = '?';
                term_screen_put(s, r, gut + col, &hat[0], 1, 1, st);
                term_screen_put(s, r, gut + col + 1, &hat[1], 1, 1, st);
            } else if (w == 0) {
                /* A combining mark joins the character before it. */
                if (col > 0) {
                    TermCell* c = &s->grid[(size_t)r * (size_t)cols + (size_t)(gut + col - 1)];
                    while (c->len == 0 && c > s->grid) c--;
                    if (c->len + (j - i) <= 4) {
                        memcpy(c->ch + c->len, ed->text + i, j - i);
                        c->len = (uint8_t)(c->len + (j - i));
                    }
                }
            } else {
                term_screen_put(s, r, gut + col, ed->text + i, j - i, w, st);
            }
        }
        col += cw;
        i = j;
    }
}

static void render(TermEditor* ed) {
    TermScreen* s = &ed->scr;
    int g0, gn;
    gutters(ed, &g0, &gn);
    if (g0 != ed->lay_g0 || gn != ed->lay_gn || s->cols != ed->lay_cols) {
        for (size_t i = 0; i < ed->nlines; i++) ed->lrows[i] = 0;
        ed->lay_g0 = g0; ed->lay_gn = gn; ed->lay_cols = s->cols;
    }
    size_t cl = line_of(ed, ed->cur);
    int total = 0, cbase = 0;
    for (size_t i = 0; i < ed->nlines; i++) {
        if (ed->lrows[i] == 0) ed->lrows[i] = walk_line(ed, i, g0, gn, SIZE_MAX, NULL, NULL);
        if (i == cl) cbase = total;
        total += (int)ed->lrows[i];
    }
    int crow = 0, ccol = 0;
    (void)walk_line(ed, cl, g0, gn, ed->cur, &crow, &ccol);
    int cvrow = cbase + crow;
    ccol += (cl == 0 && crow == 0) ? g0 : gn;
    if (ccol >= s->cols) ccol = s->cols - 1;

    int V = s->rows;
    if (total <= V) ed->vtop = 0;
    else {
        if (cvrow < ed->vtop) ed->vtop = cvrow;
        if (cvrow >= ed->vtop + V) ed->vtop = cvrow - V + 1;
        if (ed->vtop > total - V) ed->vtop = total - V;
        if (ed->vtop < 0) ed->vtop = 0;
    }
    int frame = total - ed->vtop < V ? total - ed->vtop : V;
    term_screen_begin_frame(s, frame);

    size_t ma = SIZE_MAX, mb = SIZE_MAX;
    if (!ed->final_frame) {
        if (ed->br_dirty) brackets_bad(ed);
        brackets_match(ed, &ma, &mb);
    }
    int vbase = 0;
    for (size_t i = 0; i < ed->nlines; i++) {
        int rows = (int)ed->lrows[i];
        if (vbase + rows > ed->vtop) {
            if (vbase - ed->vtop >= frame) break;
            draw_line(ed, i, vbase, g0, gn, ma, mb);
        }
        vbase += rows;
    }
    term_screen_render(s, cvrow - ed->vtop, ccol);
}

/* Write the whole buffer as a transcript (prompt, gutters, colours), for an
 * accepted buffer taller than the screen: the viewport only ever showed part
 * of it, and the scrollback should hold all of it. */
static void dump_buffer(TermEditor* ed) {
    TermScreen* s = &ed->scr;
    int g0, gn;
    gutters(ed, &g0, &gn);
    term_screen_move(s, 0, 0);
    term_screen_raw(s, "\033[J", 3);
    char sp[256];
    memset(sp, ' ', sizeof(sp));
    int cur = 0;
    for (size_t l = 0; l < ed->nlines; l++) {
        if (l == 0) {
            if (s->color) { const char* q = term_hl_sgr(TERM_HL_PROMPT); term_screen_raw(s, q, strlen(q)); }
            term_screen_raw(s, ed->prompt, (size_t)g0);
            if (s->color) term_screen_raw(s, "\033[0m", 4);
        } else {
            term_screen_raw(s, "\r\n", 2);
            term_screen_raw(s, sp, (size_t)(gn < 256 ? gn : 256));
        }
        size_t a = ed->lstart[l], e = line_end(ed, l);
        for (size_t i = a; i < e; i++) {
            int st = ed->sty[i];
            if (s->color && st != cur) {
                term_screen_raw(s, "\033[0m", 4);
                const char* q = term_hl_sgr(st);
                term_screen_raw(s, q, strlen(q));
                cur = st;
            }
            unsigned char c = (unsigned char)ed->text[i];
            if (c < 0x20 || c == 0x7F) {
                char hat[2] = { '^', (char)(c ^ 0x40) };
                term_screen_raw(s, hat, 2);
            } else {
                term_screen_raw(s, ed->text + i, 1);
            }
        }
        if (s->color && cur != 0) { term_screen_raw(s, "\033[0m", 4); cur = 0; }
    }
}

/* End the read: the final frame (no bracket highlight), then the cursor goes
 * below the buffer. */
static void finish(TermEditor* ed, TermEdStatus st) {
    TermScreen* s = &ed->scr;
    ed->status = st;
    ed->final_frame = true;
    int g0, gn;
    gutters(ed, &g0, &gn);
    int total = 0;
    for (size_t i = 0; i < ed->nlines; i++) {
        if (ed->lrows[i] == 0) ed->lrows[i] = walk_line(ed, i, g0, gn, SIZE_MAX, NULL, NULL);
        total += (int)ed->lrows[i];
    }
    if (total > s->rows) {
        dump_buffer(ed);
    } else {
        ed->cur = ed->len;
        render(ed);
    }
    if (st == TERM_ED_INTERRUPT) term_screen_raw(s, "^C", 2);
    if (st != TERM_ED_EOF) term_screen_raw(s, "\r\n", 2);
    term_screen_raw(s, "\033[J", 3);
}

/* ---- completion (inline; the popup replaces this in a later phase) ----------- */

static void do_complete(TermEditor* ed) {
    TermCompletions c;
    if (ed->cfg.complete) ed->cfg.complete(ed->text, ed->cur, &c);
    else term_complete(ed->text, ed->cur, &c);
    if (c.count == 0) {
        /* Nothing to complete at the start of a line: indent. */
        size_t a = ed->lstart[line_of(ed, ed->cur)];
        bool blank = true;
        for (size_t i = a; i < ed->cur; i++) if (ed->text[i] != ' ') blank = false;
        if (blank) insert(ed, "  ", 2, EK_OTHER);
        term_completions_free(&c);
        return;
    }
    size_t lcp = strlen(c.items[0]);
    for (size_t i = 1; i < c.count; i++) {
        size_t k = 0;
        while (k < lcp && c.items[i][k] == c.items[0][k]) k++;
        lcp = k;
    }
    size_t wlen = ed->cur - c.word_start;
    if (lcp > wlen && ed_replace(ed, c.word_start, wlen, c.items[0], lcp, EK_OTHER))
        ed->cur = c.word_start + lcp;
    term_completions_free(&c);
}

/* ---- key handling --------------------------------------------------------------- */

static void handle_ctrl(TermEditor* ed, char k) {
    size_t l = line_of(ed, ed->cur);
    switch (k) {
    case 'a': ed->cur = ed->lstart[l]; break;
    case 'e': ed->cur = line_end(ed, l); break;
    case 'b': ed->cur = prev_char(ed, ed->cur); break;
    case 'f': ed->cur = next_char(ed, ed->cur); break;
    case 'd':
        if (ed->len == 0) { finish(ed, TERM_ED_EOF); return; }
        if (ed->cur < ed->len) ed_replace(ed, ed->cur, next_char(ed, ed->cur) - ed->cur, "", 0, EK_DEL_FWD);
        break;
    case 'k': {
        size_t e = line_end(ed, l);
        if (e == ed->cur && e < ed->len) e++;           /* at the end: join lines */
        kill_range(ed, ed->cur, e, +1);
        break;
    }
    case 'u': kill_range(ed, ed->lstart[l], ed->cur, -1); break;
    case 'w': {
        size_t a = ed->cur;
        while (a > 0 && (ed->text[a - 1] == ' ' || ed->text[a - 1] == '\t')) a--;
        while (a > 0 && ed->text[a - 1] != ' ' && ed->text[a - 1] != '\t' && ed->text[a - 1] != '\n') a--;
        kill_range(ed, a, ed->cur, -1);
        break;
    }
    case 'y':
        if (ed->kill.len) insert(ed, ed->kill.p, ed->kill.len, EK_OTHER);
        break;
    case 't': {
        if (ed->cur == 0 || ed->len < 2) break;
        size_t b = ed->cur == line_end(ed, l) ? prev_char(ed, ed->cur) : ed->cur;
        size_t a = prev_char(ed, b), e = next_char(ed, b);
        if (a == b || ed->text[a] == '\n' || ed->text[b] == '\n') break;
        char tmp[8];
        size_t n1 = b - a, n2 = e - b;
        memcpy(tmp, ed->text + b, n2);
        memcpy(tmp + n2, ed->text + a, n1);
        if (ed_replace(ed, a, n1 + n2, tmp, n1 + n2, EK_OTHER)) ed->cur = e;
        break;
    }
    case 'l':                                       /* drawn by the feed's render */
        term_screen_raw(&ed->scr, "\033[H\033[2J", 7);
        term_screen_reset_area(&ed->scr, 0, true);
        break;
    case 'c': finish(ed, TERM_ED_INTERRUPT); break;
    case 'z': finish(ed, TERM_ED_SUSPEND); break;
    case '_': undo(ed); break;
    case 'p': move_vertical(ed, -1); break;
    case 'n': move_vertical(ed, +1); break;
    default: break;
    }
}

static void handle_key(TermEditor* ed, const Key* k) {
    ed->edited = false;
    ed->goal_keep = false;
    ed->this_kill = 0;
    size_t l;
    switch (k->code) {
    case K_CHAR:
        insert(ed, k->b, k->n, EK_INSERT);
        break;
    case K_ENTER:
        if (term_lex_input_complete(ed->text)) finish(ed, TERM_ED_ACCEPT);
        else insert(ed, "\n", 1, EK_OTHER);
        break;
    case K_FORCE_ENTER: finish(ed, TERM_ED_ACCEPT); break;
    case K_TAB: do_complete(ed); break;
    case K_BS:
        if (ed->cur > 0) {
            size_t p = prev_char(ed, ed->cur);
            if (ed_replace(ed, p, ed->cur - p, "", 0, EK_DEL_BACK)) ed->cur = p;
        }
        break;
    case K_DEL:
        if (ed->cur < ed->len) ed_replace(ed, ed->cur, next_char(ed, ed->cur) - ed->cur, "", 0, EK_DEL_FWD);
        break;
    case K_LEFT:  ed->cur = prev_char(ed, ed->cur); break;
    case K_RIGHT: ed->cur = next_char(ed, ed->cur); break;
    case K_UP:    move_vertical(ed, -1); break;
    case K_DOWN:  move_vertical(ed, +1); break;
    case K_HOME:  l = line_of(ed, ed->cur); ed->cur = ed->lstart[l]; break;
    case K_END:   l = line_of(ed, ed->cur); ed->cur = line_end(ed, l); break;
    case K_PGUP:  ed->cur = 0; break;
    case K_PGDN:  ed->cur = ed->len; break;
    case K_WLEFT: ed->cur = word_left(ed, ed->cur); break;
    case K_WRIGHT: ed->cur = word_right(ed, ed->cur); break;
    case K_ALT_BS: kill_range(ed, word_left(ed, ed->cur), ed->cur, -1); break;
    case K_PASTE:
        ed->in_paste = true;
        ed->paste_cr = false;
        ed->last_kind = EK_NONE;
        break;
    case K_ALT:
        switch (k->b[0]) {
        case 'b': case 'B': ed->cur = word_left(ed, ed->cur); break;
        case 'f': case 'F': ed->cur = word_right(ed, ed->cur); break;
        case 'd': case 'D': kill_range(ed, ed->cur, word_right(ed, ed->cur), +1); break;
        case '<':
            if (hist_count(ed)) {
                if (!ed->hbrowse) { bytes_set(&ed->draft, ed->text, ed->len); ed->hbrowse = true; }
                ed->hidx = 0;
                history_show(ed, 0, true);
            }
            break;
        case '>':
            if (ed->hbrowse) { ed->hidx = hist_count(ed); history_show(ed, ed->hidx, false); ed->hbrowse = false; }
            break;
        default: break;
        }
        break;
    case K_CTRL: handle_ctrl(ed, k->b[0]); break;
    default: break;
    }
    if (!ed->edited) ed->last_kind = EK_NONE;
    if (!ed->goal_keep) ed->goal_ok = false;
    ed->last_kill = ed->this_kill;
}

/* ---- input parsing ---------------------------------------------------------------- */

static size_t utf8_seq_len(unsigned char c) {
    if (c >= 0xC2 && c <= 0xDF) return 2;
    if (c >= 0xE0 && c <= 0xEF) return 3;
    if (c >= 0xF0 && c <= 0xF4) return 4;
    return 1;
}

static int csi_letter(char f) {
    switch (f) {
    case 'A': return K_UP;   case 'B': return K_DOWN;
    case 'C': return K_RIGHT; case 'D': return K_LEFT;
    case 'H': return K_HOME; case 'F': return K_END;
    default: return K_NONE;
    }
}

/* Parse one key from p[0, n). Returns the bytes it used, or 0 when the bytes
 * are a prefix of a longer sequence and !final. */
static size_t parse_key(const char* ps, size_t n, bool final, Key* k) {
    const unsigned char* p = (const unsigned char*)ps;
    k->code = K_NONE;
    k->n = 0;
    unsigned char c = p[0];
    if (c == 0x1b) {
        if (n == 1) { if (!final) return 0; k->code = K_ESC; return 1; }
        unsigned char d = p[1];
        if (d == '[') {
            size_t i = 2;
            while (i < n && p[i] >= 0x30 && p[i] <= 0x3F) i++;
            while (i < n && p[i] >= 0x20 && p[i] <= 0x2F) i++;
            if (i >= n) {
                if (!final && n < 32) return 0;
                k->code = K_ESC;
                return 1;
            }
            char f = (char)p[i];
            /* Parameters: the first number and the modifier after ';'. */
            int num = 0, mod = 0;
            size_t j = 2;
            while (j < i && p[j] >= '0' && p[j] <= '9') num = num * 10 + (p[j++] - '0');
            if (j < i && p[j] == ';') {
                j++;
                while (j < i && p[j] >= '0' && p[j] <= '9') mod = mod * 10 + (p[j++] - '0');
            }
            int code = csi_letter(f);
            if (code != K_NONE) {
                bool ctrl = mod == 5 || mod == 6 || mod == 7 || mod == 8;
                bool alt = mod == 3 || mod == 4 || mod == 7 || mod == 8;
                if ((ctrl || alt) && code == K_LEFT) code = K_WLEFT;
                else if ((ctrl || alt) && code == K_RIGHT) code = K_WRIGHT;
                k->code = code;
            } else if (f == 'Z') {
                k->code = K_BTAB;
            } else if (f == '~') {
                switch (num) {
                case 1: case 7: k->code = K_HOME; break;
                case 4: case 8: k->code = K_END; break;
                case 3: k->code = K_DEL; break;
                case 5: k->code = K_PGUP; break;
                case 6: k->code = K_PGDN; break;
                case 200: k->code = K_PASTE; break;
                default: break;                 /* 201 outside a paste, F-keys */
                }
            }
            return i + 1;
        }
        if (d == 'O') {
            if (n < 3) {
                if (!final) return 0;
                k->code = K_ALT; k->b[0] = 'O'; k->n = 1;
                return 2;
            }
            k->code = csi_letter((char)p[2]);
            return 3;
        }
        if (d == '\r' || d == '\n') { k->code = K_FORCE_ENTER; return 2; }
        if (d == 0x7F || d == 0x08) { k->code = K_ALT_BS; return 2; }
        if (d < 0x20 || d >= 0x80) { k->code = K_ESC; return 1; }
        k->code = K_ALT;
        k->b[0] = (char)d;
        k->n = 1;
        return 2;
    }
    if (c >= 0x80) {
        size_t want = utf8_seq_len(c);
        if (want > n && !final) return 0;
        uint32_t cp;
        size_t j = term_utf8_next(ps, n, 0, &cp);
        k->code = K_CHAR;
        memcpy(k->b, ps, j);
        k->n = j;
        return j;
    }
    switch (c) {
    case '\r': case '\n': k->code = K_ENTER; return 1;
    case '\t': k->code = K_TAB; return 1;
    case 0x08: case 0x7F: k->code = K_BS; return 1;
    case 0x1F: k->code = K_CTRL; k->b[0] = '_'; return 1;
    default: break;
    }
    if (c < 0x20) {
        if (c >= 1 && c <= 26) { k->code = K_CTRL; k->b[0] = (char)('a' + c - 1); }
        return 1;
    }
    k->code = K_CHAR;
    k->b[0] = (char)c;
    k->n = 1;
    return 1;
}

static const char PASTE_END[] = "\033[201~";

/* Insert pasted text verbatim (CR LF and CR become LF, a tab becomes spaces,
 * other control characters are dropped) up to the end marker. Returns the
 * bytes used; 0 when only a possible prefix of the marker is left. */
static size_t paste_step(TermEditor* ed, const char* p, size_t n, bool final) {
    char chunk[256];
    size_t m = 0, i = 0;
    while (i < n) {
        if (p[i] == '\033') {
            size_t left = n - i, ml = sizeof(PASTE_END) - 1;
            size_t cmp = left < ml ? left : ml;
            if (memcmp(p + i, PASTE_END, cmp) == 0) {
                if (left >= ml) {
                    if (m) insert(ed, chunk, m, EK_INSERT);
                    ed->in_paste = false;
                    ed->last_kind = EK_NONE;
                    return i + ml;
                }
                if (!final) break;                /* a partial marker: wait */
            }
            i++;
            continue;                             /* other escapes are dropped */
        }
        unsigned char c = (unsigned char)p[i++];
        if (c == '\n' && ed->paste_cr) { ed->paste_cr = false; continue; }
        ed->paste_cr = (c == '\r');
        if (c == '\r') c = '\n';
        if (c == '\t') {
            for (int t = 0; t < 4; t++) chunk[m++] = ' ';
        } else if (c >= 0x20 || c == '\n') {
            chunk[m++] = (char)c;
        }
        if (m + 8 > sizeof(chunk)) { insert(ed, chunk, m, EK_INSERT); m = 0; }
    }
    if (m) insert(ed, chunk, m, EK_INSERT);
    return i;
}

static void process(TermEditor* ed, bool final) {
    size_t i = 0;
    while (i < ed->inq.len && ed->status == TERM_ED_CONTINUE) {
        size_t k;
        if (ed->in_paste) {
            k = paste_step(ed, ed->inq.p + i, ed->inq.len - i, final);
        } else {
            Key key;
            k = parse_key(ed->inq.p + i, ed->inq.len - i, final, &key);
            if (k) handle_key(ed, &key);
        }
        if (k == 0) break;
        i += k;
    }
    memmove(ed->inq.p, ed->inq.p + i, ed->inq.len - i);
    ed->inq.len -= i;
}

/* ---- public API ------------------------------------------------------------------- */

static void sink_write(void* ctx, const char* s, size_t n) {
    TermEditor* ed = (TermEditor*)ctx;
    if (ed->cfg.write) ed->cfg.write(ed->cfg.write_ctx, s, n);
}

TermEditor* term_editor_new(const TermEditorConfig* cfg) {
    TermEditor* ed = calloc(1, sizeof(*ed));
    if (!ed) return NULL;
    ed->cfg = *cfg;
    if (!term_screen_init(&ed->scr, cfg->cols > 0 ? cfg->cols : 80, cfg->rows > 0 ? cfg->rows : 24,
                          cfg->color, term_hl_sgr, sink_write, ed, &ED_ALLOC)) {
        free(ed);
        return NULL;
    }
    ed->cap = 1024;
    ed->text = ed_realloc(NULL, ed->cap);
    ed->sty = ed_realloc(NULL, ed->cap);
    if (!ed->text || !ed->sty || !lines_reserve(ed, 64) || !bytes_reserve(&ed->inq, 4096)
        || !bytes_reserve(&ed->arena, 16384)) {
        term_editor_free(ed);
        return NULL;
    }
    ed->opcap = 256;
    ed->ops = ed_realloc(NULL, ed->opcap * sizeof(UndoOp));
    if (!ed->ops) { term_editor_free(ed); return NULL; }
    ed->text[0] = '\0';
    ed->nlines = 1;
    ed->lstart[0] = 0;
    ed->lstate[0] = term_hl_state_initial();
    ed->lrows[0] = 0;
    return ed;
}

void term_editor_free(TermEditor* ed) {
    if (!ed) return;
    term_screen_free(&ed->scr);
    free(ed->text); free(ed->sty);
    free(ed->lstart); free(ed->lstate); free(ed->lrows);
    free(ed->inq.p); free(ed->draft.p); free(ed->kill.p);
    free(ed->ops); free(ed->arena.p);
    free(ed);
}

void term_editor_begin(TermEditor* ed, const char* prompt, int cursor_abs_row, int start_col) {
    ed->len = 0;
    ed->text[0] = '\0';
    ed->cur = 0;
    ed->nlines = 1;
    ed->lstart[0] = 0;
    ed->lstate[0] = term_hl_state_initial();
    ed->lrows[0] = 0;
    ed->nops = 0;
    ed->arena.len = 0;
    ed->last_kind = EK_NONE;
    ed->hbrowse = false;
    ed->goal_ok = false;
    ed->last_kill = 0;
    ed->vtop = 0;
    ed->final_frame = false;
    ed->br_dirty = true;
    ed->in_paste = false;
    ed->status = TERM_ED_CONTINUE;
    size_t pl = strlen(prompt);
    if (pl >= sizeof(ed->prompt)) pl = sizeof(ed->prompt) - 1;
    memcpy(ed->prompt, prompt, pl);
    ed->prompt[pl] = '\0';
    ed->prompt_w = (int)pl;
    ed->lay_cols = -1;
    TermScreen* s = &ed->scr;
    if (start_col > 0) {
        term_screen_reset_area(s, cursor_abs_row, cursor_abs_row >= 0);
        term_screen_raw(s, "\r\n", 2);
        if (cursor_abs_row >= 0 && cursor_abs_row < s->rows - 1) cursor_abs_row++;
    }
    term_screen_reset_area(s, cursor_abs_row, cursor_abs_row >= 0);
    s->cur_col = s->cols;                 /* start with \r: the column is not trusted */
    render(ed);
    term_screen_commit(s);
}

TermEdStatus term_editor_feed(TermEditor* ed, const char* bytes, size_t n) {
    if (ed->status != TERM_ED_CONTINUE) return ed->status;
    if (n && !bytes_append(&ed->inq, bytes, n)) return ed->status;
    process(ed, false);
    if (ed->status == TERM_ED_CONTINUE) render(ed);
    term_screen_commit(&ed->scr);
    return ed->status;
}

bool term_editor_pending(const TermEditor* ed) {
    return ed->status == TERM_ED_CONTINUE && ed->inq.len > 0;
}

TermEdStatus term_editor_timeout(TermEditor* ed) {
    if (ed->status != TERM_ED_CONTINUE) return ed->status;
    process(ed, true);
    if (ed->status == TERM_ED_CONTINUE) render(ed);
    term_screen_commit(&ed->scr);
    return ed->status;
}

void term_editor_resize(TermEditor* ed, int cols, int rows) {
    TermScreen* s = &ed->scr;
    term_screen_move(s, 0, 0);
    term_screen_raw(s, "\r\033[J", 4);
    int top = s->top_abs;
    if (!term_screen_resize(s, cols, rows)) return;
    term_screen_reset_area(s, top, false);
    ed->lay_cols = -1;
    if (ed->status == TERM_ED_CONTINUE) render(ed);
    term_screen_commit(s);
}

void term_editor_set_size(TermEditor* ed, int cols, int rows) {
    if (cols == ed->scr.cols && rows == ed->scr.rows) return;
    if (term_screen_resize(&ed->scr, cols, rows)) ed->lay_cols = -1;
}

void term_editor_redraw(TermEditor* ed, bool clear) {
    TermScreen* s = &ed->scr;
    if (clear) {
        term_screen_raw(s, "\033[H\033[2J", 7);
        term_screen_reset_area(s, 0, true);
    } else {
        term_screen_reset_area(s, -1, false);
        s->cur_col = s->cols;
        ed->status = TERM_ED_CONTINUE;
        ed->final_frame = false;
    }
    render(ed);
    term_screen_commit(s);
}

const char* term_editor_text(const TermEditor* ed, size_t* len) {
    if (len) *len = ed->len;
    return ed->text;
}

size_t term_editor_write_count(const TermEditor* ed) { return ed->scr.writes; }
size_t term_editor_cursor(const TermEditor* ed) { return ed->cur; }
