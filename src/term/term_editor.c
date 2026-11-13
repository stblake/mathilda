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
    K_WLEFT, K_WRIGHT, K_ESC, K_PASTE, K_ALT_BS, K_ALT, K_CTRL, K_CPR
};

typedef struct {
    int    code;
    char   b[4];     /* K_CHAR: the UTF-8 bytes; K_ALT / K_CTRL: the letter */
    size_t n;
    int    p1, p2;   /* K_CPR: the reported row and column (1-based)      */
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

    /* The completion popup (Tab). */
    bool            pop_open;
    bool            pop_inline;              /* too small to draw: Tab cycles */
    TermCompletions pop_c;
    size_t*         pop_vis;                 /* indices of pop_c.items shown  */
    size_t          pop_nvis, pop_viscap;
    size_t          pop_sel, pop_top;
    size_t          pop_word;                /* word start in the buffer      */
    size_t          pop_qlen;                /* word length when queried      */

    /* Ctrl-R history search. */
    bool            srch_open;
    Bytes           srch_q;
    size_t*         srch_hits;               /* history indices, newest first */
    size_t          srch_nhits, srch_hitcap;
    size_t          srch_sel, srch_top;
    Bytes           srch_saved;              /* the buffer before the search  */
    size_t          srch_saved_cur;
    uint32_t        srch_hash[1024];         /* dedupe table: hash, 0 = free  */
    size_t          srch_hidx[1024];

    TermEdStatus status;

    bool            pop_slash;               /* the popup lists slash commands */
    size_t          nredo;                   /* ops[nops, nredo) can be redone */
    bool            cpr_wait, cpr_seen;      /* cursor-position reply          */

    /* The picker (term_editor_pick_begin). */
    bool            pick_open;
    TermPick        pick;
    size_t          pick_sel, pick_top;
};

static const char SRCH_LABEL[] = "search: ";

/* Defined with the popup, used by rendering and finish(). */
static void pop_close(TermEditor* ed);
static void pop_draw(TermEditor* ed, int crow, int wcol);
static void srch_draw(TermEditor* ed);
static int  srch_rows(const TermEditor* ed);
static void srch_cancel(TermEditor* ed);
static int  put_text(TermScreen* s, int row, int c, const char* t, size_t n, int maxw, int style);

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
    ed->nredo = 0;
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
    ed->nredo = 0;                       /* a new edit ends the redo chain */
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

/* Defined with the layout code. */
static void gutters(const TermEditor* ed, int* g0, int* gn);
static uint32_t walk_line(const TermEditor* ed, size_t l, int g0, int gn, size_t upto,
                          int* crow, int* ccol);

/* The byte offset on visual row `vrow` of line l nearest to screen column
 * `goal` (not past it): the rule of pos_at_col, per wrapped row. */
static size_t pos_at_visual(const TermEditor* ed, size_t l, int vrow, int goal, int g0, int gn) {
    int cols = ed->scr.cols;
    int avail = (l == 0 ? cols - g0 : cols - gn), availn = cols - gn;
    if (avail < 1) avail = 1;
    if (availn < 1) availn = 1;
    int gut = (l == 0 ? g0 : gn);
    size_t i = ed->lstart[l], e = line_end(ed, l), last = i;
    int row = 0, col = 0;
    while (i < e) {
        uint32_t cp;
        size_t j = term_utf8_next(ed->text, e, i, &cp);
        int w = char_cols(cp);
        if (col + w > avail) { row++; col = 0; avail = availn; gut = gn; }
        if (row > vrow) return last;
        if (row == vrow) {
            if (gut + col + w > goal) return i;
            last = i;
        }
        col += w;
        i = j;
    }
    if (col >= avail) row++;                  /* the end sits on a row of its own */
    return row <= vrow ? e : last;
}

/* Up/Down: to the neighbouring visual row of a wrapped line first, then to
 * the next or previous line (its last / first row), then history. The goal
 * is a screen column, kept across a run of vertical moves. */
static void move_vertical(TermEditor* ed, int dir) {
    size_t l = line_of(ed, ed->cur);
    int g0, gn;
    gutters(ed, &g0, &gn);
    int crow = 0, ccol = 0;
    int nrows = (int)walk_line(ed, l, g0, gn, ed->cur, &crow, &ccol);
    int scol = ccol + ((l == 0 && crow == 0) ? g0 : gn);
    if (!ed->goal_ok) { ed->goal = scol; ed->goal_ok = true; }
    if (dir < 0 && crow > 0) {
        ed->cur = pos_at_visual(ed, l, crow - 1, ed->goal, g0, gn);
    } else if (dir > 0 && crow + 1 < nrows) {
        ed->cur = pos_at_visual(ed, l, crow + 1, ed->goal, g0, gn);
    } else if (dir < 0) {
        if (l == 0) { history_up(ed); return; }
        int pr = (int)walk_line(ed, l - 1, g0, gn, SIZE_MAX, NULL, NULL);
        ed->cur = pos_at_visual(ed, l - 1, pr - 1, ed->goal, g0, gn);
    } else {
        if (l + 1 >= ed->nlines) { history_down(ed); return; }
        ed->cur = pos_at_visual(ed, l + 1, 0, ed->goal, g0, gn);
    }
    ed->goal_keep = true;
}

/* ---- undo -------------------------------------------------------------------- */

static void undo(TermEditor* ed) {
    if (ed->nops == 0) return;
    if (ed->nredo < ed->nops) ed->nredo = ed->nops;
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

/* Re-apply the group undo() last took back. The undone operations and their
 * text are still in ops[] and the arena (undo only shortened them), until a
 * new edit records over them. */
static void redo(TermEditor* ed) {
    if (ed->nredo <= ed->nops) return;
    unsigned g = ed->ops[ed->nops].group;
    size_t k = ed->nops;
    ed->undoing = true;
    while (k < ed->nredo && ed->ops[k].group == g) {
        const UndoOp* op = &ed->ops[k];
        if (ed_replace(ed, op->pos, op->del_len, ed->arena.p + op->ins_off, op->ins_len, EK_OTHER))
            ed->cur = op->pos + op->ins_len;
        k++;
    }
    ed->undoing = false;
    ed->nops = k;
    ed->arena.len = ed->ops[k - 1].ins_off + ed->ops[k - 1].ins_len;
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

static void pick_render(TermEditor* ed);

static void render(TermEditor* ed) {
    TermScreen* s = &ed->scr;
    if (ed->pick_open) { pick_render(ed); return; }
    if (ed->srch_open && !ed->final_frame) {
        ed->vtop = 0;
        term_screen_begin_frame(s, srch_rows(ed));
        srch_draw(ed);
        int qc = (int)(sizeof(SRCH_LABEL) - 1)
               + (int)term_text_width(ed->srch_q.p ? ed->srch_q.p : "", ed->srch_q.len);
        term_screen_render(s, 0, qc < s->cols ? qc : s->cols - 1);
        return;
    }
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
    if (ed->pop_open && !ed->pop_inline && !ed->final_frame) {
        /* The column of the word's start, when it is on the cursor's row. */
        int wrow = 0, wcol = 0;
        size_t wpos = ed->pop_word < ed->lstart[cl] ? ed->lstart[cl] : ed->pop_word;
        (void)walk_line(ed, cl, g0, gn, wpos, &wrow, &wcol);
        wcol += (cl == 0 && wrow == 0) ? g0 : gn;
        if (wrow != crow) wcol = ccol;
        pop_draw(ed, cvrow - ed->vtop, wcol);
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
    pop_close(ed);
    if (ed->srch_open) srch_cancel(ed);
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

/* ---- Ctrl-R: history search ------------------------------------------------------- */

#define SRCH_MAX_ROWS 12

static uint32_t str_hash(const char* s) {
    uint32_t h = 2166136261u;
    for (; *s; s++) h = (h ^ (unsigned char)*s) * 16777619u;
    return h ? h : 1;
}

/* Matching entries, newest first, each distinct text once. */
static void srch_run(TermEditor* ed) {
    size_t n = hist_count(ed);
    ed->srch_nhits = 0;
    ed->srch_sel = ed->srch_top = 0;
    if (n > ed->srch_hitcap) {
        size_t nc = ed->srch_hitcap ? ed->srch_hitcap : 256;
        while (nc < n) nc *= 2;
        size_t* nh = ed_realloc(ed->srch_hits, nc * sizeof(size_t));
        if (!nh) return;
        ed->srch_hits = nh;
        ed->srch_hitcap = nc;
    }
    memset(ed->srch_hash, 0, sizeof(ed->srch_hash));
    const char* q = ed->srch_q.p ? ed->srch_q.p : "";
    for (size_t k = n; k > 0; k--) {
        const char* e = term_history_get(ed->cfg.history, k - 1);
        if (*q && !strstr(e, q)) continue;
        uint32_t h = str_hash(e);
        size_t slot = h & 1023u, probes = 0;
        bool dup = false;
        while (ed->srch_hash[slot] && probes < 1024) {
            if (ed->srch_hash[slot] == h
                && strcmp(term_history_get(ed->cfg.history, ed->srch_hidx[slot]), e) == 0) { dup = true; break; }
            slot = (slot + 1) & 1023u;
            probes++;
        }
        if (dup) continue;
        if (probes < 1024) { ed->srch_hash[slot] = h; ed->srch_hidx[slot] = k - 1; }
        ed->srch_hits[ed->srch_nhits++] = k - 1;
    }
}

static void srch_start(TermEditor* ed) {
    if (!ed->cfg.history) return;
    pop_close(ed);
    bytes_set(&ed->srch_saved, ed->text, ed->len);
    ed->srch_saved_cur = ed->cur;
    bytes_set(&ed->srch_q, "", 0);
    ed->srch_open = true;
    srch_run(ed);
}

static void srch_close(TermEditor* ed) { ed->srch_open = false; }

static void srch_cancel(TermEditor* ed) {
    if (!ed->srch_open) return;
    srch_close(ed);
    if (ed->len != ed->srch_saved.len || memcmp(ed->text, ed->srch_saved.p, ed->len) != 0)
        set_buffer(ed, ed->srch_saved.p, ed->srch_saved.len);
    ed->cur = ed->srch_saved_cur <= ed->len ? ed->srch_saved_cur : ed->len;
}

static void srch_accept(TermEditor* ed) {
    srch_close(ed);
    if (ed->srch_nhits == 0) return;
    const char* e = term_history_get(ed->cfg.history, ed->srch_hits[ed->srch_sel]);
    set_buffer(ed, e, strlen(e));
    ed->last_kind = EK_NONE;
    ed->hbrowse = false;
}

/* Keys while searching. Returns true when the key was used. */
static bool srch_key(TermEditor* ed, const Key* k) {
    switch (k->code) {
    case K_CHAR:
        bytes_append(&ed->srch_q, k->b, k->n);
        srch_run(ed);
        return true;
    case K_BS:
        if (ed->srch_q.len) {
            ed->srch_q.len = term_utf8_prev(ed->srch_q.p, ed->srch_q.len);
            ed->srch_q.p[ed->srch_q.len] = '\0';
            srch_run(ed);
        }
        return true;
    case K_UP:
        if (ed->srch_sel + 1 < ed->srch_nhits) ed->srch_sel++;
        return true;
    case K_DOWN:
        if (ed->srch_sel > 0) ed->srch_sel--;
        return true;
    case K_ENTER: case K_TAB:
        srch_accept(ed);
        return true;
    case K_ESC:
        srch_cancel(ed);
        return true;
    case K_CTRL:
        switch (k->b[0]) {
        case 'r': case 'p': if (ed->srch_sel + 1 < ed->srch_nhits) ed->srch_sel++; return true;
        case 's': case 'n': if (ed->srch_sel > 0) ed->srch_sel--; return true;
        case 'g': srch_cancel(ed); return true;
        case 'c': srch_cancel(ed); return false;      /* then Ctrl-C clears */
        default: break;
        }
        break;
    default: break;
    }
    srch_accept(ed);            /* anything else: take the entry, then act */
    return false;
}

/* Rows of entry `e` in the list. */
static int entry_rows(const char* e) {
    return 1 + (int)count_nl(e, strlen(e));
}

static int srch_rows(const TermEditor* ed) {
    int maxr = ed->scr.rows - 1 < SRCH_MAX_ROWS ? ed->scr.rows - 1 : SRCH_MAX_ROWS;
    if (maxr < 1) maxr = 1;
    int r = 0;
    for (size_t i = 0; i < ed->srch_nhits && r < maxr; i++)
        r += entry_rows(term_history_get(ed->cfg.history, ed->srch_hits[i]));
    if (ed->srch_nhits == 0) r = 1;
    return 1 + (r < maxr ? r : maxr);
}

/* One line of an entry, coloured, clipped with an ellipsis. */
static void srch_line(TermEditor* ed, int row, int col, const char* t, size_t n,
                      TermHlState* st, size_t qpos, size_t qlen) {
    TermScreen* s = &ed->scr;
    uint8_t sty[512];
    size_t ln = n < sizeof(sty) ? n : sizeof(sty);
    *st = term_hl_line(t, ln, *st, sty, ed->cfg.classify, ed->cfg.classify_ctx);
    int end = s->cols;
    size_t i = 0;
    while (i < ln) {
        uint32_t cp;
        size_t j = term_utf8_next(t, ln, i, &cp);
        int w = term_cp_width(cp);
        if (w < 0) { cp = '?'; w = 1; }
        if (w == 0) { i = j; continue; }
        if (col + w > end - 1 && j < ln) { term_screen_put(s, row, col, "\xE2\x80\xA6", 3, 1, TERM_HL_DIM); return; }
        int stl = sty[i] == TERM_HL_BRACKET ? TERM_HL_PLAIN : sty[i];
        if (i >= qpos && i < qpos + qlen) stl = TERM_HL_BRACKET_MATCH;
        if (term_cp_width(cp) < 0) term_screen_put(s, row, col, "?", 1, 1, stl);
        else term_screen_put(s, row, col, t + i, j - i, w, stl);
        col += w;
        i = j;
    }
}

static void srch_draw(TermEditor* ed) {
    TermScreen* s = &ed->scr;
    int maxr = s->frame_rows - 1;
    /* Keep the selected entry in view. */
    if (ed->srch_sel < ed->srch_top) ed->srch_top = ed->srch_sel;
    for (;;) {
        int r = 0;
        for (size_t i = ed->srch_top; i <= ed->srch_sel && i < ed->srch_nhits; i++)
            r += entry_rows(term_history_get(ed->cfg.history, ed->srch_hits[i]));
        if (r <= maxr || ed->srch_top >= ed->srch_sel) break;
        ed->srch_top++;
    }
    int c = term_screen_puts(s, 0, 0, SRCH_LABEL, TERM_HL_PROMPT);
    const char* q = ed->srch_q.p ? ed->srch_q.p : "";
    put_text(s, 0, c, q, ed->srch_q.len, s->cols - c - 1, TERM_HL_PLAIN);
    if (ed->srch_nhits == 0) {
        term_screen_puts(s, 1, 2, "(no match)", TERM_HL_DIM);
        return;
    }
    int row = 1;
    for (size_t i = ed->srch_top; i < ed->srch_nhits && row <= maxr; i++) {
        bool sel = i == ed->srch_sel;
        const char* e = term_history_get(ed->cfg.history, ed->srch_hits[i]);
        const char* qp = *q ? strstr(e, q) : NULL;
        size_t qoff = qp ? (size_t)(qp - e) : SIZE_MAX;
        TermHlState st = term_hl_state_initial();
        size_t a = 0, n = strlen(e);
        bool first = true;
        while (row <= maxr) {
            const char* nl = memchr(e + a, '\n', n - a);
            size_t b = nl ? (size_t)(nl - e) : n;
            if (sel) term_screen_style_range(s, row, 0, 2, TERM_HL_POPUP_SEL);
            term_screen_puts(s, row, 0, first ? (sel ? "> " : "  ") : "  ", sel ? TERM_HL_POPUP_SEL : TERM_HL_PLAIN);
            size_t qp_line = (qoff >= a && qoff < b) ? qoff - a : SIZE_MAX;
            srch_line(ed, row, first ? 2 : 4, e + a, b - a, &st, qp_line, ed->srch_q.len);
            row++;
            first = false;
            if (!nl) break;
            a = b + 1;
            if (row > maxr) {                       /* cut: mark it */
                term_screen_puts(s, maxr, s->cols - 2, "\xE2\x80\xA6", TERM_HL_DIM);
                break;
            }
        }
    }
}

/* ---- completion popup ------------------------------------------------------------ */

#define POP_MAX_ROWS 8

static void pop_close(TermEditor* ed) {
    ed->pop_slash = false;
    if (!ed->pop_open) return;
    term_completions_free(&ed->pop_c);
    ed->pop_open = false;
    ed->pop_inline = false;
    ed->pop_nvis = 0;
}

static bool pop_query(TermEditor* ed) {
    term_completions_free(&ed->pop_c);
    if (ed->pop_slash) {
        if (ed->cfg.slash_complete) ed->cfg.slash_complete(ed->text, ed->cur, &ed->pop_c);
    } else if (ed->cfg.complete) {
        ed->cfg.complete(ed->text, ed->cur, &ed->pop_c);
    } else {
        term_complete(ed->text, ed->cur, &ed->pop_c);
    }
    ed->pop_word = ed->pop_c.word_start;
    ed->pop_qlen = ed->cur - ed->pop_c.word_start;
    if (ed->pop_c.count > ed->pop_viscap) {
        size_t nc = ed->pop_viscap ? ed->pop_viscap : 64;
        while (nc < ed->pop_c.count) nc *= 2;
        size_t* nv = ed_realloc(ed->pop_vis, nc * sizeof(size_t));
        if (!nv) return false;
        ed->pop_vis = nv;
        ed->pop_viscap = nc;
    }
    return ed->pop_c.count > 0;
}

/* Recompute the visible candidates for the word now typed, keeping the
 * selected item when it is still among them. */
static void pop_filter(TermEditor* ed) {
    size_t keep = ed->pop_nvis ? ed->pop_vis[ed->pop_sel] : SIZE_MAX;
    size_t wlen = ed->cur - ed->pop_word;
    const char* w = ed->text + ed->pop_word;
    ed->pop_nvis = 0;
    ed->pop_sel = 0;
    for (size_t i = 0; i < ed->pop_c.count; i++) {
        if (strncmp(ed->pop_c.items[i], w, wlen) != 0) continue;
        if (i == keep) ed->pop_sel = ed->pop_nvis;
        ed->pop_vis[ed->pop_nvis++] = i;
    }
    if (ed->pop_nvis == 0) { pop_close(ed); return; }
    /* A single candidate already typed out in full: nothing left to offer. */
    if (ed->pop_nvis == 1 && strlen(ed->pop_c.items[ed->pop_vis[0]]) == wlen) pop_close(ed);
}

static void pop_replace_word(TermEditor* ed, const char* item) {
    size_t n = strlen(item);
    if (ed_replace(ed, ed->pop_word, ed->cur - ed->pop_word, item, n, EK_OTHER))
        ed->cur = ed->pop_word + n;
}

static void pop_accept(TermEditor* ed) {
    if (ed->pop_nvis) pop_replace_word(ed, ed->pop_c.items[ed->pop_vis[ed->pop_sel]]);
    pop_close(ed);
}

static void pop_move(TermEditor* ed, int d) {
    if (!ed->pop_nvis) return;
    if (d > 0) ed->pop_sel = (ed->pop_sel + 1) % ed->pop_nvis;
    else ed->pop_sel = ed->pop_sel ? ed->pop_sel - 1 : ed->pop_nvis - 1;
    if (ed->pop_inline) pop_replace_word(ed, ed->pop_c.items[ed->pop_vis[ed->pop_sel]]);
}

/* Tab with the popup closed: a single candidate is inserted at once; several
 * insert their common prefix and open the popup (or, on a terminal too small
 * for one, start inline cycling); none at the start of a line indents. */
static void do_complete(TermEditor* ed) {
    if (!pop_query(ed)) {
        size_t a = ed->lstart[line_of(ed, ed->cur)];
        bool blank = true;
        for (size_t i = a; i < ed->cur; i++) if (ed->text[i] != ' ') blank = false;
        if (blank) insert(ed, "  ", 2, EK_OTHER);
        term_completions_free(&ed->pop_c);
        return;
    }
    TermCompletions* c = &ed->pop_c;
    if (c->count == 1) {
        pop_replace_word(ed, c->items[0]);
        term_completions_free(c);
        return;
    }
    size_t lcp = strlen(c->items[0]);
    for (size_t i = 1; i < c->count; i++) {
        size_t k = 0;
        while (k < lcp && c->items[i][k] == c->items[0][k]) k++;
        lcp = k;
    }
    size_t wlen = ed->cur - ed->pop_word;
    if (lcp > wlen && ed_replace(ed, ed->pop_word, wlen, c->items[0], lcp, EK_OTHER))
        ed->cur = ed->pop_word + lcp;
    ed->pop_open = true;
    ed->pop_nvis = 0;
    ed->pop_top = 0;
    pop_filter(ed);
    if (ed->pop_open && (ed->scr.cols < 16 || ed->scr.rows < 5)) {
        ed->pop_inline = true;
        ed->pop_sel = 0;
        pop_replace_word(ed, c->items[ed->pop_vis[0]]);
    }
}

/* The slash-command menu: open it over the commands matching the buffer
 * ("/" or "/na"). */
static void slash_open(TermEditor* ed) {
    pop_close(ed);
    ed->pop_slash = true;
    if (!pop_query(ed)) { term_completions_free(&ed->pop_c); ed->pop_slash = false; return; }
    ed->pop_open = true;
    ed->pop_nvis = 0;
    ed->pop_top = 0;
    ed->pop_sel = 0;
    pop_filter(ed);
}

/* Should the menu be (re)opened after this edit? A command name being typed:
 * the buffer is '/' plus letters and digits, the cursor at its end. */
static bool slash_wanted(const TermEditor* ed) {
    if (!ed->cfg.slash_complete || ed->len == 0 || ed->text[0] != '/' || ed->cur != ed->len) return false;
    for (size_t i = 1; i < ed->len; i++) {
        char c = ed->text[i];
        if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9'))) return false;
    }
    return true;
}

/* Keys while the popup is open. Returns true when the key was used. */
static bool pop_key(TermEditor* ed, const Key* k) {
    if (ed->pop_slash) {
        if (k->code == K_ENTER || k->code == K_FORCE_ENTER) {   /* run the command */
            pop_accept(ed);
            finish(ed, TERM_ED_ACCEPT);
            return true;
        }
        if (k->code == K_TAB) {                                  /* complete the name */
            pop_accept(ed);
            insert(ed, " ", 1, EK_INSERT);
            return true;
        }
    }
    switch (k->code) {
    case K_TAB:   pop_move(ed, +1); return true;
    case K_BTAB:  pop_move(ed, -1); return true;
    case K_DOWN:  if (ed->pop_inline) break; pop_move(ed, +1); return true;
    case K_UP:    if (ed->pop_inline) break; pop_move(ed, -1); return true;
    case K_ENTER:
        if (ed->pop_inline) { pop_close(ed); return true; }
        pop_accept(ed);
        return true;
    case K_ESC:   pop_close(ed); return true;
    case K_CTRL:
        if (k->b[0] == 'n' && !ed->pop_inline) { pop_move(ed, +1); return true; }
        if (k->b[0] == 'p' && !ed->pop_inline) { pop_move(ed, -1); return true; }
        if (k->b[0] == 'g') { pop_close(ed); return true; }
        break;
    case K_CHAR:
        if (!ed->pop_inline && word_byte(k->b[0])) {
            insert(ed, k->b, k->n, EK_INSERT);
            pop_filter(ed);
            return true;
        }
        break;
    case K_BS:
        if (!ed->pop_inline && ed->cur > ed->pop_word) {
            size_t p = prev_char(ed, ed->cur);
            if (ed_replace(ed, p, ed->cur - p, "", 0, EK_DEL_BACK)) ed->cur = p;
            if (ed->cur == ed->pop_word) { pop_close(ed); return true; }
            if (ed->cur - ed->pop_word < ed->pop_qlen) {
                /* Shorter than the query: candidates may have been missed. */
                if (!pop_query(ed)) { pop_close(ed); return true; }
                ed->pop_nvis = 0;
            }
            pop_filter(ed);
            return true;
        }
        break;
    default: break;
    }
    pop_close(ed);              /* any other key closes it, then acts */
    return false;
}

/* The name a candidate describes: its leading symbol characters. */
static size_t item_name_len(const char* item) {
    size_t n = 0;
    while (item[n] && term_lex_is_word_char(item[n])) n++;
    return n;
}

/* Put text into a row from column c, at most `maxw` columns, control
 * characters as spaces; returns the column after it. An ellipsis marks a
 * cut. */
static int put_text(TermScreen* s, int row, int c, const char* t, size_t n, int maxw, int style) {
    if (maxw <= 0) return c;
    bool cut = (int)term_text_width(t, n) > maxw;
    int end = c + (cut ? maxw - 1 : maxw);
    size_t i = 0;
    while (i < n) {
        uint32_t cp;
        size_t j = term_utf8_next(t, n, i, &cp);
        int w = term_cp_width(cp);
        if (w != 0) {
            if (w < 0) w = 1;
            if (c + w > end) break;
            if (term_cp_width(cp) < 0) term_screen_put(s, row, c, " ", 1, 1, style);
            else term_screen_put(s, row, c, t + i, j - i, w, style);
            c += w;
        }
        i = j;
    }
    if (cut) term_screen_put(s, row, c++, "\xE2\x80\xA6", 3, 1, style);
    return c;
}

/* Lay out and draw the popup for a cursor on frame row `crow`, word column
 * `wcol`; may grow the frame. Called with the buffer already drawn. */
static void pop_draw(TermEditor* ed, int crow, int wcol) {
    TermScreen* s = &ed->scr;
    int cols = s->cols;
    int h = (int)(ed->pop_nvis < POP_MAX_ROWS ? ed->pop_nvis : POP_MAX_ROWS);
    /* The command menu shows each command's summary beside it instead of a
     * description row. */
    bool slash = ed->pop_slash && ed->cfg.slash_doc && cols >= 30;
    bool doc_row = !slash && cols >= 40 && ed->cfg.doc;
    const char* doc = NULL;
    size_t doc_len = 0;
    if (doc_row) {
        const char* item = ed->pop_c.items[ed->pop_vis[ed->pop_sel]];
        doc = ed->cfg.doc(item, item_name_len(item), ed->cfg.doc_ctx);
        if (doc) { const char* nl = strchr(doc, '\n'); doc_len = nl ? (size_t)(nl - doc) : strlen(doc); }
    }
    bool counter = ed->pop_nvis > (size_t)h;
    if (!doc_row) doc_row = counter && cols >= 20;
    int ph = h + (doc_row ? 1 : 0);
    if (ph > s->rows - 1) { h = s->rows - 1 - (doc_row ? 1 : 0); ph = s->rows - 1; }
    if (h < 1) return;

    /* Above the cursor when the rows below are not on screen but the area
     * has room above; otherwise below (the screen scrolls if it must). */
    int below = term_screen_rows_below(s, crow);
    int top;
    if (below >= ph || crow < ph) top = crow + 1;
    else top = crow - ph;
    if (top + ph > s->rows) top = s->rows - ph;
    if (top + ph > s->frame_rows) term_screen_extend_frame(s, top + ph);

    /* Keep the selection in the visible window. */
    if (ed->pop_sel < ed->pop_top) ed->pop_top = ed->pop_sel;
    if (ed->pop_sel >= ed->pop_top + (size_t)h) ed->pop_top = ed->pop_sel - (size_t)h + 1;
    if (ed->pop_top + (size_t)h > ed->pop_nvis) ed->pop_top = ed->pop_nvis - (size_t)h;

    int iw = 0;
    for (size_t i = 0; i < ed->pop_nvis; i++) {
        int w = (int)term_text_width(ed->pop_c.items[ed->pop_vis[i]], strlen(ed->pop_c.items[ed->pop_vis[i]]));
        if (w > iw) iw = w;
    }
    if (iw > cols - 3) iw = cols - 3;
    if (iw > 48) iw = 48;
    int pw = iw + 2;
    int dw = 0;                                   /* the summary column (slash) */
    if (slash) {
        for (size_t i = 0; i < ed->pop_nvis; i++) {
            const char* it = ed->pop_c.items[ed->pop_vis[i]];
            const char* d = ed->cfg.slash_doc(it, strlen(it), ed->cfg.slash_doc_ctx);
            int w = d ? (int)term_text_width(d, strlen(d)) : 0;
            if (w > dw) dw = w;
        }
        if (dw > 0) {
            int room = (cols < 80 ? cols : 80) - (wcol > 1 ? wcol - 1 : 0) - pw - 2;
            if (dw > room) dw = room;
            if (dw < 8) dw = 0;
            if (dw) pw += 2 + dw;
        }
    }
    int x = wcol - 1;
    if (x + pw > cols) x = cols - pw;
    if (x < 0) x = 0;

    bool plain = !s->color;
    for (int r = 0; r < h; r++) {
        size_t idx = ed->pop_top + (size_t)r;
        bool sel = idx == ed->pop_sel;
        int st = sel ? TERM_HL_POPUP_SEL : TERM_HL_POPUP;
        int row = top + r;
        for (int c = x; c < x + pw; c++) term_screen_put(s, row, c, " ", 1, 1, st);
        if (plain && sel) term_screen_put(s, row, x, ">", 1, 1, st);
        const char* item = ed->pop_c.items[ed->pop_vis[idx]];
        put_text(s, row, x + 1, item, strlen(item), iw, st);
        if (dw) {
            const char* d = ed->cfg.slash_doc(item, strlen(item), ed->cfg.slash_doc_ctx);
            if (d) put_text(s, row, x + 1 + iw + 2, d, strlen(d), dw, sel ? st : TERM_HL_POPUP_DOC);
        }
    }
    if (doc_row) {
        int row = top + h;
        int dw = cols - x;
        if (dw > 72) dw = 72;
        if (dw < pw) dw = pw;
        if (x + dw > cols) dw = cols - x;
        for (int c = x; c < x + dw; c++) term_screen_put(s, row, c, " ", 1, 1, TERM_HL_POPUP_DOC);
        char cnt[32] = "";
        int cl = 0;
        if (counter) cl = snprintf(cnt, sizeof(cnt), " %zu/%zu", ed->pop_sel + 1, ed->pop_nvis);
        if (doc && doc_len) put_text(s, row, x + 1, doc, doc_len, dw - 2 - cl, TERM_HL_POPUP_DOC);
        if (cl) term_screen_puts(s, row, x + dw - cl - 1, cnt, TERM_HL_POPUP_DOC);
    }
}

/* ---- the picker ------------------------------------------------------------------ */

static void pick_select(TermEditor* ed, size_t i) {
    if (i >= ed->pick.count || i == ed->pick_sel) return;
    ed->pick_sel = i;
    if (ed->pick.on_move) {
        ed->pick.on_move(i, ed->pick.ctx);
        term_screen_invalidate(&ed->scr);    /* the styles' colours changed */
    }
}

static void pick_move(TermEditor* ed, int d) {
    size_t n = ed->pick.count;
    if (n == 0) return;
    pick_select(ed, d > 0 ? (ed->pick_sel + 1) % n : (ed->pick_sel ? ed->pick_sel - 1 : n - 1));
}

/* End the picker: erase it, show the cursor again. */
static void pick_finish(TermEditor* ed, TermEdStatus st) {
    TermScreen* s = &ed->scr;
    ed->pick_open = false;
    ed->status = st;
    term_screen_reset_style(s);
    term_screen_move(s, 0, 0);
    term_screen_raw(s, "\r\033[J\033[?25h", 10);
    term_screen_reset_area(s, s->top_abs, s->top_known);
}

static void pick_key(TermEditor* ed, const Key* k) {
    switch (k->code) {
    case K_UP: case K_BTAB: case K_LEFT:   pick_move(ed, -1); break;
    case K_DOWN: case K_TAB: case K_RIGHT: pick_move(ed, +1); break;
    case K_HOME: case K_PGUP: pick_select(ed, 0); break;
    case K_END: case K_PGDN:  if (ed->pick.count) pick_select(ed, ed->pick.count - 1); break;
    case K_ENTER: case K_FORCE_ENTER: pick_finish(ed, TERM_ED_ACCEPT); break;
    case K_ESC: pick_finish(ed, TERM_ED_INTERRUPT); break;
    case K_PASTE: ed->in_paste = true; ed->paste_cr = false; break;   /* ignored */
    case K_CTRL:
        switch (k->b[0]) {
        case 'p': pick_move(ed, -1); break;
        case 'n': pick_move(ed, +1); break;
        case 'c': case 'g': pick_finish(ed, TERM_ED_INTERRUPT); break;
        case 'd': pick_finish(ed, TERM_ED_EOF); break;
        case 'z': ed->status = TERM_ED_SUSPEND; break;
        case 'l': term_screen_raw(&ed->scr, "\033[H\033[2J", 7); term_screen_reset_area(&ed->scr, 0, true); break;
        default: break;
        }
        break;
    case K_CHAR:
        if (k->b[0] == 'k') pick_move(ed, -1);
        else if (k->b[0] == 'j') pick_move(ed, +1);
        else if (k->b[0] == 'q') pick_finish(ed, TERM_ED_INTERRUPT);
        else if (k->b[0] >= '1' && k->b[0] <= '9') pick_select(ed, (size_t)(k->b[0] - '1'));
        break;
    default: break;
    }
}

/* One preview line from column 2 of `row`. */
static void pick_preview(TermEditor* ed, int row, const TermPreviewLine* pv) {
    TermScreen* s = &ed->scr;
    int cols = s->cols;
    int c = 2;
    const char* lab = pv->label ? pv->label : "";
    int lst = pv->kind == TERM_PV_INPUT ? TERM_HL_PROMPT
            : pv->kind == TERM_PV_OUTPUT ? TERM_HL_OUT
            : pv->kind == TERM_PV_STATUS ? TERM_HL_SPIN : TERM_HL_PLAIN;
    c = put_text(s, row, c, lab, strlen(lab), cols - c, lst);
    const char* t = pv->text ? pv->text : "";
    size_t n = strlen(t);
    if (pv->kind == TERM_PV_INPUT || pv->kind == TERM_PV_OUTPUT) {
        uint8_t sty[256];
        if (n > sizeof sty) n = sizeof sty;
        (void)term_hl_line(t, n, term_hl_state_initial(), sty, ed->cfg.classify, ed->cfg.classify_ctx);
        size_t i = 0;
        while (i < n) {
            uint32_t cp;
            size_t j = term_utf8_next(t, n, i, &cp);
            int w = term_cp_width(cp);
            if (w < 1) w = 1;
            if (c + w > cols) break;
            term_screen_put(s, row, c, t + i, j - i, w, sty[i] == TERM_HL_BRACKET ? TERM_HL_PLAIN : sty[i]);
            c += w;
            i = j;
        }
    } else if (pv->kind == TERM_PV_MESSAGE) {
        const char* colon = strstr(t, ": ");
        size_t hn = colon ? (size_t)(colon - t) : 0;
        c = put_text(s, row, c, t, hn, cols - c, TERM_HL_MSG);
        put_text(s, row, c, t + hn, n - hn, cols - c, TERM_HL_PLAIN);
    } else {
        put_text(s, row, c, t, n, cols - c, pv->kind == TERM_PV_STATUS || pv->kind == TERM_PV_NOTE
                                                 ? TERM_HL_DIM : TERM_HL_PLAIN);
    }
}

static void pick_render(TermEditor* ed) {
    TermScreen* s = &ed->scr;
    /* A move can redraw every row (a new theme): synchronized output keeps
     * terminals that support it (DEC mode 2026) from showing a half frame;
     * the others ignore it. */
    term_screen_raw(s, "\033[?2026h", 8);
    int cols = s->cols;
    int maxr = s->rows - 1 > 1 ? s->rows - 1 : 1;
    int n = (int)ed->pick.count;
    int npv = (int)ed->pick.npreview;
    int h = n;
    if (1 + h + (npv ? 1 + npv : 0) > maxr) npv = 0;     /* the preview goes first */
    if (1 + h > maxr) h = maxr - 1;
    if (h < 1) h = 1;
    term_screen_begin_frame(s, 1 + h + (npv ? 1 + npv : 0));

    int c = term_screen_puts(s, 0, 0, ed->pick.title ? ed->pick.title : "", TERM_HL_LABEL);
    static const char HINT[] = "   \xE2\x86\x91\xE2\x86\x93 preview \xC2\xB7 Enter apply \xC2\xB7 Esc cancel";
    if (cols - c > 20) put_text(s, 0, c, HINT, sizeof(HINT) - 1, cols - c - 1, TERM_HL_DIM);

    if (ed->pick_sel < ed->pick_top) ed->pick_top = ed->pick_sel;
    if (ed->pick_sel >= ed->pick_top + (size_t)h) ed->pick_top = ed->pick_sel - (size_t)h + 1;

    int iw = 0, dw = 0;
    for (int i = 0; i < n; i++) {
        int w = (int)term_text_width(ed->pick.items[i], strlen(ed->pick.items[i]));
        if (w > iw) iw = w;
        if (ed->pick.docs && ed->pick.docs[i]) {
            w = (int)term_text_width(ed->pick.docs[i], strlen(ed->pick.docs[i]));
            if (w > dw) dw = w;
        }
    }
    if (iw > 24) iw = 24;
    int pw = 2 + iw + (dw ? 2 + dw : 0) + 1;
    if (pw > cols) pw = cols;
    if (pw > 76) pw = 76;
    bool plain = !s->color;
    for (int r = 0; r < h; r++) {
        size_t idx = ed->pick_top + (size_t)r;
        if (idx >= (size_t)n) break;
        bool sel = idx == ed->pick_sel;
        int st = sel ? TERM_HL_POPUP_SEL : TERM_HL_POPUP;
        int row = 1 + r;
        for (int x = 0; x < pw; x++) term_screen_put(s, row, x, " ", 1, 1, st);
        if (plain && sel) term_screen_put(s, row, 0, ">", 1, 1, st);
        const char* it = ed->pick.items[idx];
        put_text(s, row, 2, it, strlen(it), iw < pw - 2 ? iw : pw - 2, st);
        if (dw && ed->pick.docs[idx] && 2 + iw + 2 < pw - 1) {
            const char* d = ed->pick.docs[idx];
            put_text(s, row, 2 + iw + 2, d, strlen(d), pw - 1 - (2 + iw + 2), sel ? st : TERM_HL_POPUP_DOC);
        }
    }
    for (int i = 0; i < npv; i++) pick_preview(ed, 2 + h + i, &ed->pick.preview[i]);
    term_screen_render(s, 0, c < cols ? c : cols - 1);
    term_screen_raw(s, "\033[?2026l", 8);
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
    case 'r': srch_start(ed); break;
    case 'p': move_vertical(ed, -1); break;
    case 'n': move_vertical(ed, +1); break;
    default: break;
    }
}

static void pick_key(TermEditor* ed, const Key* k);

static void handle_key(TermEditor* ed, const Key* k) {
    if (k->code == K_CPR) {                 /* a cursor-position reply, not a key */
        if (ed->cpr_wait) {
            ed->cpr_wait = false;
            ed->cpr_seen = true;
            ed->scr.top_abs = k->p1 - 1;
            ed->scr.top_known = true;
        }
        return;
    }
    if (ed->pick_open) { pick_key(ed, k); return; }
    ed->edited = false;
    ed->goal_keep = false;
    ed->this_kill = 0;
    size_t l;
    if (ed->srch_open && srch_key(ed, k)) {
        if (!ed->edited) ed->last_kind = EK_NONE;
        ed->last_kill = 0;
        return;
    }
    if (ed->pop_open && pop_key(ed, k)) {
        if (!ed->edited) ed->last_kind = EK_NONE;
        ed->last_kill = 0;
        return;
    }
    switch (k->code) {
    case K_CHAR:
        insert(ed, k->b, k->n, EK_INSERT);
        break;
    case K_ENTER:
        /* A slash command is one line: Return always runs it. */
        if ((ed->cfg.slash_complete && ed->len && ed->text[0] == '/')
            || term_lex_input_complete(ed->text)) finish(ed, TERM_ED_ACCEPT);
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
        case '_': redo(ed); break;
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
    /* Typing a command name (re)opens the slash menu. */
    if (ed->edited && (k->code == K_CHAR || k->code == K_BS) && ed->status == TERM_ED_CONTINUE
        && !ed->pop_open && !ed->srch_open && slash_wanted(ed))
        slash_open(ed);
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
            if (f == 'R' && num > 0 && mod > 0) {         /* ESC [ row ; col R */
                k->code = K_CPR;
                k->p1 = num;
                k->p2 = mod;
            } else if (code != K_NONE) {
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
    term_completions_free(&ed->pop_c);
    free(ed->pop_vis);
    free(ed->srch_q.p); free(ed->srch_saved.p); free(ed->srch_hits);
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
    ed->nredo = 0;
    ed->arena.len = 0;
    ed->last_kind = EK_NONE;
    ed->hbrowse = false;
    ed->pick_open = false;
    ed->goal_ok = false;
    ed->last_kill = 0;
    ed->vtop = 0;
    ed->final_frame = false;
    ed->br_dirty = true;
    ed->in_paste = false;
    pop_close(ed);
    ed->srch_open = false;
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

void term_editor_pick_begin(TermEditor* ed, const TermPick* p, int cursor_abs_row) {
    pop_close(ed);
    ed->srch_open = false;
    ed->pick = *p;
    ed->pick_open = true;
    ed->pick_sel = p->selected < p->count ? p->selected : 0;
    ed->pick_top = 0;
    ed->status = TERM_ED_CONTINUE;
    ed->final_frame = false;
    ed->in_paste = false;
    TermScreen* s = &ed->scr;
    term_screen_reset_area(s, cursor_abs_row, cursor_abs_row >= 0);
    s->cur_col = s->cols;
    term_screen_raw(s, "\033[?25l", 6);           /* the selection is the cursor */
    render(ed);
    term_screen_commit(s);
}

size_t term_editor_pick_index(const TermEditor* ed) { return ed->pick_sel; }

void term_editor_expect_cpr(TermEditor* ed) { ed->cpr_wait = true; }
bool term_editor_cpr_seen(const TermEditor* ed) { return ed->cpr_seen; }
bool term_editor_cpr_pending(const TermEditor* ed) { return ed->cpr_wait; }

const char* term_editor_text(const TermEditor* ed, size_t* len) {
    if (len) *len = ed->len;
    return ed->text;
}

size_t term_editor_write_count(const TermEditor* ed) { return ed->scr.writes; }
size_t term_editor_cursor(const TermEditor* ed) { return ed->cur; }
