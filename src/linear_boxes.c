/* linear_boxes.c -- linear box syntax, \!\( ... \), as input.
 *
 * The reference lets a box structure be written inline: \( ... \) groups box
 * text, \* introduces a box expression, and the optional \! in front asks for
 * the boxes to be read as the expression they display.  FeynCalc's tests use
 * it for powers and fractions copied out of a notebook:
 *
 *     \!\(\*SuperscriptBox[\(mb\), \(2\)]\) - \!\(\*FractionBox[\(1\), \(2\)]\)
 *
 * is mb^2 - 1/2.  Mathilda has no box language, so a group is lowered straight
 * to InputForm text, which the ordinary parser then reads:
 *   - plain text in a group is copied, and a nested \( X \) becomes ( X );
 *   - \*Box[...] is parsed as an expression whose \( X \) arguments are box
 *     strings, and is lowered by box_to_text.
 * The lowering is recursive, so a box inside a box inside a group works.
 */

#include "linear_boxes.h"
#include "parse.h"
#include "expr.h"
#include "print.h"
#include "utf8.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* The box control escapes and their private-use code points (the reference's
 * ToCharacterCode["\!\(\)\*\@\^\&\%\+\/\_\`"]). */
static const struct { char c; uint32_t cp; const char* esc; } LB_ESC[] = {
    { '!', 0xF7C1, "\\!" }, { '(', 0xF7C9, "\\(" }, { ')', 0xF7C0, "\\)" },
    { '*', 0xF7C8, "\\*" }, { '@', 0xF7C2, "\\@" }, { '^', 0xF7C6, "\\^" },
    { '&', 0xF7C7, "\\&" }, { '%', 0xF7C5, "\\%" }, { '+', 0xF7CB, "\\+" },
    { '/', 0xF7CC, "\\/" }, { '_', 0xF7CA, "\\_" }, { '`', 0xF7CD, "\\`" },
};
#define LB_NESC (sizeof LB_ESC / sizeof LB_ESC[0])

uint32_t lb_escape_cp(char c) {
    for (size_t i = 0; i < LB_NESC; i++) if (LB_ESC[i].c == c) return LB_ESC[i].cp;
    return 0;
}

const char* lb_cp_escape(uint32_t cp) {
    for (size_t i = 0; i < LB_NESC; i++) if (LB_ESC[i].cp == cp) return LB_ESC[i].esc;
    return NULL;
}

/* The box control character at p, in either spelling: its ASCII letter
 * ('(' for \( or U+F7C9, ...) with its length in *len; 0 otherwise. */
static char lb_ctl(const char* p, size_t* len) {
    if (p[0] == '\\' && p[1]) {
        for (size_t i = 0; i < LB_NESC; i++)
            if (LB_ESC[i].c == p[1]) { *len = 2; return p[1]; }
        return 0;
    }
    if ((unsigned char)p[0] == 0xEF) {           /* U+F7xx is EF 9E/9F xx */
        size_t n = 0, i = 0;
        while (n < 3 && p[n]) n++;
        uint32_t cp = u8_decode(p, n, &i);
        for (size_t k = 0; k < LB_NESC; k++)
            if (LB_ESC[k].cp == cp) { *len = i; return LB_ESC[k].c; }
    }
    return 0;
}

size_t lb_group_open(const char* p) {
    size_t a = 0, b = 0;
    char c = lb_ctl(p, &a);
    if (c == '(') return a;
    if (c == '!' && lb_ctl(p + a, &b) == '(') return a + b;
    return 0;
}

/* ---- a growable text buffer ---- */
typedef struct { char* s; size_t n, cap; bool oom; } Buf;

static void buf_add(Buf* b, const char* s, size_t n) {
    if (b->oom) return;
    if (b->n + n + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 64;
        while (cap < b->n + n + 1) cap *= 2;
        char* t = realloc(b->s, cap);
        if (!t) { b->oom = true; return; }
        b->s = t; b->cap = cap;
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = '\0';
}
static void buf_str(Buf* b, const char* s) { buf_add(b, s, strlen(s)); }

/* Past the string literal opening at p (p[0] == '"'), or NULL if unclosed. */
static const char* skip_string(const char* p) {
    for (p++; *p; p++) {
        if (*p == '\\' && p[1]) { p++; continue; }
        if (*p == '"') return p + 1;
    }
    return NULL;
}

/* Past the \) matching the group whose content starts at p, with the content's
 * end in *cend; NULL if unbalanced. */
static const char* group_close(const char* p, const char** cend) {
    int depth = 1;
    while (*p) {
        size_t l = 0;
        char c = lb_ctl(p, &l);
        if (c == '(') { depth++; p += l; continue; }
        if (c == ')') {
            if (--depth == 0) { *cend = p; return p + l; }
            p += l; continue;
        }
        if (*p == '"') { p = skip_string(p); if (!p) return NULL; continue; }
        p++;
    }
    return NULL;
}

static bool group_text(const char* p, const char* end, Buf* out);

/* Lower a parsed box expression to InputForm text. */
static bool box_to_text(const Expr* b, Buf* out) {
    if (b->type == EXPR_STRING) { buf_str(out, b->data.string); return true; }
    if (b->type != EXPR_FUNCTION || b->data.function.head->type != EXPR_SYMBOL) {
        char* t = expr_to_string_fullform((Expr*)b);   /* a number or symbol */
        if (!t) return false;
        buf_str(out, t);
        free(t);
        return true;
    }
    const char* h = b->data.function.head->data.symbol.name;
    size_t n = b->data.function.arg_count;
    Expr* const* A = b->data.function.args;
    const char* fn = NULL;
    if (h == SYM_SuperscriptBox && n == 2) fn = "Power";
    else if (h == SYM_SubscriptBox && n == 2) fn = "Subscript";
    else if (strcmp(h, "FractionBox") == 0 && n == 2) fn = "Divide";
    else if (strcmp(h, "SqrtBox") == 0 && n == 1) fn = "Sqrt";
    else if (h == SYM_OverscriptBox && n == 2) fn = "Overscript";
    else if (strcmp(h, "UnderscriptBox") == 0 && n == 2) fn = "Underscript";
    else if (h == SYM_SubsuperscriptBox && n == 3) fn = "Subsuperscript";
    if (fn) {
        buf_str(out, fn);
        buf_str(out, "[");
        for (size_t i = 0; i < n; i++) {
            if (i) buf_str(out, ", ");
            buf_str(out, "(");
            if (!box_to_text(A[i], out)) return false;
            buf_str(out, ")");
        }
        buf_str(out, "]");
        return true;
    }
    if (h == SYM_RowBox && n == 1 && A[0]->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < A[0]->data.function.arg_count; i++)
            if (!box_to_text(A[0]->data.function.args[i], out)) return false;
        return true;
    }
    /* Wrappers that display their first argument. */
    if (n >= 1 && (strcmp(h, "FormBox") == 0 || strcmp(h, "StyleBox") == 0 ||
                   h == SYM_TagBox || strcmp(h, "AdjustmentBox") == 0))
        return box_to_text(A[0], out);
    /* InterpretationBox[boxes, expr, ...] stands for expr. */
    if (strcmp(h, "InterpretationBox") == 0 && n >= 2) {
        char* t = expr_to_string_fullform(A[1]);
        if (!t) return false;
        buf_str(out, "(");
        buf_str(out, t);
        buf_str(out, ")");
        free(t);
        return true;
    }
    return false;
}

/* \*Box[...] at p (just past the \*): parse the box expression, its \( X \)
 * arguments read as the box strings "X", lower it, and return the position
 * after it; NULL on error. */
static const char* star_box(const char* p, Buf* out) {
    Buf src = {0};
    int depth = 0;
    const char* q = p;
    /* the head, then a balanced [...] */
    while (*q && *q != '[' && (*q == '`' || *q == '$' || (*q >= '0' && *q <= '9') ||
                               (*q >= 'A' && *q <= 'Z') || (*q >= 'a' && *q <= 'z'))) q++;
    if (*q != '[' || q == p) { free(src.s); return NULL; }
    buf_add(&src, p, (size_t)(q - p));
    while (*q) {
        size_t l = 0;
        char c = lb_ctl(q, &l);
        if (c == '(') {
            const char* cend = NULL;
            const char* after = group_close(q + l, &cend);
            if (!after) { free(src.s); return NULL; }
            Buf inner = {0};
            if (!group_text(q + l, cend, &inner)) { free(inner.s); free(src.s); return NULL; }
            buf_str(&src, "\"");
            for (size_t i = 0; i < inner.n; i++) {
                if (inner.s[i] == '"' || inner.s[i] == '\\') buf_add(&src, "\\", 1);
                buf_add(&src, inner.s + i, 1);
            }
            buf_str(&src, "\"");
            free(inner.s);
            q = after;
            continue;
        }
        if (c) { free(src.s); return NULL; }      /* a stray \) or \* */
        if (*q == '"') {
            const char* e = skip_string(q);
            if (!e) { free(src.s); return NULL; }
            buf_add(&src, q, (size_t)(e - q));
            q = e;
            continue;
        }
        if (*q == '[' || *q == '{') depth++;
        if (*q == ']' || *q == '}') depth--;
        buf_add(&src, q, 1);
        q++;
        if (depth == 0) break;
    }
    if (depth != 0 || src.oom) { free(src.s); return NULL; }
    Expr* box = parse_expression(src.s);
    free(src.s);
    if (!box) return NULL;
    buf_str(out, "(");
    bool ok = box_to_text(box, out);
    buf_str(out, ")");
    expr_free(box);
    return ok ? q : NULL;
}

/* The InputForm text of the group content [p, end). */
static bool group_text(const char* p, const char* end, Buf* out) {
    while (p < end) {
        size_t l = 0;
        char c = lb_ctl(p, &l);
        if (c == '!') { p += l; continue; }
        if (c == '(') {
            const char* cend = NULL;
            const char* after = group_close(p + l, &cend);
            if (!after || after > end) return false;
            buf_str(out, "(");
            if (!group_text(p + l, cend, out)) return false;
            buf_str(out, ")");
            p = after;
            continue;
        }
        if (c == '*') {
            p = star_box(p + l, out);
            if (!p || p > end) return false;
            continue;
        }
        if (c) return false;
        /* "\ " is a space token in box text (p^\ 2 is RowBox[{p, ^, " ", 2}],
         * which reads as p^2): lower it to a plain space. */
        if (*p == '\\' && p + 1 < end && p[1] == ' ') {
            buf_add(out, " ", 1);
            p += 2;
            continue;
        }
        if (*p == '"') {
            const char* e = skip_string(p);
            if (!e || e > end) return false;
            buf_add(out, p, (size_t)(e - p));
            p = e;
            continue;
        }
        buf_add(out, p, 1);
        p++;
    }
    return !out->oom;
}

char* lb_group_to_input(const char* p, const char** end) {
    size_t open = lb_group_open(p);
    if (!open) return NULL;
    const char* cend = NULL;
    const char* after = group_close(p + open, &cend);
    if (!after) return NULL;
    Buf out = {0};
    buf_str(&out, "(");
    if (!group_text(p + open, cend, &out)) { free(out.s); return NULL; }
    buf_str(&out, ")");
    if (out.oom) { free(out.s); return NULL; }
    *end = after;
    return out.s;
}
