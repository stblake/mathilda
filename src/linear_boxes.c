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
 *     strings, and is lowered by box_to_text;
 *   - the infix box escapes build the same boxes inline and are lowered to
 *     the expressions those boxes stand for:
 *         a\^b     SuperscriptBox     a^b
 *         a\_b     SubscriptBox       Subscript[a, b]
 *         a\_b\%c  SubsuperscriptBox  Subscript[a, b]^c   (also a\^c\%b)
 *         a\&b     OverscriptBox      Overscript[a, b]
 *         a\+b     UnderscriptBox     Underscript[a, b]
 *         a\+b\%c  UnderoverscriptBox Underoverscript[a, b, c]  (a\&c\%b)
 *         a\/b     FractionBox        a/b
 *         \@a      SqrtBox            Sqrt[a]
 *         \@a\%n   RadicalBox         a^(1/n)
 * The lowering is recursive, so a box inside a box inside a group works.
 *
 * The escapes are operators over box OPERANDS: one primary token (a symbol,
 * number, string, slot, bracketed or \( \) group, \*Box) together with its
 * postfix forms (f[x], x[[1]], x', x!, f@x).  Tightest first:
 *     \& \+   then  \_   then  \^   then  \@ (prefix)   then  \/
 * The script operators group to the right (a\^b\^c is a^(b^c)), \/ to the
 * left (a\/b\/c is (a/b)/c), and every one is tighter than the ordinary
 * operators around it: a b\/c d is a (b/c) d, -x\^2 is -(x^2).  A right
 * operand may carry a sign (x\^-1, x\/-y, \@-x).  An operand that is a bare
 * \( \) group contributes its contents (x\_\(i,j\) is Subscript[x, i, j]).
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
    else if (strcmp(h, "SqrtBox") == 0 && n == 1) fn = "Sqrt";
    else if (h == SYM_OverscriptBox && n == 2) fn = "Overscript";
    else if (strcmp(h, "UnderscriptBox") == 0 && n == 2) fn = "Underscript";
    else if (strcmp(h, "UnderoverscriptBox") == 0 && n == 3) fn = "Underoverscript";
    /* SubsuperscriptBox[a, b, c] displays (and so stands for) Subscript[a, b]^c,
     * RadicalBox[a, n] the n-th root a^(1/n). */
    if (h == SYM_SubsuperscriptBox && n == 3) {
        buf_str(out, "Power[Subscript[(");
        if (!box_to_text(A[0], out)) return false;
        buf_str(out, "), (");
        if (!box_to_text(A[1], out)) return false;
        buf_str(out, ")], (");
        if (!box_to_text(A[2], out)) return false;
        buf_str(out, ")]");
        return true;
    }
    if (strcmp(h, "FractionBox") == 0 && n == 2) {   /* (a)/(b) */
        buf_str(out, "((");
        if (!box_to_text(A[0], out)) return false;
        buf_str(out, ")/(");
        if (!box_to_text(A[1], out)) return false;
        buf_str(out, "))");
        return true;
    }
    if (strcmp(h, "RadicalBox") == 0 && n == 2) {
        buf_str(out, "Power[(");
        if (!box_to_text(A[0], out)) return false;
        buf_str(out, "), 1/(");
        if (!box_to_text(A[1], out)) return false;
        buf_str(out, ")]");
        return true;
    }
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

/* ---- the infix box escapes ------------------------------------------------
 *
 * group_text tokenises the group content into box OPERANDS (already lowered to
 * InputForm text), box OPERATORS (the escapes), whitespace and the other
 * characters, then parses each run of operands and operators by the tiers in
 * the header comment and emits the expression they stand for as function-call
 * text, which carries its own grouping. Everything else is copied. */

enum { LT_UNIT, LT_OP, LT_WS, LT_CH };

typedef struct {
    int kind;
    char op;            /* LT_OP: '^' '_' '%' '&' '+' '/' '@' */
    char* text;         /* LT_UNIT / LT_WS / LT_CH: owned */
    char* inner;        /* LT_UNIT of a bare \( \) group: its contents, owned */
} LTok;

typedef struct { LTok* t; size_t n, cap; bool oom; } LToks;

static void ltoks_free(LToks* v) {
    for (size_t i = 0; i < v->n; i++) { free(v->t[i].text); free(v->t[i].inner); }
    free(v->t);
}

static void ltoks_push(LToks* v, int kind, char op, char* text, char* inner) {
    if (v->n == v->cap) {
        size_t cap = v->cap ? 2 * v->cap : 16;
        LTok* t = realloc(v->t, cap * sizeof(LTok));
        if (!t) { v->oom = true; free(text); free(inner); return; }
        v->t = t; v->cap = cap;
    }
    v->t[v->n].kind = kind; v->t[v->n].op = op;
    v->t[v->n].text = text; v->t[v->n].inner = inner;
    v->n++;
}

static char* dup_range(const char* p, size_t n) {
    char* d = malloc(n + 1);
    if (!d) return NULL;
    memcpy(d, p, n);
    d[n] = '\0';
    return d;
}

static bool lb_ident_char(unsigned char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
           c == '$' || c == '`' || c == '_' || c >= 0x80;
}

/* Past the bracket opening at p (one of ( [ {), its contents lowered between
 * the brackets into out; NULL if unbalanced. Strings and \( \) groups inside
 * are skipped over when matching. */
static const char* bracket_unit(const char* p, const char* end, Buf* out) {
    const char* q = p + 1;
    int depth = 1;
    while (q < end) {
        size_t l = 0;
        char c = lb_ctl(q, &l);
        if (c == '(') {
            const char* cend = NULL;
            q = group_close(q + l, &cend);
            if (!q || q > end) return NULL;
            continue;
        }
        if (c) { q += l; continue; }
        if (*q == '\\' && q + 1 < end) { q += 2; continue; }
        if (*q == '"') { q = skip_string(q); if (!q || q > end) return NULL; continue; }
        if (*q == '(' || *q == '[' || *q == '{') depth++;
        else if (*q == ')' || *q == ']' || *q == '}') {
            if (--depth == 0) break;
        }
        q++;
    }
    if (q >= end) return NULL;
    buf_add(out, p, 1);
    if (!group_text(p + 1, q, out)) return NULL;
    buf_add(out, q, 1);
    return q + 1;
}

/* One box operand at p, lowered into out: a primary and its postfix forms.
 * NULL when no operand starts at p (and nothing is written); *bad is set on a
 * malformed one. A bare \( \) group's contents go to *inner when it is the
 * whole operand. */
static const char* scan_unit(const char* p, const char* end, Buf* out, Buf* inner, bool* bad) {
    size_t l = 0;
    char c = lb_ctl(p, &l);
    const char* q = NULL;
    bool group = false;
    if (c == '(') {
        const char* cend = NULL;
        const char* after = group_close(p + l, &cend);
        if (!after || after > end) { *bad = true; return NULL; }
        Buf g = {0};
        if (!group_text(p + l, cend, &g)) { free(g.s); *bad = true; return NULL; }
        buf_str(out, "(");
        if (g.s) buf_str(out, g.s);
        buf_str(out, ")");
        if (inner && g.s) buf_str(inner, g.s);
        free(g.s);
        q = after;
        group = true;
    } else if (c == '*') {
        q = star_box(p + l, out);
        if (!q || q > end) { *bad = true; return NULL; }
    } else if (c) {
        return NULL;                              /* an operator escape */
    } else if (*p == '"') {
        q = skip_string(p);
        if (!q || q > end) { *bad = true; return NULL; }
        buf_add(out, p, (size_t)(q - p));
    } else if (*p == '(' || *p == '[' || *p == '{') {
        q = bracket_unit(p, end, out);
        if (!q) { *bad = true; return NULL; }
    } else if ((*p >= '0' && *p <= '9') || (*p == '.' && p + 1 < end && p[1] >= '0' && p[1] <= '9')) {
        q = p;                                    /* 12, 1.5, 2.5`20, 1.*^-3 */
        while (q < end && ((*q >= '0' && *q <= '9') || *q == '.')) q++;
        while (q < end && *q == '`') {
            q++;
            while (q < end && ((*q >= '0' && *q <= '9') || *q == '.' || *q == '`')) q++;
        }
        if (q + 1 < end && q[0] == '*' && q[1] == '^') {
            const char* r = q + 2;
            if (r < end && (*r == '-' || *r == '+')) r++;
            if (r < end && *r >= '0' && *r <= '9') {
                while (r < end && *r >= '0' && *r <= '9') r++;
                q = r;
            }
        }
        buf_add(out, p, (size_t)(q - p));
    } else if (*p == '#' || *p == '%') {          /* slots, Out */
        q = p + 1;
        while (q < end && (*q == *p || lb_ident_char((unsigned char)*q))) q++;
        buf_add(out, p, (size_t)(q - p));
    } else if (lb_ident_char((unsigned char)*p) || (*p == '\\' && p + 1 < end && p[1] == '[')) {
        q = p;                                    /* symbols, x_h, \[Alpha] */
        while (q < end) {
            if (*q == '\\' && q + 1 < end && q[1] == '[') {
                const char* r = q + 2;
                while (r < end && *r != ']') r++;
                if (r >= end) break;
                q = r + 1;
            } else if (lb_ident_char((unsigned char)*q) && !lb_ctl(q, &l)) {
                q++;
            } else break;
        }
        buf_add(out, p, (size_t)(q - p));
    } else {
        return NULL;
    }
    /* Postfix forms bind tighter than every box operator. */
    for (;;) {
        if (q < end && *q == '[') {
            q = bracket_unit(q, end, out);
            if (!q) { *bad = true; return NULL; }
        } else if (q < end && *q == '\'') {
            buf_add(out, q++, 1);
        } else if (q < end && *q == '!' && !(q + 1 < end && q[1] == '=')) {
            buf_add(out, q++, 1);
        } else if (q + 1 < end && *q == '@' && q[1] != '@' && q[1] != '*') {
            Buf rhs = {0};
            const char* r = scan_unit(q + 1, end, &rhs, NULL, bad);
            if (!r) { free(rhs.s); if (*bad) return NULL; break; }
            buf_str(out, "@");
            if (rhs.s) buf_str(out, rhs.s);
            free(rhs.s);
            q = r;
        } else break;
        group = false;
    }
    if (!group && inner) { free(inner->s); inner->s = NULL; inner->n = inner->cap = 0; }
    return q;
}

static bool lb_tokenize(const char* p, const char* end, LToks* v) {
    while (p < end && !v->oom) {
        size_t l = 0;
        char c = lb_ctl(p, &l);
        if (c == '!') { p += l; continue; }
        if (c == '^' || c == '_' || c == '%' || c == '&' || c == '+' || c == '/' || c == '@') {
            ltoks_push(v, LT_OP, c, NULL, NULL);
            p += l;
            continue;
        }
        if (c && c != '(' && c != '*') return false;   /* a stray \) or \` */
        /* "\ " is a space token in box text (p^\ 2 is RowBox[{p, ^, " ", 2}],
         * which reads as p^2): lower it to a plain space. */
        if (*p == '\\' && p + 1 < end && p[1] == ' ') {
            ltoks_push(v, LT_WS, 0, dup_range(" ", 1), NULL);
            p += 2;
            continue;
        }
        if (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') {
            const char* q = p;
            while (q < end && (*q == ' ' || *q == '\t' || *q == '\n' || *q == '\r')) q++;
            ltoks_push(v, LT_WS, 0, dup_range(p, (size_t)(q - p)), NULL);
            p = q;
            continue;
        }
        Buf u = {0}, in = {0};
        bool bad = false;
        const char* q = scan_unit(p, end, &u, &in, &bad);
        if (bad || u.oom || in.oom) { free(u.s); free(in.s); return false; }
        if (q) {
            ltoks_push(v, LT_UNIT, 0, u.s ? u.s : dup_range("", 0), in.s);
            p = q;
            continue;
        }
        free(u.s); free(in.s);
        size_t n = (*p == '\\' && p + 1 < end) ? 2 : 1;   /* other text */
        ltoks_push(v, LT_CH, 0, dup_range(p, n), NULL);
        p += n;
    }
    return !v->oom;
}

typedef struct { const LToks* v; size_t i; bool fail; } LParse;

static size_t lp_skip_ws(const LParse* P, size_t i) {
    while (i < P->v->n && P->v->t[i].kind == LT_WS) i++;
    return i;
}

/* Consume the operator `op` if it is next (across whitespace). */
static bool lp_take(LParse* P, char op) {
    size_t j = lp_skip_ws(P, P->i);
    if (j < P->v->n && P->v->t[j].kind == LT_OP && P->v->t[j].op == op) { P->i = j + 1; return true; }
    return false;
}

typedef char* (*LTier)(LParse*);

static char* lp_frac(LParse* P);
static char* lp_sqrt(LParse* P);

/* A right operand: optional signs, then the tier. */
static char* lp_signed(LParse* P, LTier tier) {
    size_t j = lp_skip_ws(P, P->i);
    if (j < P->v->n && P->v->t[j].kind == LT_CH &&
        (strcmp(P->v->t[j].text, "-") == 0 || strcmp(P->v->t[j].text, "+") == 0)) {
        P->i = j + 1;
        char* r = lp_signed(P, tier);
        if (!r) return NULL;
        Buf b = {0};
        buf_str(&b, P->v->t[j].text[0] == '-' ? "(-" : "(+");
        buf_str(&b, r);
        buf_str(&b, ")");
        free(r);
        return b.oom ? (free(b.s), (char*)NULL) : b.s;
    }
    return tier(P);
}

/* "head[a, b]" or "head[a, b, c]" (c may be NULL); frees the arguments. */
static char* lp_call(const char* head, char* a, char* b, char* c) {
    Buf o = {0};
    if (a && b) {
        buf_str(&o, head); buf_str(&o, "[");
        buf_str(&o, a); buf_str(&o, ", "); buf_str(&o, b);
        if (c) { buf_str(&o, ", "); buf_str(&o, c); }
        buf_str(&o, "]");
    }
    bool ok = a && b && !o.oom;
    free(a); free(b); free(c);
    if (!ok) { free(o.s); return NULL; }
    return o.s;
}

/* The innermost operand: a unit (its contents, for a bare group). */
static char* lp_unit(LParse* P) {
    size_t j = lp_skip_ws(P, P->i);
    if (j < P->v->n && P->v->t[j].kind == LT_OP && P->v->t[j].op == '@') return lp_sqrt(P);
    if (j >= P->v->n || P->v->t[j].kind != LT_UNIT) { P->fail = true; return NULL; }
    P->i = j + 1;
    const LTok* t = &P->v->t[j];
    const char* txt = t->inner ? t->inner : t->text;
    return dup_range(txt, strlen(txt));
}

/* \& \+ (and their \% second script), grouping to the right. */
static char* lp_over(LParse* P) {
    char* b = lp_unit(P);
    if (!b) return NULL;
    char op = lp_take(P, '&') ? '&' : lp_take(P, '+') ? '+' : 0;
    if (!op) return b;
    char* r = lp_signed(P, lp_over);
    if (!r) { free(b); return NULL; }
    if (lp_take(P, '%')) {
        char* r2 = lp_signed(P, lp_over);
        if (!r2) { free(b); free(r); return NULL; }
        /* UnderoverscriptBox[base, under, over] */
        return op == '+' ? lp_call("Underoverscript", b, r, r2)
                         : lp_call("Underoverscript", b, r2, r);
    }
    return lp_call(op == '&' ? "Overscript" : "Underscript", b, r, NULL);
}

/* Subscript[base, sub]^sup: the SubsuperscriptBox reading. */
static char* lp_subsup(char* base, char* sub, char* sup) {
    char* s = lp_call("Subscript", base, sub, NULL);
    if (!s) { free(sup); return NULL; }
    return lp_call("Power", s, sup, NULL);
}

/* \_ (and \_ \% sup), grouping to the right. */
static char* lp_sub(LParse* P) {
    char* b = lp_over(P);
    if (!b || !lp_take(P, '_')) return b;
    char* r = lp_signed(P, lp_sub);
    if (!r) { free(b); return NULL; }
    if (lp_take(P, '%')) {
        char* r2 = lp_signed(P, lp_sub);
        if (!r2) { free(b); free(r); return NULL; }
        return lp_subsup(b, r, r2);
    }
    return lp_call("Subscript", b, r, NULL);
}

/* \^ (and \^ \% sub), grouping to the right; its operand may be a \@ root. */
static char* lp_sup(LParse* P) {
    char* b = lp_sub(P);
    if (!b || !lp_take(P, '^')) return b;
    char* r = lp_signed(P, lp_sqrt);
    if (!r) { free(b); return NULL; }
    if (lp_take(P, '%')) {
        char* r2 = lp_signed(P, lp_sub);
        if (!r2) { free(b); free(r); return NULL; }
        return lp_subsup(b, r2, r);
    }
    return lp_call("Power", b, r, NULL);
}

/* \@a (and \@a \% n), prefix; otherwise a script expression. */
static char* lp_sqrt(LParse* P) {
    if (!lp_take(P, '@')) return lp_sup(P);
    char* a = lp_signed(P, lp_sqrt);
    if (!a) return NULL;
    if (lp_take(P, '%')) {
        char* n = lp_signed(P, lp_sqrt);
        if (!n) { free(a); return NULL; }
        Buf inv = {0};
        buf_str(&inv, "1/(");
        buf_str(&inv, n);
        buf_str(&inv, ")");
        free(n);
        if (inv.oom) { free(a); free(inv.s); return NULL; }
        return lp_call("Power", a, inv.s, NULL);
    }
    Buf o = {0};
    buf_str(&o, "Sqrt[");
    buf_str(&o, a);
    buf_str(&o, "]");
    free(a);
    if (o.oom) { free(o.s); return NULL; }
    return o.s;
}

/* \/ , grouping to the left. */
static char* lp_frac(LParse* P) {
    char* l = lp_sqrt(P);
    while (l && lp_take(P, '/')) {
        char* r = lp_signed(P, lp_sqrt);
        if (!r) { free(l); return NULL; }
        /* (a)/(b), as a FractionBox reads (not a held Divide[a, b]). */
        Buf o = {0};
        buf_str(&o, "(("); buf_str(&o, l); buf_str(&o, ")/("); buf_str(&o, r); buf_str(&o, "))");
        free(l); free(r);
        if (o.oom) { free(o.s); return NULL; }
        l = o.s;
    }
    return l;
}

/* The InputForm text of the group content [p, end). */
static bool group_text(const char* p, const char* end, Buf* out) {
    LToks v = {0};
    if (!lb_tokenize(p, end, &v)) { ltoks_free(&v); return false; }
    LParse P = { &v, 0, false };
    bool ok = true;
    while (ok && P.i < v.n) {
        const LTok* t = &v.t[P.i];
        bool box = false;
        if (t->kind == LT_OP) {
            if (t->op != '@') { ok = false; break; }   /* no left operand */
            box = true;
        } else if (t->kind == LT_UNIT) {
            size_t j = lp_skip_ws(&P, P.i + 1);
            box = j < v.n && v.t[j].kind == LT_OP && v.t[j].op != '@';
        }
        if (!box) { buf_str(out, t->text); P.i++; continue; }
        char* e = lp_frac(&P);
        if (!e || P.fail) { free(e); ok = false; break; }
        /* Parenthesised, so the call text cannot fuse with its neighbours
         * (2x\^2 is 2 (x^2), a.b\^2 is a . (b^2)). */
        buf_str(out, "(");
        buf_str(out, e);
        buf_str(out, ")");
        free(e);
    }
    ltoks_free(&v);
    return ok && !out->oom;
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
