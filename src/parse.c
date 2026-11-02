#include "parse.h"
#include "utf8.h"   /* \[Name] / \:hhhh escapes in string literals */
#include "linear_boxes.h"   /* \!\( ... \) */
#include "common.h"
#include "expr.h"
#include "context.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <gmp.h>
#include "numeric.h"

// Private parser state
typedef struct {
    const char* input;  // Original input string
    const char* pos;    // Current position
    int bracket_depth;  // Nesting depth of (), [], {}, <|...|>, [[...]].
                        // Newlines are insignificant inside brackets, so a
                        // multi-line list/call stays one expression.
    int saw_newline;    // Set by skip_whitespace: did the most recent skip
                        // cross a genuine (non-comment) newline? Used to make a
                        // top-level line break terminate a statement instead of
                        // acting as implicit multiplication.
    MthSpanSink* spans; // When non-NULL (mth_parse_spans only), every
                        // subexpression's source range is recorded here for the
                        // notebook's structural selection. NULL on the normal
                        // parse/eval path -- all existing {..} initializers omit
                        // it, so C99 zero-fills it and that path is unchanged.
    int colon_hold;     // bracket_depth+1 while parsing the right operand of
                        // `sym:p` (0 = none): at that bracket level a further
                        // `:` ends the pattern and becomes its Optional default
                        // (`a:b:c` is Optional[a:b, c]), so parse_colon_suffix
                        // must not claim it for an inner token.
    int line_statements; // Set only by parse_statement_sequence (ToExpression
                        // in InputForm): a top-level line break ends the
                        // statement before ANY operator, not just before an
                        // implicit Times, and a `;` at the end of a line takes
                        // an implicit Null instead of continuing onto the next
                        // line -- so `a\n// f` is a syntax error and `a;\nb` is
                        // two statements, as in the reference language.
    const char* ws_start; // [ws_start, ws_end) is the most recent run of
    const char* ws_end;   // whitespace/comments skip_whitespace crossed. After
                        // a top-level parse ws_start is the tight end of the
                        // statement, before the whitespace and line breaks the
                        // infix loop skipped while looking for an operator.
    Expr* seed_left;    // When non-NULL, parse_expression_prec takes it (owned)
                        // as its already-parsed left operand and goes straight
                        // to the infix loop -- how parse_next_expression
                        // continues `lhs := rhs;/; cond` past a statement `;`.
    const char* chain_head; // Set by parse_expression_prec on return: the
                        // flat operator chain (Plus/Times/...) it built its
                        // result with, or NULL -- so a caller can tell a bare
                        // `b c` operand from a parenthesised `(b c)`.
} ParserState;

/* Quiet mode: syntax errors are not printed, only signalled by the NULL result,
 * so a caller that reports through the message funnel (Get, ToExpression,
 * SyntaxQ) is not doubled by -- or, under Quiet[], leaked through -- the
 * parser's own diagnostics.  Out-of-memory notes are always printed. */
static int g_parse_quiet = 0;

int parse_set_quiet(int quiet) {
    int prev = g_parse_quiet;
    g_parse_quiet = quiet;
    return prev;
}

static void parse_err(const char* fmt, ...) {
    if (g_parse_quiet) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

/* ------------------------------------------------------------------ */
/* Non-ASCII characters in source code                                */
/* ------------------------------------------------------------------ */

static long char_escape_cp(const char* p, size_t* esc_len);

/* The non-ASCII character at p -- a raw UTF-8 sequence or a character escape
 * (\[Name], \:hhhh, \|hhhhhh, ...) -- as a code point.  Returns its length in
 * the source, or 0 when p starts an ASCII character (an escape that names an
 * ASCII character, \[RawPlus], counts as ASCII: it is not lexed here). */
static size_t src_nonascii(const char* p, uint32_t* cp) {
    unsigned char c = (unsigned char)p[0];
    if (c >= 0x80) {
        size_t n = 0, i = 0;
        while (n < 4 && p[n]) n++;
        *cp = u8_decode(p, n, &i);
        return i;
    }
    if (c == '\\') {
        size_t len = 0;
        long v = char_escape_cp(p, &len);
        if (v >= 0x80) { *cp = (uint32_t)v; return len; }
    }
    return 0;
}

/* Length of the identifier letter at p that is not ASCII (\[Mu], a raw "μ"),
 * or 0.  Every non-ASCII character is a letter unless the class table
 * (u8_char_class) says it is an operator, whitespace or not source at all. */
static size_t src_letter(const char* p) {
    uint32_t cp;
    size_t n = src_nonascii(p, &cp);
    if (!n) return 0;
    if (lb_cp_escape(cp)) return 0;              /* a box control character */
    const U8CharInfo* ci = u8_char_class(cp);
    return (!ci || ci->cls == U8C_LETTER) ? n : 0;
}

/* The class entry of the non-letter, non-ASCII character at p (an operator,
 * whitespace, ...), with its source length in *len; NULL otherwise. */
static const U8CharInfo* src_special(const char* p, size_t* len) {
    uint32_t cp;
    size_t n = src_nonascii(p, &cp);
    if (!n) return NULL;
    const U8CharInfo* ci = u8_char_class(cp);
    if (!ci || ci->cls == U8C_LETTER) return NULL;
    *len = n;
    return ci;
}

/* Record the half-open byte range [start, end) of one subexpression into the
 * sink. No-op when no sink is attached (the normal path) or the range is empty.
 * `start`/`end` are pointers into s->input, so offsets are absolute even across
 * ;-separated statements. */
static void record_span(ParserState* s, const char* start, const char* end) {
    if (!s->spans) return;
    /* Trim trailing whitespace: a subexpression's end can sit past a space when
     * its last sub-parse skipped trailing whitespace before finding no further
     * operator (e.g. the `x` in `2 x + 1` leaves pos at the `+`, so the Times
     * span would otherwise read "2 x "). */
    while (end > start && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n' ||
                           end[-1] == '\r' || end[-1] == '\f' || end[-1] == '\v'))
        end--;
    if (end <= start) return;
    MthSpanSink* k = s->spans;
    if (k->n == k->cap) {
        size_t nc = k->cap ? k->cap * 2 : 32;
        MthSpan* nv = realloc(k->v, nc * sizeof(MthSpan));
        if (!nv) return; /* drop this span rather than crash on OOM */
        k->v = nv;
        k->cap = nc;
    }
    k->v[k->n].start = (int)(start - s->input);
    k->v[k->n].end = (int)(end - s->input);
    k->n++;
}

// Forward declarations
static Expr* parse_expression_state(ParserState* s);
static void skip_whitespace(ParserState* s);

/* If `left` is an n-ary call whose head is the symbol `head_name`, append
 * `right` as a new argument in-place (taking ownership of `right`) and return
 * true. Otherwise return false. Used during parsing so that repeated
 * left-associative uses of Flat operators (Plus, Times) produce an n-ary
 * call (Plus[a, b, c]) instead of a nested binary tree (Plus[Plus[a, b], c]).
 * This matters for held expressions: without this, Length[Hold[a+b+c]] would
 * misleadingly report 2. */
static bool extend_flat_head(Expr* left, const char* head_name, Expr* right) {
    if (!head_is(left, intern_symbol(head_name))) return false;

    size_t old_count = left->data.function.arg_count;
    Expr** new_args = realloc(left->data.function.args, sizeof(Expr*) * (old_count + 1));
    if (!new_args) return false;
    new_args[old_count] = right;
    left->data.function.args = new_args;
    left->data.function.arg_count = old_count + 1;
    return true;
}

/* Heads that participate in Mathematica-style chained comparisons. WL
 * collapses `a < b <= c == d > e` into one Inequality[...] node. Unequal
 * (`!=`) is handled separately: a run of `!=` alone is the variadic
 * Unequal[a, b, c] ("all distinct"), and it joins an Inequality only when
 * mixed with these heads (see fold_into_inequality). */
static bool is_chain_compare_head(const char* head_name) {
    return strcmp(head_name, "Less") == 0
        || strcmp(head_name, "LessEqual") == 0
        || strcmp(head_name, "Greater") == 0
        || strcmp(head_name, "GreaterEqual") == 0
        || strcmp(head_name, "Equal") == 0;
}

/* Extend an existing Inequality[v0, op0, v1, op1, ..., vk] in place by
 * appending the operator symbol new_op_name and the new value new_val.
 * Mirrors the structure of extend_flat_head but for Inequality. */
static bool extend_inequality(Expr* ineq, const char* new_op_name, Expr* new_val) {
    if (!head_is(ineq, intern_symbol("Inequality"))) return false;
    size_t old_count = ineq->data.function.arg_count;
    Expr** new_args = realloc(ineq->data.function.args,
                              sizeof(Expr*) * (old_count + 2));
    if (!new_args) return false;
    new_args[old_count]     = expr_new_symbol(new_op_name);
    new_args[old_count + 1] = new_val;
    ineq->data.function.args = new_args;
    ineq->data.function.arg_count = old_count + 2;
    return true;
}

/* Comparison heads that are variadic predicates rather than pairwise chains:
 * a run of the SAME operator collects into one call, SameQ[a, b, c]. */
static bool is_nary_compare_head(const char* head_name) {
    return strcmp(head_name, "SameQ") == 0
        || strcmp(head_name, "UnsameQ") == 0
        || strcmp(head_name, "Unequal") == 0;
}

/* Fold `left op right` into an Inequality, taking ownership of `right` and
 * possibly replacing *left. `left` is either an Inequality (extended in place)
 * or a bare comparison H[x1, ..., xk] with H a chain head or Unequal, which
 * becomes Inequality[x1, H, x2, ..., H, xk, op, right]. Returns false, leaving
 * both untouched, when `left` is neither. */
static bool fold_into_inequality(Expr** left, const char* op_name, Expr* right) {
    Expr* l = *left;
    if (extend_inequality(l, op_name, right)) return true;
    if (l->type != EXPR_FUNCTION || l->data.function.head->type != EXPR_SYMBOL) return false;
    const char* inner = l->data.function.head->data.symbol.name;
    if (!is_chain_compare_head(inner) && inner != SYM_Unequal) return false;
    size_t k = l->data.function.arg_count;
    if (k < 2) return false;
    size_t n = 2 * k + 1;
    Expr** args = malloc(sizeof(Expr*) * n);
    if (!args) return false;
    size_t j = 0;
    for (size_t i = 0; i < k; i++) {
        if (i > 0) args[j++] = expr_new_symbol(inner);
        args[j++] = l->data.function.args[i];
        l->data.function.args[i] = NULL;    /* stolen */
    }
    args[j++] = expr_new_symbol(op_name);
    args[j++] = right;
    expr_free(l);
    *left = expr_new_function(expr_new_symbol(SYM_Inequality), args, n);
    free(args);
    return true;
}

/* ------------------- Basic Token Parsers ------------------- */

/* Leave one bracket level. A line break seen inside the brackets is not a
 * statement boundary, so the sticky saw_newline flag is cleared with it:
 * otherwise `f[x\n] b` would look as though a newline separated `f[x]` from
 * `b` and split into two statements. */
static void close_bracket(ParserState* s) {
    s->bracket_depth--;
    s->saw_newline = 0;
}

// Skips whitespace and comments
static void skip_whitespace(ParserState* s) {
    /* saw_newline is sticky: set here when a genuine line break is skipped, and
     * cleared only when a primary token is actually consumed (see
     * parse_primary). This survives the nested skip_whitespace calls made while
     * parsing a sub-expression's right operand, so the infix loop can still see
     * that a newline separated two juxtaposed statements. */
    /* A skip that starts where the previous one ended is still inside that
     * run (nested parse levels each skip after the same token), so it keeps
     * the run's start. */
    if (s->pos != s->ws_end) s->ws_start = s->pos;
    while (1) {
        if (s->pos[0] == '\\' && (s->pos[1] == '\n' || s->pos[1] == '\r')) {
            /* Line continuation: a backslash immediately before a newline
             * joins the two physical lines into one logical line
             * (Mathematica syntax). Consume the backslash and the newline
             * WITHOUT setting saw_newline — the continued line is a
             * continuation of the current statement, not a new one. Handle
             * LF, CR, and CRLF endings. */
            s->pos++;                              /* skip '\\' */
            if (s->pos[0] == '\r') s->pos++;       /* CR (or CR of CRLF) */
            if (s->pos[0] == '\n') s->pos++;       /* LF (or LF of CRLF) */
        } else if (isspace(*s->pos)) {
            if (*s->pos == '\n' || *s->pos == '\r') s->saw_newline = 1;
            s->pos++;
        } else if ((unsigned char)s->pos[0] >= 0x80 || s->pos[0] == '\\') {
            /* \[InvisibleSpace], \[VeryThinSpace], ... are whitespace and
             * \[IndentingNewLine] is a line break, raw or escaped. */
            size_t n = 0;
            const U8CharInfo* ci = src_special(s->pos, &n);
            if (!ci || (ci->cls != U8C_SPACE && ci->cls != U8C_NEWLINE)) break;
            if (ci->cls == U8C_NEWLINE) s->saw_newline = 1;
            s->pos += n;
        } else if (s->pos[0] == '(' && s->pos[1] == '*') {
            int depth = 1;
            s->pos += 2;
            while (*s->pos != '\0' && depth > 0) {
                if (s->pos[0] == '(' && s->pos[1] == '*') {
                    depth++;
                    s->pos += 2;
                } else if (s->pos[0] == '*' && s->pos[1] == ')') {
                    depth--;
                    s->pos += 2;
                } else {
                    s->pos++;
                }
            }
        } else {
            break;
        }
    }
    s->ws_end = s->pos;
}

/* Can the text at p only continue an expression -- does it start a binary or
 * postfix operator that no expression can begin with (`a; ::x` is a MessageName
 * of the compound, too)?  `+`, `-`, `!`, `.` (a number), `?`, `'`, `;;`, `<|`
 * and `<<` are excluded: each may begin a new statement.  Shared by the
 * top-level statement reader and the infix `;` operator. */
static bool starts_binary_operator(const char* p) {
    switch (p[0]) {
        case '*': case '/': case '^': case '=': case '>': case '|':
        case '&': case '@': case '~':
            return true;
        case ':':                               /* := :> :: (MessageName) */
            return true;
        case '<':                               /* <| association, << Get */
            return p[1] != '|' && p[1] != '<';
        case '!':
            return p[1] == '=';                 /* != ; !x is a new statement */
        case '-':
            return p[1] == '>' || p[1] == '=';  /* -> and -= ; -x is a new statement */
        case '+':
            return p[1] == '=';                 /* += ; +x is a new statement */
        default:
            return false;
    }
}

/* Skip the whitespace after an operator whose right operand may be omitted
 * (`a;`, `a;;`). In line-statement mode a top-level line break there ends the
 * statement: rewind to just after the operator, so the infix loop still sees
 * the break, and return 1 -- the operand is omitted. Otherwise return 0 with
 * the whitespace skipped. */
static int operand_omitted_at_line_end(ParserState* s) {
    if (!s->line_statements || s->bracket_depth != 0) {
        skip_whitespace(s);
        return 0;
    }
    const char* after = s->pos;
    s->saw_newline = 0;
    skip_whitespace(s);
    if (!s->saw_newline) return 0;
    s->pos = s->ws_start = s->ws_end = after;
    s->saw_newline = 0;
    return 1;
}

// Parses a symbol (x, `name`, $var) and resolves it through the context
// system, producing a canonical (possibly qualified) symbol name.
/* The buffer grows for the same reason parse_string's does: dropping the tail
 * of an over-long name silently ALIASES two distinct symbols onto one, which is
 * a wrong answer rather than a slow one.  Context-qualified names
 * (Long`Nested`Context`name) are the realistic way to exceed 255. */
static Expr* parse_symbol(ParserState* s) {
    char stackbuf[256];
    char* buffer = stackbuf;
    size_t cap = sizeof(stackbuf);
    size_t i = 0;

    /* A non-ASCII letter (\[Mu], a raw "μ") is stored as its UTF-8 encoding,
     * so every spelling of a character names the same symbol. */
    size_t letter_len = 0;
    while (isalnum((unsigned char)*s->pos) || *s->pos == '`' || *s->pos == '$' ||
           (letter_len = src_letter(s->pos)) > 0) {
        if (i + 5 >= cap) {
            size_t ncap = cap * 2;
            char* nb = (buffer == stackbuf) ? (char*)malloc(ncap)
                                            : (char*)realloc(buffer, ncap);
            if (!nb) {
                if (buffer != stackbuf) free(buffer);
                fprintf(stderr, "Out of memory parsing symbol name\n");
                return NULL;
            }
            if (buffer == stackbuf) memcpy(nb, stackbuf, i);
            buffer = nb;
            cap = ncap;
        }
        if (letter_len) {
            uint32_t cp = 0;
            src_nonascii(s->pos, &cp);
            i += u8_encode(cp, buffer + i);
            s->pos += letter_len;
            letter_len = 0;
            continue;
        }
        buffer[i++] = *s->pos;
        s->pos++;
    }
    buffer[i] = '\0';

    /* A one-character name that is a constant's letter is that constant:
     * \[Pi] is Pi and \[ImaginaryI] is I, while a\[Pi] stays one symbol. */
    {
        uint32_t cp = 0;
        size_t k = 0;
        if ((unsigned char)buffer[0] >= 0x80 && (cp = u8_decode(buffer, i, &k), k == i)) {
            const U8CharInfo* ci = u8_char_class(cp);
            if (ci && ci->cls == U8C_LETTER && ci->head) {
                Expr* out = expr_new_symbol(ci->head);
                if (buffer != stackbuf) free(buffer);
                return out;
            }
        }
    }

    char* resolved = context_resolve_name(buffer);
    Expr* out = expr_new_symbol(resolved ? resolved : buffer);
    free(resolved);
    if (buffer != stackbuf) free(buffer);
    return out;
}

/* Parse an optional Mathematica-style precision / accuracy suffix
 * appended to a numeric literal. Returns the kind of suffix found:
 *    0 = no suffix             (leave `value` as-is)
 *    1 = single backtick       (`n or bare `)
 *    2 = double backtick       (``n)
 * and, when a numeric precision follows, stores its value in *out_digits.
 * *has_digits is set to true if a digit value was consumed.
 *
 * Disambiguation from the context-separator use of backtick:
 * a precision suffix is only recognized when the character after the
 * backtick(s) is a digit, `.`, `+`, or `-`, or (for a single bare-`
 * with no digit) a non-identifier character. A backtick immediately
 * followed by a letter / underscore / `$` is treated as a context
 * separator and left for the identifier parser to handle — we do NOT
 * consume it here. */
static int parse_precision_suffix(ParserState* s, double* out_digits, bool* has_digits) {
    *has_digits = false;
    if (*s->pos != '`') return 0;

    /* Peek: is this a precision suffix or a context separator? */
    const char* look = s->pos + 1;
    int bt_count = 1;
    if (*look == '`') { bt_count = 2; look++; }

    /* A context-separator backtick would be followed by an identifier
     * character. Precision markers are followed by a digit/./+/-, or (for
     * the bare single-backtick form) end-of-number. */
    bool looks_like_number = (isdigit((unsigned char)*look) || *look == '.'
                              || *look == '+' || *look == '-');
    bool looks_like_context = (isalpha((unsigned char)*look) || *look == '_'
                               || *look == '$');

    if (looks_like_context) return 0;  /* leave for identifier parser */

    if (bt_count == 2 && !looks_like_number) {
        /* Double backtick MUST be followed by an accuracy value. */
        return 0;
    }

    /* Consume the backtick(s). */
    s->pos += bt_count;

    if (looks_like_number) {
        char* end;
        double v = strtod(s->pos, &end);
        /* strtod would swallow the first dot of a following `..`/`...`
         * (`2.5`10..` is Repeated[2.5`10]); give it back, as parse_number
         * does for the mantissa. */
        if (end > s->pos && end[-1] == '.' && *end == '.') end--;
        if (end != s->pos) {
            s->pos = end;
            *out_digits = v;
            *has_digits = true;
        }
    }
    return bt_count;
}

/* Build an MPFR Expr from a mantissa string `mantissa_str` at the
 * requested `digits` of Mathematica-style precision. `accuracy_mode`
 * means `digits` is accuracy instead of precision — bits must be
 * derived from both the requested accuracy and the value's magnitude.
 *
 * In a USE_MPFR=0 build this function emits a one-shot warning and
 * returns an EXPR_REAL constructed via strtod. */
static Expr* build_precision_literal(const char* mantissa_str, double digits,
                                     bool accuracy_mode) {
#ifdef USE_MPFR
    /* Convert digits → bits. For precision mode, bits = ceil(d * log2(10)).
     * For accuracy mode the precision needed is approximately
     *     bits ≈ ceil((d + log10|value|) * log2(10))
     * We first parse the value at a generous precision to measure its
     * magnitude, then re-round at the final precision. */
    long bits_prec;
    if (accuracy_mode) {
        /* Parse once to find magnitude. */
        mpfr_t tmp; mpfr_init2(tmp, 128);
        mpfr_set_str(tmp, mantissa_str, 10, MPFR_RNDN);
        double log10_abs = 0.0;
        if (!mpfr_zero_p(tmp)) {
            mpfr_t abs_t; mpfr_init2(abs_t, 128);
            mpfr_abs(abs_t, tmp, MPFR_RNDN);
            log10_abs = floor(log10(mpfr_get_d(abs_t, MPFR_RNDN)));
            mpfr_clear(abs_t);
        }
        double total_digits = digits + log10_abs + 1.0;
        if (total_digits < 1.0) total_digits = 1.0;
        bits_prec = (long)ceil(total_digits * 3.3219280948873626);
        mpfr_clear(tmp);
    } else {
        if (digits < 1.0) digits = 1.0;
        bits_prec = (long)ceil(digits * 3.3219280948873626);
    }
    if (bits_prec < 2) bits_prec = 2;
    return expr_new_mpfr_from_str(mantissa_str, (mpfr_prec_t)bits_prec);
#else
    (void)accuracy_mode; (void)digits;
    static bool warned = false;
    if (!warned) {
        parse_err(
                "parser: precision literal ignored (USE_MPFR=0); using "
                "machine precision.\n");
        warned = true;
    }
    char* end;
    double dval = strtod(mantissa_str, &end);
    return expr_new_real(dval);
#endif
}

/* Scale `val` by 10^exp, consuming `val` and returning a new Expr. Used to
 * implement Mathematica's `*^` scaled-scientific notation (e.g. `1.23*^4`
 * → 12300., `123*^-4` → 123/10000). Preserves the mantissa's type when
 * possible: Real stays Real, Integer with non-negative exp stays Integer
 * (possibly promoting to BigInt), Integer with negative exp becomes a
 * reduced Rational. MPFR mantissas multiply/divide by 10^|exp| at the
 * mantissa's own precision. */
static Expr* apply_pow10_scale(Expr* val, long exp) {
    if (val == NULL) return NULL;

    if (val->type == EXPR_REAL) {
        double d = val->data.real * pow(10.0, (double)exp);
        expr_free(val);
        return expr_new_real(d);
    }

#ifdef USE_MPFR
    if (val->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(val->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        unsigned long abs_exp = (unsigned long)(exp >= 0 ? exp : -exp);
        mpfr_t pow10;
        mpfr_init2(pow10, prec);
        mpfr_ui_pow_ui(pow10, 10, abs_exp, MPFR_RNDN);
        if (exp >= 0) {
            mpfr_mul(out->data.mpfr, val->data.mpfr, pow10, MPFR_RNDN);
        } else {
            mpfr_div(out->data.mpfr, val->data.mpfr, pow10, MPFR_RNDN);
        }
        mpfr_clear(pow10);
        expr_free(val);
        return out;
    }
#endif

    if (val->type == EXPR_INTEGER || val->type == EXPR_BIGINT) {
        mpz_t base;
        mpz_init(base);
        expr_to_mpz(val, base);
        expr_free(val);

        mpz_t pow10;
        mpz_init(pow10);
        mpz_ui_pow_ui(pow10, 10, (unsigned long)(exp >= 0 ? exp : -exp));

        if (exp >= 0) {
            mpz_mul(base, base, pow10);
            Expr* r = expr_new_bigint_from_mpz(base);
            mpz_clear(base);
            mpz_clear(pow10);
            return expr_bigint_normalize(r);
        }

        /* Negative exp: build a reduced Rational[num, den]. */
        mpz_t g; mpz_init(g);
        mpz_gcd(g, base, pow10);
        if (mpz_cmp_ui(g, 1) > 0) {
            mpz_divexact(base, base, g);
            mpz_divexact(pow10, pow10, g);
        }
        mpz_clear(g);

        if (mpz_cmp_ui(pow10, 1) == 0) {
            Expr* r = expr_new_bigint_from_mpz(base);
            mpz_clear(base);
            mpz_clear(pow10);
            return expr_bigint_normalize(r);
        }

        Expr* num = expr_bigint_normalize(expr_new_bigint_from_mpz(base));
        Expr* den = expr_bigint_normalize(expr_new_bigint_from_mpz(pow10));
        mpz_clear(base);
        mpz_clear(pow10);

        Expr* rargs[2] = { num, den };
        return expr_new_function(expr_new_symbol(SYM_Rational), rargs, 2);
    }

    /* Fallback for any other numeric atom: leave it to the evaluator. */
    Expr* p_args[2] = { expr_new_integer(10), expr_new_integer((int64_t)exp) };
    Expr* p = expr_new_function(expr_new_symbol(SYM_Power), p_args, 2);
    Expr* t_args[2] = { val, p };
    return expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
}

// Parses numbers (integers, reals, scientific notation), plus optional
// Mathematica-style precision/accuracy suffix (`n, ``n) and optional
// `*^exponent` scaled-scientific-notation suffix.
static Expr* parse_number(ParserState* s) {
    char* end;

    // Check if it's potentially a real number
    int is_real = 0;
    const char* p = s->pos;
    if (*p == '-' || *p == '+') p++;
    while (isdigit(*p)) p++;
    /* `..` / `...` are single tokens (Repeated / RepeatedNull) and win by
     * longest munch, as in the Wolfram Language: `0..` is Repeated[0] and
     * `0...` is RepeatedNull[0], so a '.' that is itself followed by '.' is
     * never the literal's decimal point (issue #100). */
    /* An `e`/`E` exponent (a Mathilda extension: the reference has only
     * `*^`) counts only when digits follow it: `2e5` is 200000., while
     * `2ep` is 2 ep and `3E` is 3 E, as the reference reads them -- the 'e'
     * starts a symbol. It used to make `2ep` the Real 2. times ep. */
    if ((*p == '.' && p[1] != '.') ||
        ((*p == 'e' || *p == 'E') &&
         (isdigit((unsigned char)p[1]) ||
          ((p[1] == '+' || p[1] == '-') && isdigit((unsigned char)p[2]))))) {
        is_real = 1;
    }

    /* Snapshot the mantissa text — used if a precision suffix follows,
     * so we can feed an exact string to MPFR rather than losing bits via
     * a double round-trip. */
    const char* mantissa_start = s->pos;
    Expr* result;

    if (is_real) {
        /* Scan the literal once, identifying:
         *   - the *mantissa* span (sign + integer digits + '.' + fraction);
         *   - the *full* span (mantissa + optional 'e[+-]?digits' exponent);
         *   - frac_digits = digits after the decimal point in the mantissa.
         *
         * The mantissa portion (without the exponent) drives the implied-
         * precision decision: a literal like `1.0e22` has mantissa "1.0"
         * with one fractional digit and magnitude 1, so implied precision
         * is just 1 — well below MachinePrecision — and the result stays
         * a machine double. A literal like `1.234567890123456789012345`
         * has implied precision ≈ 24 and is promoted to MPFR. */
        const char* lp = mantissa_start;
        if (*lp == '+' || *lp == '-') lp++;
        int int_digits = 0;
        while (isdigit((unsigned char)*lp)) { lp++; int_digits++; }
        bool has_dot = false;
        int frac_digits = 0;
        if (*lp == '.' && lp[1] != '.') {   /* `1..` / `1...`: see above */
            has_dot = true;
            lp++;
            while (isdigit((unsigned char)*lp)) { lp++; frac_digits++; }
        }
        const char* mantissa_end = lp;
        /* Optional exponent — only consume it if at least one digit follows.
         * `1.5e` (trailing bare 'e') must not eat the 'e'. */
        if (*lp == 'e' || *lp == 'E') {
            const char* save = lp;
            lp++;
            if (*lp == '+' || *lp == '-') lp++;
            if (isdigit((unsigned char)*lp)) {
                while (isdigit((unsigned char)*lp)) lp++;
            } else {
                lp = save;
            }
        }
        const char* full_end = lp;

        /* No digits at all (e.g. bare ".") → not a valid number. */
        if (int_digits == 0 && frac_digits == 0) return NULL;

        s->pos = full_end;

        /* Build a NUL-terminated copy of the full literal text. */
        size_t flen = (size_t)(full_end - mantissa_start);
        char fbuf[128];
        char* fheap = NULL;
        char* fstr = fbuf;
        if (flen + 1 > sizeof(fbuf)) {
            fheap = malloc(flen + 1);
            fstr = fheap;
        }
        memcpy(fstr, mantissa_start, flen);
        fstr[flen] = '\0';

        /* Compute implied precision in decimal digits from the mantissa:
         *     implied = frac_digits + log10(|mantissa_value|)
         * If implied ≤ MachinePrecision, keep a machine double; otherwise
         * build an MPFR at ceil(implied * log2(10)) bits, fed directly
         * from the literal text so we don't lose precision via a double
         * round-trip. Pure dot-less reals (no fractional part typed) and
         * zero mantissas always stay machine precision. */
        long mpfr_bits = 0;
        if (has_dot) {
            size_t mlen = (size_t)(mantissa_end - mantissa_start);
            char mbuf[128];
            char* mheap = NULL;
            char* mstr = mbuf;
            if (mlen + 1 > sizeof(mbuf)) {
                mheap = malloc(mlen + 1);
                mstr = mheap;
            }
            memcpy(mstr, mantissa_start, mlen);
            mstr[mlen] = '\0';
            double mval = strtod(mstr, NULL);
            /* significant digits of the literal text (leading zeros do not
             * count): strtod may round 9.9999999999999999 up to 10. */
            int sig = 0;
            bool lead = true;
            for (const char* c = mstr; *c; c++) {
                if (!isdigit((unsigned char)*c)) continue;
                if (lead && *c == '0') continue;
                lead = false;
                sig++;
            }
            if (mheap) free(mheap);
            if (isfinite(mval) && mval != 0.0 && sig >= 18) {
                double implied = (double)frac_digits + log10(fabs(mval));
                /* A literal of up to 17 significant digits -- enough to
                 * name every double -- is a machine real, as in the
                 * reference (Precision[1.2887955878011237] is
                 * MachinePrecision); 18 or more make it arbitrary precision
                 * (Precision[1.00000000000000001] is 17.), at the implied
                 * precision. */
                if (implied > 0.0) {
                    mpfr_bits = (long)ceil(implied * 3.3219280948873626);
                    if (mpfr_bits < 2) mpfr_bits = 2;
                }
            }
        }

#ifdef USE_MPFR
        if (mpfr_bits > 0) {
            result = expr_new_mpfr_from_str(fstr, (mpfr_prec_t)mpfr_bits);
            if (!result) {
                /* mpfr_set_str rejected the text — should not happen for a
                 * well-formed decimal literal. Fall back to a machine
                 * double so the parse still succeeds. */
                result = expr_new_real(strtod(fstr, NULL));
            }
        } else {
            result = expr_new_real(strtod(fstr, NULL));
        }
#else
        (void)mpfr_bits;
        result = expr_new_real(strtod(fstr, NULL));
#endif
        if (fheap) free(fheap);
        (void)end;
    } else {
        const char* start = s->pos;
        int negative = 0;
        if (*start == '-') {
            negative = 1;
            start++;
        } else if (*start == '+') {
            start++;
        }

        const char* current = start;
        while (isdigit(*current)) {
            current++;
        }

        if (current == start) return NULL;

        // Collect digit string into buffer (with optional sign)
        size_t digit_len = (size_t)(current - start);
        size_t buf_len = digit_len + (negative ? 1 : 0) + 1;
        char buf[256];
        char* heap_buf = NULL;
        char* num_str = buf;
        if (buf_len > sizeof(buf)) {
            heap_buf = malloc(buf_len);
            num_str = heap_buf;
        }
        size_t idx = 0;
        if (negative) num_str[idx++] = '-';
        memcpy(num_str + idx, start, digit_len);
        num_str[idx + digit_len] = '\0';

        s->pos = current;

        if (digit_len > 19) {
            // More digits than INT64_MAX (19 digits) — skip strtoll, go straight to bigint
            result = expr_new_bigint_from_str(num_str);
        } else {
            errno = 0;
            char* endptr;
            long long llval = strtoll(num_str, &endptr, 10);
            result = (errno == ERANGE) ? expr_new_bigint_from_str(num_str) : expr_new_integer((int64_t)llval);
        }
        if (heap_buf) free(heap_buf);
    }

    /* Optional precision / accuracy suffix. `3.14`50` → 50-digit MPFR,
     * `3.14``49` → 49-digit accuracy, `3`30` → integer 3 at 30 digits. */
    {
        double suffix_digits = 0.0;
        bool has_digits = false;
        int kind = parse_precision_suffix(s, &suffix_digits, &has_digits);
        if (kind != 0) {
            const char* mantissa_end = mantissa_start;
            /* Walk forward to end of the mantissa text we already consumed.
             * We recognize: optional sign, digits, optional '.', digits,
             * optional exponent. This mirrors the pre-suffix scan above. */
            if (*mantissa_end == '+' || *mantissa_end == '-') mantissa_end++;
            while (isdigit((unsigned char)*mantissa_end)) mantissa_end++;
            if (*mantissa_end == '.' && mantissa_end[1] != '.') {
                mantissa_end++;
                while (isdigit((unsigned char)*mantissa_end)) mantissa_end++;
            }
            if (*mantissa_end == 'e' || *mantissa_end == 'E') {
                mantissa_end++;
                if (*mantissa_end == '+' || *mantissa_end == '-') mantissa_end++;
                while (isdigit((unsigned char)*mantissa_end)) mantissa_end++;
            }
            size_t mlen = (size_t)(mantissa_end - mantissa_start);
            char mbuf[128];
            char* mheap = NULL;
            char* mstr = mbuf;
            if (mlen + 1 > sizeof(mbuf)) {
                mheap = malloc(mlen + 1);
                mstr = mheap;
            }
            memcpy(mstr, mantissa_start, mlen);
            mstr[mlen] = '\0';

            /* Default precision when no number follows a single backtick:
             * treat as MachinePrecision (leave the existing real/integer
             * untouched — this is the `3.14` form). */
            Expr* replacement = NULL;
            if (kind == 1 && !has_digits) {
                /* Bare `: force machine precision. If we parsed an integer,
                 * promote to EXPR_REAL so the printer shows a dot. */
                if (result && result->type != EXPR_REAL) {
                    double dval = strtod(mstr, &end);
                    replacement = expr_new_real(dval);
                }
            } else if (has_digits) {
                replacement = build_precision_literal(mstr, suffix_digits,
                                                     /*accuracy_mode=*/(kind == 2));
            }
            if (mheap) free(mheap);
            if (replacement) {
                expr_free(result);
                result = replacement;
            }
        }
    }

    /* Optional Mathematica scaled-scientific notation suffix:
     *   mantissa *^ signed-integer-exponent
     * e.g. `1.23*^4` → 12300., `123*^-4` → 123/10000. Must immediately
     * follow the mantissa (and any precision suffix), with no whitespace. */
    if (s->pos[0] == '*' && s->pos[1] == '^') {
        const char* save = s->pos;
        s->pos += 2;
        int exp_neg = 0;
        if (*s->pos == '+') s->pos++;
        else if (*s->pos == '-') { exp_neg = 1; s->pos++; }
        const char* exp_start = s->pos;
        while (isdigit((unsigned char)*s->pos)) s->pos++;
        if (s->pos == exp_start) {
            parse_err( "Expected exponent digits after '*^'\n");
            s->pos = save;
            expr_free(result);
            return NULL;
        }
        errno = 0;
        long exp_val = strtol(exp_start, NULL, 10);
        if (errno == ERANGE) {
            parse_err( "Exponent in '*^' suffix out of range\n");
            expr_free(result);
            return NULL;
        }
        if (exp_neg) exp_val = -exp_val;
        /* A machine-real mantissa is rescaled from its decimal text in one
         * correctly rounded strtod: d * 10^e rounds twice, so `1.5*^-5`
         * would not be the double nearest 1.5e-5 and printed forms of
         * reals would not re-read to the same value. */
        size_t mlen = 0;
        while (mantissa_start + mlen < save
               && (isdigit((unsigned char)mantissa_start[mlen])
                   || mantissa_start[mlen] == '.' || mantissa_start[mlen] == '-'
                   || mantissa_start[mlen] == '+'))
            mlen++;
        if (result && result->type == EXPR_REAL && mlen > 0 && mlen < 96
            && (mantissa_start[mlen] == '`' || mantissa_start + mlen == save)) {
            char tbuf[128];
            snprintf(tbuf, sizeof tbuf, "%.*se%ld", (int)mlen, mantissa_start, exp_val);
            double d = strtod(tbuf, NULL);
            expr_free(result);
            result = expr_new_real(d);
        } else {
            result = apply_pow10_scale(result, exp_val);
        }
    }

    return result;
}

/* Value of hex digit c, or -1. */
static int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

/* `ndig` hex digits at p as a number, or -1. */
static long hexnum(const char* p, int ndig) {
    long v = 0;
    for (int k = 0; k < ndig; k++) {
        int d = hexval(p[k]);
        if (d < 0) return -1;
        v = v * 16 + d;
    }
    return v;
}

/* The code point of the character escape at p (p[0] == '\\'): \[Name],
 * \:hhhh, \|hhhhhh, \.hh or \ooo, with the escape's source length in
 * *esc_len; -1 if p does not start such an escape.  Shared by string literals
 * and by escapes in source code (src_nonascii). */
static long char_escape_cp(const char* p, size_t* esc_len) {
    long cp = -1;
    size_t len = 0;
    if (p[1] == '[') {
        const char* close = p + 2;
        while (isalnum((unsigned char)*close)) close++;
        if (*close == ']' && close > p + 2) {
            cp = u8_named_char(p + 2, (size_t)(close - (p + 2)));
            len = (size_t)(close - p) + 1;
        }
    } else if (p[1] == ':') {
        cp = hexnum(p + 2, 4); len = 6;
    } else if (p[1] == '|') {
        cp = hexnum(p + 2, 6); len = 8;
    } else if (p[1] == '.') {
        cp = hexnum(p + 2, 2); len = 4;
    } else if (p[1] >= '0' && p[1] <= '7' && p[2] >= '0' && p[2] <= '7' &&
               p[3] >= '0' && p[3] <= '7') {
        cp = (p[1] - '0') * 64 + (p[2] - '0') * 8 + (p[3] - '0'); len = 4;
    }
    if (cp >= 0) *esc_len = len;
    return cp;
}

/* Decode the character escape at p to UTF-8 in enc (>= 4 bytes): writes the
 * escape's source length to *esc_len and returns the encoded length, or 0 if
 * p does not start such an escape. */
static size_t decode_char_escape(const char* p, size_t* esc_len, char* enc) {
    size_t len = 0;
    long cp = char_escape_cp(p, &len);
    if (cp < 0) return 0;
    size_t n = u8_encode((uint32_t)cp, enc);
    if (n) *esc_len = len;
    return n;
}

// Parses quoted strings ("text")
/* The buffer GROWS, like parse_list's: a string literal has no length limit.
 * It was a fixed char[256] whose overflow was silently DROPPED, so every
 * literal longer than 255 characters lost its tail with no diagnostic --
 * `StringLength["<300 c's>"]` answered 255.  Silent truncation of a literal is
 * a wrong answer by itself, and for the common case of a literal that holds an
 * expression (ToExpression[...], a generated .m file, a docstring) the cut
 * lands mid-token and surfaces instead as a bogus syntax error at a point the
 * source does not contain.  A string built at run time (StringJoin, <>) was
 * never affected, which is why this survived: only source literals truncate.
 * The stack buffer is kept as the fast path -- .m loading parses many short
 * strings and should not pay a malloc for each. */
static Expr* parse_string(ParserState* s) {
    s->pos++;  // Skip opening quote
    char stackbuf[256];
    char* buffer = stackbuf;
    size_t cap = sizeof(stackbuf);
    size_t i = 0;

    while (*s->pos && *s->pos != '"') {
        char c = *s->pos;
        /* Mathematica character escapes, decoded to UTF-8:
         *   \[Name]   named character (\[Alpha] -> U+03B1)
         *   \:hhhh    4-hex-digit code point
         *   \|hhhhhh  6-hex-digit code point
         *   \.hh      2-hex-digit code point
         *   \ooo      3-octal-digit code point
         * An escape that does not decode is kept verbatim. */
        char enc[4];
        size_t esc_len = 0, enc_len = 0;
        if (c == '\\' && (enc_len = decode_char_escape(s->pos, &esc_len, enc)) > 0) {
            if (i + enc_len >= cap) {
                size_t ncap = cap * 2;
                char* nb = (buffer == stackbuf) ? (char*)malloc(ncap)
                                                : (char*)realloc(buffer, ncap);
                if (!nb) {
                    if (buffer != stackbuf) free(buffer);
                    fprintf(stderr, "Out of memory parsing string literal\n");
                    return NULL;
                }
                if (buffer == stackbuf) memcpy(nb, stackbuf, i);
                buffer = nb;
                cap = ncap;
            }
            memcpy(buffer + i, enc, enc_len);
            i += enc_len;
            s->pos += esc_len;
            continue;
        }
        uint32_t box_cp;
        if (c == '\\' && (box_cp = lb_escape_cp(s->pos[1])) != 0) {
            /* \! \( \) \* ...: the box control characters, kept as the
             * reference's private-use code points so a string holding linear
             * box syntax still reads as boxes (ToExpression["\!\(x\)"]). */
            enc_len = u8_encode(box_cp, enc);
            if (i + enc_len >= cap) {
                size_t ncap = cap * 2;
                char* nb = (buffer == stackbuf) ? (char*)malloc(ncap)
                                                : (char*)realloc(buffer, ncap);
                if (!nb) {
                    if (buffer != stackbuf) free(buffer);
                    fprintf(stderr, "Out of memory parsing string literal\n");
                    return NULL;
                }
                if (buffer == stackbuf) memcpy(nb, stackbuf, i);
                buffer = nb;
                cap = ncap;
            }
            memcpy(buffer + i, enc, enc_len);
            i += enc_len;
            s->pos += 2;
            continue;
        }
        if (c == '\\' && (s->pos[1] == '\n' || s->pos[1] == '\r')) {
            /* Line continuation: a backslash at the end of a line drops both
             * characters, so "ab\<newline>cd" is "abcd".  Before a CRLF only
             * the CR goes with it, and the LF stays. */
            s->pos += 2;
            continue;
        }
        if (c == '\\' && s->pos[1] != '\0') {
            /* Decode standard C-style escapes. `\\` and `\"` decode to the
             * literal character (unchanged from the historical behaviour);
             * `\n`, `\t`, `\r` decode to the actual control characters
             * (previously mis-handled as the bare letter). Unknown escapes
             * keep the following character verbatim (backslash dropped). */
            s->pos++;
            switch (*s->pos) {
                case 'n': c = '\n'; break;
                case 't': c = '\t'; break;
                case 'r': c = '\r'; break;
                default:  c = *s->pos; break;
            }
        }
        if (i + 1 >= cap) {                       /* room for c and the NUL */
            size_t ncap = cap * 2;
            char* nb = (buffer == stackbuf) ? (char*)malloc(ncap)
                                            : (char*)realloc(buffer, ncap);
            if (!nb) {
                if (buffer != stackbuf) free(buffer);
                fprintf(stderr, "Out of memory parsing string literal\n");
                return NULL;
            }
            if (buffer == stackbuf) memcpy(nb, stackbuf, i);
            buffer = nb;
            cap = ncap;
        }
        buffer[i++] = c;
        s->pos++;
    }

    if (*s->pos != '"') {
        if (buffer != stackbuf) free(buffer);
        parse_err( "Unterminated string\n");
        return NULL;
    }
    s->pos++;  // Skip closing quote

    buffer[i] = '\0';
    Expr* out = expr_new_string(buffer);
    if (buffer != stackbuf) free(buffer);
    return out;
}

/* ------------------- Compound Expressions ------------------- */

// Parses lists: {1,2,3}
/* Parses the comma-separated elements of `{...}` / `f[...]` up to (not past)
 * the closing character `close`; the opening bracket has been consumed.
 * Mathematica's rule for empty slots: an element omitted between, before or
 * after a comma is Null -- `f[,x]` is f[Null, x], `f[x,]` is f[x, Null],
 * `{a,,b}` is {a, Null, b} -- while an empty pair `f[]` has no elements.
 * Returns a malloc'd array (count in *out_count) with s->pos on `close`, or
 * NULL on a syntax error with every parsed element freed. A missing closer
 * is reported here, once, as "Expected '<close>' to close <what>"; a failure
 * inside an element has already been reported by the inner parse. */
static Expr** parse_sequence(ParserState* s, char close, const char* what,
                             size_t* out_count) {
    size_t cap = 16, count = 0;
    Expr** elems = malloc(cap * sizeof(Expr*));

    skip_whitespace(s);
    if (*s->pos != close) {
        for (;;) {
            skip_whitespace(s);
            Expr* e;
            if (*s->pos == ',' || *s->pos == close) {
                e = expr_new_symbol(SYM_Null);          /* empty slot */
            } else {
                e = parse_expression_state(s);
                if (!e) goto fail;
                skip_whitespace(s);
            }
            if (count >= cap) {
                cap *= 2;
                elems = realloc(elems, cap * sizeof(Expr*));
            }
            elems[count++] = e;

            if (*s->pos == ',') { s->pos++; continue; }
            if (*s->pos == close) break;
            parse_err( "Expected '%c' to close %s\n", close, what);
            goto fail;                                  /* junk or EOF */
        }
    }
    *out_count = count;
    return elems;

fail:
    while (count--) expr_free(elems[count]);
    free(elems);
    return NULL;
}

static Expr* parse_list(ParserState* s) {
    s->pos++;  // Skip '{'
    s->bracket_depth++;
    size_t count = 0;
    Expr** elements = parse_sequence(s, '}', "list", &count);
    if (!elements) return NULL;
    s->pos++;  // Skip '}'
    close_bracket(s);

    // Create List[...] expression
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elements, count);
    free(elements);
    return list;
}

/* Parses associations: <|a -> 1, b -> 2|>  ->  Association[Rule[a,1], Rule[b,2]].
 * Each entry is parsed as an ordinary expression; because `->`/`:>` are normal
 * operators (precedence 1500) the `k -> v` shape becomes a Rule automatically.
 * The `<|` and `|>` delimiters are recognised as OP_NONE by get_operator, so
 * value parsing halts exactly at `|>`. Semantic validation (entries must be
 * rules, duplicate-key collapse) is left to builtin_association. */
static Expr* parse_association(ParserState* s) {
    s->pos += 2;  // Skip '<|'
    s->bracket_depth++;
    size_t cap = 16;
    Expr** entries = malloc(cap * sizeof(Expr*));
    size_t count = 0;

    skip_whitespace(s);
    while (*s->pos && strncmp(s->pos, "|>", 2) != 0) {
        if (*s->pos == ',') s->pos++;
        skip_whitespace(s);
        if (strncmp(s->pos, "|>", 2) == 0) break;  /* tolerate trailing comma */

        Expr* entry = parse_expression_state(s);
        if (!entry) {
            while (count--) expr_free(entries[count]);
            free(entries);
            return NULL;
        }
        if (count >= cap) {
            cap *= 2;
            entries = realloc(entries, cap * sizeof(Expr*));
        }
        entries[count++] = entry;
        skip_whitespace(s);
    }

    if (strncmp(s->pos, "|>", 2) != 0) {
        parse_err( "Expected '|>' to close association\n");
        while (count--) expr_free(entries[count]);
        free(entries);
        return NULL;
    }
    s->pos += 2;  // Skip '|>'
    close_bracket(s);

    Expr* assoc = expr_new_function(expr_new_symbol(SYM_Association), entries, count);
    free(entries);
    return assoc;
}

// Parses functions: f[x,y]
static Expr* parse_function(ParserState* s, Expr* head) {
    s->pos++;  // Skip '['
    s->bracket_depth++;
    size_t count = 0;
    Expr** args = parse_sequence(s, ']', "function", &count);
    if (!args) {
        expr_free(head);  /* we own the head; the caller drops it on NULL */
        return NULL;
    }
    s->pos++;  // Skip ']'
    close_bracket(s);

    if (head && head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Sqrt && count == 1) {
        expr_free(head);
        Expr* rat_args[2] = { expr_new_integer(1), expr_new_integer(2) };
        Expr* half = expr_new_function(expr_new_symbol(SYM_Rational), rat_args, 2);
        Expr* pow_args[2] = { args[0], half };
        Expr* result = expr_new_function(expr_new_symbol(SYM_Power), pow_args, 2);
        free(args);
        return result;
    }

    Expr* func = expr_new_function(head, args, count);
    free(args);
    return func;
}

/* ------------------- Main Parser ------------------- */

/* `:` precedences, in the reference's relative order
 * (Pattern 150 > Optional 140 > StringExpression 135 > Condition 130):
 *   sym : p   Pattern  -- above `~~` (2100), below `|` (2300) and `..`
 *                         (2500), so `x:a|b` and `x:p..` bind into the
 *                         pattern while `x:a ~~ b` splits at `~~`;
 *   p : v     Optional -- below Pattern, so `x:_:1` is
 *                         Optional[Pattern[x,_],1].
 * The right operand is read at prec+1 (left-associative). */
#define PREC_PATTERN  2200
#define PREC_OPTIONAL 2150

static Expr* parse_expression_prec(ParserState* s, int min_prec);

typedef enum {
    OP_NONE = 0,
    OP_PLUS,
    OP_MINUS,
    OP_TIMES,
    OP_DIVIDE,
    OP_DOT,
    OP_POWER,
    OP_SET,
    OP_SETDELAYED,
    OP_UNSET,
    OP_EQUAL,
    OP_UNEQUAL,
    OP_LESS,
    OP_GREATER,
    OP_LESSEQUAL,
    OP_GREATEREQUAL,
    OP_SAMEQ,
    OP_UNSAMEQ,
    OP_APPLY,
    OP_APPLY1,
    OP_RULE,
    OP_RULEDELAYED,
    OP_TWOWAYRULE,
    OP_CONDITION,
    OP_ALTERNATIVES,
    OP_STRINGEXPRESSION,
    OP_MAP,
    OP_MAPALL,
    OP_REPLACEALL,
    OP_REPLACEREPEATED,
    OP_SPAN,
    OP_PART,
    OP_CALL,
    OP_AND,
    OP_OR,
    OP_COMPOUND,
    OP_FUNCTION,
    OP_POSTFIX,
    OP_PREFIX,
    OP_PATTERNTEST,
    OP_COLON,
    OP_MESSAGENAME,
    OP_INFORMATION,
    OP_FACTORIAL,
    OP_FACTORIAL2,
    OP_REPEATED,
    OP_REPEATEDNULL,
    OP_STRINGJOIN,
    OP_INCREMENT,
    OP_DECREMENT,
    OP_ADDTO,
    OP_SUBTRACTFROM,
    OP_TIMESBY,
    OP_DIVIDEBY,
    OP_DERIVATIVE,
    OP_COMPOSITION,
    OP_PUT,
    OP_PUTAPPEND,
    OP_NONCOMMULT,    /* a ** b ** c -> NonCommutativeMultiply[a, b, c] */
    OP_NAMED_INFIX,   /* a \[CenterDot] b: a named-character operator with no
                         ASCII spelling; head_name is its head */
    OP_NAMED_POSTFIX, /* a \[Transpose] -> Transpose[a] */
#if UP_VALUES
    OP_UPSET,         /* lhs ^= rhs   */
    OP_UPSETDELAYED,  /* lhs ^:= rhs  */
    OP_TAGSET         /* f /: lhs = / := / =. rhs  (ternary, dedicated led branch) */
#endif
} OperatorType;

typedef struct {
    OperatorType type;
    int prec;
    int right_assoc;
    const char* head_name;
    int len;
    int flat;          /* OP_NAMED_INFIX: a<op>b<op>c is one n-ary head */
} OperatorDef;

/* The ASCII spelling of the operators a named character can also be written
 * with, so \[Rule], \[LessEqual], \[Times], ... reuse their ASCII forms'
 * parsing (Inequality chains, Times flattening, ...). */
static const char* named_operator_ascii(const char* head) {
    static const struct { const char* head; const char* ascii; } k[] = {
        {"Rule", "->"}, {"RuleDelayed", ":>"}, {"Equal", "=="}, {"Unequal", "!="},
        {"Less", "<"}, {"Greater", ">"}, {"LessEqual", "<="}, {"GreaterEqual", ">="},
        {"Times", "*"}, {"Plus", "+"}, {"And", "&&"}, {"Or", "||"}, {"Dot", "."},
        {"Power", "^"}, {"Alternatives", "|"}, {"TwoWayRule", "<->"},
    };
    for (size_t i = 0; i < sizeof k / sizeof k[0]; i++)
        if (strcmp(k[i].head, head) == 0) return k[i].ascii;
    return NULL;
}

static OperatorDef get_operator(const char* pos) {
    OperatorDef def = {OP_NONE, -1, 0, NULL, 0, 0};

    /* A named-character operator, raw (U+00D7) or escaped (\[Times]).  The
     * raw U+2192 is excluded: it stays the Rule synonym matched below. */
    if (((unsigned char)pos[0] >= 0x80 || pos[0] == '\\') &&
        !((unsigned char)pos[0] == 0xE2 && (unsigned char)pos[1] == 0x86 &&
          (unsigned char)pos[2] == 0x92)) {
        size_t n = 0;
        const U8CharInfo* ci = src_special(pos, &n);
        if (ci && (ci->cls == U8C_INFIX || ci->cls == U8C_POSTFIX) && ci->head) {
            const char* ascii = ci->cls == U8C_INFIX ? named_operator_ascii(ci->head) : NULL;
            if (ascii) {
                def = get_operator(ascii);
            } else {
                def.type = ci->cls == U8C_INFIX ? OP_NAMED_INFIX : OP_NAMED_POSTFIX;
                def.prec = ci->prec;
                def.right_assoc = ci->assoc == U8A_RIGHT;
                def.flat = ci->assoc == U8A_FLAT;
                def.head_name = ci->head;
            }
            def.len = (int)n;
        }
        return def;
    }

    /* Association delimiters `<|` and `|>` are never operators: they are
     * handled structurally by parse_association. Returning OP_NONE here stops
     * the Pratt binary loop cleanly at a `|>` (rather than lexing `|` as
     * Alternatives / `>` as Greater) and prevents `<|` being read as Less. */
    if (strncmp(pos, "<|", 2) == 0 || strncmp(pos, "|>", 2) == 0) return def;
    /* `<<file` is the prefix Get form, read in parse_primary; never Less Less. */
    if (strncmp(pos, "<<", 2) == 0) return def;

    if (strncmp(pos, ">>>", 3) == 0) {
        /* PutAppend (expr >>> file). Lower than Set (500) so `a = 5 >>> "f"`
         * parses as Set[a, PutAppend[5, "f"]] — matches Mathematica. */
        def.type = OP_PUTAPPEND; def.prec = 300; def.head_name = "PutAppend"; def.len = 3;
    } else if (strncmp(pos, "===", 3) == 0) {
        def.type = OP_SAMEQ; def.prec = 3200; def.head_name = "SameQ"; def.len = 3;
    } else if (strncmp(pos, "=!=", 3) == 0) {
        def.type = OP_UNSAMEQ; def.prec = 3200; def.head_name = "UnsameQ"; def.len = 3;
    } else if (strncmp(pos, "@@@", 3) == 0) {
        def.type = OP_APPLY1; def.prec = 7000; def.right_assoc = 1; def.head_name = "Apply1"; def.len = 3;
    } else if (strncmp(pos, "//.", 3) == 0 && !isdigit(pos[3])) {
        def.type = OP_REPLACEREPEATED; def.prec = 1200; def.right_assoc = 0; def.head_name = "ReplaceRepeated"; def.len = 3;
    } else if (strncmp(pos, "//@", 3) == 0) {
        def.type = OP_MAPALL; def.prec = 7000; def.right_assoc = 1; def.head_name = "MapAll"; def.len = 3;
    } else if (strncmp(pos, "<->", 3) == 0) {
        /* TwoWayRule (u <-> v). Same precedence/associativity as Rule; the
         * graph subsystem normalizes it to UndirectedEdge on construction.
         * Checked before "<>", "<=", "<" so the 3-char form wins. */
        def.type = OP_TWOWAYRULE; def.prec = 1500; def.right_assoc = 1; def.head_name = "TwoWayRule"; def.len = 3;
    } else if (strncmp(pos, "//", 2) == 0) {
        def.type = OP_POSTFIX; def.prec = 800; def.head_name = "Postfix"; def.len = 2;
    } else if (strncmp(pos, "/.", 2) == 0 && !isdigit(pos[2])) {
        def.type = OP_REPLACEALL; def.prec = 1200; def.right_assoc = 0; def.head_name = "ReplaceAll"; def.len = 2;
    } else if (strncmp(pos, "@@", 2) == 0) {
        def.type = OP_APPLY; def.prec = 7000; def.right_assoc = 1; def.head_name = "Apply"; def.len = 2;
    } else if (strncmp(pos, "@*", 2) == 0) {
        def.type = OP_COMPOSITION; def.prec = 7200; def.right_assoc = 0; def.head_name = "Composition"; def.len = 2;
    } else if (strncmp(pos, "/@", 2) == 0) {
        def.type = OP_MAP; def.prec = 7000; def.right_assoc = 1; def.head_name = "Map"; def.len = 2;
    } else if (strncmp(pos, ":>", 2) == 0) {
        def.type = OP_RULEDELAYED; def.prec = 1500; def.right_assoc = 1; def.head_name = "RuleDelayed"; def.len = 2;
    } else if (strncmp(pos, "->", 2) == 0) {
        def.type = OP_RULE; def.prec = 1500; def.right_assoc = 1; def.head_name = "Rule"; def.len = 2;
    } else if ((unsigned char)pos[0] == 0xE2 && (unsigned char)pos[1] == 0x86 &&
               (unsigned char)pos[2] == 0x92) {
        /* Unicode → (U+2192, Wolfram \[Rule]) is a synonym for `->`, so pasted
         * Wolfram-Language rules and associations parse directly. Three UTF-8
         * bytes E2 86 92. */
        def.type = OP_RULE; def.prec = 1500; def.right_assoc = 1; def.head_name = "Rule"; def.len = 3;
    } else if (strncmp(pos, "/;", 2) == 0) {
        def.type = OP_CONDITION; def.prec = 1700; def.right_assoc = 0; def.head_name = "Condition"; def.len = 2;
    } else if (strncmp(pos, ";;", 2) == 0) {
        def.type = OP_SPAN; def.prec = 3200; def.right_assoc = 0; def.head_name = "Span"; def.len = 2;
    } else if (*pos == ';') {
        def.type = OP_COMPOUND; def.prec = 100; def.right_assoc = 0; def.head_name = "CompoundExpression"; def.len = 1;
    } else if (strncmp(pos, "~~", 2) == 0) {
        /* StringExpression (a ~~ b ~~ c): concatenation of string patterns.
         * Precedence just below Alternatives (2300) so `a ~~ b | c` groups as
         * StringExpression[a, Alternatives[b, c]], matching Mathematica. The
         * StringExpression head is Flat, so nested nodes flatten at eval. */
        def.type = OP_STRINGEXPRESSION; def.prec = 2100; def.head_name = "StringExpression"; def.len = 2;
    } else if (strncmp(pos, "||", 2) == 0) {
        /* Or binds LOOSER than And (2800): `a || b && c` is `a || (b && c)`. */
        def.type = OP_OR; def.prec = 2700; def.head_name = "Or"; def.len = 2;
    } else if (*pos == '|') {
        def.type = OP_ALTERNATIVES; def.prec = 2300; def.head_name = "Alternatives"; def.len = 1;
    } else if (strncmp(pos, "==", 2) == 0) {
        def.type = OP_EQUAL; def.prec = 3200; def.head_name = "Equal"; def.len = 2;
    } else if (strncmp(pos, "&&", 2) == 0) {
        def.type = OP_AND; def.prec = 2800; def.head_name = "And"; def.len = 2;
    } else if (*pos == '&') {
        def.type = OP_FUNCTION; def.prec = 1000; def.head_name = "Function"; def.len = 1;
    } else if (strncmp(pos, "!=", 2) == 0) {
        def.type = OP_UNEQUAL; def.prec = 3200; def.head_name = "Unequal"; def.len = 2;
    } else if (strncmp(pos, "<>", 2) == 0) {
        def.type = OP_STRINGJOIN; def.prec = 6700; def.head_name = "StringJoin"; def.len = 2;
    } else if (strncmp(pos, "<=", 2) == 0) {
        def.type = OP_LESSEQUAL; def.prec = 3200; def.head_name = "LessEqual"; def.len = 2;
    } else if (strncmp(pos, ">=", 2) == 0) {
        def.type = OP_GREATEREQUAL; def.prec = 3200; def.head_name = "GreaterEqual"; def.len = 2;
    } else if (strncmp(pos, ">>", 2) == 0) {
        /* Put (expr >> file). Same precedence as PutAppend; right side is
         * a filename token (bare identifier or quoted string), parsed
         * specially in the operator loop below. */
        def.type = OP_PUT; def.prec = 300; def.head_name = "Put"; def.len = 2;
    } else if (*pos == '<') {
        def.type = OP_LESS; def.prec = 3200; def.head_name = "Less"; def.len = 1;
    } else if (*pos == '>') {
        def.type = OP_GREATER; def.prec = 3200; def.head_name = "Greater"; def.len = 1;
    } else if (strncmp(pos, ":=", 2) == 0) {
        def.type = OP_SETDELAYED; def.prec = 500; def.right_assoc = 1; def.head_name = "SetDelayed"; def.len = 2;
    } else if (strncmp(pos, "[[", 2) == 0) {
        def.type = OP_PART; def.prec = 10000; def.head_name = "Part"; def.len = 2;
    } else if (strncmp(pos, "::", 2) == 0) {
        /* MessageName: f::tag. High precedence so it binds tighter than
         * function application's operands; the tag is read literally (as a
         * string) by the binary loop, not parsed as an expression. */
        def.type = OP_MESSAGENAME; def.prec = 8600; def.head_name = "MessageName"; def.len = 2;
    } else if (*pos == ':') {
        /* `:` normally never reaches the binary loop: a symbol or blank
         * token followed by `:` is consumed at atom level by
         * parse_colon_suffix. What remains is the lenient fallback for a
         * left operand the reference rejects (`(a):_g`, `f[x]:0`, `x_.:0`); the binary
         * loop re-rates it once `left` is known (see the OP_COLON fix-up
         * there). The value here is only the Pattern rating. */
        def.type = OP_COLON; def.prec = PREC_PATTERN; def.right_assoc = 0; def.head_name = "Pattern"; def.len = 1;
    } else if (strncmp(pos, "=.", 2) == 0 && !isdigit((unsigned char)pos[2])) {
        /* Unset (postfix `lhs =.`). Guarded against `=.5` etc. (a real
         * literal on the RHS of Set) by requiring the char after `.` not be
         * a digit. Low precedence, like Set, so it captures the whole
         * preceding expression: `a b =.` -> Unset[a b]. */
        def.type = OP_UNSET; def.prec = 500; def.head_name = "Unset"; def.len = 2;
    } else if (*pos == '=') {
        def.type = OP_SET; def.prec = 500; def.right_assoc = 1; def.head_name = "Set"; def.len = 1;
    } else if (strncmp(pos, "++", 2) == 0) {
        def.type = OP_INCREMENT; def.prec = 7500; def.head_name = "Increment"; def.len = 2;
    } else if (strncmp(pos, "--", 2) == 0) {
        def.type = OP_DECREMENT; def.prec = 7500; def.head_name = "Decrement"; def.len = 2;
    } else if (strncmp(pos, "+=", 2) == 0) {
        def.type = OP_ADDTO; def.prec = 500; def.right_assoc = 1; def.head_name = "AddTo"; def.len = 2;
    } else if (strncmp(pos, "-=", 2) == 0) {
        def.type = OP_SUBTRACTFROM; def.prec = 500; def.right_assoc = 1; def.head_name = "SubtractFrom"; def.len = 2;
    } else if (strncmp(pos, "**", 2) == 0) {
        /* NonCommutativeMultiply. Associative (Flat) but non-commutative.
         * Precedence 5900 sits in the empty gap between Dot (5300) and
         * Power (6500), reproducing Mathematica's ordering
         * Plus < Times < Dot < NonCommutativeMultiply < Power.
         * Must precede the bare `*` arm below; differs from `*=` at offset 1
         * so order among the `*` forms is otherwise free. Left-associative. */
        def.type = OP_NONCOMMULT; def.prec = 5900; def.head_name = "NonCommutativeMultiply"; def.len = 2;
    } else if (strncmp(pos, "*=", 2) == 0) {
        def.type = OP_TIMESBY; def.prec = 500; def.right_assoc = 1; def.head_name = "TimesBy"; def.len = 2;
    } else if (strncmp(pos, "/=", 2) == 0) {
        def.type = OP_DIVIDEBY; def.prec = 500; def.right_assoc = 1; def.head_name = "DivideBy"; def.len = 2;
    } else if (*pos == '+') {
        def.type = OP_PLUS; def.prec = 3500; def.head_name = "Plus"; def.len = 1;
    } else if (*pos == '-') {
        def.type = OP_MINUS; def.prec = 3500; def.head_name = "Plus"; def.len = 1;
    } else if (*pos == '*') {
        def.type = OP_TIMES; def.prec = 4500; def.head_name = "Times"; def.len = 1;
#if UP_VALUES
    } else if (strncmp(pos, "/:", 2) == 0) {
        /* TagSet family: f /: lhs = / := / =. rhs. A loose, statement-level
         * ternary handled by a dedicated led branch below; prec 500 only needs to
         * be low enough to take the tag on its left and stop the binary loop here.
         * Differs from /., /;, /@, /=, //, //. at offset 1, so order among the
         * `/` forms is free -- it just has to precede the bare `/` (Divide). */
        def.type = OP_TAGSET; def.prec = 500; def.head_name = "TagSet"; def.len = 2;
#endif
    } else if (*pos == '/') {
        def.type = OP_DIVIDE; def.prec = 5000; def.head_name = "Divide"; def.len = 1;
    } else if (strncmp(pos, "...", 3) == 0) {
        def.type = OP_REPEATEDNULL; def.prec = 2500; def.head_name = "RepeatedNull"; def.len = 3;
    } else if (strncmp(pos, "..", 2) == 0) {
        def.type = OP_REPEATED; def.prec = 2500; def.head_name = "Repeated"; def.len = 2;
    } else if (*pos == '.' && !isdigit(pos[1])) {
        def.type = OP_DOT; def.prec = 5300; def.head_name = "Dot"; def.len = 1;
#if UP_VALUES
    } else if (strncmp(pos, "^:=", 3) == 0) {
        /* UpSetDelayed. Same tier as SetDelayed (500). MUST precede the bare `^`
         * test below, or `^:=` would lex as Power then SetDelayed. */
        def.type = OP_UPSETDELAYED; def.prec = 500; def.right_assoc = 1; def.head_name = "UpSetDelayed"; def.len = 3;
    } else if (strncmp(pos, "^=", 2) == 0) {
        /* UpSet. Same tier as Set (500). Precede the bare `^`. */
        def.type = OP_UPSET; def.prec = 500; def.right_assoc = 1; def.head_name = "UpSet"; def.len = 2;
#endif
    } else if (*pos == '^') {
        def.type = OP_POWER; def.prec = 6500; def.right_assoc = 1; def.head_name = "Power"; def.len = 1;
    } else if (*pos == '?') {
        def.type = OP_PATTERNTEST; def.prec = 7900; def.right_assoc = 0; def.head_name = "PatternTest"; def.len = 1;
    } else if (*pos == '@') {
        def.type = OP_PREFIX; def.prec = 7000; def.right_assoc = 1; def.head_name = "Prefix"; def.len = 1;
    } else if (*pos == '[') {
        def.type = OP_CALL; def.prec = 9500; def.len = 1;
    } else if (*pos == '!' && pos[1] == '!') {
        def.type = OP_FACTORIAL2; def.prec = 8200; def.head_name = "Factorial2"; def.len = 2;
    } else if (*pos == '!' && pos[1] != '=') {
        def.type = OP_FACTORIAL; def.prec = 8200; def.head_name = "Factorial"; def.len = 1;
    } else if (*pos == '\'') {
        def.type = OP_DERIVATIVE; def.prec = 7700; def.head_name = "Derivative"; def.len = 1;
    }

    return def;
}

// Parses primary expressions (atoms, lists, functions, parens)
/* True when `p` is a lone `:` -- not `::` (MessageName), `:=` (SetDelayed)
 * or `:>` (RuleDelayed). */
static bool is_single_colon(const char* p) {
    return p[0] == ':' && p[1] != ':' && p[1] != '=' && p[1] != '>';
}

typedef enum { COLON_LHS_SYMBOL, COLON_LHS_BLANK } ColonLhsKind;

/* Attach a `:` suffix to the token just read, at ATOM level.
 *
 * The left operand of `:` is the single token
 * immediately before it -- a bare symbol (`a:p` is Pattern) or a blank
 * token (`x_:v` is Optional) -- not the whole expression built so far. So
 * `n_.*a:_g` is Times[n_., a:_g], `n x_:0` is Times[n, x_:0], `x+a:_g==b`
 * is Plus[x, a:(_g==b)] and `x?f:_` is PatternTest[x, f:_]. Reading `:` as
 * an ordinary low-precedence infix operator instead built
 * Optional[n_.*a, _g] -- a silently wrong pattern (issue #102).
 *
 * While `p` is read, a `:` at the same bracket level is not claimed by an
 * inner token (see ParserState.colon_hold): it ends `p`. So after `sym:p`
 * one further `:` supplies a default (`a:_g:1` is
 * Optional[a:_g, 1]); the default itself may again start with `sym:`
 * (`x:a:b:c:d` is Optional[x:a, Optional[b:c, d]], read recursively through
 * parse_primary). Consumes `tok`; returns NULL (with `tok` freed) on a
 * malformed right operand. If no `:` follows, `tok` is returned unchanged
 * and the position is left where it was. */
static Expr* parse_colon_suffix(ParserState* s, Expr* tok, ColonLhsKind kind) {
    if (s->colon_hold == s->bracket_depth + 1) return tok;
    const char* save_pos = s->pos;
    int save_newline = s->saw_newline;
    skip_whitespace(s);
    if (!is_single_colon(s->pos)) {
        s->pos = save_pos;
        s->saw_newline = save_newline;
        return tok;
    }
    s->pos++;
    if (kind == COLON_LHS_SYMBOL) {
        int save_hold = s->colon_hold;
        s->colon_hold = s->bracket_depth + 1;
        Expr* patt = parse_expression_prec(s, PREC_PATTERN + 1);
        s->colon_hold = save_hold;
        if (!patt) { expr_free(tok); return NULL; }
        Expr* args[2] = { tok, patt };
        tok = expr_new_function(expr_new_symbol(SYM_Pattern), args, 2);

        save_pos = s->pos;
        save_newline = s->saw_newline;
        skip_whitespace(s);
        if (!is_single_colon(s->pos)) {
            s->pos = save_pos;
            s->saw_newline = save_newline;
            return tok;
        }
        s->pos++;
    }
    Expr* def = parse_expression_prec(s, PREC_OPTIONAL + 1);
    if (!def) { expr_free(tok); return NULL; }
    Expr* args[2] = { tok, def };
    return expr_new_function(expr_new_symbol(SYM_Optional), args, 2);
}

/* The file name after `<<`, `>>` or `>>>`: a string literal, or a bare token
 * of letters, digits, `` !$*-./:?\_`~ `` and balanced `[...]` (so `<<f[x]` is
 * Get["f[x]"]).  A bare token ends at whitespace, `;`, `,`, `)`, `]`, `}` or
 * end of input; any other character there is a syntax error (`<<f+1`).
 * Returns a String, or NULL after reporting the error. */
static bool file_token_char(char c) {
    return isalnum((unsigned char)c) || (unsigned char)c >= 0x80 ||
           strchr("!$*-./:?\\_`~", c) != NULL;
}

static Expr* parse_file_token(ParserState* s, const char* who) {
    skip_whitespace(s);
    if (*s->pos == '"') return parse_string(s);
    const char* start = s->pos;
    const char* p = s->pos;
    int depth = 0;
    for (;;) {
        char c = *p;
        if (c && file_token_char(c)) { p++; continue; }
        if (c == '[') { depth++; p++; continue; }
        if (c == ']' && depth > 0) { depth--; p++; continue; }
        break;
    }
    bool ends_cleanly = *p == '\0' || isspace((unsigned char)*p) || strchr(";,)]}", *p) != NULL;
    if (p == start || depth != 0 || !ends_cleanly) {
        parse_err("Expected a file name after %s\n", who);
        return NULL;
    }
    size_t n = (size_t)(p - start);
    char* buf = malloc(n + 1);
    if (!buf) return NULL;
    memcpy(buf, start, n);
    buf[n] = '\0';
    s->pos = p;
    Expr* out = expr_new_string(buf);
    free(buf);
    return out;
}

static Expr* parse_primary(ParserState* s) {
    skip_whitespace(s);
    if (!*s->pos) return NULL;
    /* A primary token starts here: reset the sticky newline flag so that only
     * line breaks appearing *after* this token can terminate the statement. */
    s->saw_newline = 0;

    /* \!\( ... \) / \( ... \): linear box syntax, lowered to InputForm text
     * and read as the expression the boxes display (linear_boxes.c). */
    if (lb_group_open(s->pos)) {
        const char* end = NULL;
        char* text = lb_group_to_input(s->pos, &end);
        Expr* e = text ? parse_expression(text) : NULL;
        free(text);
        if (!e) {
            parse_err("Invalid linear box syntax\n");
            return NULL;
        }
        s->pos = end;
        return e;
    }

    switch (*s->pos) {
        case '"': return parse_string(s);
        case '{': return parse_list(s);
        case '<':
            /* `<<file`: Get["file"].  The name is a token, not an expression. */
            if (s->pos[1] == '<') {
                s->pos += 2;
                Expr* name = parse_file_token(s, "<<");
                if (!name) return NULL;
                Expr* a[1] = { name };
                return expr_new_function(expr_new_symbol(SYM_Get), a, 1);
            }
            parse_err("Unexpected character: '<'\n");
            return NULL;
        /* A bare '[' cannot start an expression: there is no head to apply.
         * Mathematica reports `f[x,[x]]` as a syntax error, and so do we.
         * (This used to build an EXPR_FUNCTION with a NULL head, which the
         * evaluator then dereferenced -- a segfault.) `g[[x]]` (Part) and
         * `f[x]` never reach here: both are consumed as postfix operators on
         * an already-parsed head. */
        case '[': return NULL;
        case '(': {
            s->pos++;
            s->bracket_depth++;
            Expr* inner = parse_expression_prec(s, 0);
            skip_whitespace(s);
            if (*s->pos == ')') {
                s->pos++;
            } else {
                /* An unclosed "(" is a syntax error, not an expression that
                 * happens to end where the input does: "(2 + x" must not read
                 * as 2 + x. */
                parse_err( "Expected ')'\n");
                close_bracket(s);
                expr_free(inner);
                return NULL;
            }
            close_bracket(s);
            return inner;
        }

        case '%': {
            s->pos++;
            int count = 1;
            while (*s->pos == '%') {
                count++;
                s->pos++;
            }
            if (isdigit(*s->pos)) {
                char* end;
                int64_t n = strtoll(s->pos, &end, 10);
                s->pos = end;
                Expr* out_args[1] = { expr_new_integer(n) };
                return expr_new_function(expr_new_symbol(SYM_Out), out_args, 1);
            } else {
                Expr* out_args[1] = { expr_new_integer(-count) };
                return expr_new_function(expr_new_symbol(SYM_Out), out_args, 1);
            }
        }
        
        case '#': {
            s->pos++;
            int count = 1;
            while (*s->pos == '#') {
                count++;
                s->pos++;
            }
            /* Named slots (single '#' only, as in Mathematica):
             *   #name      -> Slot["name"]   name = [A-Za-z$][A-Za-z0-9$`]*
             *   #"string"  -> Slot["string"] any string literal
             * `##name` stays SlotSequence[1]*name and `#1a` stays Slot[1]*a,
             * matching the Wolfram Language lexer. */
            if (count == 1 && (isalpha((unsigned char)*s->pos) || *s->pos == '$')) {
                const char* start = s->pos;
                while (isalnum((unsigned char)*s->pos) || *s->pos == '$' ||
                       (*s->pos == '`' && (isalpha((unsigned char)s->pos[1]) || s->pos[1] == '$')))
                    s->pos++;
                size_t len = (size_t)(s->pos - start);
                char* name = malloc(len + 1);
                memcpy(name, start, len);
                name[len] = '\0';
                Expr* args[1] = { expr_new_string(name) };
                free(name);
                return expr_new_function(expr_new_symbol(SYM_Slot), args, 1);
            }
            if (count == 1 && *s->pos == '"') {
                Expr* str = parse_string(s);
                if (!str) return NULL;
                Expr* args[1] = { str };
                return expr_new_function(expr_new_symbol(SYM_Slot), args, 1);
            }
            int64_t n = 1; // Default
            if (isdigit(*s->pos)) {
                char* end;
                n = strtoll(s->pos, &end, 10);
                s->pos = end;
            }
            if (count == 1) {
                Expr* args[1] = { expr_new_integer(n) };
                return expr_new_function(expr_new_symbol(SYM_Slot), args, 1);
            } else {
                Expr* args[1] = { expr_new_integer(n) };
                return expr_new_function(expr_new_symbol(SYM_SlotSequence), args, 1);
            }
        }
        
        default:
            if (isalpha((unsigned char)*s->pos) || *s->pos == '`' || *s->pos == '$' ||
                *s->pos == '_' || src_letter(s->pos)) {
                Expr* head = NULL;
                if (*s->pos != '_') {
                    head = parse_symbol(s);
                }
                
                int underscores = 0;
                while (*s->pos == '_') {
                    underscores++;
                    s->pos++;
                }
                
                if (underscores > 0) {
                    Expr* blank = NULL;
                    Expr* blank_head = NULL;
                    if (isalpha((unsigned char)*s->pos) || *s->pos == '`' || *s->pos == '$' ||
                        src_letter(s->pos)) {
                        blank_head = parse_symbol(s);
                    }
                    
                    int n_args = blank_head ? 1 : 0;
                    Expr* args_b[1];
                    if (blank_head) args_b[0] = blank_head;
                    
                    if (underscores == 1) {
                        blank = expr_new_function(expr_new_symbol(SYM_Blank), blank_head ? args_b : NULL, n_args);
                    } else if (underscores == 2) {
                        blank = expr_new_function(expr_new_symbol(SYM_BlankSequence), blank_head ? args_b : NULL, n_args);
                    } else {
                        blank = expr_new_function(expr_new_symbol(SYM_BlankNullSequence), blank_head ? args_b : NULL, n_args);
                    }
                    
                    if (head) {
                        Expr* args[2] = { head, blank };
                        head = expr_new_function(expr_new_symbol(SYM_Pattern), args, 2);
                    } else {
                        head = blank;
                    }
                    
                    /* `x_.` is Optional. As in the Wolfram Language the
                     * `_.` token is only formed by a BARE single blank, so
                     * `_...` is (_.).. and `_..` is a syntax error; after a
                     * head or a sequence blank the dots are the `..` / `...`
                     * operator: `_Integer..` is Repeated[_Integer], `x__..`
                     * is Repeated[x__].  (A lone `_h.` / `x__.` is still read
                     * as Optional, which WL rejects -- kept for leniency.) */
                    bool bare_single = (underscores == 1 && !blank_head);
                    if (*s->pos == '.' && (bare_single || s->pos[1] != '.')) {
                        s->pos++;
                        Expr* opt_args[1] = { head };
                        head = expr_new_function(expr_new_symbol(SYM_Optional), opt_args, 1);
                        /* `x_.` already carries its default; a following `:`
                         * is left to the binary loop's lenient fallback. */
                        return head;
                    }
                    return parse_colon_suffix(s, head, COLON_LHS_BLANK);
                }

                return parse_colon_suffix(s, head, COLON_LHS_SYMBOL);

            }
            if (isdigit(*s->pos) || ((*s->pos == '-' || *s->pos == '+') && (isdigit(s->pos[1]) || (s->pos[1] == '.' && isdigit(s->pos[2])))) || (*s->pos == '.' && isdigit(s->pos[1]))) {
                return parse_number(s);
            }
            parse_err( "Unexpected character: '%c'\n", *s->pos);
            return NULL;
    }
}

/* Does the parser position start a named-character prefix operator? */
static bool named_prefix(const ParserState* s) {
    size_t n = 0;
    const U8CharInfo* ci = src_special(s->pos, &n);
    return ci && ci->cls == U8C_PREFIX && ci->head;
}

static bool can_start_primary(const char* p) {
    unsigned char c = (unsigned char)*p;
    return isalpha(c) || isdigit(c) || c == '.' || c == '{' || c == '(' || c == '"' || c == '_' || c == '`' || c == '$' || c == '#' || c == '%' ||
           src_letter(p) > 0 || lb_group_open(p) > 0;
}

/* Scan a numeric literal starting at `p` (pointing at the first digit, or a
 * leading '.' of a `.5`-style literal) and report whether it is immediately
 * followed by a '^' (Power) operator, allowing intervening spaces/tabs.
 * Mirrors the integer/real shape parse_number accepts: digits, optional
 * '.'fraction, optional 'e[+-]?digits' exponent. */
static bool number_followed_by_power(const char* p) {
    while (isdigit((unsigned char)*p)) p++;
    /* A '.' starting `..`/`...` is not a decimal point (see parse_number). */
    if (*p == '.' && p[1] != '.') { p++; while (isdigit((unsigned char)*p)) p++; }
    if (*p == 'e' || *p == 'E') {
        const char* save = p;
        p++;
        if (*p == '+' || *p == '-') p++;
        if (isdigit((unsigned char)*p)) { while (isdigit((unsigned char)*p)) p++; }
        else p = save;
    }
    while (*p == ' ' || *p == '\t') p++;
    return *p == '^';
}

/* Decide whether a leading '-' at `s` binds as a prefix operator
 * (Times[-1, ...]) rather than folding into a numeric literal. Symbolic
 * operands (`-x`, `-(...)`, `-Sin[x]`) are always prefix. A numeric operand
 * is prefix only when the literal is immediately followed by '^', so that
 * `-3^2` parses as `-(3^2)` = Times[-1, Power[3, 2]] (matching Mathematica)
 * while plain `-3` stays a negative integer literal. */
static bool minus_is_prefix(const char* s) {
    const char* n = s + 1;
    if (!isdigit((unsigned char)*n) && !(*n == '.' && isdigit((unsigned char)n[1])))
        return true;
    return number_followed_by_power(n);
}

/* True when applying `op` to `left` would merely append to the SAME flat head
 * (Plus or Times) that `left` already has — i.e. an a+b+... or a*b*... chain
 * continuing in place. Used only during structural-selection span recording to
 * suppress partial-sum/product spans so the ladder matches Mathematica's flat
 * n-ary grouping. */
static bool flat_continuation(const Expr* left, const OperatorDef* op) {
    if (!left || left->type != EXPR_FUNCTION) return false;
    const Expr* h = left->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* lh = h->data.symbol.name;
    const char* oh = NULL;
    if (op->type == OP_MINUS) oh = "Plus";          /* a - b extends Plus */
    else if (op->type == OP_TIMES) oh = "Times";    /* explicit or implicit * */
    else if (op->head_name &&
             (strcmp(op->head_name, "Plus") == 0 || strcmp(op->head_name, "Times") == 0))
        oh = op->head_name;
    return oh && lh && strcmp(lh, oh) == 0;
}

// Pratt parser for operator precedence
static Expr* parse_expression_prec(ParserState* s, int min_prec) {
    s->chain_head = NULL;
    const char* left_start;
    Expr* left = NULL;
    bool left_bare_compare = false;
    /* Head of the flat named-character operator chain `left` was built by
     * (a \[CenterDot] b), so a further \[CenterDot] appends to it, while a
     * parenthesized (a \[CenterDot] b) \[CenterDot] c nests, as in the
     * reference language. NULL when `left` is anything else. */
    const char* left_named_flat = NULL;
    /* Head of the Plus / Times / NonCommutativeMultiply / StringJoin chain this
     * loop built `left` with, so the next operator of that kind extends it
     * (a b c is Times[a, b, c]). Anything else stays one operand, as the
     * reference reads it: (a b) c and Times[a, b] c are Times[Times[a, b], c],
     * and a prefix minus -(x + y) z is Times[Times[-1, x + y], z], so the
     * negated sum distributes before the product forms. NULL otherwise. */
    const char* left_flat_head = NULL;
    /* `left` is a `;`-chain built at this level (not a parenthesised one), so
     * a further `;` appends to it: a;b;c is flat, (a;b);c stays nested. */
    bool left_bare_compound = false;
    if (s->seed_left) {
        /* A caller-supplied left operand: straight to the infix loop. */
        left = s->seed_left;
        s->seed_left = NULL;
        left_start = s->pos;
        goto infix_loop;
    }
    skip_whitespace(s);
    if (!*s->pos) return NULL;
    /* No expression token follows a closing delimiter, a comma, or a
     * binary-only operator character. Return NULL silently so callers
     * (notably the OP_COMPOUND branch that turns `expr; // f` or a trailing
     * `expr;` into CompoundExpression[expr, Null]) can decide what to do,
     * without parse_primary printing a spurious "Unexpected character"
     * message. The chars below never start a primary or a prefix operator;
     * `+`, `-`, `.`, `!`, `?`, `'` are intentionally excluded because they
     * either begin numbers (`+5`, `-5`, `.5`) or prefix forms (`!x`, `-x`, `+x`,
     * `?sym`, `++x`, `--x`); a leading `;;` is the implicit-LHS Span
     * prefix and is also left alone. */
    if (*s->pos == ']' || *s->pos == '}' || *s->pos == ')' || *s->pos == ',' ||
        *s->pos == '*' || *s->pos == '/' || *s->pos == '^' || *s->pos == '=' ||
        (*s->pos == '<' && s->pos[1] != '|' && s->pos[1] != '<') || /* '<|' association, '<<' Get */
        *s->pos == '>' || *s->pos == '|' || *s->pos == '&' ||
        *s->pos == ':' || *s->pos == '@' ||
        (*s->pos == ';' && s->pos[1] != ';')) return NULL;

    /* Start of this (sub)expression, for structural-selection span recording. */
    left_start = s->pos;

    /* Did `left` come from an UNPARENTHESISED chainable comparison at this level?
     *
     * Chained comparisons are folded into one variadic Inequality, and the fold used to decide by
     * INSPECTING the built subtree -- which cannot tell `a >= b == c` (a chain) from
     * `(a >= b) == c` (an Equal whose left operand happens to be a comparison). Parentheses were
     * therefore ignored on the left operand: `(2.0 >= 2.) == (1.0 > 0.)` parsed as
     * Inequality[2.0, GreaterEqual, 2.0, Equal, Greater[1.0, 0.0]] and evaluated to `2.0 == True`.
     * The RIGHT operand was unaffected, since it is parsed as its own subexpression, which is why
     * only one side looked wrong.
     *
     * Plus/Times flattening a few branches below has the same shape and is harmless there, because
     * Plus is associative so `(a + b) + c` and `a + b + c` mean the same thing. Chaining changes
     * MEANING, so comparisons need the provenance tracked rather than guessed. */
    if (*s->pos == '?') {
        s->pos++;
        /* `?Name` looks up one docstring; `?Pat*` searches the symbol table.
         * A wildcard cannot occur inside a symbol, so scan the whole run and
         * decide before handing anything to parse_symbol -- that stops at the
         * `*`, leaving it for the caller to report as "Extra characters after
         * expression", which is how `?Find*` failed. The pattern is passed on
         * as a String, which Information already accepts. */
        const char* start = s->pos;
        const char* p = start;
        bool wild = false;
        while (*p && (isalnum((unsigned char)*p) || *p == '_' || *p == '$'
                      || *p == '`' || *p == '*' || *p == '@')) {
            if (*p == '*' || *p == '@') wild = true;
            p++;
        }
        if (wild && p > start) {
            size_t n = (size_t)(p - start);
            char* pat = malloc(n + 1);
            if (!pat) return NULL;
            memcpy(pat, start, n);
            pat[n] = '\0';
            s->pos = p;
            Expr* args[1] = { expr_new_string(pat) };
            free(pat);
            return expr_new_function(expr_new_symbol(SYM_Information), args, 1);
        }
        Expr* sym = parse_symbol(s);
        if (!sym) return NULL;
        Expr* args[1] = { sym };
        left = expr_new_function(expr_new_symbol(SYM_Information), args, 1);
        return left; // ?symbol is top-level only usually
    }

    if (strncmp(s->pos, "<|", 2) == 0) {
        /* Association literal. Assigning to `left` (rather than returning)
         * lets postfix forms like <|...|>[[key]] and <|...|>[args] parse. */
        left = parse_association(s);
    } else if (strncmp(s->pos, ";;", 2) == 0) {
        left = expr_new_integer(1);
    } else if (strncmp(s->pos, "++", 2) == 0) {
        s->pos += 2;
        Expr* right = parse_expression_prec(s, 7500);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol(SYM_PreIncrement), args, 1);
    } else if (strncmp(s->pos, "--", 2) == 0) {
        s->pos += 2;
        Expr* right = parse_expression_prec(s, 7500);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol(SYM_PreDecrement), args, 1);
    } else if (*s->pos == '-' && minus_is_prefix(s->pos)) {
        s->pos++;
        // Use a precedence higher than Plus (3500) and Times (4500)
        Expr* right = parse_expression_prec(s, 5100);
        if (!right) return NULL;
        Expr* minus_one = expr_new_integer(-1);
        Expr* args[2] = { minus_one, right };
        left = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    } else if (*s->pos == '+' && !isdigit((unsigned char)s->pos[1]) &&
               !(s->pos[1] == '.' && isdigit((unsigned char)s->pos[2]))) {
        /* Prefix `+` (the "++" PreIncrement case is consumed above; a signed
         * literal `+5` / `+.5` is left to parse_number). As in the Wolfram
         * Language it parses to a one-argument Plus -- `Hold[+x]` is
         * `Hold[Plus[x]]` -- which the evaluator then reduces to `x`. The
         * operand binds tighter than binary Plus, so `+a b` is
         * Plus[Times[a, b]] and `+a == b` is Equal[Plus[a], b]. */
        s->pos++;
        Expr* right = parse_expression_prec(s, 3501);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol(SYM_Plus), args, 1);
    } else if (*s->pos == '!' && s->pos[1] != '=') {
        s->pos++;
        Expr* right = parse_expression_prec(s, 3000); // Not precedence is 3000
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol(SYM_Not), args, 1);
    } else if (((unsigned char)*s->pos >= 0x80 || *s->pos == '\\') && named_prefix(s)) {
        /* A named-character prefix operator: \[Not] x, \[Del] f, \[Sqrt] x. */
        size_t n = 0;
        const U8CharInfo* ci = src_special(s->pos, &n);
        s->pos += n;
        Expr* right = parse_expression_prec(s, ci->prec);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol(ci->head), args, 1);
    } else {
        left = parse_primary(s);
    }
    
    if (!left) return NULL;

infix_loop:
    while (1) {
        /* End of `left` at its current extent, captured BEFORE the trailing
         * whitespace is skipped so a recorded span's end is tight. */
        const char* tight_end = s->pos;
        skip_whitespace(s);
        /* Line-statement mode: a complete expression at the end of a line is
         * a complete statement, whatever the next line starts with. */
        if (s->line_statements && s->bracket_depth == 0 && s->saw_newline) break;
        OperatorDef op_def = get_operator(s->pos);

        /* Lenient `:` fallback (the standard forms are handled at atom level by
         * parse_colon_suffix): Pattern when `left` is still a bare symbol
         * (e.g. `(a):p`), Optional otherwise. See PREC_PATTERN. */
        if (op_def.type == OP_COLON && left->type != EXPR_SYMBOL) {
            op_def.prec = PREC_OPTIONAL;
            op_def.head_name = "Optional";
        }
        
        // Handle implicit multiplication: if no explicit operator but another expression follows
        // (an association literal `<|...|>` and a `<<file` Get begin with `<`,
        // which is not a can_start_primary char, so admit them explicitly).
        if (op_def.type == OP_NONE && (can_start_primary(s->pos) || strncmp(s->pos, "<|", 2) == 0 ||
                                       strncmp(s->pos, "<<", 2) == 0)) {
            // A top-level newline terminates the statement rather than acting as
            // implicit multiplication: in a file `a\nb` is two statements (as in
            // Mathematica), not `a b`. Inside brackets newlines stay
            // insignificant, so `{a\nb}` and multi-line calls remain one
            // expression. An explicit operator (e.g. a trailing `+`) is handled
            // above this branch, so operator-continued lines still join.
            if (s->bracket_depth == 0 && s->saw_newline) break;
            // Implicit Times has same precedence as explicit Times (4500)
            if (4500 < min_prec) break;
            op_def.type = OP_TIMES;
            op_def.prec = 4500;
            op_def.head_name = "Times";
            op_def.len = 0; // Don't advance pos
        }

        /* Structural-selection span recording. Record `left` as one selectable
         * stage, EXCEPT when the fold about to run merely appends to the same
         * flat head (Plus/Times) `left` already is: those partial sums/products
         * are not selectable units in Mathematica's flat n-ary model (a+b+c
         * selects a..c then the whole Plus, with no a+b step). The complete flat
         * node is still recorded on the pass that finds no continuing operator. */
        {
            bool stop = (op_def.type == OP_NONE || op_def.prec < min_prec);
            if (s->spans && (stop || !flat_continuation(left, &op_def)))
                record_span(s, left_start, tight_end);
            if (stop) break;
        }

        /* Recomputed each iteration: only the branch that actually builds a bare chainable
         * comparison sets it back to true, so any other operator in between breaks the chain. */
        bool prev_bare_compare = left_bare_compare;
        left_bare_compare = false;
        const char* prev_named_flat = left_named_flat;
        left_named_flat = NULL;
        const char* prev_flat_head = left_flat_head;
        left_flat_head = NULL;
        bool prev_bare_compound = left_bare_compound;
        left_bare_compound = false;
        
        s->pos += op_def.len;

        if (op_def.type == OP_MESSAGENAME) {
            /* f::tag -> MessageName[f, "tag"]. The tag is read literally: a
             * bare identifier becomes a string, or a quoted "tag" is taken as
             * given. This mirrors Mathematica, where `f::usage` has the string
             * "usage" as its second argument. */
            skip_whitespace(s);
            Expr* tag = NULL;
            if (*s->pos == '"') {
                tag = parse_string(s);
            } else {
                char buf[256];
                size_t bi = 0;
                while (*s->pos && bi + 1 < sizeof(buf) &&
                       (isalnum((unsigned char)*s->pos) || *s->pos == '$')) {
                    buf[bi++] = *s->pos++;
                }
                buf[bi] = '\0';
                if (bi == 0) {
                    parse_err( "Expected message tag after '::'\n");
                    expr_free(left);
                    return NULL;
                }
                tag = expr_new_string(buf);
            }
            if (!tag) { expr_free(left); return NULL; }
            Expr* args[2] = { left, tag };
            left = expr_new_function(expr_new_symbol(SYM_MessageName), args, 2);
            continue;
        }

        if (op_def.type == OP_PART) {
            Expr* args[64];
            size_t count = 0;
            args[count++] = left;
            s->bracket_depth++;

            while (*s->pos && strncmp(s->pos, "]]", 2) != 0) {
                skip_whitespace(s);
                if (count > 1 && *s->pos == ',') s->pos++;
                
                Expr* arg = parse_expression_state(s);
                if (!arg) {
                    while (count--) expr_free(args[count]);
                    return NULL;
                }
                args[count++] = arg;
                skip_whitespace(s);
            }
            if (strncmp(s->pos, "]]", 2) != 0) {
                parse_err( "Expected ']]'\n");
                while (count--) expr_free(args[count]);
                return NULL;
            }
            s->pos += 2; // skip ']]'
            close_bracket(s);
            left = expr_new_function(expr_new_symbol(SYM_Part), args, count);
            continue;
        } else if (op_def.type == OP_CALL) {
            s->pos--; // parse_function expects '[' to be there
            left = parse_function(s, left);
            if (!left) return NULL;
            continue;
        } else if (op_def.type == OP_FUNCTION) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Function), args, 1);
            continue;
        } else if (op_def.type == OP_FACTORIAL) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Factorial), args, 1);
            continue;
        } else if (op_def.type == OP_FACTORIAL2) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Factorial2), args, 1);
            continue;
        } else if (op_def.type == OP_DERIVATIVE) {
            /* Collapse consecutive apostrophes: f'''  ->  Derivative[3][f]. */
            int n = 1;
            while (*s->pos == '\'') { n++; s->pos++; }
            Expr* order_args[1] = { expr_new_integer(n) };
            Expr* deriv_head = expr_new_function(expr_new_symbol(SYM_Derivative), order_args, 1);
            Expr* f_args[1] = { left };
            left = expr_new_function(deriv_head, f_args, 1);
            continue;
        } else if (op_def.type == OP_INCREMENT) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Increment), args, 1);
            continue;
        } else if (op_def.type == OP_DECREMENT) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Decrement), args, 1);
            continue;
        } else if (op_def.type == OP_NAMED_POSTFIX) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(op_def.head_name), args, 1);
            continue;
        } else if (op_def.type == OP_UNSET) {
            /* Postfix `lhs =.` -> Unset[lhs]. */
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Unset), args, 1);
            continue;
        } else if (op_def.type == OP_REPEATED) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_Repeated), args, 1);
            continue;
        } else if (op_def.type == OP_REPEATEDNULL) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol(SYM_RepeatedNull), args, 1);
            continue;
        } else if (op_def.type == OP_PUT || op_def.type == OP_PUTAPPEND) {
            /* Right side of `>>` / `>>>` is a file name, read literally
             * (parse_file_token), not as an expression. */
            Expr* filename = parse_file_token(s, op_def.head_name);
            if (!filename) {
                expr_free(left);
                return NULL;
            }
            Expr* put_args[2] = { left, filename };
            left = expr_new_function(expr_new_symbol(op_def.head_name),
                                     put_args, 2);
            continue;
        } else if (op_def.type == OP_SPAN) {
            Expr* span_args[3];
            span_args[0] = left;
            int next_prec = op_def.prec + 1; // 3201
            Expr* right = NULL;
            if (operand_omitted_at_line_end(s) || *s->pos == ']' || *s->pos == ',' || *s->pos == '}' || *s->pos == '\0' || *s->pos == ';' || strncmp(s->pos, ";;", 2) == 0) {
                right = expr_new_symbol(SYM_All);
            } else {
                right = parse_expression_prec(s, next_prec);
                if (!right) right = expr_new_symbol(SYM_All);
            }
            span_args[1] = right;

            skip_whitespace(s);
            if (!(s->line_statements && s->bracket_depth == 0 && s->saw_newline)
                && strncmp(s->pos, ";;", 2) == 0) {
                s->pos += 2;
                Expr* step = NULL;
                if (operand_omitted_at_line_end(s) || *s->pos == ']' || *s->pos == ',' || *s->pos == '}' || *s->pos == '\0' || *s->pos == ';') {
                    step = expr_new_symbol(SYM_All);
                } else {
                    step = parse_expression_prec(s, next_prec);
                    if (!step) step = expr_new_symbol(SYM_All);
                }
                span_args[2] = step;
                left = expr_new_function(expr_new_symbol(SYM_Span), span_args, 3);
            } else {
                left = expr_new_function(expr_new_symbol(SYM_Span), span_args, 2);
            }
            continue;
        }
#if UP_VALUES
        else if (op_def.type == OP_TAGSET) {
            /* f /: lhs = rhs   -> TagSet[f, lhs, rhs]
             * f /: lhs := rhs  -> TagSetDelayed[f, lhs, rhs]
             * f /: lhs =.      -> TagUnset[f, lhs]
             * `left` is the tag f. Parse the middle lhs stopping just before the
             * trailing assignment operator (min_prec 501 > the 500 of Set /
             * SetDelayed / Unset, so those are NOT consumed into lhs), then
             * dispatch on what follows. The rhs then parses at Set's own tier. */
            Expr* tag = left;
            skip_whitespace(s);
            Expr* tlhs = parse_expression_prec(s, 501);
            if (!tlhs) { expr_free(tag); return NULL; }
            skip_whitespace(s);
            if (strncmp(s->pos, "=.", 2) == 0 && !isdigit((unsigned char)s->pos[2])) {
                s->pos += 2;
                Expr* ta[2] = { tag, tlhs };
                left = expr_new_function(expr_new_symbol(SYM_TagUnset), ta, 2);
                continue;
            } else if (strncmp(s->pos, ":=", 2) == 0) {
                s->pos += 2;
                Expr* trhs = parse_expression_prec(s, 500);
                if (!trhs) { expr_free(tag); expr_free(tlhs); return NULL; }
                Expr* ta[3] = { tag, tlhs, trhs };
                left = expr_new_function(expr_new_symbol(SYM_TagSetDelayed), ta, 3);
                continue;
            } else if (*s->pos == '=' && s->pos[1] != '=') {
                s->pos += 1;
                Expr* trhs = parse_expression_prec(s, 500);
                if (!trhs) { expr_free(tag); expr_free(tlhs); return NULL; }
                Expr* ta[3] = { tag, tlhs, trhs };
                left = expr_new_function(expr_new_symbol(SYM_TagSet), ta, 3);
                continue;
            } else {
                parse_err( "Expected '=', ':=' or '=.' after /: (TagSet)\n");
                expr_free(tag); expr_free(tlhs);
                return NULL;
            }
        }
#endif

        int next_min_prec = op_def.right_assoc ? op_def.prec : op_def.prec + 1;

        /* `a;` ending a line: the omitted statement is an implicit Null and
         * the line break ends the CompoundExpression. */
        Expr* right = NULL;
        /* `a; -> c`: an operator that can only continue an expression
         * means the operand after `;` is an implicit Null and the compound is
         * that operator's left operand.  Decide BEFORE parsing an operand: a
         * failed operand parse would already have consumed the `-` of `->`
         * or `-=` as a prefix minus, leaving `> c` / `= c` behind. */
        if (op_def.type == OP_COMPOUND &&
            (operand_omitted_at_line_end(s) || starts_binary_operator(s->pos)))
            right = expr_new_symbol(SYM_Null);
        if (!right) right = parse_expression_prec(s, next_min_prec);
        const char* right_chain = right ? s->chain_head : NULL;
        if (!right) {
            skip_whitespace(s);
            if (op_def.type == OP_COMPOUND && *s->pos != '[') {
                /* `a;` -- the omitted last statement is an implicit Null.
                 * A bare '[' here (`a; [y]`) is a syntax error, not an
                 * implicit Null followed by a call on the compound. */
                right = expr_new_symbol(SYM_Null);
            } else {
                expr_free(left);
                return NULL;
            }
        }
        
        if (op_def.type == OP_MINUS) {
            /* a - b c is Plus[a, Times[-1, b, c]]: the -1 joins a bare product
             * operand, while a - (b c) keeps Times[-1, Times[b, c]]. */
            Expr* neg_right;
            if (right_chain && strcmp(right_chain, "Times") == 0 && right->type == EXPR_FUNCTION &&
                right->data.function.head->type == EXPR_SYMBOL &&
                right->data.function.head->data.symbol.name == SYM_Times) {
                size_t rn = right->data.function.arg_count;
                Expr** ta = malloc(sizeof(Expr*) * (rn + 1));
                ta[0] = expr_new_integer(-1);
                for (size_t i = 0; i < rn; i++) ta[i + 1] = expr_copy(right->data.function.args[i]);
                expr_free(right);
                neg_right = expr_new_function(expr_new_symbol(SYM_Times), ta, rn + 1);
                free(ta);
            } else {
                Expr* args_times[2] = { expr_new_integer(-1), right };
                neg_right = expr_new_function(expr_new_symbol(SYM_Times), args_times, 2);
            }
            if (prev_flat_head && strcmp(prev_flat_head, "Plus") == 0 &&
                extend_flat_head(left, "Plus", neg_right)) {
                /* left is now the extended Plus; keep it */
            } else {
                Expr* args_plus[2] = { left, neg_right };
                left = expr_new_function(expr_new_symbol(SYM_Plus), args_plus, 2);
            }
            left_flat_head = "Plus";
        } else if (op_def.type == OP_DIVIDE) {
            if (left->type == EXPR_INTEGER && right->type == EXPR_INTEGER) {
                Expr* rat_args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(SYM_Rational), rat_args, 2);
            } else {
                Expr* minus_one = expr_new_integer(-1);
                Expr* args_power[2] = { right, minus_one };
                Expr* inv_right = expr_new_function(expr_new_symbol(SYM_Power), args_power, 2);
                
                if (left->type == EXPR_INTEGER && left->data.integer == 1) {
                    expr_free(left);
                    left = inv_right;
                } else {
                    Expr* args_times[2] = { left, inv_right };
                    left = expr_new_function(expr_new_symbol(SYM_Times), args_times, 2);
                }
            }
        } else if (op_def.type == OP_APPLY) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_Apply), args, 2);
        } else if (op_def.type == OP_APPLY1) {
            Expr* level_args[1] = { expr_new_integer(1) };
            Expr* level = expr_new_function(expr_new_symbol(SYM_List), level_args, 1);
            Expr* args[3] = { left, right, level };
            left = expr_new_function(expr_new_symbol(SYM_Apply), args, 3);
        } else if (op_def.type == OP_MAP) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_Map), args, 2);
        } else if (op_def.type == OP_MAPALL) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_MapAll), args, 2);
        } else if (op_def.type == OP_REPLACEALL) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_ReplaceAll), args, 2);
        } else if (op_def.type == OP_REPLACEREPEATED) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_ReplaceRepeated), args, 2);
        } else if (op_def.type == OP_ALTERNATIVES || op_def.type == OP_STRINGEXPRESSION) {
            /* a | b | c and a ~~ b ~~ c parse n-ary, as in Mathematica, so a
             * held pattern is Alternatives[a, b, c] rather than a nested pair. */
            if (!extend_flat_head(left, op_def.head_name, right)) {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(op_def.head_name), args, 2);
            }
        } else if (op_def.type == OP_NAMED_INFIX) {
            if (op_def.flat && prev_named_flat && strcmp(prev_named_flat, op_def.head_name) == 0 &&
                extend_flat_head(left, op_def.head_name, right)) {
                /* left is the extended chain */
            } else {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(op_def.head_name), args, 2);
            }
            if (op_def.flat) left_named_flat = op_def.head_name;
        } else if (op_def.type == OP_CONDITION) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol(SYM_Condition), args, 2);
        } else if (op_def.type == OP_POSTFIX) {
            Expr* args[1] = { left };
            left = expr_new_function(right, args, 1);
        } else if (op_def.type == OP_PREFIX) {
            Expr* args[1] = { right };
            left = expr_new_function(left, args, 1);
        } else if (op_def.type == OP_COLON) {
            if (left->type == EXPR_SYMBOL) {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(SYM_Pattern), args, 2);
            } else {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(SYM_Optional), args, 2);
            }
        } else if (op_def.type == OP_COMPOUND) {
            /* Flatten `;` chains at parse time so a;b;c parses as the n-ary
             * CompoundExpression[a, b, c] rather than a nested binary tree.
             * `;` associates LEFT, so the chain is complete before a later
             * continuation operator binds: `a; b; -> c` is
             * Rule[CompoundExpression[a, b, Null], c], not a Rule on `(b;)`.
             * Only a bare chain is extended: `(a; b); c` keeps the inner
             * compound as its first statement. */
            left_bare_compound = true;
            if (!(prev_bare_compound && extend_flat_head(left, op_def.head_name, right))) {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(op_def.head_name), args, 2);
            }
        } else {
            /* Comparison chains. The parser associates left, so each operator
             * sees the already-built chain as `left`; the node is rewritten at
             * construction time into the shape Mathematica produces:
             *
             *   - SameQ, UnsameQ and Unequal are VARIADIC predicates, not
             *     pairwise chains: `a === b === c` is SameQ[a, b, c] ("all
             *     identical"), `a =!= b =!= c` is UnsameQ[a, b, c] and
             *     `a != b != c` is Unequal[a, b, c] ("all distinct"). Nesting
             *     them instead (SameQ[SameQ[a, b], c]) compared a Boolean with
             *     `c`, so `1 === 1 === 1` answered False.
             *   - Less, LessEqual, Greater, GreaterEqual and Equal fold into
             *     one Inequality[a, Less, b, LessEqual, c, ...]. Unequal joins
             *     such a chain when mixed with them (`a < b != c` is
             *     Inequality[a, Less, b, Unequal, c]).
             *
             * Only an operand built BARE at this level is extended
             * (prev_bare_compare): a parenthesised `(a === b) === c` stays
             * SameQ[SameQ[a, b], c]. */
            const char* hn = op_def.head_name;
            if (hn && prev_bare_compare && is_nary_compare_head(hn)
                && extend_flat_head(left, hn, right)) {
                left_bare_compare = true;      /* still a bare chain */
            } else if (hn && prev_bare_compare
                       && (is_chain_compare_head(hn) || strcmp(hn, "Unequal") == 0)
                       && fold_into_inequality(&left, hn, right)) {
                left_bare_compare = true;
            }
            /* Flatten repeated Plus/Times/NonCommutativeMultiply at parse time
             * so that held expressions reflect the n-ary form
             * (a+b+c -> Plus[a,b,c]; a**b**c -> NonCommutativeMultiply[a,b,c]). */
            else if (op_def.head_name && prev_flat_head &&
                strcmp(prev_flat_head, op_def.head_name) == 0 &&
                extend_flat_head(left, op_def.head_name, right)) {
                left_flat_head = prev_flat_head;      /* left was extended in place */
            } else {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(op_def.head_name), args, 2);
                /* A BARE binary chainable comparison, built here rather than parenthesised, so the
                 * next comparison at this level may fold it into a chain. This is the provenance the
                 * fold above needs and could not previously obtain, since by the time it looks the
                 * node is indistinguishable from a parenthesised one. */
                if (op_def.head_name && (is_chain_compare_head(op_def.head_name)
                                         || is_nary_compare_head(op_def.head_name)))
                    left_bare_compare = true;
                if (op_def.head_name &&
                    (strcmp(op_def.head_name, "Plus") == 0 || strcmp(op_def.head_name, "Times") == 0 ||
                     strcmp(op_def.head_name, "NonCommutativeMultiply") == 0 ||
                     strcmp(op_def.head_name, "StringJoin") == 0))
                    left_flat_head = op_def.head_name;
            }
        }
    }

    s->chain_head = left_flat_head;
    return left;
}

static Expr* parse_expression_state(ParserState* s) {
    return parse_expression_prec(s, 0);
}

// Public interface
Expr* parse_expression(const char* input) {
    ParserState state = {input, input, 0, 0, NULL, 0, 0, NULL, NULL, NULL, NULL};
    Expr* result = parse_expression_state(&state);
    
    // Check for trailing garbage
    skip_whitespace(&state);
    if (*state.pos != '\0') {
        parse_err( "Extra characters after expression: %s\n", state.pos);
        expr_free(result);
        return NULL;
    }
    
    return result;
}

/* Reading user input CREATES its symbols, as the reference reader does: after
 * `Global`v = Names["Global`*"]` the list already contains v, because the
 * statement was read (creating v) before Names ran.  The interner alone only
 * canonicalises a name; materialising it makes it a symbol Names lists.  Not
 * done while a kernel module loads (its scratch names must not appear in a
 * fresh session's Names["Global`*"]), nor by parse_expression, which C code
 * uses for internal templates. */
static void materialize_symbols(const Expr* e) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) { (void)symtab_get_def(e->data.symbol.name); return; }
    if (e->type != EXPR_FUNCTION) return;
    materialize_symbols(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        materialize_symbols(e->data.function.args[i]);
}

void parse_materialize_user_symbols(const Expr* e) {
    if (intern_system_scope_depth() == 0) materialize_symbols(e);
}

Expr* parse_next_expression(const char** input_ptr) {
    if (!input_ptr || !*input_ptr) return NULL;
    ParserState state = {*input_ptr, *input_ptr, 0, 0, NULL, 0, 0, NULL, NULL, NULL, NULL};

    /* Skip leading whitespace, comments, and empty ';' separators so that a
     * stray or doubled separator never yields a spurious empty statement. */
    for (;;) {
        skip_whitespace(&state);
        if (*state.pos == ';' && state.pos[1] != ';') { state.pos++; continue; }
        break;
    }
    if (*state.pos == '\0') { *input_ptr = state.pos; return NULL; }

    /* Parse a SINGLE top-level statement. We require a minimum precedence just
     * above CompoundExpression (prec 10) so the ';' separator is NOT consumed:
     * each ';'-terminated statement in a file is parsed and evaluated on its
     * own. This is essential for context-changing prologues (BeginPackage /
     * Begin["`Private`"]) to affect the symbols parsed in *subsequent*
     * statements -- if the whole ';'-chain were parsed as one
     * CompoundExpression up front (prec 0), every symbol would be resolved
     * under the context in force before BeginPackage had a chance to run. */
    Expr* result = parse_expression_prec(&state, 101);

    /* A statement normally ends at a top-level ';'.  But when what follows the
     * ';' can only CONTINUE an expression -- a binary or postfix operator such
     * as `/;`, `//`, `/.`, `->`, `:=`, `==`, `&&`, `&` -- the ';' is the
     * CompoundExpression operator with an implicit Null, and that compound is
     * the operator's left operand: `f[x_] := g[x];/; x > 0` is
     * Condition[CompoundExpression[f[x_] := g[x], Null], x > 0], as the
     * reference reader parses it.  Repeat for a chain of such continuations.  The
 * operator must be on the same line as the `;`. */
    while (result) {
        skip_whitespace(&state);
        if (!(*state.pos == ';' && state.pos[1] != ';')) break;
        ParserState peek = state;
        peek.pos++;
        peek.saw_newline = 0;
        skip_whitespace(&peek);
        /* Only on the same line: across a line break the `;` has ended a
         * complete statement, and a leading operator is a syntax error. */
        if (peek.saw_newline || !starts_binary_operator(peek.pos)) break;
        Expr* cargs[2] = { result, expr_new_symbol(SYM_Null) };
        state.seed_left = expr_new_function(expr_new_symbol(SYM_CompoundExpression), cargs, 2);
        state.pos = peek.pos;
        result = parse_expression_prec(&state, 101);
    }

    /* Consume one trailing top-level ';' separator so the next call begins at
     * the following statement. */
    skip_whitespace(&state);
    if (*state.pos == ';' && state.pos[1] != ';') state.pos++;

    parse_materialize_user_symbols(result);
    *input_ptr = state.pos;
    return result;
}


/* Skip blanks, comments and backslash line continuations on the current line,
 * stopping at a line terminator (LF or CR) or the end of input. A comment may
 * itself span lines; its line breaks are not statement boundaries. */
/* Length of the line break at p: 1 for '\n' / '\r', the source length of a
 * named-character line break (\[IndentingNewLine]), else 0. */
static size_t line_break_len(const char* p) {
    if (*p == '\n' || *p == '\r') return 1;
    size_t n = 0;
    const U8CharInfo* ci = src_special(p, &n);
    return (ci && ci->cls == U8C_NEWLINE) ? n : 0;
}

static const char* skip_line_blanks(const char* p) {
    size_t n = 0;
    const U8CharInfo* ci;
    for (;;) {
        if (*p == ' ' || *p == '\t' || *p == '\f' || *p == '\v') {
            p++;
        } else if ((ci = src_special(p, &n)) != NULL && ci->cls == U8C_SPACE) {
            p += n;                /* \[InvisibleSpace], \[VeryThinSpace], ... */
        } else if (p[0] == '\\' && (p[1] == '\n' || p[1] == '\r')) {
            p += (p[1] == '\r' && p[2] == '\n') ? 3 : 2;
        } else if (p[0] == '(' && p[1] == '*') {
            int depth = 1;
            p += 2;
            while (*p && depth > 0) {
                if (p[0] == '(' && p[1] == '*') { depth++; p += 2; }
                else if (p[0] == '*' && p[1] == ')') { depth--; p += 2; }
                else p++;
            }
        } else {
            return p;
        }
    }
}

int parse_statement_each(const char* input, int line_statements, const char** err_pos,
                         ParseStatementFn fn, void* ctx) {
    if (err_pos) *err_pos = NULL;
    if (!input) return 0;

    const char* p = input;
    int any = 0;
    for (;;) {
        Expr* stmt = NULL;
        p = skip_line_blanks(p);
        if (*p == '\0') break;
        size_t brk = line_break_len(p);
        if (brk && !line_statements) {
            p += brk;              /* StandardForm: a line break is whitespace */
            continue;
        }
        if (brk) {
            /* An empty (or blank/comment-only) line is a Null statement. */
            stmt = expr_new_symbol(SYM_Null);
            p += brk;
        } else {
            /* Without line_statements (StandardForm) line breaks are plain
             * whitespace: start one bracket level deep so no newline can
             * terminate the expression, and require it to fill the input. */
            ParserState st = {input, p, line_statements ? 0 : 1, 0, NULL, 0,
                              line_statements, NULL, NULL, NULL, NULL};
            stmt = parse_expression_prec(&st, 0);
            int ok = stmt != NULL;
            if (ok && line_statements) {
                /* Resume at the statement's tight end: the rest of its line
                 * must be blank, and exactly one line break is consumed so a
                 * following empty line still reads as a Null statement. */
                p = st.ws_start;
                if (p < input || p > st.pos) p = st.pos;
                p = skip_line_blanks(p);
                if ((brk = line_break_len(p)) > 0) p += brk;
                else ok = *p == '\0';
            } else if (ok) {
                p = st.pos;        /* already past trailing whitespace */
                ok = *p == '\0';
            }
            if (!ok) {
                if (err_pos) *err_pos = p;
                expr_free(stmt);
                return 0;
            }
        }
        any = 1;
        /* The statement is handed over before the next is parsed, so a caller
         * that evaluates it (ToExpression) lets Begin[...] / BeginPackage[...]
         * govern the symbols of the statements after it. */
        if (!fn(stmt, ctx)) return 2;
    }
    /* No statement at all (empty, blank or comment-only input) is Null. */
    if (!any && !fn(expr_new_symbol(SYM_Null), ctx)) return 2;
    return 1;
}

typedef struct { Expr** v; size_t n, cap; } StmtList;

static int stmt_collect(Expr* stmt, void* ctx) {
    StmtList* l = (StmtList*)ctx;
    if (l->n == l->cap) {
        size_t nc = l->cap ? 2 * l->cap : 4;
        Expr** g = realloc(l->v, nc * sizeof(Expr*));
        if (!g) { expr_free(stmt); return 0; }
        l->v = g; l->cap = nc;
    }
    l->v[l->n++] = stmt;
    return 1;
}

Expr* parse_statement_sequence(const char* input, int line_statements,
                               const char** err_pos) {
    StmtList l = { NULL, 0, 0 };
    int rc = parse_statement_each(input, line_statements, err_pos, stmt_collect, &l);
    if (rc != 1) {
        for (size_t i = 0; i < l.n; i++) expr_free(l.v[i]);
        free(l.v);
        return NULL;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), l.v, l.n);
    free(l.v);
    return list;
}

/* Parse `input` for structural selection, recording every subexpression's source
 * range without evaluating. Multi-statement input (a;b;c) is walked statement by
 * statement with absolute offsets, mirroring parse_next_expression; a hard parse
 * error advances one byte so the walk always terminates. Tolerant: the spans of
 * the parts that did parse are returned even if the tail is unparseable. */
MthSpanSink* mth_parse_spans(const char* input) {
    MthSpanSink* sink = calloc(1, sizeof *sink);
    if (!sink || !input) return sink;
    const char* p = input;
    while (*p) {
        ParserState state = {input, p, 0, 0, sink, 0, 0, NULL, NULL, NULL, NULL};
        /* Skip leading whitespace, comments, and empty ';' separators. */
        for (;;) {
            skip_whitespace(&state);
            if (*state.pos == ';' && state.pos[1] != ';') { state.pos++; continue; }
            break;
        }
        if (*state.pos == '\0') break;
        const char* before = state.pos;
        Expr* result = parse_expression_prec(&state, 101);
        if (result) expr_free(result);
        skip_whitespace(&state);
        if (*state.pos == ';' && state.pos[1] != ';') state.pos++;
        if (state.pos <= before) state.pos = before + 1; /* guarantee progress */
        p = state.pos;
    }
    return sink;
}

void mth_span_sink_free(MthSpanSink* sink) {
    if (!sink) return;
    free(sink->v);
    free(sink);
}
