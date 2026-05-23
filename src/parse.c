#include "parse.h"
#include "common.h"
#include "expr.h"
#include "context.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <gmp.h>
#include "numeric.h"

// Private parser state
typedef struct {
    const char* input;  // Original input string
    const char* pos;    // Current position
} ParserState;

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

/* ------------------- Basic Token Parsers ------------------- */

// Skips whitespace and comments
static void skip_whitespace(ParserState* s) {
    while (1) {
        if (isspace(*s->pos)) {
            s->pos++;
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
}

// Parses a symbol (x, `name`, $var) and resolves it through the context
// system, producing a canonical (possibly qualified) symbol name.
static Expr* parse_symbol(ParserState* s) {
    char buffer[256];
    size_t i = 0;

    while (isalnum(*s->pos) || *s->pos == '`' || *s->pos == '$') {
        if (i < sizeof(buffer)-1) buffer[i++] = *s->pos;
        s->pos++;
    }
    buffer[i] = '\0';

    char* resolved = context_resolve_name(buffer);
    Expr* out = expr_new_symbol(resolved ? resolved : buffer);
    free(resolved);
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
        fprintf(stderr,
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
        return expr_new_function(expr_new_symbol("Rational"), rargs, 2);
    }

    /* Fallback for any other numeric atom: leave it to the evaluator. */
    Expr* p_args[2] = { expr_new_integer(10), expr_new_integer((int64_t)exp) };
    Expr* p = expr_new_function(expr_new_symbol("Power"), p_args, 2);
    Expr* t_args[2] = { val, p };
    return expr_new_function(expr_new_symbol("Times"), t_args, 2);
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
    if (*p == '.' || *p == 'e' || *p == 'E') {
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
        if (*lp == '.') {
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
            if (mheap) free(mheap);
            if (isfinite(mval) && mval != 0.0) {
                double implied = (double)frac_digits + log10(fabs(mval));
                /* MachinePrecision threshold; derived from DBL_MANT_DIG
                 * via numeric.h so the parser tracks the platform's float
                 * representation. */
                if (implied > NUMERIC_MACHINE_PRECISION_DIGITS) {
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
            if (*mantissa_end == '.') {
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
            fprintf(stderr, "Expected exponent digits after '*^'\n");
            s->pos = save;
            expr_free(result);
            return NULL;
        }
        errno = 0;
        long exp_val = strtol(exp_start, NULL, 10);
        if (errno == ERANGE) {
            fprintf(stderr, "Exponent in '*^' suffix out of range\n");
            expr_free(result);
            return NULL;
        }
        if (exp_neg) exp_val = -exp_val;
        result = apply_pow10_scale(result, exp_val);
    }

    return result;
}

// Parses quoted strings ("text")
static Expr* parse_string(ParserState* s) {
    s->pos++;  // Skip opening quote
    char buffer[256];
    size_t i = 0;
    
    while (*s->pos && *s->pos != '"') {
        if (*s->pos == '\\') s->pos++;  // Handle escapes
        if (i < sizeof(buffer)-1) buffer[i++] = *s->pos;
        s->pos++;
    }
    
    if (*s->pos != '"') {
        fprintf(stderr, "Unterminated string\n");
        return NULL;
    }
    s->pos++;  // Skip closing quote
    
    buffer[i] = '\0';
    return expr_new_string(buffer);
}

/* ------------------- Compound Expressions ------------------- */

// Parses lists: {1,2,3}
static Expr* parse_list(ParserState* s) {
    s->pos++;  // Skip '{'
    /* Dynamic buffer: lists may have arbitrarily many elements (e.g. the
     * trigsimp rule lists exceed 60 entries). Grow as needed rather than
     * relying on a fixed stack-allocated array. */
    size_t cap = 16;
    Expr** elements = malloc(cap * sizeof(Expr*));
    size_t count = 0;

    while (*s->pos && *s->pos != '}') {
        skip_whitespace(s);
        if (*s->pos == ',') s->pos++;

        Expr* elem = parse_expression_state(s);
        if (!elem) {
            while (count--) expr_free(elements[count]);
            free(elements);
            return NULL;
        }
        if (count >= cap) {
            cap *= 2;
            elements = realloc(elements, cap * sizeof(Expr*));
        }
        elements[count++] = elem;
    }

    if (*s->pos != '}') {
        fprintf(stderr, "Expected '}' to close list\n");
        while (count--) expr_free(elements[count]);
        free(elements);
        return NULL;
    }
    s->pos++;  // Skip '}'

    // Create List[...] expression
    Expr* list = expr_new_function(expr_new_symbol("List"), elements, count);
    free(elements);
    return list;
}

// Parses functions: f[x,y]
static Expr* parse_function(ParserState* s, Expr* head) {
    s->pos++;  // Skip '['
    /* Dynamic buffer: function calls may have arbitrarily many args. */
    size_t cap = 16;
    Expr** args = malloc(cap * sizeof(Expr*));
    size_t count = 0;

    while (*s->pos && *s->pos != ']') {
        skip_whitespace(s);
        if (*s->pos == ',') s->pos++;

        Expr* arg = parse_expression_state(s);
        if (!arg) {
            while (count--) expr_free(args[count]);
            free(args);
            return NULL;
        }
        if (count >= cap) {
            cap *= 2;
            args = realloc(args, cap * sizeof(Expr*));
        }
        args[count++] = arg;
    }

    if (*s->pos != ']') {
        fprintf(stderr, "Expected ']' to close function\n");
        while (count--) expr_free(args[count]);
        free(args);
        return NULL;
    }
    s->pos++;  // Skip ']'

    if (head && head->type == EXPR_SYMBOL && head->data.symbol == SYM_Sqrt && count == 1) {
        expr_free(head);
        Expr* rat_args[2] = { expr_new_integer(1), expr_new_integer(2) };
        Expr* half = expr_new_function(expr_new_symbol("Rational"), rat_args, 2);
        Expr* pow_args[2] = { args[0], half };
        Expr* result = expr_new_function(expr_new_symbol("Power"), pow_args, 2);
        free(args);
        return result;
    }

    Expr* func = expr_new_function(head, args, count);
    free(args);
    return func;
}

/* ------------------- Main Parser ------------------- */

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
    OP_CONDITION,
    OP_ALTERNATIVES,
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
    OP_DERIVATIVE,
    OP_COMPOSITION,
    OP_PUT,
    OP_PUTAPPEND
} OperatorType;

typedef struct {
    OperatorType type;
    int prec;
    int right_assoc;
    const char* head_name;
    int len;
} OperatorDef;

static OperatorDef get_operator(const char* pos) {
    OperatorDef def = {OP_NONE, -1, 0, NULL, 0};
    
    if (strncmp(pos, ">>>", 3) == 0) {
        /* PutAppend (expr >>> file). Lower than Set (40) so `a = 5 >>> "f"`
         * parses as Set[a, PutAppend[5, "f"]] — matches Mathematica. */
        def.type = OP_PUTAPPEND; def.prec = 30; def.head_name = "PutAppend"; def.len = 3;
    } else if (strncmp(pos, "===", 3) == 0) {
        def.type = OP_SAMEQ; def.prec = 290; def.head_name = "SameQ"; def.len = 3;
    } else if (strncmp(pos, "=!=", 3) == 0) {
        def.type = OP_UNSAMEQ; def.prec = 290; def.head_name = "UnsameQ"; def.len = 3;
    } else if (strncmp(pos, "@@@", 3) == 0) {
        def.type = OP_APPLY1; def.prec = 620; def.right_assoc = 1; def.head_name = "Apply1"; def.len = 3;
    } else if (strncmp(pos, "//.", 3) == 0 && !isdigit(pos[3])) {
        def.type = OP_REPLACEREPEATED; def.prec = 110; def.right_assoc = 0; def.head_name = "ReplaceRepeated"; def.len = 3;
    } else if (strncmp(pos, "//@", 3) == 0) {
        def.type = OP_MAPALL; def.prec = 620; def.right_assoc = 1; def.head_name = "MapAll"; def.len = 3;
    } else if (strncmp(pos, "//", 2) == 0) {
        def.type = OP_POSTFIX; def.prec = 70; def.head_name = "Postfix"; def.len = 2;
    } else if (strncmp(pos, "/.", 2) == 0 && !isdigit(pos[2])) {
        def.type = OP_REPLACEALL; def.prec = 110; def.right_assoc = 0; def.head_name = "ReplaceAll"; def.len = 2;
    } else if (strncmp(pos, "@@", 2) == 0) {
        def.type = OP_APPLY; def.prec = 620; def.right_assoc = 1; def.head_name = "Apply"; def.len = 2;
    } else if (strncmp(pos, "@*", 2) == 0) {
        def.type = OP_COMPOSITION; def.prec = 625; def.right_assoc = 0; def.head_name = "Composition"; def.len = 2;
    } else if (strncmp(pos, "/@", 2) == 0) {
        def.type = OP_MAP; def.prec = 620; def.right_assoc = 1; def.head_name = "Map"; def.len = 2;
    } else if (strncmp(pos, ":>", 2) == 0) {
        def.type = OP_RULEDELAYED; def.prec = 120; def.right_assoc = 1; def.head_name = "RuleDelayed"; def.len = 2;
    } else if (strncmp(pos, "->", 2) == 0) {
        def.type = OP_RULE; def.prec = 120; def.right_assoc = 1; def.head_name = "Rule"; def.len = 2;
    } else if (strncmp(pos, "/;", 2) == 0) {
        def.type = OP_CONDITION; def.prec = 130; def.right_assoc = 1; def.head_name = "Condition"; def.len = 2;
    } else if (strncmp(pos, ";;", 2) == 0) {
        def.type = OP_SPAN; def.prec = 290; def.right_assoc = 0; def.head_name = "Span"; def.len = 2;
    } else if (*pos == ';') {
        def.type = OP_COMPOUND; def.prec = 10; def.right_assoc = 1; def.head_name = "CompoundExpression"; def.len = 1;
    } else if (strncmp(pos, "||", 2) == 0) {
        def.type = OP_OR; def.prec = 215; def.head_name = "Or"; def.len = 2;
    } else if (*pos == '|') {
        def.type = OP_ALTERNATIVES; def.prec = 160; def.head_name = "Alternatives"; def.len = 1;
    } else if (strncmp(pos, "==", 2) == 0) {
        def.type = OP_EQUAL; def.prec = 290; def.head_name = "Equal"; def.len = 2;
    } else if (strncmp(pos, "&&", 2) == 0) {
        def.type = OP_AND; def.prec = 215; def.head_name = "And"; def.len = 2;
    } else if (*pos == '&') {
        def.type = OP_FUNCTION; def.prec = 90; def.head_name = "Function"; def.len = 1;
    } else if (strncmp(pos, "!=", 2) == 0) {
        def.type = OP_UNEQUAL; def.prec = 290; def.head_name = "Unequal"; def.len = 2;
    } else if (strncmp(pos, "<>", 2) == 0) {
        def.type = OP_STRINGJOIN; def.prec = 600; def.head_name = "StringJoin"; def.len = 2;
    } else if (strncmp(pos, "<=", 2) == 0) {
        def.type = OP_LESSEQUAL; def.prec = 290; def.head_name = "LessEqual"; def.len = 2;
    } else if (strncmp(pos, ">=", 2) == 0) {
        def.type = OP_GREATEREQUAL; def.prec = 290; def.head_name = "GreaterEqual"; def.len = 2;
    } else if (strncmp(pos, ">>", 2) == 0) {
        /* Put (expr >> file). Same precedence as PutAppend; right side is
         * a filename token (bare identifier or quoted string), parsed
         * specially in the operator loop below. */
        def.type = OP_PUT; def.prec = 30; def.head_name = "Put"; def.len = 2;
    } else if (*pos == '<') {
        def.type = OP_LESS; def.prec = 290; def.head_name = "Less"; def.len = 1;
    } else if (*pos == '>') {
        def.type = OP_GREATER; def.prec = 290; def.head_name = "Greater"; def.len = 1;
    } else if (strncmp(pos, ":=", 2) == 0) {
        def.type = OP_SETDELAYED; def.prec = 40; def.right_assoc = 1; def.head_name = "SetDelayed"; def.len = 2;
    } else if (strncmp(pos, "[[", 2) == 0) {
        def.type = OP_PART; def.prec = 1100; def.head_name = "Part"; def.len = 2;
    } else if (*pos == ':') {
        def.type = OP_COLON; def.prec = 140; def.right_assoc = 1; def.head_name = "Optional"; def.len = 1;
    } else if (*pos == '=') {
        def.type = OP_SET; def.prec = 40; def.right_assoc = 1; def.head_name = "Set"; def.len = 1;
    } else if (strncmp(pos, "++", 2) == 0) {
        def.type = OP_INCREMENT; def.prec = 660; def.head_name = "Increment"; def.len = 2;
    } else if (strncmp(pos, "--", 2) == 0) {
        def.type = OP_DECREMENT; def.prec = 660; def.head_name = "Decrement"; def.len = 2;
    } else if (strncmp(pos, "+=", 2) == 0) {
        def.type = OP_ADDTO; def.prec = 40; def.right_assoc = 1; def.head_name = "AddTo"; def.len = 2;
    } else if (strncmp(pos, "-=", 2) == 0) {
        def.type = OP_SUBTRACTFROM; def.prec = 40; def.right_assoc = 1; def.head_name = "SubtractFrom"; def.len = 2;
    } else if (*pos == '+') {
        def.type = OP_PLUS; def.prec = 310; def.head_name = "Plus"; def.len = 1;
    } else if (*pos == '-') {
        def.type = OP_MINUS; def.prec = 310; def.head_name = "Plus"; def.len = 1;
    } else if (*pos == '*') {
        def.type = OP_TIMES; def.prec = 400; def.head_name = "Times"; def.len = 1;
    } else if (*pos == '/') {
        def.type = OP_DIVIDE; def.prec = 470; def.head_name = "Divide"; def.len = 1;
    } else if (strncmp(pos, "...", 3) == 0) {
        def.type = OP_REPEATEDNULL; def.prec = 170; def.head_name = "RepeatedNull"; def.len = 3;
    } else if (strncmp(pos, "..", 2) == 0) {
        def.type = OP_REPEATED; def.prec = 170; def.head_name = "Repeated"; def.len = 2;
    } else if (*pos == '.' && !isdigit(pos[1])) {
        def.type = OP_DOT; def.prec = 490; def.head_name = "Dot"; def.len = 1;
    } else if (*pos == '^') {
        def.type = OP_POWER; def.prec = 590; def.right_assoc = 1; def.head_name = "Power"; def.len = 1;
    } else if (*pos == '?') {
        def.type = OP_PATTERNTEST; def.prec = 680; def.right_assoc = 0; def.head_name = "PatternTest"; def.len = 1;
    } else if (*pos == '@') {
        def.type = OP_PREFIX; def.prec = 620; def.right_assoc = 1; def.head_name = "Prefix"; def.len = 1;
    } else if (*pos == '[') {
        def.type = OP_CALL; def.prec = 1000; def.len = 1;
    } else if (*pos == '!' && pos[1] == '!') {
        def.type = OP_FACTORIAL2; def.prec = 710; def.head_name = "Factorial2"; def.len = 2;
    } else if (*pos == '!' && pos[1] != '=') {
        def.type = OP_FACTORIAL; def.prec = 710; def.head_name = "Factorial"; def.len = 1;
    } else if (*pos == '\'') {
        def.type = OP_DERIVATIVE; def.prec = 670; def.head_name = "Derivative"; def.len = 1;
    }

    return def;
}

// Parses primary expressions (atoms, lists, functions, parens)
static Expr* parse_primary(ParserState* s) {
    skip_whitespace(s);
    if (!*s->pos) return NULL;
    
    switch (*s->pos) {
        case '"': return parse_string(s);
        case '{': return parse_list(s);
        case '[': return parse_function(s, NULL);  // Anonymous function
        case '(': {
            s->pos++;
            Expr* inner = parse_expression_prec(s, 0);
            skip_whitespace(s);
            if (*s->pos == ')') {
                s->pos++;
            } else {
                fprintf(stderr, "Expected ')'\n");
            }
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
                return expr_new_function(expr_new_symbol("Out"), out_args, 1);
            } else {
                Expr* out_args[1] = { expr_new_integer(-count) };
                return expr_new_function(expr_new_symbol("Out"), out_args, 1);
            }
        }
        
        case '#': {
            s->pos++;
            int count = 1;
            while (*s->pos == '#') {
                count++;
                s->pos++;
            }
            int64_t n = 1; // Default
            if (isdigit(*s->pos)) {
                char* end;
                n = strtoll(s->pos, &end, 10);
                s->pos = end;
            }
            if (count == 1) {
                Expr* args[1] = { expr_new_integer(n) };
                return expr_new_function(expr_new_symbol("Slot"), args, 1);
            } else {
                Expr* args[1] = { expr_new_integer(n) };
                return expr_new_function(expr_new_symbol("SlotSequence"), args, 1);
            }
        }
        
        default:
            if (isalpha(*s->pos) || *s->pos == '`' || *s->pos == '$' || *s->pos == '_') {
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
                    if (isalpha(*s->pos) || *s->pos == '`' || *s->pos == '$') {
                        blank_head = parse_symbol(s);
                    }
                    
                    int n_args = blank_head ? 1 : 0;
                    Expr* args_b[1];
                    if (blank_head) args_b[0] = blank_head;
                    
                    if (underscores == 1) {
                        blank = expr_new_function(expr_new_symbol("Blank"), blank_head ? args_b : NULL, n_args);
                    } else if (underscores == 2) {
                        blank = expr_new_function(expr_new_symbol("BlankSequence"), blank_head ? args_b : NULL, n_args);
                    } else {
                        blank = expr_new_function(expr_new_symbol("BlankNullSequence"), blank_head ? args_b : NULL, n_args);
                    }
                    
                    if (head) {
                        Expr* args[2] = { head, blank };
                        head = expr_new_function(expr_new_symbol("Pattern"), args, 2);
                    } else {
                        head = blank;
                    }
                    
                    if (*s->pos == '.') {
                        s->pos++;
                        Expr* opt_args[1] = { head };
                        head = expr_new_function(expr_new_symbol("Optional"), opt_args, 1);
                    }
                }
                
                return head;

            }
            if (isdigit(*s->pos) || ((*s->pos == '-' || *s->pos == '+') && (isdigit(s->pos[1]) || (s->pos[1] == '.' && isdigit(s->pos[2])))) || (*s->pos == '.' && isdigit(s->pos[1]))) {
                return parse_number(s);
            }
            fprintf(stderr, "Unexpected character: '%c'\n", *s->pos);
            return NULL;
    }
}

static bool can_start_primary(char c) {
    return isalpha(c) || isdigit(c) || c == '.' || c == '{' || c == '(' || c == '"' || c == '_' || c == '`' || c == '$' || c == '#' || c == '%';
}

// Pratt parser for operator precedence
static Expr* parse_expression_prec(ParserState* s, int min_prec) {
    skip_whitespace(s);
    if (!*s->pos) return NULL;
    /* No expression token follows a closing delimiter, a comma, or a
     * binary-only operator character. Return NULL silently so callers
     * (notably the OP_COMPOUND branch that turns `expr; // f` or a trailing
     * `expr;` into CompoundExpression[expr, Null]) can decide what to do,
     * without parse_primary printing a spurious "Unexpected character"
     * message. The chars below never start a primary or a prefix operator;
     * `+`, `-`, `.`, `!`, `?`, `'` are intentionally excluded because they
     * either begin numbers (`+5`, `-5`, `.5`) or prefix forms (`!x`, `-x`,
     * `?sym`, `++x`, `--x`); a leading `;;` is the implicit-LHS Span
     * prefix and is also left alone. */
    if (*s->pos == ']' || *s->pos == '}' || *s->pos == ')' || *s->pos == ',' ||
        *s->pos == '*' || *s->pos == '/' || *s->pos == '^' || *s->pos == '=' ||
        *s->pos == '<' || *s->pos == '>' || *s->pos == '|' || *s->pos == '&' ||
        *s->pos == ':' || *s->pos == '@' ||
        (*s->pos == ';' && s->pos[1] != ';')) return NULL;

    Expr* left = NULL;
    
    if (*s->pos == '?') {
        s->pos++;
        Expr* sym = parse_symbol(s);
        if (!sym) return NULL;
        Expr* args[1] = { sym };
        left = expr_new_function(expr_new_symbol("Information"), args, 1);
        return left; // ?symbol is top-level only usually
    }

    if (strncmp(s->pos, ";;", 2) == 0) {
        left = expr_new_integer(1);
    } else if (strncmp(s->pos, "++", 2) == 0) {
        s->pos += 2;
        Expr* right = parse_expression_prec(s, 660);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol("PreIncrement"), args, 1);
    } else if (strncmp(s->pos, "--", 2) == 0) {
        s->pos += 2;
        Expr* right = parse_expression_prec(s, 660);
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol("PreDecrement"), args, 1);
    } else if (*s->pos == '-' && !isdigit(s->pos[1]) && !(s->pos[1] == '.' && isdigit(s->pos[2]))) {
        s->pos++;
        // Use a precedence higher than Plus (310) and Times (400)
        Expr* right = parse_expression_prec(s, 480);
        if (!right) return NULL;
        Expr* minus_one = expr_new_integer(-1);
        Expr* args[2] = { minus_one, right };
        left = expr_new_function(expr_new_symbol("Times"), args, 2);
    } else if (*s->pos == '!' && s->pos[1] != '=') {
        s->pos++;
        Expr* right = parse_expression_prec(s, 230); // Not precedence is 230
        if (!right) return NULL;
        Expr* args[1] = { right };
        left = expr_new_function(expr_new_symbol("Not"), args, 1);
    } else {
        left = parse_primary(s);
    }
    
    if (!left) return NULL;
    
    while (1) {
        skip_whitespace(s);
        OperatorDef op_def = get_operator(s->pos);
        
        // Handle implicit multiplication: if no explicit operator but another expression follows
        if (op_def.type == OP_NONE && can_start_primary(*s->pos)) {
            // Implicit Times has same precedence as explicit Times (400)
            if (400 < min_prec) break;
            op_def.type = OP_TIMES;
            op_def.prec = 400;
            op_def.head_name = "Times";
            op_def.len = 0; // Don't advance pos
        }

        if (op_def.type == OP_NONE || op_def.prec < min_prec) break;
        
        s->pos += op_def.len;
        
        if (op_def.type == OP_PART) {
            Expr* args[64];
            size_t count = 0;
            args[count++] = left;
            
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
                fprintf(stderr, "Expected ']]'\n");
                while (count--) expr_free(args[count]);
                return NULL;
            }
            s->pos += 2; // skip ']]'
            left = expr_new_function(expr_new_symbol("Part"), args, count);
            continue;
        } else if (op_def.type == OP_CALL) {
            s->pos--; // parse_function expects '[' to be there
            left = parse_function(s, left);
            if (!left) return NULL;
            continue;
        } else if (op_def.type == OP_FUNCTION) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Function"), args, 1);
            continue;
        } else if (op_def.type == OP_FACTORIAL) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Factorial"), args, 1);
            continue;
        } else if (op_def.type == OP_FACTORIAL2) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Factorial2"), args, 1);
            continue;
        } else if (op_def.type == OP_DERIVATIVE) {
            /* Collapse consecutive apostrophes: f'''  ->  Derivative[3][f]. */
            int n = 1;
            while (*s->pos == '\'') { n++; s->pos++; }
            Expr* order_args[1] = { expr_new_integer(n) };
            Expr* deriv_head = expr_new_function(expr_new_symbol("Derivative"), order_args, 1);
            Expr* f_args[1] = { left };
            left = expr_new_function(deriv_head, f_args, 1);
            continue;
        } else if (op_def.type == OP_INCREMENT) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Increment"), args, 1);
            continue;
        } else if (op_def.type == OP_DECREMENT) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Decrement"), args, 1);
            continue;
        } else if (op_def.type == OP_REPEATED) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("Repeated"), args, 1);
            continue;
        } else if (op_def.type == OP_REPEATEDNULL) {
            Expr* args[1] = { left };
            left = expr_new_function(expr_new_symbol("RepeatedNull"), args, 1);
            continue;
        } else if (op_def.type == OP_PUT || op_def.type == OP_PUTAPPEND) {
            /* Right side of `>>` / `>>>` is a filename — either a quoted
             * string or a bare token. Mathematica treats the right side
             * literally rather than as an expression. */
            skip_whitespace(s);
            Expr* filename = NULL;
            if (*s->pos == '"') {
                filename = parse_string(s);
            } else {
                /* Scan identifier-ish chars plus a few path-safe extras. */
                char buf[1024];
                size_t bi = 0;
                while (*s->pos && bi + 1 < sizeof(buf) &&
                       (isalnum((unsigned char)*s->pos) ||
                        *s->pos == '_' || *s->pos == '.' ||
                        *s->pos == '/' || *s->pos == '\\' ||
                        *s->pos == '-' || *s->pos == '~' ||
                        *s->pos == '$')) {
                    buf[bi++] = *s->pos++;
                }
                buf[bi] = '\0';
                if (bi == 0) {
                    fprintf(stderr, "Expected filename after %s\n",
                            op_def.head_name);
                    expr_free(left);
                    return NULL;
                }
                filename = expr_new_string(buf);
            }
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
            int next_prec = op_def.prec + 1; // 291
            skip_whitespace(s);
            Expr* right = NULL;
            if (*s->pos == ']' || *s->pos == ',' || *s->pos == '}' || *s->pos == '\0' || *s->pos == ';' || strncmp(s->pos, ";;", 2) == 0) {
                right = expr_new_symbol("All");
            } else {
                right = parse_expression_prec(s, next_prec);
                if (!right) right = expr_new_symbol("All");
            }
            span_args[1] = right;
            
            skip_whitespace(s);
            if (strncmp(s->pos, ";;", 2) == 0) {
                s->pos += 2;
                skip_whitespace(s);
                Expr* step = NULL;
                if (*s->pos == ']' || *s->pos == ',' || *s->pos == '}' || *s->pos == '\0' || *s->pos == ';') {
                    step = expr_new_symbol("All");
                } else {
                    step = parse_expression_prec(s, next_prec);
                    if (!step) step = expr_new_symbol("All");
                }
                span_args[2] = step;
                left = expr_new_function(expr_new_symbol("Span"), span_args, 3);
            } else {
                left = expr_new_function(expr_new_symbol("Span"), span_args, 2);
            }
            continue;
        }
        
        int next_min_prec = op_def.right_assoc ? op_def.prec : op_def.prec + 1;
        
        Expr* right = parse_expression_prec(s, next_min_prec);
        if (!right) {
            if (op_def.type == OP_COMPOUND) {
                right = expr_new_symbol("Null");
            } else {
                expr_free(left);
                return NULL;
            }
        }
        
        if (op_def.type == OP_MINUS) {
            Expr* minus_one = expr_new_integer(-1);
            Expr* args_times[2] = { minus_one, right };
            Expr* neg_right = expr_new_function(expr_new_symbol("Times"), args_times, 2);
            if (extend_flat_head(left, "Plus", neg_right)) {
                /* left is now the extended Plus; keep it */
            } else {
                Expr* args_plus[2] = { left, neg_right };
                left = expr_new_function(expr_new_symbol("Plus"), args_plus, 2);
            }
        } else if (op_def.type == OP_DIVIDE) {
            if (left->type == EXPR_INTEGER && right->type == EXPR_INTEGER) {
                Expr* rat_args[2] = { left, right };
                left = expr_new_function(expr_new_symbol("Rational"), rat_args, 2);
            } else {
                Expr* minus_one = expr_new_integer(-1);
                Expr* args_power[2] = { right, minus_one };
                Expr* inv_right = expr_new_function(expr_new_symbol("Power"), args_power, 2);
                
                if (left->type == EXPR_INTEGER && left->data.integer == 1) {
                    expr_free(left);
                    left = inv_right;
                } else {
                    Expr* args_times[2] = { left, inv_right };
                    left = expr_new_function(expr_new_symbol("Times"), args_times, 2);
                }
            }
        } else if (op_def.type == OP_APPLY) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("Apply"), args, 2);
        } else if (op_def.type == OP_APPLY1) {
            Expr* level_args[1] = { expr_new_integer(1) };
            Expr* level = expr_new_function(expr_new_symbol("List"), level_args, 1);
            Expr* args[3] = { left, right, level };
            left = expr_new_function(expr_new_symbol("Apply"), args, 3);
        } else if (op_def.type == OP_MAP) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("Map"), args, 2);
        } else if (op_def.type == OP_MAPALL) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("MapAll"), args, 2);
        } else if (op_def.type == OP_REPLACEALL) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("ReplaceAll"), args, 2);
        } else if (op_def.type == OP_REPLACEREPEATED) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("ReplaceRepeated"), args, 2);
        } else if (op_def.type == OP_ALTERNATIVES) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("Alternatives"), args, 2);
        } else if (op_def.type == OP_CONDITION) {
            Expr* args[2] = { left, right };
            left = expr_new_function(expr_new_symbol("Condition"), args, 2);
        } else if (op_def.type == OP_POSTFIX) {
            Expr* args[1] = { left };
            left = expr_new_function(right, args, 1);
        } else if (op_def.type == OP_PREFIX) {
            Expr* args[1] = { right };
            left = expr_new_function(left, args, 1);
        } else if (op_def.type == OP_COLON) {
            if (left->type == EXPR_SYMBOL) {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol("Pattern"), args, 2);
            } else {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol("Optional"), args, 2);
            }
        } else {
            /* Flatten repeated Plus/Times at parse time so that held
             * expressions reflect the n-ary form (a+b+c -> Plus[a,b,c]). */
            if (op_def.head_name &&
                (strcmp(op_def.head_name, "Plus") == 0 ||
                 strcmp(op_def.head_name, "Times") == 0) &&
                extend_flat_head(left, op_def.head_name, right)) {
                /* left was extended in place */
            } else {
                Expr* args[2] = { left, right };
                left = expr_new_function(expr_new_symbol(op_def.head_name), args, 2);
            }
        }
    }

    return left;
}

static Expr* parse_expression_state(ParserState* s) {
    return parse_expression_prec(s, 0);
}

// Public interface
Expr* parse_expression(const char* input) {
    ParserState state = {input, input};
    Expr* result = parse_expression_state(&state);
    
    // Check for trailing garbage
    skip_whitespace(&state);
    if (*state.pos != '\0') {
        fprintf(stderr, "Extra characters after expression: %s\n", state.pos);
        expr_free(result);
        return NULL;
    }
    
    return result;
}

Expr* parse_next_expression(const char** input_ptr) {
    if (!input_ptr || !*input_ptr) return NULL;
    ParserState state = {*input_ptr, *input_ptr};
    skip_whitespace(&state);
    if (*state.pos == '\0') return NULL;
    
    Expr* result = parse_expression_state(&state);
    
    skip_whitespace(&state);
    *input_ptr = state.pos;
    
    return result;
}

