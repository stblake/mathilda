/* Mathilda — NumberForm and a minimal Row (numeric-display formatting).
 *
 * NumberForm[expr], NumberForm[expr, n], NumberForm[expr, {n, f}] is a PRINT
 * WRAPPER: the head survives in the expression tree (it is inert — the builtin
 * returns NULL) and only changes how the wrapped expression is *rendered*. The
 * printer (print.c) recognises SYM_NumberForm in print_standard, builds a
 * NumberFormCtx from the spec and options, installs it as the active context,
 * and recurses; every numeric leaf then routes through numberform_render_number
 * instead of the default real/integer formatter.
 *
 * Row[{a, b, c}] concatenates the OutputForm rendering of its list elements,
 * needed by NumberForm's NumberFormat option (which assembles the mantissa,
 * base, and exponent). It, too, is an inert head handled in print.c.
 *
 * See docs/spec/builtins for the surface behaviour and every option.
 */
#ifndef NUMBERFORM_H
#define NUMBERFORM_H

#include "expr.h"
#include <stdbool.h>

/* Active formatting context for one NumberForm[...] wrapper. Built by
 * numberform_build_ctx from the wrapper's spec (n / {n,f}) and options, then
 * consulted by numberform_render_number for each numeric leaf. All string
 * fields are BORROWED — they point into static defaults, the symbol's
 * registered default-options list, or the wrapper's own argument exprs, every
 * one of which outlives the print. The ctx is a plain stack value with no
 * cleanup. */
typedef struct NumberFormCtx {
    /* precision spec */
    long n;        /* requested significant digits; -1 ⇒ default_print_precision */
    long f;        /* fractional-digit count for {n,f}; -1 if absent */
    bool have_f;

    /* options (resolved) */
    const char* sign_neg;   /* NumberSigns[[1]] */
    const char* sign_pos;   /* NumberSigns[[2]] */
    const char* point;      /* NumberPoint */
    const char* mult;       /* NumberMultiplier */
    long block_left;        /* DigitBlock int-side (0 = Infinity / no grouping) */
    long block_right;       /* DigitBlock frac-side */
    const char* sep_int;    /* NumberSeparator[[1]] */
    const char* sep_frac;   /* NumberSeparator[[2]] */
    const char* pad_left;   /* NumberPadding[[1]] */
    const char* pad_right;  /* NumberPadding[[2]] */
    bool sign_padding;      /* SignPadding */
    long exp_step;          /* ExponentStep (>= 1) */
    const Expr* exp_func;   /* ExponentFunction; NULL ⇒ Automatic */
    const Expr* num_format; /* NumberFormat;      NULL ⇒ Automatic */
    long sci_lo, sci_hi;    /* ScientificNotationThreshold {lo, hi} */
    long default_print_precision;

    /* padding: true iff either pad string is non-empty (measure + align run
     * only then). */
    bool has_padding;

    /* measure-pass maxima across every numeric leaf of the wrapped expr */
    int m_maxsignn;
    int m_maxintdig;
    int m_maxfrac;
    int m_maxscin;
    bool m_anypoint;

    /* right-justify field metrics, derived from the maxima after the measure
     * pass (valid only when has_padding). */
    int field_sign;
    int field_int;
    int field_point;
    int field_frac;
} NumberFormCtx;

/* Build a context from a NumberForm[...] call node `res`
 * (NumberForm[expr, spec?, opts...]). Fills *ctx and, when padding is active,
 * runs the measure pass over the wrapped expression. Returns true when the
 * call has at least one argument (the printer then installs the ctx and
 * recurses into args[0]); false for the degenerate 0-argument NumberForm[],
 * which the printer renders generically. */
bool numberform_build_ctx(const Expr* res, NumberFormCtx* ctx);

/* Render one numeric atom (Integer / BigInt / Real / MPFR) to stdout under
 * `ctx`. Returns true when handled; false for a value the engine declines
 * (non-finite real / MPFR), so the caller falls back to the default printer. */
bool numberform_render_number(const Expr* e, const NumberFormCtx* ctx);

/* Register NumberForm and Row (builtins, attributes, docstrings). */
void numberform_init(void);

#endif /* NUMBERFORM_H */
