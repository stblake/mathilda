/* groebnerbasis.c
 *
 * `GroebnerBasis[polys, vars]` and `GroebnerBasis[polys, mainVars,
 * elimVars]` -- the user-facing entry point that drives the Buchberger
 * core in groebner.c.  See SPEC.md and docs/spec/builtins/arithmetic-
 * and-algebra.md for the surface this implements.
 *
 * Supported options (others emit `GroebnerBasis::nimpl`):
 *   MonomialOrder      -> Lexicographic            (default)
 *                       | DegreeReverseLexicographic
 *                       | EliminationOrder         (forced by the 3-arg form)
 *   CoefficientDomain  -> Rationals | Automatic    (default)
 *   Method             -> "Buchberger" | Automatic (default)
 */

#include "groebnerbasis.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "groebner.h"
#include "internal.h"
#include "rationalize.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/*  Diagnostic helpers                                                 */
/* ------------------------------------------------------------------ */

static void warn_once(const char* tag, const char* msg) {
    fprintf(stderr, "GroebnerBasis::%s: %s\n", tag, msg);
}

/* ------------------------------------------------------------------ */
/*  Argument shape detection                                           */
/* ------------------------------------------------------------------ */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
         || e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

static bool is_equal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Equal;
}

/* ------------------------------------------------------------------ */
/*  Option extraction                                                  */
/* ------------------------------------------------------------------ */

/* Selectable Gröbner engine (Method -> ...). */
typedef enum {
    GB_METHOD_BUCHBERGER = 0,       /* "Buchberger" / Automatic (default) */
    GB_METHOD_WALK                  /* "GroebnerWalk" */
} GBMethod;

/* Coefficient domain (CoefficientDomain -> ...). */
typedef enum {
    GB_DOM_RATIONALS = 0,           /* Rationals / Automatic (default) */
    GB_DOM_INTEGERS,                /* Integers: strong GB over Z */
    GB_DOM_RATIONAL_FUNCTIONS,      /* RationalFunctions: GB over Q(params) */
    GB_DOM_POLYNOMIALS,             /* Polynomials[x,...]: x is a coeff-ring var */
    GB_DOM_INEXACT                  /* InexactNumbers[prec]: approximate output */
} GBCoeffDomain;

/* Bundle of user-visible options recognised by `GroebnerBasis[]`.  The
 * fields whose `*_set` companion is false mean "use the engine default";
 * the only one that is consulted directly is `order` (see below). */
typedef struct {
    GBOrder order;                  /* MonomialOrder -> symbolic order */
    Expr*   matrix_order;           /* MonomialOrder -> {{...},...}: borrowed
                                       weight-matrix value, NULL otherwise.
                                       Validated at dispatch (needs n_main). */
    GBMethod method;                /* Method -> ... */
    bool    sort_desc;              /* Sort -> True reverses the default
                                       LM-ascending output ordering */
    bool    parameter_vars_given;   /* ParameterVariables option set */
    Expr*   parameter_vars;         /* borrowed; List[...] or symbol or NULL */
    GBCoeffDomain domain;           /* CoefficientDomain -> ... */
    int     inexact_prec;           /* decimal digits for GB_DOM_INEXACT */
    Expr*   poly_domain_vars;       /* borrowed; args of Polynomials[...] */
} GBOptions;

/* Forward declaration (defined with the weight-matrix helpers below). */
static bool expr_to_i64(const Expr* e, int64_t* out);

/* In-out: `*n_pos` starts at the full argument count; the function
 * decrements it for each trailing Rule[] arg.  Recognised options are
 * folded into `*opt`; unknown options emit a `nimpl` diagnostic but
 * never abort evaluation.  Options whose value we cannot honour
 * (e.g. CoefficientDomain -> Reals, Modulus -> 7) fall back to the
 * default behaviour with a once-per-call note. */
static void extract_options(const Expr* res, size_t* n_pos,
                            GBOptions* opt) {
    size_t argc = res->data.function.arg_count;
    size_t cut = argc;
    while (cut > 0 && is_rule(res->data.function.args[cut - 1])) cut--;
    *n_pos = cut;

    /* Defaults. */
    opt->order = GB_ORDER_LEX;
    opt->matrix_order = NULL;
    opt->method = GB_METHOD_BUCHBERGER;
    opt->sort_desc = false;
    opt->parameter_vars_given = false;
    opt->parameter_vars = NULL;
    opt->domain = GB_DOM_RATIONALS;
    opt->inexact_prec = 16;
    opt->poly_domain_vars = NULL;

    for (size_t i = cut; i < argc; i++) {
        Expr* rule = res->data.function.args[i];
        if (rule->data.function.arg_count != 2) continue;
        Expr* key = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (key->type != EXPR_SYMBOL) continue;

        if (key->data.symbol.name == SYM_MonomialOrder) {
            if (val->type == EXPR_SYMBOL) {
                if (val->data.symbol.name == SYM_Lexicographic
                 || val->data.symbol.name == SYM_Automatic) {
                    opt->order = GB_ORDER_LEX;
                } else if (val->data.symbol.name == SYM_DegreeReverseLexicographic) {
                    opt->order = GB_ORDER_GREVLEX;
                } else if (val->data.symbol.name == SYM_EliminationOrder) {
                    opt->order = GB_ORDER_ELIM;
                } else {
                    warn_once("nimpl", "unsupported MonomialOrder value; "
                                       "falling back to Lexicographic");
                    opt->order = GB_ORDER_LEX;
                }
            } else {
                /* A non-symbol value is a candidate weight matrix.  It is
                 * validated at dispatch (the column count must match the
                 * number of main variables, which is not known here). */
                opt->matrix_order = val;
            }
        } else if (key->data.symbol.name == SYM_CoefficientDomain) {
            if (val->type == EXPR_SYMBOL) {
                if (val->data.symbol.name == SYM_Rationals
                 || val->data.symbol.name == SYM_Automatic) {
                    opt->domain = GB_DOM_RATIONALS;
                } else if (val->data.symbol.name == SYM_Integers) {
                    opt->domain = GB_DOM_INTEGERS;
                } else if (val->data.symbol.name == SYM_RationalFunctions) {
                    opt->domain = GB_DOM_RATIONAL_FUNCTIONS;
                } else if (val->data.symbol.name == SYM_InexactNumbers) {
                    /* Bare InexactNumbers -> machine precision (~16 digits). */
                    opt->domain = GB_DOM_INEXACT;
                    opt->inexact_prec = 16;
                } else {
                    warn_once("nimpl", "unsupported CoefficientDomain value; "
                                       "falling back to Rationals");
                    opt->domain = GB_DOM_RATIONALS;
                }
            } else if (val->type == EXPR_FUNCTION
                       && val->data.function.head
                       && val->data.function.head->type == EXPR_SYMBOL) {
                const char* h = val->data.function.head->data.symbol.name;
                if (h == SYM_InexactNumbers) {
                    opt->domain = GB_DOM_INEXACT;
                    int64_t p = 16;
                    if (val->data.function.arg_count >= 1)
                        (void)expr_to_i64(val->data.function.args[0], &p);
                    /* Clamp to a sane positive range. */
                    if (p < 1)   p = 1;
                    if (p > 100000) p = 100000;
                    opt->inexact_prec = (int)p;
                } else if (h == SYM_Polynomials) {
                    /* Polynomials[x,...]: the listed symbols become
                     * coefficient-ring variables, i.e. parameters. */
                    opt->domain = GB_DOM_POLYNOMIALS;
                    opt->poly_domain_vars = val;
                } else {
                    warn_once("nimpl", "unsupported CoefficientDomain value; "
                                       "falling back to Rationals");
                    opt->domain = GB_DOM_RATIONALS;
                }
            } else {
                warn_once("nimpl", "unsupported CoefficientDomain value; "
                                   "falling back to Rationals");
                opt->domain = GB_DOM_RATIONALS;
            }
        } else if (key->data.symbol.name == SYM_Method) {
            bool ok = false;
            if (val->type == EXPR_SYMBOL
                && val->data.symbol.name == SYM_Automatic) ok = true;
            if (val->type == EXPR_SYMBOL
                && val->data.symbol.name == SYM_GroebnerWalk) {
                opt->method = GB_METHOD_WALK; ok = true;
            }
            if (val->type == EXPR_STRING && val->data.string
                && strcmp(val->data.string, "Buchberger") == 0) ok = true;
            if (val->type == EXPR_STRING && val->data.string
                && strcmp(val->data.string, "GroebnerWalk") == 0) {
                opt->method = GB_METHOD_WALK; ok = true;
            }
            if (!ok) {
                warn_once("nimpl", "only Method -> \"Buchberger\" or "
                                   "\"GroebnerWalk\" is implemented; "
                                   "falling back");
            }
        } else if (key->data.symbol.name == SYM_Modulus) {
            /* Honest "not implemented" instead of silently computing the
             * characteristic-0 basis as if Modulus weren't there.  The
             * surface accepts any value and ignores it; the basis comes
             * back over the rationals. */
            warn_once("modnotimpl", "Modulus option is not yet supported; "
                                    "the basis is being computed over the "
                                    "rationals.  Use Mathematica for "
                                    "modular bases.");
        } else if (key->data.symbol.name == SYM_Sort) {
            /* Mathematica's `Sort -> True` reverses the user-supplied
             * main-variable list before computing the basis, so the
             * "smallest" (rightmost-in-original-vars) variable becomes
             * the lex-leading one.  Empirically: with vars = {x, y, z}
             * and Sort -> True the basis is the same as with
             * vars = {z, y, x}.  We rebind the variable list at the
             * partition step below. */
            if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True) {
                opt->sort_desc = true;
            } else if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_False) {
                opt->sort_desc = false;
            } else {
                warn_once("nimpl", "Sort option accepts only True or False");
            }
        } else if (key->data.symbol.name == SYM_ParameterVariables) {
            /* The variables named here are treated as parameters: they
             * survive in coefficients of the basis polynomials but are
             * not allowed to appear as leading variables of any
             * surviving polynomial.  See the per-element discovery
             * pass below for the realisation. */
            opt->parameter_vars_given = true;
            opt->parameter_vars = val;
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Polynomial-equation preprocessing                                  */
/* ------------------------------------------------------------------ */

/* Rewrite `Equal[a, b]` as `a - b`, then run `Expand` so products of
 * polynomials (e.g. `(x - 1) (x - 2)`) reach a Plus-of-monomials shape
 * that `gb_from_expr` can ingest.  Other shapes are expanded as well
 * since input polynomials may carry factored sub-expressions.  Caller
 * owns the return. */
static Expr* normalise_polynomial(const Expr* p) {
    Expr* base;
    if (is_equal(p) && p->data.function.arg_count == 2) {
        Expr* a = expr_copy(p->data.function.args[0]);
        Expr* b = expr_copy(p->data.function.args[1]);
        base = internal_subtract((Expr*[]){ a, b }, 2);
    } else {
        base = expr_copy((Expr*)p);
    }
    /* The Gröbner engine ingests only exact (Integer/BigInt/Rational)
     * coefficients.  If the input carries inexact leaves (Real/MPFR) --
     * common with CoefficientDomain -> InexactNumbers, but possible in
     * any domain -- force-rationalise them first so gb_from_expr can
     * accept the polynomial. */
    if (internal_contains_inexact(base)) {
        Expr* rat = internal_force_rationalize(base);
        expr_free(base);
        base = rat;
    }
    /* Expand factored products / collect like terms.  Expand is also a
     * cheap no-op on already-expanded inputs. */
    Expr* expanded = internal_expand((Expr*[]){ base }, 1);
    /* Finally evaluate to fix-point so Plus/Times normalise. */
    Expr* normalised = evaluate(expanded);
    expr_free(expanded);
    return normalised;
}

/* ------------------------------------------------------------------ */
/*  Parameter discovery                                                */
/* ------------------------------------------------------------------ */

/* Names that look like free symbols but are really mathematical
 * constants (or evaluator sentinels).  They must NOT be auto-promoted
 * to parameter variables since they have a numeric value. */
static bool is_known_constant_symbol(const Expr* e) {
    if (!e || e->type != EXPR_SYMBOL) return false;
    const char* s = e->data.symbol.name;
    return  s == SYM_Pi
         || s == SYM_E
         || s == SYM_EulerGamma
         || s == SYM_Catalan
         || s == SYM_Degree
         || s == SYM_True
         || s == SYM_False
         || s == SYM_Infinity
         || s == SYM_DirectedInfinity
         || s == SYM_Null;
}

/* Set-like list of Expr* (interned-symbol identity, since symbol names
 * live in the global intern table).  Used to collect parameter symbols
 * during the polynomial walk. */
typedef struct {
    Expr** items;
    size_t n, cap;
} ExprSet;

static void exprset_init(ExprSet* s) { s->items = NULL; s->n = s->cap = 0; }

static void exprset_free(ExprSet* s, bool owned) {
    if (owned) for (size_t i = 0; i < s->n; i++) expr_free(s->items[i]);
    free(s->items);
}

static bool exprset_contains(const ExprSet* s, const Expr* e) {
    for (size_t i = 0; i < s->n; i++) if (expr_eq(s->items[i], (Expr*)e)) return true;
    return false;
}

static void exprset_push_borrowed(ExprSet* s, Expr* e) {
    if (exprset_contains(s, e)) return;
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 8;
        s->items = (Expr**)realloc(s->items, sizeof(Expr*) * s->cap);
    }
    s->items[s->n++] = e;
}

/* Walk `e`, collecting any EXPR_SYMBOL leaves that are
 *   - not a known mathematical constant, and
 *   - not present in any of the `excluded[*]` sets.
 * Results are pushed (in first-appearance order) into `out` as
 * BORROWED pointers — caller must outlive the source tree. */
static void collect_free_symbols(Expr* e,
                                 const ExprSet* const* excluded,
                                 size_t n_excluded,
                                 ExprSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_known_constant_symbol(e)) return;
        for (size_t i = 0; i < n_excluded; i++) {
            if (exprset_contains(excluded[i], e)) return;
        }
        exprset_push_borrowed(out, e);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    /* Do not descend into the head: a function head like `Sin` is not a
     * parameter, it's an operator name.  We assume `gb_from_expr` will
     * reject any unrecognised head shape later. */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        collect_free_symbols(e->data.function.args[i], excluded, n_excluded, out);
    }
}

/* Push every symbol element of a List[...] (or a bare symbol) into
 * `dst` as a borrowed pointer.  Non-symbol elements are silently
 * skipped — the caller is responsible for noticing shape errors. */
static void push_var_list(Expr* lst_or_sym, ExprSet* dst) {
    if (!lst_or_sym) return;
    if (lst_or_sym->type == EXPR_SYMBOL) {
        exprset_push_borrowed(dst, lst_or_sym);
        return;
    }
    if (!is_list(lst_or_sym)) return;
    for (size_t i = 0; i < lst_or_sym->data.function.arg_count; i++) {
        Expr* v = lst_or_sym->data.function.args[i];
        if (v->type == EXPR_SYMBOL) exprset_push_borrowed(dst, v);
    }
}

/* ------------------------------------------------------------------ */
/*  Weight-matrix parsing                                              */
/* ------------------------------------------------------------------ */

/* Extract a signed 64-bit integer from an Integer/BigInt Expr.  Returns
 * false for any other shape or a BigInt too large for int64. */
static bool expr_to_i64(const Expr* e, int64_t* out) {
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = (int64_t)mpz_get_si(e->data.bigint);
        return true;
    }
    return false;
}

/* Parse a `{{...}, {...}, ...}` value into a fresh row-major int64 matrix
 * (caller frees with free()).  Requires a non-empty rectangular list of
 * integer rows.  Returns NULL (and leaves rows/cols untouched) on any
 * shape mismatch or non-integer entry. */
static int64_t* parse_weight_matrix(const Expr* v, int* rows, int* cols) {
    if (!is_list(v)) return NULL;
    size_t r = v->data.function.arg_count;
    if (r == 0) return NULL;

    int ncols = -1;
    for (size_t i = 0; i < r; i++) {
        const Expr* row = v->data.function.args[i];
        if (!is_list(row)) return NULL;
        int c = (int)row->data.function.arg_count;
        if (c == 0) return NULL;
        if (ncols < 0) ncols = c;
        else if (c != ncols) return NULL;       /* ragged */
    }

    int64_t* m = (int64_t*)malloc(sizeof(int64_t) * r * (size_t)ncols);
    for (size_t i = 0; i < r; i++) {
        const Expr* row = v->data.function.args[i];
        for (int j = 0; j < ncols; j++) {
            if (!expr_to_i64(row->data.function.args[j],
                             &m[i * (size_t)ncols + j])) {
                free(m);
                return NULL;
            }
        }
    }
    *rows = (int)r;
    *cols = ncols;
    return m;
}

/* ------------------------------------------------------------------ */
/*  Builtin entry                                                      */
/* ------------------------------------------------------------------ */

Expr* builtin_groebner_basis(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) {
        fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis called with "
                        "0 arguments; 2 or 3 expected.\n");
        return NULL;
    }

    size_t n_pos;
    GBOptions opt;
    extract_options(res, &n_pos, &opt);

    /* The 1-arg form `GroebnerBasis[polys, ParameterVariables -> ...]`
     * is allowed: the main-variable list is auto-discovered from the
     * polynomials below.  Otherwise we want at least the (polys, vars)
     * pair.  3-arg = (polys, vars, elim_vars). */
    if (n_pos < 1 || n_pos > 3) {
        fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis takes 2 or 3 "
                        "positional arguments.\n");
        return NULL;
    }
    if (n_pos == 1 && !opt.parameter_vars_given) {
        fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis called with "
                        "1 argument; 2 or 3 expected (or pass "
                        "ParameterVariables -> ...).\n");
        return NULL;
    }

    Expr* polys_list = res->data.function.args[0];
    Expr* vars_list  = (n_pos >= 2) ? res->data.function.args[1] : NULL;
    Expr* elim_list  = (n_pos == 3) ? res->data.function.args[2] : NULL;

    if (!is_list(polys_list)) return NULL;

    /* Accept a single-symbol shorthand for the variable list:
     *   GroebnerBasis[polys, x]    -> GroebnerBasis[polys, {x}]
     *   GroebnerBasis[polys, m, x] -> GroebnerBasis[polys, m, {x}]
     * We wrap symbol args in a synthetic List for uniform downstream
     * processing.  The wrappers are freed at the end of the function. */
    Expr* wrap_vars = NULL;
    Expr* wrap_elim = NULL;
    if (vars_list && !is_list(vars_list)) {
        if (vars_list->type != EXPR_SYMBOL) return NULL;
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        wrapped[0] = expr_copy(vars_list);
        wrap_vars = expr_new_function(expr_new_symbol(SYM_List), wrapped, 1);
        free(wrapped);
        vars_list = wrap_vars;
    }
    if (n_pos == 3 && !is_list(elim_list)) {
        if (elim_list->type != EXPR_SYMBOL) {
            if (wrap_vars) expr_free(wrap_vars);
            return NULL;
        }
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        wrapped[0] = expr_copy(elim_list);
        wrap_elim = expr_new_function(expr_new_symbol(SYM_List), wrapped, 1);
        free(wrapped);
        elim_list = wrap_elim;
    }

    /* Empty polynomial list -> empty basis. */
    size_t n_polys = polys_list->data.function.arg_count;
    if (n_polys == 0) {
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    /* ---- Parameter and main-variable discovery. -------------------- */
    /*
     * Two complementary inputs steer the partitioning:
     *
     *   - explicit main vars (`vars_list`) and elim vars (`elim_list`).
     *   - the `ParameterVariables -> ...` option (`opt.parameter_vars`).
     *
     * Any symbol mentioned by the polynomials that is not in any of
     * those sets and not a known mathematical constant is auto-
     * promoted to a parameter ("issue 1" / "issue 6" behaviour:
     * `GroebnerBasis[{a x^2 + 5 x - 1, ...}, {x, y}]` treats `a` as a
     * parameter instead of refusing to evaluate).  Parameters become
     * additional joint-array variables that the Buchberger engine
     * threads through; after the run, polynomials that mention no
     * main- or elim-block variable are filtered out as
     * annihilator-trivial.  The remaining polynomials carry the
     * parameters naturally as factor terms.
     */
    ExprSet param_explicit;  exprset_init(&param_explicit);
    ExprSet main_set;        exprset_init(&main_set);
    ExprSet elim_set;        exprset_init(&elim_set);

    if (opt.parameter_vars_given) push_var_list(opt.parameter_vars, &param_explicit);
    if (vars_list)                push_var_list(vars_list,           &main_set);
    if (elim_list)                push_var_list(elim_list,           &elim_set);

    /* CoefficientDomain -> Polynomials[x,...]: the named symbols are
     * coefficient-ring variables, which is exactly the parameter role.
     * Fold them into the explicit-parameter set (the Polynomials head is
     * a plain function whose args are the symbols). */
    if (opt.domain == GB_DOM_POLYNOMIALS && opt.poly_domain_vars) {
        for (size_t i = 0; i < opt.poly_domain_vars->data.function.arg_count; i++) {
            Expr* v = opt.poly_domain_vars->data.function.args[i];
            if (v->type == EXPR_SYMBOL) exprset_push_borrowed(&param_explicit, v);
        }
    }

    /* Symbols the user has named as parameters take priority over any
     * accidental presence in the main- or elim-var list.  Drop them
     * from those two sets so the joint var array stays partitioned. */
    if (param_explicit.n > 0) {
        size_t w = 0;
        for (size_t i = 0; i < main_set.n; i++) {
            if (!exprset_contains(&param_explicit, main_set.items[i])) {
                main_set.items[w++] = main_set.items[i];
            }
        }
        main_set.n = w;
        w = 0;
        for (size_t i = 0; i < elim_set.n; i++) {
            if (!exprset_contains(&param_explicit, elim_set.items[i])) {
                elim_set.items[w++] = elim_set.items[i];
            }
        }
        elim_set.n = w;
    }

    /* Auto-discover parameters: any free symbol in polys that is not
     * in main, elim, explicit-params, or a known constant.  Also
     * auto-discover main vars when none were given (the 1-arg
     * `ParameterVariables -> ...` shorthand). */
    ExprSet param_auto;  exprset_init(&param_auto);
    {
        const ExprSet* excl[3] = { &main_set, &elim_set, &param_explicit };
        for (size_t i = 0; i < n_polys; i++) {
            collect_free_symbols(polys_list->data.function.args[i],
                                 excl, 3, &param_auto);
        }
    }

    /* The combined parameter set, in first-appearance order
     * (explicit first, then auto-discovered). */
    ExprSet params;  exprset_init(&params);
    for (size_t i = 0; i < param_explicit.n; i++)
        exprset_push_borrowed(&params, param_explicit.items[i]);
    for (size_t i = 0; i < param_auto.n; i++)
        exprset_push_borrowed(&params, param_auto.items[i]);

    /* `GroebnerBasis[polys, ParameterVariables -> p]`: the main vars
     * are auto-derived as "every free symbol in polys that is not p". */
    if (main_set.n == 0 && !vars_list && param_explicit.n > 0) {
        const ExprSet* excl[3] = { &param_explicit, &elim_set, NULL };
        for (size_t i = 0; i < n_polys; i++) {
            collect_free_symbols(polys_list->data.function.args[i],
                                 excl, 2, &main_set);
        }
    }

    /* Sort -> True reverses the main-variable list before the joint
     * array is built, matching Mathematica's empirical behaviour
     * (`GroebnerBasis[polys, {x, y, z}, Sort -> True]` returns the
     * same basis as `GroebnerBasis[polys, {z, y, x}]`). */
    if (opt.sort_desc && main_set.n > 1) {
        for (size_t i = 0, j = main_set.n - 1; i < j; i++, j--) {
            Expr* tmp = main_set.items[i];
            main_set.items[i] = main_set.items[j];
            main_set.items[j] = tmp;
        }
    }

    size_t n_main = main_set.n;
    size_t n_elim = elim_set.n;
    size_t n_params = params.n;
    size_t n_vars = n_main + n_elim + n_params;
    if (n_vars == 0) {
        exprset_free(&param_explicit, false);
        exprset_free(&main_set, false);
        exprset_free(&elim_set, false);
        exprset_free(&param_auto, false);
        exprset_free(&params, false);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return NULL;
    }

    /* Joint var array: [elim..., main..., params...].  Elim is the
     * lex-leading block (preserves the elimination-theorem semantics);
     * params go last so they act as the lowest-priority symbols and the
     * Buchberger run preserves their role as "coefficient carriers". */
    Expr** all_vars = (Expr**)malloc(sizeof(Expr*) * n_vars);
    for (size_t i = 0; i < n_elim; i++)   all_vars[i] = elim_set.items[i];
    for (size_t i = 0; i < n_main; i++)   all_vars[n_elim + i] = main_set.items[i];
    for (size_t i = 0; i < n_params; i++) all_vars[n_elim + n_main + i] = params.items[i];
    int elim_pivot = (int)n_elim;

    /* The 3-arg form forces Lex (with elim variables placed first) so
     * the elimination theorem applies: any monomial involving any elim
     * variable is lex-larger than any monomial without one.  Likewise
     * the parametric path forces Lex so the params-last placement
     * orders parametric coefficients below main-variable monomials and
     * the post-run "drop pure-parameter polynomials" filter is
     * mathematically the elimination-ideal contraction. */
    GBOrder use_order = (n_elim > 0 || n_params > 0) ? GB_ORDER_LEX : opt.order;

    /* Weight-matrix MonomialOrder: honoured only in the plain form (no
     * elimination, no parameters) where the joint var array equals the
     * main-variable list, so the user matrix maps directly.  The matrix
     * must be rectangular, integer, have exactly n_main columns, and
     * define a valid term order. */
    GBWeightMatrix wmat_storage;
    const GBWeightMatrix* wmat_ptr = NULL;
    int64_t* wmat_buf = NULL;
    if (opt.matrix_order) {
        if (n_elim == 0 && n_params == 0) {
            int wr = 0, wc = 0;
            int64_t* buf = parse_weight_matrix(opt.matrix_order, &wr, &wc);
            if (buf && wc == (int)n_main && gb_wmat_validate(buf, wr, wc)) {
                wmat_buf = buf;
                wmat_storage.n_rows = wr;
                wmat_storage.n_vars = wc;
                wmat_storage.w = wmat_buf;
                wmat_ptr = &wmat_storage;
                use_order = GB_ORDER_MATRIX;
            } else {
                free(buf);
                warn_once("nimpl", "MonomialOrder matrix is not a valid term "
                                   "order for these variables; falling back "
                                   "to Lexicographic");
            }
        } else {
            warn_once("nimpl", "weight-matrix MonomialOrder is not supported "
                               "with elimination or parameter variables; "
                               "falling back to Lexicographic");
        }
    }

    /* Convert each input polynomial. */
    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * n_polys);
    size_t nF = 0;
    bool failed = false;
    for (size_t i = 0; i < n_polys; i++) {
        Expr* norm = normalise_polynomial(polys_list->data.function.args[i]);
        GBPoly* p = gb_from_expr(norm, all_vars, (int)n_vars,
                                 use_order, elim_pivot, wmat_ptr);
        expr_free(norm);
        if (!p) { failed = true; break; }
        if (p->n_terms == 0) { gb_poly_free(p); continue; }
        F[nF++] = p;
    }
    if (failed) {
        for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        free(all_vars);
        free(wmat_buf);
        exprset_free(&param_explicit, false);
        exprset_free(&main_set, false);
        exprset_free(&elim_set, false);
        exprset_free(&param_auto, false);
        exprset_free(&params, false);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return NULL;
    }

    if (nF == 0) {
        /* All inputs reduced to zero -> empty basis. */
        free(F);
        free(all_vars);
        free(wmat_buf);
        exprset_free(&param_explicit, false);
        exprset_free(&main_set, false);
        exprset_free(&elim_set, false);
        exprset_free(&param_auto, false);
        exprset_free(&params, false);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    /* Fast path: any input is a non-zero constant -> ideal = <1>.  Not
     * valid over the integers, where a constant c generates only c*Z[x]
     * (the strong basis keeps the constant), so skip it for Integers. */
    for (size_t i = 0; i < nF && opt.domain != GB_DOM_INTEGERS; i++) {
        if (gb_poly_is_constant(F[i])) {
            for (size_t j = 0; j < nF; j++) gb_poly_free(F[j]);
            free(F);
            free(all_vars);
            free(wmat_buf);
            exprset_free(&param_explicit, false);
            exprset_free(&main_set, false);
            exprset_free(&elim_set, false);
            exprset_free(&param_auto, false);
            exprset_free(&params, false);
            if (wrap_vars) expr_free(wrap_vars);
            if (wrap_elim) expr_free(wrap_elim);
            Expr** one = (Expr**)malloc(sizeof(Expr*));
            one[0] = expr_new_integer(1);
            Expr* lst = expr_new_function(expr_new_symbol(SYM_List), one, 1);
            free(one);
            return lst;
        }
    }

    /* Dispatch to the selected engine.  GroebnerWalk is honoured only in
     * the plain form (no elimination, no parameters); otherwise the
     * params-last / elim-first joint order required by the contraction
     * filter does not match the walk's source order, so we run Buchberger
     * directly (the reduced basis is identical either way). */
    size_t out_n = 0;
    GBPoly** G;
    if (opt.domain == GB_DOM_INTEGERS) {
        /* Strong Gröbner basis over Z (parameters are ordinary variables
         * in the PID polynomial ring Z[vars]). */
        G = gb_strong_buchberger(F, nF, &out_n);
    } else if (opt.method == GB_METHOD_WALK && n_elim == 0 && n_params == 0) {
        G = gb_groebner_walk(F, nF, use_order, wmat_ptr, &out_n);
    } else {
        G = gb_buchberger(F, nF, &out_n);
    }

    /* The weight matrix is no longer needed: the returned basis only feeds
     * gb_to_expr (which reads exponents/coefficients, not the order). */
    free(wmat_buf);

    /* Free input working set. */
    for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);

    /* Elimination filter: drop polynomials that still mention any elim
     * variable.  By the elimination theorem, the surviving polynomials
     * form a Gröbner basis of the elimination ideal. */
    if (n_elim > 0) {
        int* elim_idx = (int*)malloc(sizeof(int) * n_elim);
        for (size_t i = 0; i < n_elim; i++) elim_idx[i] = (int)i;
        size_t k = 0;
        for (size_t i = 0; i < out_n; i++) {
            if (gb_poly_free_of_vars(G[i], elim_idx, (int)n_elim)) {
                G[k++] = G[i];
            } else {
                gb_poly_free(G[i]);
            }
        }
        out_n = k;
        free(elim_idx);
    }

    /* Note: we INTENTIONALLY do NOT drop polynomials whose every term
     * has exponent zero in every main / elim variable.  Mathematica
     * keeps such "pure-parameter" polynomials in the basis when they
     * arise -- they are precisely the consistency conditions on the
     * parameter values that an over-determined system imposes (see the
     * `1625 x^8 + ... + 2` element returned by the issue-6 case
     * `GroebnerBasis[{...}, {y, z}]`).  The Buchberger run over
     * Z[params][main_vars] under params-last lex already produces them
     * as constants in the (y, z) coefficient ring. */

    /* CoefficientDomain -> RationalFunctions: the parameters live in the
     * coefficient FIELD Q(params), so they are units.  Autoreduce the ring
     * basis into the (smaller) Gröbner basis over Q(params)[main_vars],
     * dropping generators whose main leading monomial is redundant.  Only
     * meaningful when parameters are present (Q() over no params is just
     * Q) and in the plain form (no elimination block). */
    if (opt.domain == GB_DOM_RATIONAL_FUNCTIONS && n_params > 0
        && n_elim == 0) {
        int n_field_main = (int)(n_vars - n_params);
        gb_rational_function_reduce(G, &out_n, n_field_main);
    }

    /* If the basis collapses to {<non-zero constant>} -> {1}.  Over the
     * integers a constant generator c means the ideal is c*Z[x] (not <1>
     * unless |c| == 1), so the strong basis keeps the constant verbatim. */
    bool has_const = false;
    for (size_t i = 0; i < out_n && opt.domain != GB_DOM_INTEGERS; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            has_const = true; break;
        }
    }
    if (has_const) {
        gb_basis_free(G, out_n);
        free(all_vars);
        exprset_free(&param_explicit, false);
        exprset_free(&main_set, false);
        exprset_free(&elim_set, false);
        exprset_free(&param_auto, false);
        exprset_free(&params, false);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        Expr** one = (Expr**)malloc(sizeof(Expr*));
        one[0] = expr_new_integer(1);
        Expr* lst = expr_new_function(expr_new_symbol(SYM_List), one, 1);
        free(one);
        return lst;
    }

    /* Convert each basis polynomial back to an Expr.  Parameters
     * naturally appear as factor terms because they live in the joint
     * `all_vars` array and gb_to_expr emits one factor per non-zero
     * exponent slot. */
    Expr** items = (out_n > 0) ? (Expr**)malloc(sizeof(Expr*) * out_n) : NULL;
    if (opt.domain == GB_DOM_INEXACT) {
        /* Approximate-arithmetic output: monic, MPFR coefficients at the
         * requested decimal precision.  bits = ceil(prec * log2(10)). */
        mpfr_prec_t bits = (mpfr_prec_t)((opt.inexact_prec * 33219L) / 10000L) + 1;
        for (size_t i = 0; i < out_n; i++) {
            items[i] = gb_to_expr_inexact(G[i], all_vars, bits);
        }
    } else {
        for (size_t i = 0; i < out_n; i++) {
            items[i] = gb_to_expr(G[i], all_vars);
        }
    }
    gb_basis_free(G, out_n);
    free(all_vars);
    exprset_free(&param_explicit, false);
    exprset_free(&main_set, false);
    exprset_free(&elim_set, false);
    exprset_free(&param_auto, false);
    exprset_free(&params, false);
    if (wrap_vars) expr_free(wrap_vars);
    if (wrap_elim) expr_free(wrap_elim);

    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), items, out_n);
    free(items);
    return lst;
}

/* ------------------------------------------------------------------ */
/*  Initialisation                                                     */
/* ------------------------------------------------------------------ */

void groebner_init(void) {
    symtab_add_builtin("GroebnerBasis", builtin_groebner_basis);
    SymbolDef* def = symtab_get_def("GroebnerBasis");
    if (def) def->attributes |= ATTR_PROTECTED;
}
