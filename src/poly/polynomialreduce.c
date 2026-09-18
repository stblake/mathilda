/* polynomialreduce.c
 *
 * `PolynomialReduce[poly, {p1, ..., pn}, {x1, ..., xk}]` -- multivariate
 * polynomial division with cofactors.  Returns {{a1, ..., an}, b} such that
 *
 *     a1 p1 + a2 p2 + ... + an pn + b == poly   (exactly)
 *
 * and b is fully reduced: no term of b is divisible by any leading term of the
 * pi under the chosen monomial order.  This is the normal-form division that
 * sits at the heart of Buchberger's algorithm, exposed as a user builtin; it is
 * the multivariate-division sibling of GroebnerBasis and shares its options.
 *
 * Three coefficient domains are implemented:
 *
 *   - Rationals over Q (and RationalFunctions with no parameters, which is the
 *     same field): the exact GBPoly engine gb_divmod (groebner.c), with a FLINT
 *     fast path (fmpq_mpoly_divrem_ideal) for the Lexicographic case on a
 *     FLINT-enabled build.  gb_divmod is the authoritative, Wolfram-matching
 *     engine.
 *   - RationalFunctions over Q(params): when free symbols outside the variable
 *     list appear, they are coefficient-field elements, so the cofactors can
 *     carry denominators in the parameters (e.g. y/(4 a)).  The GBPoly engine
 *     (rational coefficients only) cannot express those; this case runs an
 *     Expr-coefficient division whose coefficient arithmetic is Together/Cancel
 *     over Q(params) -- exact, no numeric oracle.
 *   - Modulus -> p (GF(p)): via the gbmod.c gfp_divmod engine.
 *
 * CoefficientDomain -> Integers and InexactNumbers (+ Tolerance) are distinct
 * ring/numeric engines and are DEFERRED: the option is accepted but the head
 * declines (returns NULL) with a note, so behaviour is never silently wrong.
 *
 * Attributes: Protected.  Options: MonomialOrder, CoefficientDomain, Modulus,
 * ParameterVariables, Tolerance (defaults registered in options_builtin.c).
 *
 * PolynomialReduce is a symbolic/structural head (it returns lists of symbolic
 * polynomials, not element-wise machine numbers over a buffer), so it is
 * genuinely exempt from the packed/NDArray and Compile[] numeric surfaces --
 * there is nothing element-wise to lower.
 */

#include "poly.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "internal.h"
#include "groebner.h"
#include "gbmod.h"
#include "flint_bridge.h"
#include "rationalize.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/*  Small structural predicates                                        */
/* ------------------------------------------------------------------ */

static bool pr_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

static bool pr_is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_Rule
         || e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

static void pr_warn(const char* tag, const char* msg) {
    fprintf(stderr, "PolynomialReduce::%s: %s\n", tag, msg);
}

/* ------------------------------------------------------------------ */
/*  Options                                                            */
/* ------------------------------------------------------------------ */

typedef enum {
    PR_DOM_RATIONALS = 0,        /* Rationals */
    PR_DOM_RATIONAL_FUNCTIONS,   /* RationalFunctions / Automatic (default) */
    PR_DOM_DEFERRED              /* Integers / InexactNumbers / Polynomials */
} PRDomain;

typedef struct {
    const Expr* monomial_order;  /* borrowed; MonomialOrder value or NULL */
    PRDomain    domain;
    int64_t     modulus;         /* Modulus -> p (prime); 0 = char 0 */
    const Expr* parameter_vars;  /* borrowed; ParameterVariables value or NULL */
    bool        tolerance_nonzero;
} PROptions;

/* Extract options; sets *n_pos to the positional-argument count. */
static void pr_extract_options(const Expr* res, size_t* n_pos, PROptions* opt) {
    size_t argc = res->data.function.arg_count;
    size_t cut = argc;
    while (cut > 0 && pr_is_rule(res->data.function.args[cut - 1])) cut--;
    *n_pos = cut;

    opt->monomial_order    = NULL;
    opt->domain            = PR_DOM_RATIONAL_FUNCTIONS;
    opt->modulus           = 0;
    opt->parameter_vars    = NULL;
    opt->tolerance_nonzero = false;

    for (size_t i = cut; i < argc; i++) {
        Expr* rule = res->data.function.args[i];
        if (rule->data.function.arg_count != 2) continue;
        Expr* key = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (key->type != EXPR_SYMBOL) continue;

        if (key->data.symbol.name == SYM_MonomialOrder) {
            opt->monomial_order = val;
        } else if (key->data.symbol.name == SYM_CoefficientDomain) {
            if (val->type == EXPR_SYMBOL) {
                if (val->data.symbol.name == SYM_Rationals) {
                    opt->domain = PR_DOM_RATIONALS;
                } else if (val->data.symbol.name == SYM_RationalFunctions
                        || val->data.symbol.name == SYM_Automatic) {
                    opt->domain = PR_DOM_RATIONAL_FUNCTIONS;
                } else {
                    /* Integers, InexactNumbers, ... : deferred. */
                    opt->domain = PR_DOM_DEFERRED;
                }
            } else {
                /* InexactNumbers[prec], Polynomials[...], ... : deferred. */
                opt->domain = PR_DOM_DEFERRED;
            }
        } else if (key->data.symbol.name == SYM_Modulus) {
            if (val->type == EXPR_INTEGER && val->data.integer >= 2
                && val->data.integer < (int64_t)1 << 31
                && gfp_is_prime((uint64_t)val->data.integer)) {
                opt->modulus = val->data.integer;
            } else if (!(val->type == EXPR_INTEGER && val->data.integer == 0)) {
                pr_warn("modnotimpl", "Modulus -> p requires a prime p in "
                                      "[2, 2^31); computing over the rationals.");
            }
        } else if (key->data.symbol.name == SYM_ParameterVariables) {
            opt->parameter_vars = val;
        } else if (key->data.symbol.name == SYM_Tolerance) {
            if (!(val->type == EXPR_INTEGER && val->data.integer == 0))
                opt->tolerance_nonzero = true;   /* nonzero tolerance -> deferred */
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Monomial-order resolution                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    bool     ok;          /* recognised / valid */
    bool     is_matrix;   /* explicit user weight matrix */
    int      ds, rev, neg;/* named-order flags (when !is_matrix) */
    GBOrder  native;      /* LEX / GREVLEX / MATRIX hint (when !is_matrix) */
    int      flint_kind;  /* 0 Lex, 1 DegLex, 2 DegRevLex, -1 none */
    int64_t* mat;         /* explicit matrix rows*k (owned; when is_matrix) */
    int      mat_rows;
} PROrder;

/* Parse a `{{...}, ...}` value into a fresh rows*k row-major int64 matrix. */
static int64_t* pr_parse_matrix(const Expr* v, int k, int* rows_out) {
    if (!pr_is_list(v)) return NULL;
    size_t r = v->data.function.arg_count;
    if (r == 0) return NULL;
    int64_t* m = (int64_t*)malloc(sizeof(int64_t) * r * (size_t)(k > 0 ? k : 1));
    for (size_t i = 0; i < r; i++) {
        const Expr* row = v->data.function.args[i];
        if (!pr_is_list(row) || (int)row->data.function.arg_count != k) {
            free(m); return NULL;
        }
        for (int j = 0; j < k; j++) {
            const Expr* e = row->data.function.args[j];
            if (e->type == EXPR_INTEGER) {
                m[i * (size_t)k + j] = e->data.integer;
            } else {
                free(m); return NULL;
            }
        }
    }
    *rows_out = (int)r;
    return m;
}

static void pr_resolve_order(const Expr* mo, int k, PROrder* out) {
    out->ok = false; out->is_matrix = false;
    out->ds = out->rev = out->neg = 0;
    out->native = GB_ORDER_LEX; out->flint_kind = 0;
    out->mat = NULL; out->mat_rows = 0;

    if (!mo) { out->ok = true; return; }            /* default Lexicographic */

    if (mo->type == EXPR_SYMBOL || mo->type == EXPR_STRING) {
        const char* name = (mo->type == EXPR_SYMBOL) ? mo->data.symbol.name
                                                     : mo->data.string;
        /* Automatic -> Lexicographic. */
        if (mo->type == EXPR_SYMBOL && mo->data.symbol.name == SYM_Automatic) {
            out->ok = true; out->native = GB_ORDER_LEX; out->flint_kind = 0;
            return;
        }
        int ds = 0, rv = 0, nl = 0; GBOrder hint = GB_ORDER_LEX;
        if (gb_classify_named_order(name, &ds, &rv, &nl, &hint)) {
            out->ok = true;
            out->ds = ds; out->rev = rv; out->neg = nl; out->native = hint;
            if (hint == GB_ORDER_LEX) {
                out->flint_kind = 0;                                 /* Lex */
            } else if (hint == GB_ORDER_GREVLEX) {
                out->flint_kind = 2;                                 /* DegRevLex */
            } else {
                /* DegreeLexicographic / Negative* need a weight matrix.  Build
                 * and validate it: the Negative family is not a well-founded
                 * (global) term order, so -- like the GB engine -- downgrade it
                 * to Lexicographic rather than risk a non-terminating reduction. */
                int rows = 0;
                int64_t* W = gb_build_order_matrix(k, ds, rv, nl, &rows);
                if (W && gb_wmat_validate(W, rows, k)) {
                    out->is_matrix = true; out->mat = W; out->mat_rows = rows;
                    out->native = GB_ORDER_MATRIX;
                    out->flint_kind = (ds == 1 && rv == 0 && nl == 0) ? 1 : -1;
                } else {
                    free(W);
                    pr_warn("nimpl", "MonomialOrder is not a well-founded term "
                                     "order for these variables; using "
                                     "Lexicographic.");
                    out->ds = out->rev = out->neg = 0;
                    out->native = GB_ORDER_LEX; out->flint_kind = 0;
                }
            }
        }
        return;
    }

    if (pr_is_list(mo)) {
        int rows = 0;
        int64_t* m = pr_parse_matrix(mo, k, &rows);
        if (m && gb_wmat_validate(m, rows, k)) {
            out->ok = true; out->is_matrix = true;
            out->mat = m; out->mat_rows = rows; out->native = GB_ORDER_MATRIX;
            out->flint_kind = -1;
        } else {
            free(m);
        }
        return;
    }
}

/* ------------------------------------------------------------------ */
/*  Polynomial input normalisation                                     */
/* ------------------------------------------------------------------ */

/* Expand, rationalise inexact leaves, evaluate to fixpoint -- so gb_from_expr /
 * CoefficientRules ingest a clean Plus-of-monomials tree.  Caller owns result. */
static Expr* pr_normalise(const Expr* p) {
    Expr* base = expr_copy((Expr*)p);
    if (internal_contains_inexact(base)) {
        Expr* rat = internal_force_rationalize(base);
        expr_free(base);
        base = rat;
    }
    Expr* expanded = internal_expand((Expr*[]){ base }, 1);   /* consumes base */
    Expr* norm = evaluate(expanded);
    expr_free(expanded);
    return norm;
}

/* ------------------------------------------------------------------ */
/*  Variable / parameter discovery                                     */
/* ------------------------------------------------------------------ */

/* All free symbols (constants excluded) of {poly, divisors...}, canonical
 * order, via Variables[].  Returns an owned List Expr; caller frees. */
static Expr* pr_all_variables(const Expr* poly, Expr* divisors_list) {
    size_t nd = divisors_list->data.function.arg_count;
    Expr** items = (Expr**)malloc(sizeof(Expr*) * (nd + 1));
    items[0] = expr_copy((Expr*)poly);
    for (size_t i = 0; i < nd; i++)
        items[i + 1] = expr_copy(divisors_list->data.function.args[i]);
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), items, nd + 1);
    free(items);
    Expr* vars = internal_variables((Expr*[]){ lst }, 1);    /* consumes lst */
    return vars;   /* a List of symbols */
}

static bool sym_in_list(const Expr* sym, Expr* const* arr, int n) {
    for (int i = 0; i < n; i++) if (expr_eq((Expr*)sym, arr[i])) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/*  Field arithmetic over Q(params) (coefficients are Expr)            */
/* ------------------------------------------------------------------ */

/* Together[e] -> a single reduced fraction.  Consumes e; returns owned. */
static Expr* field_norm(Expr* e) {
    Expr* call = expr_new_function(expr_new_symbol("Together"),
                                   (Expr*[]){ e }, 1);
    return eval_and_free(call);
}

static Expr* field_add(const Expr* a, const Expr* b) {
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                    (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2);
    return field_norm(sum);
}

static Expr* field_mul(const Expr* a, const Expr* b) {
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2);
    return field_norm(prod);
}

/* a / b in the field, reduced to lowest terms. */
static Expr* field_div(const Expr* a, const Expr* b) {
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy((Expr*)b), expr_new_integer(-1) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy((Expr*)a), inv }, 2);
    return field_norm(prod);
}

/* Negate; CONSUMES `a`. */
static Expr* field_neg(Expr* a) {
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), a }, 2);
    return field_norm(prod);   /* Together canonicalises (distributes the -1) */
}

static bool coeff_is_zero(const Expr* c) {
    if (c->type == EXPR_INTEGER) return c->data.integer == 0;
    Expr* t = field_norm(expr_copy((Expr*)c));
    bool z = (t->type == EXPR_INTEGER && t->data.integer == 0);
    expr_free(t);
    return z;
}

/* ------------------------------------------------------------------ */
/*  Sparse polynomial with field coefficients (RationalFunctions path) */
/* ------------------------------------------------------------------ */

typedef struct { int* e; Expr* c; } RTerm;      /* e: k exps (owned); c owned */
typedef struct { RTerm* t; size_t n, cap; int k; } RPoly;

static void rpoly_init(RPoly* p, int k) { p->t = NULL; p->n = p->cap = 0; p->k = k; }

static void rpoly_free(RPoly* p) {
    for (size_t i = 0; i < p->n; i++) { free(p->t[i].e); expr_free(p->t[i].c); }
    free(p->t);
    p->t = NULL; p->n = p->cap = 0;
}

/* Add c * x^exps into p (merging like monomials).  Takes ownership of `c`;
 * copies `exps`.  Drops a term whose coefficient reduces to 0. */
static void rpoly_add_term(RPoly* p, const int* exps, Expr* c) {
    if (coeff_is_zero(c)) { expr_free(c); return; }
    for (size_t i = 0; i < p->n; i++) {
        if (memcmp(p->t[i].e, exps, sizeof(int) * (size_t)p->k) == 0) {
            Expr* sum = field_add(p->t[i].c, c);
            expr_free(p->t[i].c);
            expr_free(c);
            if (coeff_is_zero(sum)) {
                expr_free(sum);
                free(p->t[i].e);
                p->t[i] = p->t[p->n - 1];   /* swap-remove */
                p->n--;
            } else {
                p->t[i].c = sum;
            }
            return;
        }
    }
    if (p->n == p->cap) {
        p->cap = p->cap ? p->cap * 2 : 8;
        p->t = (RTerm*)realloc(p->t, sizeof(RTerm) * p->cap);
    }
    p->t[p->n].e = (int*)malloc(sizeof(int) * (size_t)(p->k > 0 ? p->k : 1));
    memcpy(p->t[p->n].e, exps, sizeof(int) * (size_t)p->k);
    p->t[p->n].c = c;
    p->n++;
}

/* Descending compare of exponent vectors under weight matrix W (rows*k). */
static int pr_wcompare(const int* a, const int* b, const int64_t* W,
                       int rows, int k) {
    for (int r = 0; r < rows; r++) {
        long long sa = 0, sb = 0;
        for (int j = 0; j < k; j++) {
            sa += (long long)W[(size_t)r * k + j] * a[j];
            sb += (long long)W[(size_t)r * k + j] * b[j];
        }
        if (sa != sb) return (sa > sb) ? -1 : 1;   /* larger monomial first */
    }
    return 0;
}

static const int64_t* g_sort_W; static int g_sort_rows, g_sort_k;
static int rterm_cmp(const void* a, const void* b) {
    const RTerm* x = (const RTerm*)a;
    const RTerm* y = (const RTerm*)b;
    return pr_wcompare(x->e, y->e, g_sort_W, g_sort_rows, g_sort_k);
}

static void rpoly_sort(RPoly* p, const int64_t* W, int rows) {
    g_sort_W = W; g_sort_rows = rows; g_sort_k = p->k;
    if (p->n > 1) qsort(p->t, p->n, sizeof(RTerm), rterm_cmp);
}

/* Build an RPoly from an Expr polynomial-in-mainvars, coefficients in the
 * params field, via CoefficientRules[expr, {mainvars}].  Returns false on a
 * non-polynomial input (the evaluated CoefficientRules stayed a function). */
static bool rpoly_from_expr(const Expr* poly, Expr* const* mainvars, int k,
                            RPoly* out) {
    rpoly_init(out, k);
    /* {mainvars} */
    Expr** vv = (Expr**)malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
    for (int i = 0; i < k; i++) vv[i] = expr_copy(mainvars[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vv, (size_t)k);
    free(vv);
    Expr* call = expr_new_function(expr_new_symbol("CoefficientRules"),
                    (Expr*[]){ expr_copy((Expr*)poly), vlist }, 2);
    Expr* rules = eval_and_free(call);
    if (!pr_is_list(rules)) { expr_free(rules); return false; }
    bool ok = true;
    for (size_t i = 0; i < rules->data.function.arg_count && ok; i++) {
        Expr* rule = rules->data.function.args[i];
        if (!pr_is_rule(rule) || rule->data.function.arg_count != 2) { ok = false; break; }
        Expr* ev = rule->data.function.args[0];   /* List of k exponents */
        Expr* cf = rule->data.function.args[1];
        if (!pr_is_list(ev) || (int)ev->data.function.arg_count != k) { ok = false; break; }
        int* exps = (int*)malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
        for (int j = 0; j < k; j++) {
            Expr* ee = ev->data.function.args[j];
            exps[j] = (ee->type == EXPR_INTEGER) ? (int)ee->data.integer : 0;
        }
        rpoly_add_term(out, exps, field_norm(expr_copy(cf)));
        free(exps);
    }
    expr_free(rules);
    if (!ok) { rpoly_free(out); return false; }
    return true;
}

/* Render an RPoly back to an Expr (Plus of Times[coeff, mono]); 0 if empty. */
static Expr* rpoly_to_expr(const RPoly* p, Expr* const* mainvars) {
    if (p->n == 0) return expr_new_integer(0);
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * p->n);
    for (size_t i = 0; i < p->n; i++) {
        Expr** fac = (Expr**)malloc(sizeof(Expr*) * (size_t)(p->k + 1));
        size_t nf = 0;
        fac[nf++] = expr_copy(p->t[i].c);
        for (int v = 0; v < p->k; v++) {
            int ev = p->t[i].e[v];
            if (ev == 0) continue;
            if (ev == 1) fac[nf++] = expr_copy(mainvars[v]);
            else fac[nf++] = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(mainvars[v]), expr_new_integer(ev) }, 2);
        }
        terms[i] = (nf == 1) ? fac[0]
                             : expr_new_function(expr_new_symbol(SYM_Times), fac, nf);
        free(fac);
    }
    Expr* out = (p->n == 1) ? terms[0]
                            : expr_new_function(expr_new_symbol(SYM_Plus), terms, p->n);
    free(terms);
    return out;
}

static bool exp_divides(const int* lm, const int* e, int* q, int k) {
    for (int i = 0; i < k; i++) {
        if (e[i] < lm[i]) return false;
        q[i] = e[i] - lm[i];
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  RationalFunctions engine: division with cofactors over Q(params)   */
/* ------------------------------------------------------------------ */

/* On success returns a freshly-built List {List{a1,...,an}, b}.  divisors and
 * poly are the (normalised) Exprs; mainvars are the k division variables. */
static Expr* pr_reduce_field(const Expr* poly, Expr* const* divs, int ndiv,
                             Expr* const* mainvars, int k,
                             const int64_t* W, int rows) {
    RPoly r;
    if (!rpoly_from_expr(poly, mainvars, k, &r)) return NULL;

    RPoly* D = (RPoly*)malloc(sizeof(RPoly) * (size_t)ndiv);
    RPoly* Q = (RPoly*)malloc(sizeof(RPoly) * (size_t)ndiv);
    bool ok = true;
    int built = 0;
    for (int i = 0; i < ndiv; i++) {
        if (!rpoly_from_expr(divs[i], mainvars, k, &D[i])) { ok = false; break; }
        rpoly_sort(&D[i], W, rows);
        rpoly_init(&Q[i], k);
        built++;
    }
    if (!ok) {
        for (int i = 0; i < built; i++) rpoly_free(&D[i]);
        free(D); free(Q); rpoly_free(&r);
        return NULL;
    }

    int* m = (int*)malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
    bool reduced;
    do {
        reduced = false;
        rpoly_sort(&r, W, rows);
        for (size_t t = 0; t < r.n && !reduced; t++) {
            for (int i = 0; i < ndiv; i++) {
                if (D[i].n == 0) continue;
                if (!exp_divides(D[i].t[0].e, r.t[t].e, m, k)) continue;
                /* c = r_t.coeff / LC(D[i]) */
                Expr* c = field_div(r.t[t].c, D[i].t[0].c);
                /* quot[i] += c x^m */
                rpoly_add_term(&Q[i], m, expr_copy(c));
                /* r -= c x^m D[i] */
                int* ne = (int*)malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
                for (size_t j = 0; j < D[i].n; j++) {
                    for (int v = 0; v < k; v++) ne[v] = m[v] + D[i].t[j].e[v];
                    Expr* nc = field_neg(field_mul(c, D[i].t[j].c));
                    rpoly_add_term(&r, ne, nc);
                }
                free(ne);
                expr_free(c);
                reduced = true;
                break;
            }
        }
    } while (reduced);
    free(m);

    /* Assemble {List{a1..an}, b}. */
    Expr** qs = (Expr**)malloc(sizeof(Expr*) * (size_t)ndiv);
    for (int i = 0; i < ndiv; i++) qs[i] = rpoly_to_expr(&Q[i], mainvars);
    Expr* qlist = expr_new_function(expr_new_symbol(SYM_List), qs, (size_t)ndiv);
    free(qs);
    Expr* rem = rpoly_to_expr(&r, mainvars);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ qlist, rem }, 2);

    for (int i = 0; i < ndiv; i++) { rpoly_free(&D[i]); rpoly_free(&Q[i]); }
    free(D); free(Q); rpoly_free(&r);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Pure-Q engine (Rationals / parameter-free RationalFunctions)       */
/* ------------------------------------------------------------------ */

static Expr* pr_reduce_rational(const Expr* polyN, Expr** divsN, int ndiv,
                                Expr** mainvars, int k, const PROrder* ord) {
    /* FLINT fast path: Lexicographic only (its variable/order convention is
     * unambiguous; DegLex/DegRevLex are left to gb_divmod pending verification
     * on a FLINT-enabled build). */
    if (flint_bridge_available() && !ord->is_matrix && ord->flint_kind == 0) {
        Expr** quots = NULL; Expr* rem = NULL;
        if (flint_polynomial_reduce(polyN, (const Expr* const*)divsN, ndiv,
                                    (const Expr* const*)mainvars, k,
                                    ord->flint_kind, &quots, &rem)) {
            Expr* qlist = expr_new_function(expr_new_symbol(SYM_List),
                                            quots, (size_t)ndiv);
            free(quots);
            Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                            (Expr*[]){ qlist, rem }, 2);
            return out;
        }
    }

    /* gb_divmod engine.  Resolve the order to (use_order, weight matrix). */
    GBOrder use_order;
    int64_t* wmat_buf = NULL;
    GBWeightMatrix wmat_storage;
    const GBWeightMatrix* wmat_ptr = NULL;
    if (ord->is_matrix) {
        use_order = GB_ORDER_MATRIX;
        wmat_buf = (int64_t*)malloc(sizeof(int64_t) * (size_t)ord->mat_rows * (size_t)k);
        memcpy(wmat_buf, ord->mat, sizeof(int64_t) * (size_t)ord->mat_rows * (size_t)k);
        wmat_storage.n_rows = ord->mat_rows; wmat_storage.n_vars = k;
        wmat_storage.w = wmat_buf; wmat_ptr = &wmat_storage;
    } else if (ord->native == GB_ORDER_LEX) {
        use_order = GB_ORDER_LEX;
    } else if (ord->native == GB_ORDER_GREVLEX) {
        use_order = GB_ORDER_GREVLEX;
    } else {
        int wr = 0;
        wmat_buf = gb_build_order_matrix(k, ord->ds, ord->rev, ord->neg, &wr);
        use_order = GB_ORDER_MATRIX;
        wmat_storage.n_rows = wr; wmat_storage.n_vars = k;
        wmat_storage.w = wmat_buf; wmat_ptr = &wmat_storage;
    }

    GBPoly* P = gb_from_expr((Expr*)polyN, mainvars, k, use_order, 0, wmat_ptr);
    if (!P) { free(wmat_buf); return NULL; }

    GBPoly** basis = (GBPoly**)malloc(sizeof(GBPoly*) * (size_t)ndiv);
    bool ok = true; int nb = 0;
    for (int i = 0; i < ndiv; i++) {
        basis[i] = gb_from_expr(divsN[i], mainvars, k, use_order, 0, wmat_ptr);
        if (!basis[i]) { ok = false; break; }
        nb++;
    }
    if (!ok) {
        for (int i = 0; i < nb; i++) gb_poly_free(basis[i]);
        free(basis); gb_poly_free(P); free(wmat_buf);
        return NULL;
    }

    GBPoly** quot = NULL;
    GBPoly* rem = gb_divmod(P, basis, (size_t)ndiv, &quot);

    Expr** qs = (Expr**)malloc(sizeof(Expr*) * (size_t)ndiv);
    for (int i = 0; i < ndiv; i++) qs[i] = gb_to_expr(quot[i], mainvars);
    Expr* qlist = expr_new_function(expr_new_symbol(SYM_List), qs, (size_t)ndiv);
    free(qs);
    Expr* remE = gb_to_expr(rem, mainvars);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ qlist, remE }, 2);

    for (int i = 0; i < ndiv; i++) { gb_poly_free(basis[i]); gb_poly_free(quot[i]); }
    free(basis); free(quot); gb_poly_free(P); gb_poly_free(rem); free(wmat_buf);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Modulus engine (GF(p))                                             */
/* ------------------------------------------------------------------ */

static Expr* pr_reduce_modular(const Expr* polyN, Expr** divsN, int ndiv,
                               Expr** mainvars, int k, const PROrder* ord,
                               int64_t modulus) {
    GBOrder morder = (!ord->is_matrix && ord->native == GB_ORDER_GREVLEX)
                     ? GB_ORDER_GREVLEX : GB_ORDER_LEX;
    if (ord->is_matrix || (ord->native != GB_ORDER_LEX
                           && ord->native != GB_ORDER_GREVLEX)) {
        pr_warn("nimpl", "Modulus supports only Lexicographic and "
                         "DegreeReverseLexicographic; using Lexicographic.");
        morder = GB_ORDER_LEX;
    }

    GBPoly* Pq = gb_from_expr((Expr*)polyN, mainvars, k, morder, 0, NULL);
    if (!Pq) return NULL;
    GFpPoly* P = gfp_from_gbpoly(Pq, morder, (uint64_t)modulus);
    gb_poly_free(Pq);
    if (!P) { pr_warn("modpole", "a coefficient has no image in GF(p)."); return NULL; }

    GFpPoly** basis = (GFpPoly**)malloc(sizeof(GFpPoly*) * (size_t)ndiv);
    bool ok = true; int nb = 0;
    for (int i = 0; i < ndiv; i++) {
        GBPoly* dq = gb_from_expr(divsN[i], mainvars, k, morder, 0, NULL);
        if (!dq) { ok = false; break; }
        basis[i] = gfp_from_gbpoly(dq, morder, (uint64_t)modulus);
        gb_poly_free(dq);
        if (!basis[i]) { ok = false; break; }
        nb++;
    }
    if (!ok) {
        for (int i = 0; i < nb; i++) gfp_poly_free(basis[i]);
        free(basis); gfp_poly_free(P);
        return NULL;
    }

    GFpPoly** quot = NULL;
    GFpPoly* rem = gfp_divmod(P, basis, (size_t)ndiv, &quot);

    Expr** qs = (Expr**)malloc(sizeof(Expr*) * (size_t)ndiv);
    for (int i = 0; i < ndiv; i++) {
        GBPoly* qg = gbpoly_from_gfp(quot[i]);
        qs[i] = gb_to_expr(qg, mainvars);
        gb_poly_free(qg);
    }
    Expr* qlist = expr_new_function(expr_new_symbol(SYM_List), qs, (size_t)ndiv);
    free(qs);
    GBPoly* rg = gbpoly_from_gfp(rem);
    Expr* remE = gb_to_expr(rg, mainvars);
    gb_poly_free(rg);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ qlist, remE }, 2);

    for (int i = 0; i < ndiv; i++) { gfp_poly_free(basis[i]); gfp_poly_free(quot[i]); }
    free(basis); free(quot); gfp_poly_free(P); gfp_poly_free(rem);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Builtin entry                                                      */
/* ------------------------------------------------------------------ */

Expr* builtin_polynomialreduce(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t n_pos;
    PROptions opt;
    pr_extract_options(res, &n_pos, &opt);
    if (n_pos < 2 || n_pos > 3) return NULL;

    Expr* poly      = res->data.function.args[0];
    Expr* divisors  = res->data.function.args[1];
    Expr* vars_arg  = (n_pos == 3) ? res->data.function.args[2] : NULL;

    if (!pr_is_list(divisors)) return NULL;
    int ndiv = (int)divisors->data.function.arg_count;
    if (ndiv == 0) return NULL;

    /* Deferred domains / options: decline (leave unevaluated) with a note. */
    if (opt.domain == PR_DOM_DEFERRED || opt.tolerance_nonzero) {
        pr_warn("nimpl", "CoefficientDomain -> Integers / InexactNumbers and "
                         "nonzero Tolerance are not yet supported.");
        return NULL;
    }

    /* ---- Resolve main variables and parameters. ---- */
    Expr* allvars = pr_all_variables(poly, divisors);   /* owned List of symbols */
    if (!pr_is_list(allvars)) { expr_free(allvars); return NULL; }
    int nall = (int)allvars->data.function.arg_count;

    /* main-variable array (borrowed pointers into allvars or vars_arg). */
    Expr** mainvars = NULL; int k = 0;
    Expr* wrap = NULL;   /* synthetic single-var list to free later */
    if (vars_arg) {
        if (pr_is_list(vars_arg)) {
            k = (int)vars_arg->data.function.arg_count;
            mainvars = (Expr**)malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
            for (int i = 0; i < k; i++) mainvars[i] = vars_arg->data.function.args[i];
        } else if (vars_arg->type == EXPR_SYMBOL) {
            k = 1;
            mainvars = (Expr**)malloc(sizeof(Expr*));
            mainvars[0] = vars_arg;
        } else {
            expr_free(allvars); return NULL;
        }
    } else if (opt.parameter_vars) {
        /* vars omitted, ParameterVariables given: mainvars = allvars \ params. */
        Expr** params = NULL; int np = 0;
        if (pr_is_list(opt.parameter_vars)) {
            np = (int)opt.parameter_vars->data.function.arg_count;
            params = (Expr**)malloc(sizeof(Expr*) * (size_t)(np > 0 ? np : 1));
            for (int i = 0; i < np; i++) params[i] = opt.parameter_vars->data.function.args[i];
        } else if (opt.parameter_vars->type == EXPR_SYMBOL) {
            np = 1; params = (Expr**)malloc(sizeof(Expr*)); params[0] = (Expr*)opt.parameter_vars;
        }
        mainvars = (Expr**)malloc(sizeof(Expr*) * (size_t)(nall > 0 ? nall : 1));
        for (int i = 0; i < nall; i++) {
            Expr* v = allvars->data.function.args[i];
            if (!sym_in_list(v, params, np)) mainvars[k++] = v;
        }
        free(params);
    } else {
        /* vars omitted: all variables are main variables. */
        k = nall;
        mainvars = (Expr**)malloc(sizeof(Expr*) * (size_t)(nall > 0 ? nall : 1));
        for (int i = 0; i < nall; i++) mainvars[i] = allvars->data.function.args[i];
    }
    if (k == 0) { free(mainvars); expr_free(allvars); if (wrap) expr_free(wrap); return NULL; }

    /* Parameters = free symbols not among the main variables. */
    bool has_params = false;
    for (int i = 0; i < nall; i++) {
        if (!sym_in_list(allvars->data.function.args[i], mainvars, k)) { has_params = true; break; }
    }

    /* ---- Resolve the monomial order. ---- */
    PROrder ord;
    pr_resolve_order(opt.monomial_order, k, &ord);
    if (!ord.ok) {
        pr_warn("nimpl", "unrecognised MonomialOrder; leaving unevaluated.");
        free(mainvars); expr_free(allvars); if (wrap) expr_free(wrap);
        return NULL;
    }

    /* ---- Normalise inputs (expand / rationalise / evaluate). ---- */
    Expr* polyN = pr_normalise(poly);
    Expr** divsN = (Expr**)malloc(sizeof(Expr*) * (size_t)ndiv);
    for (int i = 0; i < ndiv; i++) divsN[i] = pr_normalise(divisors->data.function.args[i]);

    /* ---- Dispatch. ---- */
    Expr* out = NULL;
    if (opt.modulus != 0) {
        if (has_params) {
            pr_warn("nimpl", "Modulus is not supported with parameter variables.");
        } else {
            out = pr_reduce_modular(polyN, divsN, ndiv, mainvars, k, &ord, opt.modulus);
        }
    } else if (has_params) {
        if (opt.domain == PR_DOM_RATIONALS) {
            pr_warn("nimpl", "CoefficientDomain -> Rationals with a free "
                             "parameter is not supported (use the default "
                             "RationalFunctions).");
        } else {
            const int64_t* W; int rows; int64_t* Wbuf = NULL;
            if (ord.is_matrix) { W = ord.mat; rows = ord.mat_rows; }
            else { Wbuf = gb_build_order_matrix(k, ord.ds, ord.rev, ord.neg, &rows); W = Wbuf; }
            out = pr_reduce_field(polyN, divsN, ndiv, mainvars, k, W, rows);
            free(Wbuf);
        }
    } else {
        out = pr_reduce_rational(polyN, divsN, ndiv, mainvars, k, &ord);
    }

    /* ---- Cleanup. ---- */
    expr_free(polyN);
    for (int i = 0; i < ndiv; i++) expr_free(divsN[i]);
    free(divsN);
    free(ord.mat);
    free(mainvars);
    expr_free(allvars);
    if (wrap) expr_free(wrap);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Initialisation                                                     */
/* ------------------------------------------------------------------ */

void polynomialreduce_init(void) {
    symtab_add_builtin("PolynomialReduce", builtin_polynomialreduce);
    SymbolDef* def = symtab_get_def("PolynomialReduce");
    if (def) def->attributes |= ATTR_PROTECTED;
}
