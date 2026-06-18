/*
 * product.c -- Product dispatcher for Mathilda.
 *
 * The multiplicative analogue of Sum (src/sum/sum.c).  Product is HoldAll: the
 * product variable and bounds must be held so that the iterator is not
 * prematurely evaluated against an outer binding (exactly as Sum/Table/Do hold
 * their iterator specs).
 *
 * Responsibilities of this file (Stage 0):
 *   - strip trailing options (Method -> "...", VerifyConvergence -> ..., etc.);
 *   - rewrite multiple iterators Product[f, s1, ..., sk] into nested single-spec
 *     products (outer-depends-on-inner bounds come for free);
 *   - finite explicit expansion: when a range resolves to a finite span of
 *     integers, or the spec iterates an explicit list, bind the variable and
 *     fold the evaluated terms with Times (an empty product is 1);
 *   - otherwise (symbolic bounds, Infinity, or the indefinite form Product[f,i])
 *     run a Method cascade over the context-qualified sub-algorithms
 *     Product`Telescoping, Product`Rational, Product`Geometric, Product`QProduct.
 *     Each sub-builtin returns the closed form (definite:
 *     Product`M[f,i,imin,imax]; indefinite: Product`M[f,i]) or comes back
 *     unevaluated to signal "fall through".  When all stages fall through the
 *     Product[...] is returned unevaluated (held).
 *
 * Adding a later stage is purely additive: a new src/product/product_*.c file,
 * one try_* line in the cascade, and one *_init() call in product_init().
 *
 * Memory contract: builtin_product takes ownership of res but must not free it
 * (the evaluator owns it).  Every Expr* allocated here is freed on all paths.
 */

#include "product.h"
#include "product_internal.h"   /* g_product_verify_convergence */
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "iter.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational */
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/* Upper bound on the number of factors a finite explicit expansion will
 * multiply before giving up (a runaway guard, not a correctness limit). */
#define PRODUCT_MAX_FINITE_TERMS 100000000LL

/* ------------------------------------------------------------------ */
/*  Small helpers                                                      */
/* ------------------------------------------------------------------ */

/* True if e is a function with the given symbol head. */
static bool head_is_name(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* A result "came back unevaluated" if its head is still the stage we called
 * (no rule fired). */
static bool result_is_unresolved(const Expr* result, const char* head_name) {
    return head_is_name(result, head_name);
}

/* True if the step is the integer 1.  The closed-form stages take no step
 * argument and assume a unit step, so they may only be consulted when di == 1;
 * a non-unit step (e.g. {i, 1, n, 2}) must not be dispatched to them or it
 * would silently produce the wrong, step-1 closed form. */
static bool is_unit_step(const Expr* di) {
    return di && di->type == EXPR_INTEGER && di->data.integer == 1;
}

/* ------------------------------------------------------------------ */
/*  Method option                                                      */
/* ------------------------------------------------------------------ */

typedef enum {
    PROD_METHOD_AUTOMATIC = 0,
    PROD_METHOD_TELESCOPING,
    PROD_METHOD_RATIONAL,      /* aka "Hypergeometric"            */
    PROD_METHOD_GEOMETRIC,     /* aka "PolynomialExponential"     */
    PROD_METHOD_QPRODUCT,
    PROD_METHOD_INVALID
} ProdMethod;

/* Is arg an option of the form sym -> rhs (Rule or RuleDelayed)? */
static bool is_option(const Expr* arg) {
    if (arg->type != EXPR_FUNCTION) return false;
    if (arg->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = arg->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (arg->data.function.arg_count != 2) return false;
    return arg->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Map a Method -> "Name" option to a ProdMethod.  Returns PROD_METHOD_AUTOMATIC
 * for any non-Method option (those are ignored in this pass). */
static ProdMethod parse_method_option(const Expr* opt) {
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->data.symbol != SYM_Method) return PROD_METHOD_AUTOMATIC;
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)
        return PROD_METHOD_AUTOMATIC;
    if (rhs->type != EXPR_STRING) return PROD_METHOD_INVALID;
    const char* m = rhs->data.string;
    if (strcmp(m, "Automatic") == 0)              return PROD_METHOD_AUTOMATIC;
    if (strcmp(m, "Telescoping") == 0)            return PROD_METHOD_TELESCOPING;
    if (strcmp(m, "Rational") == 0)               return PROD_METHOD_RATIONAL;
    if (strcmp(m, "Hypergeometric") == 0)         return PROD_METHOD_RATIONAL;
    if (strcmp(m, "Geometric") == 0)              return PROD_METHOD_GEOMETRIC;
    if (strcmp(m, "PolynomialExponential") == 0)  return PROD_METHOD_GEOMETRIC;
    if (strcmp(m, "QProduct") == 0)               return PROD_METHOD_QPRODUCT;
    return PROD_METHOD_INVALID;
}

/* ------------------------------------------------------------------ */
/*  Stage calls                                                        */
/* ------------------------------------------------------------------ */

/* Build head[f, i, imin, imax] (definite) and evaluate it. */
static Expr* call_stage_def(const char* head, Expr* f, Expr* var,
                            Expr* imin, Expr* imax) {
    Expr* args[4] = { expr_copy(f), expr_copy(var),
                      expr_copy(imin), expr_copy(imax) };
    Expr* call = expr_new_function(expr_new_symbol(head), args, 4);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Build head[f, i] (indefinite) and evaluate it. */
static Expr* call_stage_indef(const char* head, Expr* f, Expr* var) {
    Expr* args[2] = { expr_copy(f), expr_copy(var) };
    Expr* call = expr_new_function(expr_new_symbol(head), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Try one definite stage; NULL means "fall through". */
static Expr* try_def(const char* head, Expr* f, Expr* var,
                     Expr* imin, Expr* imax) {
    Expr* r = call_stage_def(head, f, var, imin, imax);
    if (!r) return NULL;
    if (result_is_unresolved(r, head)) { expr_free(r); return NULL; }
    return r;
}

/* Try one indefinite stage; NULL means "fall through". */
static Expr* try_indef(const char* head, Expr* f, Expr* var) {
    Expr* r = call_stage_indef(head, f, var);
    if (!r) return NULL;
    if (result_is_unresolved(r, head)) { expr_free(r); return NULL; }
    return r;
}

/* Run the Automatic cascade (or a single strict stage) for a definite product.
 *
 * Cascade order is cheapest/most-specific first so the nicest closed form wins:
 * Telescoping (Gamma-free rational) before Rational (Pochhammer/Gamma) before
 * Geometric (base^k) before QProduct (q-rational). */
static Expr* dispatch_def(ProdMethod method, Expr* f, Expr* var,
                          Expr* imin, Expr* imax) {
    switch (method) {
        case PROD_METHOD_AUTOMATIC: {
            Expr* r = try_def("Product`Telescoping", f, var, imin, imax);
            if (!r) r = try_def("Product`Rational",  f, var, imin, imax);
            if (!r) r = try_def("Product`Geometric", f, var, imin, imax);
            if (!r) r = try_def("Product`QProduct",  f, var, imin, imax);
            if (!r) r = try_def("Product`Special",   f, var, imin, imax);
            if (!r) r = try_def("Product`Infinite",  f, var, imin, imax);
            return r;
        }
        case PROD_METHOD_TELESCOPING: return try_def("Product`Telescoping", f, var, imin, imax);
        case PROD_METHOD_RATIONAL:    return try_def("Product`Rational",    f, var, imin, imax);
        case PROD_METHOD_GEOMETRIC:   return try_def("Product`Geometric",   f, var, imin, imax);
        case PROD_METHOD_QPRODUCT:    return try_def("Product`QProduct",    f, var, imin, imax);
        default:                      return NULL;
    }
}

/* Run the Automatic cascade (or a single strict stage) for an indefinite
 * product.  (Product`Infinite is definite-only, so it is not in this cascade.) */
static Expr* dispatch_indef(ProdMethod method, Expr* f, Expr* var) {
    switch (method) {
        case PROD_METHOD_AUTOMATIC: {
            Expr* r = try_indef("Product`Telescoping", f, var);
            if (!r) r = try_indef("Product`Rational",  f, var);
            if (!r) r = try_indef("Product`Geometric", f, var);
            if (!r) r = try_indef("Product`QProduct",  f, var);
            return r;
        }
        case PROD_METHOD_TELESCOPING: return try_indef("Product`Telescoping", f, var);
        case PROD_METHOD_RATIONAL:    return try_indef("Product`Rational",    f, var);
        case PROD_METHOD_GEOMETRIC:   return try_indef("Product`Geometric",   f, var);
        case PROD_METHOD_QPRODUCT:    return try_indef("Product`QProduct",    f, var);
        default:                      return NULL;
    }
}

/* ------------------------------------------------------------------ */
/*  Finite explicit expansion                                          */
/* ------------------------------------------------------------------ */

/* Build Times[terms...] from a count of owned term Expr*, evaluate, and free
 * the array.  An empty product is 1 (the multiplicative identity). */
static Expr* fold_times(Expr** terms, size_t count) {
    if (count == 0) { free(terms); return expr_new_integer(1); }
    Expr* times = expr_new_function(expr_new_symbol(SYM_Times), terms, count);
    free(terms);
    Expr* r = evaluate(times);
    expr_free(times);
    return r;
}

/* Iterate an explicit list spec, binding var to each element and collecting
 * the evaluated body.  Returns the folded Times.  Never fails. */
static Expr* expand_list(Expr* f, Expr* var, Expr* list) {
    size_t n = list->data.function.arg_count;
    Expr** terms = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t count = 0;

    Rule* saved = iter_spec_shadow(var);
    for (size_t i = 0; i < n; i++) {
        symtab_add_own_value(var->data.symbol, var, list->data.function.args[i]);
        terms[count++] = evaluate(f);
    }
    iter_spec_restore(var, saved);

    return fold_times(terms, count);
}

/*
 * Iterate a numeric arithmetic-progression range [imin : imax : di], binding
 * var to the exact value each step and collecting the evaluated body.  Returns
 * the folded Times, or NULL if the span exceeds PRODUCT_MAX_FINITE_TERMS.
 *
 * Exact value stepping mirrors Do/Sum: curr_e is advanced via Plus[curr_e, di]
 * while a synchronised double drives the termination test.
 */
static Expr* expand_range(Expr* f, Expr* var, Expr* imin, Expr* imax, Expr* di,
                          double min_val, double max_val, double di_val,
                          bool is_real) {
    (void)imax;
    size_t cap = 16, count = 0;
    Expr** terms = malloc(sizeof(Expr*) * cap);

    Rule* saved = iter_spec_shadow(var);
    double val = min_val;
    Expr* curr_e = expr_copy(imin);
    bool overflow = false;

    while ((di_val > 0 && val <= max_val + 1e-14) ||
           (di_val < 0 && val >= max_val - 1e-14)) {
        if ((int64_t)count >= PRODUCT_MAX_FINITE_TERMS) { overflow = true; break; }

        Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
        symtab_add_own_value(var->data.symbol, var, i_val);
        if (count == cap) { cap *= 2; terms = realloc(terms, sizeof(Expr*) * cap); }
        terms[count++] = evaluate(f);
        expr_free(i_val);

        /* advance curr_e and val */
        Expr* nargs[2] = { expr_copy(curr_e), expr_copy(di) };
        Expr* nexpr = expr_new_function(expr_new_symbol(SYM_Plus), nargs, 2);
        Expr* next_e = evaluate(nexpr);
        expr_free(nexpr);
        expr_free(curr_e);
        curr_e = next_e;
        if (!is_real) {
            int64_t n, d;
            if (curr_e->type == EXPR_INTEGER)      val = (double)curr_e->data.integer;
            else if (curr_e->type == EXPR_REAL)    val = curr_e->data.real;
            else if (is_rational(curr_e, &n, &d))  val = (double)n / d;
        } else {
            val += di_val;
        }
    }
    expr_free(curr_e);
    iter_spec_restore(var, saved);

    if (overflow) {
        for (size_t i = 0; i < count; i++) expr_free(terms[i]);
        free(terms);
        return NULL;
    }
    return fold_times(terms, count);
}

/* ------------------------------------------------------------------ */
/*  Single-spec handler                                                */
/* ------------------------------------------------------------------ */

/* Handle Product[f, spec] for one already-isolated spec, with the parsed
 * method.  Returns the product value, or NULL to leave Product[...] held. */
static Expr* product_one_spec(Expr* f, Expr* spec, ProdMethod method) {
    /* Indefinite form Product[f, i]: spec is a bare symbol. */
    if (spec->type == EXPR_SYMBOL) {
        Expr* var = spec;
        return dispatch_indef(method, f, var);
    }

    IterSpec s;
    if (!iter_spec_parse(spec, &s)) return NULL;

    /* Product needs a named iterator; a bare-count {n} is not a valid spec. */
    if (s.kind == ITER_KIND_COUNT) { iter_spec_free(&s); return NULL; }

    if (s.kind == ITER_KIND_LIST) {
        Expr* r = expand_list(f, s.var, s.list);
        iter_spec_free(&s);
        return r;
    }

    /* RANGE: finite numeric span -> explicit expansion; otherwise closed form. */
    double min_val, max_val, di_val;
    bool is_real, is_inf;
    bool numeric = iter_spec_resolve_numeric(&s, /*allow_inf=*/true,
                                             &min_val, &max_val, &di_val,
                                             &is_real, &is_inf);
    if (numeric && !is_inf) {
        /* Closed form first: for a unit-step, non-empty, integer-bounded range
         * a rational/geometric body collapses to g(imax+1)/g(imin), whose cost
         * is independent of the span width -- far cheaper than evaluating the
         * body once per factor.  Guards (all required):
         *   - !is_real:  integer bounds/step, no float edge cases;
         *   - di_val==1: the unit-step anti-quotient is invalid for di != 1;
         *   - min<=max:  empty ranges must fold to 1 via expansion, not the
         *                telescoping form (which would give a wrong value).
         * The iterator is shadowed because Product is HoldAll. */
        if (!is_real && di_val == 1.0 && min_val <= max_val) {
            Rule* saved = iter_spec_shadow(s.var);
            Expr* cf = dispatch_def(method, f, s.var, s.imin, s.imax);
            iter_spec_restore(s.var, saved);
            if (cf) { iter_spec_free(&s); return cf; }
        }
        Expr* r = expand_range(f, s.var, s.imin, s.imax, s.di,
                               min_val, max_val, di_val, is_real);
        if (r) { iter_spec_free(&s); return r; }
        /* span too large: fall through to closed form */
    }

    /* Symbolic / Infinity / over-wide ranges: only the closed-form stages
     * remain, and they assume a unit step.  A non-unit step has no step-aware
     * closed form here, so leave the Product held rather than emit the wrong
     * step-1 result. */
    if (!is_unit_step(s.di)) { iter_spec_free(&s); return NULL; }

    Expr* r = dispatch_def(method, f, s.var, s.imin, s.imax);
    iter_spec_free(&s);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_product(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2)
        return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    /* Split trailing options from the leading body+specs. */
    ProdMethod method = PROD_METHOD_AUTOMATIC;
    int verify = 1;                      /* VerifyConvergence (default True) */
    size_t nspecs_end = argc;            /* one past the last non-option arg */
    for (size_t i = argc; i > 1; i--) {
        if (is_option(a[i - 1])) {
            ProdMethod m = parse_method_option(a[i - 1]);
            if (m == PROD_METHOD_INVALID) return NULL;
            if (m != PROD_METHOD_AUTOMATIC) method = m;
            /* VerifyConvergence -> False disables the infinite-product gate. */
            Expr* opt = a[i - 1];
            if (opt->data.function.args[0]->data.symbol == SYM_VerifyConvergence
                    && opt->data.function.args[1]->type == EXPR_SYMBOL
                    && opt->data.function.args[1]->data.symbol == SYM_False)
                verify = 0;
            nspecs_end = i - 1;
        } else {
            break;
        }
    }
    /* a[0] = body; a[1 .. nspecs_end-1] = iterator specs. */
    if (nspecs_end < 2) return NULL;
    size_t nspecs = nspecs_end - 1;

    Expr* f = a[0];

    /* Apply VerifyConvergence for the duration of this evaluation; nested
     * Product re-evaluations save/restore their own value. */
    int saved_vc = g_product_verify_convergence;
    g_product_verify_convergence = verify;

    /* Multiple iterators: rewrite Product[f, s1, ..., sk, opts] as
     *   Product[Product[f, sk, opts], s1, ..., s_{k-1}, opts]
     * and re-evaluate.  The inner product is fully evaluated under each binding
     * of the outer variables, so outer bounds may depend on inner variables. */
    if (nspecs > 1) {
        size_t nopts = argc - nspecs_end;

        /* inner = Product[f, sk, opts...] */
        Expr** in = malloc(sizeof(Expr*) * (2 + nopts));
        in[0] = expr_copy(f);
        in[1] = expr_copy(a[nspecs_end - 1]);
        for (size_t k = 0; k < nopts; k++) in[2 + k] = expr_copy(a[nspecs_end + k]);
        Expr* inner = expr_new_function(expr_new_symbol(SYM_Product), in, 2 + nopts);
        free(in);

        /* outer = Product[inner, s1, ..., s_{k-1}, opts...] */
        size_t outc = 1 + (nspecs - 1) + nopts;
        Expr** out = malloc(sizeof(Expr*) * outc);
        size_t w = 0;
        out[w++] = inner;
        for (size_t i = 1; i < nspecs; i++) out[w++] = expr_copy(a[i]);
        for (size_t k = 0; k < nopts; k++) out[w++] = expr_copy(a[nspecs_end + k]);
        Expr* outer = expr_new_function(expr_new_symbol(SYM_Product), out, outc);
        free(out);

        Expr* r = evaluate(outer);
        expr_free(outer);
        g_product_verify_convergence = saved_vc;
        return r;
    }

    Expr* r = product_one_spec(f, a[1], method);
    g_product_verify_convergence = saved_vc;
    return r;
}

void product_init(void) {
    symtab_add_builtin("Product", builtin_product);
    symtab_get_def("Product")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init), per the contributor brief. */

    /* Sub-algorithm stages (registered as Product`Telescoping,
     * Product`Rational, Product`Geometric, ...).  Added incrementally; the
     * cascade tolerates absent stages (Product`QProduct and Product`Infinite
     * are wired in their own phases and simply fall through until then). */
    void product_telescoping_init(void);
    product_telescoping_init();
    void product_rational_init(void);
    product_rational_init();
    void product_geometric_init(void);
    product_geometric_init();
    void product_qproduct_init(void);
    product_qproduct_init();
    void product_special_init(void);
    product_special_init();
    void product_infinite_init(void);
    product_infinite_init();
}
