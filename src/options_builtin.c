/* options_builtin.c — Options, SetOptions, OptionValue and the registry of
 * default option settings for option-accepting builtins.
 *
 * Mathilda stores a symbol's default options as a List[Rule[name, val], ...]
 * on SymbolDef.default_options (the DefaultValues-equivalent), reached through
 * symtab_set_options / symtab_get_options. This file implements:
 *
 *   Options[...]      query a symbol's defaults or an expression's explicit
 *                     options, optionally selected by name.
 *   SetOptions[...]   redefine individual default options of a symbol.
 *   OptionValue[...]  resolve a single option value from explicit options plus
 *                     defaults; the bare/2-arg forms are resolved inside a rule
 *                     by optionvalue_inject_context (see apply_down_values).
 *   options_register_defaults  the comprehensive table wired from core_init.
 *
 * Memory: every result is freshly built. Sub-expressions taken from `res` or
 * from the stored options are duplicated with expr_copy (a refcount bump), so
 * nothing aliased is mutated in place and the evaluator remains free to release
 * `res` after a builtin returns.
 */
#include "options.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "eval.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---------- small option helpers ---------- */

/* Textual name of an option key (symbol or string), or NULL. */
static const char* opt_name_text(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name;
    if (e->type == EXPR_STRING) return e->data.string;
    return NULL;
}

/* Context-insensitive comparison: drop any `Context`` prefix so a symbol s and
 * the string "SymbolName[s]" compare equal, matching OptionValue semantics. */
static const char* strip_context(const char* s) {
    const char* bt = strrchr(s, '`');
    return bt ? bt + 1 : s;
}
static bool opt_name_eq(const Expr* a, const Expr* b) {
    const char* sa = opt_name_text(a);
    const char* sb = opt_name_text(b);
    if (!sa || !sb) return false;
    return strcmp(strip_context(sa), strip_context(sb)) == 0;
}

/* True for Rule[name, val] / RuleDelayed[name, val] with a symbol/string name. */
static bool is_option_rule(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    return lhs && (lhs->type == EXPR_SYMBOL || lhs->type == EXPR_STRING);
}

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

static Expr* empty_list(void) {
    return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
}

/* Build List[items...]; adopts the element pointers, copies the array. */
static Expr* list_from(Expr** items, size_t n) {
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), items, n);
    return l;
}

/* First matching rule's RHS (copied), else NULL. */
static Expr* lookup_value(const Expr* list, const Expr* name) {
    if (!is_list(list)) return NULL;
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        Expr* r = list->data.function.args[i];
        if (is_option_rule(r) && opt_name_eq(r->data.function.args[0], name))
            return expr_copy(r->data.function.args[1]);
    }
    return NULL;
}

/* First matching whole rule (copied), else NULL. */
static Expr* lookup_rule(const Expr* list, const Expr* name) {
    if (!is_list(list)) return NULL;
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        Expr* r = list->data.function.args[i];
        if (is_option_rule(r) && opt_name_eq(r->data.function.args[0], name))
            return expr_copy(r);
    }
    return NULL;
}

/* The option list for an object as a NEW List (never NULL):
 *   symbol   -> a copy of its registered defaults (or {});
 *   compound -> the option rules explicitly present in its arguments;
 *   other    -> {}. */
static Expr* options_of_object(const Expr* obj) {
    if (!obj) return empty_list();
    if (obj->type == EXPR_SYMBOL) {
        Expr* stored = symtab_get_options(obj->data.symbol.name);
        return is_list(stored) ? expr_copy(stored) : empty_list();
    }
    if (obj->type == EXPR_FUNCTION) {
        size_t n = obj->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* a = obj->data.function.args[i];
            if (is_option_rule(a)) items[k++] = expr_copy(a);
        }
        Expr* l = list_from(items, k);
        free(items);
        return l;
    }
    return empty_list();
}

/* ---------- Options ---------- */

Expr* builtin_options(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n < 1 || n > 2) return NULL;

    Expr* obj = res->data.function.args[0];
    Expr* full = options_of_object(obj);
    if (n == 1) return full;

    /* Options[obj, name] / Options[obj, {names}] -> list of selected rules. */
    Expr* sel = res->data.function.args[1];
    Expr** items;
    size_t k = 0;
    if (is_list(sel)) {
        size_t m = sel->data.function.arg_count;
        items = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* r = lookup_rule(full, sel->data.function.args[i]);
            if (r) items[k++] = r;
        }
    } else {
        items = malloc(sizeof(Expr*) * 1);
        Expr* r = lookup_rule(full, sel);
        if (r) items[k++] = r;
    }
    expr_free(full);
    Expr* out = list_from(items, k);
    free(items);
    return out;
}

/* ---------- SetOptions ---------- */

Expr* builtin_setoptions(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n < 1) return NULL;

    Expr* sobj = res->data.function.args[0];
    if (sobj->type != EXPR_SYMBOL) return NULL;
    const char* sym = sobj->data.symbol.name;

    if (get_attributes(sym) & ATTR_LOCKED) {
        fprintf(stderr, "SetOptions::locked: Symbol %s is locked and cannot be modified.\n", sym);
        return NULL;
    }

    /* Working copy of the current options as a flat vector of rule copies. */
    Expr* cur = symtab_get_options(sym);
    size_t cn = is_list(cur) ? cur->data.function.arg_count : 0;
    Expr** work = malloc(sizeof(Expr*) * (cn + n)); /* cn existing + up to n-1 new */
    size_t wn = 0;
    for (size_t i = 0; i < cn; i++) work[wn++] = expr_copy(cur->data.function.args[i]);

    for (size_t i = 1; i < n; i++) {
        Expr* r = res->data.function.args[i];
        if (!is_option_rule(r)) {                 /* malformed argument */
            for (size_t j = 0; j < wn; j++) expr_free(work[j]);
            free(work);
            return NULL;
        }
        Expr* nm = r->data.function.args[0];
        size_t idx = (size_t)-1;
        for (size_t j = 0; j < wn; j++) {
            if (opt_name_eq(work[j]->data.function.args[0], nm)) { idx = j; break; }
        }
        if (idx == (size_t)-1) {                  /* SetOptions can't add options */
            const char* nms = opt_name_text(nm);
            fprintf(stderr, "SetOptions::optnf: %s is not a known option for %s.\n",
                    nms ? nms : "?", sym);
            for (size_t j = 0; j < wn; j++) expr_free(work[j]);
            free(work);
            return NULL;
        }
        expr_free(work[idx]);
        work[idx] = expr_copy(r);                 /* preserve position, take new value */
    }

    Expr* newlist = list_from(work, wn);
    free(work);
    symtab_set_options(sym, newlist);             /* store (ownership transferred) */
    eval_clock_bump();                            /* options changed: drop eval cache */
    return expr_copy(symtab_get_options(sym));    /* return the new Options[s] */
}

/* ---------- OptionValue ---------- */

/* Look up `name` in the defaults derived from `f`, where f may be Automatic
 * (no defaults), a symbol (-> its Options), a single rule, or a list of
 * symbols/rules. First specification wins. Returns a copied RHS, or NULL. */
static Expr* defaults_lookup(const Expr* f, const Expr* name) {
    if (!f) return NULL;
    if (f->type == EXPR_SYMBOL) {
        if (f->data.symbol.name == SYM_Automatic) return NULL;
        return lookup_value(symtab_get_options(f->data.symbol.name), name);
    }
    if (is_option_rule(f)) {
        if (opt_name_eq(f->data.function.args[0], name))
            return expr_copy(f->data.function.args[1]);
        return NULL;
    }
    if (is_list(f)) {
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            Expr* v = defaults_lookup(f->data.function.args[i], name);
            if (v) return v;
        }
    }
    return NULL;
}

/* Resolve `name`: explicit options first, then defaults derived from f. */
static Expr* ov_resolve(const Expr* f, const Expr* opts, const Expr* name) {
    if (opts) {
        if (is_list(opts)) {
            Expr* v = lookup_value(opts, name);
            if (v) return v;
        } else if (is_option_rule(opts)) {
            if (opt_name_eq(opts->data.function.args[0], name))
                return expr_copy(opts->data.function.args[1]);
        }
    }
    return defaults_lookup(f, name);
}

Expr* builtin_optionvalue(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    const Expr* f = NULL;
    const Expr* opts = NULL;
    const Expr* name = NULL;
    bool hold = false;

    if (n == 1) {
        /* OptionValue[name] only resolves inside a rule (handled by the
         * inject pass). At top level there is no context -> leave it alone. */
        return NULL;
    } else if (n == 2) {
        f = a[0]; name = a[1];
    } else if (n == 3) {
        f = a[0]; opts = a[1]; name = a[2];
    } else if (n == 4) {
        f = a[0]; opts = a[1]; name = a[2];
        if (!(a[3]->type == EXPR_SYMBOL && a[3]->data.symbol.name == SYM_Hold)) return NULL;
        hold = true;
    } else {
        return NULL;
    }

    if (!(name->type == EXPR_SYMBOL || name->type == EXPR_STRING)) return NULL;

    Expr* val = ov_resolve(f, opts, name);
    if (!val) return NULL;                         /* unknown option: leave unevaluated */

    if (hold) {
        Expr** h = malloc(sizeof(Expr*));
        h[0] = val;
        Expr* r = expr_new_function(expr_new_symbol(SYM_Hold), h, 1);
        free(h);
        return r;
    }
    return val;
}

/* Functional rewrite of OptionValue nodes inside a fired rule's RHS. Returns a
 * NEW tree; the input is never mutated (replace_bindings output may share
 * refcounted nodes with the stored rule template). */
Expr* optionvalue_inject_context(const Expr* e, const char* head_sym,
                                 const Expr* opts) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;

    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_OptionValue
        && (n == 1 || n == 2)) {
        Expr** na = malloc(sizeof(Expr*) * 3);
        if (n == 1) {                              /* [name] -> [head, opts, name] */
            na[0] = expr_new_symbol(head_sym);
            na[1] = expr_copy((Expr*)opts);
            na[2] = expr_copy(e->data.function.args[0]);
        } else {                                   /* [f, name] -> [f, opts, name] */
            na[0] = expr_copy(e->data.function.args[0]);
            na[1] = expr_copy((Expr*)opts);
            na[2] = expr_copy(e->data.function.args[1]);
        }
        Expr* r = expr_new_function(expr_new_symbol(SYM_OptionValue), na, 3);
        free(na);
        return r;
    }

    /* Rebuild this node, recursing into head and every argument. */
    Expr* new_head = optionvalue_inject_context(h, head_sym, opts);
    Expr** na = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++)
        na[i] = optionvalue_inject_context(e->data.function.args[i], head_sym, opts);
    Expr* r = expr_new_function(new_head, na, n);
    if (na) free(na);
    return r;
}

/* ====================================================================
 * Default-options registry (comprehensive sweep).
 *
 * One readable table mirroring the hardcoded C defaults of every
 * option-accepting builtin. Values are written in their Wolfram surface
 * form (sentinels like -1 become Automatic, machine precision becomes
 * MachinePrecision, ...). Keep in sync with each builtin's parser.
 * ==================================================================== */

typedef struct { Expr** items; size_t n, cap; } OptBuf;
static void ob_init(OptBuf* b) { b->cap = 8; b->n = 0; b->items = malloc(b->cap * sizeof(Expr*)); }
static void ob_add(OptBuf* b, Expr* rule) {
    if (b->n == b->cap) { b->cap *= 2; b->items = realloc(b->items, b->cap * sizeof(Expr*)); }
    b->items[b->n++] = rule;
}
static void ob_commit(OptBuf* b, const char* fname) {
    Expr* l = list_from(b->items, b->n);
    free(b->items);
    symtab_set_options(fname, l);
}

static Expr* rule2(Expr* name, Expr* val) {
    Expr** ra = malloc(2 * sizeof(Expr*));
    ra[0] = name; ra[1] = val;
    Expr* r = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    free(ra);
    return r;
}
static Expr* r_sym(const char* nm, const char* val) { return rule2(expr_new_symbol(nm), expr_new_symbol(val)); }
static Expr* r_int(const char* nm, long v)          { return rule2(expr_new_symbol(nm), expr_new_integer(v)); }
static Expr* r_real(const char* nm, double v)       { return rule2(expr_new_symbol(nm), expr_new_real(v)); }
static Expr* r_str(const char* nm, const char* val) { return rule2(expr_new_symbol(nm), expr_new_string(val)); }
/* name -> {} (e.g. D's NonConstants default). */
static Expr* r_list0(const char* nm) { return rule2(expr_new_symbol(nm), empty_list()); }

void options_register_defaults(void) {
    OptBuf b;

    /* ---- Numerical calculus ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_sym("AccuracyGoal", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_sym("MaxRecursion", "Automatic"));
    ob_add(&b, r_int("MinRecursion", 0));
    ob_add(&b, r_sym("MaxPoints", "Automatic"));
    ob_add(&b, r_sym("Exclusions", "None"));
    ob_commit(&b, "NIntegrate");

    /* ---- Number theory ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "PrimePi");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_int("NSumTerms", 15));
    ob_add(&b, r_sym("NSumExtraTerms", "Automatic"));
    ob_add(&b, r_int("WynnDegree", 1));
    ob_add(&b, r_sym("VerifyConvergence", "True"));
    ob_add(&b, r_sym("AccuracyGoal", "Infinity"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_commit(&b, "NSum");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_int("NProductFactors", 15));
    ob_add(&b, r_sym("NProductExtraFactors", "Automatic"));
    ob_add(&b, r_int("WynnDegree", 1));
    ob_add(&b, r_sym("VerifyConvergence", "True"));
    ob_commit(&b, "NProduct");

    ob_init(&b);
    ob_add(&b, r_real("Radius", 0.01));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_int("MaxRecursion", 10));
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "NResidue");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "EulerSum"));
    ob_add(&b, r_sym("Direction", "Automatic"));
    ob_add(&b, r_sym("Scale", "Automatic"));
    ob_add(&b, r_int("Terms", 7));
    ob_add(&b, r_int("WynnDegree", 1));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_commit(&b, "NLimit");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "EulerSum"));
    ob_add(&b, r_sym("Scale", "Automatic"));
    ob_add(&b, r_int("Terms", 7));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_commit(&b, "ND");

    ob_init(&b);
    ob_add(&b, r_real("Radius", 1.0));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_commit(&b, "NSeries");

    /* ---- Numerical root-finding / solving ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_int("MaxIterations", 100));
    ob_add(&b, r_sym("AccuracyGoal", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_real("DampingFactor", 1.0));
    ob_add(&b, r_sym("Jacobian", "Automatic"));
    ob_add(&b, r_sym("EvaluationMonitor", "None"));
    ob_add(&b, r_sym("StepMonitor", "None"));
    ob_commit(&b, "FindRoot");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_int("MaxIterations", 500));
    ob_add(&b, r_sym("AccuracyGoal", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_sym("Gradient", "Automatic"));
    ob_add(&b, r_sym("EvaluationMonitor", "None"));
    ob_add(&b, r_sym("StepMonitor", "None"));
    ob_commit(&b, "FindMinimum");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_int("MaxIterations", 500));
    ob_add(&b, r_sym("AccuracyGoal", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_sym("Gradient", "Automatic"));
    ob_add(&b, r_sym("EvaluationMonitor", "None"));
    ob_add(&b, r_sym("StepMonitor", "None"));
    ob_commit(&b, "FindMaximum");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("MaxIterations", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_commit(&b, "NRoots");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("MaxRoots", "Automatic"));
    ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
    ob_add(&b, r_sym("VerifySolutions", "Automatic"));
    ob_add(&b, r_sym("PrecisionGoal", "Automatic"));
    ob_commit(&b, "NSolve");

    /* ---- Symbolic solving ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Cubics", "False"));
    ob_add(&b, r_sym("Quartics", "False"));
    ob_add(&b, r_sym("GeneratedParameters", "C"));
    ob_add(&b, r_sym("InverseFunctions", "Automatic"));
    ob_add(&b, r_sym("VerifySolutions", "Automatic"));
    ob_add(&b, r_sym("Assumptions", "Automatic"));
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_int("Modulus", 0));
    ob_commit(&b, "Solve");

    /* ---- Polynomial / rational ---- */
    const char* poly_fns[] = {
        "PolynomialGCD", "PolynomialLCM", "PolynomialQuotient",
        "PolynomialRemainder", "PolynomialQuotientRemainder"
    };
    for (size_t i = 0; i < sizeof(poly_fns) / sizeof(poly_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("Extension", "None"));
        ob_add(&b, r_int("Modulus", 0));
        ob_commit(&b, poly_fns[i]);
    }
    const char* rat_fns[] = { "Together", "Cancel", "Apart", "Factor" };
    for (size_t i = 0; i < sizeof(rat_fns) / sizeof(rat_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("Extension", "None"));
        ob_add(&b, r_int("Modulus", 0));
        ob_commit(&b, rat_fns[i]);
    }

    /* ---- Fitting / linear algebra ---- */
    ob_init(&b);
    ob_add(&b, r_sym("WorkingPrecision", "Automatic"));
    ob_add(&b, r_sym("FitRegularization", "None"));
    ob_add(&b, r_sym("NormFunction", "Automatic"));
    ob_commit(&b, "Fit");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_int("Modulus", 0));
    ob_add(&b, r_sym("ZeroTest", "Automatic"));
    ob_commit(&b, "LinearSolve");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("Tolerance", "Automatic"));
    ob_commit(&b, "LeastSquares");

    /* ---- Graphics ---- */
    ob_init(&b);
    ob_add(&b, r_int("PlotPoints", 50));
    ob_add(&b, r_int("MaxRecursion", 6));
    ob_add(&b, r_sym("MaxPlotPoints", "Automatic"));
    ob_add(&b, r_sym("Mesh", "False"));
    ob_add(&b, r_sym("AspectRatio", "Automatic"));
    ob_add(&b, r_sym("Axes", "True"));
    ob_add(&b, r_sym("Frame", "False"));
    ob_add(&b, r_sym("PlotRange", "Automatic"));
    ob_commit(&b, "Plot");

    ob_init(&b);
    ob_add(&b, r_sym("Joined", "False"));
    ob_add(&b, r_sym("DataRange", "Automatic"));
    ob_add(&b, r_sym("Filling", "None"));
    ob_add(&b, r_sym("FillingStyle", "Automatic"));
    ob_add(&b, r_sym("PlotMarkers", "None"));
    ob_add(&b, r_sym("PlotStyle", "Automatic"));
    ob_add(&b, r_sym("PlotLegends", "None"));
    ob_add(&b, r_sym("Axes", "True"));
    /* AspectRatio -> 1/GoldenRatio (its Wolfram surface default). */
    {
        Expr* gr[2] = { expr_new_symbol("GoldenRatio"), expr_new_integer(-1) };
        Expr* inv = expr_new_function(expr_new_symbol("Power"), gr, 2);
        ob_add(&b, rule2(expr_new_symbol("AspectRatio"), inv));
    }
    ob_add(&b, r_sym("Frame", "False"));
    ob_add(&b, r_sym("PlotRange", "Automatic"));
    ob_commit(&b, "ListPlot");

    /* Graphics[] honours these via render.c's gfx_options_parse(); the set
     * mirrors that parser exactly. Alphabetical, in their Wolfram surface
     * defaults. Prolog/Epilog/LabelStyle default to {} (nothing extra to
     * draw / no overriding style). */
    ob_init(&b);
    ob_add(&b, r_sym("AspectRatio", "Automatic"));
    ob_add(&b, r_sym("Axes", "False"));
    ob_add(&b, r_sym("AxesLabel", "None"));
    ob_add(&b, r_sym("AxesOrigin", "Automatic"));
    ob_add(&b, r_sym("AxesStyle", "Automatic"));
    ob_add(&b, r_sym("Background", "Automatic"));
    ob_add(&b, r_list0("Epilog"));
    ob_add(&b, r_sym("Frame", "False"));
    ob_add(&b, r_sym("FrameLabel", "None"));
    ob_add(&b, r_sym("FrameStyle", "Automatic"));
    ob_add(&b, r_sym("FrameTicks", "Automatic"));
    ob_add(&b, r_sym("GridLines", "None"));
    ob_add(&b, r_sym("GridLinesStyle", "Automatic"));
    ob_add(&b, r_sym("ImageSize", "Automatic"));
    ob_add(&b, r_list0("LabelStyle"));
    ob_add(&b, r_sym("PlotLabel", "None"));
    ob_add(&b, r_sym("PlotRange", "Automatic"));
    ob_add(&b, r_sym("PlotRangePadding", "Automatic"));
    ob_add(&b, r_sym("PlotStyle", "Automatic"));
    ob_add(&b, r_list0("Prolog"));
    ob_add(&b, r_sym("RotateLabel", "True"));
    ob_add(&b, r_sym("TicksStyle", "Automatic"));
    ob_commit(&b, "Graphics");

    /* ---- Symbolic calculus ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));     /* "Automatic" | "BronsteinRational" | ... */
    ob_commit(&b, "Integrate");

    ob_init(&b);
    ob_add(&b, r_sym("Direction", "Automatic"));
    ob_add(&b, r_sym("Assumptions", "Automatic"));
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "Limit");

    ob_init(&b);
    ob_add(&b, r_sym("Assumptions", "Automatic"));
    ob_commit(&b, "Series");

    ob_init(&b);
    ob_add(&b, r_list0("NonConstants"));
    ob_commit(&b, "D");

    ob_init(&b);
    ob_add(&b, r_list0("NonConstants"));
    ob_commit(&b, "Dt");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "Sum");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "Product");

    /* ---- Simplification ---- */
    const char* simp_fns[] = { "Simplify", "FullSimplify" };
    for (size_t i = 0; i < sizeof(simp_fns) / sizeof(simp_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("Assumptions", "Automatic"));
        ob_add(&b, r_sym("ComplexityFunction", "Automatic"));
        ob_add(&b, r_sym("TransformationFunctions", "Automatic"));
        ob_commit(&b, simp_fns[i]);
    }

    ob_init(&b);
    ob_add(&b, r_sym("Assumptions", "Automatic"));
    ob_commit(&b, "PowerExpand");

    /* ---- Polynomial algebra ---- */
    ob_init(&b);
    ob_add(&b, r_sym("MonomialOrder", "Lexicographic"));
    ob_add(&b, r_sym("CoefficientDomain", "Rationals"));
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("Sort", "False"));
    ob_add(&b, r_int("Modulus", 0));
    ob_commit(&b, "GroebnerBasis");

    ob_init(&b);
    ob_add(&b, r_sym("Extension", "None"));
    ob_commit(&b, "MinimalPolynomial");

    ob_init(&b);
    ob_add(&b, r_sym("Extension", "None"));
    ob_add(&b, r_sym("GaussianIntegers", "Automatic"));
    ob_add(&b, r_int("Modulus", 0));
    ob_commit(&b, "IrreduciblePolynomialQ");

    ob_init(&b);
    ob_add(&b, r_sym("GaussianIntegers", "Automatic"));
    ob_add(&b, r_int("Modulus", 0));
    ob_commit(&b, "SquareFreeQ");

    ob_init(&b);
    ob_add(&b, r_int("Modulus", 0));
    ob_commit(&b, "PolynomialExtendedGCD");

    /* ---- Linear algebra ---- */
    ob_init(&b);
    ob_add(&b, r_sym("Cubics", "True"));
    ob_add(&b, r_sym("Quartics", "True"));
    ob_commit(&b, "Eigenvalues");

    ob_init(&b);
    ob_add(&b, r_sym("Cubics", "True"));
    ob_add(&b, r_sym("Quartics", "True"));
    ob_commit(&b, "Eigenvectors");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_commit(&b, "NullSpace");

    ob_init(&b);
    ob_add(&b, r_sym("Method", "Automatic"));
    ob_add(&b, r_sym("Tolerance", "Automatic"));
    ob_commit(&b, "MatrixRank");

    ob_init(&b);
    ob_add(&b, r_sym("Tolerance", "Automatic"));
    ob_commit(&b, "PseudoInverse");

    ob_init(&b);
    ob_add(&b, r_sym("Tolerance", "Automatic"));
    ob_add(&b, r_str("TargetStructure", "Dense"));
    ob_commit(&b, "SingularValueDecomposition");

    /* ---- Number theory ---- */
    ob_init(&b);
    ob_add(&b, r_sym("GaussianIntegers", "False"));
    ob_commit(&b, "PrimeQ");

    ob_init(&b);
    ob_add(&b, r_sym("GaussianIntegers", "False"));
    ob_commit(&b, "CoprimeQ");

    ob_init(&b);
    ob_add(&b, r_sym("GaussianIntegers", "False"));
    ob_commit(&b, "FactorInteger");

    ob_init(&b);
    ob_add(&b, r_sym("GaussianIntegers", "False"));
    ob_commit(&b, "Divisors");

    /* ---- Random number generation ---- */
    const char* rand_fns[] = { "RandomInteger", "RandomReal", "RandomComplex" };
    for (size_t i = 0; i < sizeof(rand_fns) / sizeof(rand_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("WorkingPrecision", "MachinePrecision"));
        ob_commit(&b, rand_fns[i]);
    }

    /* ---- Structural functions reading Heads -> True|False ---- */
    struct { const char* fn; const char* dflt; } heads_fns[] = {
        { "Cases", "False" }, { "Count", "False" }, { "DeleteCases", "False" },
        { "MemberQ", "False" }, { "Position", "True" }, { "FreeQ", "True" },
        { "Level", "False" }, { "Depth", "False" }, { "LeafCount", "False" },
        { "Map", "False" }, { "Apply", "False" }, { "MapAll", "False" },
    };
    for (size_t i = 0; i < sizeof(heads_fns) / sizeof(heads_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("Heads", heads_fns[i].dflt));
        ob_commit(&b, heads_fns[i].fn);
    }

    /* ---- Fixed-point iteration ---- */
    const char* fp_fns[] = { "FixedPoint", "FixedPointList" };
    for (size_t i = 0; i < sizeof(fp_fns) / sizeof(fp_fns[0]); i++) {
        ob_init(&b);
        ob_add(&b, r_sym("SameTest", "Automatic"));
        ob_commit(&b, fp_fns[i]);
    }
}

void options_builtin_init(void) {
    symtab_add_builtin("Options", builtin_options);
    symtab_get_def("Options")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("SetOptions", builtin_setoptions);
    symtab_get_def("SetOptions")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("OptionValue", builtin_optionvalue);
    symtab_get_def("OptionValue")->attributes |= ATTR_PROTECTED;

    options_register_defaults();
}
