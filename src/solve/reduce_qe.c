/*
 * reduce_qe.c
 *
 * Quantifier elimination for `Reduce` (REDUCE_PLAN.md, Phase 7): the front-end
 * for the `Exists`, `ForAll` and `Resolve` heads.  See reduce_qe.h for the shape
 * of the method and the three-case (by free-variable count) routing.
 *
 * This file owns the front-end only -- quantifier normalisation (flatten a
 * same-kind chain, fold a 3-argument condition), free-variable collection, the
 * fully-quantified DECISION path (Case A, which reuses the whole Reduce engine),
 * the routing to reduce_cad_qe for the parametric (>=1 free variable) path, and
 * the recursive composition that eliminates an alternating quantifier prefix
 * inner-block-first.  The CAD projection/lifting/fold and the multi-free-variable
 * emission live in reduce_cad.c.
 *
 * Hard invariant: any decline (a malformed node, a non-Reals domain, or an
 * undecidable/unsolvable sub-problem -- including an alternating sub-block whose
 * inner elimination declines) returns NULL, leaving the input unevaluated --
 * never a wrong formula.
 */
#include "reduce_qe.h"
#include "reduce_form.h"
#include "reduce_cad.h"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "attr.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool is_quantifier(const Expr* e) {
    return is_head(e, SYM_Exists) || is_head(e, SYM_ForAll);
}

/* The known constant symbols -- excluded from the free-variable set (they are
 * numeric coefficients, not variables).  Mirrors nsolve.c's is_constant_symbol. */
static bool qe_is_constant_symbol(const char* s) {
    return s == SYM_Pi || s == SYM_E || s == SYM_I || s == SYM_Degree
        || s == SYM_EulerGamma || s == SYM_GoldenRatio || s == SYM_Catalan
        || s == SYM_Glaisher || s == SYM_Khinchin || s == SYM_Indeterminate
        || s == SYM_Infinity || s == SYM_ComplexInfinity
        || s == SYM_True || s == SYM_False;
}

/* Node builders (each CONSUMES its Expr* arguments). */
static Expr* mkfun1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mkfun2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mkfun3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}

/* Append a distinct interned name to (*arr,*n,*cap). */
static void name_push(const char*** arr, int* n, int* cap, const char* s) {
    for (int i = 0; i < *n; i++) if ((*arr)[i] == s) return;
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *arr = realloc(*arr, (size_t)*cap * sizeof(char*)); }
    (*arr)[(*n)++] = s;
}

/* Collect the distinct non-constant symbols appearing as leaves (in argument
 * position, never a head) of `e`.  Names are interned pointers, borrowed. */
static void qe_collect_symbols(const Expr* e, const char*** out, int* n, int* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (!qe_is_constant_symbol(e->data.symbol.name))
            name_push(out, n, cap, e->data.symbol.name);
        return;
    }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            qe_collect_symbols(e->data.function.args[i], out, n, cap);
}

/* Append the variable names of a bound spec (a symbol or a List of symbols) to
 * (*B,*nb,*bcap).  An empty List is accepted (adds nothing -> nbound==0 strip).
 * Returns false on a malformed spec (a non-symbol member). */
static bool qe_add_boundvars(const Expr* spec, const char*** B, int* nb, int* bcap) {
    if (spec->type == EXPR_SYMBOL) { name_push(B, nb, bcap, spec->data.symbol.name); return true; }
    if (is_head(spec, SYM_List)) {
        for (size_t i = 0; i < spec->data.function.arg_count; i++) {
            const Expr* v = spec->data.function.args[i];
            if (v->type != EXPR_SYMBOL) return false;
            name_push(B, nb, bcap, v->data.symbol.name);
        }
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Quantifier normalisation                                           *
 * ------------------------------------------------------------------ */

/* Peel a maximal chain of SAME-kind quantifiers off `q` (an Exists/ForAll expr).
 * Fills (*B,*nb,*bcap) with the bound-variable names, sets *quant (0 = Exists,
 * 1 = ForAll), and returns a freshly-owned body Expr (with any 3-argument
 * condition folded in: Exists[x,c,g] -> c && g, ForAll[x,c,g] -> !c || g).  Sets
 * *ok=false on a malformed node.  If the remaining body is itself a
 * different-kind quantifier, *alternating is set (the caller declines). */
static Expr* qe_normalize(const Expr* q, int* quant, const char*** B, int* nb,
                          int* bcap, bool* alternating, bool* ok) {
    *ok = true; *alternating = false; *quant = -1;
    const Expr* cur = q;
    Expr* body_owned = NULL;   /* set when a 3-argument fold builds a new node */

    while (cur && is_quantifier(cur)) {
        int kind = is_head(cur, SYM_Exists) ? 0 : 1;
        if (*quant == -1) *quant = kind;
        else if (kind != *quant) break;                 /* alternating: cur is body */
        size_t ac = cur->data.function.arg_count;
        if (ac != 2 && ac != 3) { *ok = false; break; }
        if (!qe_add_boundvars(cur->data.function.args[0], B, nb, bcap)) { *ok = false; break; }
        const Expr* inner = cur->data.function.args[ac - 1];
        if (ac == 3) {
            const Expr* cond = cur->data.function.args[1];
            Expr* folded = (kind == 0)
                ? mkfun2(SYM_And, expr_copy((Expr*)cond), expr_copy((Expr*)inner))
                : mkfun2(SYM_Or, mkfun1(SYM_Not, expr_copy((Expr*)cond)), expr_copy((Expr*)inner));
            if (body_owned) expr_free(body_owned);
            body_owned = folded;
            cur = body_owned;                            /* And/Or: loop ends */
        } else {
            cur = inner;
        }
    }
    if (!*ok || *quant < 0) { if (body_owned) expr_free(body_owned); return NULL; }
    if (is_quantifier(cur)) *alternating = true;         /* different-kind remnant */

    if (cur == body_owned) return body_owned;            /* transfer ownership */
    Expr* result = expr_copy((Expr*)cur);
    if (body_owned) expr_free(body_owned);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Case A -- fully-quantified decision (reuse the whole engine)        *
 * ------------------------------------------------------------------ */

/* Build and evaluate `Reduce[body, varlist, Reals]`.  `names[0..nv-1]` are the
 * variable names (a single symbol when nv==1, else a List).  Returns the
 * evaluated Expr (owned): True / False / a formula, or an unevaluated
 * `Reduce[...]` when the engine declined.  Requires nv >= 1. */
static Expr* qe_call_reduce(const Expr* body, const char** names, int nv) {
    Expr* vlist;
    if (nv == 1) {
        vlist = expr_new_symbol(names[0]);
    } else {
        Expr** vs = malloc((size_t)nv * sizeof(Expr*));
        for (int i = 0; i < nv; i++) vs[i] = expr_new_symbol(names[i]);
        vlist = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)nv);
        free(vs);
    }
    Expr* call = mkfun3(SYM_Reduce, expr_copy((Expr*)body), vlist, expr_new_symbol(SYM_Reals));
    return eval_and_free(call);
}

/* Fully-quantified decision (nfree==0).  Exists[{B},g] is True unless the
 * solution set of g over the reals is empty; ForAll[{B},g] is True only when it
 * is all of R^|B|.  Declines (NULL) when Reduce leaves the sub-problem
 * unevaluated. */
static Expr* qe_decide(const Expr* body, int quant, const char** B, int nb) {
    Expr* r = qe_call_reduce(body, B, nb);
    if (!r) return NULL;
    if (is_head(r, SYM_Reduce)) { expr_free(r); return NULL; }   /* engine declined */
    bool r_true  = is_sym(r, SYM_True);
    bool r_false = is_sym(r, SYM_False);
    expr_free(r);
    if (quant == 0) return expr_new_symbol(r_false ? SYM_False : SYM_True);  /* Exists */
    return expr_new_symbol(r_true ? SYM_True : SYM_False);                    /* ForAll */
}

/* ------------------------------------------------------------------ *
 *  Case B -- parametric single-free-variable QE (via reduce_cad_qe)    *
 * ------------------------------------------------------------------ */

static Expr* qe_parametric(const Expr* body, int quant, const char** FREE, int nfree,
                           const char** B, int nb) {
    int nvall = nfree + nb;
    Expr** vall = malloc((size_t)nvall * sizeof(Expr*));
    for (int i = 0; i < nfree; i++) vall[i]         = expr_new_symbol(FREE[i]);
    for (int i = 0; i < nb; i++)    vall[nfree + i] = expr_new_symbol(B[i]);

    bool ok = true;
    RForm* F = reduce_form_from_expr(body, vall, nvall, &ok);
    Expr* out = NULL;
    if (ok) {
        rform_simplify(F, vall, nvall);
        out = reduce_cad_qe(F, vall, nfree, &vall[nfree], nb, quant);
    }
    rform_free(F);
    for (int i = 0; i < nvall; i++) expr_free(vall[i]);
    free(vall);
    return out;
}

/* Rebuild a single quantifier node `quant[{B...}, body]` (0 = Exists, 1 = ForAll).
 * CONSUMES `body`; the bound spec is a List of the names B[0..nb-1]. */
static Expr* qe_rebuild_quant(int quant, const char** B, int nb, Expr* body) {
    Expr** vs = malloc((size_t)(nb > 0 ? nb : 1) * sizeof(Expr*));
    for (int i = 0; i < nb; i++) vs[i] = expr_new_symbol(B[i]);
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)nb);
    free(vs);
    return mkfun2(quant == 0 ? SYM_Exists : SYM_ForAll, spec, body);
}

/* strcmp-order a free-variable name array in place, so the outermost CAD level is
 * canonical (alphabetical) regardless of the order the names were discovered in
 * the body.  A no-op for nfree <= 1.  Interned names compare by their string. */
static void qe_sort_names(const char** a, int n) {
    for (int i = 1; i < n; i++) {
        const char* key = a[i];
        int j = i - 1;
        while (j >= 0 && strcmp(a[j], key) > 0) { a[j + 1] = a[j]; j--; }
        a[j + 1] = key;
    }
}

/* ------------------------------------------------------------------ *
 *  Dispatch                                                           *
 * ------------------------------------------------------------------ */

Expr* reduce_qe_dispatch(const Expr* qexpr, const Expr* dom) {
    /* v1: quantified problems are over the Reals; an explicit non-Reals domain
     * (Complexes / Integers / Rationals / ...) declines. */
    if (dom && !is_sym(dom, SYM_Reals)) return NULL;

    int quant = -1; const char** B = NULL; int nb = 0, bcap = 0;
    bool alternating = false, ok = true;
    Expr* body = qe_normalize(qexpr, &quant, &B, &nb, &bcap, &alternating, &ok);
    if (!ok || !body) { expr_free(body); free(B); return NULL; }

    /* Alternating prefix: `body` is the inner different-kind quantifier.  Eliminate
     * it first to a quantifier-free formula psi, then re-eliminate the peeled outer
     * block over psi (now non-alternating).  This composes to arbitrary alternation
     * depth.  Any inner decline -- or an inner result still carrying a quantifier /
     * unevaluated Reduce/Resolve -- declines the whole thing. */
    if (alternating) {
        Expr* psi = reduce_qe_dispatch(body, dom);
        expr_free(body);
        if (!psi || is_head(psi, SYM_Exists) || is_head(psi, SYM_ForAll)
                 || is_head(psi, SYM_Reduce) || is_head(psi, SYM_Resolve)) {
            if (psi) expr_free(psi);
            free(B);
            return NULL;
        }
        Expr* q2 = qe_rebuild_quant(quant, B, nb, psi);   /* consumes psi */
        Expr* r = reduce_qe_dispatch(q2, dom);            /* borrows q2 */
        expr_free(q2);
        free(B);
        return r;
    }

    /* nbound==0: Exists[{},g] == ForAll[{},g] == g -- reduce g over its own
     * variables (or evaluate it when it is a constant statement). */
    if (nb == 0) {
        const char** FV = NULL; int nfv = 0, fcap = 0;
        qe_collect_symbols(body, &FV, &nfv, &fcap);
        Expr* r = (nfv == 0) ? eval_and_free(expr_copy(body))
                             : qe_call_reduce(body, FV, nfv);
        free(FV); free(B); expr_free(body);
        if (r && is_head(r, SYM_Reduce)) { expr_free(r); return NULL; }
        return r;
    }

    /* Free vars = leaf symbols of the body minus the bound vars (a bound var
     * that also appears free is shadowed by the binding). */
    const char** allsyms = NULL; int nas = 0, acap = 0;
    qe_collect_symbols(body, &allsyms, &nas, &acap);
    const char** FREE = malloc((size_t)(nas > 0 ? nas : 1) * sizeof(char*));
    int nfree = 0;
    for (int i = 0; i < nas; i++) {
        bool bound = false;
        for (int j = 0; j < nb; j++) if (allsyms[i] == B[j]) { bound = true; break; }
        if (!bound) FREE[nfree++] = allsyms[i];
    }
    free(allsyms);

    Expr* out;
    if (nfree == 0) {
        out = qe_decide(body, quant, B, nb);                 /* Case A: decision */
    } else {
        /* Cases B (nfree==1) & C (nfree>=2): parametric QE over the free vars,
         * outermost in canonical (alphabetical) order so the answer is stable. */
        qe_sort_names(FREE, nfree);
        out = qe_parametric(body, quant, FREE, nfree, B, nb);
    }

    free(FREE); free(B); expr_free(body);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Builtins                                                           *
 * ------------------------------------------------------------------ */

/* Exists / ForAll are inert: HoldAll keeps the bound variables symbolic and the
 * head carries a binding for Reduce / Resolve to eliminate.  They do not
 * evaluate on their own. */
Expr* builtin_exists(Expr* res) { (void)res; return NULL; }
Expr* builtin_forall(Expr* res) { (void)res; return NULL; }

Expr* builtin_resolve(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* qexpr = res->data.function.args[0];
    Expr* dom   = (argc >= 2) ? res->data.function.args[1] : NULL;
    if (!is_quantifier(qexpr)) return NULL;      /* only quantified input */
    return reduce_qe_dispatch(qexpr, dom);
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void reduce_qe_init(void) {
    symtab_add_builtin("Exists", builtin_exists);
    symtab_add_builtin("ForAll", builtin_forall);
    symtab_add_builtin("Resolve", builtin_resolve);

    SymbolDef* d;
    d = symtab_get_def("Exists"); if (d) d->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    d = symtab_get_def("ForAll"); if (d) d->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    d = symtab_get_def("Resolve"); if (d) d->attributes |= ATTR_PROTECTED;

    symtab_set_docstring("Exists",
        "Exists[x, expr]\n"
        "\tThe quantified statement that there exists a value of x for which\n"
        "\texpr is True.  Exists[{x1, x2, ...}, expr] binds several variables\n"
        "\tand Exists[x, cond, expr] restricts to values satisfying cond.\n"
        "\tExists is inert on its own (HoldAll); it is eliminated by Reduce\n"
        "\tor Resolve over the reals.");
    symtab_set_docstring("ForAll",
        "ForAll[x, expr]\n"
        "\tThe quantified statement that expr is True for all values of x.\n"
        "\tForAll[{x1, x2, ...}, expr] binds several variables and\n"
        "\tForAll[x, cond, expr] quantifies over values satisfying cond.\n"
        "\tForAll is inert on its own (HoldAll); it is eliminated by Reduce\n"
        "\tor Resolve over the reals.");
    symtab_set_docstring("Resolve",
        "Resolve[expr]\n"
        "Resolve[expr, dom]\n"
        "\tEliminates the quantifiers (Exists, ForAll) from expr over the\n"
        "\tdomain dom (Reals; the default and only supported domain), returning\n"
        "\tan equivalent quantifier-free statement -- True or False for a fully\n"
        "\tquantified sentence, or a condition on the remaining free variables\n"
        "\t(one or more). Alternating quantifier prefixes are eliminated\n"
        "\tinner-block-first. An undecidable sign, a non-Reals domain, or an\n"
        "\tunsupported case is left unevaluated rather than guessed.");
}
