/* Mathilda — autocompile adapter.  See autocompile.h. */

#include "autocompile.h"
#include "compile.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "../expr.h"
#include "../arithmetic.h"   /* make_complex */
#include "../sym_intern.h"   /* intern_symbol */
#include "../print.h"        /* expr_to_string — the diagnostic below */

#define AC_MAX_VARS 16       /* boxed fallback path caps here; real path is unbounded */

/* When a numeric builtin's body does not compile, that builtin silently keeps
 * its interpreter path: same answer, 10-40x slower, no symptom anywhere.  Set
 * MATHILDA_COMPILE_DIAG=1 to have each such fallback name itself and, more
 * usefully, name the single subexpression that caused it — the compilable subset
 * is a cliff, so one head outside it costs the whole body.
 *
 * The env var is read once: this runs on a path that is already slow, but a
 * getenv per fallback would still be per-Plot-call noise for no reason. */
static void ac_report_bail(const Expr* body) {
    static int enabled = -1;
    if (enabled < 0) {
        const char* s = getenv("MATHILDA_COMPILE_DIAG");
        enabled = (s && *s && *s != '0') ? 1 : 0;
    }
    if (!enabled) return;
    const char* why  = compiled_bail_reason();
    const char* what = compiled_bail_expr();
    char* full = expr_to_string((Expr*)body);
    fprintf(stderr, "[compile] interpreting: %s\n           cause: %s%s%s\n",
            full ? full : "?", what ? what : "(none)",
            why ? " — " : "", why ? why : "");
    free(full);
}

struct AutoCompiled {
    CompiledProgram* prog;
    size_t           nvars;
    bool             real_result;   /* result type is CT_REAL → all-real fast path */
    CompileType      arg_type;      /* CT_REAL, or CT_COMPLEX for an _z program */
};

/* -1 = not yet read from the environment. Same shape as src/numloop.c's and
 * src/pack.c's switches. */
static int g_autocompile_enabled = -1;

bool autocompile_enabled(void) {
    if (g_autocompile_enabled < 0)
        g_autocompile_enabled = getenv("MATHILDA_NO_AUTOCOMPILE") ? 0 : 1;
    return g_autocompile_enabled != 0;
}

void autocompile_set_enabled(bool on) { g_autocompile_enabled = on ? 1 : 0; }

static AutoCompiled* ac_make(const Expr* body, const Expr* const* vars, size_t nvars,
                             CompileType argt, long prec_bits) {
    if (!body || nvars == 0) return NULL;
    /* Off: every caller sees the same NULL it sees for a body outside the
     * compilable subset, so each keeps its interpreter path with no other
     * change. */
    if (!autocompile_enabled()) return NULL;
    const char** names = malloc(nvars * sizeof(*names));
    CompileType* types = malloc(nvars * sizeof(*types));
    if (!names || !types) { free(names); free(types); return NULL; }
    for (size_t i = 0; i < nvars; i++) {
        if (!vars[i] || vars[i]->type != EXPR_SYMBOL) { free(names); free(types); return NULL; }
        names[i] = intern_symbol(vars[i]->data.symbol.name);
        types[i] = argt;
    }
    /* FOLD_GLOBALS is safe here and nowhere else: an AutoCompiled is built and
     * freed inside one builtin call, so a folded symbol (e.g. the outer
     * iteration variable of a nested Table) cannot be reassigned while the
     * program lives.  It is what lets the inner Table of
     * Table[f[x,y], {y,..}, {x,..}] compile at all.
     *
     * COMPILE_WRAP_INT IS DELIBERATELY NOT SET, and must never be.  The user did
     * not ask for any of this — Plot, Table, NIntegrate and friends compile
     * behind their backs, purely as an optimisation — so an auto-compiled path
     * that wrapped a machine integer would turn a transparent speed-up into a
     * silently wrong answer in code that never mentioned Compile.  Checking on
     * means the worst an overflow can cost here is the interpreter re-running
     * that call.
     *
     * The mask makes that structural rather than a property of this literal:
     * `SetOptions[Compile, RuntimeOptions -> ...]` moves the default for the
     * USER-facing Compile[] only, and if the two are ever wired together this
     * line still refuses to let the auto-compiled path inherit it. */
    unsigned ac_flags = COMPILE_FOLD_GLOBALS & ~COMPILE_WRAP_INT;
    CompiledProgram* prog = compile_expr_prec(body, names, types, nvars, ac_flags, prec_bits);
    free(names); free(types);
    if (!prog) { ac_report_bail(body); return NULL; }

    AutoCompiled* ac = calloc(1, sizeof *ac);
    if (!ac) { compiled_free(prog); return NULL; }
    ac->prog = prog;
    ac->nvars = nvars;
    ac->arg_type = argt;
    /* The unboxed entry point needs an all-Real SIGNATURE, not just a real
     * result — a complex-argument program never qualifies. */
    ac->real_result = (argt == CT_REAL) && (compiled_result_type(prog) == CT_REAL);
    return ac;
}

AutoCompiled* autocompile_new(const Expr* body, const Expr* const* vars, size_t nvars) {
    return ac_make(body, vars, nvars, CT_REAL, 0);
}

AutoCompiled* autocompile_new_z(const Expr* body, const Expr* const* vars, size_t nvars) {
    return ac_make(body, vars, nvars, CT_COMPLEX, 0);
}

AutoCompiled* autocompile_new_prec(const Expr* body, const Expr* const* vars,
                                   size_t nvars, long prec_bits) {
    if (prec_bits <= 0) return NULL;   /* use autocompile_new for the machine path */
    return ac_make(body, vars, nvars, CT_BIGREAL, prec_bits);
}

AutoCompiled* autocompile_new_prec_z(const Expr* body, const Expr* const* vars,
                                     size_t nvars, long prec_bits) {
    if (prec_bits <= 0) return NULL;
    return ac_make(body, vars, nvars, CT_BIGCOMPLEX, prec_bits);
}

size_t autocompiled_num_vars(const AutoCompiled* ac) { return ac->nvars; }

bool autocompiled_result_is_real(const AutoCompiled* ac) {
    return ac && ac->real_result;
}

bool autocompiled_eval_real(const AutoCompiled* ac, const double* xs, double* out) {
    if (ac->real_result)
        return compiled_eval_real(ac->prog, xs, out);   /* all-real: no boxing, self-guards finite */

    /* Non-real result type (INT/COMPLEX/BOOL): box the inputs, then accept only a
     * real-valued result. */
    if (ac->nvars > AC_MAX_VARS) return false;
    CompileValue args[AC_MAX_VARS], o;
    for (size_t i = 0; i < ac->nvars; i++) { args[i].type = CT_REAL; args[i].v.r = xs[i]; }
    if (!compiled_eval(ac->prog, args, &o)) return false;
    switch (o.type) {
        case CT_INT:  *out = (double)o.v.i; return isfinite(*out);
        case CT_REAL: *out = o.v.r;         return isfinite(*out);
        case CT_COMPLEX:
            if (cimag(o.v.z) == 0.0) { *out = creal(o.v.z); return isfinite(*out); }
            return false;   /* genuinely complex → no real value here */
        default: return false;
    }
}

bool autocompiled_eval_complex(const AutoCompiled* ac, const double* xs, double _Complex* out) {
    if (ac->real_result) {
        double y;
        if (!compiled_eval_real(ac->prog, xs, &y)) return false;
        *out = y;   /* real result, zero imaginary part */
        return true;
    }
    if (ac->nvars > AC_MAX_VARS) return false;
    CompileValue args[AC_MAX_VARS], o;
    for (size_t i = 0; i < ac->nvars; i++) { args[i].type = CT_REAL; args[i].v.r = xs[i]; }
    if (!compiled_eval(ac->prog, args, &o)) return false;   /* false ⇒ non-finite */
    switch (o.type) {
        case CT_INT:     *out = (double)o.v.i; return true;
        case CT_REAL:    *out = o.v.r;         return true;
        case CT_COMPLEX: *out = o.v.z;         return true;
        default:         return false;         /* BOOL: not a number */
    }
}

Expr* autocompiled_eval_boxed(const AutoCompiled* ac, const double* xs) {
    if (ac->real_result) {                      /* all-real: no boxing, self-guards finite */
        double y;
        return compiled_eval_real(ac->prog, xs, &y) ? expr_new_real(y) : NULL;
    }
    if (ac->nvars > AC_MAX_VARS) return NULL;
    CompileValue args[AC_MAX_VARS], o;
    for (size_t i = 0; i < ac->nvars; i++) { args[i].type = CT_REAL; args[i].v.r = xs[i]; }
    if (!compiled_eval(ac->prog, args, &o)) return NULL;
    switch (o.type) {
        case CT_INT:  return expr_new_integer(o.v.i);
        case CT_REAL: return isfinite(o.v.r) ? expr_new_real(o.v.r) : NULL;
        case CT_COMPLEX:
            if (!isfinite(creal(o.v.z)) || !isfinite(cimag(o.v.z))) return NULL;
            /* A zero imaginary part is reported as a plain real, matching how the
             * interpreter's arithmetic collapses Complex[r, 0.]. */
            return cimag(o.v.z) == 0.0
                 ? expr_new_real(creal(o.v.z))
                 : make_complex(expr_new_real(creal(o.v.z)), expr_new_real(cimag(o.v.z)));
        default: return NULL;                   /* BOOL / array: interpreter handles it */
    }
}

Expr* autocompile_eval_closed(const Expr* expr) {
    if (!expr) return NULL;
    if (!autocompile_enabled()) return NULL;
    /* Same flag mask as ac_make: fold machine-number globals (so a bound such as
     * an outer iteration variable folds to a constant), but NEVER wrap a machine
     * integer — an invisible overflow must bail to the interpreter, not silently
     * wrap.  See the long comment in ac_make. */
    unsigned ac_flags = COMPILE_FOLD_GLOBALS & ~COMPILE_WRAP_INT;
    CompiledProgram* prog = compile_expr_prec(expr, NULL, NULL, 0, ac_flags, 0);
    if (!prog) { ac_report_bail(expr); return NULL; }
    /* Only a genuine machine-number total may leave the compiled path.  A CT_INT
     * result is an exact integer sum/product and must stay on the interpreter
     * (exact bignum); CT_BOOL / array / managed cannot be a scalar total. */
    CompileType rt = compiled_result_type(prog);
    if (rt != CT_REAL && rt != CT_COMPLEX) { compiled_free(prog); return NULL; }

    CompileValue o;
    Expr* out = NULL;
    if (compiled_eval(prog, NULL, &o)) {   /* zero-arg program: no args to load */
        switch (o.type) {
            case CT_REAL:
                if (isfinite(o.v.r)) out = expr_new_real(o.v.r);
                break;
            case CT_COMPLEX:
                if (isfinite(creal(o.v.z)) && isfinite(cimag(o.v.z)))
                    /* A zero imaginary part collapses to a plain real, matching
                     * how the interpreter's arithmetic reports Complex[r, 0.]. */
                    out = cimag(o.v.z) == 0.0
                        ? expr_new_real(creal(o.v.z))
                        : make_complex(expr_new_real(creal(o.v.z)),
                                       expr_new_real(cimag(o.v.z)));
                break;
            default: break;   /* result type is fixed to REAL/COMPLEX above */
        }
    }
    compiled_free(prog);
    return out;
}

bool autocompiled_eval_z(const AutoCompiled* ac, const double _Complex* zs,
                         double _Complex* out) {
    if (!ac || ac->arg_type != CT_COMPLEX || ac->nvars > AC_MAX_VARS) return false;
    CompileValue args[AC_MAX_VARS], o;
    for (size_t i = 0; i < ac->nvars; i++) { args[i].type = CT_COMPLEX; args[i].v.z = zs[i]; }
    if (!compiled_eval(ac->prog, args, &o)) return false;   /* false ⇒ non-finite */
    switch (o.type) {
        case CT_INT:     *out = (double)o.v.i; return true;
        case CT_REAL:    *out = o.v.r;         return true;
        case CT_COMPLEX: *out = o.v.z;         return true;
        default:         return false;         /* BOOL / array: not a sample value */
    }
}

Expr* autocompiled_eval_mpfr(const AutoCompiled* ac, const Expr* const* xs) {
    if (!ac || !CT_IS_MANAGED(ac->arg_type) || ac->nvars > AC_MAX_VARS) return NULL;
    CompileValue args[AC_MAX_VARS], o;
    for (size_t i = 0; i < ac->nvars; i++) {
        if (!xs[i]) return NULL;
        args[i].type = ac->arg_type;
        args[i].v.a = (Expr*)xs[i];   /* borrowed numeric Expr; load_arg sets the container */
    }
    if (!compiled_eval(ac->prog, args, &o)) return NULL;   /* false ⇒ non-finite / non-numeric */
    /* A managed result is a NEW Expr (EXPR_MPFR / Complex[MPFR,MPFR]) the caller
     * owns; on any non-managed result the caller should interpret. */
    return CT_IS_MANAGED(o.type) ? o.v.a : NULL;
}

/* ------------------------------------------------------------------------
 * autocompile_map_unary — f applied to every element of a machine vector, as
 * ONE compiled loop.  See autocompile.h for the contract.
 * ------------------------------------------------------------------------ */

/* The element placeholder: `Part[V, I]` inside the loop body. */
static const char* ac_map_v(void) { return intern_symbol("AutoCompile`MapVector"); }
static const char* ac_map_i(void) { return intern_symbol("AutoCompile`MapIndex"); }

static bool ac_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static CompileType ac_join(CompileType a, CompileType b) {
    if (a == CT_ERR || b == CT_ERR) return CT_ERR;
    if (a == CT_BOOL || b == CT_BOOL) return (a == b) ? CT_BOOL : CT_ERR;
    return (a == CT_REAL || b == CT_REAL) ? CT_REAL : CT_INT;
}

/* The EXACTNESS GATE.  A tiny type inference over the substituted body that
 * accepts only constructs whose machine lowering computes exactly what the
 * interpreter would: integer ops that stay integers (checked for overflow by the
 * VM), real ops on reals, and the int -> int roundings/tests.  It rejects what
 * the compiler would silently turn Real — `x/2` or `x^-1` on integers is a
 * Rational to the interpreter and a double to the VM — plus any symbol it cannot
 * type (a global may hold a Rational) and every head it does not know.
 * Returns CT_INT / CT_REAL / CT_BOOL, or CT_ERR to refuse. */
static CompileType ac_exact_type(const Expr* e, CompileType elem) {
    if (!e) return CT_ERR;
    switch (e->type) {
        case EXPR_INTEGER: return CT_INT;
        case EXPR_REAL:    return CT_REAL;
        case EXPR_SYMBOL: {
            const char* s = e->data.symbol.name;
            if (strcmp(s, "True") == 0 || strcmp(s, "False") == 0) return CT_BOOL;
            return CT_ERR;
        }
        case EXPR_FUNCTION: break;
        default: return CT_ERR;
    }
    const Expr* hd = e->data.function.head;
    if (hd->type != EXPR_SYMBOL) return CT_ERR;
    const char* h = hd->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr** A = e->data.function.args;

    if (strcmp(h, "Part") == 0)      /* the element placeholder Part[V, I] */
        return (n == 2 && A[0]->type == EXPR_SYMBOL && A[0]->data.symbol.name == ac_map_v())
               ? elem : CT_ERR;

    CompileType t[8];
    if (n > 8) return CT_ERR;
    for (size_t k = 0; k < n; k++) if ((t[k] = ac_exact_type(A[k], elem)) == CT_ERR) return CT_ERR;

    /* Arithmetic closed over Int and Real.  A Real sum or product of more than
     * two terms is refused: Plus/Times are Orderless, so the interpreter adds
     * the terms in canonical order while the VM adds them as written, and
     * floating-point addition is not associative. */
    static const char* const ARITH[] = { "Plus", "Times", "Subtract", "Max", "Min", NULL };
    for (int k = 0; ARITH[k]; k++)
        if (strcmp(h, ARITH[k]) == 0) {
            if (n == 0) return CT_ERR;
            CompileType r = t[0];
            for (size_t j = 1; j < n; j++) r = ac_join(r, t[j]);
            if (r == CT_REAL && n > 2) return CT_ERR;
            return r == CT_BOOL ? CT_ERR : r;
        }
    /* Mod and Equal/Unequal on Reals are not bit-for-bit the VM's: the
     * interpreter's real Equal is tolerant, and its real Mod is x - m Floor[x/m]
     * evaluated as a tree.  Integers only. */
    if (strcmp(h, "Mod") == 0)
        return (n == 2 && t[0] == CT_INT && t[1] == CT_INT) ? CT_INT : CT_ERR;
    if (strcmp(h, "Equal") == 0 || strcmp(h, "Unequal") == 0) {
        if (n < 2) return CT_ERR;
        for (size_t j = 0; j < n; j++) if (t[j] != CT_INT) return CT_ERR;
        return CT_BOOL;
    }
    if (strcmp(h, "Minus") == 0 || strcmp(h, "Abs") == 0)
        return (n == 1 && t[0] != CT_BOOL) ? t[0] : CT_ERR;
    /* Numeric -> Integer on both sides. */
    if (strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0 ||
        strcmp(h, "Sign") == 0 || strcmp(h, "IntegerPart") == 0)
        return (n == 1 && t[0] != CT_BOOL) ? CT_INT : CT_ERR;
    if (strcmp(h, "Quotient") == 0)
        return (n == 2 && t[0] == CT_INT && t[1] == CT_INT) ? CT_INT : CT_ERR;
    if (strcmp(h, "Boole") == 0) return (n == 1 && t[0] == CT_BOOL) ? CT_INT : CT_ERR;
    if (strcmp(h, "EvenQ") == 0 || strcmp(h, "OddQ") == 0)
        return (n == 1 && t[0] == CT_INT) ? CT_BOOL : CT_ERR;
    if (strcmp(h, "Positive") == 0 || strcmp(h, "Negative") == 0 ||
        strcmp(h, "NonNegative") == 0 || strcmp(h, "NonPositive") == 0)
        return (n == 1 && t[0] != CT_BOOL) ? CT_BOOL : CT_ERR;
    static const char* const CMP[] = { "Less", "LessEqual", "Greater", "GreaterEqual", NULL };
    for (int k = 0; CMP[k]; k++)
        if (strcmp(h, CMP[k]) == 0) {
            if (n < 2) return CT_ERR;
            for (size_t j = 0; j < n; j++) if (t[j] == CT_BOOL) return CT_ERR;
            return CT_BOOL;
        }
    if (strcmp(h, "And") == 0 || strcmp(h, "Or") == 0 || strcmp(h, "Xor") == 0 ||
        strcmp(h, "Not") == 0) {
        if (n == 0) return CT_ERR;
        for (size_t j = 0; j < n; j++) if (t[j] != CT_BOOL) return CT_ERR;
        return CT_BOOL;
    }
    if (strcmp(h, "If") == 0)
        return (n == 3 && t[0] == CT_BOOL) ? ac_join(t[1], t[2]) : CT_ERR;
    /* Divide / Power stay exact only when a Real is involved, or for an integer
     * raised to a non-negative integer literal. */
    if (strcmp(h, "Divide") == 0)
        return (n == 2 && t[0] != CT_BOOL && t[1] != CT_BOOL &&
                (t[0] == CT_REAL || t[1] == CT_REAL)) ? CT_REAL : CT_ERR;
    if (strcmp(h, "Power") == 0) {
        if (n != 2 || t[0] == CT_BOOL || t[1] == CT_BOOL) return CT_ERR;
        if (t[0] == CT_REAL || t[1] == CT_REAL) return CT_REAL;
        return (A[1]->type == EXPR_INTEGER && A[1]->data.integer >= 0) ? CT_INT : CT_ERR;
    }
    return CT_ERR;
}

/* Copy of `body` with the element variable replaced by `elt` (borrowed): Slot[1]
 * for a slot-form Function (param == NULL), else the named symbol `param`. */
static Expr* ac_subst(const Expr* body, const char* param, const Expr* elt) {
    if (param) {
        if (body->type == EXPR_SYMBOL && body->data.symbol.name == param) return expr_copy((Expr*)elt);
    } else if (ac_head_is(body, "Slot") && body->data.function.arg_count == 1 &&
               body->data.function.args[0]->type == EXPR_INTEGER &&
               body->data.function.args[0]->data.integer == 1) {
        return expr_copy((Expr*)elt);
    }
    if (body->type != EXPR_FUNCTION) return expr_copy((Expr*)body);
    size_t n = body->data.function.arg_count;
    Expr* out = expr_new_function(ac_subst(body->data.function.head, param, elt), NULL, n);
    for (size_t k = 0; k < n; k++)
        out->data.function.args[k] = ac_subst(body->data.function.args[k], param, elt);
    return out;
}

Expr* autocompile_map_unary(const Expr* f, const Expr* arr) {
    if (!f || !arr || !autocompile_enabled()) return NULL;
    if (arr->type != EXPR_NDARRAY || arr->data.ndarray.rank != 1) return NULL;
    NDType dt = arr->data.ndarray.dtype;
    if (dt != NDT_INT64 && dt != NDT_FLOAT64) return NULL;
    int64_t len = arr->data.ndarray.dims[0];
    if (len < 1) return NULL;
    CompileType elem = (dt == NDT_INT64) ? CT_INT : CT_REAL;

    const char* vname = ac_map_v();
    const char* iname = ac_map_i();
    Expr* pargs[2] = { expr_new_symbol(vname), expr_new_symbol(iname) };
    Expr* elt = expr_new_function(expr_new_symbol("Part"), pargs, 2);

    /* The scalar body, element variable replaced by Part[V, I]. */
    Expr* body = NULL;
    if (f->type == EXPR_SYMBOL) {
        Expr* a1 = expr_copy(elt);
        body = expr_new_function(expr_copy((Expr*)f), &a1, 1);
    } else if (ac_head_is(f, "Function")) {
        size_t fa = f->data.function.arg_count;
        if (fa == 1) {
            body = ac_subst(f->data.function.args[0], NULL, elt);
        } else if (fa == 2) {
            const Expr* p = f->data.function.args[0];
            if (ac_head_is(p, "List") && p->data.function.arg_count == 1) p = p->data.function.args[0];
            if (p->type == EXPR_SYMBOL)
                body = ac_subst(f->data.function.args[1], p->data.symbol.name, elt);
        }
    }
    expr_free(elt);
    if (!body) return NULL;
    if (ac_exact_type(body, elem) == CT_ERR) { expr_free(body); return NULL; }

    /* Table[body, {I, Length[V]}] over the declared vector V: one VM call. */
    Expr* lenargs[1] = { expr_new_symbol(vname) };
    Expr* iter[2] = { expr_new_symbol(iname),
                      expr_new_function(expr_new_symbol("Length"), lenargs, 1) };
    Expr* targs[2] = { body, expr_new_function(expr_new_symbol("List"), iter, 2) };
    Expr* table = expr_new_function(expr_new_symbol("Table"), targs, 2);

    const char* names[1] = { vname };
    CompileType types[1] = { CT_ARRAY(elem, 1) };
    CompiledProgram* prog = compile_expr_prec(table, names, types, 1, 0u, 0);
    expr_free(table);
    if (!prog) return NULL;

    CompileValue args[1], o;
    args[0].type = types[0];
    args[0].v.a = (Expr*)arr;                     /* borrowed by the program */
    Expr* out = NULL;
    if (compiled_eval(prog, args, &o)) {
        if (CT_IS_ARRAY(o.type) && o.v.a && o.v.a->type == EXPR_NDARRAY &&
            o.v.a->data.ndarray.rank == 1 && o.v.a->data.ndarray.dims[0] == len &&
            (o.v.a->data.ndarray.dtype == NDT_INT64 || o.v.a->data.ndarray.dtype == NDT_FLOAT64 ||
             o.v.a->data.ndarray.dtype == NDT_BOOL))
            out = o.v.a;
        else if (CT_IS_ARRAY(o.type) && o.v.a)
            expr_free(o.v.a);
    }
    compiled_free(prog);
    return out;
}

void autocompiled_free(AutoCompiled* ac) {
    if (!ac) return;
    if (ac->prog) compiled_free(ac->prog);
    free(ac);
}
