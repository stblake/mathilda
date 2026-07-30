/* Mathilda — autocompile adapter.  See autocompile.h. */

#include "autocompile.h"
#include "compile.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
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
                             CompileType argt) {
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
    CompiledProgram* prog = compile_expr_ex(body, names, types, nvars, ac_flags);
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
    return ac_make(body, vars, nvars, CT_REAL);
}

AutoCompiled* autocompile_new_z(const Expr* body, const Expr* const* vars, size_t nvars) {
    return ac_make(body, vars, nvars, CT_COMPLEX);
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

void autocompiled_free(AutoCompiled* ac) {
    if (!ac) return;
    if (ac->prog) compiled_free(ac->prog);
    free(ac);
}
