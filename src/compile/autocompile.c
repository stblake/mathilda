/* Mathilda — autocompile adapter.  See autocompile.h. */

#include "autocompile.h"
#include "compile.h"

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "../expr.h"
#include "../sym_intern.h"   /* intern_symbol */

#define AC_MAX_VARS 16       /* boxed fallback path caps here; real path is unbounded */

struct AutoCompiled {
    CompiledProgram* prog;
    size_t           nvars;
    bool             real_result;   /* result type is CT_REAL → all-real fast path */
};

AutoCompiled* autocompile_new(const Expr* body, const Expr* const* vars, size_t nvars) {
    if (!body || nvars == 0) return NULL;
    const char** names = malloc(nvars * sizeof(*names));
    CompileType* types = malloc(nvars * sizeof(*types));
    if (!names || !types) { free(names); free(types); return NULL; }
    for (size_t i = 0; i < nvars; i++) {
        if (!vars[i] || vars[i]->type != EXPR_SYMBOL) { free(names); free(types); return NULL; }
        names[i] = intern_symbol(vars[i]->data.symbol.name);
        types[i] = CT_REAL;
    }
    CompiledProgram* prog = compile_expr(body, names, types, nvars);
    free(names); free(types);
    if (!prog) return NULL;

    AutoCompiled* ac = calloc(1, sizeof *ac);
    if (!ac) { compiled_free(prog); return NULL; }
    ac->prog = prog;
    ac->nvars = nvars;
    ac->real_result = (compiled_result_type(prog) == CT_REAL);
    return ac;
}

size_t autocompiled_num_vars(const AutoCompiled* ac) { return ac->nvars; }

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

void autocompiled_free(AutoCompiled* ac) {
    if (!ac) return;
    if (ac->prog) compiled_free(ac->prog);
    free(ac);
}
