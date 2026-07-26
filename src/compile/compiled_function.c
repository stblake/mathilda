/* Mathilda — CompiledFunction: user-facing Compile[] wrapper.  See
 * compiled_function.h and docs/design/compile.md (milestone M1b). */

#include "compiled_function.h"
#include "compile.h"

#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "../expr.h"
#include "../arithmetic.h"     /* is_complex, make_complex, is_rational */
#include "../sym_intern.h"     /* intern_symbol */
#include "../sym_names.h"      /* SYM_Real / SYM_Integer / SYM_Complex / ... */
#include "../symtab.h"         /* symtab_add_builtin / _set_docstring / _get_def */
#include "../attr.h"           /* ATTR_HOLDALL / ATTR_PROTECTED */
#include "../match.h"          /* env_new / env_set / replace_bindings */
#include "../eval.h"           /* evaluate / eval_and_free */

struct CompiledFunction {
    unsigned         refcount;   /* Expr copies share the payload */
    size_t           nargs;
    const char**     arg_names;  /* interned pointers (interner owns the chars) */
    CompileType*     arg_types;
    Expr*            body;       /* owned; interpreter fallback + printing */
    CompiledProgram* prog;       /* NULL if body is outside the compilable subset */
};

/* ------------------------------------------------------------------ *
 *  Numeric boxing / unboxing                                          *
 * ------------------------------------------------------------------ */

/* Concrete real value of a numeric atom (Integer/BigInt/Real/MPFR/Rational). */
static bool cf_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;          return true;
        case EXPR_REAL:    *out = e->data.real;                     return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);        return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den; expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        double dd = mpz_get_d(den);
        bool ok = (dd != 0.0);
        if (ok) *out = mpz_get_d(num) / dd;
        mpz_clears(num, den, NULL);
        return ok;
    }
    return false;
}

static bool cf_to_ll(const Expr* e, long long* out) {
    if (e->type == EXPR_INTEGER) { *out = (long long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (!mpz_fits_slong_p(e->data.bigint)) return false;
        *out = (long long)mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

static bool cf_to_complex(const Expr* e, double* re, double* im) {
    Expr *r, *i;
    if (is_complex((Expr*)e, &r, &i)) {
        double a, b;
        if (!cf_to_double(r, &a) || !cf_to_double(i, &b)) return false;
        *re = a; *im = b; return true;
    }
    if (cf_to_double(e, re)) { *im = 0.0; return true; }
    return false;
}

/* Box an argument Expr into a CompileValue of the declared type.  Returns false
 * if the argument is not a concrete number of that type (→ interpreter
 * fallback). */
static bool cf_box(const Expr* e, CompileType t, CompileValue* out) {
    out->type = t;
    switch (t) {
        case CT_BOOL:
            if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True)  { out->v.b = 1; return true; }
            if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False) { out->v.b = 0; return true; }
            return false;
        case CT_INT:     return cf_to_ll(e, &out->v.i);
        case CT_REAL:    return cf_to_double(e, &out->v.r);
        case CT_COMPLEX: {
            double re, im;
            if (!cf_to_complex(e, &re, &im)) return false;
            out->v.z = re + im * I; return true;
        }
    }
    return false;
}

static Expr* cf_unbox(const CompileValue* v) {
    switch (v->type) {
        case CT_BOOL: return expr_new_symbol(v->v.b ? "True" : "False");
        case CT_INT:  return expr_new_integer((int64_t)v->v.i);
        case CT_REAL: return expr_new_real(v->v.r);
        case CT_COMPLEX: {
            double re = creal(v->v.z), im = cimag(v->v.z);
            if (im == 0.0) return expr_new_real(re);
            return make_complex(expr_new_real(re), expr_new_real(im));
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Argument-spec parsing                                              *
 * ------------------------------------------------------------------ */

static bool typesym_to_ct(const char* nm, CompileType* out) {
    if (nm == SYM_Real)    { *out = CT_REAL;    return true; }
    if (nm == SYM_Integer) { *out = CT_INT;     return true; }
    if (nm == SYM_Complex) { *out = CT_COMPLEX; return true; }
    return false;
}

/* `_Real` (Blank[Real]) / bare `Real` / Integer / Complex. */
static bool parse_typespec(const Expr* ts, CompileType* out) {
    if (ts->type == EXPR_SYMBOL) return typesym_to_ct(ts->data.symbol.name, out);
    if (ts->type == EXPR_FUNCTION && ts->data.function.head->type == EXPR_SYMBOL
        && ts->data.function.head->data.symbol.name == SYM_Blank
        && ts->data.function.arg_count == 1
        && ts->data.function.args[0]->type == EXPR_SYMBOL)
        return typesym_to_ct(ts->data.function.args[0]->data.symbol.name, out);
    return false;
}

/* ------------------------------------------------------------------ *
 *  Lifecycle                                                          *
 * ------------------------------------------------------------------ */

CompiledFunction* compiled_function_new(const Expr* argspec, const Expr* body) {
    if (!argspec || !body) return NULL;
    if (argspec->type != EXPR_FUNCTION || argspec->data.function.head->type != EXPR_SYMBOL
        || argspec->data.function.head->data.symbol.name != SYM_List
        || argspec->data.function.arg_count == 0)
        return NULL;

    size_t n = argspec->data.function.arg_count;
    const char** names = malloc(n * sizeof(*names));
    CompileType* types = malloc(n * sizeof(*types));
    if (!names || !types) { free(names); free(types); return NULL; }

    for (size_t i = 0; i < n; i++) {
        const Expr* el = argspec->data.function.args[i];
        const char* nm = NULL; CompileType ty = CT_REAL;
        if (el->type == EXPR_SYMBOL) {
            nm = el->data.symbol.name;
        } else if (el->type == EXPR_FUNCTION && el->data.function.head->type == EXPR_SYMBOL
                   && el->data.function.head->data.symbol.name == SYM_List
                   && el->data.function.arg_count == 2
                   && el->data.function.args[0]->type == EXPR_SYMBOL
                   && parse_typespec(el->data.function.args[1], &ty)) {
            nm = el->data.function.args[0]->data.symbol.name;
        } else { free(names); free(types); return NULL; }

        const char* in = intern_symbol(nm);
        for (size_t j = 0; j < i; j++)
            if (names[j] == in) { free(names); free(types); return NULL; }  /* duplicate param */
        names[i] = in;
        types[i] = ty;
    }

    CompiledFunction* cf = calloc(1, sizeof *cf);
    if (!cf) { free(names); free(types); return NULL; }
    cf->refcount  = 1;
    cf->nargs     = n;
    cf->arg_names = names;
    cf->arg_types = types;
    cf->body      = expr_copy((Expr*)body);
    cf->prog      = compile_expr(cf->body, names, types, n);  /* NULL ⇒ fallback only */
    return cf;
}

CompiledFunction* compiled_function_ref(CompiledFunction* cf) {
    if (cf) cf->refcount++;
    return cf;
}

void compiled_function_free(CompiledFunction* cf) {
    if (!cf) return;
    if (cf->refcount > 1) { cf->refcount--; return; }
    if (cf->prog) compiled_free(cf->prog);
    if (cf->body) expr_free(cf->body);
    free(cf->arg_names);
    free(cf->arg_types);
    free(cf);
}

uint64_t compiled_function_identity(const CompiledFunction* cf) {
    return (uint64_t)(uintptr_t)cf;
}

/* ------------------------------------------------------------------ *
 *  Application                                                        *
 * ------------------------------------------------------------------ */

/* Interpreter fallback: substitute the args for the parameter symbols in the
 * original body and evaluate.  replace_bindings shares the arg nodes by
 * refcount, so eval_and_free just dec-refs them — the caller keeps ownership. */
static Expr* cf_fallback(const CompiledFunction* cf, Expr* const* args, size_t nargs) {
    MatchEnv* env = env_new();
    for (size_t i = 0; i < nargs; i++) env_set(env, cf->arg_names[i], args[i]);
    Expr* sub = replace_bindings(cf->body, env);
    env_free(env);
    return eval_and_free(sub);
}

Expr* compiled_function_apply(const CompiledFunction* cf, Expr* const* args, size_t nargs) {
    if (!cf || nargs != cf->nargs) return NULL;   /* arity mismatch ⇒ leave unevaluated */

    if (cf->prog) {
        CompileValue* cv = malloc((nargs ? nargs : 1) * sizeof(*cv));
        if (cv) {
            bool all_numeric = true;
            for (size_t i = 0; i < nargs; i++)
                if (!cf_box(args[i], cf->arg_types[i], &cv[i])) { all_numeric = false; break; }
            if (all_numeric) {
                CompileValue out;
                if (compiled_eval(cf->prog, cv, &out)) { free(cv); return cf_unbox(&out); }
            }
            free(cv);
        }
    }
    return cf_fallback(cf, args, nargs);
}

/* ------------------------------------------------------------------ *
 *  Accessors + builtin                                                *
 * ------------------------------------------------------------------ */

size_t             compiled_function_num_args(const CompiledFunction* cf) { return cf ? cf->nargs : 0; }
bool               compiled_function_is_compiled(const CompiledFunction* cf) { return cf && cf->prog; }
const Expr*        compiled_function_body(const CompiledFunction* cf) { return cf ? cf->body : NULL; }
const char* const* compiled_function_arg_names(const CompiledFunction* cf) { return cf ? cf->arg_names : NULL; }
const CompileType* compiled_function_arg_types(const CompiledFunction* cf) { return cf ? cf->arg_types : NULL; }

/* Compile[argspec, body] (HoldAll).  Never evaluates the body; the raw held
 * body is compiled, and any non-arg symbol it references simply routes that
 * call through the interpreter fallback. */
static Expr* builtin_compile(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    CompiledFunction* cf = compiled_function_new(res->data.function.args[0],
                                                 res->data.function.args[1]);
    if (!cf) return NULL;   /* malformed argspec ⇒ leave Compile[...] unevaluated */
    return expr_new_compiled(cf);
}

void compiled_function_init(void) {
    symtab_add_builtin("Compile", builtin_compile);
    SymbolDef* d = symtab_get_def("Compile");
    if (d) d->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Compile",
        "Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a "
        "CompiledFunction that evaluates expr over machine numbers (types _Real, "
        "_Integer, _Complex; default _Real), falling back to the interpreter for "
        "symbolic arguments or non-compilable bodies.");
}
