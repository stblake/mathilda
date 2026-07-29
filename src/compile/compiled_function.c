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
#include "../ndarray.h"        /* ndarray_from_nested_list — packing a List argument */
#include "../symtab.h"         /* symtab_add_builtin / _set_docstring / _get_def */
#include "../attr.h"           /* ATTR_HOLDALL / ATTR_PROTECTED */
#include "../match.h"          /* env_new / env_set / replace_bindings */
#include "../eval.h"           /* evaluate / eval_and_free */
#include "../print.h"          /* expr_to_string — Thread::tdlen diagnostic */

struct CompiledFunction {
    unsigned         refcount;   /* Expr copies share the payload */
    size_t           nargs;
    const char**     arg_names;  /* interned pointers (interner owns the chars) */
    CompileType*     arg_types;
    Expr*            body;       /* owned; interpreter fallback + printing */
    CompiledProgram* prog;       /* NULL if body is outside the compilable subset */
    uint32_t         runtime_attrs; /* RuntimeAttributes: ATTR_LISTABLE or none */
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
static bool cf_box(const Expr* e, CompileType t, CompileValue* out, bool* packed) {
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
        default: break;
    }
    if (CT_IS_ARRAY(t)) {
        /* An NDArray argument is BORROWED — the program never frees it, and the
         * caller still owns the node it passed in.  A plain nested List is
         * packed here instead, and that temporary IS ours to free after the
         * call; `packed` tells the caller which is which. */
        if (e->type == EXPR_NDARRAY) {
            if (e->data.ndarray.rank != CT_RANK(t)) return false;
            out->v.a = (Expr*)e;
            return true;
        }
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_List) {
            Expr* nd = ndarray_from_nested_list(e, CT_ELEM(t) == CT_COMPLEX
                                                   ? NDT_COMPLEX64 : NDT_FLOAT64);
            if (!nd) return false;                       /* ragged or symbolic */
            if (nd->data.ndarray.rank != CT_RANK(t)) { expr_free(nd); return false; }
            out->v.a = nd;
            if (packed) *packed = true;
            return true;
        }
        return false;
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
        default: break;
    }
    /* compiled_eval hands back a NEW EXPR_NDARRAY that the caller owns, so it
     * becomes the result directly — no copy, no conversion to a nested List. */
    if (CT_IS_ARRAY(v->type)) return v->v.a;
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

/* Parse `{x, {y, _Real}, {v, _Real, 2}, ...}` into interned names + types.
 * Shared by Compile[] and CompileDiagnostics[] so the two can never disagree
 * about which signatures are legal.  On success the caller owns both arrays. */
static bool cf_parse_argspec(const Expr* argspec,
                             const char*** names_out, CompileType** types_out,
                             size_t* n_out) {
    if (!argspec) return false;
    if (argspec->type != EXPR_FUNCTION || argspec->data.function.head->type != EXPR_SYMBOL
        || argspec->data.function.head->data.symbol.name != SYM_List
        || argspec->data.function.arg_count == 0)
        return false;

    size_t n = argspec->data.function.arg_count;
    const char** names = malloc(n * sizeof(*names));
    CompileType* types = malloc(n * sizeof(*types));
    if (!names || !types) { free(names); free(types); return false; }

    for (size_t i = 0; i < n; i++) {
        const Expr* el = argspec->data.function.args[i];
        const char* nm = NULL; CompileType ty = CT_REAL;
        if (el->type == EXPR_SYMBOL) {
            nm = el->data.symbol.name;
        } else if (el->type == EXPR_FUNCTION && el->data.function.head->type == EXPR_SYMBOL
                   && el->data.function.head->data.symbol.name == SYM_List
                   && (el->data.function.arg_count == 2 || el->data.function.arg_count == 3)
                   && el->data.function.args[0]->type == EXPR_SYMBOL
                   && parse_typespec(el->data.function.args[1], &ty)) {
            /* `{v, _Real, r}` — a rank-r machine array, the same third-element
             * rank spelling Compile[] uses in the Wolfram Language.  The element
             * type must be Real or Complex: the packed buffer formats are
             * float/complex only, there is no integer dtype to hold a rank-r
             * Integer array. */
            if (el->data.function.arg_count == 3) {
                const Expr* rk = el->data.function.args[2];
                if (rk->type != EXPR_INTEGER) { free(names); free(types); return false; }
                long long r = (long long)rk->data.integer;
                if (r == 0) {
                    /* rank 0 is just the scalar */
                } else if (r < 1 || r > CT_MAX_RANK
                           || (ty != CT_REAL && ty != CT_COMPLEX)) {
                    free(names); free(types); return false;
                } else {
                    ty = CT_ARRAY(ty, (int)r);
                }
            }
            nm = el->data.function.args[0]->data.symbol.name;
        } else { free(names); free(types); return false; }

        const char* in = intern_symbol(nm);
        for (size_t j = 0; j < i; j++)
            if (names[j] == in) { free(names); free(types); return false; }  /* duplicate param */
        names[i] = in;
        types[i] = ty;
    }

    *names_out = names; *types_out = types; *n_out = n;
    return true;
}

CompiledFunction* compiled_function_new(const Expr* argspec, const Expr* body,
                                        uint32_t runtime_attrs) {
    if (!body) return NULL;
    const char** names; CompileType* types; size_t n;
    if (!cf_parse_argspec(argspec, &names, &types, &n)) return NULL;

    CompiledFunction* cf = calloc(1, sizeof *cf);
    if (!cf) { free(names); free(types); return NULL; }
    cf->refcount  = 1;
    cf->nargs     = n;
    cf->arg_names = names;
    cf->arg_types = types;
    cf->body      = expr_copy((Expr*)body);
    cf->prog      = compile_expr(cf->body, names, types, n);  /* NULL ⇒ fallback only */
    cf->runtime_attrs = runtime_attrs;
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

/* Argument count that boxes without touching the heap.  Compile[] signatures
 * are small; anything wider falls back to one allocation. */
#define CF_APPLY_STACK_ARGS 8

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

/* ---------------- RuntimeAttributes -> Listable threading ---------------- */

static bool cf_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
           && e->data.function.head->data.symbol.name == SYM_List;
}

/* Nesting depth of a List: 1 for `{…}`, 2 for `{{…}}`, … following the first
 * element at each level.  Only ever compared against a declared rank, so a
 * ragged argument needs no special case: it threads on the outer level and the
 * per-element call decides for itself what its own element is. */
static size_t cf_list_depth(const Expr* e) {
    size_t d = 0;
    while (cf_is_list(e)) {
        d++;
        if (e->data.function.arg_count == 0) break;
        e = e->data.function.args[0];
    }
    return d;
}

/* Does slot `i` hold an argument that this signature does NOT consume whole?  A
 * scalar parameter has rank 0, so any list threads; a rank-r array parameter
 * threads only over the levels ABOVE its own, so a Listable
 * `Compile[{{v, _Real, 1}}, …]` maps over the rows of a matrix and takes a plain
 * vector as one argument.
 *
 * A packed NDArray is measured by its rank for exactly the same reason it is
 * accepted as an argument at all: the object must answer the same for a List and
 * for the array it packs to. */
static bool cf_slot_threads(const CompiledFunction* cf, const Expr* a, size_t i) {
    size_t rank = CT_IS_ARRAY(cf->arg_types[i]) ? (size_t)CT_RANK(cf->arg_types[i]) : 0;
    if (cf_is_list(a)) return cf_list_depth(a) > rank;
    if (a && a->type == EXPR_NDARRAY) return (size_t)a->data.ndarray.rank > rank;
    return false;
}

/* Outer length of a threading argument. */
static size_t cf_thread_len(const Expr* a) {
    return cf_is_list(a) ? a->data.function.arg_count : (size_t)a->data.ndarray.dims[0];
}

/* Element `j` of a threading argument.  A List element is borrowed; an NDArray
 * row is a fresh sub-array (or scalar leaf) that the caller must free, which is
 * what `*owned` reports.  Returns NULL only if the array cannot be subscripted. */
static Expr* cf_thread_elem(const Expr* a, size_t j, bool* owned) {
    if (cf_is_list(a)) { *owned = false; return a->data.function.args[j]; }
    Expr* idx = expr_new_integer((int64_t)(j + 1));
    bool degrade = false;
    Expr* row = ndarray_part(a, &idx, 1, &degrade);
    expr_free(idx);
    *owned = true;
    return row;
}

/* Print Thread::tdlen against the real call, exactly as the evaluator's generic
 * Listable threading does for a symbol head. */
static void cf_report_tdlen(const CompiledFunction* cf, Expr* const* args, size_t nargs) {
    Expr** ca = malloc(nargs * sizeof(*ca));
    if (!ca) return;
    for (size_t i = 0; i < nargs; i++) ca[i] = expr_copy(args[i]);
    Expr* call = expr_new_function(
        expr_new_compiled(compiled_function_ref((CompiledFunction*)cf)), ca, nargs);
    free(ca);
    char* s = expr_to_string(call);
    printf("Thread::tdlen: Objects of unequal length in %s cannot be combined.\n",
           s ? s : "the compiled function");
    free(s);
    expr_free(call);
}

/* Thread a Listable object over its list-like arguments.  Returns false when
 * nothing threads (the caller proceeds with the ordinary application); returns
 * true with `*out` set to the threaded result, or to NULL when the threading
 * arguments disagree on length — the same "leave the application unevaluated
 * after a Thread::tdlen message" answer the evaluator gives a Listable symbol. */
static bool cf_thread(const CompiledFunction* cf, Expr* const* args, size_t nargs,
                      Expr** out) {
    bool any = false, from_nd = false;
    size_t len = 0;
    for (size_t i = 0; i < nargs; i++) {
        if (!cf_slot_threads(cf, args[i], i)) continue;
        size_t n = cf_thread_len(args[i]);
        if (args[i]->type == EXPR_NDARRAY) from_nd = true;
        if (!any) { any = true; len = n; }
        else if (n != len) { cf_report_tdlen(cf, args, nargs); *out = NULL; return true; }
    }
    if (!any) return false;

    /* Threading over empty lists yields an empty list, with no per-element work
     * (matching f[{}] -> {} for a Listable f). */
    Expr** items = len ? malloc(len * sizeof(*items)) : NULL;
    Expr** sub   = malloc(nargs * sizeof(*sub));       /* per-element arguments */
    bool*  owned = malloc(nargs * sizeof(*owned));     /* ... and which are ours */
    if ((len && !items) || !sub || !owned) {
        free(items); free(sub); free(owned); *out = NULL; return true;
    }

    for (size_t j = 0; j < len; j++) {
        bool ok = true;
        /* Cleared for EVERY slot before any is filled: the release loop below
         * walks all nargs entries, so a slot left over from the previous element
         * would be freed twice if this element gave up part-way through. */
        for (size_t i = 0; i < nargs; i++) { sub[i] = NULL; owned[i] = false; }
        for (size_t i = 0; i < nargs && ok; i++) {
            if (cf_slot_threads(cf, args[i], i)) {
                sub[i] = cf_thread_elem(args[i], j, &owned[i]);
                if (!sub[i]) ok = false;
            } else {
                sub[i] = args[i];
            }
        }
        /* Recursing through apply, not through the evaluator: the arity is
         * unchanged so this cannot fail, nested lists thread again on the way
         * down, and each leaf still picks its own bytecode-or-interpreter path. */
        items[j] = ok ? compiled_function_apply(cf, sub, nargs) : NULL;
        for (size_t i = 0; i < nargs; i++) if (owned[i] && sub[i]) expr_free(sub[i]);
        if (!items[j]) {
            for (size_t k = 0; k < j; k++) expr_free(items[k]);
            free(items); free(sub); free(owned);
            *out = NULL; return true;
        }
    }
    free(sub); free(owned);

    /* An NDArray went in, so an NDArray comes back — the same result-kind rule
     * the scalar path follows, and the reason a caller keeps state packed across
     * calls in the first place.  Elements that threaded a level further down are
     * already packed, so flatten them back to Lists first and pack once at the
     * top: that is what makes a rank-2 array over a scalar parameter come back as
     * one rank-2 array rather than a List of rows.  A result that will not pack
     * (Booleans, symbolic leaves from an interpreted element) stays a List. */
    if (from_nd) {
        for (size_t j = 0; j < len; j++) {
            if (items[j]->type != EXPR_NDARRAY) continue;
            Expr* nested = ndarray_to_nested_list(items[j]);
            if (nested) { expr_free(items[j]); items[j] = nested; }
        }
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), items, len);
    free(items);
    if (from_nd) {
        Expr* packed = ndarray_from_nested_list(r, NDT_FLOAT64);
        if (!packed) packed = ndarray_from_nested_list(r, NDT_COMPLEX64);
        if (packed) { expr_free(r); r = packed; }
    }
    *out = r;
    return true;
}

Expr* compiled_function_apply(const CompiledFunction* cf, Expr* const* args, size_t nargs) {
    if (!cf || nargs != cf->nargs) return NULL;   /* arity mismatch ⇒ leave unevaluated */

    /* Listable is a property of the OBJECT, not of its program: it threads
     * whether the body compiled or falls back, and it must run before boxing so
     * a List argument is never handed to a scalar parameter. */
    if (cf->runtime_attrs & ATTR_LISTABLE) {
        Expr* threaded;
        if (cf_thread(cf, args, nargs, &threaded)) return threaded;
    }

    if (cf->prog) {
        /* Small arities — every realistic one — box into stack storage.  This
         * path runs once per call of `cf[x]`, so a malloc/free pair here was
         * pure overhead on the hot path it exists to make fast. */
        CompileValue stackv[CF_APPLY_STACK_ARGS];
        CompileValue* cv = stackv;
        CompileValue* heap = NULL;
        if (nargs > CF_APPLY_STACK_ARGS) {
            heap = malloc(nargs * sizeof(*heap));
            if (!heap) return cf_fallback(cf, args, nargs);
            cv = heap;
        }
        /* A List argument packed into a temporary NDArray is ours to release;
         * an NDArray passed in directly is borrowed and must be left alone. */
        bool packed[CF_APPLY_STACK_ARGS];
        bool* packed_p = packed;
        bool* packed_heap = NULL;
        if (nargs > CF_APPLY_STACK_ARGS) {
            packed_heap = calloc(nargs, sizeof(bool));
            if (!packed_heap) { free(heap); return cf_fallback(cf, args, nargs); }
            packed_p = packed_heap;
        } else {
            for (size_t i = 0; i < nargs; i++) packed[i] = false;
        }

        bool all_numeric = true;
        for (size_t i = 0; i < nargs; i++)
            if (!cf_box(args[i], cf->arg_types[i], &cv[i], &packed_p[i])) { all_numeric = false; break; }
        #define CF_RELEASE_PACKED() do { \
            for (size_t q = 0; q < nargs; q++) \
                if (packed_p[q] && CT_IS_ARRAY(cf->arg_types[q])) expr_free(cv[q].v.a); \
            free(packed_heap); \
        } while (0)
        if (all_numeric) {
            /* An all-Real signature takes the unboxed entry point: no per-argument
             * type switch on the way in and no CompileValue on the way out. */
            if (compiled_program_all_real(cf->prog)) {
                double av[CF_APPLY_STACK_ARGS], o;
                if (nargs <= CF_APPLY_STACK_ARGS) {
                    for (size_t i = 0; i < nargs; i++) av[i] = cv[i].v.r;
                    bool ok = compiled_eval_real(cf->prog, av, &o);
                    CF_RELEASE_PACKED();
                    free(heap);
                    return ok ? expr_new_real(o) : cf_fallback(cf, args, nargs);
                }
            }
            CompileValue out;
            if (compiled_eval(cf->prog, cv, &out)) {
                Expr* r = cf_unbox(&out);
                /* The result KIND must follow the argument kind, or the compiled
                 * path would answer with a different head from the interpreter
                 * fallback for the very same input: given Lists the interpreter
                 * threads and returns a List, so an array result is unpacked
                 * back to one.  Given NDArrays it returns an NDArray, and so do
                 * we — no conversion, no copy.
                 *
                 * A body that BUILDS its array (ConstantArray, Table, NestList)
                 * has no kind to inherit: the interpreter returns a List however
                 * the arguments were spelled, because the construct has no
                 * packed form.  Asking only the arguments got that wrong for a
                 * body that takes an NDArray and builds its result from
                 * something else, so the PROGRAM is asked as well. */
                if (r && CT_IS_ARRAY(out.type)) {
                    bool from_nd = false;
                    for (size_t q = 0; q < nargs; q++)
                        if (CT_IS_ARRAY(cf->arg_types[q]) && !packed_p[q]) from_nd = true;
                    if (!from_nd || compiled_result_built(cf->prog)) {
                        Expr* lst = ndarray_to_nested_list(r);
                        if (lst) { expr_free(r); r = lst; }
                    }
                }
                CF_RELEASE_PACKED();
                free(heap);
                if (r) return r;
                return cf_fallback(cf, args, nargs);
            }
        }
        CF_RELEASE_PACKED();
        #undef CF_RELEASE_PACKED
        free(heap);
    }
    return cf_fallback(cf, args, nargs);
}

/* ------------------------------------------------------------------ *
 *  Accessors + builtin                                                *
 * ------------------------------------------------------------------ */

size_t             compiled_function_num_args(const CompiledFunction* cf) { return cf ? cf->nargs : 0; }
bool               compiled_function_is_compiled(const CompiledFunction* cf) { return cf && cf->prog; }
uint32_t           compiled_function_attributes(const CompiledFunction* cf) { return cf ? cf->runtime_attrs : ATTR_NONE; }
const Expr*        compiled_function_body(const CompiledFunction* cf) { return cf ? cf->body : NULL; }
const char* const* compiled_function_arg_names(const CompiledFunction* cf) { return cf ? cf->arg_names : NULL; }
const CompileType* compiled_function_arg_types(const CompiledFunction* cf) { return cf ? cf->arg_types : NULL; }
const CompiledProgram* compiled_function_program(const CompiledFunction* cf) { return cf ? cf->prog : NULL; }

/* A `RuntimeAttributes` setting as an attribute bitmask.  The accepted spellings
 * are `Listable`, `{Listable}` and the default `{}`; anything else returns false
 * so the caller can leave the call unevaluated instead of silently dropping a
 * setting it cannot honour. */
static bool cf_parse_runtime_attrs(const Expr* v, uint32_t* out) {
    if (!v) return false;
    if (v->type == EXPR_SYMBOL) {
        if (v->data.symbol.name != SYM_Listable) return false;
        *out = ATTR_LISTABLE;
        return true;
    }
    if (!cf_is_list(v)) return false;
    uint32_t a = ATTR_NONE;
    for (size_t i = 0; i < v->data.function.arg_count; i++) {
        const Expr* el = v->data.function.args[i];
        if (el->type != EXPR_SYMBOL || el->data.symbol.name != SYM_Listable) return false;
        a |= ATTR_LISTABLE;
    }
    *out = a;
    return true;
}

/* Is `e` a `RuntimeAttributes -> value` (or `:>`) option?  Returns true for the
 * option name regardless of the value; `*ok` says whether the value parsed. */
static bool cf_match_runtime_attrs_opt(const Expr* e, uint32_t* out, bool* ok) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.arg_count != 2)
        return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_RuntimeAttributes)
        return false;
    *ok = cf_parse_runtime_attrs(e->data.function.args[1], out);
    return true;
}

/* Compile[argspec, body, opts] (HoldAll).  Never evaluates the body; the raw
 * held body is compiled, and any non-arg symbol it references simply routes that
 * call through the interpreter fallback.
 *
 * The one option is `RuntimeAttributes`, whose only setting is `Listable`: the
 * resulting object then threads over List arguments (default `{}`). */
static Expr* builtin_compile(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    /* Seed from the registered defaults so SetOptions[Compile, ...] takes
     * effect, then let an explicit trailing option override below. */
    uint32_t rattrs = ATTR_NONE;
    Expr* defs = symtab_get_options("Compile");        /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            uint32_t v; bool ok = false;
            if (cf_match_runtime_attrs_opt(defs->data.function.args[i], &v, &ok) && ok)
                rattrs = v;
        }
    }

    /* Strip trailing options right to left, so the rightmost setting wins. An
     * option this function has no meaning for — an unknown name, or a
     * RuntimeAttributes value that is not a supported attribute — leaves
     * Compile[...] unevaluated rather than being quietly ignored. */
    bool explicit_attrs = false;
    while (argc > 2) {
        uint32_t v; bool ok = false;
        if (!cf_match_runtime_attrs_opt(a[argc - 1], &v, &ok) || !ok) return NULL;
        if (!explicit_attrs) { rattrs = v; explicit_attrs = true; }
        argc--;
    }
    if (argc != 2) return NULL;

    CompiledFunction* cf = compiled_function_new(a[0], a[1], rattrs);
    if (!cf) return NULL;   /* malformed argspec ⇒ leave Compile[...] unevaluated */
    return expr_new_compiled(cf);
}

/* ------------------------------------------------------------------ *
 *  CompileDiagnostics                                                  *
 * ------------------------------------------------------------------ */

static const char* ct_name(CompileType t) {
    switch (t) {
        case CT_BOOL:    return "Boolean";
        case CT_INT:     return "Integer";
        case CT_REAL:    return "Real";
        case CT_COMPLEX: return "Complex";
        default: break;
    }
    return CT_IS_ARRAY(t) ? "Array" : "Unknown";
}

static Expr* diag_rule(const char* key, Expr* value) {
    Expr* a[2] = { expr_new_string(key), value };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* CompileDiagnostics[argspec, body] (HoldAll) — did this body compile, and if
 * not, what stopped it?
 *
 * The engine's central usability problem is that a bail is silent: the caller
 * interprets, the answer stays correct, and the only symptom is being an order
 * of magnitude slower.  Worse, the compilable subset is a cliff — one head
 * outside it costs the whole body — so "why is this slow" is almost always
 * "which single subexpression bailed", a question nothing could answer before.
 *
 * On success it also reports the instruction count with and without the
 * optimiser, which is the cheapest way to see what code generation did. */
static Expr* builtin_compile_diagnostics(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* argspec = res->data.function.args[0];
    const Expr* body    = res->data.function.args[1];

    const char** names; CompileType* types; size_t n;
    if (!cf_parse_argspec(argspec, &names, &types, &n)) return NULL;

    CompiledProgram* p = compile_expr(body, names, types, n);

    Expr* items[6];
    size_t ni = 0;
    if (p) {
        items[ni++] = diag_rule("Compiled", expr_new_symbol("True"));
        items[ni++] = diag_rule("ResultType",
                                expr_new_string(ct_name(compiled_result_type(p))));
        items[ni++] = diag_rule("Instructions",
                                expr_new_integer((int64_t)compiled_num_instructions(p)));
        items[ni++] = diag_rule("CommonSubexpressions",
                                expr_new_integer((int64_t)compiled_num_cse(p)));
        /* Recompiling with the optimiser off is the honest way to report what it
         * removed: the passes rewrite in place, so the pre-pass count is gone by
         * the time anyone can ask for it. */
        CompiledProgram* raw = compile_expr_ex(body, names, types, n, COMPILE_NO_OPT);
        if (raw) {
            items[ni++] = diag_rule("InstructionsUnoptimized",
                                    expr_new_integer((int64_t)compiled_num_instructions(raw)));
            compiled_free(raw);
        }
        compiled_free(p);
    } else {
        const char* why  = compiled_bail_reason();
        const char* what = compiled_bail_expr();
        items[ni++] = diag_rule("Compiled", expr_new_symbol("False"));
        items[ni++] = diag_rule("Reason", expr_new_string(why ? why : "unknown"));
        if (what) items[ni++] = diag_rule("Subexpression", expr_new_string(what));
    }
    free(names); free(types);
    return expr_new_function(expr_new_symbol(SYM_List), items, ni);
}

/* ------------------------------------------------------------------ *
 *  CompilePrint                                                        *
 * ------------------------------------------------------------------ */

/* CompilePrint[cf] — print the bytecode of a CompiledFunction.
 *
 * NOT HoldAll, unlike Compile and CompileDiagnostics: the argument has to
 * evaluate down to the EXPR_COMPILED atom, so that both CompilePrint[Compile[
 * ...]] and `f = Compile[...]; CompilePrint[f]` work.  Anything else is left
 * unevaluated rather than reported as an error, the same way the rest of the
 * system treats an argument it has no meaning for. */
static Expr* builtin_compile_print(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* a = res->data.function.args[0];
    if (a->type != EXPR_COMPILED) return NULL;
    char* s = compiled_function_disassemble(a->data.compiled);
    if (!s) return NULL;
    printf("%s", s);
    free(s);
    return expr_new_symbol(SYM_Null);
}

void compiled_function_init(void) {
    symtab_add_builtin("Compile", builtin_compile);
    SymbolDef* d = symtab_get_def("Compile");
    if (d) d->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Compile",
        "Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a "
        "CompiledFunction that evaluates expr over machine numbers (types _Real, "
        "_Integer, _Complex; default _Real), falling back to the interpreter for "
        "symbolic arguments or non-compilable bodies. With "
        "RuntimeAttributes -> Listable the object threads over List arguments; "
        "the default is RuntimeAttributes -> {}.");

    symtab_add_builtin("CompileDiagnostics", builtin_compile_diagnostics);
    SymbolDef* dd = symtab_get_def("CompileDiagnostics");
    if (dd) dd->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("CompileDiagnostics",
        "CompileDiagnostics[argspec, expr] reports whether expr compiles for the "
        "given Compile[] argument specification, and if not, the innermost "
        "subexpression that could not be lowered. For a compiled body it also "
        "gives the result type and the instruction count with and without the "
        "optimiser.");

    symtab_add_builtin("CompilePrint", builtin_compile_print);
    SymbolDef* dp = symtab_get_def("CompilePrint");
    if (dp) dp->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CompilePrint",
        "CompilePrint[cf] prints the bytecode of the CompiledFunction cf: its "
        "argument and result registers with their types, the scalar/array/tile "
        "register banks, and one line per instruction giving both the raw "
        "operands and a readable rendering. For an object whose body did not "
        "compile it reports the bail reason instead. Returns Null.");
}
