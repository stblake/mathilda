/*
 * numloop.c -- Automatic numeric loop fast-path (see numloop.h).
 *
 * Compiles a numeric-closed arithmetic Expr body into a small stack-machine
 * bytecode over `double`, then runs the enclosing loop entirely in doubles with
 * no per-iteration Expr allocation. Gated hard on the correctness contract in
 * numloop.h: only machine-real, inexact-result computations take this path.
 */
#include "numloop.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "numeric.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Master enable switch. Defaults on; the env var MATHILDA_NO_NUMLOOP (read
 * once) or numloop_set_enabled() force every entry point to fall through to
 * the interpreter -- used by the differential test to compare the two paths
 * bit-for-bit. */
static int g_numloop_enabled = -1;   /* -1 = uninitialised */
void numloop_set_enabled(bool on) { g_numloop_enabled = on ? 1 : 0; }
static bool numloop_off(void) {
    if (g_numloop_enabled < 0)
        g_numloop_enabled = getenv("MATHILDA_NO_NUMLOOP") ? 0 : 1;
    return g_numloop_enabled == 0;
}

/* ------------------------------------------------------------------------
 *  Bytecode
 * ---------------------------------------------------------------------- */
typedef enum {
    OP_CONST,   /* push consts[a]            */
    OP_VAR,     /* push regs[a]              */
    OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_NEG,
    OP_POW,     /* pop e,b -> push pow(b,e)  */
    OP_SIN, OP_COS, OP_TAN,
    OP_SINH, OP_COSH, OP_TANH,
    OP_EXP, OP_LOG, OP_SQRT, OP_ABS, OP_ARCTAN,
    OP_LOAD     /* pop 1-based index -> push arr[a][index] (NaN if out of range) */
} NumOp;

typedef struct { uint8_t op; int32_t a; } NInsn;

/* Array context for OP_LOAD (indexed reads of NDArray buffers inside a Part
 * loop). NULL for scalar programs, which never emit OP_LOAD. OP_LOAD pops
 * `rank` 1-based indices (i1..iR, deepest on top), validates each against
 * `dims`, and pushes the row-major element (NaN if any index is out of range). */
typedef struct {
    double*        buf[4];    /* flat float64 buffers */
    int            rank[4];   /* rank of each array */
    const int64_t* dims[4];   /* dims of each array */
    size_t         count;
} ArrCtx;

typedef struct {
    NInsn*  code;   size_t ncode, ccode;
    double* consts; size_t nconsts, cconsts;
    size_t  nvars;
    size_t  max_stack;
    size_t  cur_stack;   /* compile-time bookkeeping */
    bool    ok;
} NumProg;

static void prog_init(NumProg* p, size_t nvars) {
    memset(p, 0, sizeof(*p));
    p->nvars = nvars;
    p->ok = true;
}

static void prog_free(NumProg* p) {
    free(p->code);
    free(p->consts);
}

static void prog_push_depth(NumProg* p, int delta) {
    if (delta > 0) {
        p->cur_stack += (size_t)delta;
        if (p->cur_stack > p->max_stack) p->max_stack = p->cur_stack;
    } else {
        p->cur_stack -= (size_t)(-delta);
    }
}

static void emit(NumProg* p, NumOp op, int32_t a, int stack_delta) {
    if (!p->ok) return;
    if (p->ncode == p->ccode) {
        size_t nc = p->ccode ? p->ccode * 2 : 16;
        NInsn* g = realloc(p->code, nc * sizeof(NInsn));
        if (!g) { p->ok = false; return; }
        p->code = g; p->ccode = nc;
    }
    p->code[p->ncode].op = (uint8_t)op;
    p->code[p->ncode].a  = a;
    p->ncode++;
    prog_push_depth(p, stack_delta);
}

static int32_t add_const(NumProg* p, double v) {
    if (p->nconsts == p->cconsts) {
        size_t nc = p->cconsts ? p->cconsts * 2 : 8;
        double* g = realloc(p->consts, nc * sizeof(double));
        if (!g) { p->ok = false; return 0; }
        p->consts = g; p->cconsts = nc;
    }
    p->consts[p->nconsts] = v;
    return (int32_t)p->nconsts++;
}

/* ------------------------------------------------------------------------
 *  Variable resolution + leaf helpers
 * ---------------------------------------------------------------------- */

/* When slot_var is set the single variable is Slot[1] (index 0); otherwise the
 * variables are named symbols whose interned names sit in var_names[]. */
typedef struct {
    const char** var_names;   /* interned names, length nvars (named mode) */
    size_t       nvars;
    bool         slot_var;    /* Slot[1] -> var 0 (pure-function mode) */
    const bool*  defined;     /* optional length-nvars mask: only vars marked
                                 defined resolve; an undefined read then falls to
                                 const-fold and bails. NULL = all defined. */
    const char** arr_names;   /* array-variable names: Part[name, i1..iR] on one
                                 of these compiles to R index pushes + an OP_LOAD
                                 instead of const-folding. NULL = none. */
    const int*   arr_rank;    /* rank of each array var (indices expected) */
    size_t       narr;
} VarCtx;

/* Resolve a bare symbol to an array-variable index, or -1. */
static int resolve_arr(const VarCtx* vc, const Expr* e) {
    if (!vc->arr_names || e->type != EXPR_SYMBOL) return -1;
    for (size_t i = 0; i < vc->narr; i++)
        if (e->data.symbol.name == vc->arr_names[i]) return (int)i;
    return -1;
}

/* Return k for a Slot[k] node (k >= 1), or -1 if `e` is not a numbered Slot. */
static int slot_index(const Expr* e) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Slot
        && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER)
        return (int)e->data.function.args[0]->data.integer;
    return -1;
}

/* Resolve `e` to a variable index, or -1 if it is not one of the loop vars.
 * In slot mode Slot[k] maps to var k-1 (k in 1..nvars). */
static int resolve_var(const VarCtx* vc, const Expr* e) {
    if (vc->slot_var) {
        int k = slot_index(e);
        if (k >= 1 && (size_t)k <= vc->nvars) return k - 1;
        return -1;
    }
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < vc->nvars; i++)
            if (e->data.symbol.name == vc->var_names[i]) {
                if (vc->defined && !vc->defined[i]) return -1;   /* read-before-def */
                return (int)i;
            }
    }
    return -1;
}

/* Does `e` reference any loop variable (scalar or indexed array read) anywhere
 * in its tree? An array element load is loop-varying, so it counts even when its
 * index is constant -- otherwise a constant-index read would be frozen. */
static bool contains_var(const VarCtx* vc, const Expr* e) {
    if (resolve_var(vc, e) >= 0) return true;
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol.name == SYM_Part &&
            e->data.function.arg_count >= 2 &&
            resolve_arr(vc, e->data.function.args[0]) >= 0)
            return true;
        if (contains_var(vc, e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_var(vc, e->data.function.args[i])) return true;
    }
    return false;
}

/* Coerce an already-numericalized value to a machine double. Accepts Integer,
 * Real, BigInt, and Rational[num,den]; rejects Complex/MPFR/symbolic. */
static bool to_machine_double(const Expr* v, double* out) {
    if (!v) return false;
    switch (v->type) {
        case EXPR_INTEGER: *out = (double)v->data.integer; return true;
        case EXPR_REAL:    *out = v->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(v->data.bigint); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational((Expr*)v, &n, &d) && d != 0) {
        *out = (double)n / (double)d;
        return true;
    }
    return false;
}

/* Numericalize a variable-free subexpression to a machine double. */
static bool const_fold(const Expr* e, double* out) {
    Expr* v = numericalize(e, numeric_machine_spec());
    bool ok = to_machine_double(v, out);
    expr_free(v);
    return ok;
}

/* Above this element count the interpreter's vectorized NDArray kernels
 * (tight C / BLAS loops) beat per-element bytecode interpretation, so the fused
 * fast-path only earns its keep for small arrays in tight loops -- where the
 * per-iteration evaluator + allocation overhead the fast-path removes dominates.
 * Larger arrays decline it and use the interpreter's vectorized path. */
#define NUMLOOP_ARRAY_MAX_ELEMS 256

/* True for a dense real (float64) NDArray -- the dtype whose flat buffer is a
 * plain double[] the scalar VM can iterate element-by-element. Complex/float32
 * arrays are left to the interpreter. */
static bool is_f64_ndarray(const Expr* e) {
    return e && e->type == EXPR_NDARRAY && e->data.ndarray.dtype == NDT_FLOAT64;
}

/* Element count of an NDArray (product of its dims). */
static size_t nd_elem_count(const Expr* e) {
    size_t n = 1;
    for (int i = 0; i < e->data.ndarray.rank; i++)
        n *= (size_t)e->data.ndarray.dims[i];
    return n;
}

/* True if `e` carries a machine-inexact leaf (Real / MPFR, incl. inside
 * Complex[...]). Presence of one anywhere in the body guarantees the
 * interpreter's result is inexact, which is what makes the double result
 * authoritative even from an exact seed. */
static bool has_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION) {
        if (has_inexact_leaf(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_inexact_leaf(e->data.function.args[i])) return true;
    }
    return false;
}

/* ------------------------------------------------------------------------
 *  Compiler
 * ---------------------------------------------------------------------- */

static void compile_walk(NumProg* p, const VarCtx* vc, const Expr* e);

/* Compile a variadic Plus/Times as a left fold of binary ops. */
static void compile_fold(NumProg* p, const VarCtx* vc, const Expr* e,
                         NumOp binop, double identity) {
    size_t n = e->data.function.arg_count;
    if (n == 0) { emit(p, OP_CONST, add_const(p, identity), +1); return; }
    compile_walk(p, vc, e->data.function.args[0]);
    for (size_t i = 1; i < n; i++) {
        compile_walk(p, vc, e->data.function.args[i]);
        emit(p, binop, 0, -1);   /* pop 2, push 1 */
    }
}

static NumOp unary_op_for(const char* head) {
    if (head == SYM_Sin)    return OP_SIN;
    if (head == SYM_Cos)    return OP_COS;
    if (head == SYM_Tan)    return OP_TAN;
    if (head == SYM_Sinh)   return OP_SINH;
    if (head == SYM_Cosh)   return OP_COSH;
    if (head == SYM_Tanh)   return OP_TANH;
    if (head == SYM_Exp)    return OP_EXP;
    if (head == SYM_Log)    return OP_LOG;
    if (head == SYM_Sqrt)   return OP_SQRT;
    if (head == SYM_Abs)    return OP_ABS;
    if (head == SYM_ArcTan) return OP_ARCTAN;
    return (NumOp)255;   /* not a supported unary head */
}

static void compile_walk(NumProg* p, const VarCtx* vc, const Expr* e) {
    if (!p->ok) return;

    /* Indexed array read Part[arr, i1..iR]. Handled before the const-fold check
     * because the array's elements change across iterations, so it must NOT be
     * frozen to a constant. The number of indices must match the array's rank. */
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Part &&
        e->data.function.arg_count >= 2) {
        int ai = resolve_arr(vc, e->data.function.args[0]);
        if (ai >= 0) {
            size_t nidx = e->data.function.arg_count - 1;
            if ((int)nidx != vc->arr_rank[ai]) { p->ok = false; return; }
            for (size_t j = 0; j < nidx; j++)
                compile_walk(p, vc, e->data.function.args[j + 1]);  /* push i_{j+1} */
            emit(p, OP_LOAD, ai, -(int)(nidx - 1));   /* pop R, push 1 */
            return;
        }
    }

    /* Variable reference (bare symbol or Slot[1]). */
    int vi = resolve_var(vc, e);
    if (vi >= 0) { emit(p, OP_VAR, vi, +1); return; }

    /* Any variable-free subexpression collapses to a single constant. This
     * folds Pi, E, Rational[p,q], Sqrt[2], and plain literals uniformly. */
    if (!contains_var(vc, e)) {
        double c;
        if (!const_fold(e, &c)) { p->ok = false; return; }
        emit(p, OP_CONST, add_const(p, c), +1);
        return;
    }

    /* Structural nodes that mix in a variable. */
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) {
        p->ok = false;
        return;
    }
    const char* head = e->data.function.head->data.symbol.name;
    size_t argc = e->data.function.arg_count;

    if (head == SYM_Plus)  { compile_fold(p, vc, e, OP_ADD, 0.0); return; }
    if (head == SYM_Times) { compile_fold(p, vc, e, OP_MUL, 1.0); return; }

    if (head == SYM_Power && argc == 2) {
        compile_walk(p, vc, e->data.function.args[0]);   /* base */
        compile_walk(p, vc, e->data.function.args[1]);   /* exponent */
        emit(p, OP_POW, 0, -1);                          /* pop 2, push 1 */
        return;
    }

    if (argc == 1) {
        /* N[x] is the identity on a value already carried as a machine double. */
        if (head == SYM_N) { compile_walk(p, vc, e->data.function.args[0]); return; }
        NumOp uop = unary_op_for(head);
        if (uop != (NumOp)255) {
            compile_walk(p, vc, e->data.function.args[0]);
            emit(p, uop, 0, 0);   /* pop 1, push 1 */
            return;
        }
    }

    p->ok = false;   /* unsupported head -> not compilable */
}

/* Compile `body` over the given variables. On success returns true and the
 * program is ready to run; the caller must prog_free it. */
static bool numprog_compile(NumProg* p, const Expr* body, const VarCtx* vc) {
    prog_init(p, vc->nvars);
    compile_walk(p, vc, body);
    if (!p->ok || p->cur_stack != 1) {   /* must net exactly one result */
        prog_free(p);
        return false;
    }
    return true;
}

/* Compile a callable `f` applied to `arity` machine-double arguments (register
 * indices 0..arity-1) into `p`. Handles a pure Function (slot body, named
 * single parameter, or named parameter list) and a bare unary numeric function
 * head such as Cos. *body_inexact reports whether the body carries a Real
 * literal (always false for a bare head -- its inexactness must come from the
 * seed). Returns true on success; the caller must prog_free(p). */
static bool compile_function(NumProg* p, const Expr* f, size_t arity,
                             bool* body_inexact) {
    *body_inexact = false;

    /* Bare numeric head, e.g. Nest[Cos, x0, n] / FixedPoint[Cos, x0]. */
    if (f->type == EXPR_SYMBOL) {
        if (arity != 1) return false;
        NumOp uop = unary_op_for(f->data.symbol.name);
        if (uop == (NumOp)255) return false;
        prog_init(p, 1);
        emit(p, OP_VAR, 0, +1);
        emit(p, uop, 0, 0);
        if (!p->ok || p->cur_stack != 1) { prog_free(p); return false; }
        return true;
    }

    if (f->type != EXPR_FUNCTION ||
        f->data.function.head->type != EXPR_SYMBOL ||
        f->data.function.head->data.symbol.name != SYM_Function)
        return false;

    size_t fargc = f->data.function.arg_count;
    const Expr* body = NULL;
    const char* names[8];
    VarCtx vc;

    if (fargc == 1) {                                   /* Function[body] */
        body = f->data.function.args[0];
        vc.var_names = NULL; vc.nvars = arity; vc.slot_var = true;
    } else if (fargc >= 2 && f->data.function.args[0]->type == EXPR_SYMBOL &&
               f->data.function.args[0]->data.symbol.name == SYM_Null) {
        body = f->data.function.args[1];                /* Function[Null, body, ...] */
        vc.var_names = NULL; vc.nvars = arity; vc.slot_var = true;
    } else if (fargc == 2 && f->data.function.args[0]->type == EXPR_SYMBOL) {
        if (arity != 1) return false;                   /* Function[x, body] */
        names[0] = f->data.function.args[0]->data.symbol.name;
        body = f->data.function.args[1];
        vc.var_names = names; vc.nvars = 1; vc.slot_var = false;
    } else if (fargc == 2 && f->data.function.args[0]->type == EXPR_FUNCTION &&
               f->data.function.args[0]->data.function.head->type == EXPR_SYMBOL &&
               f->data.function.args[0]->data.function.head->data.symbol.name == SYM_List) {
        const Expr* pl = f->data.function.args[0];      /* Function[{x1,...}, body] */
        size_t k = pl->data.function.arg_count;
        if (k != arity || k > 8) return false;
        for (size_t i = 0; i < k; i++) {
            if (pl->data.function.args[i]->type != EXPR_SYMBOL) return false;
            names[i] = pl->data.function.args[i]->data.symbol.name;
        }
        body = f->data.function.args[1];
        vc.var_names = names; vc.nvars = k; vc.slot_var = false;
    } else {
        return false;
    }

    *body_inexact = has_inexact_leaf(body);
    return numprog_compile(p, body, &vc);
}

/* ------------------------------------------------------------------------
 *  Runtime VM
 * ---------------------------------------------------------------------- */
static double numprog_run_ac(const NumProg* p, const double* regs,
                             double* stack, const ArrCtx* ac) {
    size_t sp = 0;
    const NInsn* c = p->code;
    size_t n = p->ncode;
    for (size_t i = 0; i < n; i++) {
        switch ((NumOp)c[i].op) {
            case OP_CONST: stack[sp++] = p->consts[c[i].a]; break;
            case OP_VAR:   stack[sp++] = regs[c[i].a];      break;
            case OP_ADD:   stack[sp-2] = stack[sp-2] + stack[sp-1]; sp--; break;
            case OP_SUB:   stack[sp-2] = stack[sp-2] - stack[sp-1]; sp--; break;
            case OP_MUL:   stack[sp-2] = stack[sp-2] * stack[sp-1]; sp--; break;
            case OP_DIV:   stack[sp-2] = stack[sp-2] / stack[sp-1]; sp--; break;
            case OP_POW:   stack[sp-2] = pow(stack[sp-2], stack[sp-1]); sp--; break;
            case OP_NEG:   stack[sp-1] = -stack[sp-1];      break;
            case OP_SIN:   stack[sp-1] = sin(stack[sp-1]);  break;
            case OP_COS:   stack[sp-1] = cos(stack[sp-1]);  break;
            case OP_TAN:   stack[sp-1] = tan(stack[sp-1]);  break;
            case OP_SINH:  stack[sp-1] = sinh(stack[sp-1]); break;
            case OP_COSH:  stack[sp-1] = cosh(stack[sp-1]); break;
            case OP_TANH:  stack[sp-1] = tanh(stack[sp-1]); break;
            case OP_EXP:   stack[sp-1] = exp(stack[sp-1]);  break;
            case OP_LOG:   stack[sp-1] = log(stack[sp-1]);  break;
            case OP_SQRT:  stack[sp-1] = sqrt(stack[sp-1]); break;
            case OP_ABS:   stack[sp-1] = fabs(stack[sp-1]); break;
            case OP_ARCTAN:stack[sp-1] = atan(stack[sp-1]); break;
            case OP_LOAD: {
                /* pop `rank` 1-based indices (i1..iR, iR on top), validate each
                 * against dims, push the row-major element; any out-of-range or
                 * non-integer index -> NaN, which the caller's isfinite check
                 * turns into a fall-back to the interpreter. */
                size_t av = (size_t)c[i].a;
                int R = ac->rank[av];
                int64_t off = 0;
                bool oob = false;
                for (int d = 0; d < R; d++) {
                    double xd = stack[sp - R + d];
                    int64_t id = (int64_t)xd;
                    if ((double)id != xd || id < 1 || id > ac->dims[av][d]) { oob = true; break; }
                    off = off * ac->dims[av][d] + (id - 1);
                }
                sp -= (size_t)R;
                stack[sp++] = oob ? (double)NAN : ac->buf[av][off];
                break;
            }
        }
    }
    return stack[0];
}

static double numprog_run(const NumProg* p, const double* regs, double* stack) {
    return numprog_run_ac(p, regs, stack, NULL);
}

/* ------------------------------------------------------------------------
 *  Seed / writeback helpers
 * ---------------------------------------------------------------------- */

/* Read the current OwnValue of a bare symbol and coerce it to a machine
 * double. Returns false if unbound or non-real. *inexact is set true when the
 * bound value is itself machine-inexact (Real/MPFR), which -- like a Real
 * literal in the body -- guarantees the interpreter's result is inexact. */
static bool seed_from_symbol(const Expr* sym, double* out, bool* inexact) {
    Expr* cur = evaluate((Expr*)sym);   /* evaluate borrows, returns fresh */
    bool ok = false;
    if (cur && cur != sym && cur->type != EXPR_SYMBOL) {
        Expr* v = numericalize(cur, numeric_machine_spec());
        ok = to_machine_double(v, out);
        if (ok) *inexact = has_inexact_leaf(cur);
        expr_free(v);
    }
    expr_free(cur);
    return ok;
}

/* Set symbol `sym`'s OwnValue to `value` (adopts nothing: add_rule copies). */
static void writeback_symbol(const Expr* sym, Expr* value) {
    symtab_add_own_value(sym->data.symbol.name, (Expr*)sym, value);
    expr_free(value);   /* add_rule kept its own copy */
}

/* Evaluate a held numeric expression to a machine double (used for For bounds
 * and start values). Returns false if not a machine number. */
static bool eval_to_double(const Expr* e, double* out) {
    Expr* v = eval_and_free(expr_copy((Expr*)e));
    Expr* nv = numericalize(v, numeric_machine_spec());
    bool ok = to_machine_double(nv, out);
    expr_free(nv);
    expr_free(v);
    return ok;
}

/* Evaluate a held expression to an exact int64, if it is one. */
static bool eval_to_int(const Expr* e, int64_t* out) {
    Expr* v = eval_and_free(expr_copy((Expr*)e));
    bool ok = (v && v->type == EXPR_INTEGER);
    if (ok) *out = v->data.integer;
    expr_free(v);
    return ok;
}

/* ------------------------------------------------------------------------
 *  Nest[f, x0, n]
 * ---------------------------------------------------------------------- */
/* True when a seed value forces the interpreter's result inexact by itself. */
static bool value_is_inexact(const Expr* v) {
    if (v->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (v->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Nest[f, arr, n] over a dense float64 NDArray: the compiled scalar body is run
 * per element (element-local, so the update is safe in place), fusing the whole
 * map with no intermediate array temporaries and no per-iteration Expr/NDArray
 * allocation. The body may reference Slot[1] and scalar constants only. */
static Expr* numloop_nest_array(const Expr* f, const Expr* x0, int64_t n) {
    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;
    if (p.nvars != 1) { prog_free(&p); return NULL; }   /* single register per element */

    size_t N = nd_elem_count(x0);
    if (N == 0 || N > NUMLOOP_ARRAY_MAX_ELEMS) { prog_free(&p); return NULL; }
    double* buf = malloc(N * sizeof(double));
    if (!buf) { prog_free(&p); return NULL; }
    memcpy(buf, x0->data.ndarray.data, N * sizeof(double));

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { free(buf); prog_free(&p); return NULL; }

    bool bail = false;
    for (int64_t i = 0; i < n && !bail; i++) {
        for (size_t k = 0; k < N; k++) {
            double v = numprog_run(&p, &buf[k], stack);
            if (!isfinite(v)) { bail = true; break; }
            buf[k] = v;
        }
    }
    free(stack);
    prog_free(&p);
    if (bail) { free(buf); return NULL; }   /* interpreter re-runs the whole Nest */
    return expr_new_ndarray(x0->data.ndarray.rank, x0->data.ndarray.dims,
                            buf, NDT_FLOAT64);   /* takes ownership of buf */
}

Expr* numloop_nest(const Expr* f, const Expr* x0, int64_t n) {
    if (numloop_off()) return NULL;
    if (n < 0) return NULL;

    if (is_f64_ndarray(x0)) return numloop_nest_array(f, x0, n);

    double x;
    if (!to_machine_double(x0, &x)) return NULL;

    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;

    /* Result must be inexact: seed already Real/MPFR, or the body carries a
     * Real literal. Otherwise the interpreter would stay exact/symbolic. */
    if (!value_is_inexact(x0) && !body_inexact) { prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { prog_free(&p); return NULL; }

    bool bail = false;
    for (int64_t i = 0; i < n; i++) {
        x = numprog_run(&p, &x, stack);
        if (!isfinite(x)) { bail = true; break; }
    }
    free(stack);
    prog_free(&p);
    if (bail) return NULL;   /* interpreter re-runs the whole Nest */
    return expr_new_real(x);
}

/* ========================================================================
 *  Imperative loop bodies as numeric assignment blocks
 *
 *  Do / For / While bodies compile as a *block* of assignments -- a single
 *  Set[v, e], or a CompoundExpression Set[v1,e1]; Set[v2,e2]; ... All mutated
 *  variables share one double register file and statements run in order each
 *  iteration, so later statements observe earlier updates (e.g.
 *  `x = 3.5 x (1-x); y = 4 x; x = y/4`). Read-only symbols with a fixed real
 *  value fold to constants; a symbolic operand, a non-Set statement, or a
 *  variable read before it is assigned makes compilation bail to the
 *  interpreter.
 * ==================================================================== */

#define NUMBLOCK_MAXVARS  16
#define NUMBLOCK_MAXSTMTS 32

typedef struct {
    const char* names[NUMBLOCK_MAXVARS];
    const Expr* syms[NUMBLOCK_MAXVARS];   /* bare-symbol Expr, for writeback */
    bool   assigned[NUMBLOCK_MAXVARS];    /* written by some statement */
    bool   defined[NUMBLOCK_MAXVARS];     /* has a value at the current point */
    bool   seeded[NUMBLOCK_MAXVARS];      /* had a machine-real initial value */
    double regs[NUMBLOCK_MAXVARS];
    size_t nvars;
    int    counter_idx;                   /* loop-counter var index, or -1 */

    NumProg progs[NUMBLOCK_MAXSTMTS];
    int     lhs[NUMBLOCK_MAXSTMTS];
    size_t  nstmts;
    size_t  max_stack;
    bool    forces_real;

    /* Array mode (all block variables are same-shape float64 NDArrays): the
     * compiled scalar bytecode is run per element over these flat buffers,
     * fusing the whole element-wise map with no intermediate array temporaries. */
    bool     is_array;
    double*  abuf[NUMBLOCK_MAXVARS];       /* owned float64 buffer per var */
    size_t   nelem;                        /* element count (all vars) */
    int      arr_rank;
    int64_t  arr_dims[8];                  /* shape for writeback */
} NumBlock;

static void numblock_free(NumBlock* b) {
    for (size_t i = 0; i < b->nstmts; i++) prog_free(&b->progs[i]);
    b->nstmts = 0;
    if (b->is_array)
        for (size_t i = 0; i < b->nvars; i++) { free(b->abuf[i]); b->abuf[i] = NULL; }
}

static int nb_var(NumBlock* b, const char* name, const Expr* sym) {
    for (size_t i = 0; i < b->nvars; i++)
        if (b->names[i] == name) return (int)i;
    if (b->nvars >= NUMBLOCK_MAXVARS) return -1;
    size_t i = b->nvars++;
    b->names[i] = name; b->syms[i] = sym;
    b->assigned[i] = b->defined[i] = b->seeded[i] = false;
    b->regs[i] = 0.0;
    return (int)i;
}

/* Recognise Set[sym, rhs] with sym a bare symbol. */
static bool stmt_is_set(const Expr* s, const Expr** sym, const Expr** rhs) {
    if (s->type != EXPR_FUNCTION ||
        s->data.function.head->type != EXPR_SYMBOL ||
        s->data.function.head->data.symbol.name != SYM_Set ||
        s->data.function.arg_count != 2 ||
        s->data.function.args[0]->type != EXPR_SYMBOL)
        return false;
    *sym = s->data.function.args[0];
    *rhs = s->data.function.args[1];
    return true;
}

/* Build a compiled block from an imperative body. counter_sym (or NULL) names a
 * loop-managed variable seeded to counter_seed; it may be read but not assigned
 * by the block. Returns true on success (caller must numblock_free). */
static bool numblock_build(NumBlock* b, const Expr* body,
                           const Expr* counter_sym, double counter_seed) {
    memset(b, 0, sizeof(*b));
    b->counter_idx = -1;

    if (counter_sym) {
        b->counter_idx = nb_var(b, counter_sym->data.symbol.name, counter_sym);
        b->defined[b->counter_idx] = b->seeded[b->counter_idx] = true;
        b->regs[b->counter_idx] = counter_seed;
    }

    /* Flatten the body into statements (a CompoundExpression, else a single). */
    const Expr* stmts[NUMBLOCK_MAXSTMTS];
    size_t ns = 0;
    if (body->type == EXPR_FUNCTION &&
        body->data.function.head->type == EXPR_SYMBOL &&
        body->data.function.head->data.symbol.name == SYM_CompoundExpression) {
        for (size_t i = 0; i < body->data.function.arg_count; i++) {
            const Expr* a = body->data.function.args[i];
            if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Null) continue;
            if (ns >= NUMBLOCK_MAXSTMTS) return false;
            stmts[ns++] = a;
        }
    } else {
        stmts[ns++] = body;
    }
    if (ns == 0) return false;

    /* Pass 1: every statement must be an assignment; register its LHS variable. */
    const Expr* rhs_of[NUMBLOCK_MAXSTMTS];
    int         lhs_of[NUMBLOCK_MAXSTMTS];
    for (size_t i = 0; i < ns; i++) {
        const Expr *sym, *rhs;
        if (!stmt_is_set(stmts[i], &sym, &rhs)) return false;
        int vi = nb_var(b, sym->data.symbol.name, sym);
        if (vi < 0 || vi == b->counter_idx) return false;  /* table full / assigns counter */
        b->assigned[vi] = true;
        lhs_of[i] = vi;
        rhs_of[i] = rhs;
    }

    /* Seed assigned variables from their current OwnValues (the loop-carried
     * initial value). A var that fails to seed stays undefined until it is first
     * assigned; reading it earlier makes pass 2 bail. */
    for (size_t i = 0; i < b->nvars; i++) {
        if ((int)i == b->counter_idx) continue;
        double v; bool inexact = false;
        if (seed_from_symbol(b->syms[i], &v, &inexact)) {
            b->regs[i] = v;
            b->defined[i] = b->seeded[i] = true;
            if (inexact) b->forces_real = true;
        }
    }

    /* Pass 2: compile each statement in order over the currently-defined vars.
     * A read of a not-yet-defined variable does not resolve, so it const-folds
     * and bails. After compiling, the LHS becomes defined for later statements. */
    VarCtx vc = { .var_names = b->names, .nvars = b->nvars,
                  .slot_var = false, .defined = b->defined };
    for (size_t i = 0; i < ns; i++) {
        if (has_inexact_leaf(rhs_of[i])) b->forces_real = true;
        if (!numprog_compile(&b->progs[i], rhs_of[i], &vc)) {
            b->nstmts = i;             /* only [0,i) were compiled */
            numblock_free(b);
            return false;
        }
        if (b->progs[i].max_stack > b->max_stack) b->max_stack = b->progs[i].max_stack;
        b->lhs[i] = lhs_of[i];
        b->defined[lhs_of[i]] = true;
    }
    b->nstmts = ns;
    if (b->max_stack == 0) b->max_stack = 1;
    return true;
}

/* Run one pass of the block, updating regs in place. Returns false on a
 * non-finite result (the caller then bails to the interpreter). */
static bool numblock_step(NumBlock* b, double* stack) {
    for (size_t i = 0; i < b->nstmts; i++) {
        double v = numprog_run(&b->progs[i], b->regs, stack);
        if (!isfinite(v)) return false;
        b->regs[b->lhs[i]] = v;
    }
    return true;
}

/* Write every assigned variable's final value back to its OwnValue as a Real. */
static void numblock_writeback(NumBlock* b) {
    for (size_t i = 0; i < b->nvars; i++)
        if (b->assigned[i])
            writeback_symbol(b->syms[i], expr_new_real(b->regs[i]));
}

/* ---- Array-mode block: all variables are same-shape float64 NDArrays ---- */

/* Build an array block from an imperative body (no loop counter). Every LHS
 * variable must currently hold a float64 NDArray of one common shape (or be an
 * assigned-before-read temporary, which gets a zero buffer of that shape).
 * Read-only operands must be scalar constants. Returns false (so the caller can
 * try the scalar path or the interpreter) when the body is not such a loop. */
static bool numblock_build_array(NumBlock* b, const Expr* body) {
    memset(b, 0, sizeof(*b));
    b->counter_idx = -1;

    const Expr* stmts[NUMBLOCK_MAXSTMTS];
    size_t ns = 0;
    if (body->type == EXPR_FUNCTION &&
        body->data.function.head->type == EXPR_SYMBOL &&
        body->data.function.head->data.symbol.name == SYM_CompoundExpression) {
        for (size_t i = 0; i < body->data.function.arg_count; i++) {
            const Expr* a = body->data.function.args[i];
            if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Null) continue;
            if (ns >= NUMBLOCK_MAXSTMTS) return false;
            stmts[ns++] = a;
        }
    } else {
        stmts[ns++] = body;
    }
    if (ns == 0) return false;

    const Expr* rhs_of[NUMBLOCK_MAXSTMTS];
    int         lhs_of[NUMBLOCK_MAXSTMTS];
    for (size_t i = 0; i < ns; i++) {
        const Expr *sym, *rhs;
        if (!stmt_is_set(stmts[i], &sym, &rhs)) return false;
        int vi = nb_var(b, sym->data.symbol.name, sym);
        if (vi < 0) return false;
        b->assigned[vi] = true;
        lhs_of[i] = vi; rhs_of[i] = rhs;
    }

    /* Seed variables from their current values, establishing the common shape.
     * At least one variable must be a small float64 NDArray for this to be an
     * array loop; the rest must match its shape (or be undefined temporaries). */
    bool have_shape = false;
    for (size_t i = 0; i < b->nvars; i++) {
        Expr* cur = evaluate((Expr*)b->syms[i]);
        if (is_f64_ndarray(cur)) {
            size_t Ni = nd_elem_count(cur);
            if (!have_shape) {
                if (Ni == 0 || Ni > NUMLOOP_ARRAY_MAX_ELEMS ||
                    cur->data.ndarray.rank > 8) { expr_free(cur); goto fail; }
                b->nelem = Ni;
                b->arr_rank = cur->data.ndarray.rank;
                memcpy(b->arr_dims, cur->data.ndarray.dims,
                       sizeof(int64_t) * (size_t)b->arr_rank);
                have_shape = true;
            } else if (Ni != b->nelem) { expr_free(cur); goto fail; }
            b->abuf[i] = malloc(b->nelem * sizeof(double));
            if (!b->abuf[i]) { expr_free(cur); goto fail; }
            memcpy(b->abuf[i], cur->data.ndarray.data, b->nelem * sizeof(double));
            b->defined[i] = b->seeded[i] = true;
        } else if (cur && cur != b->syms[i] && cur->type != EXPR_SYMBOL) {
            /* a bound non-array value (e.g. a scalar) cannot mix into an array
             * block -- decline so the scalar path / interpreter handles it. */
            expr_free(cur);
            goto fail;
        }
        /* else: unbound -> an assigned-before-read temporary (buffer below). */
        expr_free(cur);
    }
    if (!have_shape) goto fail;   /* no array variable -> not an array loop */

    /* Allocate zero buffers for undefined temporaries now that the shape is known. */
    for (size_t i = 0; i < b->nvars; i++) {
        if (!b->abuf[i]) {
            b->abuf[i] = calloc(b->nelem, sizeof(double));
            if (!b->abuf[i]) goto fail;
        }
    }

    b->is_array = true;
    b->forces_real = true;   /* NDArrays are always inexact */

    /* Compile each statement over the array variables (a read-only array operand
     * doesn't resolve as a var and won't const-fold, so it bails here). */
    VarCtx vc = { .var_names = b->names, .nvars = b->nvars,
                  .slot_var = false, .defined = b->defined };
    for (size_t i = 0; i < ns; i++) {
        if (!numprog_compile(&b->progs[i], rhs_of[i], &vc)) {
            b->nstmts = i;
            numblock_free(b);
            return false;
        }
        if (b->progs[i].max_stack > b->max_stack) b->max_stack = b->progs[i].max_stack;
        b->lhs[i] = lhs_of[i];
        b->defined[lhs_of[i]] = true;
    }
    b->nstmts = ns;
    if (b->max_stack == 0) b->max_stack = 1;
    return true;

fail:
    b->is_array = true;   /* so numblock_free releases any abuf already taken */
    numblock_free(b);
    return false;
}

/* Run one pass of the array block: fuse all statements over one traversal of the
 * element index, updating buffers in place (element-local, so safe). Returns
 * false on a non-finite result. `elem` is a scratch register file of nvars. */
static bool numblock_step_array(NumBlock* b, double* elem, double* stack) {
    for (size_t k = 0; k < b->nelem; k++) {
        for (size_t i = 0; i < b->nvars; i++) elem[i] = b->abuf[i][k];
        for (size_t s = 0; s < b->nstmts; s++) {
            double v = numprog_run(&b->progs[s], elem, stack);
            if (!isfinite(v)) return false;
            elem[b->lhs[s]] = v;
            b->abuf[b->lhs[s]][k] = v;
        }
    }
    return true;
}

/* Write each assigned variable's final buffer back as a float64 NDArray. The
 * buffer's ownership transfers to the new NDArray, so it is not freed here. */
static void numblock_writeback_array(NumBlock* b) {
    for (size_t i = 0; i < b->nvars; i++) {
        if (!b->assigned[i]) continue;
        Expr* nd = expr_new_ndarray(b->arr_rank, b->arr_dims, b->abuf[i], NDT_FLOAT64);
        b->abuf[i] = NULL;   /* ownership moved into nd */
        writeback_symbol(b->syms[i], nd);
    }
}

/* ========================================================================
 *  Part-assignment loops:  Do[a[[idx]] = rhs, {i, ...}] / For[...]
 *
 *  A counter-driven loop whose body writes one element of a 1-D float64 NDArray
 *  per iteration. The buffer is mutated *in place* (O(iterations)), avoiding the
 *  interpreter's whole-array copy per Part-set (which is O(iterations * N)). The
 *  rhs may read the counter, other elements a[[jexpr]] (OP_LOAD), and scalar
 *  constants. Single Set[Part[a, idx], rhs] statement.
 * ==================================================================== */
typedef struct {
    const Expr* arr_sym;        /* the array's bare symbol, for writeback */
    double*     buf;            /* owned float64 buffer, mutated in place */
    int64_t     nelem;
    int         rank;
    int64_t     dims[8];
    NumProg     idx_prog[8];    /* one LHS index expression per axis */
    NumProg     rhs_prog;       /* rhs, over {counter} + array reads */
    size_t      max_stack;
    bool        built;
} PartLoop;

static void partloop_free(PartLoop* pl) {
    if (!pl->built) return;
    for (int k = 0; k < pl->rank; k++) prog_free(&pl->idx_prog[k]);
    prog_free(&pl->rhs_prog);
    free(pl->buf);
    pl->built = false;
}

/* Build from a single-statement body Set[Part[a, i1..iR], rhs] where a is a bare
 * symbol bound to a rank-R float64 NDArray; counter_name is the loop variable
 * (scalar register 0). The rhs may read a[[j1..jR]] and the counter. */
static bool partloop_build(PartLoop* pl, const Expr* body, const char* counter_name) {
    memset(pl, 0, sizeof(*pl));
    if (body->type != EXPR_FUNCTION ||
        body->data.function.head->type != EXPR_SYMBOL ||
        body->data.function.head->data.symbol.name != SYM_Set ||
        body->data.function.arg_count != 2)
        return false;
    const Expr* lhs = body->data.function.args[0];
    const Expr* rhs = body->data.function.args[1];
    if (lhs->type != EXPR_FUNCTION ||
        lhs->data.function.head->type != EXPR_SYMBOL ||
        lhs->data.function.head->data.symbol.name != SYM_Part ||
        lhs->data.function.arg_count < 2 ||
        lhs->data.function.args[0]->type != EXPR_SYMBOL)
        return false;
    const Expr* asym = lhs->data.function.args[0];
    size_t nidx = lhs->data.function.arg_count - 1;

    Expr* cur = evaluate((Expr*)asym);
    if (!is_f64_ndarray(cur) || cur->data.ndarray.rank < 1 ||
        cur->data.ndarray.rank > 8 || (size_t)cur->data.ndarray.rank != nidx) {
        expr_free(cur); return false;
    }
    pl->arr_sym = asym;
    pl->rank = cur->data.ndarray.rank;
    pl->nelem = 1;
    for (int k = 0; k < pl->rank; k++) {
        pl->dims[k] = cur->data.ndarray.dims[k];
        pl->nelem *= pl->dims[k];
    }
    pl->buf = malloc((size_t)pl->nelem * sizeof(double));
    if (!pl->buf) { expr_free(cur); return false; }
    memcpy(pl->buf, cur->data.ndarray.data, (size_t)pl->nelem * sizeof(double));
    expr_free(cur);

    const char* cvn = counter_name;
    const char* avn = asym->data.symbol.name;
    int arank = pl->rank;
    VarCtx vc = { .var_names = &cvn, .nvars = 1, .slot_var = false, .defined = NULL,
                  .arr_names = &avn, .arr_rank = &arank, .narr = 1 };

    pl->max_stack = 0;
    for (size_t k = 0; k < nidx; k++) {
        if (!numprog_compile(&pl->idx_prog[k], lhs->data.function.args[k + 1], &vc)) {
            for (size_t j = 0; j < k; j++) prog_free(&pl->idx_prog[j]);
            free(pl->buf); return false;
        }
        if (pl->idx_prog[k].max_stack > pl->max_stack) pl->max_stack = pl->idx_prog[k].max_stack;
    }
    if (!numprog_compile(&pl->rhs_prog, rhs, &vc)) {
        for (size_t j = 0; j < nidx; j++) prog_free(&pl->idx_prog[j]);
        free(pl->buf); return false;
    }
    if (pl->rhs_prog.max_stack > pl->max_stack) pl->max_stack = pl->rhs_prog.max_stack;
    pl->built = true;
    return true;
}

/* Run one iteration at counter value `i`. Evaluates each LHS axis index,
 * validates it against the array's shape, and stores the rhs at the row-major
 * offset. Returns false (bail) on an out-of-range index or a non-finite rhs. */
static bool partloop_step(PartLoop* pl, int64_t i, double* regs, double* stack,
                          const ArrCtx* ac) {
    regs[0] = (double)i;
    int64_t off = 0;
    for (int k = 0; k < pl->rank; k++) {
        double xk = numprog_run_ac(&pl->idx_prog[k], regs, stack, ac);
        int64_t ik = (int64_t)xk;
        if ((double)ik != xk || ik < 1 || ik > pl->dims[k]) return false;
        off = off * pl->dims[k] + (ik - 1);
    }
    double v = numprog_run_ac(&pl->rhs_prog, regs, stack, ac);
    if (!isfinite(v)) return false;
    pl->buf[off] = v;
    return true;
}

/* Write the mutated buffer back as a float64 NDArray (ownership transfers). */
static void partloop_writeback(PartLoop* pl) {
    Expr* nd = expr_new_ndarray(pl->rank, pl->dims, pl->buf, NDT_FLOAT64);
    pl->buf = NULL;   /* ownership moved into nd */
    for (int k = 0; k < pl->rank; k++) prog_free(&pl->idx_prog[k]);
    prog_free(&pl->rhs_prog);
    pl->built = false;
    writeback_symbol(pl->arr_sym, nd);
}

/* ------------------------------------------------------------------------
 *  Do[body, {n}] count form  /  Do[body, {i, imin, imax, di}] range form
 * ---------------------------------------------------------------------- */
Expr* numloop_do_count(const Expr* body, int64_t n) {
    if (numloop_off()) return NULL;
    if (n < 1) return NULL;

    NumBlock b;

    /* Small-float64-NDArray body: fuse the element-wise map over the flat
     * buffers (no per-iteration Expr/NDArray allocation). */
    if (numblock_build_array(&b, body)) {
        double* stack = malloc(b.max_stack * sizeof(double));
        double* elem  = malloc(b.nvars * sizeof(double));
        if (!stack || !elem) { free(stack); free(elem); numblock_free(&b); return NULL; }
        bool abail = false;
        for (int64_t k = 0; k < n; k++)
            if (!numblock_step_array(&b, elem, stack)) { abail = true; break; }
        free(stack); free(elem);
        if (abail) { numblock_free(&b); return NULL; }
        numblock_writeback_array(&b);
        numblock_free(&b);
        return expr_new_symbol(SYM_Null);
    }

    if (!numblock_build(&b, body, NULL, 0.0)) return NULL;
    if (!b.forces_real) { numblock_free(&b); return NULL; }

    double* stack = malloc(b.max_stack * sizeof(double));
    if (!stack) { numblock_free(&b); return NULL; }

    bool bail = false;
    for (int64_t k = 0; k < n; k++)
        if (!numblock_step(&b, stack)) { bail = true; break; }
    free(stack);
    if (bail) { numblock_free(&b); return NULL; }   /* vars untouched; interp re-runs */

    numblock_writeback(&b);
    numblock_free(&b);
    return expr_new_symbol(SYM_Null);
}

Expr* numloop_do_range(const Expr* body, const Expr* var,
                       int64_t imin, int64_t imax, int64_t di) {
    if (numloop_off()) return NULL;
    if (di == 0) return NULL;

    /* In-place Part-assignment loop: Do[a[[idx]] = rhs, {i, imin, imax, di}]. */
    {
        PartLoop pl;
        if (partloop_build(&pl, body, var->data.symbol.name)) {
            double* stack = malloc(pl.max_stack * sizeof(double));
            if (!stack) { partloop_free(&pl); return NULL; }
            ArrCtx ac = { .count = 1 };
            ac.buf[0] = pl.buf; ac.rank[0] = pl.rank; ac.dims[0] = pl.dims;
            double regs[1];
            bool bail = false;
            for (int64_t i = imin; (di > 0) ? (i <= imax) : (i >= imax); i += di)
                if (!partloop_step(&pl, i, regs, stack, &ac)) { bail = true; break; }
            free(stack);
            if (bail) { partloop_free(&pl); return NULL; }
            partloop_writeback(&pl);   /* iterator stays localised; only `a` persists */
            return expr_new_symbol(SYM_Null);
        }
    }

    NumBlock b;
    if (!numblock_build(&b, body, var, (double)imin)) return NULL;
    if (!b.forces_real) { numblock_free(&b); return NULL; }

    double* stack = malloc(b.max_stack * sizeof(double));
    if (!stack) { numblock_free(&b); return NULL; }

    bool bail = false;
    for (int64_t i = imin; (di > 0) ? (i <= imax) : (i >= imax); i += di) {
        b.regs[b.counter_idx] = (double)i;
        if (!numblock_step(&b, stack)) { bail = true; break; }
    }
    free(stack);
    if (bail) { numblock_free(&b); return NULL; }

    /* Do localises its iterator: we never touched var's OwnValue, so it is
     * already restored. Only the block's assigned variables persist. */
    numblock_writeback(&b);
    numblock_free(&b);
    return expr_new_symbol(SYM_Null);
}

/* ------------------------------------------------------------------------
 *  For[start, test, incr, body]
 * ---------------------------------------------------------------------- */

/* Classify a two-argument comparison head; op: 0 '<', 1 '<=', 2 '>', 3 '>='. */
static bool cmp_op(const char* head, int* op) {
    if (head == SYM_Less)         { *op = 0; return true; }
    if (head == SYM_LessEqual)    { *op = 1; return true; }
    if (head == SYM_Greater)      { *op = 2; return true; }
    if (head == SYM_GreaterEqual) { *op = 3; return true; }
    return false;
}

static bool cmp_eval(double a, double b, int op) {
    switch (op) {
        case 0: return a <  b;
        case 1: return a <= b;
        case 2: return a >  b;
        default:return a >= b;
    }
}

Expr* numloop_for(const Expr* start, const Expr* test,
                  const Expr* incr, const Expr* body) {
    if (numloop_off()) return NULL;

    /* start = Set[i, <int>] */
    const Expr *ivar, *istart_rhs;
    if (!stmt_is_set(start, &ivar, &istart_rhs)) return NULL;
    int64_t i0;
    if (!eval_to_int(istart_rhs, &i0)) return NULL;

    /* incr = Increment[i] on the same counter */
    if (incr->type != EXPR_FUNCTION ||
        incr->data.function.head->type != EXPR_SYMBOL ||
        incr->data.function.head->data.symbol.name != SYM_Increment ||
        incr->data.function.arg_count != 1 ||
        incr->data.function.args[0]->type != EXPR_SYMBOL ||
        incr->data.function.args[0]->data.symbol.name != ivar->data.symbol.name)
        return NULL;

    /* test = i <cmp> <bound>, counter on the left. */
    if (test->type != EXPR_FUNCTION ||
        test->data.function.head->type != EXPR_SYMBOL ||
        test->data.function.arg_count != 2)
        return NULL;
    int op;
    if (!cmp_op(test->data.function.head->data.symbol.name, &op)) return NULL;
    if (test->data.function.args[0]->type != EXPR_SYMBOL ||
        test->data.function.args[0]->data.symbol.name != ivar->data.symbol.name)
        return NULL;
    double bound;
    if (!eval_to_double(test->data.function.args[1], &bound)) return NULL;

    /* In-place Part-assignment loop: For[i=i0, i<n, i++, a[[idx]] = rhs]. */
    {
        PartLoop pl;
        if (partloop_build(&pl, body, ivar->data.symbol.name)) {
            double* stack = malloc(pl.max_stack * sizeof(double));
            if (!stack) { partloop_free(&pl); return NULL; }
            ArrCtx ac = { .count = 1 };
            ac.buf[0] = pl.buf; ac.rank[0] = pl.rank; ac.dims[0] = pl.dims;
            double regs[1];
            bool pbail = false;
            int64_t pi = i0;
            while (cmp_eval((double)pi, bound, op)) {
                if (!partloop_step(&pl, pi, regs, stack, &ac)) { pbail = true; break; }
                pi++;
            }
            free(stack);
            if (pbail) { partloop_free(&pl); return NULL; }
            partloop_writeback(&pl);
            writeback_symbol(ivar, expr_new_integer(pi));   /* For keeps its counter */
            return expr_new_symbol(SYM_Null);
        }
    }

    NumBlock b;
    if (!numblock_build(&b, body, ivar, (double)i0)) return NULL;
    if (!b.forces_real) { numblock_free(&b); return NULL; }

    double* stack = malloc(b.max_stack * sizeof(double));
    if (!stack) { numblock_free(&b); return NULL; }

    bool bail = false;
    int64_t i = i0;
    while (cmp_eval((double)i, bound, op)) {
        b.regs[b.counter_idx] = (double)i;
        if (!numblock_step(&b, stack)) { bail = true; break; }
        i++;   /* Increment[i] */
    }
    free(stack);
    if (bail) { numblock_free(&b); return NULL; }

    /* For does not localise its counter: leave i at its final integer value and
     * the block's assigned variables at their final real values. */
    numblock_writeback(&b);
    writeback_symbol(ivar, expr_new_integer(i));
    numblock_free(&b);
    return expr_new_symbol(SYM_Null);
}

/* ------------------------------------------------------------------------
 *  While[test, body]
 * ---------------------------------------------------------------------- */
Expr* numloop_while(const Expr* test, const Expr* body) {
    if (numloop_off()) return NULL;

    /* test = <lhs> <cmp> <rhs> */
    if (test->type != EXPR_FUNCTION ||
        test->data.function.head->type != EXPR_SYMBOL ||
        test->data.function.arg_count != 2)
        return NULL;
    int op;
    if (!cmp_op(test->data.function.head->data.symbol.name, &op)) return NULL;

    NumBlock b;
    if (!numblock_build(&b, body, NULL, 0.0)) return NULL;
    if (!b.forces_real && !has_inexact_leaf(test)) { numblock_free(&b); return NULL; }

    /* The test runs before the body each iteration, so it may read only
     * variables that already hold a value (seeded) -- gate its compile on the
     * seeded mask, not the post-body defined mask. */
    VarCtx tvc = { .var_names = b.names, .nvars = b.nvars,
                   .slot_var = false, .defined = b.seeded };
    NumProg tl, tr;
    if (!numprog_compile(&tl, test->data.function.args[0], &tvc)) { numblock_free(&b); return NULL; }
    if (!numprog_compile(&tr, test->data.function.args[1], &tvc)) {
        prog_free(&tl); numblock_free(&b); return NULL;
    }

    size_t ms = b.max_stack;
    if (tl.max_stack > ms) ms = tl.max_stack;
    if (tr.max_stack > ms) ms = tr.max_stack;
    double* stack = malloc(ms * sizeof(double));
    if (!stack) { prog_free(&tl); prog_free(&tr); numblock_free(&b); return NULL; }

    bool bail = false;
    int64_t guard = 0;
    const int64_t GUARD_CAP = 100000000;   /* runaway backstop */
    while (cmp_eval(numprog_run(&tl, b.regs, stack),
                    numprog_run(&tr, b.regs, stack), op)) {
        if (!numblock_step(&b, stack)) { bail = true; break; }
        if (++guard > GUARD_CAP) { bail = true; break; }
    }
    free(stack);
    prog_free(&tl); prog_free(&tr);
    if (bail) { numblock_free(&b); return NULL; }

    numblock_writeback(&b);
    numblock_free(&b);
    return expr_new_symbol(SYM_Null);
}

/* ------------------------------------------------------------------------
 *  Fold[f, x0, list]
 * ---------------------------------------------------------------------- */
Expr* numloop_fold(const Expr* f, const Expr* x0, const Expr* list) {
    if (numloop_off()) return NULL;
    if (list->type != EXPR_FUNCTION) return NULL;
    size_t m = list->data.function.arg_count;
    if (m == 0) return NULL;   /* Fold[f,x0,{}] = x0 (possibly exact); let interp */

    double acc;
    if (!to_machine_double(x0, &acc)) return NULL;

    /* f is binary: Slot[1]=accumulator, Slot[2]=list element. */
    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 2, &body_inexact)) return NULL;

    /* Every list element must already be a machine number; gather them and note
     * whether any is inexact (which -- like a Real seed or Real body literal --
     * forces the interpreter's result inexact). */
    double* elems = malloc(m * sizeof(double));
    if (!elems) { prog_free(&p); return NULL; }
    bool elem_inexact = false;
    bool ok = true;
    for (size_t i = 0; i < m; i++) {
        const Expr* el = list->data.function.args[i];
        if (!to_machine_double(el, &elems[i])) { ok = false; break; }
        if (value_is_inexact(el)) elem_inexact = true;
    }
    if (!ok) { free(elems); prog_free(&p); return NULL; }

    if (!value_is_inexact(x0) && !body_inexact && !elem_inexact) {
        free(elems); prog_free(&p); return NULL;
    }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { free(elems); prog_free(&p); return NULL; }

    bool bail = false;
    for (size_t i = 0; i < m; i++) {
        double regs[2] = { acc, elems[i] };
        acc = numprog_run(&p, regs, stack);
        if (!isfinite(acc)) { bail = true; break; }
    }
    free(stack);
    free(elems);
    prog_free(&p);
    if (bail) return NULL;
    return expr_new_real(acc);
}

/* ------------------------------------------------------------------------
 *  FixedPoint[f, x0]  (default SameTest, no application cap)
 * ---------------------------------------------------------------------- */
Expr* numloop_fixedpoint(const Expr* f, const Expr* x0) {
    if (numloop_off()) return NULL;

    double x;
    if (!to_machine_double(x0, &x)) return NULL;

    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;
    if (!value_is_inexact(x0) && !body_inexact) { prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { prog_free(&p); return NULL; }

    /* Iterate x := f(x) until it stops changing (SameQ on machine reals is
     * exact double equality, matching the interpreter's expr_eq). A non-
     * converging orbit hits the safety cap and bails to the interpreter. */
    bool bail = false;
    int64_t guard = 0;
    const int64_t CAP = 1000000;   /* == ITER_SAFETY_CAP */
    double cur = x;
    while (1) {
        double next = numprog_run(&p, &cur, stack);
        if (!isfinite(next)) { bail = true; break; }
        bool same = (next == cur);
        cur = next;
        if (same) break;
        if (++guard >= CAP) { bail = true; break; }
    }
    free(stack);
    prog_free(&p);
    if (bail) return NULL;
    return expr_new_real(cur);
}

/* ------------------------------------------------------------------------
 *  NestWhile[f, x0, test]  (m = 1, default max / n; test a unary predicate)
 * ---------------------------------------------------------------------- */
Expr* numloop_nestwhile(const Expr* f, const Expr* x0, const Expr* test) {
    if (numloop_off()) return NULL;

    double x;
    if (!to_machine_double(x0, &x)) return NULL;

    /* test must be a pure Function whose body is a comparison of two
     * numeric-closed operands in Slot[1], e.g. (# < 100 &). */
    const Expr* tbody = NULL;
    if (test->type == EXPR_FUNCTION &&
        test->data.function.head->type == EXPR_SYMBOL &&
        test->data.function.head->data.symbol.name == SYM_Function) {
        if (test->data.function.arg_count == 1)
            tbody = test->data.function.args[0];
        else if (test->data.function.arg_count >= 2 &&
                 test->data.function.args[0]->type == EXPR_SYMBOL &&
                 test->data.function.args[0]->data.symbol.name == SYM_Null)
            tbody = test->data.function.args[1];
    }
    if (!tbody || tbody->type != EXPR_FUNCTION ||
        tbody->data.function.head->type != EXPR_SYMBOL ||
        tbody->data.function.arg_count != 2)
        return NULL;
    int op;
    if (!cmp_op(tbody->data.function.head->data.symbol.name, &op)) return NULL;

    NumProg pf, tl, tr;
    bool body_inexact;
    if (!compile_function(&pf, f, 1, &body_inexact)) return NULL;
    if (!value_is_inexact(x0) && !body_inexact) { prog_free(&pf); return NULL; }

    VarCtx vc = { .var_names = NULL, .nvars = 1, .slot_var = true };
    if (!numprog_compile(&tl, tbody->data.function.args[0], &vc)) { prog_free(&pf); return NULL; }
    if (!numprog_compile(&tr, tbody->data.function.args[1], &vc)) {
        prog_free(&pf); prog_free(&tl); return NULL;
    }

    size_t ms = pf.max_stack;
    if (tl.max_stack > ms) ms = tl.max_stack;
    if (tr.max_stack > ms) ms = tr.max_stack;
    double* stack = malloc(ms * sizeof(double));
    if (!stack) { prog_free(&pf); prog_free(&tl); prog_free(&tr); return NULL; }

    /* NestWhile: while test(current) holds, apply f; stop at the first value
     * that fails the test and return it. */
    bool bail = false;
    int64_t guard = 0;
    const int64_t CAP = 1000000;
    double cur = x;
    while (cmp_eval(numprog_run(&tl, &cur, stack),
                    numprog_run(&tr, &cur, stack), op)) {
        cur = numprog_run(&pf, &cur, stack);
        if (!isfinite(cur)) { bail = true; break; }
        if (++guard >= CAP) { bail = true; break; }
    }
    free(stack);
    prog_free(&pf); prog_free(&tl); prog_free(&tr);
    if (bail) return NULL;
    return expr_new_real(cur);
}
