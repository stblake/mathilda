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
#include "pack.h"
#include "ndarray.h"   /* is_packed_list — Map over a packed List */
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

/* Gate for the variable-free const-fold below.
 *
 * const_fold() *evaluates* whatever it is handed (numericalize walks the tree
 * and applies the evaluator to each node), so handing it an arbitrary
 * subexpression is a speculative evaluation of user code. That is wrong twice
 * over:
 *
 *   - Side effects fire an extra time. `f[] := (c = c + 1; 1.5)` in
 *     `Do[y = f[], {5}]` incremented `c` six times: once for the probe, then
 *     five times for the loop the probe failed to replace.
 *   - The probe is not even the same computation. numericalize() rewrites
 *     exact integers to machine reals *before* evaluating, so a call like
 *     `istep[grid, grid, lam, 41]` was evaluated with 41. in the argument
 *     that ends up as a `Table` bound and a `Part` subscript. Real subscripts
 *     do not resolve, so the body collapsed into one huge symbolic `Plus`
 *     whose like-term hashing is O(grid) *per element* -- the probe ran ~4x
 *     slower than the real evaluation it was trying to pre-empt, and its
 *     result was then discarded.
 *
 * So fold only what is syntactically closed over numeric literals, numeric
 * constants and the arithmetic/elementary heads this block already lowers.
 * Anything else bails to the interpreter, which is always a correct outcome
 * for a fast path -- numloop is an optimization, never a semantic. */
static bool const_foldable(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        /* A bare symbol is foldable when reading it is a pure lookup: either it
         * has no OwnValue at all (Pi, E, Degree -- numericalize resolves the
         * constant, and const_fold rejects anything else), or its OwnValue is
         * already a number, as `x = 2.5` leaves it.
         *
         * A *delayed* OwnValue (`x := (c = c + 1; 1.5)`) is a held expression
         * re-run on every read, so folding it both fires its side effects once
         * and freezes the first value as a loop constant -- `Do[y = x + 1., {5}]`
         * evaluated the body once instead of five times. Reject anything whose
         * stored value is not already a number. */
        case EXPR_SYMBOL: {
            SymbolDef* def = symtab_lookup(e->data.symbol.name);
            if (!def || !def->own_values) return true;
            const Expr* v = def->own_values->replacement;
            if (!v) return false;
            if (v->type == EXPR_INTEGER || v->type == EXPR_REAL ||
                v->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
            if (v->type == EXPR_MPFR) return true;
#endif
            {
                int64_t rn, rd;
                return is_rational((Expr*)v, &rn, &rd);
            }
        }
        case EXPR_FUNCTION:
            break;
        default:
            return false;
    }

    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* nm = h->data.symbol.name;
    if (!(nm == SYM_Plus || nm == SYM_Times || nm == SYM_Power ||
          nm == SYM_Subtract || nm == SYM_Divide || nm == SYM_Rational ||
          nm == SYM_N || unary_op_for(nm) != (NumOp)255))
        return false;

    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!const_foldable(e->data.function.args[i])) return false;
    return true;
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
     * folds Pi, E, Rational[p,q], Sqrt[2], and plain literals uniformly --
     * but only once const_foldable() has confirmed that folding it cannot
     * run a user-defined rule (see the note on that predicate). */
    if (!contains_var(vc, e)) {
        double c;
        if (!const_foldable(e))    { p->ok = false; return; }
        if (!const_fold(e, &c))    { p->ok = false; return; }
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
    /* Zero-initialised, not merely assigned field-by-field below: the branches
     * set var_names/nvars/slot_var and nothing else, so `defined`, `arr_names`,
     * `arr_rank` and `narr` would keep whatever was on the stack. Both optional
     * pointers are read as "NULL means absent" -- resolve_var dereferences
     * `defined[i]` the moment `defined` is non-NULL -- so stack litter there is a
     * wild read, not a wrong answer. It stayed invisible for as long as every
     * caller happened to leave zeros at those offsets; routing one more head
     * (FixedPointList) through here changed the stack shape and it faulted. */
    VarCtx vc = {0};

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

/* ------------------------------------------------------------------------
 *  Collecting the iterates of a *List head
 *
 *  NestList / FoldList / NestWhileList / FixedPointList differ from their
 *  scalar twins only in keeping every iterate instead of the last. Both stages
 *  below collect raw doubles and box to Expr only once the loop has finished
 *  successfully, so a bail -- the common case for a body that leaves the
 *  compilable subset or goes non-finite midway -- costs no Expr allocation at
 *  all and has nothing to unwind. (numloop_map predates this and boxes as it
 *  goes, which is why it carries a `done` counter to free the prefix.)
 * ---------------------------------------------------------------------- */

/* Box n doubles as List[Real...]. Every value here is a machine real by
 * construction, so a big enough result goes straight into a packed buffer -- one
 * memcpy instead of n Expr allocations. The four callers (NestList, FoldList and
 * the two growable-accumulator loops) all reach this with a plain array. */
static Expr* reals_to_list(const double* v, size_t n) {
    double* buf = NULL;
    Expr* packed = ndbuild_open_f64((int64_t)n, &buf);
    if (packed) { memcpy(buf, v, n * sizeof(double)); return packed; }

    Expr** items = malloc((n ? n : 1) * sizeof(Expr*));
    if (!items) return NULL;
    for (size_t i = 0; i < n; i++) items[i] = expr_new_real(v[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return r;
}

/* Growable double vector, for the two loops whose length is not known before
 * they run (NestWhileList, FixedPointList). */
typedef struct { double* v; size_t n, cap; } DVec;

static bool dvec_push(DVec* d, double x) {
    if (d->n == d->cap) {
        size_t cap = d->cap ? d->cap * 2 : 64;
        double* nv = realloc(d->v, cap * sizeof(double));
        if (!nv) return false;
        d->v = nv; d->cap = cap;
    }
    d->v[d->n++] = x;
    return true;
}

static void dvec_free(DVec* d) { free(d->v); d->v = NULL; d->n = d->cap = 0; }

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
    return expr_new_ndarray_like(x0, x0->data.ndarray.rank, x0->data.ndarray.dims,
                            buf, NDT_FLOAT64);   /* takes ownership of buf */
}

static Expr* numloop_nest_impl(const Expr* f, const Expr* x0, int64_t n, bool as_list) {
    if (numloop_off()) return NULL;
    if (n < 0) return NULL;

    /* The fused in-place array map answers with the FINAL array; NestList over
     * an array seed would have to keep n+1 whole arrays, which is the
     * interpreter's own allocation problem and none of the loop overhead this
     * path exists to remove. */
    if (is_f64_ndarray(x0)) return as_list ? NULL : numloop_nest_array(f, x0, n);

    double x;
    if (!to_machine_double(x0, &x)) return NULL;

    /* The iterate buffer is n+1 doubles; refuse a count whose byte size would
     * not compute, rather than wrapping into a short allocation. */
    if (as_list && (uint64_t)n + 1 > (uint64_t)(SIZE_MAX / sizeof(double))) return NULL;

    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;

    /* Result must be inexact: seed already Real/MPFR, or the body carries a
     * Real literal. Otherwise the interpreter would stay exact/symbolic. */
    if (!value_is_inexact(x0) && !body_inexact) { prog_free(&p); return NULL; }

    /* ...but a Real literal in the body only licenses positions the body
     * COMPUTED. NestList hands back the seed unevaluated at out[[1]], and so
     * does Nest at n = 0, so those need the seed itself to be Real:
     * NestList[# + 0. &, 1, 3] is {1, 1., 1., 1.}, not {1., 1., 1., 1.}. */
    if ((as_list || n == 0) && !value_is_inexact(x0)) { prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    double* out   = as_list ? malloc(((size_t)n + 1) * sizeof(double)) : NULL;
    if (!stack || (as_list && !out)) { free(stack); free(out); prog_free(&p); return NULL; }

    bool bail = false;
    if (as_list) out[0] = x;
    for (int64_t i = 0; i < n; i++) {
        x = numprog_run(&p, &x, stack);
        if (!isfinite(x)) { bail = true; break; }
        if (as_list) out[i + 1] = x;
    }
    free(stack);
    prog_free(&p);
    if (bail) { free(out); return NULL; }   /* interpreter re-runs the whole Nest */
    Expr* r = as_list ? reals_to_list(out, (size_t)n + 1) : expr_new_real(x);
    free(out);
    return r;
}

Expr* numloop_nest(const Expr* f, const Expr* x0, int64_t n) {
    return numloop_nest_impl(f, x0, n, false);
}

Expr* numloop_nestlist(const Expr* f, const Expr* x0, int64_t n) {
    return numloop_nest_impl(f, x0, n, true);
}

/* ------------------------------------------------------------------------
 *  Map[f, expr] at level {1}
 * ---------------------------------------------------------------------- */
static Expr* numloop_map_impl(const Expr* f, const Expr* expr, bool as_scan) {
    if (numloop_off()) return NULL;

    /* A PACKED List is accepted alongside a plain one, and that is not a
     * convenience: without it Map over a packed list fell past this fast path to
     * the ndarray leading-axis walk, which applies f through the INTERPRETER per
     * element. Measured at 10^6 elements, Map[#^2 &, x] was 424 ms packed against
     * 222 ms plain and Map[Sin[#] Exp[-#] &, x] 1120 ms against 180 ms -- so
     * automatic packing made the most-used functional head up to 6x SLOWER. Only
     * rank 1: Map over a matrix maps over its ROWS, and this body takes a scalar.
     */
    size_t n;
    bool   packed_src = false;
    NDType sdt = NDT_FLOAT64;
    if (is_packed_list(expr) && expr->data.ndarray.rank == 1) {
        sdt = expr->data.ndarray.dtype;
        if (sdt != NDT_FLOAT64 && sdt != NDT_INT64) return NULL;   /* not ours */
        n = (size_t)expr->data.ndarray.dims[0];
        packed_src = true;
    } else if (expr->type == EXPR_FUNCTION) {
        n = expr->data.function.arg_count;
    } else {
        return NULL;
    }
    if (n == 0) return NULL;   /* trivial; let the interpreter handle it */

    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;

    /* Every element must be a machine number; note whether they are ALL inexact. */
    double* vals = malloc(n * sizeof(double));
    if (!vals) { prog_free(&p); return NULL; }
    bool all_inexact = true, ok = true;
    if (packed_src) {
        /* THE O(1) TYPE DECISION. A buffer's dtype answers for every element at
         * once -- float64 means all of them are inexact machine reals, int64 that
         * all are exact Integers -- so the per-element probe the plain path has to
         * run disappears entirely. int64 widens through `double` exactly as
         * to_machine_double does for an Integer, rounding past 2^53 the same way,
         * because matching the plain List is the contract here. */
        all_inexact = (sdt == NDT_FLOAT64);
        if (sdt == NDT_FLOAT64) {
            memcpy(vals, expr->data.ndarray.data, n * sizeof(double));
        } else {
            const int64_t* src = (const int64_t*)expr->data.ndarray.data;
            for (size_t i = 0; i < n; i++) vals[i] = (double)src[i];
        }
    } else {
    for (size_t i = 0; i < n; i++) {
        if (!to_machine_double(expr->data.function.args[i], &vals[i])) { ok = false; break; }
        if (!value_is_inexact(expr->data.function.args[i])) all_inexact = false;
    }
    }
    /* Every result element must be inexact. Element i's is, if the body carries
     * a Real literal (so every element computes to a Real) or element i is
     * itself Real. ALL, not ANY: a body that just returns its argument passes
     * the exact ones straight through, so Map[# &, {1., 2, 3}] is {1., 2, 3} --
     * one Real element does not make the whole result inexact. */
    if (!ok || (!all_inexact && !body_inexact)) { free(vals); prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { free(vals); prog_free(&p); return NULL; }

    /* Every result element is a machine real (checked just above), so a List
     * result of any size worth packing is written straight into a buffer. Only
     * for a List head: Map[f, g[1., 2.]] keeps the head g, which no buffer can
     * carry. A bail discards the whole result either way -- the interpreter
     * re-runs the Map -- so there is nothing to abandon partially. */
    double* pbuf = NULL;
    Expr* packed = NULL;
    if (!as_scan && packed_src) {
        /* DERIVED from a buffer, so no threshold and no switch: Sin[packed] is
         * packed at any length and Map must be too. See ndbuild_open_like. */
        int64_t d1 = (int64_t)n;
        void* vb = NULL;
        packed = ndbuild_open_like(expr, 1, &d1, NDT_FLOAT64, &vb);
        pbuf = (double*)vb;
    } else if (!as_scan &&
               expr->data.function.head->type == EXPR_SYMBOL &&
               expr->data.function.head->data.symbol.name == SYM_List) {
        /* PRODUCED from a plain List: the producer rule, threshold and all. */
        packed = ndbuild_open_f64((int64_t)n, &pbuf);
    }

    /* Scan discards every f[element]; only Map builds a result. */
    Expr** out = (as_scan || packed) ? NULL : malloc(n * sizeof(Expr*));
    if (!as_scan && !packed && !out) { free(stack); free(vals); prog_free(&p); return NULL; }

    bool bail = false;
    size_t done = 0;
    for (size_t i = 0; i < n; i++) {
        double r = numprog_run(&p, &vals[i], stack);
        if (!isfinite(r)) { bail = true; break; }
        if (packed) pbuf[i] = r;
        else if (!as_scan) { out[i] = expr_new_real(r); done++; }
    }
    free(stack);
    free(vals);
    prog_free(&p);
    if (bail) {
        for (size_t i = 0; i < done; i++) expr_free(out[i]);
        free(out);
        if (packed) expr_free(packed);
        return NULL;
    }
    if (as_scan) return expr_new_symbol(SYM_Null);
    if (packed) return packed;
    /* A packed source always has the List head, and always takes the buffer path
     * above unless ndbuild_open declined (packing off, or under the threshold) --
     * in which case the result is a plain List of that head. */
    Expr* head = packed_src ? expr_new_symbol(SYM_List)
                            : expr_copy(expr->data.function.head);
    Expr* result = expr_new_function(head, out, n);
    free(out);
    return result;
}

Expr* numloop_map(const Expr* f, const Expr* expr) {
    return numloop_map_impl(f, expr, false);
}

/* Scan[f, expr] answers Null and exists for f's side effects -- of which a
 * numeric-closed body has none, so the loop below is observably a no-op and
 * could in principle be skipped outright. It is run anyway, because the ONE
 * thing the interpreter would still do is stop at a non-finite element (where
 * it goes complex, or divides by zero); running the compiled body keeps that
 * boundary exactly where the interpreter puts it, and hands those cases back. */
Expr* numloop_scan(const Expr* f, const Expr* expr) {
    return numloop_map_impl(f, expr, true);
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
        Expr* nd = expr_new_ndarray_raw(b->arr_rank, b->arr_dims, b->abuf[i], NDT_FLOAT64);
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
    Expr* nd = expr_new_ndarray_raw(pl->rank, pl->dims, pl->buf, NDT_FLOAT64);
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
static Expr* numloop_fold_impl(const Expr* f, const Expr* x0, const Expr* list,
                               bool as_list) {
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
     * whether they are ALL inexact. */
    double* elems = malloc(m * sizeof(double));
    if (!elems) { prog_free(&p); return NULL; }
    bool elems_inexact = true;
    bool ok = true;
    for (size_t i = 0; i < m; i++) {
        const Expr* el = list->data.function.args[i];
        if (!to_machine_double(el, &elems[i])) { ok = false; break; }
        if (!value_is_inexact(el)) elems_inexact = false;
    }
    if (!ok) { free(elems); prog_free(&p); return NULL; }

    /* Elements must be inexact throughout, not merely somewhere: a body that
     * returns its second argument (#2 &) hands each element back untouched, so
     * Fold[#2 &, 1., {1, 2, 3}] is the exact Integer 3.
     *
     * Given that, the accumulator is inexact if the seed is, or if the body
     * carries a Real literal -- except for FoldList, which emits the seed itself
     * as out[[1]] and so needs it Real however inexact the body is. */
    bool seed_ok = value_is_inexact(x0) || (!as_list && body_inexact);
    if (!elems_inexact || !seed_ok) { free(elems); prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { free(elems); prog_free(&p); return NULL; }
    /* FoldList keeps the seed and every partial result: m+1 values. */
    double* out = as_list ? malloc((m + 1) * sizeof(double)) : NULL;
    if (as_list && !out) { free(stack); free(elems); prog_free(&p); return NULL; }

    bool bail = false;
    if (as_list) out[0] = acc;
    for (size_t i = 0; i < m; i++) {
        double regs[2] = { acc, elems[i] };
        acc = numprog_run(&p, regs, stack);
        if (!isfinite(acc)) { bail = true; break; }
        if (as_list) out[i + 1] = acc;
    }
    free(stack);
    free(elems);
    prog_free(&p);
    if (bail) { free(out); return NULL; }
    Expr* r = as_list ? reals_to_list(out, m + 1) : expr_new_real(acc);
    free(out);
    return r;
}

Expr* numloop_fold(const Expr* f, const Expr* x0, const Expr* list) {
    return numloop_fold_impl(f, x0, list, false);
}

Expr* numloop_foldlist(const Expr* f, const Expr* x0, const Expr* list) {
    return numloop_fold_impl(f, x0, list, true);
}

/* ------------------------------------------------------------------------
 *  FixedPoint[f, x0]  (default SameTest, no application cap)
 * ---------------------------------------------------------------------- */
static Expr* numloop_fixedpoint_impl(const Expr* f, const Expr* x0, bool as_list) {
    if (numloop_off()) return NULL;

    double x;
    if (!to_machine_double(x0, &x)) return NULL;

    NumProg p;
    bool body_inexact;
    if (!compile_function(&p, f, 1, &body_inexact)) return NULL;
    if (!value_is_inexact(x0) && !body_inexact) { prog_free(&p); return NULL; }

    /* FixedPointList emits the seed unevaluated at out[[1]], so it needs a Real
     * seed even when the body would make every later iterate Real. The exact
     * seed also changes the LENGTH -- the interpreter's SameQ separates 1 from
     * 1., so FixedPointList[# + 0. &, 1] takes one extra step. */
    if (as_list && !value_is_inexact(x0)) { prog_free(&p); return NULL; }

    double* stack = malloc(p.max_stack * sizeof(double));
    if (!stack) { prog_free(&p); return NULL; }

    /* Iterate x := f(x) until it stops changing (SameQ on machine reals is
     * exact double equality, matching the interpreter's expr_eq). A non-
     * converging orbit hits the safety cap and bails to the interpreter.
     *
     * FixedPointList keeps the seed and every iterate INCLUDING the repeated
     * final value, so the list ends ..., fp, fp -- that trailing duplicate is
     * what the interpreted path produces (FixedPointList[Cos, 1.0] has length
     * 92 with [[-1]] === [[-2]]) and the fast path has to reproduce it. */
    bool bail = false;
    int64_t guard = 0;
    const int64_t CAP = 1000000;   /* == ITER_SAFETY_CAP */
    double cur = x;
    DVec acc = { NULL, 0, 0 };
    if (as_list && !dvec_push(&acc, cur)) bail = true;
    while (!bail) {
        double next = numprog_run(&p, &cur, stack);
        if (!isfinite(next)) { bail = true; break; }
        if (as_list && !dvec_push(&acc, next)) { bail = true; break; }
        bool same = (next == cur);
        cur = next;
        if (same) break;
        if (++guard >= CAP) { bail = true; break; }
    }
    free(stack);
    prog_free(&p);
    if (bail) { dvec_free(&acc); return NULL; }
    if (!as_list) return expr_new_real(cur);
    Expr* r = reals_to_list(acc.v, acc.n);
    dvec_free(&acc);
    return r;
}

Expr* numloop_fixedpoint(const Expr* f, const Expr* x0) {
    return numloop_fixedpoint_impl(f, x0, false);
}

Expr* numloop_fixedpointlist(const Expr* f, const Expr* x0) {
    return numloop_fixedpoint_impl(f, x0, true);
}

/* ------------------------------------------------------------------------
 *  NestWhile[f, x0, test]  (m = 1, default max / n; test a unary predicate)
 * ---------------------------------------------------------------------- */
static Expr* numloop_nestwhile_impl(const Expr* f, const Expr* x0, const Expr* test,
                                    bool as_list) {
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

    /* NestWhileList emits the seed unevaluated at out[[1]]. The scalar form
     * hands it back whole when the test fails on the first look, which is not
     * known until the loop runs -- so that case is caught after the loop. */
    if (as_list && !value_is_inexact(x0)) { prog_free(&pf); return NULL; }

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
     * that fails the test and return it. NestWhileList collects the seed and
     * every iterate up to and including that first failing value. */
    bool bail = false;
    int64_t guard = 0;
    const int64_t CAP = 1000000;
    double cur = x;
    bool applied = false;
    DVec acc = { NULL, 0, 0 };
    if (as_list && !dvec_push(&acc, cur)) bail = true;
    while (!bail && cmp_eval(numprog_run(&tl, &cur, stack),
                             numprog_run(&tr, &cur, stack), op)) {
        cur = numprog_run(&pf, &cur, stack);
        applied = true;
        if (!isfinite(cur)) { bail = true; break; }
        if (as_list && !dvec_push(&acc, cur)) { bail = true; break; }
        if (++guard >= CAP) { bail = true; break; }
    }
    free(stack);
    prog_free(&pf); prog_free(&tl); prog_free(&tr);
    /* The test failed immediately: the answer is the seed itself, unevaluated,
     * so an exact seed must stay exact -- NestWhile[# + 0. &, 1, # < 0 &] is 1,
     * not 1.. Nothing was computed, so there is no work to lose by handing it
     * back to the interpreter. */
    if (!applied && !value_is_inexact(x0)) bail = true;
    if (bail) { dvec_free(&acc); return NULL; }
    if (!as_list) return expr_new_real(cur);
    Expr* r = reals_to_list(acc.v, acc.n);
    dvec_free(&acc);
    return r;
}

Expr* numloop_nestwhile(const Expr* f, const Expr* x0, const Expr* test) {
    return numloop_nestwhile_impl(f, x0, test, false);
}

Expr* numloop_nestwhilelist(const Expr* f, const Expr* x0, const Expr* test) {
    return numloop_nestwhile_impl(f, x0, test, true);
}

/* ------------------------------------------------------------------------
 *  Accumulate[list]
 *
 *  The one head here with no function argument to compile: its interpreted
 *  path builds a Plus[running, element] node and calls evaluate() on it once
 *  per element, which is the entire cost. The running sum below is that same
 *  binary double addition in the same order, so the two agree bit-for-bit.
 *
 *  The gate is every element a machine Real -- not merely convertible. A list
 *  that mixes exact and inexact keeps an EXACT prefix (Accumulate[{1, 2., 3}]
 *  is {1, 3., 6.}, whose first element is an Integer), and a uniform buffer of
 *  doubles cannot express that; requiring Real throughout makes the case
 *  unreachable instead of approximating it. Rational and BigInt elements, which
 *  to_machine_double would happily narrow, are excluded for the same reason.
 * ---------------------------------------------------------------------- */
Expr* numloop_accumulate(const Expr* list) {
    if (numloop_off()) return NULL;
    if (!list || list->type != EXPR_FUNCTION) return NULL;
    size_t n = list->data.function.arg_count;
    if (n == 0) return NULL;   /* trivial; let the interpreter copy it */

    for (size_t i = 0; i < n; i++)
        if (list->data.function.args[i]->type != EXPR_REAL) return NULL;

    double* out = malloc(n * sizeof(double));
    if (!out) return NULL;

    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum = i ? sum + list->data.function.args[i]->data.real
                : list->data.function.args[i]->data.real;
        if (!isfinite(sum)) { free(out); return NULL; }   /* interpreter re-runs */
        out[i] = sum;
    }

    /* Accumulate preserves the argument's head, so this is not always a List. */
    Expr** items = malloc(n * sizeof(Expr*));
    if (!items) { free(out); return NULL; }
    for (size_t i = 0; i < n; i++) items[i] = expr_new_real(out[i]);
    free(out);
    Expr* r = expr_new_function(expr_copy(list->data.function.head), items, n);
    free(items);
    return r;
}
