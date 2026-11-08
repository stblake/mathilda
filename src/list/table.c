#include "list_common.h"
#include "table.h"
#include "compile/autocompile.h"
#include "iter.h"
#include "../pack.h"
#include "numeric.h"       /* numeric_mpq_to_double */
#include "checked_int.h"   /* ci_*_i64 */
#include <gmp.h>
#include <complex.h>   /* csin & co.: the builtins' machine formula */
#include <float.h>     /* DBL_MIN */


/* Element backstop for a finite Table range, matching Sum's SUM_MAX_FINITE_TERMS
 * and Product's PRODUCT_MAX_FINITE_TERMS. Reaching it returns Table[...]
 * unevaluated rather than silently truncating (the old 1000000 cap did the
 * latter, and at a hundredth of this size — so Table[i,{i,1,2000000}] came back
 * with 1000001 elements). With the exact termination test in iter_range_continue
 * a well-formed finite range terminates on its own; this only guards a range
 * whose element count genuinely exceeds the limit. */
#define TABLE_MAX_FINITE_ELEMENTS 100000000LL

/* ===========================================================================
 * Exact-rational body kernel for an integer iterator.
 *
 * Table[N[1/((i + j - 1) (i + j - 2)/2 + j)], {i, n}, {j, n}] spent ~2 us per
 * element in the evaluator -- rebuilding and re-evaluating the Plus / Times /
 * Power tree, making a Rational, then N -- against a few ns of arithmetic.
 * The machine auto-compiler cannot take it: the iterator is an exact Integer,
 * so the interpreter's answer is the EXACT value (correctly rounded by N), and
 * a double evaluation with intermediate roundings would differ in the last
 * bits.
 *
 * So this kernel evaluates the body exactly, in int64 rationals with
 * overflow-checked arithmetic. Exact rational arithmetic has one answer, so
 * every element it completes equals the interpreter's; an element it cannot
 * complete (an overflow, a zero base to a negative power, a non-integer
 * exponent, Factorial out of range) is handed to the interpreter on its own,
 * which then produces the authoritative value and any message. A top-level
 * N[g] / N[g, MachinePrecision] rounds the exact quotient once, the same
 * correctly-rounded conversion numericalize applies to a Rational leaf.
 *
 * Accepted body: integer and int64-Rational literals, the iterator, any other
 * symbol whose OwnValue is already an exact int64 number (an enclosing
 * iterator, a constant n), combined by Plus, Times, Power, Abs, Min, Max and
 * Factorial. Nothing else -- no user rule ever runs, so skipping the evaluator
 * is not observable.
 *
 * Under a top-level N, also Sin, Cos, Tan, Sinh and Cosh of an exact INTEGER
 * (Table[N[Cos[i j]/(i + j)], ...]).  The interpreter keeps Cos[6] symbolic,
 * folds every exact factor of a product into one Rational (Times is Flat) and
 * every exact term of a sum into one number, and N then numericalises the
 * leaves and multiplies / adds the machine numbers.  So an inexact value is
 * carried as a DEFERRED product coef * atom -- coef the exact Rational, atom
 * the machine value of the one transcendental factor -- and materialised
 * exactly where N would form it: as N[coef] * atom, or N[c] + N[coef] * atom
 * inside a sum.  A product or sum with two inexact operands, a transcendental
 * of a non-integer, a power or Abs/Min/Max of an inexact value, or a
 * non-normal result hands that ELEMENT to the interpreter: there the exact
 * evaluator's term merging or a different numericalisation could change the
 * rounding.  The machine value of f[k] is the builtin's own formula,
 * creal(f(k + 0 I)) (csin, ccos, ctan, csinh, ccosh), and f[0] is the exact
 * 0 or 1 the builtin returns. Compiled once per Table call (a nested Table recompiles
 * its inner body per outer value, which is what folds the outer iterator in
 * as a constant).
 * ======================================================================== */
typedef struct { int64_t n, d; } TxRat;          /* d > 0, gcd(|n|, d) = 1 */

typedef enum {
    TX_CONST, TX_VAR, TX_PLUS, TX_TIMES, TX_POWER, TX_ABS, TX_MIN, TX_MAX, TX_FACT,
    TX_FN                       /* machine mode only: a = TxFn */
} TxOp;

typedef enum { TXF_SIN, TXF_COS, TXF_TAN, TXF_SINH, TXF_COSH } TxFn;

typedef struct { uint8_t op; int32_t a; } TxInsn;

typedef struct {
    TxInsn* code;   size_t ncode, ccode;
    TxRat*  consts; size_t nconst, cconst;
    long    depth, max_depth;
    bool    ok;
    bool    machine;            /* body was N[g]: answer the rounded double */
    bool    has_fn;             /* contains TX_FN: run by tx_run_machine */
    const char* var;            /* interned iterator name */
} TxProg;

/* Binary gcd of a, b >= 0: shifts and subtractions only. The Euclidean
 * loop's 64-bit divisions dominated the kernel's cost (~40 cycles each). */
static int tx_ctz(uint64_t v) {                  /* v != 0 */
#if defined(__GNUC__)
    return __builtin_ctzll(v);
#else
    int k = 0;
    while (!(v & 1)) { v >>= 1; k++; }
    return k;
#endif
}
static int64_t tx_gcd(int64_t a0, int64_t b0) {
    uint64_t a = (uint64_t)a0, b = (uint64_t)b0;
    if (a == 0) return (int64_t)b;
    if (b == 0 || a == 1 || b == 1) return a == 0 ? (int64_t)b : (b == 0 ? (int64_t)a : 1);
    int sh = tx_ctz(a | b);
    a >>= tx_ctz(a);
    do {
        b >>= tx_ctz(b);
        if (a > b) { uint64_t t = a; a = b; b = t; }
        if (a == 1) break;                       /* coprime odd parts */
        b -= a;
    } while (b);
    return (int64_t)(a << sh);
}

static void tx_emit(TxProg* p, TxOp op, int32_t a, long delta) {
    if (!p->ok) return;
    if (p->ncode == p->ccode) {
        size_t nc = p->ccode ? 2 * p->ccode : 32;
        TxInsn* g = realloc(p->code, nc * sizeof(TxInsn));
        if (!g) { p->ok = false; return; }
        p->code = g; p->ccode = nc;
    }
    p->code[p->ncode].op = (uint8_t)op;
    p->code[p->ncode].a = a;
    p->ncode++;
    p->depth += delta;
    if (p->depth > p->max_depth) p->max_depth = p->depth;
}

static void tx_const(TxProg* p, int64_t n, int64_t d) {
    if (!p->ok) return;
    if (n == INT64_MIN || d <= 0) { p->ok = false; return; }
    if (p->nconst == p->cconst) {
        size_t nc = p->cconst ? 2 * p->cconst : 8;
        TxRat* g = realloc(p->consts, nc * sizeof(TxRat));
        if (!g) { p->ok = false; return; }
        p->consts = g; p->cconst = nc;
    }
    p->consts[p->nconst].n = n;
    p->consts[p->nconst].d = d;
    tx_emit(p, TX_CONST, (int32_t)p->nconst++, +1);
}

/* Parameter bindings while compiling an inlined pure function: a parameter
 * (named, or Slot[k] when `slots`) compiles as its argument expression, which
 * is compiled in the CALLER's environment -- exactly beta-reduction. */
typedef struct TxEnv {
    const char* const* names;    /* NULL for a slot function */
    const Expr* const* args;
    size_t n;
    const struct TxEnv* parent;
} TxEnv;

#define TX_MAX_INLINE_DEPTH 16

static void tx_walk(TxProg* p, const Expr* e, const TxEnv* env, int depth);

/* The pure Function a call head names, when finding it is a pure lookup:
 *   - a symbol whose OwnValue is a Function (f = Function[...]; f[i, j]);
 *   - sym[lits...] matching a literal DownValue of sym exactly, when EVERY
 *     DownValue of sym is literal and sym has no OwnValue (ent[1] =
 *     Function[...]; ent[1][i, j]). Literal rules are tried first and an
 *     exact match is unique, so this is the rule the interpreter applies.
 * NULL otherwise. Nothing is evaluated. */
static const Expr* tx_lookup_function(const Expr* h) {
    const Expr* v = NULL;
    if (h->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_lookup(h->data.symbol.name);
        if (def && def->own_values && !def->own_values->next)
            v = def->own_values->replacement;
    } else if (h->type == EXPR_FUNCTION && h->data.function.head->type == EXPR_SYMBOL) {
        SymbolDef* def = symtab_lookup(h->data.function.head->data.symbol.name);
        if (!def || def->own_values || !def->down_values) return NULL;
        for (const Rule* r = def->down_values; r; r = r->next) {
            if (!r->literal || r->scoped_cond) return NULL;
            if (!v && expr_eq(r->pattern, (Expr*)h)) v = r->replacement;
        }
    }
    if (!v || v->type != EXPR_FUNCTION || v->data.function.head->type != EXPR_SYMBOL ||
        v->data.function.head->data.symbol.name != SYM_Function) return NULL;
    size_t fa = v->data.function.arg_count;
    if (fa != 1 && fa != 2) return NULL;          /* no attribute argument */
    return v;
}

/* Compile f[args...] for a pure function f by inlining its body. */
static void tx_inline_call(TxProg* p, const Expr* fn, const Expr* call,
                           const TxEnv* env, int depth) {
    if (depth >= TX_MAX_INLINE_DEPTH) { p->ok = false; return; }
    size_t argc = call->data.function.arg_count;
    const char* names[16];
    TxEnv inner = { NULL, (const Expr* const*)call->data.function.args, argc, env };
    const Expr* body;
    if (fn->data.function.arg_count == 1) {
        body = fn->data.function.args[0];         /* Function[body]: slots */
    } else {
        const Expr* ps = fn->data.function.args[0];
        body = fn->data.function.args[1];
        if (ps->type == EXPR_SYMBOL) {
            if (argc != 1) { p->ok = false; return; }
            names[0] = ps->data.symbol.name;
        } else if (ps->type == EXPR_FUNCTION && ps->data.function.head->type == EXPR_SYMBOL &&
                   ps->data.function.head->data.symbol.name == SYM_List &&
                   ps->data.function.arg_count == argc && argc <= 16) {
            for (size_t i = 0; i < argc; i++) {
                if (ps->data.function.args[i]->type != EXPR_SYMBOL) { p->ok = false; return; }
                names[i] = ps->data.function.args[i]->data.symbol.name;
            }
        } else {
            p->ok = false;
            return;
        }
        inner.names = names;
    }
    tx_walk(p, body, &inner, depth + 1);
}

static void tx_walk(TxProg* p, const Expr* e, const TxEnv* env, int depth) {
    if (!p->ok) return;
    int64_t rn, rd;
    if (e->type == EXPR_INTEGER) { tx_const(p, e->data.integer, 1); return; }
    if (e->type == EXPR_SYMBOL) {
        if (env && env->names)
            for (size_t i = 0; i < env->n; i++)
                if (env->names[i] == e->data.symbol.name) {
                    tx_walk(p, env->args[i], env->parent, depth);
                    return;
                }
        if (e->data.symbol.name == p->var) { tx_emit(p, TX_VAR, 0, +1); return; }
        /* A symbol already holding an exact int64 number is a constant for the
         * whole run: the accepted body cannot assign anything. */
        SymbolDef* def = symtab_lookup(e->data.symbol.name);
        const Expr* v = (def && def->own_values) ? def->own_values->replacement : NULL;
        if (v && v->type == EXPR_INTEGER) { tx_const(p, v->data.integer, 1); return; }
        if (v && is_rational(v, &rn, &rd)) { tx_const(p, rn, rd); return; }
        p->ok = false;
        return;
    }
    if (e->type != EXPR_FUNCTION) { p->ok = false; return; }
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head->type == EXPR_SYMBOL && head->data.symbol.name == SYM_Slot) {
        /* Slot[k] of the innermost slot function. */
        if (!env || env->names || argc != 1 || e->data.function.args[0]->type != EXPR_INTEGER) {
            p->ok = false;
            return;
        }
        int64_t k = e->data.function.args[0]->data.integer;
        if (k < 1 || (uint64_t)k > env->n) { p->ok = false; return; }
        tx_walk(p, env->args[k - 1], env->parent, depth);
        return;
    }
    if (is_rational(e, &rn, &rd)) { tx_const(p, rn, rd); return; }
    TxOp op = TX_CONST;
    bool known = false;
    if (head->type == EXPR_SYMBOL && argc == 1 && p->machine) {
        const char* h = head->data.symbol.name;
        int fn = h == SYM_Sin ? TXF_SIN : h == SYM_Cos ? TXF_COS : h == SYM_Tan ? TXF_TAN :
                 h == SYM_Sinh ? TXF_SINH : h == SYM_Cosh ? TXF_COSH : -1;
        if (fn >= 0) {
            tx_walk(p, e->data.function.args[0], env, depth);
            tx_emit(p, TX_FN, fn, 0);
            p->has_fn = true;
            return;
        }
    }
    if (head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol.name;
        known = true;
        if      (h == SYM_Plus  && argc >= 1) op = TX_PLUS;
        else if (h == SYM_Times && argc >= 1) op = TX_TIMES;
        else if (h == SYM_Power && argc == 2) op = TX_POWER;
        else if (h == SYM_Abs   && argc == 1) op = TX_ABS;
        else if (h == SYM_Min   && argc >= 1) op = TX_MIN;
        else if (h == SYM_Max   && argc >= 1) op = TX_MAX;
        else if (h == SYM_Factorial && argc == 1) op = TX_FACT;
        else known = false;
    }
    if (!known) {
        /* A call of a pure function the head names by lookup: inline it. */
        const Expr* fn = tx_lookup_function(head);
        if (!fn) { p->ok = false; return; }
        tx_inline_call(p, fn, e, env, depth);
        return;
    }
    if (argc > 1000000) { p->ok = false; return; }
    for (size_t i = 0; i < argc; i++) tx_walk(p, e->data.function.args[i], env, depth);
    tx_emit(p, op, (int32_t)argc, 1 - (long)argc);
}

static void tx_free(TxProg* p) { free(p->code); free(p->consts); }

/* Compile `body` over iterator `var`. False (nothing to free) to decline. */
static bool tx_compile(TxProg* p, const Expr* body, const Expr* var) {
    memset(p, 0, sizeof(*p));
    p->ok = true;
    p->var = var->data.symbol.name;
    if (body->type == EXPR_FUNCTION && body->data.function.head->type == EXPR_SYMBOL &&
        body->data.function.head->data.symbol.name == SYM_N &&
        (body->data.function.arg_count == 1 ||
         (body->data.function.arg_count == 2 &&
          body->data.function.args[1]->type == EXPR_SYMBOL &&
          body->data.function.args[1]->data.symbol.name == SYM_MachinePrecision))) {
        p->machine = true;
        body = body->data.function.args[0];
    }
    tx_walk(p, body, NULL, 0);
    if (!p->ok || p->depth != 1) { tx_free(p); return false; }
    return true;
}

/* r = a + b; false on overflow. */
static bool tx_add(TxRat a, TxRat b, TxRat* r) {
    int64_t x, y, n, d;
    if (a.d == 1 && b.d == 1) {
        if (ci_add_i64(a.n, b.n, &n) || n == INT64_MIN) return false;
        r->n = n; r->d = 1;
        return true;
    }
    if (a.d == 1 || b.d == 1) {
        /* integer + p/q = (integer q + p)/q, already reduced: gcd(kq + p, q)
         * = gcd(p, q) = 1. No gcd loop. */
        if (b.d == 1) { TxRat t = a; a = b; b = t; }      /* a is the integer */
        if (ci_mul_i64(a.n, b.d, &x) || ci_add_i64(x, b.n, &n) || n == INT64_MIN)
            return false;
        r->n = n; r->d = b.d;
        return true;
    }
    if (ci_mul_i64(a.n, b.d, &x) || ci_mul_i64(b.n, a.d, &y) ||
        ci_add_i64(x, y, &n) || ci_mul_i64(a.d, b.d, &d) || n == INT64_MIN) return false;
    int64_t g = tx_gcd(n < 0 ? -n : n, d);
    if (g > 1) { n /= g; d /= g; }
    r->n = n; r->d = d;
    return true;
}

/* r = a * b, cross-reduced first so the products stay small. */
static bool tx_mul(TxRat a, TxRat b, TxRat* r) {
    if (a.d == 1 && b.d == 1) {
        int64_t n;
        if (ci_mul_i64(a.n, b.n, &n) || n == INT64_MIN) return false;
        r->n = n; r->d = 1;
        return true;
    }
    int64_t g1 = tx_gcd(a.n < 0 ? -a.n : a.n, b.d);   /* >= 1: b.d >= 1 */
    int64_t g2 = tx_gcd(b.n < 0 ? -b.n : b.n, a.d);
    int64_t n, d;
    if (g1 > 1) { a.n /= g1; b.d /= g1; }
    if (g2 > 1) { b.n /= g2; a.d /= g2; }
    if (ci_mul_i64(a.n, b.n, &n) || ci_mul_i64(a.d, b.d, &d) ||
        n == INT64_MIN) return false;
    if (n == 0) d = 1;
    r->n = n; r->d = d;
    return true;
}

/* *c = sign(a - b); false on overflow. */
static bool tx_cmp(TxRat a, TxRat b, int* c) {
    int64_t x, y;
    if (ci_mul_i64(a.n, b.d, &x) || ci_mul_i64(b.n, a.d, &y)) return false;
    *c = (x < y) ? -1 : (x > y);
    return true;
}

static const int64_t TX_FACTORIALS[21] = {
    INT64_C(1), INT64_C(1), INT64_C(2), INT64_C(6), INT64_C(24), INT64_C(120),
    INT64_C(720), INT64_C(5040), INT64_C(40320), INT64_C(362880), INT64_C(3628800),
    INT64_C(39916800), INT64_C(479001600), INT64_C(6227020800), INT64_C(87178291200),
    INT64_C(1307674368000), INT64_C(20922789888000), INT64_C(355687428096000),
    INT64_C(6402373705728000), INT64_C(121645100408832000), INT64_C(2432902008176640000)
};

/* Evaluate at iterator value `x`. False: this element must be interpreted. */
static bool tx_run(const TxProg* p, int64_t x, TxRat* stack, TxRat* out) {
    size_t sp = 0;
    for (size_t i = 0; i < p->ncode; i++) {
        const TxInsn* in = &p->code[i];
        switch ((TxOp)in->op) {
            case TX_CONST: stack[sp++] = p->consts[in->a]; break;
            case TX_VAR:   /* INT64_MIN has no negation: let the interpreter */
                           if (x == INT64_MIN) return false;
                           stack[sp].n = x; stack[sp].d = 1; sp++; break;
            case TX_PLUS: case TX_TIMES: {
                size_t k = (size_t)in->a, base = sp - k;
                TxRat acc = stack[base];
                for (size_t j = 1; j < k; j++) {
                    bool ok = (in->op == TX_PLUS) ? tx_add(acc, stack[base + j], &acc)
                                                  : tx_mul(acc, stack[base + j], &acc);
                    if (!ok) return false;
                }
                sp = base;
                stack[sp++] = acc;
                break;
            }
            case TX_MIN: case TX_MAX: {
                size_t k = (size_t)in->a, base = sp - k;
                TxRat acc = stack[base];
                for (size_t j = 1; j < k; j++) {
                    int c;
                    if (!tx_cmp(stack[base + j], acc, &c)) return false;
                    if ((in->op == TX_MIN) ? (c < 0) : (c > 0)) acc = stack[base + j];
                }
                sp = base;
                stack[sp++] = acc;
                break;
            }
            case TX_POWER: {
                TxRat b = stack[sp - 2], e = stack[sp - 1];
                sp -= 2;
                if (e.d != 1) return false;                 /* a radical */
                int64_t k = e.n, n = b.n, d = b.d;
                if (k == 0) {
                    if (n == 0) return false;               /* 0^0: Indeterminate */
                    n = 1; d = 1;
                } else {
                    if (k < 0) {
                        if (n == 0) return false;           /* 1/0: ComplexInfinity */
                        int64_t t = n; n = d; d = t;        /* invert, sign to n */
                        if (d < 0) { n = -n; d = -d; }
                        k = -k;
                    }
                    if (ci_powi_i64(n, k, &n) || ci_powi_i64(d, k, &d) ||
                        n == INT64_MIN || d == INT64_MIN) return false;
                }
                stack[sp].n = n; stack[sp].d = d; sp++;
                break;
            }
            case TX_ABS:
                if (stack[sp - 1].n < 0) stack[sp - 1].n = -stack[sp - 1].n;
                break;
            case TX_FACT: {
                TxRat a = stack[sp - 1];
                if (a.d != 1 || a.n < 0 || a.n > 20) return false;
                stack[sp - 1].n = TX_FACTORIALS[a.n];
                break;
            }
            case TX_FN: return false;               /* tx_run_machine's */
        }
    }
    *out = stack[0];
    return true;
}

/* int64 -> mpz without assuming long is 64 bits (LLP64). */
static void tx_mpz_set_i64(mpz_ptr z, int64_t v) {
    uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, 1, sizeof(uint64_t), 0, 0, &m);
    if (v < 0) mpz_neg(z, z);
}

static void tx_mpq_set(mpq_ptr q, TxRat r) {       /* r already canonical */
    tx_mpz_set_i64(mpq_numref(q), r.n);
    tx_mpz_set_i64(mpq_denref(q), r.d);
}

/* The double N gives a Rational/Integer leaf: one correctly-rounded division
 * when both parts are exact doubles, else the exact rounding numericalize uses. */
static double tx_to_double(TxRat r) {
    const int64_t lim = INT64_C(1) << 53;
    if (r.n <= lim && r.n >= -lim && r.d <= lim)
        return (double)r.n / (double)r.d;
    mpq_t q;
    mpq_init(q);
    tx_mpq_set(q, r);
    double v = numeric_mpq_to_double(q);
    mpq_clear(q);
    return v;
}


/* ---- machine tier: transcendental leaves under N -----------------------
 * A stack value is exact (r) or the deferred product r * atom (see the kernel
 * comment above).  Mirrors the interpreter's exact phase -- Flat products and
 * sums fold their exact parts -- and then N's numericalisation of the result. */
typedef struct { TxRat r; double atom; bool inexact; } TxM;

/* N of a deferred product: N[coef] * atom, or atom itself for coef 1. */
static double txm_value(TxM v) {
    if (!v.inexact) return tx_to_double(v.r);
    if (v.r.n == 1 && v.r.d == 1) return v.atom;
    return tx_to_double(v.r) * v.atom;
}

static bool txm_normal(double v) { return isfinite(v) && (v == 0.0 || fabs(v) >= DBL_MIN); }

static int tx_bitlen(uint64_t a) { int b = 0; for (; a; a >>= 1) b++; return b; }

/* Would this exact leaf, left in the expression N sees, make N raise the
 * evaluation to MPFR?  numeric_plan_working_spec does that when the expression
 * has a call and an exact leaf whose double is not its value with error worse
 * than the output can absorb: an integer beyond 53 bits that is no double, or
 * a non-dyadic rational with more numerator than denominator bits (N[-12/7 +
 * Sinh[12]/50] is evaluated at 128 bits and rounded, not leaf by leaf).  Such
 * an element is the interpreter's. */
static bool tx_leaf_raises(TxRat r) {
    if (r.n == 0) return false;
    uint64_t a = r.n < 0 ? (uint64_t)0 - (uint64_t)r.n : (uint64_t)r.n;
    uint64_t b = (uint64_t)r.d;
    int la = tx_bitlen(a);
    bool a_exact = la <= 53 || (a & ((((uint64_t)1) << (la - 53)) - 1u)) == 0;
    if (b == 1) return !a_exact;
    if ((b & (b - 1)) == 0 && a_exact) return false;      /* dyadic: its own double */
    return la - tx_bitlen(b) > 0;
}

/* Scratch for tx_run_machine's exact ops, each sized max_depth + 1. */
typedef struct { TxInsn* code; TxRat* consts; TxRat* stack; } TxScratch;

static bool tx_run_machine(const TxProg* p, int64_t x, TxM* st, const TxScratch* sc,
                           double* out) {
    size_t sp = 0;
    for (size_t i = 0; i < p->ncode; i++) {
        const TxInsn* in = &p->code[i];
        switch ((TxOp)in->op) {
            case TX_CONST: st[sp].r = p->consts[in->a]; st[sp].inexact = false; sp++; break;
            case TX_VAR:
                if (x == INT64_MIN) return false;
                st[sp].r.n = x; st[sp].r.d = 1; st[sp].inexact = false; sp++;
                break;
            case TX_PLUS: case TX_TIMES: {
                size_t k = (size_t)in->a, base = sp - k;
                bool plus = in->op == TX_PLUS;
                TxRat c = { plus ? 0 : 1, 1 };
                const TxM* one = NULL;               /* the single inexact operand */
                for (size_t j = 0; j < k; j++) {
                    const TxM* v = &st[base + j];
                    if (v->inexact) {
                        if (one) return false;       /* two inexact operands */
                        one = v;
                        continue;
                    }
                    if (!(plus ? tx_add(c, v->r, &c) : tx_mul(c, v->r, &c))) return false;
                }
                TxM res;
                if (!one) {
                    res.r = c; res.inexact = false;
                } else if (plus) {
                    if (c.n == 0) res = *one;        /* Plus[0, t] is t */
                    else {
                        if (tx_leaf_raises(c) || tx_leaf_raises(one->r)) return false;
                        res.atom = tx_to_double(c) + txm_value(*one);
                        if (!txm_normal(res.atom)) return false;
                        res.r.n = 1; res.r.d = 1; res.inexact = true;
                    }
                } else {
                    TxRat cr;
                    if (!tx_mul(c, one->r, &cr)) return false;
                    if (cr.n == 0) { res.r = cr; res.inexact = false; }   /* 0 * t */
                    else { res.r = cr; res.atom = one->atom; res.inexact = true; }
                }
                sp = base;
                st[sp++] = res;
                break;
            }
            case TX_FN: {
                TxM a = st[sp - 1];
                if (a.inexact || a.r.d != 1) return false;     /* integer argument only */
                const int64_t lim = INT64_C(1) << 53;
                if (a.r.n > lim || a.r.n < -lim) return false;
                if (a.r.n == 0) {                              /* the builtins' exact f[0] */
                    st[sp - 1].r.n = (in->a == TXF_COS || in->a == TXF_COSH) ? 1 : 0;
                    st[sp - 1].r.d = 1;
                    break;
                }
                double complex z = (double)a.r.n;
                double complex w;
                switch ((TxFn)in->a) {
                    case TXF_SIN:  w = csin(z);  break;
                    case TXF_COS:  w = ccos(z);  break;
                    case TXF_TAN:  w = ctan(z);  break;
                    case TXF_SINH: w = csinh(z); break;
                    default:       w = ccosh(z); break;
                }
                double v = creal(w);
                if (!txm_normal(v) || v == 0.0) return false;
                st[sp - 1].r.n = 1; st[sp - 1].r.d = 1;
                st[sp - 1].atom = v; st[sp - 1].inexact = true;
                break;
            }
            default: {
                /* Power, Abs, Min, Max, Factorial: exact operands only, run by
                 * tx_run itself as a program that loads the k operands as
                 * constants and applies this one instruction -- the exact
                 * kernel's own arithmetic, with no second copy of it. */
                size_t k = (in->op == TX_POWER) ? 2
                         : (in->op == TX_ABS || in->op == TX_FACT) ? 1 : (size_t)in->a;
                size_t base = sp - k;
                for (size_t j = 0; j < k; j++) {
                    if (st[base + j].inexact) return false;
                    sc->consts[j] = st[base + j].r;
                    sc->code[j].op = (uint8_t)TX_CONST;
                    sc->code[j].a = (int32_t)j;
                }
                sc->code[k] = *in;
                TxProg one = *p;
                one.code = sc->code;   one.ncode = k + 1;
                one.consts = sc->consts; one.nconst = k;
                TxRat r;
                if (!tx_run(&one, x, sc->stack, &r)) return false;
                sp = base;
                st[sp].r = r; st[sp].inexact = false; sp++;
                break;
            }
        }
    }
    if (st[0].inexact && tx_leaf_raises(st[0].r)) return false;
    double v = txm_value(st[0]);
    if (!txm_normal(v)) return false;
    *out = v;
    return true;
}

/* ---- GMP tier ---------------------------------------------------------
 * An element whose int64 run overflowed (2^-i 3^-j, 1/(i! j!) at i = 30) is
 * rerun here in mpq_t before falling back to the interpreter: still exact,
 * so still the interpreter's value, at a fraction of its cost. Sizes are
 * capped so a runaway power goes to the interpreter rather than being built
 * here. */
#define TX_MPQ_MAX_BITS 1000000

typedef struct { mpq_t* s; size_t n; bool init; } TxMpqStack;

static bool tx_mpq_too_big(mpq_srcptr q) {
    return mpz_sizeinbase(mpq_numref(q), 2) > TX_MPQ_MAX_BITS ||
           mpz_sizeinbase(mpq_denref(q), 2) > TX_MPQ_MAX_BITS;
}

static bool tx_run_mpq(const TxProg* p, int64_t x, TxMpqStack* st, mpq_ptr out) {
    if (!st->init) {
        st->s = malloc(st->n * sizeof(mpq_t));
        if (!st->s) return false;
        for (size_t i = 0; i < st->n; i++) mpq_init(st->s[i]);
        st->init = true;
    }
    mpq_t* S = st->s;
    size_t sp = 0;
    for (size_t i = 0; i < p->ncode; i++) {
        const TxInsn* in = &p->code[i];
        switch ((TxOp)in->op) {
            case TX_CONST: tx_mpq_set(S[sp++], p->consts[in->a]); break;
            case TX_VAR: {
                TxRat r = { x, 1 };
                tx_mpq_set(S[sp++], r);
                break;
            }
            case TX_PLUS: case TX_TIMES: {
                size_t k = (size_t)in->a, base = sp - k;
                for (size_t j = 1; j < k; j++) {
                    if (in->op == TX_PLUS) mpq_add(S[base], S[base], S[base + j]);
                    else                   mpq_mul(S[base], S[base], S[base + j]);
                    if (tx_mpq_too_big(S[base])) return false;
                }
                sp = base + 1;
                break;
            }
            case TX_MIN: case TX_MAX: {
                size_t k = (size_t)in->a, base = sp - k;
                for (size_t j = 1; j < k; j++) {
                    int c = mpq_cmp(S[base + j], S[base]);
                    if ((in->op == TX_MIN) ? (c < 0) : (c > 0)) mpq_set(S[base], S[base + j]);
                }
                sp = base + 1;
                break;
            }
            case TX_POWER: {
                mpq_ptr b = S[sp - 2], e = S[sp - 1];
                if (mpz_cmp_ui(mpq_denref(e), 1) != 0 ||
                    !mpz_fits_slong_p(mpq_numref(e))) return false;
                long k = mpz_get_si(mpq_numref(e));
                sp -= 1;                                   /* the result replaces b */
                if (k == 0) {
                    if (mpq_sgn(b) == 0) return false;     /* 0^0 */
                    mpq_set_ui(b, 1, 1);
                    break;
                }
                if (k < 0) {
                    if (mpq_sgn(b) == 0) return false;     /* 1/0 */
                    mpq_inv(b, b);
                    k = -k;
                }
                /* the result has ~k * bits(base) bits: refuse a runaway */
                size_t bn = mpz_sizeinbase(mpq_numref(b), 2);
                size_t bd = mpz_sizeinbase(mpq_denref(b), 2);
                size_t bb = bn > bd ? bn : bd;
                if (bb > 1 && (unsigned long)k > TX_MPQ_MAX_BITS / (bb - 1)) return false;
                /* powers of coprime parts stay coprime: still canonical */
                mpz_pow_ui(mpq_numref(b), mpq_numref(b), (unsigned long)k);
                mpz_pow_ui(mpq_denref(b), mpq_denref(b), (unsigned long)k);
                break;
            }
            case TX_ABS:
                mpq_abs(S[sp - 1], S[sp - 1]);
                break;
            case TX_FACT: {
                mpq_ptr a = S[sp - 1];
                if (mpz_cmp_ui(mpq_denref(a), 1) != 0 || mpz_sgn(mpq_numref(a)) < 0 ||
                    mpz_cmp_ui(mpq_numref(a), 20000) > 0) return false;
                mpz_fac_ui(mpq_numref(a), mpz_get_ui(mpq_numref(a)));
                break;
            }
            case TX_FN: return false;               /* tx_run_machine's */
        }
    }
    mpq_set(out, S[0]);
    return true;
}

static void tx_mpq_stack_free(TxMpqStack* st) {
    if (st->init) {
        for (size_t i = 0; i < st->n; i++) mpq_clear(st->s[i]);
        free(st->s);
    }
}

/* The element loop of table_exact_int for a program with transcendental
 * leaves (tx_run_machine), kept apart so the exact loop stays exactly as it
 * was: sharing one loop cost the plain exact kernel ~10%.  Takes ownership of
 * *p.  Every result is a machine number under N, so the list is built packed;
 * an element the machine tier declines is the interpreter's, as in the exact
 * loop. */
static Expr* table_exact_int_fn(TxProg* p, Expr* body, Expr* var, int64_t imin,
                                int64_t di, size_t n) {
    size_t depth = (size_t)p->max_depth + 1;
    TxM* mstack = malloc(depth * sizeof(TxM));
    TxScratch sc = { malloc((depth + 1) * sizeof(TxInsn)), malloc(depth * sizeof(TxRat)),
                     malloc(depth * sizeof(TxRat)) };
    double* dbuf = NULL;
    bool alloc_ok = mstack && sc.code && sc.consts && sc.stack;
    Expr* packed = alloc_ok ? ndbuild_open_f64((int64_t)n, &dbuf) : NULL;
    if (!packed) {
        free(mstack); free(sc.code); free(sc.consts); free(sc.stack); tx_free(p);
        return NULL;
    }
    Expr** results = NULL;
    for (size_t i = 0; i < n; i++) {
        int64_t x = (int64_t)((uint64_t)imin + (uint64_t)i * (uint64_t)di);
        double dv;
        if (tx_run_machine(p, x, mstack, &sc, &dv)) {
            if (packed) dbuf[i] = dv;
            else        results[i] = expr_new_real(dv);
            continue;
        }
        /* This element is the interpreter's: bind and evaluate it. */
        Expr* xv = expr_new_integer(x);
        iter_bind(var, xv);
        expr_free(xv);
        Expr* v = evaluate(body);
        if (packed) {
            if (v && v->type == EXPR_REAL) { dbuf[i] = v->data.real; expr_free(v); continue; }
            /* Not a machine real: finish as a boxed List. */
            results = malloc(sizeof(Expr*) * n);
            if (!results) {
                expr_free(v); expr_free(packed); free(mstack);
                free(sc.code); free(sc.consts); free(sc.stack); tx_free(p);
                return NULL;
            }
            ndbuild_abandon(packed, i, results);
            packed = NULL;
        }
        results[i] = v;
    }
    free(mstack);
    free(sc.code); free(sc.consts); free(sc.stack);
    tx_free(p);
    if (packed) return packed;
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, n);
    free(results);
    return pack_offer(list);
}

/* Table[body, {var, imin, imax, di}] over exact int64 imin/imax/di through the
 * kernel. Returns the finished list, or NULL (having evaluated nothing) when
 * the body is outside the kernel's subset. `var` must already be shadowed by
 * the caller; an element the kernel declines is bound and interpreted. */
static Expr* table_exact_int(Expr* body, Expr* var, int64_t imin, int64_t imax, int64_t di) {
    if (di == 0 || var->type != EXPR_SYMBOL) return NULL;
    uint64_t total;   /* element count, exactly, in unsigned arithmetic */
    if (di > 0) {
        if (imax < imin) return NULL;          /* empty: the general path answers {} */
        total = ((uint64_t)imax - (uint64_t)imin) / (uint64_t)di + 1;
    } else {
        if (imax > imin) return NULL;
        total = ((uint64_t)imin - (uint64_t)imax) / ((uint64_t)0 - (uint64_t)di) + 1;
    }
    if (total < 2 || total > (uint64_t)TABLE_MAX_FINITE_ELEMENTS) return NULL;

    TxProg p;
    if (!tx_compile(&p, body, var)) return NULL;
    if (p.has_fn) return table_exact_int_fn(&p, body, var, imin, di, (size_t)total);
    TxRat* stack = malloc((size_t)(p.max_depth + 1) * sizeof(TxRat));
    if (!stack) { tx_free(&p); return NULL; }

    size_t n = (size_t)total;
    double* dbuf = NULL;
    Expr* packed = p.machine ? ndbuild_open_f64((int64_t)n, &dbuf) : NULL;
    Expr** results = packed ? NULL : malloc(sizeof(Expr*) * n);
    if (!packed && !results) { free(stack); tx_free(&p); return NULL; }
    TxMpqStack mst = { NULL, (size_t)p.max_depth + 1, false };
    mpq_t mres;
    mpq_init(mres);

    for (size_t i = 0; i < n; i++) {
        /* imin + i*di stays within [imin, imax] (or [imax, imin]): no overflow */
        int64_t x = (int64_t)((uint64_t)imin + (uint64_t)i * (uint64_t)di);
        TxRat r;
        if (tx_run(&p, x, stack, &r)) {
            if (packed)         dbuf[i] = tx_to_double(r);
            else if (p.machine) results[i] = expr_new_real(tx_to_double(r));
            else                results[i] = (r.d == 1) ? expr_new_integer(r.n)
                                                        : make_rational(r.n, r.d);
            continue;
        }
        /* int64 overflowed (or an edge case): try exact GMP rationals. */
        if (tx_run_mpq(&p, x, &mst, mres)) {
            if (p.machine) {
                double dv = numeric_mpq_to_double(mres);
                /* Outside the double range numericalize answers an MPFR
                 * value: leave that to the interpreter. */
                if (isfinite(dv) && (dv != 0.0 || mpq_sgn(mres) == 0)) {
                    if (packed) dbuf[i] = dv;
                    else        results[i] = expr_new_real(dv);
                    continue;
                }
            } else {
                Expr* q = make_rational_mpz(mpq_numref(mres), mpq_denref(mres));
                if (q) { results[i] = q; continue; }   /* exact: never packed */
            }
        }
        /* This element is the interpreter's: bind and evaluate it. */
        Expr* xv = expr_new_integer(x);
        iter_bind(var, xv);
        expr_free(xv);
        Expr* v = evaluate(body);
        if (packed) {
            if (v && v->type == EXPR_REAL) { dbuf[i] = v->data.real; expr_free(v); continue; }
            /* Not a machine real: finish as a boxed List. */
            results = malloc(sizeof(Expr*) * n);
            if (!results) {
                expr_free(v); expr_free(packed); free(stack); tx_free(&p);
                tx_mpq_stack_free(&mst); mpq_clear(mres);
                return NULL;
            }
            ndbuild_abandon(packed, i, results);
            packed = NULL;
        }
        results[i] = v;
    }
    free(stack);
    tx_free(&p);
    tx_mpq_stack_free(&mst);
    mpq_clear(mres);
    if (packed) return packed;
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), results, n);
    free(results);
    return pack_offer(list);
}

Expr* builtin_table(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    
    if (res->data.function.arg_count > 2) {
        Expr** inner_args = malloc(sizeof(Expr*) * 2);
        inner_args[0] = expr_copy(res->data.function.args[0]);
        inner_args[1] = expr_copy(res->data.function.args[res->data.function.arg_count - 1]);
        Expr* inner_table = expr_new_function(expr_new_symbol(SYM_Table), inner_args, 2);
        free(inner_args);
        
        Expr** outer_args = malloc(sizeof(Expr*) * (res->data.function.arg_count - 1));
        outer_args[0] = inner_table;
        for (size_t i = 1; i < res->data.function.arg_count - 1; i++) {
            outer_args[i] = expr_copy(res->data.function.args[i]);
        }
        Expr* outer_table = expr_new_function(expr_new_symbol(SYM_Table), outer_args, res->data.function.arg_count - 1);
        free(outer_args);
        
        Expr* eval_outer = evaluate(outer_table);
        expr_free(outer_table);
        return eval_outer;
    }
    
    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* ---- Parse the iterator spec (shared helper) ---- */
    IterSpec s;
    if (!iter_spec_parse_lattice(spec, &s)) return NULL;

    int is_n_times   = (s.kind == ITER_KIND_COUNT);
    int is_list_iter = (s.kind == ITER_KIND_LIST);
    double min_val = 0, max_val = 0, di_val = 0;
    bool is_real = false, is_inf = false;

    /* Table does not iterate to Infinity; allow_inf = false. */
    if (!is_list_iter) {
        if (!iter_spec_resolve_numeric(&s, /*allow_inf=*/false,
                                       &min_val, &max_val, &di_val,
                                       &is_real, &is_inf)) {
            iter_spec_free(&s);
            return NULL;
        }
    }

    /* Convenience aliases into the owned IterSpec (freed via iter_spec_free). */
    Expr* var_sym = s.var;
    Expr* imin_e  = s.imin;
    Expr* imax_e  = s.imax;
    Expr* di_e    = s.di;
    Expr* list_e  = s.list;

    size_t results_cap = 16;
    size_t results_count = 0;
    Expr** results = malloc(sizeof(Expr*) * results_cap);
    bool overflow = false;   /* element count exceeded TABLE_MAX_FINITE_ELEMENTS */

    /* Exact-integer range longer than the backstop: decline up front, before
     * allocating a single element. The in-loop guard below (and Sum/Product)
     * only discover the overflow mid-fill; for an integer range the length is
     * O(1) to compute, so Table[i, {i, 1, 10^9}] returns unevaluated at once
     * instead of allocating a hundred million Exprs and then giving up. Real and
     * rational ranges fall through to the in-loop backstop. */
    if (!is_n_times && !is_list_iter && !is_real
        && expr_is_integer_like(imin_e) && expr_is_integer_like(imax_e)
        && di_e && expr_is_integer_like(di_e)) {
        mpz_t lo, hi, step, span;
        mpz_inits(lo, hi, step, span, NULL);
        expr_to_mpz(imin_e, lo);
        expr_to_mpz(imax_e, hi);
        expr_to_mpz(di_e, step);
        mpz_sub(span, hi, lo);                    /* imax - imin */
        bool decline = false;
        if (mpz_sgn(span) != 0 && mpz_sgn(span) == mpz_sgn(step)) {
            mpz_abs(span, span);
            mpz_abs(step, step);
            mpz_fdiv_q(span, span, step);         /* q = floor(|span|/|step|); count = q + 1 */
            decline = mpz_cmp_si(span, (long)TABLE_MAX_FINITE_ELEMENTS) >= 0;  /* count > cap */
        }
        mpz_clears(lo, hi, step, span, NULL);
        if (decline) { free(results); iter_spec_free(&s); return NULL; }
    }

    IterShadow old_own = iter_spec_shadow(var_sym);

    /* Exact integer iterator with a body in the exact-rational kernel's subset:
     * the kernel builds the whole list (see table_exact_int). */
    if (!is_n_times && !is_list_iter && !is_real && var_sym->type == EXPR_SYMBOL &&
        imin_e->type == EXPR_INTEGER && imax_e->type == EXPR_INTEGER &&
        di_e && di_e->type == EXPR_INTEGER) {
        Expr* fast = table_exact_int(expr, var_sym, imin_e->data.integer,
                                     imax_e->data.integer, di_e->data.integer);
        if (fast) {
            iter_spec_restore(var_sym, old_own);
            iter_spec_free(&s);
            free(results);
            return fast;
        }
    }

    if (is_n_times) {
        int64_t n = imax_e->data.integer;
        for (int64_t i = 0; i < n; i++) {
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else if (is_list_iter) {
        for (size_t i = 0; i < list_e->data.function.arg_count; i++) {
            iter_bind(var_sym, list_e->data.function.args[i]);
            Expr* eval_expr = evaluate(expr);
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;
        }
    } else {
        /* Auto-compile fast path — ONLY for an inexact (machine-real) iterator,
         * where the interpreter already produces machine reals, so exactness is
         * never at stake.  Exact (Integer/BigInt/Rational) iterators keep the
         * pure-interpreter path below and their results are bit-for-bit
         * unchanged.  A non-finite / complex compiled result falls back to the
         * interpreter for that element (which yields the complex/singular value).
         * When compiling, the exact running value `curr_e` is not needed, so its
         * per-element evaluate(Plus[...]) advance is skipped too. */
        AutoCompiled* ac = (is_real && var_sym->type == EXPR_SYMBOL)
                               ? autocompile_new(expr, (const Expr* const*)&var_sym, 1) : NULL;
        double val = min_val;
        Expr* curr_e = ac ? NULL : expr_copy(imin_e);

        /* Packed fast path: a compiled body whose result type is CT_REAL emits
         * nothing but machine reals, so the whole list can be written straight
         * into a buffer. Counted first with the loop's own recurrence -- doubles
         * only, no evaluation -- so the packed and interpreted paths produce the
         * same number of elements from the same sequence of `val`s.
         *
         * compiled_eval_real still self-guards finiteness, so an element where
         * the interpreter would give ComplexInfinity or a complex value returns
         * false. That abandons the buffer, hands back the elements already
         * written, and the ordinary loop below picks up from exactly where it
         * stopped -- same fallback the boxed path has always had, per element. */
        if (ac && autocompiled_result_is_real(ac)) {
            size_t total = 0;
            double v = min_val;
            /* Real packed path: is_real is true here, so the count uses the same
             * double recurrence as the fill below (no curr_e). */
            while (iter_range_continue(is_real, /*is_inf=*/false, NULL, NULL,
                                       v, max_val, di_val)) {
                if ((int64_t)total >= TABLE_MAX_FINITE_ELEMENTS) { overflow = true; break; }
                total++;
                v += di_val;
            }
            double* buf = NULL;
            Expr* packed = overflow ? NULL : ndbuild_open_f64((int64_t)total, &buf);
            if (packed) {
                size_t i = 0;
                while (i < total && autocompiled_eval_real(ac, &val, &buf[i])) {
                    i++;
                    val += di_val;
                }
                if (i == total) {
                    autocompiled_free(ac);
                    iter_spec_restore(var_sym, old_own);
                    iter_spec_free(&s);
                    free(results);
                    return packed;
                }
                /* Element i is not machine-real: finish on the List path. */
                if (total > results_cap) {
                    results_cap = total;
                    results = realloc(results, sizeof(Expr*) * results_cap);
                }
                ndbuild_abandon(packed, i, results);
                results_count = i;
            }
        }

        while (iter_range_continue(is_real, /*is_inf=*/false, curr_e, imax_e,
                                   val, max_val, di_val)) {
            if ((int64_t)results_count >= TABLE_MAX_FINITE_ELEMENTS) { overflow = true; break; }
            /* Boxed (not _eval_real) so an integer-valued body stays an Integer
             * — the element type is user-visible in the returned list. */
            Expr* eval_expr = ac ? autocompiled_eval_boxed(ac, &val) : NULL;
            if (!eval_expr) {   /* interpreter path (and per-element fallback) */
                Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
                iter_bind(var_sym, i_val);
                eval_expr = evaluate(expr);
                expr_free(i_val);
            }
            if (results_count == results_cap) { results_cap *= 2; results = realloc(results, sizeof(Expr*) * results_cap); }
            results[results_count++] = eval_expr;

            if (!ac) {   /* advance the exact running value (unused when compiling) */
                Expr* next_e = iter_step_add(curr_e, di_e);
                if (!next_e) {
                    Expr* next_args[2] = { expr_copy(curr_e), expr_copy(di_e) };
                    Expr* next_expr = expr_new_function(expr_new_symbol(SYM_Plus), next_args, 2);
                    next_e = evaluate(next_expr);
                    expr_free(next_expr);
                }
                expr_free(curr_e);
                curr_e = next_e;
            }

            val += di_val;
        }
        if (curr_e) expr_free(curr_e);
        autocompiled_free(ac);
    }

    /* A finite range whose element count exceeded the backstop: free what was
     * collected and leave Table[...] unevaluated, as Sum/Product do on overflow,
     * rather than returning a silently-truncated list. */
    if (overflow) {
        for (size_t i = 0; i < results_count; i++) expr_free(results[i]);
        free(results);
        iter_spec_restore(var_sym, old_own);
        iter_spec_free(&s);
        return NULL;
    }

    iter_spec_restore(var_sym, old_own);
    iter_spec_free(&s);

    Expr* result_list = expr_new_function(expr_new_symbol(SYM_List), results, results_count);
    free(results);
    /* Catches every branch the direct path above does not: the exact-iterator
     * loop, Table[expr, {n}], Table[expr, {i, list}], and a body that compiled to
     * an integer or complex result. Declines in O(1) on anything symbolic. */
    return pack_offer(result_list);
}
