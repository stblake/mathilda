/* Mathilda — Compile[] engine, scalar core (see compile.h, docs/design/compile.md).
 *
 * Front-end: bottom-up type inference over the Expr, monomorphic-opcode lowering
 * with widening coercions inserted where operand types differ, bailing (NULL) on
 * anything outside the compilable subset.  Back-end: a register machine whose
 * registers are raw 16-byte slots; the opcode carries the type, so execution
 * does no tag checks and the Real path never pays complex cost.  Temporaries use
 * a stack-discipline allocator, so a program needs O(expression-depth) registers.
 */
#include "compile.h"
#include "../arithmetic.h"
#include "../symtab.h"
#include "../ndarray.h"    /* NDUnaryKernel / NDBinaryKernel — shared kernel layer */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

/* ------------------------------------------------------------------ *
 *  Runtime slot, instruction, program                                 *
 * ------------------------------------------------------------------ */
typedef union { long long i; double r; double _Complex z; const void* p; } Slot;

/* scalar kernel signatures exposed by the shared ndkernels layer */
typedef bool (*kfn_r)(double, double*);
typedef bool (*kfn_c)(double, double, double*, double*);
typedef bool (*kfn_c2)(double, double, double, double, double*, double*);

enum {
    OP_CONST, OP_MOVE,
    OP_I2R, OP_I2C, OP_R2C,
    OP_ADD_I, OP_ADD_R, OP_ADD_C,
    OP_SUB_I, OP_SUB_R, OP_SUB_C,
    OP_MUL_I, OP_MUL_R, OP_MUL_C,
    OP_DIV_R, OP_DIV_C,
    OP_MOD_I, OP_QUOT_I,
    OP_NEG_I, OP_NEG_R, OP_NEG_C,
    OP_INV_R, OP_INV_C,
    OP_POWI_I, OP_POWI_R, OP_POWI_C, OP_POW_R, OP_POW_C,
    OP_SQRT_R, OP_SQRT_C, OP_EXP_R, OP_EXP_C, OP_LOG_R, OP_LOG_C,
    OP_SIN_R, OP_SIN_C, OP_COS_R, OP_COS_C, OP_TAN_R, OP_TAN_C,
    OP_SINH_R, OP_SINH_C, OP_COSH_R, OP_COSH_C, OP_TANH_R, OP_TANH_C,
    OP_ASIN_R, OP_ASIN_C, OP_ACOS_R, OP_ACOS_C, OP_ATAN_R, OP_ATAN_C,
    OP_ABS_I, OP_ABS_R, OP_ABS_C,   /* ABS_C -> real */
    OP_SIGN_I, OP_SIGN_R,
    OP_FLOOR_R, OP_CEIL_R, OP_ROUND_R,   /* real -> int */
    OP_RE_C, OP_IM_C, OP_ARG_C, OP_CONJ_C,   /* RE/IM/ARG -> real */
    OP_ATAN2_R, OP_MAX_I, OP_MAX_R, OP_MIN_I, OP_MIN_R,
    OP_ERF_R, OP_ERFC_R,
    /* generic special-function kernels from the shared ndkernels registry:
     * imm.p is the scalar kernel fn.  RR real->real; R2R real->real via cplx;
     * RC real->complex; CC complex->complex; CR complex->real (projection). */
    OP_KERN_RR, OP_KERN_R2R, OP_KERN_RC, OP_KERN_CC, OP_KERN_CR,
    OP_KERN2_RR, OP_KERN2_RC, OP_KERN2_CC,
    OP_LT_I, OP_LT_R, OP_LE_I, OP_LE_R, OP_GT_I, OP_GT_R, OP_GE_I, OP_GE_R,
    OP_EQ_I, OP_EQ_R, OP_EQ_C, OP_NE_I, OP_NE_R, OP_NE_C,
    OP_AND, OP_OR, OP_XOR, OP_NOT,
    OP_RET
};

typedef struct { uint16_t op; uint32_t dst, a, b; Slot imm; } Instr;

struct CompiledProgram {
    Instr*      code;
    size_t      n;
    int         nreg;
    int         result_reg;
    CompileType result_type;
    size_t      nargs;
    CompileType* arg_types;   /* [nargs] */
    unsigned char* argdep;    /* [nargs] which args are read */
    Slot*       frame;        /* reusable register file [nreg] (mutable via ptr) */
    bool        all_real;     /* every arg + result is CT_REAL */
};

/* ------------------------------------------------------------------ *
 *  Argument name -> index map (interned-pointer hash)                 *
 * ------------------------------------------------------------------ */
typedef struct { const char** key; int* val; size_t cap; } NameMap;
static bool nm_init(NameMap* m, const char* const* names, size_t n) {
    size_t cap = 8; while (cap < n * 2) cap <<= 1;
    m->cap = cap;
    m->key = calloc(cap, sizeof(const char*));
    m->val = malloc(cap * sizeof(int));
    if (!m->key || !m->val) { free(m->key); free(m->val); m->key = NULL; return false; }
    for (size_t k = 0; k < n; k++) {
        const char* nm = names[k];
        size_t h = ((uintptr_t)nm >> 4) & (cap - 1);
        while (m->key[h]) h = (h + 1) & (cap - 1);
        m->key[h] = nm; m->val[h] = (int)k;
    }
    return true;
}
static int nm_get(const NameMap* m, const char* nm) {
    size_t h = ((uintptr_t)nm >> 4) & (m->cap - 1);
    while (m->key[h]) { if (m->key[h] == nm) return m->val[h]; h = (h + 1) & (m->cap - 1); }
    return -1;
}
static void nm_free(NameMap* m) { free(m->key); free(m->val); }

/* ------------------------------------------------------------------ *
 *  Compiler context + emit                                            *
 * ------------------------------------------------------------------ */
typedef struct {
    Instr* code; size_t n, cap;
    int nlocals;        /* registers [0,nlocals) are args/locals */
    int temp_top;       /* current temp-stack height            */
    int maxreg;         /* high-water register count            */
    const CompileType* arg_types;
    NameMap map;
    unsigned char* argdep;
    bool ok;
} Ctx;

typedef struct { int reg; bool tmp; CompileType type; } Val;

static void ins(Ctx* c, uint16_t op, uint32_t dst, uint32_t a, uint32_t b, Slot imm) {
    if (!c->ok) return;
    if (c->n == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        Instr* nb = realloc(c->code, nc * sizeof(Instr));
        if (!nb) { c->ok = false; return; }
        c->code = nb; c->cap = nc;
    }
    c->code[c->n].op = op; c->code[c->n].dst = dst;
    c->code[c->n].a = a; c->code[c->n].b = b; c->code[c->n].imm = imm;
    c->n++;
}
static int alloc_temp(Ctx* c) {
    int r = c->nlocals + c->temp_top;
    c->temp_top++;
    if (c->nlocals + c->temp_top > c->maxreg) c->maxreg = c->nlocals + c->temp_top;
    return r;
}
static void free_if_tmp(Ctx* c, Val v) { if (v.tmp) c->temp_top--; }

static CompileType num_common(CompileType a, CompileType b) {
    if (a == CT_BOOL || b == CT_BOOL) return (CompileType)-1;
    return a > b ? a : b;
}

/* Widen v to `target` (numeric widening only), inserting a coercion op. */
static void coerce(Ctx* c, Val* v, CompileType target) {
    if (!c->ok || v->type == target) return;
    if (v->type == CT_BOOL || v->type > target) { c->ok = false; return; }
    uint16_t op = 0;
    if (v->type == CT_INT && target == CT_REAL) op = OP_I2R;
    else if (v->type == CT_INT && target == CT_COMPLEX) op = OP_I2C;
    else if (v->type == CT_REAL && target == CT_COMPLEX) op = OP_R2C;
    else { c->ok = false; return; }
    Slot z = { 0 };
    int dst = v->tmp ? v->reg : alloc_temp(c);
    ins(c, op, (uint32_t)dst, (uint32_t)v->reg, 0, z);
    v->reg = dst; v->tmp = true; v->type = target;
}

static bool emit(Ctx* c, const Expr* e, Val* out);

/* value in `a` (already at `type`), value in `b` (already at `type`): emit a
 * typed binary op, freeing operand temps LIFO and reusing a register. */
static Val binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype) {
    free_if_tmp(c, b);
    free_if_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z = { 0 };
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rtype };
    return r;
}
static Val unop(Ctx* c, uint16_t op, Val a, CompileType rtype) {
    free_if_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z = { 0 };
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, 0, z);
    Val r = { dst, true, rtype };
    return r;
}
static Val emit_const(Ctx* c, Slot imm, CompileType type) {
    int dst = alloc_temp(c);
    ins(c, OP_CONST, (uint32_t)dst, 0, 0, imm);
    Val r = { dst, true, type };
    return r;
}
/* unary/binary op carrying a kernel function pointer in imm.p */
static Val kern_unop(Ctx* c, uint16_t op, Val a, CompileType rt, const void* fn) {
    free_if_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, 0, z);
    Val r = { dst, true, rt };
    return r;
}
static Val kern_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, const void* fn) {
    free_if_tmp(c, b); free_if_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rt };
    return r;
}

/* numeric literal? -> value + type */
static bool literal(const Expr* e, Slot* imm, CompileType* type) {
    if (e->type == EXPR_INTEGER) { imm->i = e->data.integer; *type = CT_INT; return true; }
    if (e->type == EXPR_REAL)    { imm->r = e->data.real; *type = CT_REAL; return true; }
    if (e->type == EXPR_BIGINT)  { imm->r = mpz_get_d(e->data.bigint); *type = CT_REAL; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { imm->r = mpfr_get_d(e->data.mpfr, MPFR_RNDN); *type = CT_REAL; return true; }
#endif
    int64_t nn, dd;
    if (is_rational(e, &nn, &dd)) { imm->r = (double)nn / (double)dd; *type = CT_REAL; return true; }
    /* Complex[re, im] with numeric parts */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "Complex") == 0
        && e->data.function.arg_count == 2) {
        Slot ra, ia; CompileType ta, ti;
        if (literal(e->data.function.args[0], &ra, &ta) &&
            literal(e->data.function.args[1], &ia, &ti)) {
            double re = (ta == CT_INT) ? (double)ra.i : ra.r;
            double im = (ti == CT_INT) ? (double)ia.i : ia.r;
            imm->z = re + im * I; *type = CT_COMPLEX; return true;
        }
    }
    return false;
}

static bool named_const(const char* nm, double* out) {
    if (strcmp(nm, "Pi") == 0)          { *out = M_PI; return true; }
    if (strcmp(nm, "E") == 0)           { *out = M_E;  return true; }
    if (strcmp(nm, "EulerGamma") == 0)  { *out = 0.57721566490153286061; return true; }
    if (strcmp(nm, "Degree") == 0)      { *out = M_PI / 180.0; return true; }
    if (strcmp(nm, "GoldenRatio") == 0) { *out = 1.61803398874989484820; return true; }
    if (strcmp(nm, "Catalan") == 0)     { *out = 0.91596559417721901505; return true; }
    return false;
}

/* map a real/complex elementary unary head to its opcode pair; returns false if
 * not such a function. Sets *op_r, *op_c (op_c 0 if only real). */
static bool unary_math(const char* h, uint16_t* op_r, uint16_t* op_c) {
    #define U(name, r, cx) if (strcmp(h, name) == 0) { *op_r = r; *op_c = cx; return true; }
    U("Sqrt", OP_SQRT_R, OP_SQRT_C) U("Exp", OP_EXP_R, OP_EXP_C) U("Log", OP_LOG_R, OP_LOG_C)
    U("Sin", OP_SIN_R, OP_SIN_C) U("Cos", OP_COS_R, OP_COS_C) U("Tan", OP_TAN_R, OP_TAN_C)
    U("Sinh", OP_SINH_R, OP_SINH_C) U("Cosh", OP_COSH_R, OP_COSH_C) U("Tanh", OP_TANH_R, OP_TANH_C)
    U("ArcSin", OP_ASIN_R, OP_ASIN_C) U("ArcCos", OP_ACOS_R, OP_ACOS_C)
    U("Erf", OP_ERF_R, 0) U("Erfc", OP_ERFC_R, 0)
    #undef U
    return false;
}

/* emit a numeric unary function whose real/complex opcodes are op_r/op_c. */
static bool emit_unary_math(Ctx* c, const Expr* arg, uint16_t op_r, uint16_t op_c, Val* out) {
    Val a; if (!emit(c, arg, &a)) return false;
    if (a.type == CT_COMPLEX) {
        if (!op_c) { c->ok = false; return false; }
        *out = unop(c, op_c, a, CT_COMPLEX); return c->ok;
    }
    coerce(c, &a, CT_REAL);
    *out = unop(c, op_r, a, CT_REAL);
    return c->ok;
}

/* Generic special-function path: any numeric function registered in the shared
 * ndkernels layer (Gamma, LogGamma, Zeta, PolyGamma, Bessel*, ...) is lowered to
 * a KERNEL op that calls the machine kernel — no Expr, no opcode per function.
 * Returns true if a kernel handled the node; false otherwise (caller bails). */
static bool try_kernel(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (na == 1) {
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_unary_kernel) return false;
        const NDUnaryKernel* k = (const NDUnaryKernel*)d->ndarray_unary_kernel;
        if (!k->cplx && !k->real) return false;   /* degrade sentinel: no kernel */
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_COMPLEX) {
            if (!k->cplx) return false;           /* real-only kernel: can't do complex */
            *out = k->to_real ? kern_unop(c, OP_KERN_CR, a, CT_REAL, (const void*)k->cplx)
                              : kern_unop(c, OP_KERN_CC, a, CT_COMPLEX, (const void*)k->cplx);
            return true;
        }
        coerce(c, &a, CT_REAL);
        if (k->to_real && k->cplx)          *out = kern_unop(c, OP_KERN_R2R, a, CT_REAL, (const void*)k->cplx);
        else if (k->real_closed && k->real) *out = kern_unop(c, OP_KERN_RR, a, CT_REAL, (const void*)k->real);
        else if (k->real_closed && k->cplx) *out = kern_unop(c, OP_KERN_R2R, a, CT_REAL, (const void*)k->cplx);
        else if (k->real)                   *out = kern_unop(c, OP_KERN_RR, a, CT_REAL, (const void*)k->real);
        else if (k->cplx)                   *out = kern_unop(c, OP_KERN_RC, a, CT_COMPLEX, (const void*)k->cplx);
        else return false;
        return true;
    }
    if (na == 2) {
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_binary_kernel) return false;
        const NDBinaryKernel* k = (const NDBinaryKernel*)d->ndarray_binary_kernel;
        if (!k->cplx) return false;   /* degrade sentinel: no machine kernel -> bail */
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        if (t <= CT_REAL && k->real_closed) { coerce(c, &a, CT_REAL); coerce(c, &b, CT_REAL); *out = kern_binop(c, OP_KERN2_RR, a, b, CT_REAL, (const void*)k->cplx); }
        else if (t <= CT_REAL)              { coerce(c, &a, CT_REAL); coerce(c, &b, CT_REAL); *out = kern_binop(c, OP_KERN2_RC, a, b, CT_COMPLEX, (const void*)k->cplx); }
        else                                { coerce(c, &a, CT_COMPLEX); coerce(c, &b, CT_COMPLEX); *out = kern_binop(c, OP_KERN2_CC, a, b, CT_COMPLEX, (const void*)k->cplx); }
        return true;
    }
    return false;
}

static bool emit(Ctx* c, const Expr* e, Val* out) {
    if (!e || !c->ok) { c->ok = false; return false; }

    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) { *out = emit_const(c, imm, lt); return c->ok; }

    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        int k = nm_get(&c->map, nm);
        if (k >= 0) { c->argdep[k] = 1; out->reg = k; out->tmp = false; out->type = c->arg_types[k]; return true; }
        if (strcmp(nm, "True") == 0)  { imm.i = 1; *out = emit_const(c, imm, CT_BOOL); return c->ok; }
        if (strcmp(nm, "False") == 0) { imm.i = 0; *out = emit_const(c, imm, CT_BOOL); return c->ok; }
        if (strcmp(nm, "I") == 0)     { imm.z = I; *out = emit_const(c, imm, CT_COMPLEX); return c->ok; }
        double cv;
        if (named_const(nm, &cv)) { imm.r = cv; *out = emit_const(c, imm, CT_REAL); return c->ok; }
        c->ok = false; return false;
    }

    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) { c->ok = false; return false; }
    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args;
    size_t na = e->data.function.arg_count;

    /* n-ary Plus / Times */
    if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
        bool mul = h[0] == 'T';
        if (na == 0) { imm.i = mul ? 1 : 0; *out = emit_const(c, imm, CT_INT); return c->ok; }
        Val acc; if (!emit(c, A[0], &acc)) return false;
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return false;
            CompileType t = num_common(acc.type, b.type);
            if ((int)t < 0) { c->ok = false; return false; }
            coerce(c, &acc, t); coerce(c, &b, t);
            uint16_t op = mul ? (t == CT_INT ? OP_MUL_I : t == CT_REAL ? OP_MUL_R : OP_MUL_C)
                              : (t == CT_INT ? OP_ADD_I : t == CT_REAL ? OP_ADD_R : OP_ADD_C);
            acc = binop(c, op, acc, b, t);
        }
        *out = acc; return c->ok;
    }
    if ((strcmp(h, "Subtract") == 0) && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_INT ? OP_SUB_I : t == CT_REAL ? OP_SUB_R : OP_SUB_C, a, b, t);
        return c->ok;
    }
    if (strcmp(h, "Minus") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_BOOL) { c->ok = false; return false; }
        *out = unop(c, a.type == CT_INT ? OP_NEG_I : a.type == CT_REAL ? OP_NEG_R : OP_NEG_C, a, a.type);
        return c->ok;
    }
    if (strcmp(h, "Divide") == 0 && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        if (t < CT_REAL) t = CT_REAL;
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_DIV_C : OP_DIV_R, a, b, t);
        return c->ok;
    }
    if ((strcmp(h, "Mod") == 0 || strcmp(h, "Quotient") == 0) && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        if (a.type != CT_INT || b.type != CT_INT) { c->ok = false; return false; }
        *out = binop(c, h[0] == 'M' ? OP_MOD_I : OP_QUOT_I, a, b, CT_INT);
        return c->ok;
    }
    if (strcmp(h, "Power") == 0 && na == 2) {
        const Expr* base = A[0]; const Expr* ex = A[1];
        if (ex->type == EXPR_INTEGER) {
            long long nexp = ex->data.integer;
            Val a; if (!emit(c, base, &a)) return false;
            if (a.type == CT_INT && nexp >= 0) { Slot s; s.i = nexp; free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_I, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_INT; return c->ok; }
            if (a.type == CT_COMPLEX) { Slot s; s.i = nexp; free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_C, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_COMPLEX; return c->ok; }
            coerce(c, &a, CT_REAL);
            Slot s; s.i = nexp; free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_R, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_REAL; return c->ok;
        }
        int64_t rn, rd;
        if (is_rational(ex, &rn, &rd) && rd == 2 && (rn == 1 || rn == -1)) {
            Val a; if (!emit(c, base, &a)) return false;
            if (a.type == CT_COMPLEX) { a = unop(c, OP_SQRT_C, a, CT_COMPLEX); }
            else { coerce(c, &a, CT_REAL); a = unop(c, OP_SQRT_R, a, CT_REAL); }
            if (rn == -1) a = unop(c, a.type == CT_COMPLEX ? OP_INV_C : OP_INV_R, a, a.type);
            *out = a; return c->ok;
        }
        Val a, b; if (!emit(c, base, &a) || !emit(c, ex, &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        if (t < CT_REAL) t = CT_REAL;
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_POW_C : OP_POW_R, a, b, t);
        return c->ok;
    }

    /* elementary unary math */
    uint16_t op_r, op_c;
    if (na == 1 && unary_math(h, &op_r, &op_c)) return emit_unary_math(c, A[0], op_r, op_c, out);
    if (strcmp(h, "Log") == 0 && na == 2) {   /* Log[b,x] = Log[x]/Log[b] */
        Val x; if (!emit_unary_math(c, A[1], OP_LOG_R, OP_LOG_C, &x)) return false;
        Val b; if (!emit_unary_math(c, A[0], OP_LOG_R, OP_LOG_C, &b)) return false;
        CompileType t = num_common(x.type, b.type); coerce(c, &x, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_DIV_C : OP_DIV_R, x, b, t);
        return c->ok;
    }
    if (strcmp(h, "Tanh") == 0 && na == 1) return emit_unary_math(c, A[0], OP_TANH_R, OP_TANH_C, out);

    if (strcmp(h, "Abs") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_COMPLEX) { *out = unop(c, OP_ABS_C, a, CT_REAL); return c->ok; }
        if (a.type == CT_INT)     { *out = unop(c, OP_ABS_I, a, CT_INT); return c->ok; }
        coerce(c, &a, CT_REAL); *out = unop(c, OP_ABS_R, a, CT_REAL); return c->ok;
    }
    if (strcmp(h, "Sign") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_INT)  { *out = unop(c, OP_SIGN_I, a, CT_INT); return c->ok; }
        if (a.type == CT_REAL) { *out = unop(c, OP_SIGN_R, a, CT_REAL); return c->ok; }
        c->ok = false; return false;    /* complex Sign deferred */
    }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_INT) { *out = a; return c->ok; }
        if (a.type != CT_REAL) coerce(c, &a, CT_REAL);
        uint16_t op = h[0] == 'F' ? OP_FLOOR_R : h[0] == 'C' ? OP_CEIL_R : OP_ROUND_R;
        *out = unop(c, op, a, CT_INT); return c->ok;
    }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0 || strcmp(h, "Conjugate") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type != CT_COMPLEX) {
            if (strcmp(h, "Im") == 0) { free_if_tmp(c, a); imm.r = 0.0; *out = emit_const(c, imm, CT_REAL); return c->ok; }
            if (strcmp(h, "Arg") == 0) { /* real Arg: 0 if >=0 else Pi — deferred; require complex */ c->ok = false; return false; }
            *out = a; return c->ok;   /* Re/Conjugate of real = itself */
        }
        if (strcmp(h, "Re") == 0)   { *out = unop(c, OP_RE_C, a, CT_REAL); return c->ok; }
        if (strcmp(h, "Im") == 0)   { *out = unop(c, OP_IM_C, a, CT_REAL); return c->ok; }
        if (strcmp(h, "Arg") == 0)  { *out = unop(c, OP_ARG_C, a, CT_REAL); return c->ok; }
        *out = unop(c, OP_CONJ_C, a, CT_COMPLEX); return c->ok;
    }
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na >= 1) {
        bool mx = h[1] == 'a';
        Val acc; if (!emit(c, A[0], &acc)) return false;
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return false;
            CompileType t = num_common(acc.type, b.type);
            if ((int)t < 0 || t == CT_COMPLEX) { c->ok = false; return false; }
            coerce(c, &acc, t); coerce(c, &b, t);
            acc = binop(c, mx ? (t == CT_INT ? OP_MAX_I : OP_MAX_R) : (t == CT_INT ? OP_MIN_I : OP_MIN_R), acc, b, t);
        }
        *out = acc; return c->ok;
    }
    if (strcmp(h, "ArcTan") == 0) {
        if (na == 1) return emit_unary_math(c, A[0], OP_ATAN_R, OP_ATAN_C, out);
        if (na == 2) {
            Val x, y; if (!emit(c, A[0], &x) || !emit(c, A[1], &y)) return false;
            coerce(c, &x, CT_REAL); coerce(c, &y, CT_REAL);
            *out = binop(c, OP_ATAN2_R, x, y, CT_REAL); return c->ok;   /* atan2(y,x): b=y top */
        }
        c->ok = false; return false;
    }

    /* comparisons (2-arg) -> Bool */
    struct { const char* name; uint16_t oi, orr, oc; } CMP[] = {
        { "Less", OP_LT_I, OP_LT_R, 0 }, { "LessEqual", OP_LE_I, OP_LE_R, 0 },
        { "Greater", OP_GT_I, OP_GT_R, 0 }, { "GreaterEqual", OP_GE_I, OP_GE_R, 0 },
        { "Equal", OP_EQ_I, OP_EQ_R, OP_EQ_C }, { "Unequal", OP_NE_I, OP_NE_R, OP_NE_C },
    };
    for (size_t k = 0; k < sizeof(CMP) / sizeof(CMP[0]); k++) {
        if (strcmp(h, CMP[k].name) == 0 && na == 2) {
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
            CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
            coerce(c, &a, t); coerce(c, &b, t);
            uint16_t op = t == CT_COMPLEX ? CMP[k].oc : t == CT_INT ? CMP[k].oi : CMP[k].orr;
            if (!op) { c->ok = false; return false; }   /* ordering of complex */
            *out = binop(c, op, a, b, CT_BOOL); return c->ok;
        }
    }
    /* boolean logic */
    if ((strcmp(h, "And") == 0 || strcmp(h, "Or") == 0 || strcmp(h, "Xor") == 0) && na >= 1) {
        uint16_t op = h[0] == 'A' ? OP_AND : h[0] == 'O' ? OP_OR : OP_XOR;
        Val acc; if (!emit(c, A[0], &acc)) return false;
        if (acc.type != CT_BOOL) { c->ok = false; return false; }
        for (size_t i = 1; i < na; i++) {
            Val b; if (!emit(c, A[i], &b)) return false;
            if (b.type != CT_BOOL) { c->ok = false; return false; }
            acc = binop(c, op, acc, b, CT_BOOL);
        }
        *out = acc; return c->ok;
    }
    if (strcmp(h, "Not") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type != CT_BOOL) { c->ok = false; return false; }
        *out = unop(c, OP_NOT, a, CT_BOOL); return c->ok;
    }

    /* last resort: any numeric function with a machine kernel in ndkernels */
    Val kv;
    if (try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok; }

    c->ok = false; return false;   /* unsupported head -> bail */
}

/* ------------------------------------------------------------------ *
 *  VM                                                                 *
 * ------------------------------------------------------------------ */
static long long ipow_i(long long b, long long n) { long long r = 1; while (n > 0) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }
static double    ipow_r(double b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }
static double _Complex ipow_c(double _Complex b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double _Complex r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }

static void vm_run(const Instr* code, size_t n, Slot* R) {
    for (size_t k = 0; k < n; k++) {
        const Instr* c = &code[k];
        Slot* d = &R[c->dst]; const Slot* a = &R[c->a]; const Slot* b = &R[c->b];
        switch (c->op) {
            case OP_CONST: *d = c->imm; break;
            case OP_MOVE:  *d = *a; break;
            case OP_I2R: d->r = (double)a->i; break;
            case OP_I2C: d->z = (double)a->i; break;
            case OP_R2C: d->z = a->r; break;
            case OP_ADD_I: d->i = a->i + b->i; break;
            case OP_ADD_R: d->r = a->r + b->r; break;
            case OP_ADD_C: d->z = a->z + b->z; break;
            case OP_SUB_I: d->i = a->i - b->i; break;
            case OP_SUB_R: d->r = a->r - b->r; break;
            case OP_SUB_C: d->z = a->z - b->z; break;
            case OP_MUL_I: d->i = a->i * b->i; break;
            case OP_MUL_R: d->r = a->r * b->r; break;
            case OP_MUL_C: d->z = a->z * b->z; break;
            case OP_DIV_R: d->r = a->r / b->r; break;
            case OP_DIV_C: d->z = a->z / b->z; break;
            case OP_MOD_I: { long long m = b->i; long long q = a->i % m; if (q != 0 && ((q < 0) != (m < 0))) q += m; d->i = q; } break;
            case OP_QUOT_I: { long long m = b->i, x = a->i, q = x / m; if ((x % m != 0) && ((x < 0) != (m < 0))) q -= 1; d->i = q; } break;
            case OP_NEG_I: d->i = -a->i; break;
            case OP_NEG_R: d->r = -a->r; break;
            case OP_NEG_C: d->z = -a->z; break;
            case OP_INV_R: d->r = 1.0 / a->r; break;
            case OP_INV_C: d->z = 1.0 / a->z; break;
            case OP_POWI_I: d->i = ipow_i(a->i, c->imm.i); break;
            case OP_POWI_R: d->r = ipow_r(a->r, c->imm.i); break;
            case OP_POWI_C: d->z = ipow_c(a->z, c->imm.i); break;
            case OP_POW_R: d->r = pow(a->r, b->r); break;
            case OP_POW_C: d->z = cpow(a->z, b->z); break;
            case OP_SQRT_R: d->r = sqrt(a->r); break;
            case OP_SQRT_C: d->z = csqrt(a->z); break;
            case OP_EXP_R: d->r = exp(a->r); break;
            case OP_EXP_C: d->z = cexp(a->z); break;
            case OP_LOG_R: d->r = log(a->r); break;
            case OP_LOG_C: d->z = clog(a->z); break;
            case OP_SIN_R: d->r = sin(a->r); break;   case OP_SIN_C: d->z = csin(a->z); break;
            case OP_COS_R: d->r = cos(a->r); break;   case OP_COS_C: d->z = ccos(a->z); break;
            case OP_TAN_R: d->r = tan(a->r); break;   case OP_TAN_C: d->z = ctan(a->z); break;
            case OP_SINH_R: d->r = sinh(a->r); break; case OP_SINH_C: d->z = csinh(a->z); break;
            case OP_COSH_R: d->r = cosh(a->r); break; case OP_COSH_C: d->z = ccosh(a->z); break;
            case OP_TANH_R: d->r = tanh(a->r); break; case OP_TANH_C: d->z = ctanh(a->z); break;
            case OP_ASIN_R: d->r = asin(a->r); break; case OP_ASIN_C: d->z = casin(a->z); break;
            case OP_ACOS_R: d->r = acos(a->r); break; case OP_ACOS_C: d->z = cacos(a->z); break;
            case OP_ATAN_R: d->r = atan(a->r); break; case OP_ATAN_C: d->z = catan(a->z); break;
            case OP_ABS_I: d->i = a->i < 0 ? -a->i : a->i; break;
            case OP_ABS_R: d->r = fabs(a->r); break;
            case OP_ABS_C: d->r = cabs(a->z); break;
            case OP_SIGN_I: d->i = (a->i > 0) - (a->i < 0); break;
            case OP_SIGN_R: d->r = (a->r > 0) - (a->r < 0); break;
            case OP_FLOOR_R: d->i = (long long)floor(a->r); break;
            case OP_CEIL_R:  d->i = (long long)ceil(a->r); break;
            case OP_ROUND_R: d->i = (long long)llround(a->r); break;
            case OP_RE_C: d->r = creal(a->z); break;
            case OP_IM_C: d->r = cimag(a->z); break;
            case OP_ARG_C: d->r = carg(a->z); break;
            case OP_CONJ_C: d->z = conj(a->z); break;
            case OP_ATAN2_R: d->r = atan2(b->r, a->r); break;   /* ArcTan[x,y]=atan2(y,x); a=x,b=y */
            case OP_MAX_I: d->i = a->i > b->i ? a->i : b->i; break;
            case OP_MAX_R: d->r = a->r > b->r ? a->r : b->r; break;
            case OP_MIN_I: d->i = a->i < b->i ? a->i : b->i; break;
            case OP_MIN_R: d->r = a->r < b->r ? a->r : b->r; break;
            case OP_ERF_R: d->r = erf(a->r); break;
            case OP_ERFC_R: d->r = erfc(a->r); break;
            case OP_KERN_RR: { double o; d->r = ((kfn_r)c->imm.p)(a->r, &o) ? o : NAN; } break;
            case OP_KERN_R2R: { double orr, oi; d->r = ((kfn_c)c->imm.p)(a->r, 0.0, &orr, &oi) ? orr : NAN; } break;
            case OP_KERN_RC: { double orr, oi; if (((kfn_c)c->imm.p)(a->r, 0.0, &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } break;
            case OP_KERN_CC: { double orr, oi; if (((kfn_c)c->imm.p)(creal(a->z), cimag(a->z), &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } break;
            case OP_KERN_CR: { double orr, oi; d->r = ((kfn_c)c->imm.p)(creal(a->z), cimag(a->z), &orr, &oi) ? orr : NAN; } break;
            case OP_KERN2_RR: { double orr, oi; d->r = ((kfn_c2)c->imm.p)(a->r, 0.0, b->r, 0.0, &orr, &oi) ? orr : NAN; } break;
            case OP_KERN2_RC: { double orr, oi; if (((kfn_c2)c->imm.p)(a->r, 0.0, b->r, 0.0, &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } break;
            case OP_KERN2_CC: { double orr, oi; if (((kfn_c2)c->imm.p)(creal(a->z), cimag(a->z), creal(b->z), cimag(b->z), &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } break;
            case OP_LT_I: d->i = a->i < b->i; break;  case OP_LT_R: d->i = a->r < b->r; break;
            case OP_LE_I: d->i = a->i <= b->i; break; case OP_LE_R: d->i = a->r <= b->r; break;
            case OP_GT_I: d->i = a->i > b->i; break;  case OP_GT_R: d->i = a->r > b->r; break;
            case OP_GE_I: d->i = a->i >= b->i; break; case OP_GE_R: d->i = a->r >= b->r; break;
            case OP_EQ_I: d->i = a->i == b->i; break; case OP_EQ_R: d->i = a->r == b->r; break;
            case OP_EQ_C: d->i = a->z == b->z; break;
            case OP_NE_I: d->i = a->i != b->i; break; case OP_NE_R: d->i = a->r != b->r; break;
            case OP_NE_C: d->i = a->z != b->z; break;
            case OP_AND: d->i = a->i && b->i; break;
            case OP_OR:  d->i = a->i || b->i; break;
            case OP_XOR: d->i = (!!a->i) ^ (!!b->i); break;
            case OP_NOT: d->i = !a->i; break;
            case OP_RET: default: return;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Public API                                                         *
 * ------------------------------------------------------------------ */
CompiledProgram* compile_expr(const Expr* body, const char* const* arg_names,
                              const CompileType* arg_types, size_t nargs) {
    if (!body) return NULL;
    Ctx c; memset(&c, 0, sizeof(c));
    c.ok = true; c.nlocals = (int)nargs; c.arg_types = arg_types;
    c.argdep = calloc(nargs ? nargs : 1, 1);
    if (!c.argdep) return NULL;
    if (!nm_init(&c.map, arg_names, nargs)) { free(c.argdep); return NULL; }
    c.maxreg = (int)nargs;

    Val res;
    bool ok = emit(&c, body, &res) && c.ok;
    if (ok) { Slot z = { 0 }; ins(&c, OP_RET, (uint32_t)res.reg, 0, 0, z); ok = c.ok; }
    nm_free(&c.map);
    if (!ok) { free(c.code); free(c.argdep); return NULL; }

    CompiledProgram* p = calloc(1, sizeof(*p));
    if (!p) { free(c.code); free(c.argdep); return NULL; }
    p->code = c.code; p->n = c.n; p->nreg = c.maxreg;
    p->result_reg = res.reg; p->result_type = res.type;
    p->nargs = nargs; p->argdep = c.argdep;
    p->arg_types = malloc((nargs ? nargs : 1) * sizeof(CompileType));
    p->frame = malloc((size_t)(c.maxreg ? c.maxreg : 1) * sizeof(Slot));
    if (!p->arg_types || !p->frame) { compiled_free(p); return NULL; }
    memcpy(p->arg_types, arg_types, nargs * sizeof(CompileType));
    p->all_real = (res.type == CT_REAL);
    for (size_t k = 0; k < nargs; k++) if (arg_types[k] != CT_REAL) p->all_real = false;
    return p;
}

CompileType compiled_result_type(const CompiledProgram* p) { return p->result_type; }
size_t compiled_num_args(const CompiledProgram* p) { return p->nargs; }

size_t compiled_arg_deps(const CompiledProgram* p, int* deps, size_t cap) {
    size_t n = 0;
    for (size_t k = 0; k < p->nargs && n < cap; k++) if (p->argdep[k]) deps[n++] = (int)k;
    return n;
}

static void load_arg(Slot* s, const CompileValue* v, CompileType want) {
    /* coerce the boxed arg to the register's declared type */
    double _Complex z;
    switch (v->type) {
        case CT_BOOL: s->i = v->v.b ? 1 : 0; return;
        case CT_INT:  if (want == CT_INT) { s->i = v->v.i; } else if (want == CT_REAL) { s->r = (double)v->v.i; } else { s->z = (double)v->v.i; } return;
        case CT_REAL: if (want == CT_COMPLEX) { s->z = v->v.r; } else { s->r = v->v.r; } return;
        case CT_COMPLEX: z = v->v.z; s->z = z; return;
    }
}

static bool finite_result(const Slot* s, CompileType t) {
    if (t == CT_REAL) return isfinite(s->r);
    if (t == CT_COMPLEX) return isfinite(creal(s->z)) && isfinite(cimag(s->z));
    return true;   /* int/bool always ok */
}

bool compiled_eval(const CompiledProgram* p, const CompileValue* args, CompileValue* out) {
    for (size_t k = 0; k < p->nargs; k++) load_arg(&p->frame[k], &args[k], p->arg_types[k]);
    vm_run(p->code, p->n, p->frame);
    const Slot* r = &p->frame[p->result_reg];
    out->type = p->result_type;
    switch (p->result_type) {
        case CT_BOOL: out->v.b = (unsigned char)(r->i != 0); break;
        case CT_INT:  out->v.i = r->i; break;
        case CT_REAL: out->v.r = r->r; break;
        case CT_COMPLEX: out->v.z = r->z; break;
    }
    return finite_result(r, p->result_type);
}

bool compiled_eval_real(const CompiledProgram* p, const double* args, double* out) {
    if (!p->all_real) return false;
    for (size_t k = 0; k < p->nargs; k++) p->frame[k].r = args[k];
    vm_run(p->code, p->n, p->frame);
    *out = p->frame[p->result_reg].r;
    return isfinite(*out);
}

void compiled_free(CompiledProgram* p) {
    if (!p) return;
    free(p->code); free(p->arg_types); free(p->argdep); free(p->frame);
    free(p);
}
