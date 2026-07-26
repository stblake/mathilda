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
#include "../ndreduce.h"   /* ndred_total — array reductions (M3a) */
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
/* A register (or an instruction's immediate).  `p` carries an ndkernel function
 * pointer in an immediate; `arr` carries the OWNED EXPR_NDARRAY handle of an
 * array register (M3a).  The opcode says which member is live. */
typedef union { long long i; double r; double _Complex z; const void* p; Expr* arr; } Slot;

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
    OP_JMP,   /* pc = imm target in .b */
    OP_JZ,    /* if R[.a] (bool) is false, pc = .b; else fall through */
    OP_INC_I, /* R[.dst].i += imm.i (loop-counter step) */
    /* ---- arrays (M3a) -------------------------------------------------
     * An array register holds an owned EXPR_NDARRAY in Slot.arr and lives in a
     * dedicated bank at the top of the frame (see ARR_VREG), so a slot is
     * either always-array or never-array and teardown is one range sweep.
     * These ops delegate the buffer work to the NDArray subsystem; `flags`
     * says which operands are arrays, which the op consumes (frees), and what
     * element type the program promised. */
    OP_ARR_FREE,  /* expr_free R[dst]'s array and NULL the slot */
    OP_V_EW,      /* elementwise Plus (imm.i != 0) / Times, with scalar broadcast */
    OP_V_POW,     /* Power: array^array, array^scalar, scalar^array */
    OP_V_KERN,    /* map a unary ndkernel (imm.p) over the buffer */
    OP_V_KERN2,   /* binary ndkernel (imm.p) over one array + one scalar */
    OP_V_TOTAL,   /* full reduction of a rank-1 array -> scalar */
    OP_V_LEN,     /* leading dimension -> int */
    OP_RET
};

/* `flags` bit layout for the array opcodes.  `flags` occupies what was pure
 * padding after `op`, so Instr does not grow. */
#define AF_FREE_A     0x0001u   /* the op consumes (frees) R[a]'s array */
#define AF_FREE_B     0x0002u   /* the op consumes (frees) R[b]'s array */
#define AF_A_SHIFT    2         /* operand-a kind, 3 bits */
#define AF_B_SHIFT    5         /* operand-b kind, 3 bits */
#define AF_R_SHIFT    8         /* promised result element type, 2 bits */
#define AF_A(f)       (((f) >> AF_A_SHIFT) & 7u)
#define AF_B(f)       (((f) >> AF_B_SHIFT) & 7u)
#define AF_R(f)       (((f) >> AF_R_SHIFT) & 3u)
enum { AK_ARR = 0, AK_REAL = 1, AK_COMPLEX = 2 };   /* operand kinds */

typedef struct { uint16_t op, flags; uint32_t dst, a, b; Slot imm; } Instr;

struct CompiledProgram {
    Instr*      code;
    size_t      n;
    int         nreg;
    int         arr_base;     /* array registers are [arr_base, nreg) */
    int         result_reg;
    CompileType result_type;
    size_t      nargs;
    CompileType* arg_types;   /* [nargs] */
    unsigned char* argdep;    /* [nargs] which args are read */
    Slot*       frame;        /* reusable register file [nreg] (mutable via ptr) */
    bool        all_real;     /* every arg + result is CT_REAL, no array temps */
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
/* Array temporaries are allocated from their own stack into a *virtual*
 * register range tagged with ARR_VREG; at finalize the tags are rewritten to
 * `maxreg + k`, placing every array register in one contiguous bank above the
 * scalar registers.  That keeps array and scalar temps from ever sharing a
 * slot (so teardown can never mistake a double for a pointer) without a
 * second allocator or a liveness pass.  Jump targets, which also live in the
 * `b` field, are ordinary small integers and are left alone by the rewrite. */
#define ARR_VREG 0x40000000

typedef struct {
    Instr* code; size_t n, cap;
    int nlocals;        /* registers [0,nlocals) are args/locals */
    int temp_top;       /* current scalar temp-stack height      */
    int maxreg;         /* high-water scalar register count      */
    int arr_top;        /* current array temp-stack height       */
    int arr_max;        /* high-water array register count       */
    const CompileType* arg_types;
    NameMap map;
    unsigned char* argdep;
    /* lexically-scoped loop variables (Sum/Product/Nest), innermost last */
    struct { const char* name; int reg; CompileType type; } scope[16];
    int nscope;
    bool ok;
} Ctx;

/* resolve a symbol name to a scoped loop variable, or -1 */
static int scope_find(const Ctx* c, const char* nm, CompileType* type) {
    for (int s = c->nscope - 1; s >= 0; s--)
        if (c->scope[s].name == nm) { *type = c->scope[s].type; return c->scope[s].reg; }
    return -1;
}

typedef struct { int reg; bool tmp; CompileType type; } Val;

static void ins_f(Ctx* c, uint16_t op, uint16_t flags,
                  uint32_t dst, uint32_t a, uint32_t b, Slot imm) {
    if (!c->ok) return;
    if (c->n == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        Instr* nb = realloc(c->code, nc * sizeof(Instr));
        if (!nb) { c->ok = false; return; }
        c->code = nb; c->cap = nc;
    }
    c->code[c->n].op = op; c->code[c->n].flags = flags; c->code[c->n].dst = dst;
    c->code[c->n].a = a; c->code[c->n].b = b; c->code[c->n].imm = imm;
    c->n++;
}
static void ins(Ctx* c, uint16_t op, uint32_t dst, uint32_t a, uint32_t b, Slot imm) {
    ins_f(c, op, 0, dst, a, b, imm);
}
static int alloc_temp(Ctx* c) {
    int r = c->nlocals + c->temp_top;
    c->temp_top++;
    if (c->nlocals + c->temp_top > c->maxreg) c->maxreg = c->nlocals + c->temp_top;
    return r;
}
static int alloc_arr(Ctx* c) {
    int r = ARR_VREG + c->arr_top;
    c->arr_top++;
    if (c->arr_top > c->arr_max) c->arr_max = c->arr_top;
    return r;
}

/* Pop a temporary WITHOUT emitting anything: for operands whose array (if any)
 * the consuming instruction frees itself via its AF_FREE_* flags, so the free
 * happens after the operand has been read. */
static void pop_tmp(Ctx* c, Val v) {
    if (!v.tmp) return;
    if (CT_IS_ARRAY(v.type)) c->arr_top--; else c->temp_top--;
}
/* Pop a temporary whose value is now DEAD.  An array temp needs its handle
 * released here and now — inside a loop body the alternative (relying on
 * teardown) would accumulate one buffer per iteration. */
static void free_if_tmp(Ctx* c, Val v) {
    if (!v.tmp) return;
    if (CT_IS_ARRAY(v.type)) {
        Slot z = { 0 };
        ins(c, OP_ARR_FREE, (uint32_t)v.reg, 0, 0, z);
        c->arr_top--;
    } else c->temp_top--;
}

/* Common type of two operands.  An array absorbs a scalar (broadcast) and two
 * arrays must agree on rank; element types widen exactly as scalars do. */
static CompileType num_common(CompileType a, CompileType b) {
    if (a == CT_BOOL || b == CT_BOOL) return CT_ERR;
    if (CT_IS_ARRAY(a) || CT_IS_ARRAY(b)) {
        CompileType ea = CT_IS_ARRAY(a) ? CT_ELEM(a) : a;
        CompileType eb = CT_IS_ARRAY(b) ? CT_ELEM(b) : b;
        int ra = CT_IS_ARRAY(a) ? CT_RANK(a) : 0;
        int rb = CT_IS_ARRAY(b) ? CT_RANK(b) : 0;
        if (ra && rb && ra != rb) return CT_ERR;      /* no rank broadcasting */
        return CT_ARRAY(ea > eb ? ea : eb, ra > rb ? ra : rb);
    }
    return a > b ? a : b;
}

/* Widen v to `target` (numeric widening only), inserting a coercion op.  There
 * is no array coercion: the NDArray layer promotes element dtypes itself, so an
 * array only ever "coerces" to its own type. */
static void coerce(Ctx* c, Val* v, CompileType target) {
    if (!c->ok || v->type == target) return;
    if (CT_IS_ARRAY(v->type) || CT_IS_ARRAY(target)) { c->ok = false; return; }
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

/* Single choke point protecting every SCALAR opcode from an array operand.
 * Each scalar op is emitted through binop / unop / kern_unop / kern_binop, so
 * one guard here is enough: any head that has no array lowering (comparisons,
 * Max/Min, Mod, ...) bails automatically the moment an array reaches it,
 * instead of silently reinterpreting a handle as a double. */
static bool scalar_only(Ctx* c, CompileType a, CompileType b, CompileType r) {
    if (CT_IS_ARRAY(a) || CT_IS_ARRAY(b) || CT_IS_ARRAY(r)) { c->ok = false; return false; }
    return true;
}

/* value in `a` (already at `type`), value in `b` (already at `type`): emit a
 * typed binary op, freeing operand temps LIFO and reusing a register. */
static Val binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype) {
    scalar_only(c, a.type, b.type, rtype);
    pop_tmp(c, b);
    pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z = { 0 };
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rtype };
    return r;
}
static Val unop(Ctx* c, uint16_t op, Val a, CompileType rtype) {
    scalar_only(c, a.type, a.type, rtype);
    pop_tmp(c, a);
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
    scalar_only(c, a.type, a.type, rt);
    pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, 0, z);
    Val r = { dst, true, rt };
    return r;
}
static Val kern_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, const void* fn) {
    scalar_only(c, a.type, b.type, rt);
    pop_tmp(c, b); pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rt };
    return r;
}

/* ------------------------------------------------------------------ *
 *  Array emission helpers (M3a)                                       *
 * ------------------------------------------------------------------ */

/* Runtime operand kind of a compile-time type. */
static unsigned ak_of(CompileType t) {
    if (CT_IS_ARRAY(t)) return (unsigned)AK_ARR;
    return t == CT_COMPLEX ? (unsigned)AK_COMPLEX : (unsigned)AK_REAL;
}

/* Emit an array opcode.  The operand frees are encoded in the instruction's
 * flags rather than emitted as separate ARR_FREEs, so they happen *after* the
 * op has read the operands — which lets the result reuse an operand's register
 * exactly as the scalar binop does.  `b` is a dummy for unary ops. */
static Val arr_op(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, Slot imm) {
    CompileType relem = CT_IS_ARRAY(rt) ? CT_ELEM(rt) : rt;
    uint16_t f = (uint16_t)((ak_of(a.type) << AF_A_SHIFT) | (ak_of(b.type) << AF_B_SHIFT)
                            | (((unsigned)relem & 3u) << AF_R_SHIFT));
    if (a.tmp && CT_IS_ARRAY(a.type)) f |= AF_FREE_A;
    if (b.tmp && CT_IS_ARRAY(b.type)) f |= AF_FREE_B;
    pop_tmp(c, b); pop_tmp(c, a);
    int dst = CT_IS_ARRAY(rt) ? alloc_arr(c) : alloc_temp(c);
    ins_f(c, op, f, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, imm);
    Val r = { dst, true, rt };
    return r;
}

/* The unused second operand of a unary array op. */
static Val arr_noop_val(void) { Val v = { 0, false, CT_REAL }; return v; }

/* Prepare one operand of an array op: arrays pass through untouched (the ND
 * layer promotes element dtypes itself), scalars widen to Real/Complex so the
 * VM knows which half of the slot to read. */
static void arr_prep(Ctx* c, Val* v, CompileType elem) {
    if (CT_IS_ARRAY(v->type)) return;
    coerce(c, v, elem == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
}

/* A Real constant in a fresh scalar temp (exponents / negation factors). */
static Val arr_real_const(Ctx* c, double x) { Slot s; s.r = x; return emit_const(c, s, CT_REAL); }

/* Elementwise Plus/Times over operands of which at least one is an array. */
static Val arr_ew(Ctx* c, Val a, Val b, CompileType rt, bool is_plus) {
    arr_prep(c, &a, CT_ELEM(rt)); arr_prep(c, &b, CT_ELEM(rt));
    Slot ip; ip.i = is_plus ? 1 : 0;
    return arr_op(c, OP_V_EW, a, b, rt, ip);
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

/* f[array] for an already-emitted array operand `a`: map the head's registered
 * ndkernel over the buffer, so every elementary and special function that has a
 * machine kernel is reachable without one array opcode per function.  Sqrt has
 * no kernel (it is a Power in the ND layer) and lowers to the power op.
 *
 * Element typing mirrors the scalar rules: a real input promises a real result,
 * and a kernel that leaves the real axis (Sqrt/Log/ArcSin of a negative entry)
 * fails the call at runtime so the caller falls back to the interpreter —
 * the array-level form of the scalar "non-finite where the interpreter would be
 * complex" contract. */
static bool emit_arr_unary(Ctx* c, const char* head, Val a, Val* out) {
    CompileType ea = CT_ELEM(a.type);
    int rank = CT_RANK(a.type);
    if (strcmp(head, "Sqrt") == 0) {
        Val e = arr_real_const(c, 0.5);
        Slot z = { 0 };
        *out = arr_op(c, OP_V_POW, a, e, CT_ARRAY(ea, rank), z);
        return c->ok;
    }
    SymbolDef* d = symtab_lookup(head);
    const NDUnaryKernel* k = d ? (const NDUnaryKernel*)d->ndarray_unary_kernel : NULL;
    if (!k || (!k->cplx && !k->real)) { c->ok = false; return false; }   /* degrade sentinel */
    if (ea == CT_COMPLEX && !k->cplx) { c->ok = false; return false; }
    CompileType er = k->to_real ? CT_REAL : ea;
    Slot z; z.p = k;
    *out = arr_op(c, OP_V_KERN, a, arr_noop_val(), CT_ARRAY(er, rank), z);
    return c->ok;
}

/* emit a numeric unary function whose real/complex opcodes are op_r/op_c. */
static bool emit_unary_math(Ctx* c, const char* head, const Expr* arg,
                            uint16_t op_r, uint16_t op_c, Val* out) {
    Val a; if (!emit(c, arg, &a)) return false;
    if (CT_IS_ARRAY(a.type)) return emit_arr_unary(c, head, a, out);
    if (a.type == CT_COMPLEX) {
        if (!op_c) { c->ok = false; return false; }
        *out = unop(c, op_c, a, CT_COMPLEX); return c->ok;
    }
    coerce(c, &a, CT_REAL);
    *out = unop(c, op_r, a, CT_REAL);
    return c->ok;
}

/* Extract a single-parameter pure function's parameter name and body from
 * Function[u, body] or Function[{u}, body].  Returns false for anything else. */
static bool extract_function(const Expr* f, const char** pname, const Expr** body) {
    if (!f || f->type != EXPR_FUNCTION || f->data.function.head->type != EXPR_SYMBOL
        || strcmp(f->data.function.head->data.symbol.name, "Function") != 0
        || f->data.function.arg_count != 2) return false;
    const Expr* p = f->data.function.args[0];
    if (p->type == EXPR_SYMBOL) *pname = p->data.symbol.name;
    else if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL
             && strcmp(p->data.function.head->data.symbol.name, "List") == 0
             && p->data.function.arg_count == 1
             && p->data.function.args[0]->type == EXPR_SYMBOL)
        *pname = p->data.function.args[0]->data.symbol.name;
    else return false;
    *body = f->data.function.args[1];
    return true;
}

/* Fixed-point accumulator type for Nest[Function[u,body], x, n]: the body's
 * output feeds back as its input, so the accumulator register must hold a type
 * wide enough to absorb every iteration.  Starting from x's type, widen until
 * the body's output type no longer grows (bounded lattice → converges fast).
 * Returns the fixed-point CompileType, or -1 if it can't be compiled. */
static int nest_fixed_type(Ctx* c, const char* pname, const Expr* body, CompileType t0);

/* Pure type inference (no code emission) — needed to type an If's result
 * register before both branches are lowered.  Mirrors emit's result-type rules;
 * returns false for anything not compilable. */
static bool infer_type(Ctx* c, const Expr* e, CompileType* out) {
    if (!e) return false;
    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) { *out = lt; return true; }
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType st; if (scope_find(c, nm, &st) >= 0) { *out = st; return true; }
        int k = nm_get(&c->map, nm);
        if (k >= 0) { *out = c->arg_types[k]; return true; }
        if (strcmp(nm, "True") == 0 || strcmp(nm, "False") == 0) { *out = CT_BOOL; return true; }
        if (strcmp(nm, "I") == 0) { *out = CT_COMPLEX; return true; }
        double cv; if (named_const(nm, &cv)) { *out = CT_REAL; return true; }
        return false;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args; size_t na = e->data.function.arg_count;
    CompileType ta, tb;
    #define IT(idx, dst) do { if (!infer_type(c, A[idx], &dst)) return false; } while (0)
    if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
        if (na == 0) { *out = CT_INT; return true; }
        IT(0, ta); for (size_t i = 1; i < na; i++) { IT(i, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; }
        *out = ta; return true;
    }
    if (strcmp(h, "Subtract") == 0 && na == 2) { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; *out = ta; return true; }
    if (strcmp(h, "Minus") == 0 && na == 1)    { IT(0, ta); if (ta == CT_BOOL) return false; *out = ta; return true; }
    if (strcmp(h, "Divide") == 0 && na == 2)   { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; if (ta < CT_REAL) ta = CT_REAL; *out = ta; return true; }
    if ((strcmp(h, "Mod") == 0 || strcmp(h, "Quotient") == 0) && na == 2) { IT(0, ta); IT(1, tb); if (ta != CT_INT || tb != CT_INT) return false; *out = CT_INT; return true; }
    if (strcmp(h, "Power") == 0 && na == 2) {
        IT(0, ta);
        if (A[1]->type == EXPR_INTEGER) { *out = (ta == CT_INT && A[1]->data.integer >= 0) ? CT_INT : (ta == CT_COMPLEX ? CT_COMPLEX : CT_REAL); return true; }
        int64_t rn, rd; if (is_rational(A[1], &rn, &rd) && rd == 2 && (rn == 1 || rn == -1)) { *out = ta == CT_COMPLEX ? CT_COMPLEX : CT_REAL; return true; }
        IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false; if (ta < CT_REAL) ta = CT_REAL; *out = ta; return true;
    }
    uint16_t or_, oc_;
    if (na == 1 && (unary_math(h, &or_, &oc_) || strcmp(h, "Tanh") == 0)) { IT(0, ta); *out = ta == CT_COMPLEX ? CT_COMPLEX : CT_REAL; return true; }
    if (strcmp(h, "Log") == 0 && na == 2) { IT(0, ta); IT(1, tb); ta = num_common(ta, tb); *out = ta == CT_COMPLEX ? CT_COMPLEX : CT_REAL; return true; }
    if (strcmp(h, "Abs") == 0 && na == 1)  { IT(0, ta); *out = ta == CT_COMPLEX ? CT_REAL : ta; return true; }
    if (strcmp(h, "Sign") == 0 && na == 1) { IT(0, ta); if (ta == CT_COMPLEX) return false; *out = ta; return true; }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0) && na == 1) { IT(0, ta); *out = CT_INT; return true; }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0) && na == 1) { IT(0, ta); if (strcmp(h, "Arg") == 0 && ta != CT_COMPLEX) return false; *out = CT_REAL; return true; }
    if (strcmp(h, "Conjugate") == 0 && na == 1) { IT(0, ta); *out = ta; return true; }
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na >= 1) { IT(0, ta); for (size_t i = 1; i < na; i++) { IT(i, tb); ta = num_common(ta, tb); if ((int)ta < 0 || ta == CT_COMPLEX) return false; } *out = ta; return true; }
    if (strcmp(h, "ArcTan") == 0) { if (na == 1) { IT(0, ta); *out = ta == CT_COMPLEX ? CT_COMPLEX : CT_REAL; return true; } if (na == 2) { *out = CT_REAL; return true; } return false; }
    if (na == 2 && (!strcmp(h, "Less") || !strcmp(h, "LessEqual") || !strcmp(h, "Greater") || !strcmp(h, "GreaterEqual") || !strcmp(h, "Equal") || !strcmp(h, "Unequal"))) { *out = CT_BOOL; return true; }
    if (!strcmp(h, "And") || !strcmp(h, "Or") || !strcmp(h, "Xor") || !strcmp(h, "Not")) { *out = CT_BOOL; return true; }
    if (strcmp(h, "If") == 0 && na == 3) {
        CompileType tt, te; if (!infer_type(c, A[1], &tt) || !infer_type(c, A[2], &te)) return false;
        if (tt == te) { *out = tt; return true; }
        *out = num_common(tt, te); return (int)*out >= 0;
    }
    if ((strcmp(h, "Sum") == 0 || strcmp(h, "Product") == 0) && na == 2) {
        const Expr* spec = A[1];
        if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL
            || strcmp(spec->data.function.head->data.symbol.name, "List") != 0
            || spec->data.function.arg_count != 3
            || spec->data.function.args[0]->type != EXPR_SYMBOL || c->nscope >= 16) return false;
        CompileType t0, t1;
        if (!infer_type(c, spec->data.function.args[1], &t0) || !infer_type(c, spec->data.function.args[2], &t1)
            || t0 != CT_INT || t1 != CT_INT) return false;
        c->scope[c->nscope].name = spec->data.function.args[0]->data.symbol.name;
        c->scope[c->nscope].reg = 0; c->scope[c->nscope].type = CT_INT; c->nscope++;
        CompileType T; bool okT = infer_type(c, A[0], &T);
        c->nscope--;
        if (!okT || T == CT_BOOL) return false;
        *out = T; return true;
    }
    if (strcmp(h, "CompoundExpression") == 0 && na >= 1) return infer_type(c, A[na - 1], out);
    if ((strcmp(h, "With") == 0 || strcmp(h, "Module") == 0) && na == 2) {
        const Expr* L = A[0];
        if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
            || strcmp(L->data.function.head->data.symbol.name, "List") != 0
            || L->data.function.arg_count == 0
            || (int)(c->nscope + L->data.function.arg_count) > 16) return false;
        int pushed = 0;
        for (size_t i = 0; i < L->data.function.arg_count; i++) {
            const Expr* spec = L->data.function.args[i];
            const char* vname = NULL; CompileType vt = CT_REAL;
            if (spec->type == EXPR_SYMBOL) vname = spec->data.symbol.name;
            else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
                     && strcmp(spec->data.function.head->data.symbol.name, "Set") == 0
                     && spec->data.function.arg_count == 2 && spec->data.function.args[0]->type == EXPR_SYMBOL) {
                vname = spec->data.function.args[0]->data.symbol.name;
                if (!infer_type(c, spec->data.function.args[1], &vt)) { c->nscope -= pushed; return false; }
            } else { c->nscope -= pushed; return false; }
            c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = 0; c->scope[c->nscope].type = vt;
            c->nscope++; pushed++;
        }
        bool okb = infer_type(c, A[1], out);
        c->nscope -= pushed; return okb;
    }
    if ((!strcmp(h, "Set") || !strcmp(h, "AddTo") || !strcmp(h, "SubtractFrom") || !strcmp(h, "TimesBy")) && na == 2
        && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt) < 0) return false; *out = vt; return true;
    }
    if ((!strcmp(h, "Increment") || !strcmp(h, "Decrement")) && na == 1 && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt) < 0) return false; *out = vt; return true;
    }
    if ((!strcmp(h, "Do") && na == 2) || (!strcmp(h, "While") && na == 2) || (!strcmp(h, "For") && na == 4)) { *out = CT_INT; return true; }
    if (!strcmp(h, "Nest") && na == 3) {
        const char* pn; const Expr* bd; CompileType tn, tx;
        if (!extract_function(A[0], &pn, &bd) || !infer_type(c, A[2], &tn) || tn != CT_INT
            || !infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, pn, bd, tx);
        if (tfp < 0) return false;
        *out = (CompileType)tfp; return true;
    }
    /* array -> scalar reductions (M3a): the only way an array type re-enters
     * the scalar lattice, so these must be inferable inside If/Sum/... */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta) || CT_RANK(ta) != 1) return false;
        *out = (h[0] == 'L') ? CT_INT : CT_ELEM(ta);
        return true;
    }
    if (na == 1) { SymbolDef* d = symtab_lookup(h); if (d && d->ndarray_unary_kernel) { const NDUnaryKernel* k = d->ndarray_unary_kernel; if (k->cplx || k->real) { IT(0, ta); if (CT_IS_ARRAY(ta)) { *out = CT_ARRAY(k->to_real ? CT_REAL : CT_ELEM(ta), CT_RANK(ta)); return true; } if (k->to_real) { *out = CT_REAL; return true; } if (ta == CT_COMPLEX) { if (!k->cplx) return false; *out = CT_COMPLEX; return true; } *out = (k->real_closed || k->real) ? CT_REAL : CT_COMPLEX; return true; } } }
    if (na == 2) { SymbolDef* d = symtab_lookup(h); if (d && d->ndarray_binary_kernel) { const NDBinaryKernel* k = d->ndarray_binary_kernel; if (k->cplx) { IT(0, ta); IT(1, tb); CompileType t = num_common(ta, tb); if ((int)t < 0) return false; *out = (t <= CT_REAL && k->real_closed) ? CT_REAL : CT_COMPLEX; return true; } } }
    #undef IT
    return false;
}

static int nest_fixed_type(Ctx* c, const char* pname, const Expr* body, CompileType t0) {
    if (c->nscope >= 16) return -1;
    CompileType t = t0;
    for (int iter = 0; iter < 4; iter++) {
        c->scope[c->nscope].name = pname; c->scope[c->nscope].reg = 0;
        c->scope[c->nscope].type = t; c->nscope++;
        CompileType tb; bool okb = infer_type(c, body, &tb);
        c->nscope--;
        if (!okb || (int)tb < 0) return -1;
        if (tb == CT_BOOL && t != CT_BOOL) return -1;    /* can't fold a bool into a number */
        if (tb <= t) return (int)t;                       /* output coerces down into the accumulator */
        t = tb;                                           /* output widened it; grow and re-check */
    }
    return -1;
}

/* True when any argument is (or infers to) an array — the signal to take an
 * array lowering rather than the scalar one.  This is only a routing hint: a
 * wrong "false" still ends in a clean bail, because every scalar opcode is
 * guarded by scalar_only(). */
static bool any_array_arg(Ctx* c, Expr** A, size_t na) {
    for (size_t i = 0; i < na; i++) {
        CompileType t;
        if (infer_type(c, A[i], &t) && CT_IS_ARRAY(t)) return true;
    }
    return false;
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
        if (CT_IS_ARRAY(a.type)) return emit_arr_unary(c, h, a, out);
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
        if (CT_IS_ARRAY(t)) {
            /* one array + one broadcast scalar (BesselJ[n,v], ArcTan[v,y], ...);
             * the ND layer has no array-array binary-kernel map. */
            if (CT_IS_ARRAY(a.type) && CT_IS_ARRAY(b.type)) { c->ok = false; return false; }
            CompileType er = CT_ELEM(t);
            arr_prep(c, &a, er); arr_prep(c, &b, er);
            Slot z; z.p = k;
            *out = arr_op(c, OP_V_KERN2, a, b, CT_ARRAY(k->real_closed ? er : CT_COMPLEX,
                                                        CT_RANK(t)), z);
            return c->ok;
        }
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
        CompileType st; int sr = scope_find(c, nm, &st);
        if (sr >= 0) { out->reg = sr; out->tmp = false; out->type = st; return true; }
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
            if (CT_IS_ARRAY(t)) { acc = arr_ew(c, acc, b, t, !mul); if (!c->ok) return false; continue; }
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
        if (CT_IS_ARRAY(t)) {                       /* a - b = a + (-1) b */
            if (CT_IS_ARRAY(b.type)) {
                Val m1 = arr_real_const(c, -1.0);
                b = arr_ew(c, b, m1, b.type, false);
            } else {                                /* scalar subtrahend: negate in a register */
                arr_prep(c, &b, CT_ELEM(t));
                b = unop(c, b.type == CT_COMPLEX ? OP_NEG_C : OP_NEG_R, b, b.type);
            }
            if (!c->ok) return false;
            *out = arr_ew(c, a, b, t, true); return c->ok;
        }
        coerce(c, &a, t); coerce(c, &b, t);
        *out = binop(c, t == CT_INT ? OP_SUB_I : t == CT_REAL ? OP_SUB_R : OP_SUB_C, a, b, t);
        return c->ok;
    }
    if (strcmp(h, "Minus") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_BOOL) { c->ok = false; return false; }
        if (CT_IS_ARRAY(a.type)) {                  /* -v = (-1) v */
            Val m1 = arr_real_const(c, -1.0);
            *out = arr_ew(c, a, m1, a.type, false); return c->ok;
        }
        *out = unop(c, a.type == CT_INT ? OP_NEG_I : a.type == CT_REAL ? OP_NEG_R : OP_NEG_C, a, a.type);
        return c->ok;
    }
    if (strcmp(h, "Divide") == 0 && na == 2) {
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        if (CT_IS_ARRAY(t)) {
            /* a / b = a * b^-1: the ND layer has elementwise power and times,
             * but no divide. */
            CompileType et = CT_ELEM(t);
            if (CT_IS_ARRAY(b.type)) {
                Val em = arr_real_const(c, -1.0);
                Slot z = { 0 };
                b = arr_op(c, OP_V_POW, b, em, b.type, z);
            } else {
                arr_prep(c, &b, et);
                b = unop(c, b.type == CT_COMPLEX ? OP_INV_C : OP_INV_R, b, b.type);
            }
            if (!c->ok) return false;
            *out = arr_ew(c, a, b, t, false); return c->ok;
        }
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
        if (any_array_arg(c, A, na)) {
            /* base^exp over arrays: three shapes, all served by the NDArray
             * power fast paths (which promote to a complex dtype where a real
             * base leaves the real axis — caught by the element-type check). */
            Val a, b; if (!emit(c, base, &a) || !emit(c, ex, &b)) return false;
            CompileType t = num_common(a.type, b.type);
            if ((int)t < 0 || !CT_IS_ARRAY(t)) { c->ok = false; return false; }
            CompileType et = CT_ELEM(t);
            arr_prep(c, &a, et); arr_prep(c, &b, et);
            Slot z = { 0 };
            *out = arr_op(c, OP_V_POW, a, b, t, z);
            return c->ok;
        }
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
    if (na == 1 && unary_math(h, &op_r, &op_c)) return emit_unary_math(c, h, A[0], op_r, op_c, out);
    if (strcmp(h, "Log") == 0 && na == 2) {   /* Log[b,x] = Log[x]/Log[b] */
        /* over an array, the ND layer's registered two-arg Log kernel maps the
         * whole buffer in one pass — the scalar lowering below would need an
         * array divide the ND layer does not have. */
        Val kv2;
        if (any_array_arg(c, A, na)) {
            if (try_kernel(c, h, A, na, &kv2)) { *out = kv2; return c->ok; }
            c->ok = false; return false;
        }
        Val x; if (!emit_unary_math(c, h, A[1], OP_LOG_R, OP_LOG_C, &x)) return false;
        Val b; if (!emit_unary_math(c, h, A[0], OP_LOG_R, OP_LOG_C, &b)) return false;
        CompileType t = num_common(x.type, b.type); coerce(c, &x, t); coerce(c, &b, t);
        *out = binop(c, t == CT_COMPLEX ? OP_DIV_C : OP_DIV_R, x, b, t);
        return c->ok;
    }
    if (strcmp(h, "Tanh") == 0 && na == 1) return emit_unary_math(c, h, A[0], OP_TANH_R, OP_TANH_C, out);

    /* Projections and rounding over an array.  Each of these heads is handled
     * below by a dedicated scalar branch that would never reach try_kernel, yet
     * every one of them has a registered ND unary kernel — so intercept the
     * array case here and map the whole buffer in one pass. */
    if (na == 1 && (!strcmp(h, "Abs") || !strcmp(h, "Sign") || !strcmp(h, "Floor")
                    || !strcmp(h, "Ceiling") || !strcmp(h, "Round") || !strcmp(h, "Re")
                    || !strcmp(h, "Im") || !strcmp(h, "Arg") || !strcmp(h, "Conjugate"))
        && any_array_arg(c, A, na)) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (!CT_IS_ARRAY(a.type)) { c->ok = false; return false; }
        return emit_arr_unary(c, h, a, out);
    }

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
        if (na == 1) return emit_unary_math(c, h, A[0], OP_ATAN_R, OP_ATAN_C, out);
        if (na == 2) {
            Val kv2;
            if (any_array_arg(c, A, na)) {            /* ND two-arg ArcTan kernel */
                if (try_kernel(c, h, A, na, &kv2)) { *out = kv2; return c->ok; }
                c->ok = false; return false;
            }
            Val x, y; if (!emit(c, A[0], &x) || !emit(c, A[1], &y)) return false;
            coerce(c, &x, CT_REAL); coerce(c, &y, CT_REAL);
            *out = binop(c, OP_ATAN2_R, x, y, CT_REAL); return c->ok;   /* atan2(y,x): b=y top */
        }
        c->ok = false; return false;
    }

    /* Total[v] / Length[v]: the array -> scalar reductions.  Total delegates to
     * the NDArray reduction so its summation order — and therefore its
     * rounding — is identical to the interpreter's Total[]. */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (!CT_IS_ARRAY(a.type) || CT_RANK(a.type) != 1) { c->ok = false; return false; }
        Slot z = { 0 };
        *out = (h[0] == 'L') ? arr_op(c, OP_V_LEN, a, arr_noop_val(), CT_INT, z)
                             : arr_op(c, OP_V_TOTAL, a, arr_noop_val(), CT_ELEM(a.type), z);
        return c->ok;
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

    /* If[cond, then, else]: branch control flow.  Result lands in one register
     * that whichever branch runs writes; the other branch is jumped over. */
    if (strcmp(h, "If") == 0 && na == 3) {
        CompileType tt, te;
        if (!infer_type(c, A[1], &tt) || !infer_type(c, A[2], &te)) { c->ok = false; return false; }
        CompileType rt = (tt == te) ? tt : num_common(tt, te);
        /* An array result would have to be copied into the branch-join register;
         * a MOVE only duplicates the handle, so array-valued branches are not in
         * the M3a subset (see docs/design/compile_state.md). */
        if ((int)rt < 0 || CT_IS_ARRAY(rt)) { c->ok = false; return false; }
        int rr = alloc_temp(c);                     /* persistent result reg */
        Slot z = { 0 };
        Val cond; if (!emit(c, A[0], &cond)) return false;
        if (cond.type != CT_BOOL) { c->ok = false; return false; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)cond.reg, 0, z);
        free_if_tmp(c, cond);
        Val th; if (!emit(c, A[1], &th)) return false;
        if (CT_IS_ARRAY(th.type)) { c->ok = false; return false; }
        coerce(c, &th, rt);
        ins(c, OP_MOVE, (uint32_t)rr, (uint32_t)th.reg, 0, z);
        free_if_tmp(c, th);
        size_t jmp = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;   /* else label */
        Val el; if (!emit(c, A[2], &el)) return false;
        if (CT_IS_ARRAY(el.type)) { c->ok = false; return false; }
        coerce(c, &el, rt);
        ins(c, OP_MOVE, (uint32_t)rr, (uint32_t)el.reg, 0, z);
        free_if_tmp(c, el);
        if (c->ok) c->code[jmp].b = (uint32_t)c->n;   /* end label */
        out->reg = rr; out->tmp = true; out->type = rt;
        return c->ok;
    }

    /* Sum[body, {i, lo, hi}] / Product[...]: integer-counted accumulation loop.
     * The loop variable lives in a scoped register; the accumulator survives the
     * body's temporaries and is the result. */
    if ((strcmp(h, "Sum") == 0 || strcmp(h, "Product") == 0) && na == 2) {
        bool prod = h[0] == 'P';
        const Expr* spec = A[1];
        if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL
            || strcmp(spec->data.function.head->data.symbol.name, "List") != 0
            || spec->data.function.arg_count != 3
            || spec->data.function.args[0]->type != EXPR_SYMBOL) { c->ok = false; return false; }
        const char* iname = spec->data.function.args[0]->data.symbol.name;
        Expr* lo = spec->data.function.args[1];
        Expr* hi = spec->data.function.args[2];
        CompileType t0, t1;
        if (!infer_type(c, lo, &t0) || !infer_type(c, hi, &t1) || t0 != CT_INT || t1 != CT_INT
            || c->nscope >= 16) { c->ok = false; return false; }        /* integer iteration only */
        /* body type with i bound as INT */
        c->scope[c->nscope].name = iname; c->scope[c->nscope].reg = 0; c->scope[c->nscope].type = CT_INT; c->nscope++;
        CompileType T; bool okT = infer_type(c, A[0], &T);
        c->nscope--;
        /* An array accumulator would need a per-iteration copy, not a MOVE. */
        if (!okT || T == CT_BOOL || CT_IS_ARRAY(T)) { c->ok = false; return false; }
        Slot z = { 0 };
        int racc = alloc_temp(c), rhi = alloc_temp(c), ri = alloc_temp(c);
        Slot iz; iz.i = 0; if (T == CT_INT) iz.i = prod ? 1 : 0; else if (T == CT_REAL) iz.r = prod ? 1.0 : 0.0; else iz.z = prod ? 1.0 : 0.0;
        ins(c, OP_CONST, (uint32_t)racc, 0, 0, iz);
        Val vlo; if (!emit(c, lo, &vlo)) return false; ins(c, OP_MOVE, (uint32_t)ri, (uint32_t)vlo.reg, 0, z); free_if_tmp(c, vlo);
        Val vhi; if (!emit(c, hi, &vhi)) return false; ins(c, OP_MOVE, (uint32_t)rhi, (uint32_t)vhi.reg, 0, z); free_if_tmp(c, vhi);
        c->scope[c->nscope].name = iname; c->scope[c->nscope].reg = ri; c->scope[c->nscope].type = CT_INT; c->nscope++;
        size_t L = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;                                          /* free the guard temp */
        Val rb; if (!emit(c, A[0], &rb)) { c->nscope--; return false; }
        if (CT_IS_ARRAY(rb.type)) { c->nscope--; c->ok = false; return false; }
        coerce(c, &rb, T);
        uint16_t acc = prod ? (T == CT_INT ? OP_MUL_I : T == CT_REAL ? OP_MUL_R : OP_MUL_C)
                            : (T == CT_INT ? OP_ADD_I : T == CT_REAL ? OP_ADD_R : OP_ADD_C);
        ins(c, acc, (uint32_t)racc, (uint32_t)racc, (uint32_t)rb.reg, z);
        free_if_tmp(c, rb);
        Slot one; one.i = 1; ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)L, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;               /* loop-exit label */
        c->nscope--;
        c->temp_top -= 2;                                       /* free ri, rhi; racc is result */
        out->reg = racc; out->tmp = true; out->type = T;
        return c->ok;
    }

    /* ---- procedural constructs: local variables, mutation, loops ---- */

    /* CompoundExpression[e1,...,en]: run each for side effects, return the last. */
    if (strcmp(h, "CompoundExpression") == 0 && na >= 1) {
        for (size_t i = 0; i + 1 < na; i++) { Val v; if (!emit(c, A[i], &v)) return false; free_if_tmp(c, v); }
        return emit(c, A[na - 1], out);
    }

    /* With[{v=init,...}, body] / Module[{v or v=init,...}, body]: numeric locals
     * as mutable registers.  The body result is relocated onto the first local's
     * register so the locals free cleanly (LIFO). */
    if ((strcmp(h, "With") == 0 || strcmp(h, "Module") == 0) && na == 2) {
        const Expr* L = A[0];
        if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
            || strcmp(L->data.function.head->data.symbol.name, "List") != 0
            || L->data.function.arg_count == 0
            || (int)(c->nscope + L->data.function.arg_count) > 16) { c->ok = false; return false; }
        size_t nl = L->data.function.arg_count;
        Slot z = { 0 };
        int base_reg = -1, pushed = 0;
        for (size_t i = 0; i < nl; i++) {
            const Expr* spec = L->data.function.args[i];
            const char* vname = NULL; const Expr* init = NULL;
            if (spec->type == EXPR_SYMBOL) vname = spec->data.symbol.name;
            else if (spec->type == EXPR_FUNCTION && spec->data.function.head->type == EXPR_SYMBOL
                     && strcmp(spec->data.function.head->data.symbol.name, "Set") == 0
                     && spec->data.function.arg_count == 2
                     && spec->data.function.args[0]->type == EXPR_SYMBOL) {
                vname = spec->data.function.args[0]->data.symbol.name; init = spec->data.function.args[1];
            } else { c->nscope -= pushed; c->ok = false; return false; }
            int reg = alloc_temp(c);
            if (base_reg < 0) base_reg = reg;
            CompileType vt = CT_REAL;
            if (init) {
                Val iv; if (!emit(c, init, &iv)) { c->nscope -= pushed; return false; }
                /* an array local would alias, not own, its initialiser's handle */
                if (CT_IS_ARRAY(iv.type)) { c->nscope -= pushed; c->ok = false; return false; }
                vt = iv.type; ins(c, OP_MOVE, (uint32_t)reg, (uint32_t)iv.reg, 0, z); free_if_tmp(c, iv);
            } else { Slot s; s.r = 0.0; ins(c, OP_CONST, (uint32_t)reg, 0, 0, s); }
            c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = reg; c->scope[c->nscope].type = vt;
            c->nscope++; pushed++;
        }
        Val body; if (!emit(c, A[1], &body)) { c->nscope -= pushed; return false; }
        if (CT_IS_ARRAY(body.type)) { c->nscope -= pushed; c->ok = false; return false; }
        if (body.reg != base_reg) ins(c, OP_MOVE, (uint32_t)base_reg, (uint32_t)body.reg, 0, z);
        c->nscope -= pushed;
        c->temp_top = (base_reg - c->nlocals) + 1;   /* free above base_reg; keep result */
        out->reg = base_reg; out->tmp = true; out->type = body.type;
        return c->ok;
    }

    /* Set / AddTo / SubtractFrom / TimesBy on a scoped local (return new value);
     * Increment / Decrement (v++ / v--) return the OLD value. */
    {
        int kind = -1;
        if (strcmp(h, "Set") == 0) kind = 0; else if (strcmp(h, "AddTo") == 0) kind = 1;
        else if (strcmp(h, "SubtractFrom") == 0) kind = 2; else if (strcmp(h, "TimesBy") == 0) kind = 3;
        else if (strcmp(h, "Increment") == 0) kind = 4; else if (strcmp(h, "Decrement") == 0) kind = 5;
        if (kind >= 0) {
            size_t want = (kind >= 4) ? 1 : 2;
            if (na != want || A[0]->type != EXPR_SYMBOL) { c->ok = false; return false; }
            CompileType vt; int vreg = scope_find(c, A[0]->data.symbol.name, &vt);
            if (vreg < 0 || vt == CT_BOOL) { c->ok = false; return false; }   /* mutable numeric locals only */
            Slot z = { 0 };
            if (kind >= 4) {
                int old = alloc_temp(c); ins(c, OP_MOVE, (uint32_t)old, (uint32_t)vreg, 0, z);
                if (vt == CT_INT) { Slot s; s.i = (kind == 4) ? 1 : -1; ins(c, OP_INC_I, (uint32_t)vreg, 0, 0, s); }
                else { int one = alloc_temp(c); Slot s; s.r = (kind == 4) ? 1.0 : -1.0; ins(c, OP_CONST, (uint32_t)one, 0, 0, s);
                       ins(c, vt == CT_COMPLEX ? OP_ADD_C : OP_ADD_R, (uint32_t)vreg, (uint32_t)vreg, (uint32_t)one, z); c->temp_top--; }
                out->reg = old; out->tmp = true; out->type = vt; return c->ok;
            }
            Val val; if (!emit(c, A[1], &val)) return false;
            if (CT_IS_ARRAY(val.type)) { c->ok = false; return false; }
            coerce(c, &val, vt); if (!c->ok) return false;
            if (kind == 0) ins(c, OP_MOVE, (uint32_t)vreg, (uint32_t)val.reg, 0, z);
            else { uint16_t op = kind == 1 ? (vt == CT_INT ? OP_ADD_I : vt == CT_REAL ? OP_ADD_R : OP_ADD_C)
                              : kind == 2 ? (vt == CT_INT ? OP_SUB_I : vt == CT_REAL ? OP_SUB_R : OP_SUB_C)
                                          : (vt == CT_INT ? OP_MUL_I : vt == CT_REAL ? OP_MUL_R : OP_MUL_C);
                   ins(c, op, (uint32_t)vreg, (uint32_t)vreg, (uint32_t)val.reg, z); }
            free_if_tmp(c, val);
            out->reg = vreg; out->tmp = false; out->type = vt; return c->ok;
        }
    }

    /* Do[body, {i, lo, hi}]: counted loop for side effects; returns a dummy 0. */
    if (strcmp(h, "Do") == 0 && na == 2) {
        const Expr* spec = A[1];
        if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL
            || strcmp(spec->data.function.head->data.symbol.name, "List") != 0
            || spec->data.function.arg_count != 3
            || spec->data.function.args[0]->type != EXPR_SYMBOL || c->nscope >= 16) { c->ok = false; return false; }
        const char* iname = spec->data.function.args[0]->data.symbol.name;
        CompileType t0, t1;
        if (!infer_type(c, spec->data.function.args[1], &t0) || !infer_type(c, spec->data.function.args[2], &t1)
            || t0 != CT_INT || t1 != CT_INT) { c->ok = false; return false; }
        Slot z = { 0 };
        int rhi = alloc_temp(c), ri = alloc_temp(c);
        Val vlo; if (!emit(c, spec->data.function.args[1], &vlo)) return false; ins(c, OP_MOVE, (uint32_t)ri, (uint32_t)vlo.reg, 0, z); free_if_tmp(c, vlo);
        Val vhi; if (!emit(c, spec->data.function.args[2], &vhi)) return false; ins(c, OP_MOVE, (uint32_t)rhi, (uint32_t)vhi.reg, 0, z); free_if_tmp(c, vhi);
        c->scope[c->nscope].name = iname; c->scope[c->nscope].reg = ri; c->scope[c->nscope].type = CT_INT; c->nscope++;
        size_t Lp = c->n;
        int rc = alloc_temp(c); ins(c, OP_LE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
        Val bod; if (!emit(c, A[0], &bod)) { c->nscope--; return false; } free_if_tmp(c, bod);
        Slot one; one.i = 1; ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        c->nscope--; c->temp_top -= 2;
        int r0 = alloc_temp(c); Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT; return c->ok;
    }

    /* While[test, body] / For[start, test, incr, body]: data-dependent loops. */
    if ((strcmp(h, "While") == 0 && na == 2) || (strcmp(h, "For") == 0 && na == 4)) {
        bool isfor = h[0] == 'F';
        Slot z = { 0 };
        if (isfor) { Val s; if (!emit(c, A[0], &s)) return false; free_if_tmp(c, s); }
        size_t Lp = c->n;
        Val t; if (!emit(c, isfor ? A[1] : A[0], &t)) return false;
        if (t.type != CT_BOOL) { c->ok = false; return false; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)t.reg, 0, z); free_if_tmp(c, t);
        Val bod; if (!emit(c, isfor ? A[3] : A[1], &bod)) return false; free_if_tmp(c, bod);
        if (isfor) { Val ic; if (!emit(c, A[2], &ic)) return false; free_if_tmp(c, ic); }
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        int r0 = alloc_temp(c); Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT; return c->ok;
    }

    /* Nest[Function[u, body], x, n]: apply body n times, feeding each result
     * back in as u.  The accumulator lives in one persistent register (racc)
     * typed to the fixed-point type; the counted loop mirrors Do. */
    if (strcmp(h, "Nest") == 0 && na == 3) {
        const char* pname; const Expr* body;
        if (!extract_function(A[0], &pname, &body) || c->nscope >= 16) { c->ok = false; return false; }
        CompileType tn; if (!infer_type(c, A[2], &tn) || tn != CT_INT) { c->ok = false; return false; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
        int tfp = nest_fixed_type(c, pname, body, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return false; }
        CompileType t = (CompileType)tfp;
        Slot z = { 0 };
        /* Persistent registers FIRST, so freeing the init/count temps (which sit
         * above them) cannot clobber them. */
        int racc = alloc_temp(c), rn = alloc_temp(c), rcnt = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return false; coerce(c, &vx, t); if (!c->ok) return false;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        Val vn; if (!emit(c, A[2], &vn)) return false;
        ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
        Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, s0);
        size_t Lp = c->n;
        int rc = alloc_temp(c); ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
        c->scope[c->nscope].name = pname; c->scope[c->nscope].reg = racc; c->scope[c->nscope].type = t; c->nscope++;
        Val vb; if (!emit(c, body, &vb)) { c->nscope--; return false; }
        if (CT_IS_ARRAY(vb.type)) { c->nscope--; c->ok = false; return false; }
        coerce(c, &vb, t); if (!c->ok) { c->nscope--; return false; }
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        free_if_tmp(c, vb);
        c->nscope--;
        Slot one; one.i = 1; ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        c->temp_top = (racc - c->nlocals) + 1;   /* keep racc; free rn, rcnt */
        out->reg = racc; out->tmp = true; out->type = t; return c->ok;
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

/* ------------------------------------------------------------------ *
 *  Array opcodes (M3a): delegation to the NDArray subsystem           *
 * ------------------------------------------------------------------ */

/* Read a scalar register as a (re, im) pair, per its compile-time operand kind. */
static void vm_scalar_pair(const Slot* s, unsigned kind, double* re, double* im) {
    if (kind == AK_COMPLEX) { *re = creal(s->z); *im = cimag(s->z); }
    else                    { *re = s->r;        *im = 0.0; }
}

/* Box a scalar register as a temporary numeric Expr, because the ND helpers
 * take Expr operands.  Two allocations amortised over a whole-buffer pass is
 * nothing, and it keeps one implementation of the broadcast and dtype-promotion
 * rules instead of a second copy here. */
static Expr* vm_box_scalar(const Slot* s, unsigned kind) {
    double re, im;
    vm_scalar_pair(s, kind, &re, &im);
    if (im == 0.0) return expr_new_real(re);
    return make_complex(expr_new_real(re), expr_new_real(im));
}

/* Concrete (re, im) of a scalar Expr produced by an ND reduction. */
static bool vm_num_pair(const Expr* e, double* re, double* im) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *re = e->data.real;             *im = 0.0; return true; }
    if (e->type == EXPR_INTEGER) { *re = (double)e->data.integer;  *im = 0.0; return true; }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "Complex") == 0) {
        double ar, ai, br, bi;
        if (vm_num_pair(e->data.function.args[0], &ar, &ai)
            && vm_num_pair(e->data.function.args[1], &br, &bi)) {
            *re = ar; *im = br; return true;
        }
    }
    return false;
}

/* Write a reduction's scalar result into a register, honouring the element type
 * the program promised: a value that left the real axis where the program
 * promised a real one fails the call. */
static bool vm_write_scalar(const Expr* e, unsigned relem, Slot* d) {
    double re, im;
    if (!vm_num_pair(e, &re, &im)) return false;
    if (relem == (unsigned)CT_COMPLEX) { d->z = re + im * I; return true; }
    if (im != 0.0) return false;
    if (relem == (unsigned)CT_INT) d->i = (long long)re; else d->r = re;
    return true;
}

/* Execute one array opcode.  Returns false to abort the whole program — a shape
 * mismatch, an allocation failure, a kernel that declined an element, or a
 * result that left the promised element type — after which the caller falls
 * back to the interpreter.  Operands flagged AF_FREE_* are consumed here, AFTER
 * the op has read them, so the result may legitimately reuse an operand's
 * register; each freed slot is NULLed so an abort can never double-free. */
static bool vm_array_op(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const unsigned f = c->flags, ka = AF_A(f), kb = AF_B(f);
    Expr* r = NULL;

    switch (c->op) {
        case OP_ARR_FREE:
            expr_free(d->arr); d->arr = NULL;
            return true;

        case OP_V_LEN: {
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            long long len = (long long)x->data.ndarray.dims[0];
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            d->i = len;
            return true;
        }

        case OP_V_TOTAL: {
            Expr* s = ndred_total_all(a->arr);       /* borrows; same rounding as Total[] */
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            bool ok = s && vm_write_scalar(s, AF_R(f), d);
            expr_free(s);
            return ok;
        }

        case OP_V_EW: {
            Expr* boxa = (ka == AK_ARR) ? NULL : vm_box_scalar(a, ka);
            Expr* boxb = (kb == AK_ARR) ? NULL : vm_box_scalar(b, kb);
            Expr* ops[2];
            ops[0] = boxa ? boxa : a->arr;
            ops[1] = boxb ? boxb : b->arr;
            if (ops[0] && ops[1]) r = ndarray_elementwise(ops, 2, c->imm.i != 0);
            expr_free(boxa); expr_free(boxb);
            break;
        }

        case OP_V_POW: {
            double re, im;
            if (ka == AK_ARR && kb == AK_ARR)  r = ndarray_elementwise_power(a->arr, b->arr);
            else if (ka == AK_ARR) { vm_scalar_pair(b, kb, &re, &im); r = ndarray_scalar_power(a->arr, re, im); }
            else                   { vm_scalar_pair(a, ka, &re, &im); r = ndarray_base_scalar_power(re, im, b->arr); }
            break;
        }

        case OP_V_KERN:
            r = ndarray_map_unary(a->arr, (const NDUnaryKernel*)c->imm.p);
            break;

        case OP_V_KERN2: {
            Expr* boxa = (ka == AK_ARR) ? NULL : vm_box_scalar(a, ka);
            Expr* boxb = (kb == AK_ARR) ? NULL : vm_box_scalar(b, kb);
            const Expr* p0 = boxa ? boxa : a->arr;
            const Expr* p1 = boxb ? boxb : b->arr;
            if (p0 && p1) r = ndarray_map_binary(p0, p1, (const NDBinaryKernel*)c->imm.p);
            expr_free(boxa); expr_free(boxb);
            break;
        }

        default: return false;
    }

    /* shared tail for the array-PRODUCING ops */
    if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
    if (f & AF_FREE_B) { expr_free(b->arr); b->arr = NULL; }
    if (!r) return false;
    if (r->type != EXPR_NDARRAY
        || (AF_R(f) != (unsigned)CT_COMPLEX && ndt_is_complex(r->data.ndarray.dtype))) {
        expr_free(r); return false;
    }
    d->arr = r;
    return true;
}

/* Every opcode, for the computed-goto jump table (must cover the whole enum). */
#define OPLIST \
    X(JMP) X(JZ) X(INC_I) X(CONST) X(MOVE) X(I2R) X(I2C) X(R2C) \
    X(ADD_I) X(ADD_R) X(ADD_C) X(SUB_I) X(SUB_R) X(SUB_C) \
    X(MUL_I) X(MUL_R) X(MUL_C) X(DIV_R) X(DIV_C) X(MOD_I) X(QUOT_I) \
    X(NEG_I) X(NEG_R) X(NEG_C) X(INV_R) X(INV_C) \
    X(POWI_I) X(POWI_R) X(POWI_C) X(POW_R) X(POW_C) \
    X(SQRT_R) X(SQRT_C) X(EXP_R) X(EXP_C) X(LOG_R) X(LOG_C) \
    X(SIN_R) X(SIN_C) X(COS_R) X(COS_C) X(TAN_R) X(TAN_C) \
    X(SINH_R) X(SINH_C) X(COSH_R) X(COSH_C) X(TANH_R) X(TANH_C) \
    X(ASIN_R) X(ASIN_C) X(ACOS_R) X(ACOS_C) X(ATAN_R) X(ATAN_C) \
    X(ABS_I) X(ABS_R) X(ABS_C) X(SIGN_I) X(SIGN_R) \
    X(FLOOR_R) X(CEIL_R) X(ROUND_R) X(RE_C) X(IM_C) X(ARG_C) X(CONJ_C) \
    X(ATAN2_R) X(MAX_I) X(MAX_R) X(MIN_I) X(MIN_R) X(ERF_R) X(ERFC_R) \
    X(KERN_RR) X(KERN_R2R) X(KERN_RC) X(KERN_CC) X(KERN_CR) \
    X(KERN2_RR) X(KERN2_RC) X(KERN2_CC) \
    X(LT_I) X(LT_R) X(LE_I) X(LE_R) X(GT_I) X(GT_R) X(GE_I) X(GE_R) \
    X(EQ_I) X(EQ_R) X(EQ_C) X(NE_I) X(NE_R) X(NE_C) \
    X(AND) X(OR) X(XOR) X(NOT) \
    X(ARR_FREE) X(V_EW) X(V_POW) X(V_KERN) X(V_KERN2) X(V_TOTAL) X(V_LEN) \
    X(RET)

/* The bytecode interpreter.  A threaded (computed-goto) dispatch is used on
 * GCC/Clang — each opcode ends by jumping straight to the next, which the branch
 * predictor handles far better than a single switch; a portable switch is the
 * fallback.  Programs always end in OP_RET, so no per-op bounds check is needed. */
#if defined(__GNUC__) && !defined(VM_NO_THREADED)
#define VM_THREADED 1
#else
#define VM_THREADED 0
#endif

static void vm_run(const Instr* code, size_t n, Slot* R, bool* failed) {
    *failed = false;
    if (n == 0) return;
    size_t pc = 0;
    const Instr* c = &code[pc];
    Slot* d = &R[c->dst]; Slot* a = &R[c->a]; Slot* b = &R[c->b];
#if VM_THREADED
    #define X(name) [OP_##name] = &&L_##name,
    static const void* const tbl[] = { OPLIST };
    #undef X
    #define OP(name) L_##name
    #define NEXT() do { pc++; c = &code[pc]; d = &R[c->dst]; a = &R[c->a]; b = &R[c->b]; goto *tbl[c->op]; } while (0)
    #define JUMP() do { c = &code[pc]; d = &R[c->dst]; a = &R[c->a]; b = &R[c->b]; goto *tbl[c->op]; } while (0)
    goto *tbl[c->op];
#else
    #define OP(name) case OP_##name
    #define NEXT() break
    #define JUMP() continue
    while (pc < n) {
        c = &code[pc]; d = &R[c->dst]; a = &R[c->a]; b = &R[c->b];
        switch (c->op) {
#endif
    #define ARROP() do { if (!vm_array_op(c, d, a, b)) goto vm_fail; } while (0); NEXT()
            OP(JMP): pc = c->b; JUMP();
            OP(JZ):  pc = a->i ? pc + 1 : c->b; JUMP();   /* branch if false */
            OP(INC_I): d->i += c->imm.i; NEXT();
            OP(CONST): *d = c->imm; NEXT();
            OP(MOVE):  *d = *a; NEXT();
            OP(I2R): d->r = (double)a->i; NEXT();
            OP(I2C): d->z = (double)a->i; NEXT();
            OP(R2C): d->z = a->r; NEXT();
            OP(ADD_I): d->i = a->i + b->i; NEXT();
            OP(ADD_R): d->r = a->r + b->r; NEXT();
            OP(ADD_C): d->z = a->z + b->z; NEXT();
            OP(SUB_I): d->i = a->i - b->i; NEXT();
            OP(SUB_R): d->r = a->r - b->r; NEXT();
            OP(SUB_C): d->z = a->z - b->z; NEXT();
            OP(MUL_I): d->i = a->i * b->i; NEXT();
            OP(MUL_R): d->r = a->r * b->r; NEXT();
            OP(MUL_C): d->z = a->z * b->z; NEXT();
            OP(DIV_R): d->r = a->r / b->r; NEXT();
            OP(DIV_C): d->z = a->z / b->z; NEXT();
            OP(MOD_I): { long long m = b->i; long long q = a->i % m; if (q != 0 && ((q < 0) != (m < 0))) q += m; d->i = q; } NEXT();
            OP(QUOT_I): { long long m = b->i, x = a->i, q = x / m; if ((x % m != 0) && ((x < 0) != (m < 0))) q -= 1; d->i = q; } NEXT();
            OP(NEG_I): d->i = -a->i; NEXT();
            OP(NEG_R): d->r = -a->r; NEXT();
            OP(NEG_C): d->z = -a->z; NEXT();
            OP(INV_R): d->r = 1.0 / a->r; NEXT();
            OP(INV_C): d->z = 1.0 / a->z; NEXT();
            OP(POWI_I): d->i = ipow_i(a->i, c->imm.i); NEXT();
            OP(POWI_R): d->r = ipow_r(a->r, c->imm.i); NEXT();
            OP(POWI_C): d->z = ipow_c(a->z, c->imm.i); NEXT();
            OP(POW_R): d->r = pow(a->r, b->r); NEXT();
            OP(POW_C): d->z = cpow(a->z, b->z); NEXT();
            OP(SQRT_R): d->r = sqrt(a->r); NEXT();
            OP(SQRT_C): d->z = csqrt(a->z); NEXT();
            OP(EXP_R): d->r = exp(a->r); NEXT();
            OP(EXP_C): d->z = cexp(a->z); NEXT();
            OP(LOG_R): d->r = log(a->r); NEXT();
            OP(LOG_C): d->z = clog(a->z); NEXT();
            OP(SIN_R): d->r = sin(a->r); NEXT();   OP(SIN_C): d->z = csin(a->z); NEXT();
            OP(COS_R): d->r = cos(a->r); NEXT();   OP(COS_C): d->z = ccos(a->z); NEXT();
            OP(TAN_R): d->r = tan(a->r); NEXT();   OP(TAN_C): d->z = ctan(a->z); NEXT();
            OP(SINH_R): d->r = sinh(a->r); NEXT(); OP(SINH_C): d->z = csinh(a->z); NEXT();
            OP(COSH_R): d->r = cosh(a->r); NEXT(); OP(COSH_C): d->z = ccosh(a->z); NEXT();
            OP(TANH_R): d->r = tanh(a->r); NEXT(); OP(TANH_C): d->z = ctanh(a->z); NEXT();
            OP(ASIN_R): d->r = asin(a->r); NEXT(); OP(ASIN_C): d->z = casin(a->z); NEXT();
            OP(ACOS_R): d->r = acos(a->r); NEXT(); OP(ACOS_C): d->z = cacos(a->z); NEXT();
            OP(ATAN_R): d->r = atan(a->r); NEXT(); OP(ATAN_C): d->z = catan(a->z); NEXT();
            OP(ABS_I): d->i = a->i < 0 ? -a->i : a->i; NEXT();
            OP(ABS_R): d->r = fabs(a->r); NEXT();
            OP(ABS_C): d->r = cabs(a->z); NEXT();
            OP(SIGN_I): d->i = (a->i > 0) - (a->i < 0); NEXT();
            OP(SIGN_R): d->r = (a->r > 0) - (a->r < 0); NEXT();
            OP(FLOOR_R): d->i = (long long)floor(a->r); NEXT();
            OP(CEIL_R):  d->i = (long long)ceil(a->r); NEXT();
            OP(ROUND_R): d->i = (long long)llround(a->r); NEXT();
            OP(RE_C): d->r = creal(a->z); NEXT();
            OP(IM_C): d->r = cimag(a->z); NEXT();
            OP(ARG_C): d->r = carg(a->z); NEXT();
            OP(CONJ_C): d->z = conj(a->z); NEXT();
            OP(ATAN2_R): d->r = atan2(b->r, a->r); NEXT();   /* ArcTan[x,y]=atan2(y,x); a=x,b=y */
            OP(MAX_I): d->i = a->i > b->i ? a->i : b->i; NEXT();
            OP(MAX_R): d->r = a->r > b->r ? a->r : b->r; NEXT();
            OP(MIN_I): d->i = a->i < b->i ? a->i : b->i; NEXT();
            OP(MIN_R): d->r = a->r < b->r ? a->r : b->r; NEXT();
            OP(ERF_R): d->r = erf(a->r); NEXT();
            OP(ERFC_R): d->r = erfc(a->r); NEXT();
            OP(KERN_RR): { double o; d->r = ((kfn_r)c->imm.p)(a->r, &o) ? o : NAN; } NEXT();
            OP(KERN_R2R): { double orr, oi; d->r = ((kfn_c)c->imm.p)(a->r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN_RC): { double orr, oi; if (((kfn_c)c->imm.p)(a->r, 0.0, &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } NEXT();
            OP(KERN_CC): { double orr, oi; if (((kfn_c)c->imm.p)(creal(a->z), cimag(a->z), &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } NEXT();
            OP(KERN_CR): { double orr, oi; d->r = ((kfn_c)c->imm.p)(creal(a->z), cimag(a->z), &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RR): { double orr, oi; d->r = ((kfn_c2)c->imm.p)(a->r, 0.0, b->r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RC): { double orr, oi; if (((kfn_c2)c->imm.p)(a->r, 0.0, b->r, 0.0, &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } NEXT();
            OP(KERN2_CC): { double orr, oi; if (((kfn_c2)c->imm.p)(creal(a->z), cimag(a->z), creal(b->z), cimag(b->z), &orr, &oi)) d->z = orr + oi * I; else d->z = NAN + NAN * I; } NEXT();
            OP(LT_I): d->i = a->i < b->i; NEXT();  OP(LT_R): d->i = a->r < b->r; NEXT();
            OP(LE_I): d->i = a->i <= b->i; NEXT(); OP(LE_R): d->i = a->r <= b->r; NEXT();
            OP(GT_I): d->i = a->i > b->i; NEXT();  OP(GT_R): d->i = a->r > b->r; NEXT();
            OP(GE_I): d->i = a->i >= b->i; NEXT(); OP(GE_R): d->i = a->r >= b->r; NEXT();
            OP(EQ_I): d->i = a->i == b->i; NEXT(); OP(EQ_R): d->i = a->r == b->r; NEXT();
            OP(EQ_C): d->i = a->z == b->z; NEXT();
            OP(NE_I): d->i = a->i != b->i; NEXT(); OP(NE_R): d->i = a->r != b->r; NEXT();
            OP(NE_C): d->i = a->z != b->z; NEXT();
            OP(AND): d->i = a->i && b->i; NEXT();
            OP(OR):  d->i = a->i || b->i; NEXT();
            OP(XOR): d->i = (!!a->i) ^ (!!b->i); NEXT();
            OP(NOT): d->i = !a->i; NEXT();
            /* Array ops are out of line: they allocate, they can fail, and
             * keeping them out of the scalar cases costs the scalar path
             * nothing. */
            OP(ARR_FREE): ARROP();
            OP(V_EW):     ARROP();
            OP(V_POW):    ARROP();
            OP(V_KERN):   ARROP();
            OP(V_KERN2):  ARROP();
            OP(V_TOTAL):  ARROP();
            OP(V_LEN):    ARROP();
            OP(RET): return;
#if !VM_THREADED
            default: return;
        }
        pc++;
    }
    return;
#endif
vm_fail:
    *failed = true;   /* abort: the caller releases live array temporaries */
    return;
    #undef OP
    #undef NEXT
    #undef JUMP
    #undef ARROP
}
#undef OPLIST

/* ------------------------------------------------------------------ *
 *  Public API                                                         *
 * ------------------------------------------------------------------ */
/* Rewrite a virtual array register to its final slot in the array bank above
 * the scalar registers.  Ordinary register numbers and jump targets (also
 * carried in the `b` field) are far below ARR_VREG and pass through. */
static uint32_t patch_reg(uint32_t r, int base) {
    return r >= (uint32_t)ARR_VREG ? (uint32_t)base + (r - (uint32_t)ARR_VREG) : r;
}

CompiledProgram* compile_expr(const Expr* body, const char* const* arg_names,
                              const CompileType* arg_types, size_t nargs) {
    if (!body) return NULL;
    for (size_t k = 0; k < nargs; k++)            /* M3a: rank-1 arrays only */
        if (CT_IS_ARRAY(arg_types[k]) && CT_RANK(arg_types[k]) != 1) return NULL;
    Ctx c; memset(&c, 0, sizeof(c));
    c.ok = true; c.nlocals = (int)nargs; c.arg_types = arg_types;
    c.argdep = calloc(nargs ? nargs : 1, 1);
    if (!c.argdep) return NULL;
    if (!nm_init(&c.map, arg_names, nargs)) { free(c.argdep); return NULL; }
    c.maxreg = (int)nargs;

    Val res;
    bool ok = emit(&c, body, &res) && c.ok;
    /* A borrowed argument array cannot be the result: the caller owns what it
     * gets back, and freeing an argument would corrupt the caller's value. */
    if (ok && CT_IS_ARRAY(res.type) && !res.tmp) ok = false;
    if (ok) { Slot z = { 0 }; ins(&c, OP_RET, (uint32_t)res.reg, 0, 0, z); ok = c.ok; }
    nm_free(&c.map);
    if (!ok) { free(c.code); free(c.argdep); return NULL; }

    /* Place the array bank above the scalar registers and resolve the tags. */
    int arr_base = c.maxreg, nreg = c.maxreg + c.arr_max;
    for (size_t i = 0; i < c.n; i++) {
        c.code[i].dst = patch_reg(c.code[i].dst, arr_base);
        c.code[i].a   = patch_reg(c.code[i].a, arr_base);
        if (c.code[i].op != OP_JMP && c.code[i].op != OP_JZ)
            c.code[i].b = patch_reg(c.code[i].b, arr_base);
    }
    int result_reg = (int)patch_reg((uint32_t)res.reg, arr_base);

    CompiledProgram* p = calloc(1, sizeof(*p));
    if (!p) { free(c.code); free(c.argdep); return NULL; }
    p->code = c.code; p->n = c.n; p->nreg = nreg; p->arr_base = arr_base;
    p->result_reg = result_reg; p->result_type = res.type;
    p->nargs = nargs; p->argdep = c.argdep;
    p->arg_types = malloc((nargs ? nargs : 1) * sizeof(CompileType));
    p->frame = malloc((size_t)(nreg ? nreg : 1) * sizeof(Slot));
    if (!p->arg_types || !p->frame) { compiled_free(p); return NULL; }
    memcpy(p->arg_types, arg_types, nargs * sizeof(CompileType));
    p->all_real = (res.type == CT_REAL) && c.arr_max == 0;
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

/* Bind one boxed argument to its register.  Returns false if the value does not
 * match the declared type, so the caller can fall back to the interpreter.
 * An array argument is BORROWED: its handle is stored, never freed. */
static bool load_arg(Slot* s, const CompileValue* v, CompileType want) {
    if (CT_IS_ARRAY(want)) {
        Expr* x = CT_IS_ARRAY(v->type) ? v->v.a : NULL;
        if (!x || x->type != EXPR_NDARRAY || x->data.ndarray.rank != CT_RANK(want))
            return false;
        s->arr = x;
        return true;
    }
    if (CT_IS_ARRAY(v->type)) return false;
    /* coerce the boxed arg to the register's declared type */
    switch (v->type) {
        case CT_BOOL: s->i = v->v.b ? 1 : 0; return true;
        case CT_INT:  if (want == CT_INT) { s->i = v->v.i; } else if (want == CT_REAL) { s->r = (double)v->v.i; } else { s->z = (double)v->v.i; } return true;
        case CT_REAL: if (want == CT_COMPLEX) { s->z = v->v.r; } else { s->r = v->v.r; } return true;
        case CT_COMPLEX: s->z = v->v.z; return true;
        default: return false;
    }
}

static bool finite_result(const Slot* s, CompileType t) {
    if (t == CT_REAL) return isfinite(s->r);
    if (t == CT_COMPLEX) return isfinite(creal(s->z)) && isfinite(cimag(s->z));
    return true;   /* int/bool/array always ok */
}

/* Teardown for the array bank.  Every array register is NULLed before a call
 * (the frame is reused, so a handle left over from a previous call must never
 * be touched) and any that is still live afterwards is released here — the
 * belt to OP_ARR_FREE's braces, and the only cleanup on the abort path.
 * A result array is NULLed out by the caller first, so it survives. */
static void arr_reset(const CompiledProgram* p) {
    for (int r = p->arr_base; r < p->nreg; r++) p->frame[r].arr = NULL;
}
static void arr_sweep(const CompiledProgram* p) {
    for (int r = p->arr_base; r < p->nreg; r++)
        if (p->frame[r].arr) { expr_free(p->frame[r].arr); p->frame[r].arr = NULL; }
}

bool compiled_eval(const CompiledProgram* p, const CompileValue* args, CompileValue* out) {
    for (size_t k = 0; k < p->nargs; k++)
        if (!load_arg(&p->frame[k], &args[k], p->arg_types[k])) return false;
    arr_reset(p);
    bool failed = false;
    vm_run(p->code, p->n, p->frame, &failed);
    Slot* r = &p->frame[p->result_reg];
    out->type = p->result_type;
    if (CT_IS_ARRAY(p->result_type)) {
        out->v.a = failed ? NULL : r->arr;
        if (!failed) r->arr = NULL;    /* ownership transfers to the caller */
        arr_sweep(p);
        return !failed && out->v.a != NULL;
    }
    switch (p->result_type) {
        case CT_BOOL: out->v.b = (unsigned char)(r->i != 0); break;
        case CT_INT:  out->v.i = r->i; break;
        case CT_REAL: out->v.r = r->r; break;
        case CT_COMPLEX: out->v.z = r->z; break;
        default: break;
    }
    arr_sweep(p);
    return !failed && finite_result(r, p->result_type);
}

bool compiled_eval_real(const CompiledProgram* p, const double* args, double* out) {
    if (!p->all_real) return false;   /* implies no array registers */
    for (size_t k = 0; k < p->nargs; k++) p->frame[k].r = args[k];
    bool failed = false;
    vm_run(p->code, p->n, p->frame, &failed);
    *out = p->frame[p->result_reg].r;
    return isfinite(*out);
}

bool compiled_eval_real_batch(const CompiledProgram* const* progs, size_t nprogs,
                              const double* args, size_t nargs, double* out) {
    if (nprogs == 0) return true;
    /* one shared frame = the widest program's (big enough for every program's
     * registers); the argument region [0,nargs) is loaded once and never written
     * by any program (dst registers are always temporaries >= nargs). */
    size_t im = 0;
    for (size_t i = 1; i < nprogs; i++) if (progs[i]->nreg > progs[im]->nreg) im = i;
    Slot* F = progs[im]->frame;
    for (size_t k = 0; k < nargs; k++) F[k].r = args[k];
    for (size_t i = 0; i < nprogs; i++) {
        if (!progs[i]->all_real) return false;   /* implies no array registers */
        bool failed = false;
        vm_run(progs[i]->code, progs[i]->n, F, &failed);
        out[i] = F[(size_t)progs[i]->result_reg].r;
        if (!isfinite(out[i])) return false;
    }
    return true;
}

void compiled_free(CompiledProgram* p) {
    if (!p) return;
    free(p->code); free(p->arg_types); free(p->argdep); free(p->frame);
    free(p);
}
