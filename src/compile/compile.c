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
#include "compile_internal.h"    /* Slot / Instr / opcodes, shared with optimize.c */
#include "compiled_function.h"   /* inlining a CompiledFunction callee */
#include "../arithmetic.h"
#include "../symtab.h"
#include "../attr.h"      /* ATTR_LISTABLE — the gate on elementwise fusion */
#include "../ndarray.h"    /* NDUnaryKernel / NDBinaryKernel — shared kernel layer */
#include "../ndreduce.h"   /* ndred_total / ndred_accumulate — array reductions */
#include "../ndstruct.h"   /* ndstruct_reverse / _sort / ... — delegated structure */
#include "../ndarray_internal.h"  /* nd_parallel_for — threading the fused map loop */
#include "../print.h"      /* expr_to_string — printing the node a bail choked on */
#include "../sym_names.h" /* SYM_All / SYM_Span / SYM_List — Part subscript specs */
#include "../sym_intern.h" /* intern_symbol — the FN_HEAD placeholder parameters */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/* Thread-local storage for the VM's call-depth counter.  `_Thread_local` is C11
 * and `__thread` is a GNU extension, so both stay behind a guard; without either
 * the counter is a plain global, which is still correct for the single-threaded
 * build the project defaults to. */
#if defined(MATHILDA_THREADS) && (defined(__GNUC__) || defined(__clang__))
#define VM_TLS __thread
#else
#define VM_TLS
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

/* ------------------------------------------------------------------ *
 *  Runtime slot, instruction, program                                 *
 * ------------------------------------------------------------------ */
/* Slot, the kfn_* kernel signatures, Instr, the opcode enum, the AF_* array
 * flags, ParLoop and `struct CompiledProgram` itself all live in
 * compile_internal.h — shared with the optimiser and the disassembler. */

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
/* A compiled callee bigger than this is CALLed rather than pasted in: past this
 * size inlining costs more in code size and compile time than the call costs at
 * run time. */
#define INLINE_MAX_INSTRS 32

#define CSE_MAX       16     /* reserved CSE registers; bounds frame growth */
#define CSE_OCC_MAX  256     /* occurrence map entries */

/* Depth of the lexical scope stack.  One entry per live binding, and the
 * bindings now come from six places at once — Sum/Product iterators, Do
 * iterators, With/Module locals, inlined CompiledFunction parameters, fusion
 * leaves, and (since the functional heads landed) every lambda parameter of
 * every Nest/Fold/Map in the body.  Nesting three of those exhausted 16, and a
 * scope overflow bails the WHOLE body, so the array is sized for nesting rather
 * than for a single construct.  It lives on the stack-allocated Ctx, so the
 * cost is bytes at compile time and nothing at run time. */
#define CTX_MAX_SCOPE 32

/* Live Slot[] bindings — the parameters of a `Function[body]` being inlined.
 * Also bounds a lambda's parameter count and hence emit_apply's argument
 * count, so the three stay consistent by construction. */
#define FN_MAX_PARAMS 8

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
    struct { const char* name; int reg; CompileType type; bool built; } scope[CTX_MAX_SCOPE];
    int nscope;

    /* Live Slot[] bindings for an inlined `Function[body]` (`#`, `#1`, `#2`).
     *
     * Slot[k] is an EXPR_FUNCTION, not a symbol, so it cannot ride on scope[],
     * which matches interned symbol POINTERS.  One flat frame is not an
     * approximation: substitute_slots (src/purefunc.c:76) deliberately does not
     * recurse into a nested Function, so only the innermost frame is ever
     * visible, and emit_apply saves/restores nslot around each body. */
    struct { int reg; CompileType type; } slot[FN_MAX_PARAMS];
    int nslot;
    bool ok;
    unsigned flags;     /* COMPILE_FOLD_GLOBALS, ... */
    int inlining;       /* >0 while lowering an inlined CompiledFunction body */
    int fusing;         /* >0 while lowering the body of a fused elementwise loop */
    int tile_top, tile_max;   /* strip-mining tile registers (TILE_VREG bank) */
    bool vector_mode;   /* inside a strip-mined loop: tile-valued ops are legal */
    bool array_args;    /* any declared argument is an array: gates fusion probing */

    /* General-Part subscript lists (see PartSpec).  Held here while emitting and
     * handed to the program at finalize; freed here if the compile bails. */
    PartSpec** parts;
    int        nparts, parts_cap;

    /* Expr-level CSE (see cse_plan).  A chosen repeated subtree is computed once
     * into a register reserved below the temp stack, so no temp allocation and
     * none of the explicit temp_top resets in the loop lowerings can reach it. */
    struct { const Expr* node; int idx; } cse_occ[CSE_OCC_MAX];
    size_t      cse_nocc;
    int         cse_reg[CSE_MAX];
    CompileType cse_type[CSE_MAX];
    bool        cse_ready[CSE_MAX];
    int         ncse;

    /* Diagnostics: the INNERMOST subexpression `emit` could not lower.  A bail is
     * otherwise invisible — the caller just interprets, an order of magnitude
     * slower, and a body one head outside the subset costs the whole body.  Set
     * in the emit wrapper (first writer wins, so recursion unwinding leaves the
     * deepest cause) and rolled back with speculative lowering.  Borrowed: the
     * body outlives the compile. */
    const Expr* bail_node;
} Ctx;

/* The caller's argument symbols are NOT in scope inside an inlined
 * CompiledFunction body: the callee was compiled against its own parameters and
 * globals only, so a caller argument sharing a name with a callee global would
 * silently capture it. */
static int arg_find(const Ctx* c, const char* nm) {
    return c->inlining ? -1 : nm_get(&c->map, nm);
}

/* Under COMPILE_FOLD_GLOBALS, resolve a non-argument symbol that currently
 * holds a machine-number OwnValue (`y = 0.37`) to that constant — the same value
 * the interpreter would substitute for it right now.  Only an unconditional
 * assignment counts: the rule's pattern must be the bare symbol, so patterned or
 * conditional OwnValues are left alone and the symbol stays uncompilable.
 * See COMPILE_FOLD_GLOBALS in compile.h for why this is opt-in. */
static bool literal(const Expr* e, Slot* imm, CompileType* type);

static bool global_const(const Ctx* c, const char* nm, Slot* imm, CompileType* type) {
    if (!(c->flags & COMPILE_FOLD_GLOBALS)) return false;
    for (Rule* r = symtab_get_own_values(nm); r; r = r->next) {
        if (!r->pattern || r->pattern->type != EXPR_SYMBOL) continue;
        if (strcmp(r->pattern->data.symbol.name, nm) != 0) continue;
        return r->replacement && literal(r->replacement, imm, type);
    }
    return false;
}

/* A head symbol whose OwnValue is a CompiledFunction object — `newt[z, n]`
 * after `newt = Compile[...]`.  Applying it through the evaluator costs an
 * expression round-trip per call, so when the caller opted into folding we
 * inline the callee's body instead.  Same gate as global_const: the object could
 * be reassigned, which only matters for programs that outlive the call. */
static const CompiledFunction* compiled_callee(const Ctx* c, const char* nm) {
    if (!(c->flags & COMPILE_FOLD_GLOBALS)) return NULL;
    for (Rule* r = symtab_get_own_values(nm); r; r = r->next) {
        if (!r->pattern || r->pattern->type != EXPR_SYMBOL) continue;
        if (strcmp(r->pattern->data.symbol.name, nm) != 0) continue;
        return (r->replacement && r->replacement->type == EXPR_COMPILED)
             ? r->replacement->data.compiled : NULL;
    }
    return NULL;
}

/* resolve a symbol name to a scoped loop variable, or -1.  `built` is optional
 * and only meaningful for an array-typed binding (see Val.built). */
static int scope_find(const Ctx* c, const char* nm, CompileType* type, bool* built) {
    for (int s = c->nscope - 1; s >= 0; s--)
        if (c->scope[s].name == nm) {
            *type = c->scope[s].type;
            if (built) *built = c->scope[s].built;
            return c->scope[s].reg;
        }
    return -1;
}

/* A lowered value: its register, whether the consumer must free it, its type,
 * and — for an array — whether it was CONSTRUCTED by the body rather than
 * derived from an array argument.
 *
 * `built` exists because the result KIND has to match the interpreter's.  Given
 * a List the interpreter threads and returns a List; given an NDArray it
 * returns an NDArray; but a body that BUILDS its array (ConstantArray, Table,
 * NestList) returns a List *whatever the arguments were*, since the construct
 * itself has no packed form.  Deciding that at the boundary from the argument
 * kinds alone got it wrong for a body that takes an NDArray and builds a fresh
 * array from something else. */
typedef struct { int reg; bool tmp; CompileType type; bool built; } Val;

/* Result of an op over these operands: built unless some array operand traces
 * back to an argument.  An op with no array operand at all constructs its
 * result, so it is built. */
static bool arr_built(Val a, Val b) {
    bool r = true;
    if (CT_IS_ARRAY(a.type)) r = r && a.built;
    if (CT_IS_ARRAY(b.type)) r = r && b.built;
    return r;
}

static bool reg_is_tile(int r);

/* The tile opcodes occupy one contiguous run of OPLIST, from VSETLEN to the last
 * kernel form.  Used only by the safety net below. */
static bool is_tile_op(uint16_t op) {
    return op >= OP_VSETLEN && op <= OP_VKERN2_CC;
}

static void ins_f(Ctx* c, uint16_t op, uint16_t flags,
                  uint32_t dst, uint32_t a, uint32_t b, Slot imm) {
    if (!c->ok) return;
    /* SAFETY NET for strip mining.  A scalar opcode handed a tile register would
     * read the tile's POINTER as a double and quietly return nonsense — which is
     * exactly what happened when Power's integer-exponent path emitted POWI
     * directly instead of through unop().  Every lowering that can see a tile is
     * supposed to go through binop/unop/kern_unop/kern_binop; this catches any
     * that does not, and turns a wrong answer into a clean bail.  Control flow is
     * exempt: JMP/JZ/LOOP carry a jump target in `b`, not a register. */
    if (c->vector_mode && !is_tile_op(op)
        && op != OP_JMP && op != OP_JZ && op != OP_LOOP) {
        if (reg_is_tile((int)dst) || reg_is_tile((int)a) || reg_is_tile((int)b)) {
            c->ok = false;
            return;
        }
    }
    if (c->n == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 64;
        Instr* nb = realloc(c->code, nc * sizeof(Instr));
        if (!nb) { c->ok = false; return; }
        c->code = nb; c->cap = nc;
    }
    /* The overflow-checking choice is stamped onto each integer instruction here
     * — the one place every instruction passes through — rather than threaded
     * to the two dozen sites that emit integer arithmetic. */
    if ((c->flags & COMPILE_WRAP_INT) && op_is_checked_int(op)) flags |= IF_NOCHK;
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

/* A strip-mining tile register.  Like array registers these are allocated into
 * a virtual range and relocated into their own contiguous bank at finalize, so
 * "is this value a tile?" is answered by the register number alone — no extra
 * field on Val that a construction site could forget to initialise. */
static int alloc_tile(Ctx* c) {
    int r = TILE_VREG + c->tile_top;
    c->tile_top++;
    if (c->tile_top > c->tile_max) c->tile_max = c->tile_top;
    return r;
}
static bool reg_is_tile(int r) {
    return (uint32_t)r >= (uint32_t)TILE_VREG && (uint32_t)r < (uint32_t)ARR_VREG;
}
static bool val_is_tile(Val v) { return reg_is_tile(v.reg); }

/* Pop a temporary WITHOUT emitting anything: for operands whose array (if any)
 * the consuming instruction frees itself via its AF_FREE_* flags, so the free
 * happens after the operand has been read. */
static void pop_tmp(Ctx* c, Val v) {
    if (!v.tmp) return;
    if (CT_IS_ARRAY(v.type)) c->arr_top--;
    else if (val_is_tile(v)) c->tile_top--;
    else c->temp_top--;
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
    } else if (val_is_tile(v)) {
        c->tile_top--;      /* tile storage belongs to the frame, never freed */
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
    if (val_is_tile(*v)) {
        /* Tiles only ever hold Real or Complex elements, so the single widening
         * that can arise is real -> complex. */
        if (v->type != CT_REAL || target != CT_COMPLEX) { c->ok = false; return; }
        Slot z; memset(&z, 0, sizeof z);
        pop_tmp(c, *v);
        int dst = alloc_tile(c);
        ins(c, OP_VR2C, (uint32_t)dst, (uint32_t)v->reg, 0, z);
        v->reg = dst; v->tmp = true; v->type = CT_COMPLEX;
        return;
    }
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
static int cse_lookup(const Ctx* c, const Expr* e);

/* Single choke point protecting every SCALAR opcode from an array operand.
 * Each scalar op is emitted through binop / unop / kern_unop / kern_binop, so
 * one guard here is enough: any head that has no array lowering (comparisons,
 * Max/Min, Mod, ...) bails automatically the moment an array reaches it,
 * instead of silently reinterpreting a handle as a double. */
static bool scalar_only(Ctx* c, CompileType a, CompileType b, CompileType r) {
    if (CT_IS_ARRAY(a) || CT_IS_ARRAY(b) || CT_IS_ARRAY(r)) { c->ok = false; return false; }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Strip mining: scalar opcode -> tile opcode                         *
 * ------------------------------------------------------------------ *
 * Inside a strip-mined loop the ordinary lowering runs unchanged — Plus, Power,
 * the kernel dispatch, the coercions, all of it — and the four choke points
 * below swap in the tile form of whatever opcode it chose.  So the vectorisable
 * set is defined by this one table rather than by a second emitter that could
 * drift from the first.  An opcode with no tile form (comparisons, Floor, the
 * integer ops) returns 0 and the whole strip-mined attempt is rolled back. */
static uint16_t vec_op(uint16_t op) {
    switch (op) {
        case OP_ADD_R: return OP_VADD_R;   case OP_ADD_C: return OP_VADD_C;
        case OP_SUB_R: return OP_VSUB_R;   case OP_SUB_C: return OP_VSUB_C;
        case OP_MUL_R: return OP_VMUL_R;   case OP_MUL_C: return OP_VMUL_C;
        case OP_DIV_R: return OP_VDIV_R;   case OP_DIV_C: return OP_VDIV_C;
        case OP_NEG_R: return OP_VNEG_R;   case OP_NEG_C: return OP_VNEG_C;
        case OP_INV_R: return OP_VINV_R;   case OP_INV_C: return OP_VINV_C;
        case OP_POWI_R: return OP_VPOWI_R; case OP_POWI_C: return OP_VPOWI_C;
        case OP_POW_R: return OP_VPOW_R;   case OP_POW_C: return OP_VPOW_C;
        case OP_ATAN2_R: return OP_VATAN2_R;
        case OP_SQRT_R: return OP_VSQRT_R; case OP_SQRT_C: return OP_VSQRT_C;
        case OP_EXP_R: return OP_VEXP_R;   case OP_EXP_C: return OP_VEXP_C;
        case OP_LOG_R: return OP_VLOG_R;   case OP_LOG_C: return OP_VLOG_C;
        case OP_SIN_R: return OP_VSIN_R;   case OP_SIN_C: return OP_VSIN_C;
        case OP_COS_R: return OP_VCOS_R;   case OP_COS_C: return OP_VCOS_C;
        case OP_TAN_R: return OP_VTAN_R;   case OP_TAN_C: return OP_VTAN_C;
        case OP_SINH_R: return OP_VSINH_R; case OP_SINH_C: return OP_VSINH_C;
        case OP_COSH_R: return OP_VCOSH_R; case OP_COSH_C: return OP_VCOSH_C;
        case OP_TANH_R: return OP_VTANH_R; case OP_TANH_C: return OP_VTANH_C;
        case OP_ASIN_R: return OP_VASIN_R; case OP_ASIN_C: return OP_VASIN_C;
        case OP_ACOS_R: return OP_VACOS_R; case OP_ACOS_C: return OP_VACOS_C;
        case OP_ATAN_R: return OP_VATAN_R; case OP_ATAN_C: return OP_VATAN_C;
        case OP_ABS_R: return OP_VABS_R;   case OP_ABS_C: return OP_VABS_C;
        case OP_SIGN_R: return OP_VSIGN_R;
        case OP_RE_C: return OP_VRE_C;     case OP_IM_C: return OP_VIM_C;
        case OP_ARG_C: return OP_VARG_C;   case OP_CONJ_C: return OP_VCONJ_C;
        case OP_ERF_R: return OP_VERF_R;   case OP_ERFC_R: return OP_VERFC_R;
        case OP_R2C: return OP_VR2C;
        case OP_KERN_RR: return OP_VKERN_RR;   case OP_KERN_R2R: return OP_VKERN_R2R;
        case OP_KERN_RC: return OP_VKERN_RC;   case OP_KERN_CC: return OP_VKERN_CC;
        case OP_KERN_CR: return OP_VKERN_CR;
        case OP_KERN2_RR: return OP_VKERN2_RR; case OP_KERN2_RC: return OP_VKERN2_RC;
        case OP_KERN2_CC: return OP_VKERN2_CC;
        default: return 0;                 /* no tile form: abandon vectorisation */
    }
}

/* Broadcast a scalar operand into a tile so every tile op is tile-by-tile.  The
 * splat depends only on the scalar, so when that scalar is loop-invariant the
 * optimiser's LICM hoists it clean out of the strip loop. */
static Val vsplat(Ctx* c, Val v) {
    if (val_is_tile(v)) return v;
    if (CT_IS_ARRAY(v.type) || v.type == CT_BOOL || v.type == CT_INT) { c->ok = false; return v; }
    Slot z; memset(&z, 0, sizeof z);
    pop_tmp(c, v);
    int t = alloc_tile(c);
    ins(c, v.type == CT_COMPLEX ? OP_VSPLAT_C : OP_VSPLAT_R,
        (uint32_t)t, (uint32_t)v.reg, 0, z);
    Val r = { t, true, v.type, false };
    return r;
}

/* May a tile op write its result into an operand's register?
 *
 * Only when the element WIDTH is unchanged.  A complex-to-real op (Abs, Re, Im,
 * Arg) reads through `double _Complex*` and writes through `double*`; those are
 * different types, so the compiler is entitled to assume they cannot overlap and
 * to vectorise accordingly — which silently produces wrong results if the
 * emitter has pointed them at the same buffer.  Element-at-a-time the overlap
 * happens to be benign (the write always trails the read), which is exactly what
 * makes this the kind of bug that only appears once the loop vectorises.
 *
 * Not reusing costs one extra tile slot, and the whole tile bank is reclaimed at
 * the end of the loop body anyway. */
static bool tile_same_width(CompileType a, CompileType b) {
    return (a == CT_COMPLEX) == (b == CT_COMPLEX);
}

static Val vec_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype, Slot imm) {
    uint16_t vop = vec_op(op);
    if (!vop) { c->ok = false; return a; }
    a = vsplat(c, a); b = vsplat(c, b);
    if (tile_same_width(a.type, rtype) && tile_same_width(b.type, rtype)) {
        pop_tmp(c, b); pop_tmp(c, a);
    }
    int dst = alloc_tile(c);
    ins(c, vop, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, imm);
    Val r = { dst, true, rtype, false };
    return r;
}
static Val vec_unop(Ctx* c, uint16_t op, Val a, CompileType rtype, Slot imm) {
    uint16_t vop = vec_op(op);
    if (!vop) { c->ok = false; return a; }
    a = vsplat(c, a);
    if (tile_same_width(a.type, rtype)) pop_tmp(c, a);
    int dst = alloc_tile(c);
    ins(c, vop, (uint32_t)dst, (uint32_t)a.reg, 0, imm);
    Val r = { dst, true, rtype, false };
    return r;
}

/* value in `a` (already at `type`), value in `b` (already at `type`): emit a
 * typed binary op, freeing operand temps LIFO and reusing a register. */
static Val binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype) {
    if (c->vector_mode && (val_is_tile(a) || val_is_tile(b))) {
        Slot z; memset(&z, 0, sizeof z);
        return vec_binop(c, op, a, b, rtype, z);
    }
    scalar_only(c, a.type, b.type, rtype);
    pop_tmp(c, b);
    pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z = { 0 };
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rtype, false };
    return r;
}
static Val unop(Ctx* c, uint16_t op, Val a, CompileType rtype) {
    if (c->vector_mode && val_is_tile(a)) {
        Slot z; memset(&z, 0, sizeof z);
        return vec_unop(c, op, a, rtype, z);
    }
    scalar_only(c, a.type, a.type, rtype);
    pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z = { 0 };
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, 0, z);
    Val r = { dst, true, rtype, false };
    return r;
}
static Val emit_const(Ctx* c, Slot imm, CompileType type) {
    /* Normalise the immediate: `literal()` and friends assign only the member
     * matching the type, leaving the rest of the union indeterminate.  Zeroing
     * first makes two CONSTs of the same value bitwise identical, which is what
     * lets the optimiser's value numbering compare them. */
    Slot k; memset(&k, 0, sizeof k);
    switch (type) {
        case CT_BOOL: case CT_INT: k.i = imm.i; break;
        case CT_REAL:              k.r = imm.r; break;
        case CT_COMPLEX:           k.z = imm.z; break;
        default:                   k = imm;     break;
    }
    int dst = alloc_temp(c);
    ins(c, OP_CONST, (uint32_t)dst, 0, 0, k);
    Val r = { dst, true, type, false };
    return r;
}
/* unary/binary op carrying a kernel function pointer in imm.p */
static Val kern_unop(Ctx* c, uint16_t op, Val a, CompileType rt, const void* fn) {
    if (c->vector_mode && val_is_tile(a)) {
        Slot k; memset(&k, 0, sizeof k); k.p = fn;
        return vec_unop(c, op, a, rt, k);
    }
    scalar_only(c, a.type, a.type, rt);
    pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, 0, z);
    Val r = { dst, true, rt, false };
    return r;
}
static Val kern_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, const void* fn) {
    if (c->vector_mode && (val_is_tile(a) || val_is_tile(b))) {
        Slot k; memset(&k, 0, sizeof k); k.p = fn;
        return vec_binop(c, op, a, b, rt, k);
    }
    scalar_only(c, a.type, b.type, rt);
    pop_tmp(c, b); pop_tmp(c, a);
    int dst = alloc_temp(c);
    Slot z; z.p = fn;
    ins(c, op, (uint32_t)dst, (uint32_t)a.reg, (uint32_t)b.reg, z);
    Val r = { dst, true, rt, false };
    return r;
}

/* ------------------------------------------------------------------ *
 *  Array emission helpers (M3a)                                       *
 * ------------------------------------------------------------------ */

/* Runtime operand kind of a compile-time type. */
static unsigned ak_of(CompileType t) {
    if (CT_IS_ARRAY(t)) return (unsigned)AK_ARR;
    if (t == CT_COMPLEX) return (unsigned)AK_COMPLEX;
    /* An exact integer scalar stays exact: see AK_INT in compile_internal.h. A
     * Bool has no arithmetic meaning here and never reaches an array op. */
    if (t == CT_INT) return (unsigned)AK_INT;
    return (unsigned)AK_REAL;
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
    Val r = { dst, true, rt, arr_built(a, b) };
    return r;
}

/* The unused second operand of a unary array op. */
static Val arr_noop_val(void) { Val v = { 0, false, CT_REAL, false }; return v; }

/* Prepare one operand of an array op: arrays pass through untouched (the ND
 * layer promotes element dtypes itself), scalars widen to Real/Complex so the
 * VM knows which half of the slot to read. */
static void arr_prep(Ctx* c, Val* v, CompileType elem) {
    if (CT_IS_ARRAY(v->type)) return;
    /* An INTEGER-element result keeps an integer scalar exact rather than
     * widening it to Real -- the widening is what made u * 2 over an _Integer
     * array answer {2., 4., 6.}. Everything else widens as before: a Real
     * element type genuinely makes the whole operation inexact, which is the
     * interpreter's answer too (Range[1., 3.] * 2 is {2., 4., 6.}). */
    if (elem == CT_INT && v->type == CT_INT) return;
    coerce(c, v, elem == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
}

/* A Real constant in a fresh scalar temp (exponents / negation factors). */
static Val arr_real_const(Ctx* c, double x) { Slot s; s.r = x; return emit_const(c, s, CT_REAL); }

/* ------------------------------------------------------------------ *
 *  Indexed Part (M3c)                                                 *
 * ------------------------------------------------------------------ *
 * Part splits into two lowerings with very different costs, and which one
 * applies is decided purely by the SHAPE of the subscript list:
 *
 *   - every axis subscripted by a scalar integer expression, one per axis
 *     -> the flat index is built inline (one A_AXIS per axis) and the element
 *        is read with the same A_LOAD the fused elementwise loop uses.  No
 *        allocation, no call: this is the path a stencil runs in.
 *   - anything else — Span, All, a list of positions, or fewer subscripts than
 *     the rank — -> A_PART, which calls the interpreter's own ndarray_part.
 *     The result is a new array, so it costs an allocation, but the compiled
 *     answer is the interpreted one by construction rather than by agreement.
 *
 * The split is not a subset restriction: both together cover every spec Part
 * accepts on a dense array.
 */

static bool infer_type(Ctx* c, const Expr* e, CompileType* out);

/* A machine element type an array buffer can hold.  Bool has no buffer form and
 * a nested array is not an element; Int became one when NDT_INT64 arrived. */
static bool ct_is_elem(CompileType t) {
    return t == CT_INT || t == CT_REAL || t == CT_COMPLEX;
}

/* Is `node` (by identity) somewhere inside `root`?  Used only to keep a bail
 * diagnostic from pointing into a tree the emitter built and is about to free. */
static bool expr_subtree_of(const Expr* root, const Expr* node) {
    if (!root || !node) return false;
    if (root == node) return true;
    if (root->type != EXPR_FUNCTION) return false;
    if (expr_subtree_of(root->data.function.head, node)) return true;
    for (size_t i = 0; i < root->data.function.arg_count; i++)
        if (expr_subtree_of(root->data.function.args[i], node)) return true;
    return false;
}

/* Ownership: only an array the PROGRAM owns may be written through, because
 * A_SET writes the buffer in place.  Argument arrays are borrowed (the caller
 * still owns the node it passed in), and they live in the scalar register range
 * below nlocals, so the register number alone settles it. */
static bool reg_is_owned_arr(int r) { return (uint32_t)r >= (uint32_t)ARR_VREG; }

void compile_partspec_free(PartSpec* p) {
    if (!p) return;
    for (int i = 0; i < p->n; i++) if (p->lit) expr_free(p->lit[i]);
    free(p->lit); free(p->reg); free(p);
}

/* Hand a freshly built PartSpec to the context, which owns it from here on. */
static bool ctx_own_partspec(Ctx* c, PartSpec* p) {
    if (c->nparts == c->parts_cap) {
        int nc = c->parts_cap ? c->parts_cap * 2 : 4;
        PartSpec** np = realloc(c->parts, (size_t)nc * sizeof *np);
        if (!np) { compile_partspec_free(p); c->ok = false; return false; }
        c->parts = np; c->parts_cap = nc;
    }
    c->parts[c->nparts++] = p;
    return true;
}

/* A subscript that is a compile-time constant spec rather than a value to be
 * computed: an explicit All, a Span, or a list of positions.  Integers are
 * deliberately NOT included — an integer subscript is lowered as an expression
 * so that `u[[2]]` and `u[[i]]` take the same path. */
static bool subscript_is_literal_spec(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_All;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_Span || h == SYM_List;
}

/* True when `Part[a, A[0..na-1]]` can use the inline scalar path: one scalar
 * integer subscript per axis, none of them a slice spec. */
static bool part_is_scalar_indexed(Ctx* c, CompileType at, const Expr* const* A, size_t na) {
    if ((int)na != CT_RANK(at)) return false;           /* partial -> sub-array */
    for (size_t i = 0; i < na; i++) {
        if (subscript_is_literal_spec(A[i])) return false;
        CompileType t;
        if (!infer_type(c, A[i], &t) || t != CT_INT) return false;
    }
    return true;
}

/* Lower the subscripts of a scalar-indexed Part into ONE flat-index register.
 * `arr` must already be emitted.  Each axis costs a single A_AXIS, which does
 * the multiply by the axis length, the 1-based (or negative) resolution and the
 * range check together. */
static bool emit_flat_index(Ctx* c, Val arr, const Expr* const* A, size_t na, int* idx_out) {
    Slot z; memset(&z, 0, sizeof z);
    int ridx = alloc_temp(c);
    Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
    ins(c, OP_CONST, (uint32_t)ridx, 0, 0, k0);
    for (size_t i = 0; i < na; i++) {
        Val s;
        if (!emit(c, A[i], &s)) return false;
        if (s.type != CT_INT) { c->ok = false; return false; }
        Slot ax; memset(&ax, 0, sizeof ax); ax.i = (long long)i;
        ins(c, OP_A_AXIS, (uint32_t)ridx, (uint32_t)s.reg, (uint32_t)arr.reg, ax);
        free_if_tmp(c, s);
    }
    (void)z;
    *idx_out = ridx;
    return c->ok;
}

/* ConstantArray[v, n] / ConstantArray[v, {d1, ..., dr}] — the only way to bring
 * a new array into existence inside a compiled body.
 *
 * The RANK has to come from the source text (a bare dimension, or the length of
 * a literal dimension list), never from a runtime value, because it is part of
 * the result's compile-time type.  The dimensions themselves are ordinary
 * expressions and are evaluated per call. */
static bool const_array_shape(Ctx* c, const Expr* const* A, int* rank_out, CompileType* elem_out) {
    const Expr* d = A[1];
    int rank = 1;
    if (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL
        && d->data.function.head->data.symbol.name == SYM_List) {
        rank = (int)d->data.function.arg_count;
        if (rank < 1 || rank > CT_MAX_RANK || rank > NDARRAY_MAX_RANK) return false;
    }
    CompileType et;
    if (!infer_type(c, A[0], &et)) return false;
    /* An integer fill is fine now that NDT_INT64 exists: ConstantArray[0, n]
     * holds exact integer zeros in the interpreter, and a packed int64 buffer
     * holds exactly those.  Before the integer dtype this had to be refused,
     * because a float64 buffer would have answered differently rather than
     * merely faster. */
    if (!ct_is_elem(et)) return false;
    *rank_out = rank; *elem_out = et;
    return true;
}

/* Build the PartSpec for a general Part and emit any computed subscripts into a
 * run of consecutive registers, which is how A_PART finds them. */
static PartSpec* emit_partspec(Ctx* c, const Expr* const* A, size_t na, int* base_out) {
    PartSpec* p = calloc(1, sizeof *p);
    if (!p) { c->ok = false; return NULL; }
    p->n = (int)na;
    p->lit = calloc(na, sizeof *p->lit);
    p->reg = malloc(na * sizeof *p->reg);
    if (!p->lit || !p->reg) { compile_partspec_free(p); c->ok = false; return NULL; }
    for (size_t i = 0; i < na; i++) p->reg[i] = -1;

    /* Registers first, all of them, so the computed subscripts are contiguous:
     * a literal axis still burns a slot rather than perturbing the run. */
    int base = -1;
    for (size_t i = 0; i < na; i++) {
        int r = alloc_temp(c);
        if (base < 0) base = r;
        p->reg[i] = r;
    }
    for (size_t i = 0; i < na; i++) {
        if (subscript_is_literal_spec(A[i])) {
            p->lit[i] = expr_copy((Expr*)A[i]);
            if (!p->lit[i]) { compile_partspec_free(p); c->ok = false; return NULL; }
            continue;
        }
        Val s;
        if (!emit(c, A[i], &s) || s.type != CT_INT) {
            compile_partspec_free(p); c->ok = false; return NULL;
        }
        Slot z; memset(&z, 0, sizeof z);
        ins(c, OP_MOVE, (uint32_t)p->reg[i], (uint32_t)s.reg, 0, z);
        free_if_tmp(c, s);
    }
    *base_out = base;
    return p;
}

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

/* ---- compile-time function values ------------------------------------------
 *
 * The VM has no function value: no closure, no function register, and a lambda
 * is always INLINED at its use site.  What a functional head (Nest, Fold, Map,
 * ...) needs is therefore not a runtime representation but a compile-time one —
 * "which function is this, and how do I paste it in?".
 *
 * `fn_resolve` answers the first half, purely, so infer_type can ask it too;
 * `emit_apply` / `infer_apply` answer the second.  Keeping them in ONE place is
 * what stops two heads from disagreeing about which spellings of `f` they
 * accept, and it is why adding a functional head is a lowering rather than a
 * lowering plus a private parser for its function argument.
 *
 * The accepted spellings mirror apply_pure_function (src/purefunc.c:207).
 * Anything else bails, which is correct and not merely cautious: a form the
 * interpreter leaves symbolic is not a machine value, so answering with one
 * would diverge rather than merely go faster. */
typedef enum {
    FN_LAMBDA,    /* Function[u, body] / Function[{u, ...}, body] */
    FN_SLOTS,     /* Function[body] with #, #1, #2, ...           */
    FN_HEAD,      /* a bare head: Sin, Plus, or a CompiledFunction-valued symbol */
    FN_IDENTITY,  /* Identity                                     */
    FN_COMPOSE    /* Composition[f1, ..., fn]                     */
} FnKind;

typedef struct {
    FnKind      kind;
    int         nparams;                /* -1 == accepts any arity */
    const char* pname[FN_MAX_PARAMS];   /* FN_LAMBDA: interned parameter names */
    const Expr* body;                   /* FN_LAMBDA / FN_SLOTS */
    const char* head;                   /* FN_HEAD: interned head name */
    /* The `f` node as written.  FN_COMPOSE walks it; every kind uses it to
     * blame a bail on the user's own tree rather than on scaffolding. */
    const Expr* fexpr;
} FnSpec;

/* Resolve `f` as a function of exactly `want_arity` arguments.  Pure: it
 * inspects the tree and never emits, so both passes can call it. */
static bool fn_resolve(const Expr* f, int want_arity, FnSpec* out) {
    if (!f || want_arity < 1 || want_arity > FN_MAX_PARAMS) return false;
    memset(out, 0, sizeof *out);
    out->fexpr = f;

    if (f->type == EXPR_SYMBOL) {
        const char* nm = f->data.symbol.name;
        if (nm == SYM_Identity) {
            if (want_arity != 1) return false;
            out->kind = FN_IDENTITY; out->nparams = 1; return true;
        }
        /* A bare head.  Whether it can actually be lowered at these argument
         * types is emit_node's business — try_kernel and the CompiledFunction
         * callee path both live there — so resolution stays a syntactic test and
         * the single answer to "is this compilable" stays in one place. */
        out->kind = FN_HEAD; out->nparams = -1; out->head = nm; return true;
    }

    if (f->type != EXPR_FUNCTION || f->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = f->data.function.head->data.symbol.name;
    Expr* const* A = f->data.function.args;
    size_t na = f->data.function.arg_count;

    if (h == SYM_Composition) {
        /* Composition[] and Composition[f] never survive evaluation (they
         * normalise to Identity and f — src/core.c:1956-1957), so the general
         * form is the only one that reaches a body. */
        if (na < 1) return false;
        out->kind = FN_COMPOSE; out->nparams = want_arity; return true;
    }
    if (h != SYM_Function) return false;

    if (na == 1) { out->kind = FN_SLOTS; out->nparams = -1; out->body = A[0]; return true; }
    /* Function[params, body, attrs]: the attribute form can hold its arguments
     * (pure_function_attributes, src/purefunc.c:30), i.e. it changes evaluation
     * ORDER — which an inlined body does not reproduce. */
    if (na != 2) return false;

    const Expr* p = A[0];
    if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_Null) {   /* long-hand slot form */
        out->kind = FN_SLOTS; out->nparams = -1; out->body = A[1]; return true;
    }

    out->kind = FN_LAMBDA; out->body = A[1];
    if (p->type == EXPR_SYMBOL) {
        out->pname[0] = p->data.symbol.name; out->nparams = 1;
    } else if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL
               && p->data.function.head->data.symbol.name == SYM_List
               && p->data.function.arg_count >= 1
               && p->data.function.arg_count <= FN_MAX_PARAMS) {
        size_t k = p->data.function.arg_count;
        for (size_t i = 0; i < k; i++) {
            if (p->data.function.args[i]->type != EXPR_SYMBOL) return false;
            out->pname[i] = p->data.function.args[i]->data.symbol.name;
        }
        /* Duplicate parameter names would make the later binding shadow the
         * earlier one, where the interpreter substitutes both. */
        for (size_t i = 0; i < k; i++)
            for (size_t j = 0; j < i; j++)
                if (out->pname[i] == out->pname[j]) return false;
        out->nparams = (int)k;
    } else return false;

    /* Arity is EXACT.  A short call leaves the surplus parameter symbolic
     * (apply_pure_function, src/purefunc.c:271), so the interpreter's answer is
     * not a machine number and the compiled path must not produce one. */
    return out->nparams == want_arity;
}

/* `Slot[k]` (`#`, `#k`) resolved against the live slot frame, or -1.  Shared by
 * both passes so they cannot disagree about which slots are bound. */
static int fn_slot_index(const Ctx* c, Expr* const* A, size_t na) {
    long long k = 1;                                   /* bare `#` is Slot[1] */
    if (na == 1) {
        if (A[0]->type != EXPR_INTEGER) return -1;
        k = (long long)A[0]->data.integer;
    } else if (na != 0) return -1;
    return (k >= 1 && k <= (long long)c->nslot) ? (int)(k - 1) : -1;
}

/* Reserved parameter names for the FN_HEAD path (see emit_apply).  The context
 * mark makes them unproducible from user source, so binding them cannot shadow
 * anything a body could refer to. */
static const char* fn_placeholder(int i) {
    static const char* cache[FN_MAX_PARAMS];
    static const char* const names[FN_MAX_PARAMS] = {
        "System`Compile$fn1", "System`Compile$fn2", "System`Compile$fn3", "System`Compile$fn4",
        "System`Compile$fn5", "System`Compile$fn6", "System`Compile$fn7", "System`Compile$fn8"
    };
    if (!cache[i]) cache[i] = intern_symbol(names[i]);
    return cache[i];
}

/* Build `head[$1, ..., $n]` over the reserved placeholder symbols.  Caller frees. */
static Expr* fn_head_call(const char* head, int n) {
    Expr* args[FN_MAX_PARAMS];
    for (int i = 0; i < n; i++) {
        args[i] = expr_new_symbol(fn_placeholder(i));
        if (!args[i]) { while (i-- > 0) expr_free(args[i]); return NULL; }
    }
    Expr* hd = expr_new_symbol(head);
    if (!hd) { for (int i = 0; i < n; i++) expr_free(args[i]); return NULL; }
    return expr_new_function(hd, args, (size_t)n);   /* adopts hd and args */
}

/* ---- delegated structural heads -------------------------------------------
 *
 * These already have an NDArray entry point in the interpreter, and each of
 * those takes the WHOLE call and promises a result identical to the equivalent
 * List call (src/ndstruct.h:14, src/ndreduce.h:20).  So the compiled form is
 * that same function, called from the VM — which makes the compiled subset of
 * these heads the interpreted one by construction, exactly as A_PART does for
 * the general Part specs.  The win is not the operation (it was already a fast
 * buffer walk) but that a body CONTAINING one no longer bails wholesale.
 *
 * `rank_rule`: 0 = same rank as the operand, 1 = rank 1 (Flatten), 2 = rank 2
 * in and out (Transpose).  Every entry preserves the element type. */
typedef struct {
    const char* head;
    Expr* (*fn)(Expr*);
    int nextra;        /* trailing INTEGER arguments, passed through as written */
    int rank_rule;
} NdFnSpec;

static const NdFnSpec ND_FNS[] = {
    { "Reverse",    ndstruct_reverse,   0, 0 },
    { "Sort",       ndstruct_sort,      0, 0 },   /* 1-arg only: a comparator is a
                                                   * function value the ND path
                                                   * cannot call back into */
    { "Accumulate", ndred_accumulate,   0, 0 },
    { "Flatten",    ndstruct_flatten,   0, 1 },
    { "Transpose",  ndstruct_transpose, 0, 2 },
    { "Take",       ndstruct_take,      1, 0 },
    { "Drop",       ndstruct_drop,      1, 0 },
};

static const NdFnSpec* nd_fn_lookup(const char* h, size_t na) {
    for (size_t i = 0; i < sizeof ND_FNS / sizeof ND_FNS[0]; i++)
        if (strcmp(h, ND_FNS[i].head) == 0 && na == 1u + (size_t)ND_FNS[i].nextra)
            return &ND_FNS[i];
    return NULL;
}

/* Result type of a delegated head over an operand of type `ta`, or CT_ERR. */
static CompileType nd_fn_result(const NdFnSpec* s, CompileType ta) {
    if (!CT_IS_ARRAY(ta)) return CT_ERR;
    int rank = CT_RANK(ta);
    if (s->rank_rule == 1) return CT_ARRAY(CT_ELEM(ta), 1);
    if (s->rank_rule == 2) return rank == 2 ? ta : CT_ERR;
    return ta;
}

/* ---- counted-loop iterator specs -------------------------------------------
 * Do/Sum/Product accept the same counted spellings the interpreter does
 * (src/iter.c, iter_spec_parse):
 *
 *     n   or   {n}         repeat n times, binding no iterator symbol
 *     {i, hi}              i = 1, 2, ..., hi
 *     {i, lo, hi}          i = lo, ..., hi
 *     {i, lo, hi, di}      i = lo, lo + di, ...   (di a nonzero integer literal)
 *
 * `var == NULL` marks the two count-only forms; `lo == NULL` means the implicit
 * lower bound 1, so the common case synthesises no nodes.  List iteration
 * ({i, {a,b,c}}) is not compiled: the bound then fails to infer as CT_INT at the
 * call site, so it bails to the interpreter like any other unsupported form.
 * A non-literal step also bails — the loop test's direction has to be known when
 * the comparison opcode is chosen. */
typedef struct {
    const char* var;    /* iterator symbol, or NULL for the count-only forms */
    const Expr* lo;     /* NULL => literal 1 */
    const Expr* hi;
    long long   di;     /* step; nonzero */
} LoopSpec;

static bool loop_spec_parse(const Expr* spec, LoopSpec* out) {
    out->var = NULL; out->lo = NULL; out->hi = NULL; out->di = 1;
    if (!spec) return false;
    if (spec->type != EXPR_FUNCTION || spec->data.function.head->type != EXPR_SYMBOL
        || strcmp(spec->data.function.head->data.symbol.name, "List") != 0) {
        out->hi = spec;                                     /* bare count: Do[body, n] */
        return true;
    }
    size_t len = spec->data.function.arg_count;
    Expr* const* a = spec->data.function.args;
    if (len == 1) { out->hi = a[0]; return true; }           /* {n} */
    if (len < 2 || len > 4 || a[0]->type != EXPR_SYMBOL) return false;
    out->var = a[0]->data.symbol.name;
    if (len == 2) { out->hi = a[1]; return true; }           /* {i, hi} */
    out->lo = a[1]; out->hi = a[2];
    if (len == 4) {                                          /* {i, lo, hi, di} */
        if (a[3]->type != EXPR_INTEGER || a[3]->data.integer == 0) return false;
        out->di = a[3]->data.integer;
    }
    return true;
}

static bool infer_type(Ctx* c, const Expr* e, CompileType* out);

/* Both bounds must be integer-typed: these are integer-counted loops. */
static bool loop_spec_int_bounds(Ctx* c, const LoopSpec* s) {
    CompileType t;
    if (s->lo && (!infer_type(c, s->lo, &t) || t != CT_INT)) return false;
    return infer_type(c, s->hi, &t) && t == CT_INT;
}

/* Range[hi] / Range[lo, hi] / Range[lo, hi, di] as a LoopSpec, so Range shares
 * the element-COUNT arithmetic with Table rather than deriving its own — the
 * count is the whole difficulty, and two implementations of it would eventually
 * disagree at an endpoint.
 *
 * `var` stays NULL: Range has no iterator symbol to bind, its body IS the
 * iterator.  Bounds must infer as CT_INT (the interpreter's real Range walks by
 * repeated addition, which a closed-form step does not reproduce bit for bit)
 * and the step must be a nonzero literal, so the loop test's direction is known
 * when the comparison opcode is chosen. */
static bool range_spec(Ctx* c, Expr* const* A, size_t na, LoopSpec* out) {
    memset(out, 0, sizeof *out);
    out->di = 1;
    if (na == 1) { out->lo = NULL; out->hi = A[0]; }
    else         { out->lo = A[0]; out->hi = A[1]; }
    if (na == 3) {
        if (A[2]->type != EXPR_INTEGER) return false;
        out->di = (long long)A[2]->data.integer;
        if (out->di == 0) return false;
    }
    return loop_spec_int_bounds(c, out);
}

/* Materialise a loop bound into the persistent register `dst`.  A NULL bound is
 * the implicit lower limit 1, loaded straight as a constant so the count-only
 * and {i,hi} forms cost nothing extra. */
static bool emit_loop_bound(Ctx* c, const Expr* b, int dst) {
    Slot z = { 0 };
    if (!b) { Slot one; one.i = 1; ins(c, OP_CONST, (uint32_t)dst, 0, 0, one); return c->ok; }
    Val v; if (!emit(c, b, &v)) return false;
    ins(c, OP_MOVE, (uint32_t)dst, (uint32_t)v.reg, 0, z);
    free_if_tmp(c, v);
    return c->ok;
}

/* Defined below infer_apply, which they all use; declared here because the
 * inference branches for Nest/Fold/FixedPoint/NestWhile need them. */
static bool infer_apply(Ctx* c, const FnSpec* s, const CompileType* argt, int n,
                        CompileType* out);
static int  nest_fixed_type(Ctx* c, const FnSpec* s, CompileType t0);
static int  accum_fixed_type(Ctx* c, const FnSpec* s, CompileType t0,
                             const CompileType* rest, int nrest);
static CompileType vec_elem_type(Ctx* c, const Expr* e);
static bool fp_opts(Expr* const* A, size_t na, size_t start,
                    const Expr** max_out, const Expr** same_out);

/* Pure type inference (no code emission) — needed to type an If's result
 * register before both branches are lowered.  Mirrors emit's result-type rules;
 * returns false for anything not compilable. */
/* Lift a scalar result element type back over an array operand.
 *
 * Several inference branches compute "the element type this head produces" and
 * then returned it directly, which silently turned an array-valued node into a
 * scalar one.  That was invisible while every array operation was delegated to
 * the ND layer — which decides its own result dtype — and became two distinct
 * bugs once fusion started trusting these types: the wrong output buffer, and
 * whole bodies (anything rooted at Sin/Exp/Log/ArcTan over an array) never
 * being recognised as fusable at all. */
static CompileType arr_like(CompileType operand, CompileType elem_result) {
    return CT_IS_ARRAY(operand) ? CT_ARRAY(elem_result, CT_RANK(operand)) : elem_result;
}
static CompileType elem_of(CompileType t) { return CT_IS_ARRAY(t) ? CT_ELEM(t) : t; }

/* Element access opcode for an array's element type.  One place rather than the
 * fourteen `elem == CT_COMPLEX ? _C : _R` conditionals this replaces — adding
 * the integer element type to those by hand is exactly the kind of edit that
 * gets 13 of 14 sites. */
static uint16_t a_load_op(CompileType elem) {
    return elem == CT_COMPLEX ? OP_A_LOAD_C
         : elem == CT_INT     ? OP_A_LOAD_I
                              : OP_A_LOAD_R;
}
static uint16_t a_store_op(CompileType elem) {
    return elem == CT_COMPLEX ? OP_A_STORE_C
         : elem == CT_INT     ? OP_A_STORE_I
                              : OP_A_STORE_R;
}

/* Heads that are integer-CLOSED: given integer arguments the interpreter returns
 * an Integer, so a compiled body must too.  Every one of them also has a
 * registered real kernel — which is the trap, because reaching the kernel first
 * silently answers `35.` where the interpreter says `35`, and the numeric parity
 * tests cannot see the difference (see docs: result-HEAD parity).
 *
 * ONE table, consulted by both infer_type and emit, so the type a body is
 * declared to have and the opcode it actually runs cannot drift apart.
 *
 * Heads that look like they belong here and do NOT:
 *   LegendreP   `LegendreP[2, 2]` is 11/2 — Rational for some integer arguments.
 *   Beta        Rational.       HarmonicNumber  Rational.
 *   Divide      Rational; `Divide[7,3]` is 7/3 and no machine type holds it. */
typedef struct { const char* name; size_t arity; uint16_t op; } IntClosed;
static const IntClosed INT_CLOSED[] = {
    { "Factorial",  1, OP_FACT_I  },
    { "Gamma",      1, OP_GAMMA_I },
    { "Fibonacci",  1, OP_FIB_I   },
    { "LucasL",     1, OP_LUCAS_I },
    { "Binomial",   2, OP_BINOM_I },
    { "Pochhammer", 2, OP_POCH_I  },
};
static const IntClosed* int_closed_head(const char* h, size_t na) {
    for (size_t i = 0; i < sizeof INT_CLOSED / sizeof INT_CLOSED[0]; i++)
        if (INT_CLOSED[i].arity == na && strcmp(INT_CLOSED[i].name, h) == 0)
            return &INT_CLOSED[i];
    return NULL;
}
/* Integer-ONLY heads: ones with no real counterpart at all.  Unlike INT_CLOSED
 * above there is no kernel behind them to fall through to — non-integer
 * arguments simply bail, exactly as they did before.
 *
 * `*pred` distinguishes the three that answer True/False (CT_BOOL) from the ones
 * that answer an Integer.  GCD and LCM are Flat and n-ary in the interpreter, so
 * any arity from 2 up is accepted and folded left. */
static bool int_only_head(const char* h, size_t na, bool* pred) {
    *pred = false;
    if (na >= 2 && (strcmp(h, "GCD") == 0 || strcmp(h, "LCM") == 0)) return true;
    if ((na == 1 || na == 2)
        && (strcmp(h, "IntegerLength") == 0 || strcmp(h, "IntegerExponent") == 0)) return true;
    if (na == 3 && strcmp(h, "PowerMod") == 0) return true;
    if ((na == 1 && (strcmp(h, "EvenQ") == 0 || strcmp(h, "OddQ") == 0))
        || (na == 2 && strcmp(h, "Divisible") == 0)) { *pred = true; return true; }
    return false;
}

/* Peek at the argument types without emitting anything, so a head that is
 * integer-closed only on integer arguments can decline to the ordinary real
 * lowering for every other case. */
static bool all_args_int(Ctx* c, Expr* const* A, size_t na) {
    for (size_t i = 0; i < na; i++) {
        CompileType t;
        if (!infer_type(c, A[i], &t) || t != CT_INT) return false;
    }
    return true;
}

static bool infer_type(Ctx* c, const Expr* e, CompileType* out) {
    if (!e) return false;
    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) { *out = lt; return true; }
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType st; if (scope_find(c, nm, &st, NULL) >= 0) { *out = st; return true; }
        int k = arg_find(c, nm);
        if (k >= 0) { *out = c->arg_types[k]; return true; }
        if (strcmp(nm, "True") == 0 || strcmp(nm, "False") == 0) { *out = CT_BOOL; return true; }
        if (strcmp(nm, "I") == 0) { *out = CT_COMPLEX; return true; }
        double cv; if (named_const(nm, &cv)) { *out = CT_REAL; return true; }
        Slot gi; CompileType gt;
        if (global_const(c, nm, &gi, &gt)) { *out = gt; return true; }
        return false;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args; size_t na = e->data.function.arg_count;
    /* Slot[k] inside an inlined Function[body].  An index past the bound frame
     * falls through and fails, which is right: the interpreter would leave that
     * Slot unsubstituted, so its answer is not a machine number either. */
    if (h == SYM_Slot) { int k = fn_slot_index(c, A, na); if (k >= 0) { *out = c->slot[k].type; return true; } }
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
    if ((strcmp(h, "Mod") == 0 || strcmp(h, "Quotient") == 0) && na == 2) {
        IT(0, ta); IT(1, tb);
        if (ta == CT_INT && tb == CT_INT) { *out = CT_INT; return true; }
        if (ta <= CT_REAL && tb <= CT_REAL && ta != CT_BOOL && tb != CT_BOOL) {
            /* Real Mod/Quotient: the interpreter evaluates both (Mod[2.5,1.2] is
             * 0.1, Quotient[7.5,2.] is 3), so declining here would drop a whole
             * body to the interpreter over an operation we can do. */
            *out = (h[0] == 'M') ? CT_REAL : CT_INT;
            return true;
        }
        return false;
    }
    if (strcmp(h, "Power") == 0 && na == 2) {
        IT(0, ta);
        if (A[1]->type == EXPR_INTEGER) {
            CompileType el = elem_of(ta);
            CompileType res = (el == CT_INT && A[1]->data.integer >= 0)
                            ? CT_INT : (el == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            *out = arr_like(ta, res); return true;
        }
        int64_t rn, rd;
        if (is_rational(A[1], &rn, &rd) && rd == 2 && (rn == 1 || rn == -1)) {
            *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            return true;
        }
        IT(1, tb); ta = num_common(ta, tb); if ((int)ta < 0) return false;
        /* Integer base and integer exponent stay INTEGER, even though the
         * exponent's sign is only known per call: `2^n` is the Integer 1024 at
         * n = 10, and answering 1024. would be the wrong head.  OP_POW_II
         * abandons the call when the exponent turns out negative, because 2^-3
         * is the Rational 1/8 and no machine type holds it. */
        if (ta == CT_INT) { *out = CT_INT; return true; }
        if (ta < CT_REAL) ta = CT_REAL; *out = ta; return true;
    }
    uint16_t or_, oc_;
    if (na == 1 && (unary_math(h, &or_, &oc_) || strcmp(h, "Tanh") == 0)) {
        IT(0, ta);
        *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
        return true;
    }
    if (strcmp(h, "Log") == 0 && na == 2) {
        IT(0, ta); IT(1, tb); ta = num_common(ta, tb);
        if ((int)ta < 0) return false;
        *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
        return true;
    }
    /* Projections: the result is real even for a complex operand, and over an
     * array that means a REAL-element array, not an array of the operand's own
     * element type.  Getting this wrong is invisible while each operation is
     * delegated to the ND layer (which decides the result dtype itself), and
     * becomes a wrong answer the moment a fused loop allocates the output buffer
     * from this type instead. */
    if (strcmp(h, "Abs") == 0 && na == 1) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) {
            CompileType el = CT_ELEM(ta) == CT_COMPLEX ? CT_REAL : CT_ELEM(ta);
            *out = CT_ARRAY(el, CT_RANK(ta)); return true;
        }
        *out = ta == CT_COMPLEX ? CT_REAL : ta; return true;
    }
    /* Integer-closed heads on integer arguments (see INT_CLOSED).  Must be tested
     * BEFORE the generic kernel dispatch below: each of these HAS a registered
     * real kernel, and reaching it first is exactly how `Binomial[7,3]` came back
     * as 35. instead of 35. */
    if (int_closed_head(h, na) && all_args_int(c, A, na)) { *out = CT_INT; return true; }
    /* Integer-ONLY heads: no real counterpart, so they exist over CT_INT or not
     * at all.  The three predicates answer True/False, hence CT_BOOL. */
    {
        bool pred = false;
        if (int_only_head(h, na, &pred) && all_args_int(c, A, na)) {
            *out = pred ? CT_BOOL : CT_INT; return true;
        }
    }
    /* Sign of a real is the INTEGER -1, 0 or 1 — `Sign[-2.5]` is `-1`, not
     * `-1.` — so the result type is CT_INT, not the argument's type.  Sign of a
     * complex is z/|z|, which is genuinely complex. */
    if (strcmp(h, "Sign") == 0 && na == 1) {
        IT(0, ta);
        if (CT_IS_ARRAY(ta)) { *out = ta; return true; }
        *out = (ta == CT_COMPLEX) ? CT_COMPLEX : CT_INT;
        return true;
    }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0
         || strcmp(h, "IntegerPart") == 0) && na == 1) {
        IT(0, ta); if (CT_IS_ARRAY(ta)) return false; *out = CT_INT; return true; }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0) && na == 1) {
        IT(0, ta);
        if (ta == CT_BOOL) return false;
        /* Over the integers all three stay integral: Re is the value, Im is 0,
         * and Arg is 0 (or Pi, which ARG_I hands back to the interpreter). */
        if (ta == CT_INT) { *out = CT_INT; return true; }
        *out = CT_IS_ARRAY(ta) ? CT_ARRAY(CT_REAL, CT_RANK(ta)) : CT_REAL;
        return true;
    }
    /* FractionalPart of an integer is the Integer 0, not 0. */
    if (strcmp(h, "FractionalPart") == 0 && na == 1) {
        IT(0, ta);
        if (ta == CT_INT) { *out = CT_INT; return true; }
        if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        *out = ta; return true;
    }
    if (strcmp(h, "Conjugate") == 0 && na == 1) { IT(0, ta); *out = ta; return true; }
    if (strcmp(h, "HypergeometricPFQ") == 0 && na == 3) {
        const Expr* la = A[0];
        const Expr* lb = A[1];
        if (la->type != EXPR_FUNCTION || la->data.function.head->type != EXPR_SYMBOL
            || strcmp(la->data.function.head->data.symbol.name, "List") != 0
            || lb->type != EXPR_FUNCTION || lb->data.function.head->type != EXPR_SYMBOL
            || strcmp(lb->data.function.head->data.symbol.name, "List") != 0) return false;
        size_t np = la->data.function.arg_count, nq = lb->data.function.arg_count;
        if (np + nq + 2 > 8) return false;
        for (size_t i = 0; i < np; i++) {
            CompileType t;
            if (!infer_type(c, la->data.function.args[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
        }
        for (size_t j = 0; j < nq; j++) {
            CompileType t;
            if (!infer_type(c, lb->data.function.args[j], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
        }
        IT(2, ta); if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        *out = CT_REAL; return true;
    }
    if (strcmp(h, "UnitStep") == 0 && na >= 1) {
        for (size_t i = 0; i < na; i++) { IT(i, ta); if (CT_IS_ARRAY(ta) || ta == CT_BOOL) return false; }
        *out = CT_INT; return true;                    /* UnitStep[0.5] is 1, not 1. */
    }
    if ((strcmp(h, "Clip") == 0 || strcmp(h, "Rescale") == 0) && na == 2) {
        const Expr* bnd = A[1];
        if (bnd->type != EXPR_FUNCTION || bnd->data.function.head->type != EXPR_SYMBOL
            || strcmp(bnd->data.function.head->data.symbol.name, "List") != 0
            || bnd->data.function.arg_count != 2) return false;
        IT(0, ta);
        CompileType tl, tu;
        if (!infer_type(c, bnd->data.function.args[0], &tl)
            || !infer_type(c, bnd->data.function.args[1], &tu)) return false;
        if (tl != CT_REAL || tu != CT_REAL || CT_IS_ARRAY(ta) || ta == CT_BOOL) return false;
        /* Rescale carries a complex argument through; Clip cannot (Min/Max need
         * an order, and the interpreter leaves complex Clip unevaluated). */
        if (ta == CT_COMPLEX) {
            if (h[0] == 'C') return false;
            *out = CT_COMPLEX; return true;
        }
        *out = CT_REAL; return true;
    }
    if ((strcmp(h, "Max") == 0 || strcmp(h, "Min") == 0) && na >= 1) { IT(0, ta); for (size_t i = 1; i < na; i++) { IT(i, tb); ta = num_common(ta, tb); if ((int)ta < 0 || ta == CT_COMPLEX) return false; } *out = ta; return true; }
    if (strcmp(h, "ArcTan") == 0) {
        if (na == 1) {
            IT(0, ta);
            *out = arr_like(ta, elem_of(ta) == CT_COMPLEX ? CT_COMPLEX : CT_REAL);
            return true;
        }
        if (na == 2) {
            IT(0, ta); IT(1, tb); ta = num_common(ta, tb);
            if ((int)ta < 0) return false;
            *out = arr_like(ta, CT_REAL);
            return true;
        }
        return false;
    }
    if (na == 2 && (!strcmp(h, "Less") || !strcmp(h, "LessEqual") || !strcmp(h, "Greater") || !strcmp(h, "GreaterEqual") || !strcmp(h, "Equal") || !strcmp(h, "Unequal"))) { *out = CT_BOOL; return true; }
    if (!strcmp(h, "And") || !strcmp(h, "Or") || !strcmp(h, "Xor") || !strcmp(h, "Not")) { *out = CT_BOOL; return true; }
    if (strcmp(h, "If") == 0 && na == 3) {
        CompileType tt, te; if (!infer_type(c, A[1], &tt) || !infer_type(c, A[2], &te)) return false;
        if (tt == te) { *out = tt; return true; }
        *out = num_common(tt, te); return (int)*out >= 0;
    }
    if ((strcmp(h, "Sum") == 0 || strcmp(h, "Product") == 0) && na == 2) {
        LoopSpec s;
        /* s.var != NULL: the interpreter's Sum/Product reject a bare count
         * ({n}), so the compiled path must reject it too — only Do accepts it. */
        if (!loop_spec_parse(A[1], &s) || !s.var || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) return false;
        c->scope[c->nscope].name = s.var;
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
            || (int)(c->nscope + L->data.function.arg_count) > CTX_MAX_SCOPE) return false;
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
    if ((!strcmp(h, "Set") || !strcmp(h, "AddTo") || !strcmp(h, "SubtractFrom")
         || !strcmp(h, "TimesBy") || !strcmp(h, "DivideBy")) && na == 2
        && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt, NULL) < 0) return false; *out = vt; return true;
    }
    if ((!strcmp(h, "Increment") || !strcmp(h, "Decrement")) && na == 1 && A[0]->type == EXPR_SYMBOL) {
        CompileType vt; if (scope_find(c, A[0]->data.symbol.name, &vt, NULL) < 0) return false; *out = vt; return true;
    }
    if ((!strcmp(h, "Do") && na >= 2) || (!strcmp(h, "While") && na == 2) || (!strcmp(h, "For") && na == 4)
        || (!strcmp(h, "If") && na == 2)) { *out = CT_INT; return true; }
    if (!strcmp(h, "Nest") && na == 3) {
        FnSpec s; CompileType tn, tx;
        if (!fn_resolve(A[0], 1, &s) || !infer_type(c, A[2], &tn) || tn != CT_INT
            || !infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0) return false;
        *out = (CompileType)tfp; return true;
    }
    /* Delegated structural heads — the rank rule lives with the table. */
    {
        const NdFnSpec* nf = nd_fn_lookup(h, na);
        if (nf) {
            IT(0, ta);
            for (int i = 0; i < nf->nextra; i++) {
                CompileType te;
                if (!infer_type(c, A[1 + i], &te) || te != CT_INT) return false;
            }
            CompileType rt = nd_fn_result(nf, ta);
            if ((int)rt < 0) return false;
            *out = rt; return true;
        }
    }
    /* Select / TakeWhile / LengthWhile / All-, Any-, NoneTrue — one predicate
     * loop; see the emit-side block. */
    if (na == 2 && (!strcmp(h, "Select") || !strcmp(h, "TakeWhile")
                    || !strcmp(h, "LengthWhile") || !strcmp(h, "AllTrue")
                    || !strcmp(h, "AnyTrue") || !strcmp(h, "NoneTrue"))) {
        FnSpec s; if (!fn_resolve(A[1], 1, &s)) return false;
        CompileType el = vec_elem_type(c, A[0]);
        if ((int)el < 0) return false;
        CompileType tb;
        if (!infer_apply(c, &s, &el, 1, &tb) || tb != CT_BOOL) return false;
        if (!strcmp(h, "Select") || !strcmp(h, "TakeWhile")) *out = CT_ARRAY(el, 1);
        else if (!strcmp(h, "LengthWhile"))                  *out = CT_INT;
        else                                                 *out = CT_BOOL;
        return true;
    }
    /* First[v] / Last[v] are v[[1]] and v[[-1]]. */
    if ((!strcmp(h, "First") || !strcmp(h, "Last")) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta)) return false;
        *out = CT_RANK(ta) == 1 ? CT_ELEM(ta) : CT_ARRAY(CT_ELEM(ta), CT_RANK(ta) - 1);
        return true;
    }
    /* Map / Scan over a rank-1 array — see the emit-side block for why the
     * result element type must equal the source's. */
    if ((!strcmp(h, "Map") || !strcmp(h, "Scan")) && na == 2) {
        bool scan = h[0] == 'S';
        FnSpec s; if (!fn_resolve(A[0], 1, &s)) return false;
        CompileType el = vec_elem_type(c, A[1]);
        if ((int)el < 0) return false;
        CompileType rt;
        if (!infer_apply(c, &s, &el, 1, &rt) || CT_IS_ARRAY(rt)) return false;
        if (scan) { *out = CT_INT; return true; }
        if (rt != el) return false;
        *out = CT_ARRAY(rt, 1); return true;
    }
    /* Fold[f, x0, v] / Fold[f, v]: the accumulator's widening fixed point with
     * the element type held fixed. */
    if (!strcmp(h, "Fold") && (na == 2 || na == 3)) {
        FnSpec s; if (!fn_resolve(A[0], 2, &s)) return false;
        CompileType el = vec_elem_type(c, A[na - 1]);
        if ((int)el < 0) return false;
        CompileType t0 = el;                       /* Fold[f, v] seeds from v[[1]] */
        if (na == 3 && !infer_type(c, A[1], &t0)) return false;
        int tfp = accum_fixed_type(c, &s, t0, &el, 1);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        *out = (CompileType)tfp; return true;
    }
    /* NestList / FoldList — see the emit-side block for the element-type rule. */
    if ((!strcmp(h, "NestList") && na == 3)
        || (!strcmp(h, "FoldList") && (na == 2 || na == 3))) {
        bool nest = h[0] == 'N';
        FnSpec s; if (!fn_resolve(A[0], nest ? 1 : 2, &s)) return false;
        CompileType T;
        if (nest) {
            CompileType tn, tx;
            if (!infer_type(c, A[2], &tn) || tn != CT_INT
                || !infer_type(c, A[1], &tx)) return false;
            int tfp = nest_fixed_type(c, &s, tx);
            if (tfp < 0) return false;
            T = (CompileType)tfp;
        } else {
            CompileType el = vec_elem_type(c, A[na - 1]);
            if ((int)el < 0) return false;
            CompileType t0 = el;
            if (na == 3 && !infer_type(c, A[1], &t0)) return false;
            int tfp = accum_fixed_type(c, &s, t0, &el, 1);
            if (tfp < 0) return false;
            T = (CompileType)tfp;
        }
        if (!ct_is_elem(T)) return false;
        *out = CT_ARRAY(T, 1); return true;
    }
    /* FixedPointList / NestWhileList — a BUILT history, so Real/Complex only. */
    if ((!strcmp(h, "FixedPointList") && na >= 2 && na <= 4)
        || (!strcmp(h, "NestWhileList") && na == 3)) {
        bool fp = h[0] == 'F';
        FnSpec s, ss, ts; const Expr *mx = NULL, *st = NULL;
        if (!fn_resolve(A[0], 1, &s)) return false;
        if (fp) { if (!fp_opts(A, na, 2, &mx, &st)) return false; }
        else if (!fn_resolve(A[2], 1, &ts)) return false;
        CompileType tn;
        if (mx && (!infer_type(c, mx, &tn) || tn != CT_INT)) return false;
        CompileType tx; if (!infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0) return false;
        CompileType T = (CompileType)tfp;
        if (!ct_is_elem(T)) return false;
        CompileType tb;
        if (fp && st) {
            CompileType at[2] = { T, T };
            if (!fn_resolve(st, 2, &ss) || !infer_apply(c, &ss, at, 2, &tb) || tb != CT_BOOL)
                return false;
        }
        if (!fp && (!infer_apply(c, &ts, &T, 1, &tb) || tb != CT_BOOL)) return false;
        *out = CT_ARRAY(T, 1); return true;
    }
    if (!strcmp(h, "FixedPoint") && na >= 2 && na <= 4) {
        FnSpec s; const Expr *mx, *st;
        if (!fn_resolve(A[0], 1, &s) || !fp_opts(A, na, 2, &mx, &st)) return false;
        CompileType tn;
        if (mx && (!infer_type(c, mx, &tn) || tn != CT_INT)) return false;
        CompileType tx; if (!infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        if (st) {                                  /* SameTest must yield a Bool */
            FnSpec ss; CompileType tb, at[2] = { (CompileType)tfp, (CompileType)tfp };
            if (!fn_resolve(st, 2, &ss) || !infer_apply(c, &ss, at, 2, &tb) || tb != CT_BOOL)
                return false;
        }
        *out = (CompileType)tfp; return true;
    }
    if (!strcmp(h, "NestWhile") && na == 3) {
        FnSpec s, ts; CompileType tx, tb;
        if (!fn_resolve(A[0], 1, &s) || !fn_resolve(A[2], 1, &ts)
            || !infer_type(c, A[1], &tx)) return false;
        int tfp = nest_fixed_type(c, &s, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) return false;
        CompileType at = (CompileType)tfp;
        if (!infer_apply(c, &ts, &at, 1, &tb) || tb != CT_BOOL) return false;
        *out = at; return true;
    }
    /* Part (M3c).  A full run of scalar subscripts drops every axis and lands
     * back in the scalar lattice; anything else keeps at least one axis and
     * stays an array, with the rank counted the way build_axis_selector counts
     * it: a spec keeps its axis, a scalar subscript drops it, and axes past the
     * last subscript are implicit Alls and survive. */
    if (strcmp(h, "Part") == 0 && na >= 2) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta) || (int)(na - 1) > CT_RANK(ta)) return false;
        if (part_is_scalar_indexed(c, ta, (const Expr* const*)A + 1, na - 1)) {
            *out = CT_ELEM(ta); return true;
        }
        int keep = CT_RANK(ta) - (int)(na - 1);
        for (size_t i = 1; i < na; i++) if (subscript_is_literal_spec(A[i])) keep++;
        if (keep < 1 || keep > CT_MAX_RANK) return false;
        *out = CT_ARRAY(CT_ELEM(ta), keep);
        return true;
    }
    if (strcmp(h, "ConstantArray") == 0 && na == 2) {
        int rank; CompileType elem;
        if (!const_array_shape(c, (const Expr* const*)A, &rank, &elem)) return false;
        *out = CT_ARRAY(elem, rank); return true;
    }
    /* Range[hi] / Range[lo, hi] / Range[lo, hi, di] — see the emit-side block. */
    if (strcmp(h, "Range") == 0 && na >= 1 && na <= 3) {
        LoopSpec s;
        if (!range_spec(c, A, na, &s)) return false;
        *out = CT_ARRAY(CT_INT, 1); return true;
    }
    /* Table[body, spec...] — see the emit-side block for why the iterators must
     * be integer-bounded and the element type Real or Complex. */
    if (strcmp(h, "Table") == 0 && na >= 2 && (int)na - 1 <= CT_MAX_RANK) {
        int rank = (int)na - 1;
        LoopSpec sp[CT_MAX_RANK];
        for (int j = 0; j < rank; j++)
            if (!loop_spec_parse(A[j + 1], &sp[j]) || !loop_spec_int_bounds(c, &sp[j])) return false;
        if (c->nscope + rank > CTX_MAX_SCOPE) return false;
        int pushed = 0;
        for (int j = 0; j < rank; j++)
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++; pushed++;
            }
        CompileType et; bool okT = infer_type(c, A[0], &et);
        c->nscope -= pushed;
        if (!okT || !ct_is_elem(et)) return false;
        *out = CT_ARRAY(et, rank); return true;
    }

    /* array -> scalar reductions (M3a): the only way an array type re-enters
     * the scalar lattice, so these must be inferable inside If/Sum/... */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        IT(0, ta);
        if (!CT_IS_ARRAY(ta) || (h[0] == 'T' && CT_RANK(ta) != 1)) return false;
        *out = (h[0] == 'L') ? CT_INT : CT_ELEM(ta);
        return true;
    }
    if (na == 1) { SymbolDef* d = symtab_lookup(h); if (d && d->ndarray_unary_kernel) { const NDUnaryKernel* k = d->ndarray_unary_kernel; if (k->cplx || k->real) { IT(0, ta); if (CT_IS_ARRAY(ta)) { *out = CT_ARRAY(k->to_real ? CT_REAL : CT_ELEM(ta), CT_RANK(ta)); return true; } if (k->to_real) { *out = CT_REAL; return true; } if (ta == CT_COMPLEX) { if (!k->cplx) return false; *out = CT_COMPLEX; return true; } *out = (k->real_closed || k->real) ? CT_REAL : CT_COMPLEX; return true; } } }
    if (na >= 3 && na <= 8) {
        SymbolDef* d = symtab_lookup(h);
        if (d && d->ndarray_nary_kernel) {
            const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
            if (k->cplx && k->nargs == na) {
                for (size_t i = 0; i < na; i++) {
                    CompileType t;
                    if (!infer_type(c, A[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
                }
                *out = CT_REAL; return true;
            }
        }
    }
    if (na == 2) {
        SymbolDef* d = symtab_lookup(h);
        if (d && d->ndarray_binary_kernel) {
            const NDBinaryKernel* k = d->ndarray_binary_kernel;
            if (k->cplx) {
                IT(0, ta); IT(1, tb);
                CompileType t = num_common(ta, tb);
                if ((int)t < 0) return false;
                /* Over an array the result is an ARRAY, mirroring the unary
                 * branch above and the emit-side lowering.  Reporting a scalar
                 * here was harmless while every array op was delegated (the ND
                 * layer picks its own result dtype), and becomes a wrong output
                 * buffer as soon as a fused loop is sized from this type. */
                if (CT_IS_ARRAY(t)) {
                    CompileType er = CT_ELEM(t);
                    *out = CT_ARRAY(k->real_closed ? er : CT_COMPLEX, CT_RANK(t));
                    return true;
                }
                *out = (t <= CT_REAL && k->real_closed) ? CT_REAL : CT_COMPLEX;
                return true;
            }
        }
    }
    #undef IT
    /* Inlined CompiledFunction call — must type here too, or a call inside an
     * If branch / Sum body could not be lowered.  Mirrors emit's binding. */
    {
        const CompiledFunction* cf = compiled_callee(c, h);
        if (cf && na >= 1 && compiled_function_num_args(cf) == na
            && na + (size_t)c->nscope <= CTX_MAX_SCOPE && c->inlining < 8) {
            const char* const* pn = compiled_function_arg_names(cf);
            const CompileType* pt = compiled_function_arg_types(cf);
            /* Type every argument in the CALLER's environment first — binding a
             * parameter before the later arguments are typed would let it
             * shadow a caller variable of the same name. */
            for (size_t i = 0; i < na; i++) {
                CompileType ta;
                if (!infer_type(c, A[i], &ta) || CT_IS_ARRAY(ta) || CT_IS_ARRAY(pt[i]))
                    return false;
            }
            int saved_scope = c->nscope;
            for (size_t i = 0; i < na; i++) {
                c->scope[c->nscope].name = pn[i];
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = pt[i];
                c->nscope++;
            }
            c->inlining++;
            bool okb = infer_type(c, compiled_function_body(cf), out);
            c->inlining--;
            c->nscope = saved_scope;
            return okb && !CT_IS_ARRAY(*out);
        }
    }
    return false;
}

/* Type `f[a1, ..., an]` for an already-resolved function value.  Mirrors
 * emit_apply's binding exactly — the two must agree or a construct types as one
 * thing and lowers as another. */
static bool infer_apply(Ctx* c, const FnSpec* s, const CompileType* argt, int n,
                        CompileType* out) {
    switch (s->kind) {
        case FN_IDENTITY:
            if (n != 1) return false;
            *out = argt[0]; return true;

        case FN_LAMBDA: {
            if (n != s->nparams || c->nscope + n > CTX_MAX_SCOPE) return false;
            int saved_scope = c->nscope, saved_nslot = c->nslot;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = s->pname[i];
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = argt[i];
                c->nscope++;
            }
            c->nslot = 0;                    /* see emit_apply for why */
            bool ok = infer_type(c, s->body, out);
            c->nslot = saved_nslot; c->nscope = saved_scope;
            return ok;
        }

        case FN_SLOTS: {
            if (n > FN_MAX_PARAMS) return false;
            int saved_nslot = c->nslot;
            struct { int reg; CompileType type; } saved[FN_MAX_PARAMS];
            memcpy(saved, c->slot, sizeof saved);
            for (int i = 0; i < n; i++) { c->slot[i].reg = 0; c->slot[i].type = argt[i]; }
            c->nslot = n;
            bool ok = infer_type(c, s->body, out);
            c->nslot = saved_nslot;
            memcpy(c->slot, saved, sizeof saved);
            return ok;
        }

        case FN_COMPOSE: {
            /* Composition[f1,...,fk][a...] = f1[f2[...fk[a...]]] — the INNERMOST
             * takes every argument, each outer takes one (src/eval.c:1449). */
            Expr* const* F = s->fexpr->data.function.args;
            size_t nf = s->fexpr->data.function.arg_count;
            CompileType t; FnSpec inner;
            if (!fn_resolve(F[nf - 1], n, &inner) || !infer_apply(c, &inner, argt, n, &t))
                return false;
            for (size_t k = nf - 1; k > 0; k--) {
                FnSpec g;
                if (!fn_resolve(F[k - 1], 1, &g) || !infer_apply(c, &g, &t, 1, &t)) return false;
            }
            *out = t; return true;
        }

        case FN_HEAD: {
            if (n > FN_MAX_PARAMS || c->nscope + n > CTX_MAX_SCOPE) return false;
            Expr* call = fn_head_call(s->head, n);
            if (!call) return false;
            int saved_scope = c->nscope;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = fn_placeholder(i);
                c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = argt[i];
                c->nscope++;
            }
            bool ok = infer_type(c, call, out);
            c->nscope = saved_scope;
            expr_free(call);
            return ok;
        }
    }
    return false;
}

/* Fixed-point accumulator type for an iteration that feeds its own output back
 * in (Nest, Fold, FixedPoint, NestWhile): the register must hold a type wide
 * enough to absorb every iteration.  Starting from t0, widen until the output
 * type stops growing — the lattice is bounded, so this converges in a few
 * passes or not at all.
 *
 * Argument 0 is the accumulator; `rest` carries the types of any further
 * arguments, which do not vary (Fold's list element). */
static int accum_fixed_type(Ctx* c, const FnSpec* s, CompileType t0,
                            const CompileType* rest, int nrest) {
    if (nrest < 0 || nrest + 1 > FN_MAX_PARAMS) return -1;
    CompileType t = t0;
    for (int iter = 0; iter < 4; iter++) {
        CompileType at[FN_MAX_PARAMS], tb;
        at[0] = t;
        for (int i = 0; i < nrest; i++) at[i + 1] = rest[i];
        if (!infer_apply(c, s, at, nrest + 1, &tb) || (int)tb < 0) return -1;
        if (tb == CT_BOOL && t != CT_BOOL) return -1;    /* can't fold a bool into a number */
        if (tb <= t) return (int)t;                       /* output coerces down into the accumulator */
        t = tb;                                           /* output widened it; grow and re-check */
    }
    return -1;
}
static int nest_fixed_type(Ctx* c, const FnSpec* s, CompileType t0) {
    return accum_fixed_type(c, s, t0, NULL, 0);
}

/* The element type of a rank-1 array-valued expression, or CT_ERR.  Fold, Map
 * and the other list-consuming heads all need exactly this test. */
static CompileType vec_elem_type(Ctx* c, const Expr* e) {
    CompileType t;
    if (!infer_type(c, e, &t) || !CT_IS_ARRAY(t) || CT_RANK(t) != 1) return CT_ERR;
    return CT_ELEM(t);
}

/* FixedPoint's trailing arguments: an application bound and/or SameTest -> s,
 * in either order, exactly as parse_fp_opts accepts them (src/funcprog.c:2511).
 * Returns false for a spelling the interpreter itself refuses. */
static bool fp_opts(Expr* const* A, size_t na, size_t start,
                    const Expr** max_out, const Expr** same_out) {
    *max_out = NULL; *same_out = NULL;
    for (size_t i = start; i < na; i++) {
        const Expr* a = A[i];
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL
            && a->data.function.args[0]->data.symbol.name == SYM_SameTest) {
            if (*same_out) return false;
            *same_out = a->data.function.args[1];
        } else if (!*max_out) {
            *max_out = a;        /* the bound; must infer as CT_INT to compile */
        } else return false;
    }
    return true;
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
    if (na >= 3 && na <= 8 && !c->vector_mode) {
        /* N-ary machine kernel.  Arguments must land in CONSECUTIVE registers
         * because the instruction carries only their base and count, so the
         * block is reserved first and each argument lowered into it — the same
         * shape OP_CALL uses. */
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) return false;
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx || k->nargs != na) return false;
        for (size_t i = 0; i < na; i++) {
            CompileType t;
            if (!infer_type(c, A[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL) return false;
        }
        int base = alloc_temp(c);
        for (size_t i = 1; i < na; i++) (void)alloc_temp(c);
        int after = c->temp_top;
        Slot z; memset(&z, 0, sizeof z);
        for (size_t i = 0; i < na && c->ok; i++) {
            Val v;
            if (!emit(c, A[i], &v)) return false;
            coerce(c, &v, CT_REAL);
            if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)(base + (int)i), (uint32_t)v.reg, 0, z);
            c->temp_top = after;
        }
        if (!c->ok) return false;
        c->temp_top = (base - c->nlocals) + 1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        ins_f(c, OP_KERNN, (uint16_t)na, (uint32_t)base, (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = CT_REAL;
        return c->ok;
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

/* ------------------------------------------------------------------ *
 *  Elementwise fusion (M3b)                                           *
 * ------------------------------------------------------------------ *
 * The delegated array opcodes (V_EW, V_KERN, ...) call the same NDArray kernels
 * the interpreter does, so they only ever removed the per-operation evaluator
 * round-trip: `Total[Sin[v] Exp[-v] + Sqrt[v]]` measured 1.0x at length 4096,
 * because both paths make five full-length passes and allocate four temporary
 * buffers.  Array speed is not about removing interpretation; it is about
 * removing intermediate buffers.
 *
 * So an elementwise chain is lowered instead to ONE counted loop over the flat
 * index, whose body is ordinary scalar bytecode reading elements through
 * A_LOAD and writing through A_STORE.  Three things fall out of that choice:
 *
 *   - it is rank-agnostic (a flat index never looks at the shape),
 *   - a trailing `Total` folds into an accumulator register in the same pass,
 *   - the scalar optimiser applies INSIDE the loop, so a loop-invariant
 *     subexpression is hoisted out of the element loop for free.
 *
 * Whether a subtree is fusable is not analysed up front: the body is lowered
 * speculatively with the array leaves bound to element registers, and if the
 * ordinary scalar emitter cannot do it, the emitted code is rolled back and the
 * caller falls through to the delegated path.  That way the fusable set is
 * exactly the scalar-compilable set, by construction, and cannot drift from it.
 */

/* Save/restore point for speculative lowering. */
typedef struct {
    size_t n; int temp_top, maxreg, arr_top, arr_max, nscope; bool ok;
    int tile_top, tile_max; bool vector_mode;
    const Expr* bail_node;
} EmitMark;

static EmitMark emit_mark(const Ctx* c) {
    EmitMark m;
    m.n = c->n; m.temp_top = c->temp_top; m.maxreg = c->maxreg;
    m.arr_top = c->arr_top; m.arr_max = c->arr_max; m.nscope = c->nscope; m.ok = c->ok;
    m.tile_top = c->tile_top; m.tile_max = c->tile_max; m.vector_mode = c->vector_mode;
    /* A speculative lowering that fails is not a bail — the caller falls through
     * to another strategy — so its diagnostic must not survive the rollback. */
    m.bail_node = c->bail_node;
    return m;
}
static void emit_rollback(Ctx* c, EmitMark m) {
    c->n = m.n; c->temp_top = m.temp_top; c->maxreg = m.maxreg;
    c->arr_top = m.arr_top; c->arr_max = m.arr_max; c->nscope = m.nscope; c->ok = m.ok;
    c->tile_top = m.tile_top; c->tile_max = m.tile_max; c->vector_mode = m.vector_mode;
    c->bail_node = m.bail_node;
}

#define FUSE_MAX_LEAVES 8

typedef struct {
    const char* name[FUSE_MAX_LEAVES];   /* interned symbol pointer */
    int         reg[FUSE_MAX_LEAVES];    /* the array-valued register */
    CompileType type[FUSE_MAX_LEAVES];   /* the ARRAY type (element + rank) */
    bool        built[FUSE_MAX_LEAVES];  /* leaf constructed here, not an argument */
    int         n;
} FuseLeaves;

/* Collect the distinct array-valued symbols in `e`.  Returns false if the tree
 * contains more than FUSE_MAX_LEAVES of them. */
static bool fuse_collect(Ctx* c, const Expr* e, FuseLeaves* L) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType t; bool tb = false; int reg = scope_find(c, nm, &t, &tb);
        if (reg < 0) {
            int k = arg_find(c, nm);
            if (k < 0) return true;                 /* not a value we bind */
            /* A fused leaf is a READ of the argument like any other.  Recording
             * it here as well as in emit() matters because a fully fused body
             * never reaches emit()'s symbol case: without this the argument
             * looks unread, and compiled_arg_deps would tell a client (an FD
             * Jacobian, say) the function does not depend on it. */
            c->argdep[k] = 1;
            reg = k; t = c->arg_types[k]; tb = false;   /* borrowed from the caller */
        }
        if (!CT_IS_ARRAY(t)) return true;
        for (int i = 0; i < L->n; i++) if (L->name[i] == nm) return true;   /* seen */
        if (L->n >= FUSE_MAX_LEAVES) return false;
        L->name[L->n] = nm; L->reg[L->n] = reg; L->type[L->n] = t; L->built[L->n] = tb; L->n++;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return true;
    if (!fuse_collect(c, e->data.function.head, L)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!fuse_collect(c, e->data.function.args[i], L)) return false;
    return true;
}

/* Is it legal to thread this subtree element by element?
 *
 * ONLY if the interpreter would do the same, and the interpreter threads over a
 * list exactly when the head is Listable.  `Max[v, w]` is the trap: it is not
 * Listable, so the interpreter returns the single largest element, and an
 * elementwise fusion would quietly answer something else.  `If`, `With` and
 * `Sum` are the same story.  The compiled path must never answer where the
 * interpreter declines — nor differently from it — even when the loop it would
 * generate is perfectly well defined.
 *
 * Returns 1 if the subtree reaches an array leaf legally, 0 if it contains no
 * array leaf at all (a scalar subtree, which may use any head), -1 to reject. */
static int fuse_listable(Ctx* c, const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType t; int reg = scope_find(c, nm, &t, NULL);
        if (reg < 0) { int k = arg_find(c, nm); if (k < 0) return 0; t = c->arg_types[k]; }
        return CT_IS_ARRAY(t) ? 1 : 0;
    }
    if (e->type != EXPR_FUNCTION) return 0;
    int any = 0, narr = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int s = fuse_listable(c, e->data.function.args[i]);
        if (s < 0) return -1;
        if (s > 0) { any = 1; narr++; }
    }
    if (!any) return 0;                       /* scalar subtree: unconstrained */
    if (e->data.function.head->type != EXPR_SYMBOL) return -1;
    const char* h = e->data.function.head->data.symbol.name;
    if (!(get_attributes(h) & ATTR_LISTABLE)) return -1;

    /* Being Listable is necessary but not sufficient when BOTH operands are
     * arrays.  The NDArray layer's binary kernels are defined as "one array plus
     * one broadcast scalar", so for a head that has one the interpreter has no
     * array-by-array path at all and leaves the expression UNEVALUATED —
     * verified: `ArcTan[NDArray[...], NDArray[...]]` comes back untouched.  A
     * fused loop would happily compute it, which is precisely the divergence the
     * engine is not allowed to introduce.
     *
     * `Log[b, x]` is the exception, and not an arbitrary one: it never reaches a
     * binary kernel because its lowering is the arithmetic identity
     * Log[x]/Log[b], so the interpreter evaluates it array-by-array too (it
     * returns an NDArray). */
    if (narr >= 2 && e->data.function.arg_count == 2 && strcmp(h, "Log") != 0) {
        SymbolDef* d = symtab_lookup(h);
        if (d && d->ndarray_binary_kernel) return -1;
    }
    return 1;
}

static bool emit(Ctx* c, const Expr* e, Val* out);

static bool try_fuse(Ctx* c, const Expr* e, Val* out) {
    if (c->fusing) return false;                   /* already inside a fused loop */

    /* `Total[chain]` folds the reduction into the same pass.  Only for rank 1:
     * for higher rank Total reduces the LEADING axis, not every element, so it
     * is a different operation and must keep the delegated path. */
    const Expr* body = e;
    bool reduce = false;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "Total") == 0
        && e->data.function.arg_count == 1) {
        CompileType it;
        if (infer_type(c, e->data.function.args[0], &it)
            && CT_IS_ARRAY(it) && CT_RANK(it) == 1) {
            body = e->data.function.args[0];
            reduce = true;
        }
    }

    CompileType bt;
    if (!infer_type(c, body, &bt) || !CT_IS_ARRAY(bt)) return false;
    CompileType elem = CT_ELEM(bt);
    int rank = CT_RANK(bt);
    /* Integer element types are excluded from FUSION specifically: the tile
     * opcodes (VADD_R, VMUL_C, ...) are real and complex only, so an integer
     * chain takes the ordinary per-element loop instead.  Correct, just not
     * strip-mined. */
    if (elem != CT_REAL && elem != CT_COMPLEX) return false;

    if (fuse_listable(c, body) != 1) return false;

    FuseLeaves L; L.n = 0;
    if (!fuse_collect(c, body, &L) || L.n == 0) return false;
    for (int i = 0; i < L.n; i++)
        if (CT_RANK(L.type[i]) != rank) return false;          /* no broadcasting */
    if (c->nscope + L.n > (int)(sizeof c->scope / sizeof c->scope[0])) return false;

    EmitMark mark = emit_mark(c);

    /* Persistent loop registers FIRST, then their inits into temps above them:
     * allocating a persistent register after emitting a temp init lets
     * free_if_tmp drop temp_top back onto it, and the next alloc clobbers it. */
    int rn   = alloc_temp(c);                       /* element count */
    int ri   = alloc_temp(c);                       /* flat index */
    int racc = reduce ? alloc_temp(c) : -1;         /* folded Total accumulator */
    int rout = reduce ? -1 : alloc_arr(c);          /* result buffer */

    Slot z; memset(&z, 0, sizeof z);
    ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)L.reg[0], 0, z);
    for (int i = 1; i < L.n; i++)
        ins(c, OP_A_SHAPECHK, 0, (uint32_t)L.reg[0], (uint32_t)L.reg[i], z);
    if (!reduce)
        ins_f(c, OP_A_NEWLIKE, (uint16_t)((unsigned)elem << AF_R_SHIFT),
              (uint32_t)rout, (uint32_t)L.reg[0], 0, z);
    else {
        Slot k0; memset(&k0, 0, sizeof k0);
        ins(c, OP_CONST, (uint32_t)racc, 0, 0, k0);   /* 0.0 / 0.0+0.0i */
    }
    Slot k0; memset(&k0, 0, sizeof k0);
    ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);

    /* Guard once on entry, then close the loop with a single LOOP instruction:
     * per tile that is one control instruction instead of four. */
    int rc = alloc_temp(c);
    ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
    size_t jz = c->n;
    ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
    c->temp_top--;                                   /* rc is dead after the test */

    /* Fan-out marker for the strip loop that follows.  MAP ONLY.
     *
     * A map writes each output element from its own input element, so splitting
     * the index range across threads is bit-identical to running it serially —
     * no element's value depends on which chunk computed it, and the chunks
     * write disjoint memory.  A REDUCTION is not: floating-point addition is not
     * associative, so per-thread partial sums would give a different (not
     * merely differently-rounded) answer from the serial fold, and the compiled
     * path is required to agree with the interpreter, not merely to be close to
     * it.  Reductions therefore stay serial until there is a reason to match
     * some specific chunking, which would be a promise about thread count.
     *
     * At runtime this either fans out and jumps past the loop, or falls straight
     * through into exactly the serial loop that was always emitted here. */
    size_t apar = (size_t)-1;
    if (!reduce && !(c->flags & COMPILE_NO_PAR)) {
        apar = c->n;
        ins(c, OP_APAR, (uint32_t)ri, (uint32_t)rn, 0, z);
    }
    size_t head = c->n;

    /* Live length of this tile (short only on the final one). */
    ins(c, OP_VSETLEN, 0, (uint32_t)ri, (uint32_t)rn, z);

    /* Bind each array leaf to the TILE holding its current block of elements,
     * then lower the body with the ORDINARY scalar emitter.  The binding goes on
     * the scope stack, which emit() and infer_type() already consult ahead of the
     * argument map, and the tile-ness rides on the register number — so nothing
     * in the scalar lowering needs to know about strip mining at all. */
    c->vector_mode = true;
    int tlreg[FUSE_MAX_LEAVES];
    for (int i = 0; i < L.n; i++) {
        CompileType le = CT_ELEM(L.type[i]);
        tlreg[i] = alloc_tile(c);
        ins(c, le == CT_COMPLEX ? OP_VLOAD_C : OP_VLOAD_R,
            (uint32_t)tlreg[i], (uint32_t)L.reg[i], (uint32_t)ri, z);
        c->scope[c->nscope].name = L.name[i];
        c->scope[c->nscope].reg  = tlreg[i];
        c->scope[c->nscope].type = le;
        c->nscope++;
    }

    c->fusing++;
    Val bv;
    bool bok = emit(c, body, &bv);
    c->fusing--;
    c->nscope -= L.n;

    if (!bok || !c->ok || CT_IS_ARRAY(bv.type)) { emit_rollback(c, mark); return false; }
    coerce(c, &bv, elem);
    /* A body that never touched a leaf (all-scalar) yields a scalar; broadcast it
     * so the store below always has a tile to write. */
    if (c->ok && !val_is_tile(bv)) bv = vsplat(c, bv);
    if (!c->ok) { emit_rollback(c, mark); return false; }

    if (reduce)
        ins(c, elem == CT_COMPLEX ? OP_VSUM_C : OP_VSUM_R,
            (uint32_t)racc, (uint32_t)bv.reg, 0, z);
    else
        ins(c, elem == CT_COMPLEX ? OP_VSTORE_C : OP_VSTORE_R,
            (uint32_t)rout, (uint32_t)ri, (uint32_t)bv.reg, z);

    /* Drop every body tile AND the leaf tiles: the next iteration rewrites them,
     * so nothing may accumulate across the back edge. */
    c->temp_top = (rc - c->nlocals);
    c->tile_top = mark.tile_top;
    c->vector_mode = mark.vector_mode;

    Slot step; memset(&step, 0, sizeof step); step.i = VBLOCK;
    ins(c, OP_LOOP, (uint32_t)ri, (uint32_t)rn, (uint32_t)head, step);
    if (!c->ok) { emit_rollback(c, mark); return false; }
    c->code[jz].b = (uint32_t)c->n;
    if (apar != (size_t)-1) c->code[apar].b = (uint32_t)c->n;

    if (reduce) {
        c->temp_top = (racc - c->nlocals) + 1;       /* keep racc, drop rn/ri */
        out->reg = racc; out->tmp = true; out->type = elem;
    } else {
        c->temp_top = (rn - c->nlocals);             /* rn/ri dead; rout is an arr */
        /* The fused buffer is built only if every leaf it reads was. */
        bool blt = true;
        for (int i = 0; i < L.n; i++) blt = blt && L.built[i];
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(elem, rank); out->built = blt;
    }
    return true;
}

/* Lower `f[a1, ..., an]` for an already-resolved function value, with the
 * arguments ALREADY in registers.
 *
 * The caller owns those registers and must keep them live for the whole body —
 * a parameter is read wherever the body mentions it, while the body is free to
 * allocate temps above them.  In practice every caller binds a persistent
 * register (an accumulator, a loop element), so no copy is made here.
 *
 * `out` may share an argument's register when the body is a bare parameter
 * reference, exactly as an ordinary lowering may return an operand's register;
 * `out->tmp` says whether the caller must free it. */
static bool emit_apply(Ctx* c, const FnSpec* s, const Val* argv, int n, Val* out) {
    if (!c->ok) return false;

    switch (s->kind) {
        case FN_IDENTITY:
            if (n != 1) { c->ok = false; return false; }
            *out = argv[0]; out->tmp = false; return true;

        case FN_LAMBDA: {
            if (n != s->nparams || c->nscope + n > CTX_MAX_SCOPE) { c->ok = false; return false; }
            int saved_scope = c->nscope, saved_nslot = c->nslot;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = s->pname[i];
                c->scope[c->nscope].reg = argv[i].reg;
                c->scope[c->nscope].type = argv[i].type;
                c->scope[c->nscope].built = argv[i].built;
                c->nscope++;
            }
            /* Slots are INVISIBLE under a named lambda.  apply_pure_function's
             * named path runs substitute_names only (src/purefunc.c:247), so a
             * Slot[1] inside Function[u, # + u] survives into the interpreter's
             * answer and that answer is not a machine number.  Hiding the frame
             * is what keeps the compiled path from answering where it declines. */
            c->nslot = 0;
            bool ok = emit(c, s->body, out);
            c->nslot = saved_nslot; c->nscope = saved_scope;
            return ok;
        }

        case FN_SLOTS: {
            if (n > FN_MAX_PARAMS) { c->ok = false; return false; }
            int saved_nslot = c->nslot;
            struct { int reg; CompileType type; } saved[FN_MAX_PARAMS];
            memcpy(saved, c->slot, sizeof saved);
            for (int i = 0; i < n; i++) { c->slot[i].reg = argv[i].reg; c->slot[i].type = argv[i].type; }
            c->nslot = n;
            bool ok = emit(c, s->body, out);
            c->nslot = saved_nslot;
            memcpy(c->slot, saved, sizeof saved);
            return ok;
        }

        case FN_COMPOSE: {
            /* Composition[f1,...,fk][a...] = f1[f2[...fk[a...]]] — the INNERMOST
             * takes every argument, each outer takes one (src/eval.c:1449). */
            Expr* const* F = s->fexpr->data.function.args;
            size_t nf = s->fexpr->data.function.arg_count;
            FnSpec inner; Val v;
            if (!fn_resolve(F[nf - 1], n, &inner)) { c->ok = false; return false; }
            if (!emit_apply(c, &inner, argv, n, &v)) return false;
            for (size_t k = nf - 1; k > 0; k--) {
                FnSpec g;
                if (!fn_resolve(F[k - 1], 1, &g)) { c->ok = false; return false; }
                Val r;
                if (!emit_apply(c, &g, &v, 1, &r)) return false;
                free_if_tmp(c, v);
                v = r;
            }
            *out = v; return c->ok;
        }

        case FN_HEAD: {
            /* emit_node dispatches on a head name plus an Expr** argument list,
             * not on Vals, so applying a bare head means synthesizing the call
             * over reserved placeholder symbols bound to the argument registers.
             * The alternative — a Val-taking entry point into ~45 lowerings —
             * would be a large refactor with no behaviour change, and this way
             * EVERY head emit_node knows (all 93 machine kernels, plus a
             * CompiledFunction-valued symbol) becomes usable as a function value
             * with no per-head work.  A handful of nodes per compile, none per
             * call.  Same scaffolding pattern as the multi-iterator Do. */
            if (n > FN_MAX_PARAMS || c->nscope + n > CTX_MAX_SCOPE) { c->ok = false; return false; }
            Expr* call = fn_head_call(s->head, n);
            if (!call) { c->ok = false; return false; }
            int saved_scope = c->nscope;
            for (int i = 0; i < n; i++) {
                c->scope[c->nscope].name = fn_placeholder(i);
                c->scope[c->nscope].reg = argv[i].reg;
                c->scope[c->nscope].type = argv[i].type;
                c->scope[c->nscope].built = argv[i].built;
                c->nscope++;
            }
            bool ok = emit(c, call, out);
            c->nscope = saved_scope;
            /* Blame the user's own node: the scaffolding is freed on the next
             * line, and a diagnostic pointing into it would dangle. */
            if (!ok && expr_subtree_of(call, c->bail_node)) c->bail_node = s->fexpr;
            expr_free(call);
            return ok;
        }
    }
    c->ok = false; return false;
}

/* `if (R[reg] > lim) fail` — the interpreter silently truncates a Table at
 * 1e6 elements (src/list/table.c:114), so beyond that the two would disagree
 * about the LENGTH.  Declining hands the call back and lets it truncate. */
static void emit_max_guard(Ctx* c, int reg, long long lim) {
    Slot z = { 0 }, kl; memset(&kl, 0, sizeof kl); kl.i = lim;
    int rl = alloc_temp(c), rc = alloc_temp(c);
    ins(c, OP_CONST, (uint32_t)rl, 0, 0, kl);
    ins(c, OP_GT_I, (uint32_t)rc, (uint32_t)reg, (uint32_t)rl, z);
    size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);   /* in range -> skip */
    ins(c, OP_FAIL, 0, 0, 0, z);
    if (c->ok) c->code[jz].b = (uint32_t)c->n;
    c->temp_top -= 2;
}

/* Element count of a counted integer iterator, into `dst`.
 *
 *     n = max(0, floor((hi - lo) / di) + 1)
 *
 * OP_QUOT_I is FLOOR division (it corrects the sign, compile.c's QUOT_I case),
 * which is what the interpreter's `val <= hi` / `val >= hi` walk amounts to for
 * an integer step — including the wrong-direction case, where the floor goes
 * negative and the clamp gives the empty table.  Truncating division would
 * claim one element for `{i, 1, 0, 2}`, where the interpreter produces none. */
static bool emit_iter_count(Ctx* c, const LoopSpec* s, int rlo, int dst) {
    Slot z = { 0 }, k; memset(&k, 0, sizeof k);
    int rhi = alloc_temp(c), rt = alloc_temp(c), rk = alloc_temp(c);
    if (!emit_loop_bound(c, s->lo, rlo)) return false;
    if (!emit_loop_bound(c, s->hi, rhi)) return false;
    ins(c, OP_SUB_I, (uint32_t)rt, (uint32_t)rhi, (uint32_t)rlo, z);
    k.i = s->di;  ins(c, OP_CONST, (uint32_t)rk, 0, 0, k);
    ins(c, OP_QUOT_I, (uint32_t)rt, (uint32_t)rt, (uint32_t)rk, z);
    k.i = 1;      ins(c, OP_CONST, (uint32_t)rk, 0, 0, k);
    ins(c, OP_ADD_I, (uint32_t)rt, (uint32_t)rt, (uint32_t)rk, z);
    k.i = 0;      ins(c, OP_CONST, (uint32_t)rk, 0, 0, k);
    ins(c, OP_MAX_I, (uint32_t)dst, (uint32_t)rt, (uint32_t)rk, z);
    c->temp_top -= 3;
    return c->ok;
}

/* Does `e` evaluate to Null in the interpreter?  Do/While/For/Scan run for
 * their side effects and answer Null; the compiled forms answer the integer 0,
 * which only agrees where the value is thrown away.  See compile_expr_ex. */
static bool stmt_valued_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return strcmp(h, "Do") == 0 || strcmp(h, "While") == 0
        || strcmp(h, "For") == 0 || strcmp(h, "Scan") == 0
        /* If[c, t] answers Null when c is False; the compiled form answers 0. */
        || (strcmp(h, "If") == 0 && e->data.function.arg_count == 2);
}

/* SameQ opcode for an accumulator type.
 *
 * Deliberately NOT the Equal opcode.  The interpreter's default SameTest is
 * expr_eq (src/funcprog.c:2493), and expr_eq calls two NaNs the same
 * (src/expr.c:622) — which is what lets FixedPoint terminate on an orbit that
 * reaches NaN instead of spinning to the safety cap.  Integers have no NaN, so
 * there OP_EQ_I already IS expr_eq. */
static uint16_t emit_sameq_op(CompileType t) {
    return t == CT_INT ? OP_EQ_I : t == CT_COMPLEX ? OP_SAMEQ_C : OP_SAMEQ_R;
}

/* `if (R[reg] < 0) fail` — a negative application bound leaves Nest and
 * FixedPoint UNEVALUATED in the interpreter (src/funcprog.c:2153), where a
 * counted loop would silently run zero times and return the seed. */
static void emit_nonneg_guard(Ctx* c, int reg) {
    Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
    int rz = alloc_temp(c), rc = alloc_temp(c);
    ins(c, OP_CONST, (uint32_t)rz, 0, 0, k0);
    ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)reg, (uint32_t)rz, z);   /* reg < 0 ? */
    size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);   /* in range -> skip */
    ins(c, OP_FAIL, 0, 0, 0, z);
    if (c->ok) c->code[jz].b = (uint32_t)c->n;
    c->temp_top -= 2;
}

/* `if (R[reg] == 0) fail` — Fold[f, {}] stays unevaluated (src/funcprog.c:2282),
 * so a seedless fold over an empty vector must not answer. */
static void emit_nonzero_guard(Ctx* c, int reg) {
    Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
    int rz = alloc_temp(c), rc = alloc_temp(c);
    ins(c, OP_CONST, (uint32_t)rz, 0, 0, k0);
    ins(c, OP_LE_I, (uint32_t)rc, (uint32_t)reg, (uint32_t)rz, z);   /* reg <= 0 ? */
    size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);   /* non-empty -> skip */
    ins(c, OP_FAIL, 0, 0, 0, z);
    if (c->ok) c->code[jz].b = (uint32_t)c->n;
    c->temp_top -= 2;
}

/* The lowering proper.  Every bail is a plain `return false` from somewhere in
 * here; the `emit` wrapper below is what turns that into a diagnostic, so no
 * bail site needs to know diagnostics exist. */
static bool emit_node(Ctx* c, const Expr* e, Val* out) {
    if (!e || !c->ok) { c->ok = false; return false; }

    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) { *out = emit_const(c, imm, lt); return c->ok; }

    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        CompileType st; bool sb = false; int sr = scope_find(c, nm, &st, &sb);
        if (sr >= 0) { out->reg = sr; out->tmp = false; out->type = st; out->built = sb; return true; }
        int k = arg_find(c, nm);
        /* An argument array is borrowed from the caller, so anything derived
         * from it takes the caller's kind: NOT built. */
        if (k >= 0) { c->argdep[k] = 1; out->reg = k; out->tmp = false; out->type = c->arg_types[k]; out->built = false; return true; }
        if (strcmp(nm, "True") == 0)  { imm.i = 1; *out = emit_const(c, imm, CT_BOOL); return c->ok; }
        if (strcmp(nm, "False") == 0) { imm.i = 0; *out = emit_const(c, imm, CT_BOOL); return c->ok; }
        if (strcmp(nm, "I") == 0)     { imm.z = I; *out = emit_const(c, imm, CT_COMPLEX); return c->ok; }
        double cv;
        if (named_const(nm, &cv)) { imm.r = cv; *out = emit_const(c, imm, CT_REAL); return c->ok; }
        CompileType gt;
        if (global_const(c, nm, &imm, &gt)) { *out = emit_const(c, imm, gt); return c->ok; }
        c->ok = false; return false;
    }

    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) { c->ok = false; return false; }

    /* Already computed once at program entry?  The register is reserved below
     * the temp stack, so it is still live and is never a temp to be freed. */
    if (c->ncse) {
        int slot = cse_lookup(c, e);
        if (slot >= 0 && c->cse_ready[slot]) {
            out->reg = c->cse_reg[slot]; out->tmp = false; out->type = c->cse_type[slot];
            return true;
        }
    }

    const char* h = e->data.function.head->data.symbol.name;
    Expr** A = e->data.function.args;
    size_t na = e->data.function.arg_count;

    /* Slot[k] inside an inlined Function[body].  An index past the bound frame
     * falls through to the bail below, which is right: the interpreter would
     * leave that Slot unsubstituted, so its answer is not a machine number. */
    if (h == SYM_Slot) {
        int k = fn_slot_index(c, A, na);
        if (k >= 0) { out->reg = c->slot[k].reg; out->tmp = false; out->type = c->slot[k].type; return true; }
        c->ok = false; return false;
    }

    /* Elementwise fusion, tried before the delegated array lowering below.  The
     * `array_args` gate keeps this free for the overwhelmingly common scalar
     * program: with no array arguments there is nothing to fuse and we never
     * even run the type probe.  On failure try_fuse rolls back cleanly and the
     * delegated path handles the node exactly as before. */
    if (c->array_args && !c->fusing && !(c->flags & COMPILE_NO_FUSE)) {
        Val fv;
        if (try_fuse(c, e, &fv)) { *out = fv; return c->ok; }
    }

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
        CompileType mt, nt;
        if (!infer_type(c, A[0], &mt) || !infer_type(c, A[1], &nt)) { c->ok = false; return false; }
        if (mt == CT_INT && nt == CT_INT) {
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
            *out = binop(c, h[0] == 'M' ? OP_MOD_I : OP_QUOT_I, a, b, CT_INT);
            return c->ok;
        }
        if (mt <= CT_REAL && nt <= CT_REAL && mt != CT_BOOL && nt != CT_BOOL) {
            if (h[0] == 'M') {          /* the registered machine kernel for Mod */
                Val kv;
                if (try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok; }
                c->ok = false; return false;
            }
            /* Quotient[a,b] == Floor[a/b], and an INTEGER like the interpreter's. */
            Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
            coerce(c, &a, CT_REAL); coerce(c, &b, CT_REAL);
            Val q = binop(c, OP_DIV_R, a, b, CT_REAL);
            *out = unop(c, OP_FLOOR_R, q, CT_INT);
            return c->ok;
        }
        c->ok = false; return false;
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
            /* Integer exponents emit POWI directly rather than through unop(),
             * because the exponent rides in the immediate.  That means this site
             * has to make the tile decision itself — routing it through vec_unop
             * is what keeps a strip-mined `v^3` from reading a tile POINTER as a
             * double. */
            Slot s; memset(&s, 0, sizeof s); s.i = nexp;
            if (c->vector_mode && val_is_tile(a)) {
                CompileType rt = (a.type == CT_COMPLEX) ? CT_COMPLEX : CT_REAL;
                if (rt == CT_REAL) coerce(c, &a, CT_REAL);
                if (!c->ok) return false;
                *out = vec_unop(c, rt == CT_COMPLEX ? OP_POWI_C : OP_POWI_R, a, rt, s);
                return c->ok;
            }
            if (a.type == CT_INT && nexp >= 0) { free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_I, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_INT; return c->ok; }
            if (a.type == CT_COMPLEX) { free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_C, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_COMPLEX; return c->ok; }
            coerce(c, &a, CT_REAL);
            free_if_tmp(c, a); int d = alloc_temp(c); ins(c, OP_POWI_R, (uint32_t)d, (uint32_t)a.reg, 0, s); out->reg = d; out->tmp = true; out->type = CT_REAL; return c->ok;
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
        /* Integer base AND integer exponent, the exponent's value known only per
         * call.  `2^n` is an Integer for n >= 0 and a Rational below it, so the
         * only faithful lowering is one that computes exactly and abandons the
         * call on a negative exponent (and on 0^0, which is Indeterminate). */
        if (t == CT_INT && !val_is_tile(a) && !val_is_tile(b)) {
            *out = binop(c, OP_POW_II, a, b, CT_INT);
            return c->ok;
        }
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
    /* Integer-closed heads, before the generic kernel dispatch — see INT_CLOSED.
     * Declining here (rather than bailing) leaves every non-integer call to the
     * ordinary real kernel it already used. */
    {
        const IntClosed* ic = int_closed_head(h, na);
        if (ic && all_args_int(c, A, na)) {
            Val a; if (!emit(c, A[0], &a)) return false;
            if (na == 1) { *out = unop(c, ic->op, a, CT_INT); return c->ok; }
            Val b; if (!emit(c, A[1], &b)) return false;
            *out = binop(c, ic->op, a, b, CT_INT);
            return c->ok;
        }
    }
    /* Integer-only heads — see int_only_head. */
    {
        bool pred = false;
        if (int_only_head(h, na, &pred) && all_args_int(c, A, na)) {
            /* GCD / LCM: Flat and n-ary, folded left the way the interpreter's
             * own Flat attribute would associate them. */
            if (strcmp(h, "GCD") == 0 || strcmp(h, "LCM") == 0) {
                uint16_t op = (h[0] == 'G') ? OP_GCD_I : OP_LCM_I;
                Val acc; if (!emit(c, A[0], &acc)) return false;
                for (size_t i = 1; i < na; i++) {
                    Val v; if (!emit(c, A[i], &v)) return false;
                    acc = binop(c, op, acc, v, CT_INT);
                    if (!c->ok) return false;
                }
                *out = acc; return c->ok;
            }
            if (strcmp(h, "IntegerLength") == 0 || strcmp(h, "IntegerExponent") == 0) {
                uint16_t op = (h[7] == 'L') ? OP_ILEN_I : OP_IEXP_I;
                Val a; if (!emit(c, A[0], &a)) return false;
                Val b;
                if (na == 2) { if (!emit(c, A[1], &b)) return false; }
                else {
                    /* The default base is 10 for both. */
                    memset(&imm, 0, sizeof imm); imm.i = 10;
                    b = emit_const(c, imm, CT_INT);
                    if (!c->ok) return false;
                }
                *out = binop(c, op, a, b, CT_INT); return c->ok;
            }
            if (strcmp(h, "PowerMod") == 0) {
                /* Ternary: three CONSECUTIVE registers, the K_NARY shape. Same
                 * allocate-then-MOVE dance as the n-ary kernel path, and for the
                 * same reason — each argument's own lowering needs temps above
                 * the run, so the run has to be reserved first. */
                int base = alloc_temp(c);
                (void)alloc_temp(c); (void)alloc_temp(c);
                int after = c->temp_top;
                Slot z; memset(&z, 0, sizeof z);
                for (size_t i = 0; i < 3 && c->ok; i++) {
                    Val v;
                    if (!emit(c, A[i], &v)) return false;
                    ins(c, OP_MOVE, (uint32_t)(base + (int)i), (uint32_t)v.reg, 0, z);
                    c->temp_top = after;
                }
                if (!c->ok) return false;
                c->temp_top = (base - c->nlocals) + 1;
                ins_f(c, OP_POWMOD_I, 3, (uint32_t)base, (uint32_t)base, 0, z);
                out->reg = base; out->tmp = true; out->type = CT_INT;
                return c->ok;
            }
            /* EvenQ / OddQ / Divisible are Mod-and-compare, so they need no
             * opcode of their own — and they inherit MOD_I's guard, which hands
             * `Divisible[n, 0]` to the interpreter rather than dividing by it. */
            {
                Val a; if (!emit(c, A[0], &a)) return false;
                Val b;
                if (na == 2) { if (!emit(c, A[1], &b)) return false; }
                else {
                    memset(&imm, 0, sizeof imm); imm.i = 2;
                    b = emit_const(c, imm, CT_INT);
                    if (!c->ok) return false;
                }
                Val m = binop(c, OP_MOD_I, a, b, CT_INT);
                if (!c->ok) return false;
                memset(&imm, 0, sizeof imm); imm.i = (strcmp(h, "OddQ") == 0) ? 1 : 0;
                Val k = emit_const(c, imm, CT_INT);
                if (!c->ok) return false;
                *out = binop(c, OP_EQ_I, m, k, CT_BOOL);
                return c->ok;
            }
        }
    }
    if (strcmp(h, "Sign") == 0 && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_INT)  { *out = unop(c, OP_SIGN_I, a, CT_INT); return c->ok; }
        /* CT_INT, not CT_REAL: `Sign[-2.5]` is the Integer -1 in the
         * interpreter, and a compiled path answering -1. would differ in HEAD
         * from the one it is supposed to be interchangeable with.  Same reason
         * UnitStep is lowered by hand. */
        if (a.type == CT_REAL) { *out = unop(c, OP_SIGN_R, a, CT_INT); return c->ok; }
        if (a.type == CT_COMPLEX) {
            /* z/|z|, and 0 at the origin — already in the shared kernel
             * registry, which this branch was shadowing by bailing first. */
            SymbolDef* d = symtab_lookup("Sign");
            const NDUnaryKernel* k = d ? (const NDUnaryKernel*)d->ndarray_unary_kernel : NULL;
            if (!k || !k->cplx) { c->ok = false; return false; }
            *out = kern_unop(c, OP_KERN_CC, a, CT_COMPLEX, (const void*)k->cplx);
            return c->ok;
        }
        c->ok = false; return false;
    }
    if ((strcmp(h, "Floor") == 0 || strcmp(h, "Ceiling") == 0 || strcmp(h, "Round") == 0
         || strcmp(h, "IntegerPart") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_INT) { *out = a; return c->ok; }
        if (a.type != CT_REAL) coerce(c, &a, CT_REAL);
        /* IntegerPart is here rather than on its registered kernel because the
         * kernel returns a double and the interpreter returns an Integer —
         * `IntegerPart[2.5]` is `2`, not `2.`.  The kernel stays for the ARRAY
         * path, where a packed real buffer is the right answer. */
        uint16_t op = h[0] == 'F' ? OP_FLOOR_R : h[0] == 'C' ? OP_CEIL_R
                    : h[0] == 'R' ? OP_ROUND_R : OP_TRUNC_R;
        *out = unop(c, op, a, CT_INT); return c->ok;
    }
    /* FractionalPart of an integer is the Integer 0.  Its registered kernel
     * returns a double, so without this the head comes back as 0. — the same
     * trap IntegerPart is lowered by hand for. */
    if (strcmp(h, "FractionalPart") == 0 && na == 1) {
        CompileType at;
        if (infer_type(c, A[0], &at) && at == CT_INT) {
            Val a; if (!emit(c, A[0], &a)) return false;
            free_if_tmp(c, a); memset(&imm, 0, sizeof imm); imm.i = 0;
            *out = emit_const(c, imm, CT_INT);
            return c->ok;
        }
    }
    if ((strcmp(h, "Re") == 0 || strcmp(h, "Im") == 0 || strcmp(h, "Arg") == 0 || strcmp(h, "Conjugate") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        if (a.type == CT_INT && strcmp(h, "Conjugate") != 0) {
            /* Re[n] = n, Im[n] = the Integer 0, Arg[n] = the Integer 0 or Pi. */
            if (strcmp(h, "Re") == 0) { *out = a; return c->ok; }
            if (strcmp(h, "Im") == 0) {
                free_if_tmp(c, a); memset(&imm, 0, sizeof imm); imm.i = 0;
                *out = emit_const(c, imm, CT_INT); return c->ok;
            }
            *out = unop(c, OP_ARG_I, a, CT_INT); return c->ok;
        }
        if (a.type != CT_COMPLEX) {
            if (strcmp(h, "Im") == 0) { free_if_tmp(c, a); imm.r = 0.0; *out = emit_const(c, imm, CT_REAL); return c->ok; }
            /* Arg of a real is 0 or Pi, which the interpreter evaluates, so
             * widen and use the same opcode rather than bailing the whole body:
             * carg(x + 0i) gives exactly 0 for x > 0 and Pi for x < 0. */
            if (strcmp(h, "Arg") == 0) {
                coerce(c, &a, CT_COMPLEX);
                if (!c->ok) return false;
                *out = unop(c, OP_ARG_C, a, CT_REAL);
                return c->ok;
            }
            *out = a; return c->ok;   /* Re/Conjugate of real = itself */
        }
        if (strcmp(h, "Re") == 0)   { *out = unop(c, OP_RE_C, a, CT_REAL); return c->ok; }
        if (strcmp(h, "Im") == 0)   { *out = unop(c, OP_IM_C, a, CT_REAL); return c->ok; }
        if (strcmp(h, "Arg") == 0)  { *out = unop(c, OP_ARG_C, a, CT_REAL); return c->ok; }
        *out = unop(c, OP_CONJ_C, a, CT_COMPLEX); return c->ok;
    }
    /* UnitStep / Clip / Rescale: lowered by hand rather than registered as
     * kernels, because for these the RESULT TYPE is the whole difficulty.
     *
     *   UnitStep returns an INTEGER (UnitStep[0.5] is 1, not 1.) — a real-valued
     *   kernel would answer with a different head from the interpreter.  A
     *   comparison already leaves 0/1 in the integer half of the slot, so the
     *   lowering is just the comparison, typed CT_INT.
     *
     *   Clip and Rescale take a LIST of bounds, which no kernel signature can
     *   express; a literal two-element List is destructured here instead.  Both
     *   are gated on the bounds being REAL: with exact bounds the interpreter
     *   returns an exact value where it clips (Clip[5, {1, 3}] is the Integer 3)
     *   and a Real where it does not, which no single compiled type can match. */
    /* HypergeometricPFQ[{a1..ap}, {b1..bq}, z] — the head that actually matters,
     * because the evaluator canonicalises Hypergeometric0F1, 1F1 and 2F1 into it
     * before anything downstream sees them.  Its LIST arguments are what no
     * kernel signature can express, so the two literal Lists are destructured
     * here and flattened into the consecutive block the n-ary opcode wants,
     * with p passed as the leading element so the kernel can find the split. */
    if (strcmp(h, "HypergeometricPFQ") == 0 && na == 3) {
        const Expr* la = A[0];
        const Expr* lb = A[1];
        if (la->type != EXPR_FUNCTION || la->data.function.head->type != EXPR_SYMBOL
            || strcmp(la->data.function.head->data.symbol.name, "List") != 0
            || lb->type != EXPR_FUNCTION || lb->data.function.head->type != EXPR_SYMBOL
            || strcmp(lb->data.function.head->data.symbol.name, "List") != 0)
            { c->ok = false; return false; }
        size_t np = la->data.function.arg_count, nq = lb->data.function.arg_count;
        if (np + nq + 2 > 8) { c->ok = false; return false; }   /* KERNN operand cap */

        const Expr* parts[8];
        size_t total = 0;
        for (size_t i = 0; i < np; i++) parts[total++] = la->data.function.args[i];
        for (size_t j = 0; j < nq; j++) parts[total++] = lb->data.function.args[j];
        parts[total++] = A[2];
        for (size_t i = 0; i < total; i++) {
            CompileType t;
            if (!infer_type(c, parts[i], &t) || CT_IS_ARRAY(t) || t == CT_BOOL)
                { c->ok = false; return false; }
        }

        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_nary_kernel) { c->ok = false; return false; }
        const NDNaryKernel* k = (const NDNaryKernel*)d->ndarray_nary_kernel;
        if (!k->cplx) { c->ok = false; return false; }

        Slot z; memset(&z, 0, sizeof z);
        int base = alloc_temp(c);                     /* slot 0 holds p */
        for (size_t i = 0; i < total; i++) (void)alloc_temp(c);
        int after = c->temp_top;
        Slot pk; memset(&pk, 0, sizeof pk); pk.r = (double)np;
        ins(c, OP_CONST, (uint32_t)base, 0, 0, pk);
        for (size_t i = 0; i < total && c->ok; i++) {
            Val v;
            if (!emit(c, parts[i], &v)) return false;
            coerce(c, &v, CT_REAL);
            if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)(base + 1 + (int)i), (uint32_t)v.reg, 0, z);
            c->temp_top = after;
        }
        if (!c->ok) return false;
        c->temp_top = (base - c->nlocals) + 1;
        Slot kp; memset(&kp, 0, sizeof kp); kp.p = (const void*)k->cplx;
        ins_f(c, OP_KERNN, (uint16_t)(total + 1), (uint32_t)base, (uint32_t)base, 0, kp);
        out->reg = base; out->tmp = true; out->type = CT_REAL;
        return c->ok;
    }

    if (strcmp(h, "UnitStep") == 0 && na >= 1) {
        Val acc; bool have = false;
        for (size_t i = 0; i < na; i++) {
            CompileType at;
            if (!infer_type(c, A[i], &at) || CT_IS_ARRAY(at) || at == CT_BOOL) { c->ok = false; return false; }
            Val v; if (!emit(c, A[i], &v)) return false;
            Val z;
            if (at == CT_INT) { Slot k; memset(&k, 0, sizeof k); z = emit_const(c, k, CT_INT); }
            else { coerce(c, &v, CT_REAL); Slot k; memset(&k, 0, sizeof k); k.r = 0.0; z = emit_const(c, k, CT_REAL); }
            if (!c->ok) return false;
            Val ge = binop(c, (at == CT_INT) ? OP_GE_I : OP_GE_R, v, z, CT_INT);
            acc = have ? binop(c, OP_MUL_I, acc, ge, CT_INT) : ge;
            have = true;
        }
        if (!have) { c->ok = false; return false; }
        *out = acc;
        return c->ok;
    }
    if ((strcmp(h, "Clip") == 0 || strcmp(h, "Rescale") == 0) && na == 2) {
        const Expr* bnd = A[1];
        if (bnd->type != EXPR_FUNCTION || bnd->data.function.head->type != EXPR_SYMBOL
            || strcmp(bnd->data.function.head->data.symbol.name, "List") != 0
            || bnd->data.function.arg_count != 2) { c->ok = false; return false; }
        CompileType tx, tl, tu;
        if (!infer_type(c, A[0], &tx) || !infer_type(c, bnd->data.function.args[0], &tl)
            || !infer_type(c, bnd->data.function.args[1], &tu)) { c->ok = false; return false; }
        /* Real bounds only — see the note above. */
        if (tl != CT_REAL || tu != CT_REAL || CT_IS_ARRAY(tx) || tx == CT_BOOL) { c->ok = false; return false; }
        /* Rescale is just (x - lo)/(hi - lo), which is defined for a complex x
         * and which the interpreter evaluates (`Rescale[1. + I, {0., 2.}]` is
         * `0.5 + 0.5 I`).  Clip is NOT: the interpreter leaves `Clip[1. + I,
         * {0., 2.}]` unevaluated, because Min/Max need an order. */
        bool cx = (tx == CT_COMPLEX);
        if (cx && h[0] == 'C') { c->ok = false; return false; }
        /* Temporaries are a STACK: binop pops its two operands and allocates the
         * destination, so both operands must be the top of it at that moment.
         * Everything below is therefore emitted in exactly the order it is
         * consumed — and `lo`, which Rescale needs twice, is emitted twice
         * rather than held across an intervening allocation.  (It is a literal
         * bound; the optimiser's value numbering folds the duplicate away.) */
        const Expr* elo = bnd->data.function.args[0];
        const Expr* ehi = bnd->data.function.args[1];
        Val x, lo, hi;
        if (!emit(c, A[0], &x)) return false;
        coerce(c, &x, cx ? CT_COMPLEX : CT_REAL);
        if (!c->ok) return false;

        if (h[0] == 'C') {                       /* Clip[x, {lo, hi}] = Min[Max[x, lo], hi] */
            if (!emit(c, elo, &lo)) return false;
            Val mx = binop(c, OP_MAX_R, x, lo, CT_REAL);
            if (!emit(c, ehi, &hi)) return false;
            *out = binop(c, OP_MIN_R, mx, hi, CT_REAL);
        } else {                                 /* Rescale[x, {lo, hi}] = (x - lo)/(hi - lo) */
            CompileType t = cx ? CT_COMPLEX : CT_REAL;
            if (!emit(c, elo, &lo)) return false;
            if (cx) coerce(c, &lo, CT_COMPLEX);
            Val num = binop(c, cx ? OP_SUB_C : OP_SUB_R, x, lo, t);
            Val lo2, hi2;
            if (!emit(c, ehi, &hi2)) return false;
            if (cx) coerce(c, &hi2, CT_COMPLEX);
            if (!emit(c, elo, &lo2)) return false;
            if (cx) coerce(c, &lo2, CT_COMPLEX);
            Val den = binop(c, cx ? OP_SUB_C : OP_SUB_R, hi2, lo2, t);
            *out = binop(c, cx ? OP_DIV_C : OP_DIV_R, num, den, t);
        }
        return c->ok;
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

    /* ConstantArray[v, dims] (M3c): a zeroed buffer, then a fill loop only when
     * the value is not already zero.  Emitting the fill as ordinary bytecode
     * (rather than as a memset variant per element type) means a non-constant
     * fill costs nothing extra and the optimiser hoists it like any other loop. */
    if (strcmp(h, "ConstantArray") == 0 && na == 2) {
        int rank; CompileType elem;
        if (!const_array_shape(c, (const Expr* const*)A, &rank, &elem)) { c->ok = false; return false; }
        const Expr* d = A[1];
        const Expr* const* dexpr = &d;                 /* the bare-dimension case */
        if (rank > 1 || (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL
                         && d->data.function.head->data.symbol.name == SYM_List))
            dexpr = (const Expr* const*)d->data.function.args;

        int base = -1;
        for (int i = 0; i < rank; i++) {               /* dims into consecutive regs */
            int r = alloc_temp(c);
            if (base < 0) base = r;
            Val dv;
            if (!emit(c, dexpr[i], &dv)) return false;
            if (dv.type != CT_INT) { c->ok = false; return false; }
            Slot z; memset(&z, 0, sizeof z);
            ins(c, OP_MOVE, (uint32_t)r, (uint32_t)dv.reg, 0, z);
            free_if_tmp(c, dv);
        }
        Slot el; memset(&el, 0, sizeof el); el.i = (long long)elem;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, (uint16_t)rank, (uint32_t)rout, (uint32_t)base, (uint32_t)base, el);

        Slot fz; memset(&fz, 0, sizeof fz);
        bool zero_fill = (A[0]->type == EXPR_REAL && A[0]->data.real == 0.0
                          && !signbit(A[0]->data.real));
        if (!zero_fill) {
            /* for k = 0 .. size-1: out[k] = v   (v hoisted by LICM if invariant) */
            int rn = alloc_temp(c), rk = alloc_temp(c);
            ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)rout, 0, fz);
            Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk, (uint32_t)rn, fz);
            size_t jz = c->n;
            ins(c, OP_JZ, 0, (uint32_t)rc, 0, fz);
            c->temp_top--;                              /* rc dead after the guard */
            size_t body = c->n;
            Val fv;
            if (!emit(c, A[0], &fv)) return false;
            coerce(c, &fv, elem);
            ins(c,a_store_op(elem),
                (uint32_t)rout, (uint32_t)rk, (uint32_t)fv.reg, fz);
            free_if_tmp(c, fv);
            Slot one; memset(&one, 0, sizeof one); one.i = 1;
            ins(c, OP_LOOP, (uint32_t)rk, (uint32_t)rn, (uint32_t)body, one);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
            c->temp_top -= 2;                           /* rn, rk */
        }
        c->temp_top -= rank;                            /* the dimension registers */
        /* A construction site: the interpreter's ConstantArray returns a List
         * whatever the arguments were, so the result kind must not follow them. */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(elem, rank); out->built = true;
        return c->ok;
    }

    /* Table[body, spec1, ..., speck]: build a rank-k machine array.
     *
     * INTEGER iterators only, and that is not a convenience restriction.  The
     * interpreter walks a real iterator by repeated addition against a 1e-14
     * termination slack (src/list/table.c:90), so a closed-form `lo + k di`
     * differs from it in the last bits and, near the endpoint, in the element
     * COUNT.  Integer bounds reproduce that walk exactly.
     *
     * An INTEGER body is compiled, not refused: `Table[i, {i, 1, n}]` holds
     * exact Integers in the interpreter and a packed NDT_INT64 buffer holds
     * exactly those.  (Before the integer dtype existed this had to bail, since
     * a float64 buffer would have answered differently rather than faster.)
     *
     * One flat store index across k nested loops, innermost varying fastest —
     * which is row-major, and which is how the interpreter's nested rewrite
     * orders the elements too (see the multi-iterator Do). */
    if (strcmp(h, "Table") == 0 && na >= 2 && (int)na - 1 <= CT_MAX_RANK) {
        int rank = (int)na - 1;
        LoopSpec sp[CT_MAX_RANK];
        for (int j = 0; j < rank; j++)
            if (!loop_spec_parse(A[j + 1], &sp[j]) || !loop_spec_int_bounds(c, &sp[j]))
                { c->ok = false; return false; }
        if (c->nscope + rank > CTX_MAX_SCOPE) { c->ok = false; return false; }

        int pushed = 0;
        for (int j = 0; j < rank; j++)
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = 0;
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++; pushed++;
            }
        CompileType et; bool okT = infer_type(c, A[0], &et);
        c->nscope -= pushed;
        if (!okT || !ct_is_elem(et)) { c->ok = false; return false; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;

        /* The dimension registers must be ONE contiguous run: A_NEW reads rank
         * of them starting at `base`. */
        int rn[CT_MAX_RANK], rlo[CT_MAX_RANK], base = -1;
        for (int j = 0; j < rank; j++) { int r = alloc_temp(c); if (base < 0) base = r; rn[j] = r; }
        for (int j = 0; j < rank; j++) rlo[j] = alloc_temp(c);
        for (int j = 0; j < rank; j++) {
            if (!emit_iter_count(c, &sp[j], rlo[j], rn[j])) return false;
            emit_max_guard(c, rn[j], 1000000);
        }

        Slot el; memset(&el, 0, sizeof el); el.i = (long long)et;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, (uint16_t)rank, (uint32_t)rout, (uint32_t)base, (uint32_t)base, el);

        int rflat = alloc_temp(c);
        ins(c, OP_CONST, (uint32_t)rflat, 0, 0, k0);

        int riv[CT_MAX_RANK], rk[CT_MAX_RANK];
        size_t Lp[CT_MAX_RANK], jz[CT_MAX_RANK];
        int scope_entry = c->nscope;
        for (int j = 0; j < rank; j++) {
            riv[j] = alloc_temp(c); rk[j] = alloc_temp(c);
            ins(c, OP_MOVE, (uint32_t)riv[j], (uint32_t)rlo[j], 0, z);
            ins(c, OP_CONST, (uint32_t)rk[j], 0, 0, k0);
            if (sp[j].var) {
                c->scope[c->nscope].name = sp[j].var; c->scope[c->nscope].reg = riv[j];
                c->scope[c->nscope].type = CT_INT; c->scope[c->nscope].built = false;
                c->nscope++;
            }
            Lp[j] = c->n;
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk[j], (uint32_t)rn[j], z);
            jz[j] = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            c->temp_top--;
        }

        Val bv;
        if (!emit(c, A[0], &bv)) { c->nscope = scope_entry; return false; }
        if (CT_IS_ARRAY(bv.type)) { c->nscope = scope_entry; c->ok = false; return false; }
        coerce(c, &bv, et);
        if (!c->ok) { c->nscope = scope_entry; return false; }
        ins(c,a_store_op(et),
            (uint32_t)rout, (uint32_t)rflat, (uint32_t)bv.reg, z);
        free_if_tmp(c, bv);
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rflat, 0, 0, one);

        for (int j = rank - 1; j >= 0; j--) {
            Slot sd; memset(&sd, 0, sizeof sd); sd.i = sp[j].di;
            ins(c, OP_INC_I, (uint32_t)riv[j], 0, 0, sd);
            ins(c, OP_INC_I, (uint32_t)rk[j], 0, 0, one);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp[j], z);
            if (c->ok) c->code[jz[j]].b = (uint32_t)c->n;
        }
        c->nscope = scope_entry;
        c->temp_top = base - c->nlocals;      /* every scalar temp is dead; rout is an arr */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(et, rank); out->built = true;
        return c->ok;
    }

    /* Range[hi] / Range[lo, hi] / Range[lo, hi, di] — a rank-1 integer array.
     *
     * This is `Table[i, {i, lo, hi, di}]` with the body being the iterator
     * itself, so it shares the iteration machinery (emit_iter_count reproduces
     * the interpreter's element COUNT exactly) and simply stores the loop
     * variable instead of a body.
     *
     * INTEGER bounds only, for the reason Table gives: the interpreter walks a
     * real iterator by repeated addition against a termination slack, so a
     * closed-form step differs from it in the last bits and, at the endpoint, in
     * the number of elements. */
    if (strcmp(h, "Range") == 0 && na >= 1 && na <= 3) {
        LoopSpec s;
        if (!range_spec(c, A, na, &s)) { c->ok = false; return false; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot one; memset(&one, 0, sizeof one); one.i = 1;

        int rn = alloc_temp(c), rlo = alloc_temp(c);
        if (!emit_iter_count(c, &s, rlo, rn)) return false;
        emit_max_guard(c, rn, VM_ITER_SAFETY_CAP);

        Slot el; memset(&el, 0, sizeof el); el.i = (long long)CT_INT;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rn, (uint32_t)rn, el);

        int riv = alloc_temp(c), rk = alloc_temp(c);
        ins(c, OP_MOVE,  (uint32_t)riv, (uint32_t)rlo, 0, z);
        ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);

        size_t top = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rk, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        ins(c, OP_A_STORE_I, (uint32_t)rout, (uint32_t)rk, (uint32_t)riv, z);
        Slot sd; memset(&sd, 0, sizeof sd); sd.i = s.di;
        ins(c, OP_INC_I, (uint32_t)riv, 0, 0, sd);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)top, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        c->temp_top = rn - c->nlocals;         /* every scalar temp is dead */
        out->reg = rout; out->tmp = true;
        out->type = CT_ARRAY(CT_INT, 1); out->built = true;
        return c->ok;
    }

    /* Part[a, subscripts...] (M3c).  See the block comment above emit_flat_index
     * for why this splits in two. */
    if (strcmp(h, "Part") == 0 && na >= 2) {
        Val a;
        if (!emit(c, A[0], &a)) return false;
        if (!CT_IS_ARRAY(a.type) || (int)(na - 1) > CT_RANK(a.type)) { c->ok = false; return false; }
        CompileType elem = CT_ELEM(a.type);
        Slot z; memset(&z, 0, sizeof z);

        if (part_is_scalar_indexed(c, a.type, (const Expr* const*)A + 1, na - 1)) {
            int ridx;
            if (!emit_flat_index(c, a, (const Expr* const*)A + 1, na - 1, &ridx)) return false;
            int rd = alloc_temp(c);
            ins(c,a_load_op(elem),
                (uint32_t)rd, (uint32_t)a.reg, (uint32_t)ridx, z);
            /* rd sits above ridx, so drop ridx by relocating rd onto it. */
            ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)rd, 0, z);
            c->temp_top--;                              /* rd */
            free_if_tmp(c, a);
            out->reg = ridx; out->tmp = true; out->type = elem;
            return c->ok;
        }

        int keep = CT_RANK(a.type) - (int)(na - 1);
        for (size_t i = 1; i < na; i++) if (subscript_is_literal_spec(A[i])) keep++;
        if (keep < 1 || keep > CT_MAX_RANK) { c->ok = false; return false; }

        int base = -1;
        PartSpec* ps = emit_partspec(c, (const Expr* const*)A + 1, na - 1, &base);
        if (!ps) return false;
        if (!ctx_own_partspec(c, ps)) return false;
        Slot ip; memset(&ip, 0, sizeof ip); ip.p = ps;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_PART, (uint16_t)(na - 1), (uint32_t)rout, (uint32_t)base,
              (uint32_t)a.reg, ip);
        c->temp_top -= (int)(na - 1);                   /* the subscript registers */
        if (a.tmp) {
            /* The source was itself a temporary, and it sits BELOW the result in
             * the array stack, so it cannot simply be popped: the next alloc_arr
             * would hand out the register the result is living in.  Free it and
             * slide the result down into its slot, restoring LIFO. */
            ins(c, OP_ARR_FREE, (uint32_t)a.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)a.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = a.reg;
        }
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(elem, keep); out->built = a.built;
        return c->ok;
    }

    /* Total[v] / Length[v]: the array -> scalar reductions.  Total delegates to
     * the NDArray reduction so its summation order — and therefore its
     * rounding — is identical to the interpreter's Total[]. */
    if ((strcmp(h, "Total") == 0 || strcmp(h, "Length") == 0) && na == 1) {
        Val a; if (!emit(c, A[0], &a)) return false;
        /* Length is dims[0] at any rank — the number of ROWS of a matrix, as in
         * the interpreter.  Total stays rank-1: ndred_total_all collapses every
         * axis, but Total[] reduces only the leading one, so at rank 2 the two
         * would disagree. */
        if (!CT_IS_ARRAY(a.type)
            || (h[0] == 'T' && CT_RANK(a.type) != 1)) { c->ok = false; return false; }
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

    /*
     * If[cond, then] -- no else branch. A STATEMENT, not a value: the
     * interpreter answers Null when cond is False, so this lowers exactly like
     * While (guard, jump over the body, answer the integer 0) and joins
     * stmt_valued_head so the 0 can never be the program's result.
     *
     * Omitting this was not a small gap. The compilable subset is a cliff, not a
     * slope -- one head outside it costs the WHOLE body -- and `If[test, var =
     * val]` is how anyone writes a running maximum or a conditional store. Both
     * a Sieve of Eratosthenes and a Collatz longest-chain search failed to
     * compile on this single form, and so ran interpreted: Collatz to 10^6 took
     * 240 s instead of 0.35 s. docs/design/compile.md section 11 specified it
     * ("If[c,t] -> f = Null only valid where the value is unused"); it was simply
     * never emitted.
     */
    if (strcmp(h, "If") == 0 && na == 2) {
        Slot z = { 0 };
        Val cond; if (!emit(c, A[0], &cond)) return false;
        if (cond.type != CT_BOOL) { c->ok = false; return false; }
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)cond.reg, 0, z);
        free_if_tmp(c, cond);
        Val th; if (!emit(c, A[1], &th)) return false;
        free_if_tmp(c, th);                    /* value discarded, as While's body is */
        if (c->ok) c->code[jz].b = (uint32_t)c->n;   /* end label */
        int r0 = alloc_temp(c); Slot s0; s0.i = 0;
        ins(c, OP_CONST, (uint32_t)r0, 0, 0, s0);
        out->reg = r0; out->tmp = true; out->type = CT_INT;
        return c->ok;
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
        LoopSpec s;
        /* Named iterator required — see the matching note in infer_type. */
        if (!loop_spec_parse(A[1], &s) || !s.var || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) { c->ok = false; return false; }  /* integer iteration only */
        /* body type with i bound as INT */
        c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = 0;
        c->scope[c->nscope].type = CT_INT; c->nscope++;
        CompileType T; bool okT = infer_type(c, A[0], &T);
        c->nscope--;
        /* An array accumulator would need a per-iteration copy, not a MOVE. */
        if (!okT || T == CT_BOOL || CT_IS_ARRAY(T)) { c->ok = false; return false; }
        Slot z = { 0 };
        int racc = alloc_temp(c), rhi = alloc_temp(c), ri = alloc_temp(c);
        Slot iz; iz.i = 0; if (T == CT_INT) iz.i = prod ? 1 : 0; else if (T == CT_REAL) iz.r = prod ? 1.0 : 0.0; else iz.z = prod ? 1.0 : 0.0;
        ins(c, OP_CONST, (uint32_t)racc, 0, 0, iz);
        if (!emit_loop_bound(c, s.lo, ri)) return false;
        if (!emit_loop_bound(c, s.hi, rhi)) return false;
        c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = ri;
        c->scope[c->nscope].type = CT_INT; c->nscope++;
        size_t L = c->n;
        int rc = alloc_temp(c);
        ins(c, s.di > 0 ? OP_LE_I : OP_GE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;                                          /* free the guard temp */
        Val rb; if (!emit(c, A[0], &rb)) { c->nscope--; return false; }
        if (CT_IS_ARRAY(rb.type)) { c->nscope--; c->ok = false; return false; }
        coerce(c, &rb, T);
        uint16_t acc = prod ? (T == CT_INT ? OP_MUL_I : T == CT_REAL ? OP_MUL_R : OP_MUL_C)
                            : (T == CT_INT ? OP_ADD_I : T == CT_REAL ? OP_ADD_R : OP_ADD_C);
        ins(c, acc, (uint32_t)racc, (uint32_t)racc, (uint32_t)rb.reg, z);
        free_if_tmp(c, rb);
        Slot step; step.i = s.di; ins(c, OP_INC_I, (uint32_t)ri, 0, 0, step);
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
            || (int)(c->nscope + L->data.function.arg_count) > CTX_MAX_SCOPE) { c->ok = false; return false; }
        size_t nl = L->data.function.arg_count;
        Slot z = { 0 };
        int base_reg = -1, pushed = 0;
        /* Array locals (M3c) get a register in the ARRAY bank on top of their
         * scalar slot, which stays allocated but unused: keeping every local's
         * scalar register in one run is what makes the single temp_top reset at
         * the end correct whatever mix of kinds the locals are. */
        int arr_entry = c->arr_top, arr_regs[16], narr = 0;
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
                if (CT_IS_ARRAY(iv.type)) {
                    /* A local may be written through (u[[i]] = ...), so it has to
                     * OWN its array.  An initialiser that is already a temporary
                     * is adopted by leaving its slot allocated for the life of
                     * the scope; anything else — an argument, an enclosing local
                     * — is borrowed and is copied, which is exactly the value
                     * semantics the interpreter gives the same code. */
                    if (narr >= (int)(sizeof arr_regs / sizeof arr_regs[0])) {
                        c->nscope -= pushed; c->ok = false; return false;
                    }
                    int areg;
                    if (iv.tmp && reg_is_owned_arr(iv.reg)) areg = iv.reg;
                    else {
                        areg = alloc_arr(c);
                        uint16_t f = (uint16_t)(((unsigned)CT_ELEM(iv.type) & 3u) << AF_R_SHIFT);
                        ins_f(c, OP_A_COPY, f, (uint32_t)areg, (uint32_t)iv.reg, 0, z);
                    }
                    arr_regs[narr++] = areg;
                    c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = areg;
                    c->scope[c->nscope].type = iv.type; c->scope[c->nscope].built = iv.built;
                    c->nscope++; pushed++;
                    continue;
                }
                vt = iv.type; ins(c, OP_MOVE, (uint32_t)reg, (uint32_t)iv.reg, 0, z); free_if_tmp(c, iv);
            } else { Slot s; s.r = 0.0; ins(c, OP_CONST, (uint32_t)reg, 0, 0, s); }
            c->scope[c->nscope].name = vname; c->scope[c->nscope].reg = reg; c->scope[c->nscope].type = vt;
            c->nscope++; pushed++;
        }
        Val body; if (!emit(c, A[1], &body)) { c->nscope -= pushed; return false; }
        c->nscope -= pushed;

        if (CT_IS_ARRAY(body.type)) {
            /* The result has to outlive the frees below, so every array local
             * EXCEPT the one carrying it is released, the bank is wound back to
             * where the scope started, and the value is moved down into the
             * first slot — a handle transfer, not a copy. */
            for (int i = narr - 1; i >= 0; i--)
                if (arr_regs[i] != body.reg) ins(c, OP_ARR_FREE, (uint32_t)arr_regs[i], 0, 0, z);
            c->arr_top = arr_entry;
            int rres = alloc_arr(c);
            if (!reg_is_owned_arr(body.reg)) {
                /* borrowed (an argument array): copy, or the caller's node would
                 * be freed twice — once by them and once as our result */
                uint16_t f = (uint16_t)(((unsigned)CT_ELEM(body.type) & 3u) << AF_R_SHIFT);
                ins_f(c, OP_A_COPY, f, (uint32_t)rres, (uint32_t)body.reg, 0, z);
            } else if (rres != body.reg) {
                ins(c, OP_A_XFER, (uint32_t)rres, (uint32_t)body.reg, 0, z);
            }
            c->temp_top = (base_reg - c->nlocals);   /* no scalar result to keep */
            out->reg = rres; out->tmp = true; out->type = body.type; out->built = body.built;
            return c->ok;
        }

        for (int i = narr - 1; i >= 0; i--) ins(c, OP_ARR_FREE, (uint32_t)arr_regs[i], 0, 0, z);
        c->arr_top = arr_entry;
        if (body.reg != base_reg) ins(c, OP_MOVE, (uint32_t)base_reg, (uint32_t)body.reg, 0, z);
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
        else if (strcmp(h, "DivideBy") == 0) kind = 6;
        /* 4 and 5 are the unary forms (v++ / v--); everything else takes a value. */
        #define IC_UNARY(k) ((k) == 4 || (k) == 5)
        /* Part assignment (M3c): u[[i, j]] = v, and the compound forms.
         *
         * The target must be an array the program OWNS, because the write goes
         * into the buffer in place.  An argument array is borrowed — the caller
         * still holds the node it passed in — so writing through it would mutate
         * a value the caller never offered up, and for a List argument (packed
         * into a temporary at the boundary) the write would silently vanish.
         * Copy it into a local first; that is what the interpreter's value
         * semantics do anyway. */
        if (kind >= 0 && !IC_UNARY(kind) && na == 2 && A[0]->type == EXPR_FUNCTION
            && A[0]->data.function.head->type == EXPR_SYMBOL
            && strcmp(A[0]->data.function.head->data.symbol.name, "Part") == 0
            && A[0]->data.function.arg_count >= 2) {
            const Expr* pt = A[0];
            const Expr* const* S = (const Expr* const*)pt->data.function.args;
            size_t ns = pt->data.function.arg_count - 1;
            Val arr;
            if (!emit(c, S[0], &arr)) return false;
            if (!CT_IS_ARRAY(arr.type) || !reg_is_owned_arr(arr.reg) || arr.tmp
                || (int)ns > CT_RANK(arr.type)) { c->ok = false; return false; }
            CompileType elem = CT_ELEM(arr.type);
            Slot z; memset(&z, 0, sizeof z);

            if (part_is_scalar_indexed(c, arr.type, S + 1, ns)) {
                int ridx;
                if (!emit_flat_index(c, arr, S + 1, ns, &ridx)) return false;
                Val val;
                if (kind == 0) {                         /* plain Set */
                    if (!emit(c, A[1], &val)) return false;
                } else {                                 /* AddTo / SubtractFrom / TimesBy */
                    int rold = alloc_temp(c);
                    ins(c,a_load_op(elem),
                        (uint32_t)rold, (uint32_t)arr.reg, (uint32_t)ridx, z);
                    Val cur = { rold, true, elem, false }, rhs;
                    if (!emit(c, A[1], &rhs)) return false;
                    coerce(c, &rhs, elem);
                    uint16_t op = kind == 1 ? (elem == CT_COMPLEX ? OP_ADD_C : OP_ADD_R)
                                : kind == 2 ? (elem == CT_COMPLEX ? OP_SUB_C : OP_SUB_R)
                                : kind == 3 ? (elem == CT_COMPLEX ? OP_MUL_C : OP_MUL_R)
                                            : (elem == CT_COMPLEX ? OP_DIV_C : OP_DIV_R);
                    val = binop(c, op, cur, rhs, elem);
                }
                coerce(c, &val, elem);
                if (!c->ok) return false;
                ins(c,a_store_op(elem),
                    (uint32_t)arr.reg, (uint32_t)ridx, (uint32_t)val.reg, z);
                /* Set returns the stored value; relocate it onto the index
                 * register so the whole subscript computation is reclaimed. */
                ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)val.reg, 0, z);
                pop_tmp(c, val);
                c->temp_top = (ridx - c->nlocals) + 1;
                out->reg = ridx; out->tmp = true; out->type = elem;
                return c->ok;
            }

            /* General spec: only a plain Set, because the compound forms would
             * have to read a whole slice, combine it and write it back — which
             * is Part-as-an-lvalue on a sub-array, not this. */
            if (kind != 0) { c->ok = false; return false; }
            int base = -1;
            PartSpec* ps = emit_partspec(c, S + 1, ns, &base);
            if (!ps) return false;
            if (!ctx_own_partspec(c, ps)) return false;
            Val val;
            if (!emit(c, A[1], &val)) return false;
            if (CT_IS_ARRAY(val.type)) ps->rhs_kind = AK_ARR;
            else {
                coerce(c, &val, elem);
                if (!c->ok) return false;
                /* Same exactness rule as arr_prep's: writing 5 into an int64
                 * buffer must box as the Integer 5, or ndarray_part_set's
                 * mixed-exactness check sees an inexact right-hand side. */
                ps->rhs_kind = (elem == CT_COMPLEX) ? AK_COMPLEX
                             : (elem == CT_INT)     ? AK_INT
                                                    : AK_REAL;
            }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = ps;
            ins_f(c, OP_A_PARTSET, (uint16_t)ns, (uint32_t)arr.reg, (uint32_t)base,
                  (uint32_t)val.reg, ip);
            free_if_tmp(c, val);
            c->temp_top -= (int)ns;
            /* Set's value is the right-hand side; a general one is an array, and
             * an array result would need an owner, so this form is a statement:
             * it reports 0 and is only useful inside a CompoundExpression. */
            Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            *out = emit_const(c, k0, CT_INT);
            return c->ok;
        }

        if (kind >= 0) {
            size_t want = IC_UNARY(kind) ? 1 : 2;
            if (na != want || A[0]->type != EXPR_SYMBOL) { c->ok = false; return false; }
            CompileType vt; int vreg = scope_find(c, A[0]->data.symbol.name, &vt, NULL);
            if (vreg < 0 || vt == CT_BOOL) { c->ok = false; return false; }   /* mutable numeric locals only */
            /* DivideBy always produces a Real (or Complex), never an Int. */
            if (kind == 6 && vt == CT_INT) { c->ok = false; return false; }
            Slot z = { 0 };
            if (IC_UNARY(kind)) {
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
                              : kind == 3 ? (vt == CT_INT ? OP_MUL_I : vt == CT_REAL ? OP_MUL_R : OP_MUL_C)
                                          : (vt == CT_COMPLEX ? OP_DIV_C : OP_DIV_R);
                   ins(c, op, (uint32_t)vreg, (uint32_t)vreg, (uint32_t)val.reg, z); }
            free_if_tmp(c, val);
            out->reg = vreg; out->tmp = false; out->type = vt; return c->ok;
        }
        #undef IC_UNARY
    }

    /* Do[body, {i, lo, hi}]: counted loop for side effects; returns a dummy 0. */
    /* Do[body, spec1, spec2, ...]: the interpreter nests the iterators with the
     * LAST varying fastest, so rewrite to the nested form and lower that.  Doing
     * it here rather than teaching the loop lowering about several specs keeps
     * one implementation of the loop and makes the multi-iterator form correct
     * by the same argument as the single one. */
    if (strcmp(h, "Do") == 0 && na > 2) {
        Expr** ia = malloc((na - 1) * sizeof(Expr*));
        if (!ia) { c->ok = false; return false; }
        ia[0] = expr_copy((Expr*)A[0]);
        for (size_t i = 2; i < na; i++) ia[i - 1] = expr_copy((Expr*)A[i]);
        Expr* inner = expr_new_function(expr_new_symbol("Do"), ia, na - 1);
        free(ia);
        if (!inner) { c->ok = false; return false; }
        Expr* oa[2] = { inner, expr_copy((Expr*)A[1]) };
        Expr* outer = expr_new_function(expr_new_symbol("Do"), oa, 2);
        if (!outer) { expr_free(inner); c->ok = false; return false; }
        bool r = emit(c, outer, out);
        /* The rewrite is scaffolding, not the user's own tree: a bail inside it
         * would leave the diagnostic pointing into a node freed on the next
         * line, so blame the Do the user actually wrote. */
        if (!r && expr_subtree_of(outer, c->bail_node)) c->bail_node = e;
        expr_free(outer);
        return r;
    }
    if (strcmp(h, "Do") == 0 && na == 2) {
        LoopSpec s;
        if (!loop_spec_parse(A[1], &s) || !loop_spec_int_bounds(c, &s)
            || c->nscope >= CTX_MAX_SCOPE) { c->ok = false; return false; }
        Slot z = { 0 };
        int rhi = alloc_temp(c), ri = alloc_temp(c);
        if (!emit_loop_bound(c, s.lo, ri)) return false;
        if (!emit_loop_bound(c, s.hi, rhi)) return false;
        int pushed = 0;
        if (s.var) {
            c->scope[c->nscope].name = s.var; c->scope[c->nscope].reg = ri;
            c->scope[c->nscope].type = CT_INT; c->nscope++; pushed = 1;
        }
        /* Unit step gets the fused loop instruction: OP_LOOP increments, tests
         * and branches in ONE, where the general shape below spends four
         * instructions per iteration on control alone (test, branch, increment,
         * back-edge).  On `Do[s = s + 1. i, {i, 1, n}]` that is 8 inner
         * instructions down to 5.  OP_LOOP compares `++i < a`, so the bound
         * register holds hi + 1; the entry guard still runs once, because a
         * loop whose range is empty must execute the body zero times and
         * OP_LOOP tests only at the bottom.
         *
         * A non-unit or negative step keeps the general form — OP_LOOP steps by
         * one, and the direction of the test is baked into it. */
        if (s.di == 1) {
            Slot one; memset(&one, 0, sizeof one); one.i = 1;
            int rc = alloc_temp(c), rend = alloc_temp(c);
            /* Entry guard, ONCE: OP_LOOP tests at the bottom, so an empty range
             * would otherwise run the body a first time. */
            ins(c, OP_LE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            /* OP_LOOP compares `++i < a`, so the bound register holds hi + 1. */
            ins(c, OP_CONST, (uint32_t)rc, 0, 0, one);
            ins(c, OP_ADD_I, (uint32_t)rend, (uint32_t)rhi, (uint32_t)rc, z);
            size_t Lb = c->n;
            Val bod; if (!emit(c, A[0], &bod)) { c->nscope -= pushed; return false; }
            free_if_tmp(c, bod);
            ins(c, OP_LOOP, (uint32_t)ri, (uint32_t)rend, (uint32_t)Lb, one);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
            c->temp_top -= 2;                                /* rc, rend */
        } else {
            size_t Lp = c->n;
            int rc = alloc_temp(c);
            ins(c, s.di > 0 ? OP_LE_I : OP_GE_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rhi, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
            Val bod; if (!emit(c, A[0], &bod)) { c->nscope -= pushed; return false; }
            free_if_tmp(c, bod);
            Slot step; step.i = s.di; ins(c, OP_INC_I, (uint32_t)ri, 0, 0, step);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
            if (c->ok) c->code[jz].b = (uint32_t)c->n;
        }
        c->nscope -= pushed; c->temp_top -= 2;
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

    /* Nest[f, x, n]: apply f n times, feeding each result back in.  The
     * accumulator lives in one persistent register (racc) typed to the
     * fixed-point type; the counted loop mirrors Do.  `f` is any function value
     * fn_resolve accepts — Function[u,body], #-slots, a bare head, Composition,
     * a CompiledFunction — so this one lowering covers all of them. */
    if (strcmp(h, "Nest") == 0 && na == 3) {
        FnSpec fs;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return false; }
        CompileType tn; if (!infer_type(c, A[2], &tn) || tn != CT_INT) { c->ok = false; return false; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
        int tfp = nest_fixed_type(c, &fs, tx);
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
        /* Nest[f, x, -1] is UNEVALUATED (src/funcprog.c:2153); the counted loop
         * below would silently run zero times and hand back the seed. */
        emit_nonneg_guard(c, rn);
        Slot s0; s0.i = 0; ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, s0);
        size_t Lp = c->n;
        int rc = alloc_temp(c); ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z); c->temp_top--;
        Val acc = { racc, false, t, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, t); if (!c->ok) return false;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        free_if_tmp(c, vb);
        Slot one; one.i = 1; ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;
        c->temp_top = (racc - c->nlocals) + 1;   /* keep racc; free rn, rcnt */
        out->reg = racc; out->tmp = true; out->type = t; return c->ok;
    }

    /* Reverse / Sort / Accumulate / Flatten / Transpose / Take / Drop over a
     * machine array, delegated to the interpreter's own entry point. */
    {
        const NdFnSpec* nf = nd_fn_lookup(h, na);
        if (nf) {
            Val a;
            if (!emit(c, A[0], &a)) return false;
            CompileType rt = nd_fn_result(nf, a.type);
            if ((int)rt < 0) { c->ok = false; return false; }
            Slot z = { 0 };
            int base = 0;
            for (int i = 0; i < nf->nextra; i++) {
                int r = alloc_temp(c);
                if (i == 0) base = r;
                Val v;
                if (!emit(c, A[1 + i], &v)) return false;
                if (v.type != CT_INT) { c->ok = false; return false; }
                ins(c, OP_MOVE, (uint32_t)r, (uint32_t)v.reg, 0, z);
                free_if_tmp(c, v);
            }
            Slot ip; memset(&ip, 0, sizeof ip); ip.p = nf;
            int rout = alloc_arr(c);
            ins_f(c, OP_A_NDFN, (uint16_t)nf->nextra, (uint32_t)rout,
                  (uint32_t)base, (uint32_t)a.reg, ip);
            c->temp_top -= nf->nextra;
            if (a.tmp) {   /* restore LIFO — the same slide the Part lowering does */
                ins(c, OP_ARR_FREE, (uint32_t)a.reg, 0, 0, z);
                ins(c, OP_A_XFER, (uint32_t)a.reg, (uint32_t)rout, 0, z);
                c->arr_top--;
                rout = a.reg;
            }
            out->reg = rout; out->tmp = true; out->type = rt; out->built = a.built;
            return c->ok;
        }
    }

    /* Select / TakeWhile / LengthWhile / AllTrue / AnyTrue / NoneTrue over a
     * rank-1 array: one predicate loop, four things done with the answer.
     *
     * Select and TakeWhile produce a buffer whose length is only known once the
     * loop has run, so they allocate the upper bound the source gives them and
     * cut it to size (A_TRUNC).  An EMPTY result declines: the interpreter
     * cannot pack `{}` either, so it answers with a List, and a length-0 array
     * would be a different value.
     *
     * (These only became compilable once the interpreter itself grew NDArray
     * paths for them — every one used to return the call UNEVALUATED on a
     * packed argument, so there was nothing to be parity with.) */
    {
        int sel = -1;
        if (na == 2) {
            if      (strcmp(h, "Select") == 0)      sel = 0;
            else if (strcmp(h, "TakeWhile") == 0)   sel = 1;
            else if (strcmp(h, "LengthWhile") == 0) sel = 2;
            else if (strcmp(h, "AllTrue") == 0)     sel = 3;
            else if (strcmp(h, "AnyTrue") == 0)     sel = 4;
            else if (strcmp(h, "NoneTrue") == 0)    sel = 5;
        }
        if (sel >= 0) {
            bool keeps = (sel <= 1);            /* builds a buffer */
            bool pred_run = (sel <= 2);         /* scans a prefix / filters */
            FnSpec ps;
            if (!fn_resolve(A[1], 1, &ps)) { c->ok = false; return false; }
            CompileType el = vec_elem_type(c, A[0]);
            if ((int)el < 0) { c->ok = false; return false; }
            CompileType tb;
            if (!infer_apply(c, &ps, &el, 1, &tb) || tb != CT_BOOL) { c->ok = false; return false; }

            Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
            Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
            int rn = alloc_temp(c), ri = alloc_temp(c), rk = alloc_temp(c);
            Val va;
            if (!emit(c, A[0], &va)) return false;
            if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return false; }
            ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);

            int rout = -1;
            if (keeps) {
                Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)el;
                rout = alloc_arr(c);
                ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rn, (uint32_t)rn, el_i);
            }
            /* Boolean accumulator: All and None start True, Any starts False. */
            if (!pred_run) {
                Slot init; memset(&init, 0, sizeof init); init.i = (sel == 4) ? 0 : 1;
                ins(c, OP_CONST, (uint32_t)rk, 0, 0, init);
            } else ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);

            size_t Lp = c->n;
            int rc = alloc_temp(c);
            ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
            size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
            c->temp_top--;

            int relem = alloc_temp(c);
            int body_top = c->temp_top;
            ins(c,a_load_op(el),
                (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
            Val ev = { relem, false, el, false }, vt;
            if (!emit_apply(c, &ps, &ev, 1, &vt)) return false;
            if (vt.type != CT_BOOL) { c->ok = false; return false; }

            size_t jstop = 0; bool has_stop = false;
            if (sel == 0) {                       /* Select: keep when true */
                size_t jskip = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                ins(c,a_store_op(el),
                    (uint32_t)rout, (uint32_t)rk, (uint32_t)relem, z);
                ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);
                if (c->ok) c->code[jskip].b = (uint32_t)c->n;
            } else if (sel == 1 || sel == 2) {
                /* A prefix, so the write index and the final length are both
                 * just `ri` — no separate counter to keep in step. */
                jstop = c->n; has_stop = true;
                ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                if (sel == 1)
                    ins(c,a_store_op(el),
                        (uint32_t)rout, (uint32_t)ri, (uint32_t)relem, z);
            } else {
                /* Short-circuit.  AllTrue fires on a FALSE test, AnyTrue and
                 * NoneTrue on a TRUE one; the value they then answer with is
                 * True only for AnyTrue. */
                Slot outv; memset(&outv, 0, sizeof outv); outv.i = (sel == 4) ? 1 : 0;
                size_t jfalse = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
                size_t jcont = 0;
                if (sel == 3) {                   /* AllTrue: true -> keep going */
                    jcont = c->n; ins(c, OP_JMP, 0, 0, 0, z);
                    if (c->ok) c->code[jfalse].b = (uint32_t)c->n;   /* false -> fire */
                }
                ins(c, OP_CONST, (uint32_t)rk, 0, 0, outv);
                jstop = c->n; has_stop = true; ins(c, OP_JMP, 0, 0, 0, z);
                if (c->ok) {
                    if (sel == 3) c->code[jcont].b = (uint32_t)c->n;
                    else          c->code[jfalse].b = (uint32_t)c->n; /* false -> keep going */
                }
            }
            c->temp_top = body_top - 1;
            ins(c, OP_INC_I, (uint32_t)ri, 0, 0, k1);
            ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
            if (c->ok) {
                c->code[jz].b = (uint32_t)c->n;
                if (has_stop) c->code[jstop].b = (uint32_t)c->n;
            }

            if (keeps) {
                int rlen = (sel == 0) ? rk : ri;
                /* An empty result has no packed form; the interpreter answers
                 * with a List, so a length-0 array would not be the same value. */
                emit_nonzero_guard(c, rlen);
                ins_f(c, OP_A_TRUNC, 0, (uint32_t)rout, (uint32_t)rlen, 0, z);
                if (va.tmp) {          /* restore LIFO, as the Part lowering does */
                    ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
                    ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
                    c->arr_top--;
                    rout = va.reg;
                }
                c->temp_top = (rn - c->nlocals);
                out->reg = rout; out->tmp = true; out->type = CT_ARRAY(el, 1);
                out->built = va.built;     /* filtered, not constructed */
                return c->ok;
            }
            free_if_tmp(c, va);
            /* LengthWhile's answer is the prefix length, which IS `ri`; the
             * predicates' is the boolean in `rk`.  Relocate it down onto rn so
             * the scratch above is reclaimed. */
            ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)(sel == 2 ? ri : rk), 0, z);
            c->temp_top = (rn - c->nlocals) + 1;
            out->reg = rn; out->tmp = true; out->built = false;
            out->type = (sel == 2) ? CT_INT : CT_BOOL;
            return c->ok;
        }
    }

    /* First[v] / Last[v]: exactly v[[1]] and v[[-1]], so lower them as that
     * rather than duplicating the indexing — the range check that makes
     * First[{}] decline comes along for free. */
    if ((strcmp(h, "First") == 0 || strcmp(h, "Last") == 0) && na == 1) {
        CompileType ta;
        if (!infer_type(c, A[0], &ta) || !CT_IS_ARRAY(ta)) { c->ok = false; return false; }
        Expr* idx = expr_new_integer(h[0] == 'F' ? 1 : -1);
        Expr* args[2] = { expr_copy(A[0]), idx };
        Expr* part = expr_new_function(expr_new_symbol("Part"), args, 2);
        if (!part) { expr_free(idx); c->ok = false; return false; }
        bool r = emit(c, part, out);
        if (!r && expr_subtree_of(part, c->bail_node)) c->bail_node = e;
        expr_free(part);
        return r;
    }

    /* Map[f, v] / Scan[f, v] over a rank-1 array.
     *
     * One element loop serves both: Map stores each result into a fresh buffer,
     * Scan drops it and answers like Do.  Rank 1 only — at rank >= 2 the
     * interpreter's Map applies f to each ROW (map_ndarray_axis,
     * src/funcprog.c:272), which is a different operation from elementwise.
     *
     * The result element type must equal the SOURCE's, because Map over a packed
     * array repacks with the source dtype (src/funcprog.c:289): Map[Abs, cv] on
     * a complex vector comes back complex-typed with real values.  Anything else
     * would answer with a different element type, so it declines. */
    if ((strcmp(h, "Map") == 0 || strcmp(h, "Scan") == 0) && na == 2) {
        bool scan = h[0] == 'S';
        FnSpec fs;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return false; }
        CompileType el = vec_elem_type(c, A[1]);
        if ((int)el < 0) { c->ok = false; return false; }
        CompileType rt;
        if (!infer_apply(c, &fs, &el, 1, &rt) || CT_IS_ARRAY(rt)) { c->ok = false; return false; }
        if (!scan && rt != el) { c->ok = false; return false; }
        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;

        int rn = alloc_temp(c), ri = alloc_temp(c);
        Val va;
        if (!emit(c, A[1], &va)) return false;
        if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return false; }
        ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);
        /* An empty source makes the interpreter's Map return a plain List — the
         * result of ndarray_from_nested_list declining on List[] — where a
         * compiled buffer would be an empty NDArray. */
        if (!scan) emit_nonzero_guard(c, rn);

        /* Fast route.  When the body threads elementwise, Map[Function[u, body], v]
         * IS `body` with u bound to the whole array — and that is what the
         * existing elementwise fusion strip-mines into one tiled, optionally
         * threaded pass.  Measured 6.2x over the per-element loop below on
         * `#^2 + 1.` at 200k elements, for ten lines and no new machinery.
         *
         * Three conditions make the rewrite legal, and each corresponds to a way
         * the interpreter's Map differs from threading:
         *   - fuse_listable == 1: the interpreter threads a list exactly when the
         *     head is Listable, so `Map[Function[u, Total[u]], v]` (Total of each
         *     SCALAR element) must not become Total[v);
         *   - exactly one array leaf, and it is the parameter: a second array
         *     would make Map produce a list OF arrays, not an elementwise result;
         *   - the result is rank 1 with the promised element type.
         * Anything else rolls back and takes the general loop, which is correct
         * for every body — this is a cost split, not a subset restriction. */
        if (!scan && !(c->flags & COMPILE_NO_FUSE)
            && fs.kind == FN_LAMBDA && fs.nparams == 1 && c->nscope < CTX_MAX_SCOPE) {
            EmitMark mk = emit_mark(c);
            c->scope[c->nscope].name = fs.pname[0]; c->scope[c->nscope].reg = va.reg;
            c->scope[c->nscope].type = va.type;     c->scope[c->nscope].built = va.built;
            c->nscope++;
            FuseLeaves L; L.n = 0;
            bool legal = fuse_listable(c, fs.body) == 1
                      && fuse_collect(c, fs.body, &L) && L.n == 1
                      && L.name[0] == fs.pname[0];
            Val bv;
            bool got = legal && emit(c, fs.body, &bv) && CT_IS_ARRAY(bv.type)
                    && CT_RANK(bv.type) == 1 && CT_ELEM(bv.type) == rt;
            c->nscope--;
            if (got) {
                int res = bv.reg;
                if (va.tmp && bv.tmp && res != va.reg) {  /* restore LIFO, as Part does */
                    ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
                    ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)res, 0, z);
                    c->arr_top--;
                    res = va.reg;
                }
                c->temp_top = (rn - c->nlocals);          /* rn, ri unused on this route */
                out->reg = res; out->tmp = bv.tmp; out->type = bv.type;
                out->built = va.built;
                return c->ok;
            }
            emit_rollback(c, mk);
        }

        int rout = -1;
        if (!scan) {
            rout = alloc_arr(c);
            ins_f(c, OP_A_NEWLIKE, (uint16_t)(((unsigned)rt & 3u) << AF_R_SHIFT),
                  (uint32_t)rout, (uint32_t)va.reg, 0, z);
        }
        ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int relem = alloc_temp(c);
        int body_top = c->temp_top;
        ins(c,a_load_op(el),
            (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
        Val ev = { relem, false, el, false }, vb;
        if (!emit_apply(c, &fs, &ev, 1, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        if (!scan) {
            coerce(c, &vb, rt); if (!c->ok) return false;
            ins(c,a_store_op(rt),
                (uint32_t)rout, (uint32_t)ri, (uint32_t)vb.reg, z);
        }
        c->temp_top = body_top - 1;                     /* body temps and relem */
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        if (scan) {
            free_if_tmp(c, va);
            c->temp_top = (rn - c->nlocals);
            int r0 = alloc_temp(c); ins(c, OP_CONST, (uint32_t)r0, 0, 0, k0);
            out->reg = r0; out->tmp = true; out->type = CT_INT; out->built = false;
            return c->ok;
        }
        /* The source may itself be a temporary sitting BELOW the result in the
         * array stack, so it cannot just be popped: free it and slide the result
         * down into its slot, restoring LIFO.  Same dance as the Part lowering. */
        if (va.tmp) {
            ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = va.reg;
        }
        c->temp_top = (rn - c->nlocals);                /* rn, ri dead; rout is an arr */
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(rt, 1);
        out->built = va.built;                          /* Map follows its source's kind */
        return c->ok;
    }

    /* Fold[f, x0, v] / Fold[f, v]: reduce a rank-1 array into a scalar
     * accumulator.  Nest's loop with an index and a per-iteration element. */
    if (strcmp(h, "Fold") == 0 && (na == 2 || na == 3)) {
        FnSpec fs;
        if (!fn_resolve(A[0], 2, &fs)) { c->ok = false; return false; }
        CompileType el = vec_elem_type(c, A[na - 1]);
        if ((int)el < 0) { c->ok = false; return false; }
        CompileType t0 = el;
        if (na == 3 && !infer_type(c, A[1], &t0)) { c->ok = false; return false; }
        int tfp = accum_fixed_type(c, &fs, t0, &el, 1);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return false; }
        CompileType T = (CompileType)tfp;
        Slot z = { 0 };

        /* Persistent registers first, so the operand lowerings below allocate
         * ABOVE them and freeing those temps cannot reach them. */
        int racc = alloc_temp(c), rn = alloc_temp(c), ri = alloc_temp(c);
        Val va;
        if (!emit(c, A[na - 1], &va)) return false;
        if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return false; }
        ins(c, OP_A_SIZE, (uint32_t)rn, (uint32_t)va.reg, 0, z);

        uint16_t ldop = a_load_op(el);
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        if (na == 3) {
            Val vx; if (!emit(c, A[1], &vx)) return false;
            coerce(c, &vx, T); if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);
        } else {
            /* Fold[f, {}] stays UNEVALUATED (src/funcprog.c:2282), so an empty
             * vector must fail the call rather than answer with a seed. */
            emit_nonzero_guard(c, rn);
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k0);
            int rs = alloc_temp(c);
            ins(c, ldop, (uint32_t)rs, (uint32_t)va.reg, (uint32_t)ri, z);
            Val sv = { rs, true, el, false };
            coerce(c, &sv, T); if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)sv.reg, 0, z);
            free_if_tmp(c, sv);
            Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
            ins(c, OP_CONST, (uint32_t)ri, 0, 0, k1);   /* the seed is consumed */
        }

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)ri, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;                                  /* guard temp is dead */

        int relem = alloc_temp(c);
        int body_top = c->temp_top;
        ins(c, ldop, (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ri, z);
        Val argv[2] = { { racc, false, T, false }, { relem, false, el, false } }, vb;
        if (!emit_apply(c, &fs, argv, 2, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, T); if (!c->ok) return false;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = body_top - 1;                     /* drop the body's temps and relem */
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)ri, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        free_if_tmp(c, va);                             /* release the source buffer */
        c->temp_top = (racc - c->nlocals) + 1;          /* keep racc */
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok;
    }

    /* NestList[f, x, n] / FoldList[f, x0, v] / FoldList[f, v]: the same
     * accumulator loops as Nest and Fold, writing every iterate into a buffer
     * whose length is known before the loop starts (n + 1, or the source length
     * plus one for the seed).
     *
     * The element type must be Real or Complex — these BUILD their result, so
     * the ConstantArray rule applies: `NestList[2 # &, 1, 5]` holds exact
     * Integers in the interpreter and a packed buffer has no integer dtype. */
    if ((strcmp(h, "NestList") == 0 && na == 3)
        || (strcmp(h, "FoldList") == 0 && (na == 2 || na == 3))) {
        bool nest = h[0] == 'N';
        FnSpec fs;
        if (!fn_resolve(A[0], nest ? 1 : 2, &fs)) { c->ok = false; return false; }

        CompileType el = CT_ERR, T;
        if (nest) {
            CompileType tn;
            if (!infer_type(c, A[2], &tn) || tn != CT_INT) { c->ok = false; return false; }
            CompileType tx;
            if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
            int tfp = nest_fixed_type(c, &fs, tx);
            if (tfp < 0) { c->ok = false; return false; }
            T = (CompileType)tfp;
        } else {
            el = vec_elem_type(c, A[na - 1]);
            if ((int)el < 0) { c->ok = false; return false; }
            CompileType t0 = el;
            if (na == 3 && !infer_type(c, A[1], &t0)) { c->ok = false; return false; }
            int tfp = accum_fixed_type(c, &fs, t0, &el, 1);
            if (tfp < 0) { c->ok = false; return false; }
            T = (CompileType)tfp;
        }
        if (!ct_is_elem(T)) { c->ok = false; return false; }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
        uint16_t stop = a_store_op(T);

        /* Persistent registers first; `rlen` must be the sole dimension operand
         * of A_NEW, which reads a contiguous run starting there. */
        int rlen = alloc_temp(c), racc = alloc_temp(c);
        int rn = alloc_temp(c), rcnt = alloc_temp(c), rone = alloc_temp(c);
        ins(c, OP_CONST, (uint32_t)rone, 0, 0, k1);

        Val va = { 0, false, CT_REAL, false };
        if (nest) {
            Val vn; if (!emit(c, A[2], &vn)) return false;
            if (vn.type != CT_INT) { c->ok = false; return false; }
            ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            emit_nonneg_guard(c, rn);                    /* n < 0 is unevaluated */
            ins(c, OP_ADD_I, (uint32_t)rlen, (uint32_t)rn, (uint32_t)rone, z);
        } else {
            if (!emit(c, A[na - 1], &va)) return false;
            if (!CT_IS_ARRAY(va.type) || CT_RANK(va.type) != 1) { c->ok = false; return false; }
            ins(c, OP_A_SIZE, (uint32_t)rlen, (uint32_t)va.reg, 0, z);
            if (na == 3) {
                /* seeded: the history is the seed plus one entry per element */
                ins(c, OP_MOVE, (uint32_t)rn, (uint32_t)rlen, 0, z);
                ins(c, OP_ADD_I, (uint32_t)rlen, (uint32_t)rlen, (uint32_t)rone, z);
            } else {
                /* seedless: the first element IS the seed, so one fewer step.
                 * FoldList[f, {}] is `{}` in the interpreter, and an empty
                 * packed result is not the same value, so decline it. */
                emit_nonzero_guard(c, rlen);
                ins(c, OP_SUB_I, (uint32_t)rn, (uint32_t)rlen, (uint32_t)rone, z);
            }
        }
        emit_max_guard(c, rlen, 1 << 26);   /* bound the allocation, not the answer */

        Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)T;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rlen, (uint32_t)rlen, el_i);

        /* seed into racc, and into slot 0 of the buffer */
        if (nest) {
            Val vx; if (!emit(c, A[1], &vx)) return false;
            coerce(c, &vx, T); if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        } else if (na == 3) {
            Val vx; if (!emit(c, A[1], &vx)) return false;
            coerce(c, &vx, T); if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        } else {
            int rs = alloc_temp(c);
            ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
            ins(c,a_load_op(el),
                (uint32_t)rs, (uint32_t)va.reg, (uint32_t)rcnt, z);
            Val sv = { rs, true, el, false };
            coerce(c, &sv, T); if (!c->ok) return false;
            ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)sv.reg, 0, z);
            free_if_tmp(c, sv);
        }
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
        ins(c, stop, (uint32_t)rout, (uint32_t)rcnt, (uint32_t)racc, z);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rn, z);
        size_t jz = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int body_top = c->temp_top, relem = -1;
        Val argv[2] = { { racc, false, T, false }, { 0, false, CT_REAL, false } }, vb;
        int nargs = 1;
        if (!nest) {
            /* Element k of the source pairs with history slot k+1: seeded folds
             * consume v[[k+1]] at step k, seedless ones v[[k+2]] (the first
             * element was the seed). */
            relem = alloc_temp(c);
            body_top = c->temp_top;
            int ridx = alloc_temp(c);
            if (na == 3) ins(c, OP_MOVE, (uint32_t)ridx, (uint32_t)rcnt, 0, z);
            else         ins(c, OP_ADD_I, (uint32_t)ridx, (uint32_t)rcnt, (uint32_t)rone, z);
            ins(c,a_load_op(el),
                (uint32_t)relem, (uint32_t)va.reg, (uint32_t)ridx, z);
            c->temp_top--;                                /* ridx dead after the load */
            argv[1].reg = relem; argv[1].type = el;
            nargs = 2;
        }
        if (!emit_apply(c, &fs, argv, nargs, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, T); if (!c->ok) return false;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = nest ? body_top : body_top - 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, k1);
        ins(c, stop, (uint32_t)rout, (uint32_t)rcnt, (uint32_t)racc, z);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) c->code[jz].b = (uint32_t)c->n;

        if (!nest && va.tmp) {           /* restore LIFO, as the Part lowering does */
            ins(c, OP_ARR_FREE, (uint32_t)va.reg, 0, 0, z);
            ins(c, OP_A_XFER, (uint32_t)va.reg, (uint32_t)rout, 0, z);
            c->arr_top--;
            rout = va.reg;
        }
        c->temp_top = (rlen - c->nlocals);
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(T, 1);
        /* NestList always constructs; FoldList's history is packed by the
         * interpreter with the SOURCE dtype, so it follows the source's kind. */
        out->built = nest ? true : va.built;
        return c->ok;
    }

    /* FixedPointList[f, x, ...] / NestWhileList[f, x, test]: the same loops as
     * their scalar twins, keeping every iterate.
     *
     * The length is not known until the loop has run, so the buffer GROWS
     * (A_PUSH) and is cut to size at the end (A_TRUNC).  Allocating the safety
     * cap up front instead would be 8 MB per call, and running the body twice
     * to count first would double every side effect a `Set` in it performs.
     * The capacity a particular run reached is never observable. */
    if ((strcmp(h, "FixedPointList") == 0 && na >= 2 && na <= 4)
        || (strcmp(h, "NestWhileList") == 0 && na == 3)) {
        bool fp = h[0] == 'F';
        FnSpec fs, ss, ts;
        const Expr *mx = NULL, *st = NULL;
        if (!fn_resolve(A[0], 1, &fs)) { c->ok = false; return false; }
        if (fp) { if (!fp_opts(A, na, 2, &mx, &st)) { c->ok = false; return false; } }
        else if (!fn_resolve(A[2], 1, &ts)) { c->ok = false; return false; }

        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0) { c->ok = false; return false; }
        CompileType T = (CompileType)tfp;
        /* A built history, so the ConstantArray element-type rule applies. */
        if (!ct_is_elem(T)) { c->ok = false; return false; }
        if (fp && st && !fn_resolve(st, 2, &ss)) { c->ok = false; return false; }
        if (!fp) {   /* the while-test must yield a Bool on the accumulator */
            CompileType tb;
            if (!infer_apply(c, &ts, &T, 1, &tb) || tb != CT_BOOL) { c->ok = false; return false; }
        }

        Slot z = { 0 }, k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        Slot k1; memset(&k1, 0, sizeof k1); k1.i = 1;
        uint16_t pflags = (uint16_t)(((unsigned)T & 3u) << AF_R_SHIFT);

        int rcap = alloc_temp(c), racc = alloc_temp(c);
        int rk = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Slot cap0; memset(&cap0, 0, sizeof cap0); cap0.i = 16;   /* grown as needed */
        ins(c, OP_CONST, (uint32_t)rcap, 0, 0, cap0);

        Val vx; if (!emit(c, A[1], &vx)) return false;
        coerce(c, &vx, T); if (!c->ok) return false;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        if (mx) {
            Val vn; if (!emit(c, mx, &vn)) return false;
            if (vn.type != CT_INT) { c->ok = false; return false; }
            ins(c, OP_MOVE, (uint32_t)rlim, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            emit_nonneg_guard(c, rlim);
        } else {
            Slot cp; memset(&cp, 0, sizeof cp); cp.i = VM_ITER_SAFETY_CAP;
            ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cp);
        }

        Slot el_i; memset(&el_i, 0, sizeof el_i); el_i.i = (long long)T;
        int rout = alloc_arr(c);
        ins_f(c, OP_A_NEW, 1, (uint32_t)rout, (uint32_t)rcap, (uint32_t)rcap, el_i);

        ins(c, OP_CONST, (uint32_t)rk, 0, 0, k0);
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);
        ins_f(c, OP_A_PUSH, pflags, (uint32_t)rout, (uint32_t)rk, (uint32_t)racc, z);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);

        size_t Lp = c->n, jend = 0, jcap = 0;
        int body_top = c->temp_top;
        if (!fp) {                       /* NestWhileList tests BEFORE applying */
            Val acc0 = { racc, false, T, false }, vt;
            if (!emit_apply(c, &ts, &acc0, 1, &vt)) return false;
            if (vt.type != CT_BOOL) { c->ok = false; return false; }
            jend = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
            c->temp_top = body_top;
        }
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        Val acc = { racc, false, T, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, T); if (!c->ok) return false;
        int rsame = -1;
        if (fp) {                        /* compare before the accumulator moves */
            if (st) {
                Val sargv[2] = { { racc, false, T, false }, vb }, sv;
                sargv[1].tmp = false;
                if (!emit_apply(c, &ss, sargv, 2, &sv)) return false;
                if (sv.type != CT_BOOL) { c->ok = false; return false; }
                rsame = sv.reg;
            } else {
                rsame = alloc_temp(c);
                ins(c, emit_sameq_op(T), (uint32_t)rsame, (uint32_t)vb.reg, (uint32_t)racc, z);
            }
        }
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, k1);
        ins_f(c, OP_A_PUSH, pflags, (uint32_t)rout, (uint32_t)rk, (uint32_t)racc, z);
        ins(c, OP_INC_I, (uint32_t)rk, 0, 0, k1);
        if (fp) ins(c, OP_JZ, 0, (uint32_t)rsame, (uint32_t)Lp, z);  /* not same -> iterate */
        else    ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        c->temp_top = body_top;
        size_t jdone = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jcap].b = (uint32_t)c->n;
        /* Reaching the cap of an UNBOUNDED run is where the interpreter gives
         * up too; with a user bound, falling through keeps what was collected. */
        if (!mx) ins(c, OP_FAIL, 0, 0, 0, z);
        if (c->ok) {
            c->code[jdone].b = (uint32_t)c->n;
            if (!fp) c->code[jend].b = (uint32_t)c->n;
        }
        ins_f(c, OP_A_TRUNC, 0, (uint32_t)rout, (uint32_t)rk, 0, z);

        c->temp_top = (rcap - c->nlocals);
        out->reg = rout; out->tmp = true; out->type = CT_ARRAY(T, 1); out->built = true;
        return c->ok;
    }

    /* FixedPoint[f, x] / FixedPoint[f, x, n] / SameTest -> s: iterate until two
     * successive values are SameQ.  See emit_sameq for why that is not Equal. */
    if (strcmp(h, "FixedPoint") == 0 && na >= 2 && na <= 4) {
        FnSpec fs, ss; const Expr *mx, *st;
        if (!fn_resolve(A[0], 1, &fs) || !fp_opts(A, na, 2, &mx, &st)) { c->ok = false; return false; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return false; }
        CompileType T = (CompileType)tfp;
        if (st && !fn_resolve(st, 2, &ss)) { c->ok = false; return false; }
        if (!st && T == CT_BOOL) { c->ok = false; return false; }   /* no SameQ opcode for Bool */
        Slot z = { 0 };

        int racc = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return false;
        coerce(c, &vx, T); if (!c->ok) return false;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        if (mx) {
            Val vn; if (!emit(c, mx, &vn)) return false;
            if (vn.type != CT_INT) { c->ok = false; return false; }
            ins(c, OP_MOVE, (uint32_t)rlim, (uint32_t)vn.reg, 0, z); free_if_tmp(c, vn);
            /* A negative bound leaves the whole call unevaluated. */
            emit_nonneg_guard(c, rlim);
        } else {
            Slot cap; memset(&cap, 0, sizeof cap); cap.i = VM_ITER_SAFETY_CAP;
            ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cap);
        }
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);

        size_t Lp = c->n;
        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        size_t jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        int body_top = c->temp_top;
        Val acc = { racc, false, T, false }, vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, T); if (!c->ok) return false;
        /* Compare BEFORE the accumulator is overwritten. */
        int rsame;
        if (st) {
            Val sargv[2] = { { racc, false, T, false }, vb }, sv;
            sargv[1].tmp = false;
            if (!emit_apply(c, &ss, sargv, 2, &sv)) return false;
            if (sv.type != CT_BOOL) { c->ok = false; return false; }
            rsame = sv.reg;
        } else {
            rsame = alloc_temp(c);
            ins(c, emit_sameq_op(T), (uint32_t)rsame, (uint32_t)vb.reg, (uint32_t)racc, z);
        }
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JZ, 0, (uint32_t)rsame, (uint32_t)Lp, z);   /* not same -> iterate */
        c->temp_top = body_top;
        size_t jend = c->n; ins(c, OP_JMP, 0, 0, 0, z);
        if (c->ok) c->code[jcap].b = (uint32_t)c->n;
        /* Reaching the cap of an UNBOUNDED run is where the interpreter gives up
         * too, so fail the call rather than answer with a non-fixed point.  With
         * a user bound, falling through returns f^n(x), which is what
         * FixedPoint[f, x, n] means. */
        if (!mx) ins(c, OP_FAIL, 0, 0, 0, z);
        if (c->ok) c->code[jend].b = (uint32_t)c->n;

        c->temp_top = (racc - c->nlocals) + 1;
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok;
    }

    /* NestWhile[f, x, test]: apply f while test holds on the CURRENT value, so
     * the test runs before the first application (src/funcprog.c:2323). */
    if (strcmp(h, "NestWhile") == 0 && na == 3) {
        FnSpec fs, ts;
        if (!fn_resolve(A[0], 1, &fs) || !fn_resolve(A[2], 1, &ts)) { c->ok = false; return false; }
        CompileType tx; if (!infer_type(c, A[1], &tx)) { c->ok = false; return false; }
        int tfp = nest_fixed_type(c, &fs, tx);
        if (tfp < 0 || CT_IS_ARRAY((CompileType)tfp)) { c->ok = false; return false; }
        CompileType T = (CompileType)tfp;
        Slot z = { 0 };

        int racc = alloc_temp(c), rcnt = alloc_temp(c), rlim = alloc_temp(c);
        Val vx; if (!emit(c, A[1], &vx)) return false;
        coerce(c, &vx, T); if (!c->ok) return false;
        ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vx.reg, 0, z); free_if_tmp(c, vx);
        Slot cap; memset(&cap, 0, sizeof cap); cap.i = VM_ITER_SAFETY_CAP;
        ins(c, OP_CONST, (uint32_t)rlim, 0, 0, cap);
        Slot k0; memset(&k0, 0, sizeof k0); k0.i = 0;
        ins(c, OP_CONST, (uint32_t)rcnt, 0, 0, k0);

        size_t Lp = c->n;
        int body_top = c->temp_top;
        Val acc = { racc, false, T, false }, vt;
        if (!emit_apply(c, &ts, &acc, 1, &vt)) return false;
        if (vt.type != CT_BOOL) { c->ok = false; return false; }
        size_t jend = c->n; ins(c, OP_JZ, 0, (uint32_t)vt.reg, 0, z);
        c->temp_top = body_top;

        int rc = alloc_temp(c);
        ins(c, OP_LT_I, (uint32_t)rc, (uint32_t)rcnt, (uint32_t)rlim, z);
        size_t jcap = c->n; ins(c, OP_JZ, 0, (uint32_t)rc, 0, z);
        c->temp_top--;

        Val vb;
        if (!emit_apply(c, &fs, &acc, 1, &vb)) return false;
        if (CT_IS_ARRAY(vb.type)) { c->ok = false; return false; }
        coerce(c, &vb, T); if (!c->ok) return false;
        if (vb.reg != racc) ins(c, OP_MOVE, (uint32_t)racc, (uint32_t)vb.reg, 0, z);
        c->temp_top = body_top;
        Slot one; memset(&one, 0, sizeof one); one.i = 1;
        ins(c, OP_INC_I, (uint32_t)rcnt, 0, 0, one);
        ins(c, OP_JMP, 0, 0, (uint32_t)Lp, z);
        if (c->ok) { c->code[jcap].b = (uint32_t)c->n; }
        ins(c, OP_FAIL, 0, 0, 0, z);                    /* cap: interpreter gives up too */
        if (c->ok) c->code[jend].b = (uint32_t)c->n;

        c->temp_top = (racc - c->nlocals) + 1;
        out->reg = racc; out->tmp = true; out->type = T; out->built = false;
        return c->ok;
    }

    /* last resort: any numeric function with a machine kernel in ndkernels */
    Val kv;
    if (try_kernel(c, h, A, na, &kv)) { *out = kv; return c->ok; }

    /* g[a1,...,aN] where g is a CompiledFunction: inline its body with the
     * parameters bound to registers holding the evaluated arguments.  Each
     * argument is lowered in the CALLER's environment, then the body is lowered
     * with only the parameters in scope (see arg_find). */
    {
        const CompiledFunction* cf = compiled_callee(c, h);

        /* Inline or call?
         *
         * Inlining is the faster of the two — no frame, no argument copy — and
         * stays the default.  But pasting a LARGE callee in at every use site
         * multiplies code size and compile time, and the inliner also declines
         * outright when the depth cap or the scope budget is hit, which used to
         * bail the entire body to the interpreter.  So a big callee is CALLed
         * instead: machine values in registers, its own frame, no Expr and no
         * evaluator round-trip.  A callee with no program of its own (its body
         * did not compile standalone) can only be inlined. */
        const CompiledProgram* cp = cf ? compiled_function_program(cf) : NULL;
        bool prefer_call = cp && compiled_num_instructions(cp) > INLINE_MAX_INSTRS;

        if (cf && !prefer_call && na >= 1 && compiled_function_num_args(cf) == na
            && na + (size_t)c->nscope <= CTX_MAX_SCOPE
            && c->inlining < 8) {          /* depth cap: a self-referential body */
            const char* const* pn = compiled_function_arg_names(cf);
            const CompileType* pt = compiled_function_arg_types(cf);
            Slot z = { 0 };
            /* Sized from the scope bound, which is what the guard above tests:
             * the two must be the same constant or raising one overflows this. */
            int preg[CTX_MAX_SCOPE], res_top = 0;
            for (size_t i = 0; i < na; i++) {
                Val v;
                if (!emit(c, A[i], &v)) return false;
                if (CT_IS_ARRAY(v.type) || CT_IS_ARRAY(pt[i])) { c->ok = false; return false; }
                coerce(c, &v, pt[i]);
                if (!c->ok) return false;
                /* Copy into a dedicated register: the parameter stays live for
                 * the whole body, while v may sit in a caller temp that the body
                 * would otherwise be free to reuse. */
                preg[i] = alloc_temp(c);
                if (i == 0) res_top = c->temp_top;   /* body result lands here */
                ins(c, OP_MOVE, (uint32_t)preg[i], (uint32_t)v.reg, 0, z);
            }
            int saved_scope = c->nscope;
            for (size_t i = 0; i < na; i++) {
                c->scope[c->nscope].name = pn[i];
                c->scope[c->nscope].reg = preg[i];
                c->scope[c->nscope].type = pt[i];
                c->nscope++;
            }
            c->inlining++;
            Val bv; bool okb = emit(c, compiled_function_body(cf), &bv);
            c->inlining--;
            c->nscope = saved_scope;
            if (!okb) return false;
            if (CT_IS_ARRAY(bv.type)) { c->ok = false; return false; }
            /* Land the result in the first parameter's register and pop
             * everything the call allocated above it. */
            if (bv.reg != preg[0]) ins(c, OP_MOVE, (uint32_t)preg[0], (uint32_t)bv.reg, 0, z);
            c->temp_top = res_top;
            out->reg = preg[0]; out->tmp = true; out->type = bv.type;
            return c->ok;
        }

        /* Inlining declined — the depth cap, too many parameters, or an array
         * signature.  Emit a real CALL rather than bailing: the callee runs on
         * its own frame with machine values passed in registers, so a deep or
         * self-referential chain compiles instead of dropping the whole body to
         * the interpreter.  Arguments go in consecutive registers because the
         * instruction carries only their base and count. */
        if (cf && cp && na >= 1 && compiled_function_num_args(cf) == na && na <= 16
            && !c->vector_mode) {
            const CompileType* pt = compiled_function_arg_types(cf);
            bool ok_args = true;
            for (size_t i = 0; i < na && ok_args; i++)
                if (CT_IS_ARRAY(pt[i])) ok_args = false;
            if (ok_args && !CT_IS_ARRAY(compiled_result_type(cp))) {
                EmitMark mk = emit_mark(c);
                Slot z; memset(&z, 0, sizeof z);
                /* Reserve the argument block FIRST so it is contiguous: each
                 * argument's own lowering allocates temps ABOVE it, which would
                 * otherwise interleave with the parameter slots. */
                int base = alloc_temp(c);
                for (size_t i = 1; i < na; i++) (void)alloc_temp(c);
                int after_params = c->temp_top;
                for (size_t i = 0; i < na && c->ok; i++) {
                    Val v;
                    if (!emit(c, A[i], &v) || CT_IS_ARRAY(v.type)) { c->ok = false; break; }
                    coerce(c, &v, pt[i]);
                    if (!c->ok) break;
                    ins(c, OP_MOVE, (uint32_t)(base + (int)i), (uint32_t)v.reg, 0, z);
                    c->temp_top = after_params;         /* drop this argument's temps */
                }
                if (c->ok) {
                    Slot k; memset(&k, 0, sizeof k); k.p = cp;
                    c->temp_top = (base - c->nlocals) + 1;   /* result lands in base */
                    ins_f(c, OP_CALL, (uint16_t)na, (uint32_t)base, (uint32_t)base, 0, k);
                    out->reg = base; out->tmp = true;
                    out->type = compiled_result_type(cp);
                    return c->ok;
                }
                emit_rollback(c, mk);
            }
        }
    }

    c->ok = false; return false;   /* unsupported head -> bail */
}

/* Lower `e`, recording the innermost failure.
 *
 * Every recursive descent in emit_node goes through this wrapper, so on the way
 * back up the first writer is the deepest subexpression that could not be
 * lowered — which is the one a user needs to see.  Doing it in one place rather
 * than at each of the ~40 bail sites means a new bail is diagnosed the day it is
 * written, with no chance of a site being forgotten. */
static bool emit(Ctx* c, const Expr* e, Val* out) {
    bool r = emit_node(c, e, out);
    if (!r && !c->bail_node) c->bail_node = e;
    return r;
}

/* ------------------------------------------------------------------ *
 *  VM                                                                 *
 * ------------------------------------------------------------------ */
/* The integer form is `ci_powi` in compile_internal.h: same binary exponentiation,
 * overflow-checked at every step, and shared with the optimiser's folding. */
static double    ipow_r(double b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }

/* ------------------------------------------------------------------ *
 *  Exact integer kernels for the integer-closed heads                 *
 * ------------------------------------------------------------------ *
 * Each returns TRUE ON FAILURE — overflow, or an argument outside the domain
 * where the interpreter yields an integer — matching the ci_* convention, so a
 * VM body is the same one-line IOP() as an arithmetic opcode.  Failure hands the
 * call to the interpreter, which then answers exactly (a bigint) or symbolically
 * (ComplexInfinity, a Rational), whichever is right.
 *
 * The domains are not guesses; they are what the interpreter actually does:
 *   Factorial[-1]     ComplexInfinity     Gamma[0]        ComplexInfinity
 *   Pochhammer[7,-3]  1/120  (Rational)   7^-3            1/343  (Rational)
 *   Arg[-3]           Pi     (Symbol)     Binomial[3,7]   0      (all integer)
 * so Factorial needs n >= 0, Gamma n >= 1, Pochhammer n >= 0, Power e >= 0, and
 * Binomial no guard at all. */

static bool int_factorial(long long n, long long* out) {
    if (n < 0) return true;                       /* ComplexInfinity */
    long long r = 1;
    for (long long k = 2; k <= n; k++) if (ci_mul(r, k, &r)) return true;
    *out = r;
    return false;                                 /* 21! already overflows */
}

/* Gamma[n] = (n-1)! on the positive integers; Gamma[0] and Gamma of a negative
 * integer are ComplexInfinity, which is not a machine number. */
static bool int_gamma(long long n, long long* out) {
    if (n < 1) return true;
    return int_factorial(n - 1, out);
}

/* Binomial[n, k] for machine integers, by the multiplicative recurrence
 * C(n,k) = C(n,k-1) * (n-k+1) / k, which is exact at every step because the
 * partial product of k consecutive integers is divisible by k!.  Computing
 * n!/(k!(n-k)!) directly would overflow at n = 21 for results that fit easily. */
static bool int_binomial(long long n, long long k, long long* out) {
    if (n >= 0) {
        if (k < 0 || k > n) { *out = 0; return false; }       /* Binomial[3,7] = 0 */
        if (k > n - k) k = n - k;                             /* symmetry: fewer steps */
        long long r = 1;
        for (long long i = 1; i <= k; i++) {
            if (ci_mul(r, n - k + i, &r)) return true;
            r /= i;                                           /* exact by construction */
        }
        *out = r;
        return false;
    }
    /* Negative upper index: Binomial[-n, k] = (-1)^k Binomial[n+k-1, k].  Left to
     * the interpreter — it is rare, and the identity is easy to get subtly wrong
     * in a way no test here would catch. */
    return true;
}

/* Pochhammer[a, n] = a (a+1) ... (a+n-1); a falling/negative n gives a Rational. */
static bool int_pochhammer(long long a, long long n, long long* out) {
    if (n < 0) return true;
    long long r = 1;
    for (long long i = 0; i < n; i++) {
        long long t;
        if (ci_add(a, i, &t) || ci_mul(r, t, &r)) return true;
    }
    *out = r;
    return false;
}

/* Fibonacci / LucasL by iteration, including the negative index (the interpreter
 * defines both there: F[-n] = (-1)^(n+1) F[n], L[-n] = (-1)^n L[n]).  Iterative
 * rather than the closed form because these must be EXACT — a double loses
 * Fibonacci exactly where it starts to matter, and F[92] is the last one that
 * fits in an int64 anyway. */
static bool int_fib2(long long n, bool lucas, long long* out) {
    if (n == LLONG_MIN) return true;              /* -n would overflow; test FIRST */
    bool neg = n < 0;
    long long m = neg ? -n : n;
    long long a = lucas ? 2 : 0, b = 1;           /* (L0,L1) = (2,1); (F0,F1) = (0,1) */
    for (long long i = 0; i < m; i++) {
        long long t;
        if (ci_add(a, b, &t)) return true;
        a = b; b = t;
    }
    /* F[-n] flips sign for even n, L[-n] for odd n. */
    if (neg && ((m % 2 == 0) == !lucas)) { if (ci_neg(a, &a)) return true; }
    *out = a;
    return false;
}
static bool int_fib(long long n, long long* out)   { return int_fib2(n, false, out); }
static bool int_lucas(long long n, long long* out) { return int_fib2(n, true,  out); }

/* a^e on the integers.  A negative exponent is a Rational in the interpreter
 * (7^-3 is 1/343) and 0^0 is Indeterminate, so both abandon the call. */
static bool int_pow(long long a, long long e, long long* out) {
    if (e < 0) return true;
    if (e == 0 && a == 0) return true;            /* Indeterminate */
    return ci_powi(a, e, out);
}

/* GCD is non-negative and GCD[0, 0] is 0, matching the interpreter.  The only
 * failure is INT64_MIN, whose magnitude is not representable. */
static bool int_gcd(long long a, long long b, long long* out) {
    if (a == LLONG_MIN || b == LLONG_MIN) return true;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long long t = a % b; a = b; b = t; }
    *out = a;
    return false;
}

/* LCM[a, b] = |a b| / gcd, non-negative, and 0 when either side is 0.  Divide
 * BEFORE multiplying so the intermediate stays in range whenever the answer
 * does. */
static bool int_lcm(long long a, long long b, long long* out) {
    if (a == 0 || b == 0) { *out = 0; return false; }
    long long g;
    if (int_gcd(a, b, &g) || g == 0) return true;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    return ci_mul(a / g, b, out);
}

/* Number of digits of n in the given base; IntegerLength[0] is 0. */
static bool int_ilen(long long n, long long base, long long* out) {
    if (base < 2) return true;
    if (n == LLONG_MIN) return true;
    if (n < 0) n = -n;
    long long k = 0;
    while (n > 0) { n /= base; k++; }
    *out = k;
    return false;
}

/* Largest e with base^e dividing n.  IntegerExponent[0, b] is Infinity in the
 * interpreter, which is not a machine integer. */
static bool int_iexp(long long n, long long base, long long* out) {
    if (base < 2 || n == 0) return true;
    if (n == LLONG_MIN) return true;
    if (n < 0) n = -n;
    long long k = 0;
    while (n % base == 0) { n /= base; k++; }
    *out = k;
    return false;
}

/* PowerMod[a, e, m], including the negative exponent (a modular INVERSE, which
 * exists only when gcd(a, m) is 1 — the interpreter reports that, so we defer).
 *
 * The modulus is capped at sqrt(INT64_MAX) so every intermediate product fits an
 * int64.  Going wider needs a 128-bit multiply, and `__int128` is a GNU
 * extension this file cannot use (CLAUDE.md: strict C99); a larger modulus
 * therefore hands the call to the interpreter, which has GMP. */
#define POWMOD_MAX_M 3037000499LL             /* floor(sqrt(2^63 - 1)) */
static bool int_powmod(long long a, long long e, long long m, long long* out) {
    if (m == 0 || m > POWMOD_MAX_M || m < -POWMOD_MAX_M) return true;
    if (m < 0) return true;                   /* interpreter's sign convention */
    if (m == 1) { *out = 0; return false; }
    a %= m; if (a < 0) a += m;
    if (e < 0) {
        /* Modular inverse by the extended Euclid, then the positive power. */
        long long r0 = m, r1 = a, s0 = 0, s1 = 1;
        while (r1) {
            long long q = r0 / r1;
            long long t = r0 - q * r1; r0 = r1; r1 = t;
            t = s0 - q * s1;           s0 = s1; s1 = t;
        }
        if (r0 != 1) return true;             /* not invertible */
        a = s0 % m; if (a < 0) a += m;
        if (e == LLONG_MIN) return true;      /* -e would overflow */
        e = -e;
    }
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = (r * a) % m;
        e >>= 1;
        if (e > 0) a = (a * a) % m;
    }
    *out = r;
    return false;
}
static double _Complex ipow_c(double _Complex b, long long n) { if (n < 0) { b = 1.0 / b; n = -n; } double _Complex r = 1; while (n) { if (n & 1) r *= b; b *= b; n >>= 1; } return r; }

/* ------------------------------------------------------------------ *
 *  Array opcodes (M3a): delegation to the NDArray subsystem           *
 * ------------------------------------------------------------------ */

/* Read a scalar register as a (re, im) pair, per its compile-time operand kind. */
static void vm_scalar_pair(const Slot* s, unsigned kind, double* re, double* im) {
    if      (kind == AK_COMPLEX) { *re = creal(s->z);   *im = cimag(s->z); }
    else if (kind == AK_INT)     { *re = (double)s->i;  *im = 0.0; }
    else                         { *re = s->r;          *im = 0.0; }
}

/* Box a scalar register as a temporary numeric Expr, because the ND helpers
 * take Expr operands.  Two allocations amortised over a whole-buffer pass is
 * nothing, and it keeps one implementation of the broadcast and dtype-promotion
 * rules instead of a second copy here. */
static Expr* vm_box_scalar(const Slot* s, unsigned kind) {
    /* An integer register boxes as an Integer, so the ND layer sees an EXACT
     * scalar and keeps an int64 buffer exact instead of promoting it. Reading it
     * back through the double pair would defeat the whole point of AK_INT. */
    if (kind == AK_INT) return expr_new_integer((int64_t)s->i);
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
/* An independent copy of `x` in the program's CANONICAL dtype for `relem`.
 *
 * Every array a program owns is float64 or complex64, never float32, and that
 * is load-bearing rather than tidy: A_STORE writes the buffer at its declared
 * width, so one float32 array reaching a store would write doubles into half-
 * sized slots.  Arguments may be any dtype — they are only ever read — so the
 * narrowing is done here, at the point ownership begins. */
/* The buffer dtype a program uses for a given element type.  A program owns only
 * these three — float32 never appears, see nd_own_copy. */
static NDType ct_elem_ndt(unsigned relem) {
    return relem == (unsigned)CT_COMPLEX ? NDT_COMPLEX64
         : relem == (unsigned)CT_INT     ? NDT_INT64
                                         : NDT_FLOAT64;
}

static Expr* nd_own_copy(const Expr* x, unsigned relem) {
    if (!x || x->type != EXPR_NDARRAY) return NULL;
    NDType dt = ct_elem_ndt(relem);
    NDType sdt = x->data.ndarray.dtype;
    /* An integer-typed program will not silently take a float buffer: the two
     * are different element types to the interpreter as well (`Total` of a
     * float64 NDArray is a Real), so the call goes back rather than rounding. */
    if ((dt == NDT_INT64) != (sdt == NDT_INT64)) return NULL;
    size_t n = ndarray_size(x), esz = ndt_elem_size(dt);
    void* buf = malloc(esz * (n ? n : 1));
    if (!buf) return NULL;
    if (sdt == dt) memcpy(buf, x->data.ndarray.data, esz * n);
    else for (size_t k = 0; k < n; k++) {
        double re, im;
        ndt_get(x->data.ndarray.data, k, sdt, &re, &im);
        /* A complex source into a real program is the array form of the scalar
         * contract: the program promised real, so it fails rather than truncate. */
        if (im != 0.0 && dt != NDT_COMPLEX64) { free(buf); return NULL; }
        ndt_set(buf, k, dt, re, im);
    }
    Expr* nw = expr_new_ndarray_like(x, x->data.ndarray.rank, x->data.ndarray.dims, buf, dt);
    if (!nw) free(buf);
    return nw;
}

/* Materialise a general Part's subscript list for this call: literal specs are
 * borrowed straight from the PartSpec, computed ones are boxed from registers. */
static bool ps_build_indices(const PartSpec* ps, const Slot* R, Expr** idx) {
    for (int i = 0; i < ps->n; i++) {
        if (ps->lit[i]) { idx[i] = ps->lit[i]; continue; }
        idx[i] = expr_new_integer((int64_t)R[ps->reg[i]].i);
        if (!idx[i]) {
            for (int j = 0; j < i; j++) if (!ps->lit[j]) expr_free(idx[j]);
            return false;
        }
    }
    return true;
}
static void ps_free_indices(const PartSpec* ps, Expr** idx) {
    for (int i = 0; i < ps->n; i++) if (!ps->lit[i]) expr_free(idx[i]);
}

/* The array opcodes whose operands are a RUN of registers, so they need the
 * whole frame rather than three slots.  Same abort contract as vm_array_op. */
static bool vm_range_array_op(const Instr* c, Slot* R) {
    switch (c->op) {
        case OP_A_NEW: {                  /* ConstantArray[0, {d1, ..., dr}] */
            int rank = (int)c->flags;
            if (rank < 1 || rank > NDARRAY_MAX_RANK) return false;
            int64_t dims[NDARRAY_MAX_RANK];
            size_t n = 1;
            for (int i = 0; i < rank; i++) {
                long long d = R[c->a + (unsigned)i].i;
                if (d < 0) return false;          /* Table[..., {n}] with n < 0 */
                dims[i] = (int64_t)d;
                n *= (size_t)d;
            }
            NDType dt = ct_elem_ndt((unsigned)c->imm.i);
            void* buf = calloc(n ? n : 1, ndt_elem_size(dt));
            if (!buf) return false;
            Expr* nw = expr_new_ndarray_raw(rank, dims, buf, dt);
            if (!nw) { free(buf); return false; }
            expr_free(R[c->dst].arr);             /* register reused from a prior call */
            R[c->dst].arr = nw;
            return true;
        }

        case OP_A_PART: {                 /* Span / All / list / partial indexing */
            const PartSpec* ps = (const PartSpec*)c->imm.p;
            Expr* idx[NDARRAY_MAX_RANK];
            const Expr* src = R[c->b].arr;
            if (!src || src->type != EXPR_NDARRAY) return false;
            if (!ps_build_indices(ps, R, idx)) return false;
            bool degrade = false;
            Expr* r = ndarray_part(src, idx, (size_t)ps->n, &degrade);
            ps_free_indices(ps, idx);
            /* A spec ndarray_part cannot do natively (degrade) or an out-of-range
             * subscript (r == NULL) both abort: the interpreter re-runs the body
             * and produces whatever Part[] properly produces, including a
             * diagnostic.  The compiled path never invents an answer. */
            if (!r) return false;
            if (r->type != EXPR_NDARRAY) { expr_free(r); return false; }
            expr_free(R[c->dst].arr);
            R[c->dst].arr = r;
            return true;
        }

        case OP_A_NDFN: {                 /* Reverse / Sort / Flatten / Take / ... */
            const NdFnSpec* fn = (const NdFnSpec*)c->imm.p;
            Expr* src = R[c->b].arr;
            if (!src || src->type != EXPR_NDARRAY) return false;
            size_t nx = (size_t)c->flags;
            /* The entry point takes the whole CALL, so rebuild it.  expr_copy is
             * a refcount bump (src/expr.c), not a buffer copy, so this costs a
             * node — and the array is never mutated: these paths all allocate a
             * fresh result. */
            Expr* args[1 + NDARRAY_MAX_RANK];
            args[0] = expr_copy(src);
            for (size_t i = 0; i < nx; i++)
                args[1 + i] = expr_new_integer((int64_t)R[c->a + (unsigned)i].i);
            Expr* call = expr_new_function(expr_new_symbol(fn->head), args, 1 + nx);
            if (!call) { expr_free(args[0]); return false; }
            Expr* r = fn->fn(call);
            expr_free(call);
            /* A case the fast path does not handle comes back as a nested List
             * (ndarray_delist_and_reeval), which is precisely the signal to
             * decline: the interpreter then answers, exactly as it would have. */
            if (!r) return false;
            if (r->type != EXPR_NDARRAY) { expr_free(r); return false; }
            expr_free(R[c->dst].arr);
            R[c->dst].arr = r;
            return true;
        }

        case OP_A_PARTSET: {              /* u[[spec...]] = rhs, in place */
            const PartSpec* ps = (const PartSpec*)c->imm.p;
            Expr* idx[NDARRAY_MAX_RANK];
            Expr* tgt = R[c->dst].arr;
            if (!tgt || tgt->type != EXPR_NDARRAY) return false;
            if (!ps_build_indices(ps, R, idx)) return false;
            Expr* rhs = (ps->rhs_kind == AK_ARR)
                      ? R[c->b].arr
                      : vm_box_scalar(&R[c->b], (unsigned)ps->rhs_kind);
            bool ok = rhs && ndarray_part_set(tgt, idx, (size_t)ps->n, rhs);
            if (ps->rhs_kind != AK_ARR) expr_free(rhs);
            ps_free_indices(ps, idx);
            return ok;
        }

        default: return false;
    }
}

static bool vm_array_op(const Instr* c, Slot* d, Slot* a, Slot* b) {
    const unsigned f = c->flags, ka = AF_A(f), kb = AF_B(f);
    Expr* r = NULL;

    switch (c->op) {
        case OP_ARR_FREE:
            expr_free(d->arr); d->arr = NULL;
            return true;

        /* ---- ownership (M3c) ------------------------------------------ */
        case OP_A_COPY: {                 /* a local initialised from a borrowed array */
            Expr* nw = nd_own_copy(a->arr, AF_R(f));
            if (!nw) return false;
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            expr_free(d->arr);
            d->arr = nw;
            return true;
        }

        case OP_A_XFER:                   /* move the handle, do not duplicate it */
            if (d == a) return true;
            expr_free(d->arr);
            d->arr = a->arr;
            a->arr = NULL;                /* or teardown would free it twice */
            return true;

        /* ---- fused-loop setup (M3b) -----------------------------------
         * Out of line because each runs ONCE per call, outside the element
         * loop; only A_LOAD/A_STORE are inline in the dispatch loop. */
        case OP_A_SIZE: {                 /* total elements, any rank */
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            d->i = (long long)ndarray_size(x);
            return true;
        }

        case OP_A_SHAPECHK: {             /* all leaves of a fused loop agree */
            const Expr* x = a->arr;
            const Expr* y = b->arr;
            if (!x || !y || x->type != EXPR_NDARRAY || y->type != EXPR_NDARRAY) return false;
            if (x->data.ndarray.rank != y->data.ndarray.rank) return false;
            for (int i = 0; i < x->data.ndarray.rank; i++)
                if (x->data.ndarray.dims[i] != y->data.ndarray.dims[i]) return false;
            return true;
        }

        case OP_A_NEWLIKE: {              /* result buffer, shape of the operand */
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            NDType dt = ct_elem_ndt(AF_R(f));
            size_t nelem = ndarray_size(x);
            size_t esz = ndt_elem_size(dt);
            void* buf = calloc(nelem ? nelem : 1, esz);
            if (!buf) return false;
            Expr* nw = expr_new_ndarray_like(x, x->data.ndarray.rank, x->data.ndarray.dims, buf, dt);
            if (!nw) { free(buf); return false; }
            expr_free(d->arr);            /* reused register from a prior call */
            d->arr = nw;
            return true;
        }

        case OP_V_LEN: {
            const Expr* x = a->arr;
            if (!x || x->type != EXPR_NDARRAY) return false;
            long long len = (long long)x->data.ndarray.dims[0];
            if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
            d->i = len;
            return true;
        }

        case OP_V_TOTAL: {
            /* An INTEGER array sums here rather than through ndred_total_all,
             * which accumulates in a double and is therefore exact only to 2^53
             * — `Total[{9007199254740993, 1}]` came back one short.  Summing in
             * int64 with the same overflow rule as the scalar opcodes keeps the
             * answer identical to the interpreter's, or hands the call back. */
            if (a->arr && a->arr->type == EXPR_NDARRAY
                && a->arr->data.ndarray.dtype == NDT_INT64) {
                const NDArrayData* A_ = &a->arr->data.ndarray;
                const int64_t* p = (const int64_t*)A_->data;
                size_t n = ndarray_size(a->arr);
                long long acc = 0;
                bool ovf = false;
                for (size_t k = 0; k < n && !ovf; k++) ovf = ci_add(acc, (long long)p[k], &acc);
                if (f & AF_FREE_A) { expr_free(a->arr); a->arr = NULL; }
                if (ovf) return false;
                d->i = acc;
                return true;
            }
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
/* OPLIST now lives in compile_internal.h: one list drives the opcode enum, the
 * VM jump table below, and the optimiser's instruction-property table. */

/* The bytecode interpreter.  A threaded (computed-goto) dispatch is used on
 * GCC/Clang — each opcode ends by jumping straight to the next, which the branch
 * predictor handles far better than a single switch; a portable switch is the
 * fallback.  Programs always end in OP_RET, so no per-op bounds check is needed. */
#if defined(__GNUC__) && !defined(VM_NO_THREADED)
#define VM_THREADED 1
#else
#define VM_THREADED 0
#endif

static bool vm_call(const CompiledProgram* cp, const Slot* argv, unsigned nargs, Slot* dst);
static void vm_run(const Instr* code, size_t n, Slot* R, bool* failed);

/* ------------------------------------------------------------------ *
 *  Parallel strip loop (OP_APAR)                                      *
 * ------------------------------------------------------------------ *
 * A fused MAP is embarrassingly parallel: element i of the output depends only
 * on element i of the inputs, so a worker can be handed a sub-range of the flat
 * index space and the result is BIT-IDENTICAL to the serial pass — same
 * operations, same order, on the same values. Nothing here changes an answer;
 * it only changes which core computes it.
 *
 * Each worker gets its OWN frame, seeded by copying the parent's registers.
 * That copy is what makes this safe with no locking anywhere:
 *   - scalars and loop bounds are inherited by value;
 *   - array registers are inherited as raw pointers, which is correct because
 *     the workers only READ the inputs and write DISJOINT output elements;
 *   - tile registers are re-pointed at the worker's own tile storage, so no two
 *     workers ever touch the same VBLOCK scratch buffer.
 * No worker allocates or frees an array (the output buffer was allocated by the
 * parent before the loop), so no ownership crosses a thread boundary.
 *
 * Returns false to DECLINE — too small to be worth threading, out of memory, or
 * a worker failed. Declining is always safe: OP_APAR then falls through into the
 * serial loop, which is still sitting immediately after it and which recomputes
 * the same map from the same starting index. */
#define VM_STACK_SLOTS 512

typedef struct { const ParLoop* pl; const Slot* parent; } ParCtx;

static bool par_chunk(void* vctx, size_t lo, size_t hi) {
    const ParCtx* x = (const ParCtx*)vctx;
    const ParLoop* pl = x->pl;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (pl->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(pl->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    memcpy(R, x->parent, (size_t)pl->nreg * sizeof(Slot));
    Slot* tiles = R + pl->nreg;
    for (int k = 0; k < pl->ntiles; k++)
        R[pl->tile_base + k].p = tiles + (size_t)k * VBLOCK;

    R[pl->ri].i = (long long)lo;
    R[pl->rn].i = (long long)hi;


    bool failed = false;
    vm_run(pl->code, pl->n, R, &failed);
    free(heap);
    return !failed;
}

static bool vm_par_run(const ParLoop* pl, Slot* R) {
#ifdef MATHILDA_THREADS
    long long n = R[pl->rn].i;
    /* nd_parallel_for owns the "is this worth threading" decision and the
     * chunking, so the compiled path fans out on exactly the same terms as the
     * rest of the ND layer rather than inventing a second policy. It runs the
     * range serially below its threshold, which for us would be pure overhead —
     * hence the explicit check first. */
    if (n < (long long)NDARRAY_THREAD_THRESHOLD) return false;

    ParCtx ctx; ctx.pl = pl; ctx.parent = R;
    return nd_parallel_for((size_t)n, par_chunk, &ctx);
#else
    (void)pl; (void)R;
    return false;
#endif
}

static void vm_run(const Instr* code, size_t n, Slot* R, bool* failed) {
    *failed = false;
    if (n == 0) return;
    size_t pc = 0;
    const Instr* c = &code[pc];
    /* Live elements of the current tile, set by VSETLEN.  Only the load, store
     * and reduction read it — every arithmetic tile op covers the full VBLOCK —
     * and all three are pinned inside the strip loop because they read the loop
     * index, so this cannot be stranded by a hoist. */
    int vlen = VBLOCK;
    /* Operands are addressed lazily: an opcode pays only for the registers it
     * actually reads.  Computing all three up front in NEXT() cost three
     * shift-and-adds per instruction where the hottest ops (ADD_R, MUL_R) use
     * two and the unary ops use one. */
    #define RD (R[c->dst])
    #define RA (R[c->a])
    #define RB (R[c->b])
#if VM_THREADED
    #define X(name, kind) [OP_##name] = &&L_##name,
    static const void* const tbl[] = { OPLIST };
    #undef X
    #define OP(name) L_##name
    #define NEXT() do { c = &code[++pc]; goto *tbl[c->op]; } while (0)
    #define JUMP() do { c = &code[pc]; goto *tbl[c->op]; } while (0)
    goto *tbl[c->op];
#else
    #define OP(name) case OP_##name
    #define NEXT() break
    #define JUMP() continue
    while (pc < n) {
        c = &code[pc];
        switch (c->op) {
#endif
    #define ARROP() do { if (!vm_array_op(c, &RD, &RA, &RB)) goto vm_fail; } while (0); NEXT()
    /* An integer opcode that can overflow.  The helper returns true on overflow,
     * and the call is abandoned exactly as OP_FAIL abandons it, so the caller
     * falls back to the interpreter and gets the exact bigint answer.
     *
     * Two levels of opting out, and they measure different things:
     *
     *   IF_NOCHK in the instruction (COMPILE_WRAP_INT / Compile[]'s
     *   "CatchMachineIntegerOverflow" -> False) keeps the wrapped value instead
     *   of falling back.  Because `&&` short-circuits, the flag is read only
     *   once an overflow has actually been detected, so a program that does not
     *   overflow executes the identical instruction path either way — the option
     *   costs nothing.
     *
     *   `-DVM_NO_INT_CHECK` at BUILD time removes the detection itself.  That is
     *   the one that measures what the feature really costs, because it is the
     *   detection — not the never-taken branch — that is the price.  (Named for
     *   the VM, like VM_NO_THREADED, and deliberately NOT after the
     *   COMPILE_WRAP_INT flag: compile.h defines that as a bit value, so an
     *   `#ifdef` on it is true in every build and silently disables the
     *   checks — which is exactly the bug this comment now prevents.) */
#ifdef VM_NO_INT_CHECK
    #define IOP(chk) do { (void)(chk); } while (0); NEXT()
#else
    #define IOP(chk) do { if ((chk) && !(c->flags & IF_NOCHK)) goto vm_fail; } while (0); NEXT()
#endif
            OP(JMP): pc = c->b; JUMP();
            OP(JZ):  pc = RA.i ? pc + 1 : c->b; JUMP();   /* branch if false */
            /* Checked: this is Increment/Decrement's opcode as well as a loop
             * step (compile.c's Increment lowering), and `x = 2^63-1; x++` is a
             * perfectly ordinary thing to write. */
            OP(INC_I): IOP(ci_add(RD.i, c->imm.i, &RD.i));
            /* Increment, test and branch in one.  A counted loop otherwise spends
             * four instructions per iteration on control alone (INC, LT, JZ,
             * JMP), which is most of the body when the body is one element of a
             * fused array pass.
             *
             * Deliberately NOT overflow-checked, unlike INC_I.  The index only
             * overflows if the limit sits within one step of INT64_MAX, and
             * reaching there from the initial value takes on the order of 10^18
             * iterations — no terminating run gets close.  The check would
             * otherwise land in the innermost loop of every fused array pass. */
            OP(LOOP): { RD.i += c->imm.i; if (RD.i < RA.i) { pc = c->b; JUMP(); } } NEXT();
            /* Run the strip loop that follows across threads and skip past it;
             * on a decline, fall through and run it right here, serially. */
            OP(APAR): if (vm_par_run((const ParLoop*)c->imm.p, R)) { pc = c->b; JUMP(); }
                      NEXT();
            OP(CONST): RD = c->imm; NEXT();
            OP(MOVE):  RD = RA; NEXT();
            OP(I2R): RD.r = (double)RA.i; NEXT();
            OP(I2C): RD.z = (double)RA.i; NEXT();
            OP(R2C): RD.z = RA.r; NEXT();
            OP(ADD_I): IOP(ci_add(RA.i, RB.i, &RD.i));
            OP(ADD_R): RD.r = RA.r + RB.r; NEXT();
            OP(ADD_C): RD.z = RA.z + RB.z; NEXT();
            /* Immediate forms (K_BINK).  One register read instead of two, and
             * the CONST that would have materialised the operand is gone.  Each
             * is a SINGLE arithmetic operation, so no floating-point contraction
             * is possible and the result is bit-identical to the register form
             * the optimiser rewrote. */
            OP(ADD_RK): RD.r = RA.r + c->imm.r; NEXT();
            OP(SUB_RK): RD.r = RA.r - c->imm.r; NEXT();
            OP(SUB_KR): RD.r = c->imm.r - RA.r; NEXT();
            OP(MUL_RK): RD.r = RA.r * c->imm.r; NEXT();
            OP(DIV_RK): RD.r = RA.r / c->imm.r; NEXT();
            OP(DIV_KR): RD.r = c->imm.r / RA.r; NEXT();
            OP(ADD_IK): IOP(ci_add(RA.i, c->imm.i, &RD.i));
            OP(SUB_IK): IOP(ci_sub(RA.i, c->imm.i, &RD.i));
            OP(SUB_KI): IOP(ci_sub(c->imm.i, RA.i, &RD.i));
            OP(MUL_IK): IOP(ci_mul(RA.i, c->imm.i, &RD.i));
            OP(LT_RK): RD.i = RA.r <  c->imm.r; NEXT();
            OP(LE_RK): RD.i = RA.r <= c->imm.r; NEXT();
            OP(GT_RK): RD.i = RA.r >  c->imm.r; NEXT();
            OP(GE_RK): RD.i = RA.r >= c->imm.r; NEXT();
            OP(LT_IK): RD.i = RA.i <  c->imm.i; NEXT();
            OP(LE_IK): RD.i = RA.i <= c->imm.i; NEXT();
            OP(GT_IK): RD.i = RA.i >  c->imm.i; NEXT();
            OP(GE_IK): RD.i = RA.i >= c->imm.i; NEXT();
            OP(SUB_I): IOP(ci_sub(RA.i, RB.i, &RD.i));
            OP(SUB_R): RD.r = RA.r - RB.r; NEXT();
            OP(SUB_C): RD.z = RA.z - RB.z; NEXT();
            OP(MUL_I): IOP(ci_mul(RA.i, RB.i, &RD.i));
            OP(MUL_R): RD.r = RA.r * RB.r; NEXT();
            OP(MUL_C): RD.z = RA.z * RB.z; NEXT();
            OP(DIV_R): RD.r = RA.r / RB.r; NEXT();
            OP(DIV_C): RD.z = RA.z / RB.z; NEXT();
            /* Both integer divisions guard the two inputs the hardware traps on
             * rather than merely answering wrongly: a zero divisor, and
             * INT64_MIN / -1 whose quotient is not representable.  `Mod[5, 0]`
             * is left unevaluated by the interpreter, and the compiled path used
             * to take the whole process down with SIGFPE. */
            OP(MOD_I): { long long m = RB.i, x = RA.i;
                         if (m == 0 || (m == -1 && x == LLONG_MIN)) goto vm_fail;
                         long long q = x % m; if (q != 0 && ((q < 0) != (m < 0))) q += m;
                         RD.i = q; } NEXT();
            OP(QUOT_I): { long long m = RB.i, x = RA.i;
                          if (m == 0 || (m == -1 && x == LLONG_MIN)) goto vm_fail;
                          long long q = x / m; if ((x % m != 0) && ((x < 0) != (m < 0))) q -= 1;
                          RD.i = q; } NEXT();
            OP(NEG_I): IOP(ci_neg(RA.i, &RD.i));
            OP(NEG_R): RD.r = -RA.r; NEXT();
            OP(NEG_C): RD.z = -RA.z; NEXT();
            OP(INV_R): RD.r = 1.0 / RA.r; NEXT();
            OP(INV_C): RD.z = 1.0 / RA.z; NEXT();
            OP(POWI_I): IOP(ci_powi(RA.i, c->imm.i, &RD.i));
            OP(POWI_R): RD.r = ipow_r(RA.r, c->imm.i); NEXT();
            OP(POWI_C): RD.z = ipow_c(RA.z, c->imm.i); NEXT();
            OP(POW_R): RD.r = pow(RA.r, RB.r); NEXT();
            OP(POW_C): RD.z = cpow(RA.z, RB.z); NEXT();
            OP(SQRT_R): RD.r = sqrt(RA.r); NEXT();
            OP(SQRT_C): RD.z = csqrt(RA.z); NEXT();
            OP(EXP_R): RD.r = exp(RA.r); NEXT();
            OP(EXP_C): RD.z = cexp(RA.z); NEXT();
            OP(LOG_R): RD.r = log(RA.r); NEXT();
            OP(LOG_C): RD.z = clog(RA.z); NEXT();
            OP(SIN_R): RD.r = sin(RA.r); NEXT();   OP(SIN_C): RD.z = csin(RA.z); NEXT();
            OP(COS_R): RD.r = cos(RA.r); NEXT();   OP(COS_C): RD.z = ccos(RA.z); NEXT();
            OP(TAN_R): RD.r = tan(RA.r); NEXT();   OP(TAN_C): RD.z = ctan(RA.z); NEXT();
            OP(SINH_R): RD.r = sinh(RA.r); NEXT(); OP(SINH_C): RD.z = csinh(RA.z); NEXT();
            OP(COSH_R): RD.r = cosh(RA.r); NEXT(); OP(COSH_C): RD.z = ccosh(RA.z); NEXT();
            OP(TANH_R): RD.r = tanh(RA.r); NEXT(); OP(TANH_C): RD.z = ctanh(RA.z); NEXT();
            OP(ASIN_R): RD.r = asin(RA.r); NEXT(); OP(ASIN_C): RD.z = casin(RA.z); NEXT();
            OP(ACOS_R): RD.r = acos(RA.r); NEXT(); OP(ACOS_C): RD.z = cacos(RA.z); NEXT();
            OP(ATAN_R): RD.r = atan(RA.r); NEXT(); OP(ATAN_C): RD.z = catan(RA.z); NEXT();
            OP(ABS_I): IOP(ci_abs(RA.i, &RD.i));
            OP(ABS_R): RD.r = fabs(RA.r); NEXT();
            OP(ABS_C): RD.r = cabs(RA.z); NEXT();
            OP(SIGN_I): RD.i = (RA.i > 0) - (RA.i < 0); NEXT();
            /* Integer-closed heads.  Same IOP() shape as the arithmetic
             * opcodes: the helper reports overflow OR an out-of-domain
             * argument, and either way the interpreter takes the call. */
            OP(POW_II):  IOP(int_pow(RA.i, RB.i, &RD.i));
            OP(FACT_I):  IOP(int_factorial(RA.i, &RD.i));
            OP(GAMMA_I): IOP(int_gamma(RA.i, &RD.i));
            OP(BINOM_I): IOP(int_binomial(RA.i, RB.i, &RD.i));
            OP(POCH_I):  IOP(int_pochhammer(RA.i, RB.i, &RD.i));
            OP(FIB_I):   IOP(int_fib(RA.i, &RD.i));
            OP(LUCAS_I): IOP(int_lucas(RA.i, &RD.i));
            /* Arg[n] is 0 for n >= 0; below that it is the symbol Pi. */
            OP(ARG_I):   { if (RA.i < 0) goto vm_fail; RD.i = 0; } NEXT();
            OP(GCD_I):   IOP(int_gcd(RA.i, RB.i, &RD.i));
            OP(LCM_I):   IOP(int_lcm(RA.i, RB.i, &RD.i));
            OP(ILEN_I):  IOP(int_ilen(RA.i, RB.i, &RD.i));
            OP(IEXP_I):  IOP(int_iexp(RA.i, RB.i, &RD.i));
            /* Ternary, so the operands are a RUN of registers starting at `a`
             * (the K_NARY shape KERNN uses), not the two operand fields. */
            OP(POWMOD_I): IOP(int_powmod(R[c->a].i, R[c->a + 1].i, R[c->a + 2].i, &RD.i));
            /* Writes the INTEGER slot: Sign of a real is an Integer in the
             * interpreter.  (The tile form VSIGN_R stays real — it feeds packed
             * float arrays, where the ND kernel's result dtype is real.) */
            OP(SIGN_R): RD.i = (RA.r > 0) - (RA.r < 0); NEXT();
            OP(FLOOR_R): RD.i = (long long)floor(RA.r); NEXT();
            /* IntegerPart truncates TOWARD ZERO, which is not Floor for a
             * negative argument: IntegerPart[-1.5] is -1, Floor[-1.5] is -2. */
            OP(TRUNC_R): RD.i = (long long)trunc(RA.r); NEXT();
            OP(CEIL_R):  RD.i = (long long)ceil(RA.r); NEXT();
            OP(ROUND_R): RD.i = (long long)llround(RA.r); NEXT();
            OP(RE_C): RD.r = creal(RA.z); NEXT();
            OP(IM_C): RD.r = cimag(RA.z); NEXT();
            OP(ARG_C): RD.r = carg(RA.z); NEXT();
            OP(CONJ_C): RD.z = conj(RA.z); NEXT();
            OP(ATAN2_R): RD.r = atan2(RB.r, RA.r); NEXT();   /* ArcTan[x,y]=atan2(y,x); a=x,b=y */
            OP(MAX_I): RD.i = RA.i > RB.i ? RA.i : RB.i; NEXT();
            OP(MAX_R): RD.r = RA.r > RB.r ? RA.r : RB.r; NEXT();
            OP(MIN_I): RD.i = RA.i < RB.i ? RA.i : RB.i; NEXT();
            OP(MIN_R): RD.r = RA.r < RB.r ? RA.r : RB.r; NEXT();
            OP(ERF_R): RD.r = erf(RA.r); NEXT();
            OP(ERFC_R): RD.r = erfc(RA.r); NEXT();
            OP(KERN_RR): { double o; RD.r = ((kfn_r)c->imm.p)(RA.r, &o) ? o : NAN; } NEXT();
            OP(KERN_R2R): { double orr, oi; RD.r = ((kfn_c)c->imm.p)(RA.r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN_RC): { double orr, oi; if (((kfn_c)c->imm.p)(RA.r, 0.0, &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN_CC): { double orr, oi; if (((kfn_c)c->imm.p)(creal(RA.z), cimag(RA.z), &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN_CR): { double orr, oi; RD.r = ((kfn_c)c->imm.p)(creal(RA.z), cimag(RA.z), &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RR): { double orr, oi; RD.r = ((kfn_c2)c->imm.p)(RA.r, 0.0, RB.r, 0.0, &orr, &oi) ? orr : NAN; } NEXT();
            OP(KERN2_RC): { double orr, oi; if (((kfn_c2)c->imm.p)(RA.r, 0.0, RB.r, 0.0, &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(KERN2_CC): { double orr, oi; if (((kfn_c2)c->imm.p)(creal(RA.z), cimag(RA.z), creal(RB.z), cimag(RB.z), &orr, &oi)) RD.z = orr + oi * I; else RD.z = NAN + NAN * I; } NEXT();
            OP(LT_I): RD.i = RA.i < RB.i; NEXT();  OP(LT_R): RD.i = RA.r < RB.r; NEXT();
            OP(LE_I): RD.i = RA.i <= RB.i; NEXT(); OP(LE_R): RD.i = RA.r <= RB.r; NEXT();
            OP(GT_I): RD.i = RA.i > RB.i; NEXT();  OP(GT_R): RD.i = RA.r > RB.r; NEXT();
            OP(GE_I): RD.i = RA.i >= RB.i; NEXT(); OP(GE_R): RD.i = RA.r >= RB.r; NEXT();
            OP(EQ_I): RD.i = RA.i == RB.i; NEXT(); OP(EQ_R): RD.i = RA.r == RB.r; NEXT();
            OP(EQ_C): RD.i = RA.z == RB.z; NEXT();
            OP(NE_I): RD.i = RA.i != RB.i; NEXT(); OP(NE_R): RD.i = RA.r != RB.r; NEXT();
            OP(NE_C): RD.i = RA.z != RB.z; NEXT();
            /* SameQ on machine numbers, matching expr_eq (src/expr.c:622):
             * two NaNs ARE the same, which is what lets an iteration whose
             * orbit reaches NaN terminate instead of spinning.  Complex is
             * componentwise because expr_eq recurses into Complex[re, im]. */
            OP(SAMEQ_R): RD.i = (RA.r == RB.r) || (isnan(RA.r) && isnan(RB.r)); NEXT();
            OP(SAMEQ_C): {
                double ar = creal(RA.z), ai = cimag(RA.z);
                double br = creal(RB.z), bi = cimag(RB.z);
                RD.i = ((ar == br) || (isnan(ar) && isnan(br)))
                    && ((ai == bi) || (isnan(ai) && isnan(bi)));
            } NEXT();
            OP(FAIL): goto vm_fail;
            OP(AND): RD.i = RA.i && RB.i; NEXT();
            OP(OR):  RD.i = RA.i || RB.i; NEXT();
            OP(XOR): RD.i = (!!RA.i) ^ (!!RB.i); NEXT();
            OP(NOT): RD.i = !RA.i; NEXT();
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
            /* ---- fused elementwise loop (M3b) ----------------------------
             * A_LOAD/A_STORE are the whole point of fusion: an elementwise
             * chain becomes ONE pass over the buffers driven by ordinary scalar
             * bytecode, instead of one full-length ND pass and one temporary
             * buffer per operation.  They are inline (not routed through
             * vm_array_op) because they run once per element.
             *
             * The source dtype is checked per element rather than hoisted: it is
             * a perfectly predicted branch next to the arithmetic it feeds, and
             * hoisting it would mean either specialising the loop at runtime or
             * refusing float32 arrays that work today. */
            OP(A_LOAD_R): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t k_ = (size_t)RB.i;
                if (A_->dtype == NDT_FLOAT64)      RD.r = ((const double*)A_->data)[k_];
                else if (A_->dtype == NDT_FLOAT32) RD.r = (double)((const float*)A_->data)[k_];
                else goto vm_fail;   /* promised real, buffer is complex */
            } NEXT();
            OP(A_LOAD_C): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t k_ = (size_t)RB.i;
                switch (A_->dtype) {
                    case NDT_FLOAT64:   RD.z = ((const double*)A_->data)[k_]; break;
                    case NDT_FLOAT32:   RD.z = (double)((const float*)A_->data)[k_]; break;
                    case NDT_COMPLEX64: RD.z = ((const double*)A_->data)[2*k_]
                                             + ((const double*)A_->data)[2*k_+1] * I; break;
                    default:            RD.z = (double)((const float*)A_->data)[2*k_]
                                             + (double)((const float*)A_->data)[2*k_+1] * I; break;
                }
            } NEXT();
            /* Integer element access.  Unlike the real form there is no width
             * variation to absorb — a program's integer buffers are always
             * NDT_INT64 (nd_own_copy refuses to build one from anything else),
             * so a foreign dtype here means the promised element type was not
             * kept and the call goes back to the interpreter. */
            OP(A_LOAD_I): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                if (A_->dtype != NDT_INT64) goto vm_fail;
                RD.i = ((const int64_t*)A_->data)[(size_t)RB.i];
            } NEXT();
            OP(A_STORE_I): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                ((int64_t*)A_->data)[(size_t)RA.i] = RB.i;
            } NEXT();
            OP(A_STORE_R): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                ((double*)A_->data)[(size_t)RA.i] = RB.r;
            } NEXT();
            OP(A_STORE_C): {
                NDArrayData* A_ = &RD.arr->data.ndarray;
                size_t k_ = (size_t)RA.i;
                ((double*)A_->data)[2*k_]   = creal(RB.z);
                ((double*)A_->data)[2*k_+1] = cimag(RB.z);
            } NEXT();
            /* Append at R[a], growing the buffer when it runs out.  The array
             * was allocated by this call and nothing else holds it, so growing
             * it in place is safe; `dims[0]` is the CAPACITY until A_TRUNC sets
             * the real length. */
            OP(A_PUSH): {
                Expr* A_ = RD.arr;
                if (!A_ || A_->type != EXPR_NDARRAY || A_->data.ndarray.rank != 1) goto vm_fail;
                long long k_ = RA.i;
                if (k_ < 0) goto vm_fail;
                NDArrayData* nd_ = &A_->data.ndarray;
                if (k_ >= nd_->dims[0]) {
                    int64_t cap_ = nd_->dims[0] > 0 ? nd_->dims[0] * 2 : 16;
                    if (cap_ <= k_) cap_ = k_ + 1;
                    if (cap_ > (int64_t)VM_ITER_SAFETY_CAP + 1) goto vm_fail;
                    size_t es_ = ndt_elem_size(nd_->dtype);
                    void* nb_ = realloc(nd_->data, (size_t)cap_ * es_);
                    if (!nb_) goto vm_fail;
                    nd_->data = nb_;
                    nd_->dims[0] = cap_;
                }
                if (AF_R(c->flags) == (unsigned)CT_COMPLEX) {
                    ((double*)nd_->data)[2*(size_t)k_]     = creal(RB.z);
                    ((double*)nd_->data)[2*(size_t)k_ + 1] = cimag(RB.z);
                } else ((double*)nd_->data)[(size_t)k_] = RB.r;
            } NEXT();
            /* Final length, so a capacity the run happened to reach is never
             * visible.  Shrinking in place: a failed realloc keeps the buffer,
             * which is still large enough, so only dims[0] has to be right. */
            OP(A_TRUNC): {
                Expr* A_ = RD.arr;
                if (!A_ || A_->type != EXPR_NDARRAY || A_->data.ndarray.rank != 1) goto vm_fail;
                long long n_ = RA.i;
                NDArrayData* nd_ = &A_->data.ndarray;
                if (n_ < 0 || n_ > nd_->dims[0]) goto vm_fail;
                size_t es_ = ndt_elem_size(nd_->dtype);
                void* nb_ = realloc(nd_->data, (size_t)(n_ ? n_ : 1) * es_);
                if (nb_) nd_->data = nb_;
                nd_->dims[0] = n_;
            } NEXT();
            OP(A_SIZE):     ARROP();
            OP(A_NEWLIKE):  ARROP();
            OP(A_SHAPECHK): ARROP();
            OP(A_COPY):     ARROP();
            OP(A_XFER):     ARROP();
            /* ---- indexed Part (M3c) --------------------------------------
             * One instruction per axis, resolving the subscript against that
             * axis and folding it into the running flat index.  The range check
             * has to be per axis: u[[1, n + 5]] on an n x n array is inside the
             * buffer and reads the row below, which is the one indexing bug a
             * flat bounds check cannot catch. */
            OP(A_AXIS): {
                const Expr* A_ = RB.arr;
                long long ax = c->imm.i;
                if (!A_ || A_->type != EXPR_NDARRAY
                    || ax >= (long long)A_->data.ndarray.rank) goto vm_fail;
                long long len = (long long)A_->data.ndarray.dims[ax];
                long long k_ = RA.i;
                if (k_ < 0) k_ = len + k_ + 1;          /* Part counts from the end */
                if (k_ < 1 || k_ > len) goto vm_fail;   /* -> interpreter, which reports it */
                RD.i = RD.i * len + (k_ - 1);
            } NEXT();
            OP(A_NEW):     do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_PART):    do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_PARTSET): do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();
            OP(A_NDFN):    do { if (!vm_range_array_op(c, R)) goto vm_fail; } while (0); NEXT();

            /* ---- strip-mined tile ops (M5b) -------------------------------
             * One opcode, VBLOCK elements, in a loop shaped so the C compiler
             * can vectorise it.  Arithmetic always covers the FULL tile: the
             * load pads a short tail with 1.0, so no op reads uninitialised
             * memory and only load/store/reduce need the live length. */
            #define TD_R  ((double*)RD.p)
            #define TA_R  ((const double*)RA.p)
            #define TB_R  ((const double*)RB.p)
            #define TD_C  ((double _Complex*)RD.p)
            #define TA_C  ((const double _Complex*)RA.p)
            #define TB_C  ((const double _Complex*)RB.p)
            #define VBIN(NAME, TAG, CTY, EXPR) \
                OP(NAME): { CTY* d_ = TD_##TAG; const CTY* a_ = TA_##TAG; \
                            const CTY* b_ = TB_##TAG; (void)a_; (void)b_; \
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = (EXPR); } NEXT()
            #define VUN(NAME, TAG, CTY, EXPR) \
                OP(NAME): { CTY* d_ = TD_##TAG; const CTY* a_ = TA_##TAG; \
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = (EXPR); } NEXT()

            OP(VSETLEN): {                    /* live elements of this tile */
                long long rem = RB.i - RA.i;
                vlen = (int)(rem < VBLOCK ? rem : VBLOCK);
                if (vlen < 0) vlen = 0;
            } NEXT();

            OP(VLOAD_R): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t off = (size_t)RB.i;
                double* d_ = TD_R;
                if (A_->dtype == NDT_FLOAT64) {
                    const double* s_ = (const double*)A_->data + off;
                    for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[k_];
                } else if (A_->dtype == NDT_FLOAT32) {
                    const float* s_ = (const float*)A_->data + off;
                    for (int k_ = 0; k_ < vlen; k_++) d_[k_] = (double)s_[k_];
                } else goto vm_fail;          /* promised real, buffer is complex */
                for (int k_ = vlen; k_ < VBLOCK; k_++) d_[k_] = 1.0;   /* safe pad */
            } NEXT();
            OP(VLOAD_C): {
                const NDArrayData* A_ = &RA.arr->data.ndarray;
                size_t off = (size_t)RB.i;
                double _Complex* d_ = TD_C;
                switch (A_->dtype) {
                    case NDT_FLOAT64: { const double* s_ = (const double*)A_->data + off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[k_]; break; }
                    case NDT_FLOAT32: { const float* s_ = (const float*)A_->data + off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = (double)s_[k_]; break; }
                    case NDT_COMPLEX64: { const double* s_ = (const double*)A_->data + 2 * off;
                        for (int k_ = 0; k_ < vlen; k_++) d_[k_] = s_[2*k_] + s_[2*k_+1] * I; break; }
                    default: { const float* s_ = (const float*)A_->data + 2 * off;
                        for (int k_ = 0; k_ < vlen; k_++)
                            d_[k_] = (double)s_[2*k_] + (double)s_[2*k_+1] * I; break; }
                }
                for (int k_ = vlen; k_ < VBLOCK; k_++) d_[k_] = 1.0;
            } NEXT();
            /* The store is where a fused loop honours the element-type promise.
             * A real-typed program that produced a non-finite element is exactly
             * the case where the interpreter would have gone complex (Sqrt of a
             * negative, Log of a negative) or hit a pole, and the delegated path
             * fails there too because its kernels reject non-finite elements.
             * Returning a buffer full of NaN instead would silently hand back a
             * real array where the interpreter gives a complex one.
             *
             * The check is accumulated branchlessly and tested once per tile, so
             * it does not break the vectorisation of the copy. */
            OP(VSTORE_R): {
                double* d_ = (double*)RD.arr->data.ndarray.data + (size_t)RA.i;
                const double* s_ = TB_R;
                int bad_ = 0;
                for (int k_ = 0; k_ < vlen; k_++) {
                    double x_ = s_[k_];
                    d_[k_] = x_;
                    bad_ |= !(x_ - x_ == 0.0);          /* NaN or +-Inf */
                }
                if (bad_) goto vm_fail;
            } NEXT();
            OP(VSTORE_C): {
                double* d_ = (double*)RD.arr->data.ndarray.data + 2 * (size_t)RA.i;
                const double _Complex* s_ = TB_C;
                int bad_ = 0;
                for (int k_ = 0; k_ < vlen; k_++) {
                    double re_ = creal(s_[k_]), im_ = cimag(s_[k_]);
                    d_[2*k_] = re_; d_[2*k_+1] = im_;
                    bad_ |= !(re_ - re_ == 0.0) | !(im_ - im_ == 0.0);
                }
                if (bad_) goto vm_fail;
            } NEXT();
            /* Sequential accumulation, deliberately: it reproduces the
             * element-at-a-time fused loop's summation order exactly, so
             * strip-mining does not change any already-tested result. */
            OP(VSUM_R): { const double* s_ = TA_R; double acc_ = RD.r;
                          for (int k_ = 0; k_ < vlen; k_++) acc_ += s_[k_];
                          RD.r = acc_; } NEXT();
            OP(VSUM_C): { const double _Complex* s_ = TA_C; double _Complex acc_ = RD.z;
                          for (int k_ = 0; k_ < vlen; k_++) acc_ += s_[k_];
                          RD.z = acc_; } NEXT();

            OP(VSPLAT_R): { double* d_ = TD_R; double v_ = RA.r;
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = v_; } NEXT();
            OP(VSPLAT_C): { double _Complex* d_ = TD_C; double _Complex v_ = RA.z;
                            for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = v_; } NEXT();
            OP(VR2C): { double _Complex* d_ = TD_C; const double* a_ = TA_R;
                        for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = a_[k_]; } NEXT();

            VBIN(VADD_R, R, double, a_[k_] + b_[k_]);
            VBIN(VSUB_R, R, double, a_[k_] - b_[k_]);
            VBIN(VMUL_R, R, double, a_[k_] * b_[k_]);
            VBIN(VDIV_R, R, double, a_[k_] / b_[k_]);
            VUN(VNEG_R, R, double, -a_[k_]);
            VUN(VINV_R, R, double, 1.0 / a_[k_]);
            VBIN(VADD_C, C, double _Complex, a_[k_] + b_[k_]);
            VBIN(VSUB_C, C, double _Complex, a_[k_] - b_[k_]);
            VBIN(VMUL_C, C, double _Complex, a_[k_] * b_[k_]);
            VBIN(VDIV_C, C, double _Complex, a_[k_] / b_[k_]);
            VUN(VNEG_C, C, double _Complex, -a_[k_]);
            VUN(VINV_C, C, double _Complex, 1.0 / a_[k_]);
            /* Small exponents are unrolled with EXACTLY the association
             * ipow_* uses (x^2 = x*x, x^3 = x*(x*x), x^4 = (x*x)*(x*x)), so the
             * result is bitwise identical to the general path while the tile
             * loop becomes a vectorisable multiply instead of a per-element
             * function call with a loop inside it. */
            OP(VPOWI_R): { double* d_ = TD_R; const double* a_ = TA_R; long long e_ = c->imm.i;
                if (e_ == 2)      for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_]; d_[k_] = t_ * t_; }
                else if (e_ == 3) for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_]; d_[k_] = t_ * (t_ * t_); }
                else if (e_ == 4) for (int k_ = 0; k_ < VBLOCK; k_++) { double t_ = a_[k_], s_ = t_ * t_; d_[k_] = s_ * s_; }
                else              for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = ipow_r(a_[k_], e_);
            } NEXT();
            OP(VPOWI_C): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C; long long e_ = c->imm.i;
                if (e_ == 2)      for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_]; d_[k_] = t_ * t_; }
                else if (e_ == 3) for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_]; d_[k_] = t_ * (t_ * t_); }
                else if (e_ == 4) for (int k_ = 0; k_ < VBLOCK; k_++) { double _Complex t_ = a_[k_], s_ = t_ * t_; d_[k_] = s_ * s_; }
                else              for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = ipow_c(a_[k_], e_);
            } NEXT();
            VBIN(VPOW_R, R, double, pow(a_[k_], b_[k_]));
            VBIN(VPOW_C, C, double _Complex, cpow(a_[k_], b_[k_]));
            VBIN(VATAN2_R, R, double, atan2(b_[k_], a_[k_]));
            VUN(VSQRT_R, R, double, sqrt(a_[k_]));   VUN(VSQRT_C, C, double _Complex, csqrt(a_[k_]));
            VUN(VEXP_R, R, double, exp(a_[k_]));    VUN(VEXP_C, C, double _Complex, cexp(a_[k_]));
            VUN(VLOG_R, R, double, log(a_[k_]));    VUN(VLOG_C, C, double _Complex, clog(a_[k_]));
            VUN(VSIN_R, R, double, sin(a_[k_]));    VUN(VSIN_C, C, double _Complex, csin(a_[k_]));
            VUN(VCOS_R, R, double, cos(a_[k_]));    VUN(VCOS_C, C, double _Complex, ccos(a_[k_]));
            VUN(VTAN_R, R, double, tan(a_[k_]));    VUN(VTAN_C, C, double _Complex, ctan(a_[k_]));
            VUN(VSINH_R, R, double, sinh(a_[k_]));   VUN(VSINH_C, C, double _Complex, csinh(a_[k_]));
            VUN(VCOSH_R, R, double, cosh(a_[k_]));   VUN(VCOSH_C, C, double _Complex, ccosh(a_[k_]));
            VUN(VTANH_R, R, double, tanh(a_[k_]));   VUN(VTANH_C, C, double _Complex, ctanh(a_[k_]));
            VUN(VASIN_R, R, double, asin(a_[k_]));   VUN(VASIN_C, C, double _Complex, casin(a_[k_]));
            VUN(VACOS_R, R, double, acos(a_[k_]));   VUN(VACOS_C, C, double _Complex, cacos(a_[k_]));
            VUN(VATAN_R, R, double, atan(a_[k_]));   VUN(VATAN_C, C, double _Complex, catan(a_[k_]));
            VUN(VABS_R, R, double, fabs(a_[k_]));
            VUN(VSIGN_R, R, double, (double)((a_[k_] > 0) - (a_[k_] < 0)));
            VUN(VERF_R, R, double, erf(a_[k_]));
            VUN(VERFC_R, R, double, erfc(a_[k_]));
            /* complex tile in, real tile out */
            OP(VABS_C): { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = cabs(a_[k_]); } NEXT();
            OP(VRE_C):  { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = creal(a_[k_]); } NEXT();
            OP(VIM_C):  { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = cimag(a_[k_]); } NEXT();
            OP(VARG_C): { double* d_ = TD_R; const double _Complex* a_ = TA_C;
                          for (int k_ = 0; k_ < VBLOCK; k_++) d_[k_] = carg(a_[k_]); } NEXT();
            VUN(VCONJ_C, C, double _Complex, conj(a_[k_]));

            /* Kernel tiles: the indirect call per element blocks vectorisation,
             * but these are libm-class kernels where the call dominates anyway —
             * the win here is purely the amortised dispatch. */
            OP(VKERN_RR): { double* d_ = TD_R; const double* a_ = TA_R; double o_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_r)c->imm.p)(a_[k_], &o_) ? o_ : NAN; } NEXT();
            OP(VKERN_R2R): { double* d_ = TD_R; const double* a_ = TA_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(a_[k_], 0.0, &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN_RC): { double _Complex* d_ = TD_C; const double* a_ = TA_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(a_[k_], 0.0, &or_, &oi_) ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN_CC): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), &or_, &oi_)
                             ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN_CR): { double* d_ = TD_R; const double _Complex* a_ = TA_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN2_RR): { double* d_ = TD_R; const double* a_ = TA_R; const double* b_ = TB_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(a_[k_], 0.0, b_[k_], 0.0, &or_, &oi_) ? or_ : NAN; } NEXT();
            OP(VKERN2_RC): { double _Complex* d_ = TD_C; const double* a_ = TA_R; const double* b_ = TB_R; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(a_[k_], 0.0, b_[k_], 0.0, &or_, &oi_)
                             ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            OP(VKERN2_CC): { double _Complex* d_ = TD_C; const double _Complex* a_ = TA_C;
                             const double _Complex* b_ = TB_C; double or_, oi_;
                for (int k_ = 0; k_ < VBLOCK; k_++)
                    d_[k_] = ((kfn_c2)c->imm.p)(creal(a_[k_]), cimag(a_[k_]), creal(b_[k_]), cimag(b_[k_]),
                                                &or_, &oi_) ? or_ + oi_ * I : NAN + NAN * I; } NEXT();
            #undef VBIN
            #undef VUN
            /* A compiled callee, invoked directly: machine values in, machine
             * value out, no Expr and no evaluator round-trip.  This is what
             * lifts the inliner's depth cap — beyond it the callee is CALLED
             * rather than pasted in, so deep chains compile instead of bailing. */
            OP(CALL): { if (!vm_call((const CompiledProgram*)c->imm.p,
                                     &RA, c->flags, &RD)) goto vm_fail; } NEXT();
            /* An n-ary machine kernel: arguments in `flags` consecutive
             * registers starting at `a`, the kernel pointer in the immediate. */
            OP(KERNN): {
                double ar[8], ai[8], orr, oi;
                unsigned na_ = c->flags;
                if (na_ > 8) goto vm_fail;
                for (unsigned k_ = 0; k_ < na_; k_++) { ar[k_] = R[c->a + k_].r; ai[k_] = 0.0; }
                if (!((kfn_cn)c->imm.p)(ar, ai, (size_t)na_, &orr, &oi)) RD.r = NAN;
                else RD.r = orr;
            } NEXT();
            OP(NOP): NEXT();
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

/* ------------------------------------------------------------------ *
 *  Common-subexpression elimination, at the Expr level                *
 * ------------------------------------------------------------------ *
 * The bytecode optimiser has a value-numbering CSE, and it almost never fires.
 * The reason is structural rather than a bug in the pass: binop/unop pop their
 * operands BEFORE allocating the destination, so a computation normally writes
 * into one of its own operand registers (`SIN_R t1, t1`).  The value-number
 * entry keyed on that register is therefore invalidated by the very instruction
 * that created it, and the value is genuinely gone — no bytecode-level pass can
 * recover it without first changing how registers are allocated.
 *
 * Doing it on the Expr instead sidesteps all of that, and is where the
 * information is richest anyway: repeated subtrees are visible directly.  A
 * chosen subtree is computed ONCE into a register reserved BELOW the temp stack
 * (by raising `nlocals`), which is exactly the discipline that already protects
 * arguments — no temp allocation, and none of the explicit `temp_top` resets in
 * the loop lowerings, can reach it.
 *
 * Eligibility is deliberately narrow, and every clause earns its place:
 *   - every free symbol must be an ARGUMENT.  A subtree mentioning a loop
 *     variable or a With/Module local cannot be hoisted to program entry, where
 *     that name is not bound yet.  (Loop-INVARIANT hoisting is the optimiser's
 *     LICM pass; this is a different job.)
 *   - no control flow and no assignment anywhere inside, so evaluating it early
 *     and unconditionally is observationally identical.
 *   - scalar-typed: an array value carries ownership, and a tile belongs to the
 *     fused loop that created it.
 */
#define CSE_MAX_NODES 4000   /* skip planning on very large trees */

static bool cse_head_is_pure(const char* h) {
    static const char* impure[] = {
        "Set", "SetDelayed", "AddTo", "SubtractFrom", "TimesBy", "DivideBy",
        "Increment", "Decrement", "CompoundExpression", "Do", "While", "For",
        "Sum", "Product", "Nest", "If", "With", "Module", "Block", "Function",
        "Which", "Piecewise", "Total", "Length"
    };
    for (size_t i = 0; i < sizeof impure / sizeof impure[0]; i++)
        if (strcmp(h, impure[i]) == 0) return false;
    return true;
}

/* Every free symbol an argument, every head pure.  Also counts nodes. */
static bool cse_eligible(Ctx* c, const Expr* e, int* nodes) {
    if (!e) return false;
    (*nodes)++;
    if (*nodes > CSE_MAX_NODES) return false;
    Slot imm; CompileType lt;
    if (literal(e, &imm, &lt)) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm_get(&c->map, nm) >= 0) return true;          /* an argument */
        double cv;
        if (strcmp(nm, "I") == 0 || named_const(nm, &cv)) return true;
        return false;                    /* loop var, local, or free symbol */
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    if (!cse_head_is_pure(e->data.function.head->data.symbol.name)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!cse_eligible(c, e->data.function.args[i], nodes)) return false;
    return true;
}

typedef struct { const Expr* rep; uint64_t hash; int nodes; int count; int slot; } CseCand;

/* Collect eligible subtrees, count structural duplicates, and record where each
 * one occurs.  Occurrences are keyed on the Expr POINTER — the planner walks the
 * very tree the emitter will walk, so pointer identity is exact and the
 * emit-time lookup needs no structural comparison. */
static void cse_scan(Ctx* c, const Expr* e, CseCand* cand, int* ncand) {
    if (!e || e->type != EXPR_FUNCTION) return;
    int nodes = 0;
    if (e->data.function.head->type == EXPR_SYMBOL && cse_eligible(c, e, &nodes) && nodes >= 3) {
        CompileType t;
        if (infer_type(c, e, &t) && (int)t >= 0 && !CT_IS_ARRAY(t)) {
            uint64_t h = expr_hash(e);
            int found = -1;
            for (int i = 0; i < *ncand; i++)
                if (cand[i].hash == h && expr_eq(cand[i].rep, e)) { found = i; break; }
            if (found < 0 && *ncand < CSE_MAX * 4) {
                found = (*ncand)++;
                cand[found].rep = e; cand[found].hash = h;
                cand[found].nodes = nodes; cand[found].count = 0; cand[found].slot = -1;
            }
            if (found >= 0) {
                cand[found].count++;
                if (c->cse_nocc < CSE_OCC_MAX) {
                    c->cse_occ[c->cse_nocc].node = e;
                    c->cse_occ[c->cse_nocc].idx  = found;
                    c->cse_nocc++;
                }
            }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        cse_scan(c, e->data.function.args[i], cand, ncand);
}

/* Slot of the CSE register holding `e`, or -1. */
static int cse_lookup(const Ctx* c, const Expr* e) {
    for (size_t i = 0; i < c->cse_nocc; i++)
        if (c->cse_occ[i].node == e) return c->cse_occ[i].idx;
    return -1;
}

/* Choose repeated subtrees, reserve a register for each, and emit each once.
 * Smallest first, so a larger candidate containing a smaller one picks up the
 * smaller one's register rather than recomputing it. */
static void cse_plan(Ctx* c, const Expr* body) {
    CseCand cand[CSE_MAX * 4];
    int ncand = 0;
    cse_scan(c, body, cand, &ncand);
    if (!ncand) { c->cse_nocc = 0; return; }

    int chosen[CSE_MAX], nchosen = 0;
    while (nchosen < CSE_MAX) {
        int best = -1;
        for (int i = 0; i < ncand; i++) {
            if (cand[i].count < 2 || cand[i].slot >= 0) continue;
            if (best < 0 || cand[i].nodes < cand[best].nodes) best = i;
        }
        if (best < 0) break;
        cand[best].slot = nchosen;
        c->cse_reg[nchosen] = c->nlocals + nchosen;   /* reserved BELOW the temps */
        c->cse_ready[nchosen] = false;
        chosen[nchosen++] = best;
    }
    if (!nchosen) { c->cse_nocc = 0; return; }

    /* Keep only occurrences of chosen candidates, remapped to their slot. */
    size_t w = 0;
    for (size_t i = 0; i < c->cse_nocc; i++) {
        int slot = cand[c->cse_occ[i].idx].slot;
        if (slot < 0) continue;
        c->cse_occ[w].node = c->cse_occ[i].node;
        c->cse_occ[w].idx  = slot;
        w++;
    }
    c->cse_nocc = w;

    c->nlocals += nchosen;
    if (c->nlocals > c->maxreg) c->maxreg = c->nlocals;
    c->ncse = nchosen;

    /* Emit each chosen subtree once.  `cse_ready` is switched on as we go, so
     * candidate k reuses candidates 0..k-1 and can never reference itself. */
    for (int k = 0; k < nchosen; k++) {
        Val v;
        if (!emit(c, cand[chosen[k]].rep, &v) || !c->ok) {
            /* Not compilable in isolation after all — drop this one and every
             * later one, and carry on without them rather than failing. */
            c->ok = true;
            c->ncse = k;
            for (size_t i = 0; i < c->cse_nocc; ) {
                if (c->cse_occ[i].idx >= k) c->cse_occ[i] = c->cse_occ[--c->cse_nocc];
                else i++;
            }
            return;
        }
        Slot z; memset(&z, 0, sizeof z);
        ins(c, OP_MOVE, (uint32_t)c->cse_reg[k], (uint32_t)v.reg, 0, z);
        free_if_tmp(c, v);
        c->cse_type[k]  = v.type;
        c->cse_ready[k] = true;
    }
}

/* ------------------------------------------------------------------ *
 *  Public API                                                         *
 * ------------------------------------------------------------------ */
/* Rewrite a virtual array register to its final slot in the array bank above
 * the scalar registers.  Ordinary register numbers and jump targets (also
 * carried in the `b` field) are far below ARR_VREG and pass through. */
/* Lift every OP_APAR's strip loop into a standalone ParLoop.
 *
 * Runs AFTER the optimiser, which is the only correct place: LICM and compaction
 * both move instruction indices, so a range recorded at emit time would name the
 * wrong instructions by the time it mattered.
 *
 * An APAR whose loop cannot be lifted (a branch leaving the range — nothing
 * emits one today, but an `If` lowering that grew a jump past the loop would)
 * is turned into a NOP, and execution simply falls into the serial loop that is
 * still sitting right there. Declining to parallelise is always safe; the serial
 * path is the same code. */
static void extract_par_loops(CompiledProgram* p) {
    int cap = 0;
    for (size_t i = 0; i < p->n; i++) if (p->code[i].op == OP_APAR) cap++;
    if (!cap) return;
    p->ploops = calloc((size_t)cap, sizeof(ParLoop));
    if (!p->ploops) return;

    for (size_t i = 0; i < p->n; i++) {
        if (p->code[i].op != OP_APAR) continue;
        size_t lo = i + 1, hi = p->code[i].b;      /* body is [lo, hi) */
        if (hi <= lo || hi > p->n) { p->code[i].op = OP_NOP; continue; }

        bool liftable = true;
        for (size_t j = lo; j < hi && liftable; j++) {
            int k = compile_op_kind[p->code[j].op];
            if (k != K_JMP && k != K_JZ && k != K_LOOP && k != K_APAR) continue;
            size_t t = p->code[j].b;
            if (t < lo || t > hi) liftable = false;   /* escapes the range */
        }
        size_t len = hi - lo;
        Instr* sub = liftable ? malloc((len + 1) * sizeof(Instr)) : NULL;
        if (!sub) { p->code[i].op = OP_NOP; continue; }

        memcpy(sub, &p->code[lo], len * sizeof(Instr));
        for (size_t j = 0; j < len; j++) {
            int k = compile_op_kind[sub[j].op];
            if (k == K_JMP || k == K_JZ || k == K_LOOP || k == K_APAR)
                sub[j].b = (uint32_t)(sub[j].b - lo);
        }
        /* Falling off the end of the lifted loop must RETURN, not run into
         * whatever followed it in the parent program. */
        memset(&sub[len], 0, sizeof sub[len]);
        sub[len].op = OP_RET;

        ParLoop* pl = &p->ploops[p->nploops];
        pl->code = sub; pl->n = len + 1;
        pl->ri = p->code[i].dst; pl->rn = p->code[i].a;
        pl->frame_slots = p->frame_slots; pl->nreg = p->nreg;
        pl->tile_base = p->tile_base;     pl->ntiles = p->ntiles;
        /* The array is fully sized up front, so this pointer stays valid. */
        p->code[i].imm.p = pl;
        p->nploops++;
    }
}

static uint32_t patch_reg(uint32_t r, int arr_base, int tile_base) {
    if (r >= (uint32_t)ARR_VREG)  return (uint32_t)arr_base  + (r - (uint32_t)ARR_VREG);
    if (r >= (uint32_t)TILE_VREG) return (uint32_t)tile_base + (r - (uint32_t)TILE_VREG);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Bail diagnostics                                                    *
 *                                                                      *
 *  A bail is silent by construction: compile_expr returns NULL, the     *
 *  caller interprets, the answer is still right, and the only symptom   *
 *  is being 10-40x slower.  Worse, the compilable subset is a cliff —   *
 *  one unsupported head costs the WHOLE body.  So the emitter records   *
 *  the innermost subexpression it could not lower, and this is where    *
 *  that becomes readable.                                               *
 *                                                                      *
 *  Single-writer state: compilation happens on the calling thread only  *
 *  (workers execute finished programs, they never compile), and the     *
 *  values are valid until the next compile call.                        *
 * ------------------------------------------------------------------ */

static char*       g_bail_expr   = NULL;   /* printed innermost failing node */
static const char* g_bail_reason = NULL;   /* static classification, or NULL */

static void bail_clear(void) {
    free(g_bail_expr);
    g_bail_expr   = NULL;
    g_bail_reason = NULL;
}

/* What kind of thing the emitter choked on.  Derived from the node rather than
 * threaded out of each bail site, which keeps the ~40 sites free of diagnostic
 * bookkeeping at the cost of a coarser message. */
static const char* bail_classify(const Expr* e) {
    if (!e) return "empty expression";
    switch (e->type) {
        case EXPR_SYMBOL:
            return "symbol is not a declared argument and holds no machine value";
        case EXPR_STRING:
            return "a string is not a machine number";
        case EXPR_FUNCTION:
            if (e->data.function.head->type != EXPR_SYMBOL)
                return "the head is not a symbol";
            return "no machine lowering for this head at these argument types";
        default:
            return "not a machine number";
    }
}

static void bail_record(const Ctx* c) {
    bail_clear();
    if (!c->bail_node) { g_bail_reason = "the body has no machine result type"; return; }
    g_bail_reason = bail_classify(c->bail_node);
    g_bail_expr   = expr_to_string((Expr*)c->bail_node);
}

const char* compiled_bail_reason(void) { return g_bail_reason; }
const char* compiled_bail_expr(void)   { return g_bail_expr; }

CompiledProgram* compile_expr(const Expr* body, const char* const* arg_names,
                              const CompileType* arg_types, size_t nargs) {
    return compile_expr_ex(body, arg_names, arg_types, nargs, 0u);
}

CompiledProgram* compile_expr_ex(const Expr* body, const char* const* arg_names,
                                 const CompileType* arg_types, size_t nargs,
                                 unsigned flags) {
    bail_clear();
    if (!body) { g_bail_reason = "empty body"; return NULL; }
    Ctx c; memset(&c, 0, sizeof(c));
    for (size_t k = 0; k < nargs; k++)
        if (CT_IS_ARRAY(arg_types[k])) {
            /* Rank is bounded only by the packed type encoding: the fused
             * elementwise loop walks a flat index and never looks at the shape,
             * so rank > 1 needs no separate machinery.  Constructs that DO care
             * about the shape (Total over the leading axis, Length) still gate
             * on rank themselves. */
            if (CT_RANK(arg_types[k]) < 1 || CT_RANK(arg_types[k]) > CT_MAX_RANK) {
                g_bail_reason = "declared argument array rank is out of range";
                return NULL;
            }
            c.array_args = true;
        }
    c.flags = flags;
    c.ok = true; c.nlocals = (int)nargs; c.arg_types = arg_types;
    c.argdep = calloc(nargs ? nargs : 1, 1);
    if (!c.argdep) return NULL;
    if (!nm_init(&c.map, arg_names, nargs)) { free(c.argdep); return NULL; }
    c.maxreg = (int)nargs;

    /* Hoist repeated subtrees into reserved registers before lowering the body.
     * Must run before ANY temp is allocated: it works by raising `nlocals`, the
     * floor the temp allocator builds on. */
    if (!(flags & COMPILE_NO_CSE)) cse_plan(&c, body);

    Val res;
    bool ok = emit(&c, body, &res) && c.ok;
    /* A borrowed argument array cannot be the result: the caller owns what it
     * gets back, and freeing an argument would corrupt the caller's value. */
    if (ok && CT_IS_ARRAY(res.type) && !res.tmp) {
        ok = false;
        c.bail_node = body;   /* the whole body is a borrowed argument array */
    }
    /* A statement-shaped head lowers to the integer 0 where the interpreter
     * answers Null, so it may only appear where its VALUE is discarded.  Inside
     * a CompoundExpression that is exactly what happens (free_if_tmp drops it
     * and nothing observes it); in RESULT position the two disagree, which was
     * reachable as `Compile[{n}, Do[..., {n}]][3]` giving 0 against the
     * interpreter's Null.  Null is not worth a fifth type in the lattice — it
     * would ripple through num_common, coerce, finite_result and cf_unbox — so
     * the honest fix is to decline the one position where it shows. */
    if (ok && stmt_valued_head(body)) {
        ok = false;
        c.bail_node = body;
    }
    if (ok) { Slot z = { 0 }; ins(&c, OP_RET, (uint32_t)res.reg, 0, 0, z); ok = c.ok; }
    nm_free(&c.map);
    if (!ok) {
        bail_record(&c);
        for (int i = 0; i < c.nparts; i++) compile_partspec_free(c.parts[i]);
        free(c.parts);
        free(c.code); free(c.argdep); return NULL;
    }

    /* Three contiguous banks: scalars, then array handles, then strip-mining
     * tiles.  A slot therefore has one kind for the whole life of the program,
     * so teardown can never mistake a double for a pointer, and each bank is a
     * range sweep. */
    int arr_base  = c.maxreg;
    int tile_base = arr_base + c.arr_max;
    int nreg      = tile_base + c.tile_max;
    for (size_t i = 0; i < c.n; i++) {
        c.code[i].dst = patch_reg(c.code[i].dst, arr_base, tile_base);
        c.code[i].a   = patch_reg(c.code[i].a, arr_base, tile_base);
        /* `b` is a branch TARGET on the jumping kinds and a register everywhere
         * else.  Asking the kind table rather than listing the opcodes means a
         * new branch opcode cannot have its target silently relocated into the
         * array bank — which is a corruption with no symptom until it jumps. */
        int bk = compile_op_kind[c.code[i].op];
        if (bk != K_JMP && bk != K_JZ && bk != K_LOOP && bk != K_APAR)
            c.code[i].b = patch_reg(c.code[i].b, arr_base, tile_base);
    }
    int result_reg = (int)patch_reg((uint32_t)res.reg, arr_base, tile_base);

    /* Optimise the emitted bytecode.  Runs after patch_reg so the array bank is
     * already at its final place and `arr_base` means what the optimiser expects.
     * Register numbers are preserved, so `result_reg` stays valid.  A failure
     * here is non-fatal: the unoptimised code is still correct. */
    if (!(flags & COMPILE_NO_OPT)) compile_optimize(c.code, &c.n, nreg, arr_base, tile_base);

    CompiledProgram* p = calloc(1, sizeof(*p));
    if (!p) {
        for (int i = 0; i < c.nparts; i++) compile_partspec_free(c.parts[i]);
        free(c.parts);
        free(c.code); free(c.argdep); return NULL;
    }
    /* The general-Part subscript lists become the program's: their literal specs
     * are pointed at from instruction immediates and must outlive the body. */
    p->parts = c.parts; p->nparts = c.nparts;
    p->code = c.code; p->n = c.n; p->nreg = nreg; p->arr_base = arr_base;
    p->tile_base = tile_base;
    p->result_reg = result_reg; p->result_type = res.type; p->ncse = c.ncse;
    p->result_built = CT_IS_ARRAY(res.type) && res.built;
    p->nargs = nargs; p->argdep = c.argdep;
    p->arg_types = malloc((nargs ? nargs : 1) * sizeof(CompileType));
    p->ntiles = c.tile_max;
    p->frame_slots = (size_t)nreg + (size_t)c.tile_max * VBLOCK;
    if (p->frame_slots == 0) p->frame_slots = 1;
    if (!p->arg_types) { compiled_free(p); return NULL; }
    memcpy(p->arg_types, arg_types, nargs * sizeof(CompileType));
    p->all_real = (res.type == CT_REAL) && c.arr_max == 0;
    for (size_t k = 0; k < nargs; k++) if (arg_types[k] != CT_REAL) p->all_real = false;
    extract_par_loops(p);
    return p;
}

CompileType compiled_result_type(const CompiledProgram* p) { return p->result_type; }
bool        compiled_result_built(const CompiledProgram* p) { return p && p->result_built; }
size_t compiled_num_args(const CompiledProgram* p) { return p->nargs; }
size_t compiled_num_instructions(const CompiledProgram* p) { return p->n; }
bool compiled_program_all_real(const CompiledProgram* p) { return p && p->all_real; }
size_t compiled_num_cse(const CompiledProgram* p) { return p ? (size_t)p->ncse : 0; }

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
/* ------------------------------------------------------------------ *
 *  Call frames                                                        *
 * ------------------------------------------------------------------ *
 * The frame used to be a single buffer owned by the CompiledProgram and mutated
 * through a `const CompiledProgram*`, which meant a program was neither
 * reentrant nor thread-safe: two live calls of the same program shared one
 * register file, and (once fusion landed) one set of tile buffers.
 *
 * Frames now come from the C stack, which is per-thread and naturally nested, so
 * reentrancy and thread-safety both fall out with no arena to own, grow or free
 * — and no allocation at all for any program that fits the fixed buffer, which
 * is every realistic one.  Larger programs fall back to a single malloc.
 *
 * Layout: `nreg` register Slots, then `ntiles * VBLOCK` Slots of tile storage.
 * A Slot is exactly `sizeof(double _Complex)`, so tile storage is correctly
 * aligned for both real and complex tiles.  VM_STACK_SLOTS is defined up with
 * the parallel strip loop, which needs the same frame shape for its workers. */

/* Point each tile register at its slice of the frame's tile storage. */
static void frame_bind_tiles(const CompiledProgram* p, Slot* R) {
    Slot* tiles = R + p->nreg;
    for (int k = 0; k < p->ntiles; k++)
        R[p->tile_base + k].p = tiles + (size_t)k * VBLOCK;
}

/* Only the ARRAY bank: tile slots hold pointers into the frame's own storage,
 * so clearing them would strand the tiles for the rest of the call. */
static void arr_reset(const CompiledProgram* p, Slot* R) {
    for (int r = p->arr_base; r < p->tile_base; r++) R[r].arr = NULL;
}
static void arr_sweep(const CompiledProgram* p, Slot* R) {
    for (int r = p->arr_base; r < p->tile_base; r++)
        if (R[r].arr) { expr_free(R[r].arr); R[r].arr = NULL; }
}

/* Depth guard for OP_CALL.  Frames live on the C stack, so unbounded nesting
 * would overflow it rather than fail cleanly; a compiled program that recurses
 * past this simply fails and the caller falls back to the interpreter. */
#define VM_MAX_CALL_DEPTH 200
static VM_TLS int vm_call_depth = 0;

/* Run `cp` on `nargs` machine values taken straight from the caller's registers.
 * The caller coerced them to the callee's declared types at emit time, so the
 * copy is a raw Slot move — no boxing, no Expr, no evaluator.  The callee gets
 * its OWN frame, which is what makes a compiled program reentrant. */
static bool vm_call(const CompiledProgram* cp, const Slot* argv, unsigned nargs, Slot* dst) {
    if (!cp || cp->nargs != (size_t)nargs) return false;
    if (vm_call_depth >= VM_MAX_CALL_DEPTH) return false;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (cp->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(cp->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    if (cp->ntiles) frame_bind_tiles(cp, R);
    for (unsigned k = 0; k < nargs; k++) R[k] = argv[k];
    arr_reset(cp, R);

    vm_call_depth++;
    bool failed = false;
    vm_run(cp->code, cp->n, R, &failed);
    vm_call_depth--;

    Slot* r = &R[cp->result_reg];
    bool good = !failed && finite_result(r, cp->result_type)
                && !CT_IS_ARRAY(cp->result_type);
    if (good) *dst = *r;
    arr_sweep(cp, R);
    free(heap);
    return good;
}

bool compiled_eval(const CompiledProgram* p, const CompileValue* args, CompileValue* out) {
    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (p->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(p->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    if (p->ntiles) frame_bind_tiles(p, R);

    for (size_t k = 0; k < p->nargs; k++)
        if (!load_arg(&R[k], &args[k], p->arg_types[k])) { free(heap); return false; }
    arr_reset(p, R);
    bool failed = false;
    vm_run(p->code, p->n, R, &failed);
    Slot* r = &R[p->result_reg];
    out->type = p->result_type;
    if (CT_IS_ARRAY(p->result_type)) {
        out->v.a = failed ? NULL : r->arr;
        if (!failed) r->arr = NULL;    /* ownership transfers to the caller */
        arr_sweep(p, R);
        free(heap);
        return !failed && out->v.a != NULL;
    }
    switch (p->result_type) {
        case CT_BOOL: out->v.b = (unsigned char)(r->i != 0); break;
        case CT_INT:  out->v.i = r->i; break;
        case CT_REAL: out->v.r = r->r; break;
        case CT_COMPLEX: out->v.z = r->z; break;
        default: break;
    }
    arr_sweep(p, R);
    bool good = !failed && finite_result(r, p->result_type);
    free(heap);
    return good;
}

bool compiled_eval_real(const CompiledProgram* p, const double* args, double* out) {
    if (!p->all_real) return false;   /* implies no array registers and no tiles */
    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* R = stackframe;
    if (p->frame_slots > VM_STACK_SLOTS) {
        heap = malloc(p->frame_slots * sizeof(Slot));
        if (!heap) return false;
        R = heap;
    }
    for (size_t k = 0; k < p->nargs; k++) R[k].r = args[k];
    bool failed = false;
    vm_run(p->code, p->n, R, &failed);
    *out = R[p->result_reg].r;
    free(heap);
    /* `failed` was computed and dropped here for as long as no opcode an
     * all-Real program could contain was able to abort — array ops can, and an
     * all-Real program has no array registers.  OP_FAIL changed that: an
     * unbounded FixedPoint hitting the safety cap must drop the call to the
     * interpreter, not hand back whatever the accumulator happened to hold. */
    return !failed && isfinite(*out);
}

bool compiled_eval_real_batch(const CompiledProgram* const* progs, size_t nprogs,
                              const double* args, size_t nargs, double* out) {
    if (nprogs == 0) return true;
    /* one shared frame = the widest program's (big enough for every program's
     * registers); the argument region [0,nargs) is loaded once and never written
     * by any program (dst registers are always temporaries >= nargs). */
    size_t widest = 0;
    for (size_t i = 0; i < nprogs; i++)
        if (progs[i]->frame_slots > widest) widest = progs[i]->frame_slots;

    Slot stackframe[VM_STACK_SLOTS];
    Slot* heap = NULL;
    Slot* F = stackframe;
    if (widest > VM_STACK_SLOTS) {
        heap = malloc(widest * sizeof(Slot));
        if (!heap) return false;
        F = heap;
    }
    for (size_t k = 0; k < nargs; k++) F[k].r = args[k];
    for (size_t i = 0; i < nprogs; i++) {
        /* all_real implies no array registers and no tiles, so one shared frame
         * is enough and needs no per-program tile binding. */
        if (!progs[i]->all_real) { free(heap); return false; }
        bool failed = false;
        vm_run(progs[i]->code, progs[i]->n, F, &failed);
        out[i] = F[(size_t)progs[i]->result_reg].r;
        if (failed || !isfinite(out[i])) { free(heap); return false; }
    }
    free(heap);
    return true;
}

void compiled_free(CompiledProgram* p) {
    if (!p) return;
    for (int i = 0; i < p->nploops; i++) free(p->ploops[i].code);
    free(p->ploops);
    for (int i = 0; i < p->nparts; i++) compile_partspec_free(p->parts[i]);
    free(p->parts);
    free(p->code); free(p->arg_types); free(p->argdep);
    free(p);
}
