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
#include "compile_emit.h"        /* NameMap / Val / Ctx builder + emitter constants */
#include "compiled_function.h"   /* inlining a CompiledFunction callee */
#include "../arithmetic.h"
#include "../symtab.h"
#include "../attr.h"      /* ATTR_LISTABLE — the gate on elementwise fusion */
#include "../ndarray.h"    /* NDUnaryKernel / NDBinaryKernel — shared kernel layer */
#include "../assoc.h"      /* is_association / assoc_lookup_value / assoc_values_list — B1 */
#include "../print.h"      /* expr_to_string — printing the node a bail choked on */
#include "../sym_names.h" /* SYM_All / SYM_Span / SYM_List — Part subscript specs */
#include "../sym_intern.h" /* intern_symbol — the FN_HEAD placeholder parameters */
#ifdef USE_MPFR
#include "../numeric_complex.h"  /* ncpx — arbitrary-precision complex containers */
#endif
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
/* Slot, the kfn_* kernel signatures, Instr, the opcode enum, the AF_* array
 * flags, ParLoop and `struct CompiledProgram` itself all live in
 * compile_internal.h — shared with the optimiser and the disassembler. */

/* ------------------------------------------------------------------ *
 *  Argument name -> index map (interned-pointer hash)                 *
 * ------------------------------------------------------------------ */
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

/* The caller's argument symbols are NOT in scope inside an inlined
 * CompiledFunction body: the callee was compiled against its own parameters and
 * globals only, so a caller argument sharing a name with a callee global would
 * silently capture it. */
int arg_find(const Ctx* c, const char* nm) {
    return c->inlining ? -1 : nm_get(&c->map, nm);
}

/* Under COMPILE_FOLD_GLOBALS, resolve a non-argument symbol that currently
 * holds a machine-number OwnValue (`y = 0.37`) to that constant — the same value
 * the interpreter would substitute for it right now.  Only an unconditional
 * assignment counts: the rule's pattern must be the bare symbol, so patterned or
 * conditional OwnValues are left alone and the symbol stays uncompilable.
 * See COMPILE_FOLD_GLOBALS in compile.h for why this is opt-in. */
bool literal(const Expr* e, Slot* imm, CompileType* type);

bool global_const(const Ctx* c, const char* nm, Slot* imm, CompileType* type) {
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
const CompiledFunction* compiled_callee(const Ctx* c, const char* nm) {
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
int scope_find(const Ctx* c, const char* nm, CompileType* type, bool* built) {
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

void ins_f(Ctx* c, uint16_t op, uint16_t flags,
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
    if ((c->flags & COMPILE_WRAP_INT) && op_is_checked_int(op)
        && !(flags & IF_FORCECHK)) flags |= IF_NOCHK;
    c->code[c->n].op = op; c->code[c->n].flags = flags; c->code[c->n].dst = dst;
    c->code[c->n].a = a; c->code[c->n].b = b; c->code[c->n].imm = imm;
    c->n++;
}
void ins(Ctx* c, uint16_t op, uint32_t dst, uint32_t a, uint32_t b, Slot imm) {
    ins_f(c, op, 0, dst, a, b, imm);
}
int alloc_temp(Ctx* c) {
    int r = c->nlocals + c->temp_top;
    c->temp_top++;
    if (c->nlocals + c->temp_top > c->maxreg) c->maxreg = c->nlocals + c->temp_top;
    return r;
}
int alloc_arr(Ctx* c) {
    int r = ARR_VREG + c->arr_top;
    c->arr_top++;
    if (c->arr_top > c->arr_max) c->arr_max = c->arr_top;
    return r;
}

/* A strip-mining tile register.  Like array registers these are allocated into
 * a virtual range and relocated into their own contiguous bank at finalize, so
 * "is this value a tile?" is answered by the register number alone — no extra
 * field on Val that a construction site could forget to initialise. */
int alloc_tile(Ctx* c) {
    int r = TILE_VREG + c->tile_top;
    c->tile_top++;
    if (c->tile_top > c->tile_max) c->tile_max = c->tile_top;
    return r;
}
static bool reg_is_tile(int r) {
    return (uint32_t)r >= (uint32_t)TILE_VREG && (uint32_t)r < (uint32_t)ARR_VREG;
}
bool val_is_tile(Val v) { return reg_is_tile(v.reg); }

/* Pop a temporary WITHOUT emitting anything: for operands whose array (if any)
 * the consuming instruction frees itself via its AF_FREE_* flags, so the free
 * happens after the operand has been read. */
void pop_tmp(Ctx* c, Val v) {
    if (!v.tmp) return;
    /* An owned association temp (B3) lives in the array bank exactly like an
     * array temp — a borrowed-argument association is tmp==false and never
     * reaches here — so both pop the array bank. */
    if (CT_IS_ARRAY(v.type) || CT_IS_ASSOC(v.type)) c->arr_top--;
    else if (val_is_tile(v)) c->tile_top--;
    else c->temp_top--;
}
/* Pop a temporary whose value is now DEAD.  An array temp needs its handle
 * released here and now — inside a loop body the alternative (relying on
 * teardown) would accumulate one buffer per iteration. */
void free_if_tmp(Ctx* c, Val v) {
    if (!v.tmp) return;
    /* Owned association temps free through the same OP_ARR_FREE as array temps
     * (expr_free handles any Expr in the .arr slot). */
    if (CT_IS_ARRAY(v.type) || CT_IS_ASSOC(v.type)) {
        Slot z = { 0 };
        ins(c, OP_ARR_FREE, (uint32_t)v.reg, 0, 0, z);
        c->arr_top--;
    } else if (val_is_tile(v)) {
        c->tile_top--;      /* tile storage belongs to the frame, never freed */
    } else c->temp_top--;
}

/* Common type of two operands.  An array absorbs a scalar (broadcast) and two
 * arrays must agree on rank; element types widen exactly as scalars do. */
CompileType num_common(CompileType a, CompileType b) {
    if (a == CT_BOOL || b == CT_BOOL) return CT_ERR;
    /* An Association bag is never a numeric operand; without this it would fall
     * through to `a > b ? a : b` and be returned as a bogus "wider scalar". */
    if (CT_IS_ASSOC(a) || CT_IS_ASSOC(b)) return CT_ERR;
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
void coerce(Ctx* c, Val* v, CompileType target) {
    if (!c->ok || v->type == target) return;
    if (CT_IS_ARRAY(v->type) || CT_IS_ARRAY(target)) { c->ok = false; return; }
    if (CT_IS_ASSOC(v->type) || CT_IS_ASSOC(target)) { c->ok = false; return; }
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

bool emit(Ctx* c, const Expr* e, Val* out);
int cse_lookup(const Ctx* c, const Expr* e);

/* Single choke point protecting every SCALAR opcode from an array operand.
 * Each scalar op is emitted through binop / unop / kern_unop / kern_binop, so
 * one guard here is enough: any head that has no array lowering (comparisons,
 * Max/Min, Mod, ...) bails automatically the moment an array reaches it,
 * instead of silently reinterpreting a handle as a double. */
bool scalar_only(Ctx* c, CompileType a, CompileType b, CompileType r) {
    if (CT_IS_ARRAY(a) || CT_IS_ARRAY(b) || CT_IS_ARRAY(r)) { c->ok = false; return false; }
    /* Same guard for an Association handle: it is not a machine scalar, so any
     * scalar op that receives one must bail rather than read it as a double. */
    if (CT_IS_ASSOC(a) || CT_IS_ASSOC(b) || CT_IS_ASSOC(r)) { c->ok = false; return false; }
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
Val vsplat(Ctx* c, Val v) {
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
Val vec_unop(Ctx* c, uint16_t op, Val a, CompileType rtype, Slot imm) {
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
Val binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rtype) {
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
Val unop(Ctx* c, uint16_t op, Val a, CompileType rtype) {
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
Val emit_const(Ctx* c, Slot imm, CompileType type) {
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
Val kern_unop(Ctx* c, uint16_t op, Val a, CompileType rt, const void* fn) {
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
Val kern_binop(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, const void* fn) {
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
Val arr_op(Ctx* c, uint16_t op, Val a, Val b, CompileType rt, Slot imm) {
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
Val arr_noop_val(void) { Val v = { 0, false, CT_REAL, false }; return v; }

/* Prepare one operand of an array op: arrays pass through untouched (the ND
 * layer promotes element dtypes itself), scalars widen to Real/Complex so the
 * VM knows which half of the slot to read. */
void arr_prep(Ctx* c, Val* v, CompileType elem) {
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
Val arr_real_const(Ctx* c, double x) { Slot s; s.r = x; return emit_const(c, s, CT_REAL); }

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

bool infer_type(Ctx* c, const Expr* e, CompileType* out);

/* A machine element type an array buffer can hold (a nested array is not an
 * element).  Int became one when NDT_INT64 arrived; Bool when NDT_BOOL did, so a
 * body that computes a truth value (Table[Not[v[[i]]], ...], ConstantArray[True,
 * n]) now builds a one-byte bool buffer instead of refusing to lower. */
bool ct_is_elem(CompileType t) {
    return t == CT_INT || t == CT_REAL || t == CT_COMPLEX || t == CT_BOOL;
}

/* Is `node` (by identity) somewhere inside `root`?  Used only to keep a bail
 * diagnostic from pointing into a tree the emitter built and is about to free. */
bool expr_subtree_of(const Expr* root, const Expr* node) {
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
bool reg_is_owned_arr(int r) { return (uint32_t)r >= (uint32_t)ARR_VREG; }

void compile_partspec_free(PartSpec* p) {
    if (!p) return;
    for (int i = 0; i < p->n; i++) if (p->lit) expr_free(p->lit[i]);
    free(p->lit); free(p->reg); free(p);
}

/* Hand a freshly built PartSpec to the context, which owns it from here on. */
bool ctx_own_partspec(Ctx* c, PartSpec* p) {
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
bool subscript_is_literal_spec(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_All;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_Span || h == SYM_List;
}

/* True when `Part[a, A[0..na-1]]` can use the inline scalar path: one scalar
 * integer subscript per axis, none of them a slice spec. */
bool part_is_scalar_indexed(Ctx* c, CompileType at, const Expr* const* A, size_t na) {
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
bool emit_flat_index(Ctx* c, Val arr, const Expr* const* A, size_t na, int* idx_out) {
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
bool const_array_shape(Ctx* c, const Expr* const* A, int* rank_out, CompileType* elem_out) {
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
PartSpec* emit_partspec(Ctx* c, const Expr* const* A, size_t na, int* base_out) {
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
Val arr_ew(Ctx* c, Val a, Val b, CompileType rt, bool is_plus) {
    arr_prep(c, &a, CT_ELEM(rt)); arr_prep(c, &b, CT_ELEM(rt));
    Slot ip; ip.i = is_plus ? 1 : 0;
    return arr_op(c, OP_V_EW, a, b, rt, ip);
}

/* numeric literal? -> value + type */
bool literal(const Expr* e, Slot* imm, CompileType* type) {
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

bool named_const(const char* nm, double* out) {
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
bool unary_math(const char* h, uint16_t* op_r, uint16_t* op_c) {
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
/* The element type `ndarray_map_unary` will ACTUALLY produce for kernel `k`
 * over an operand whose element type is `ea`, or CT_ERR if it will decline.
 *
 * This has to be exact, not conservative.  The declared type is what every
 * downstream opcode reads the destination slot as, so a kernel that writes an
 * NDT_INT64 buffer under a CT_REAL declaration does not merely lose precision
 * — the consumer reinterprets the integer bits as a double.  That was live:
 *
 *     Compile[{{v, _Real, 1}}, Total[Floor[v]]][{1.5, 2.5, 3.5}]
 *
 * answered 2.96439*10^-323 (the int64 6, read as an IEEE double) where the
 * interpreter answers 6, and `Floor[v] + 1` came back as Reals where the
 * interpreter gives exact Integers.  Standalone `Floor[v]` was right, because
 * the caller re-reads the result buffer's real dtype — only a CONSUMER inside
 * the same program saw the lie.  Found by tools/compile_coverage.py on
 * 2026-08-02, which asked which heads compile over an array and made the
 * narrowing family worth looking at.
 *
 * The condition mirrors ndarray_map_unary (src/ndarray.c) line for line; the
 * two must not drift. */
CompileType nd_unary_elem(const NDUnaryKernel* k, CompileType ea) {
    if (!k) return CT_ERR;
    if (k->to_int && ea != CT_COMPLEX
        && (ea == CT_INT ? k->to_int_i != NULL : k->to_int_r != NULL))
        return CT_INT;
    /* Neither arm is the degrade sentinel: no machine kernel at all. */
    if (!k->cplx && !k->real) return CT_ERR;
    if (ea == CT_COMPLEX && !k->cplx) return CT_ERR;
    return k->to_real ? CT_REAL : ea;
}

/* The binary twin of nd_unary_elem: the element type ndarray_map_binary will
 * produce for kernel `k` over an array of element type `ea` broadcast against a
 * scalar of type `es`, or CT_ERR if it will decline.
 *
 * `es == CT_INT` is exactly the compiler's way of knowing what that function
 * calls `scal->type == EXPR_INTEGER` -- an EXACT integer, which is what decides
 * between the exact int64 arm and the double one.
 *
 * Without this the array form of every to_int-only binary kernel was refused
 * outright, because `!k->cplx` reads as the degrade sentinel and Mod, Quotient,
 * GCD, LCM, IntegerLength and IntegerExponent have no complex arm at all.  Each
 * has both a working interpreter buffer path AND a compiled SCALAR opcode, so
 * `Mod[x, 3]` compiled and `Mod[v, 3]` did not -- and being outside the subset
 * is a cliff, so the array spelling took its whole body to the interpreter.
 * Mirrors ndarray_map_binary (src/ndarray.c); the two must not drift. */
CompileType nd_binary_elem(const NDBinaryKernel* k, CompileType ea,
                                  CompileType es) {
    if (!k) return CT_ERR;
    if (k->to_int && ea != CT_COMPLEX && es != CT_COMPLEX) {
        bool use_i = (ea == CT_INT) && (es == CT_INT) && k->to_int_i != NULL;
        bool use_r = !use_i && k->to_int_r != NULL && ea != CT_INT;
        if (use_i || use_r) return CT_INT;
    }
    if (!k->cplx) return CT_ERR;                  /* sentinel: no machine kernel */
    return k->real_closed ? ea : CT_COMPLEX;
}

bool emit_arr_unary(Ctx* c, const char* head, Val a, Val* out) {
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
    CompileType er = nd_unary_elem(k, ea);
    if ((int)er < 0) { c->ok = false; return false; }
    Slot z; z.p = k;
    *out = arr_op(c, OP_V_KERN, a, arr_noop_val(), CT_ARRAY(er, rank), z);
    return c->ok;
}

/* emit a numeric unary function whose real/complex opcodes are op_r/op_c. */
bool emit_unary_math(Ctx* c, const char* head, const Expr* arg,
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

/* Resolve `f` as a function of exactly `want_arity` arguments.  Pure: it
 * inspects the tree and never emits, so both passes can call it. */
bool fn_resolve(const Expr* f, int want_arity, FnSpec* out) {
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

/* Replace the single implicit slot (# = Slot[1] / Slot[]) of a one-argument pure
 * function with the named parameter `pn`, so an FN_SLOTS body compiles as a
 * standalone function of one machine argument.  Returns NULL (bail) if the body
 * references a higher slot (#2, …) or a named slot — i.e. more than one param. */
static Expr* subst_slot1(const Expr* e, const char* pn) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Slot) {
        size_t na = e->data.function.arg_count;
        if (na == 0) return expr_new_symbol(pn);                       /* #        */
        if (na == 1 && e->data.function.args[0]->type == EXPR_INTEGER
            && e->data.function.args[0]->data.integer == 1) return expr_new_symbol(pn);  /* #1 */
        return NULL;                                                   /* #2+ / #x */
    }
    if (e->type == EXPR_FUNCTION) {
        Expr* head = subst_slot1(e->data.function.head, pn);
        if (!head) return NULL;
        size_t na = e->data.function.arg_count;
        Expr** args = malloc((na ? na : 1) * sizeof(Expr*));
        if (!args) { expr_free(head); return NULL; }
        for (size_t i = 0; i < na; i++) {
            args[i] = subst_slot1(e->data.function.args[i], pn);
            if (!args[i]) { for (size_t j = 0; j < i; j++) expr_free(args[j]); free(args); expr_free(head); return NULL; }
        }
        Expr* r = expr_new_function(head, args, na);
        free(args);
        return r;
    }
    return expr_copy((Expr*)e);
}

/* Compile an embedded one-argument function `fn` — one machine-scalar value of
 * type `in` in, one machine result out — into a self-contained callee program, so
 * a B4 kernel can invoke it per entry with vm_call (no evaluator).  Handles the
 * three function forms: `Function[u, body]` (compile body of u), `#…&`
 * (substitute the slot for a named param), and a bare head `h` (compile
 * h[param]); Identity is the trivial pass-through.  Returns NULL if the function
 * shape is unsupported or its body is outside the compilable subset (the caller
 * then bails).  Inherits the parent's global-folding flag; the callee is
 * throwaway (freed with the program), so folding is lifetime-safe. */
struct CompiledProgram* compile_value_callee(Ctx* c, const Expr* fn, CompileType in) {
    FnSpec s;
    if (!fn_resolve(fn, 1, &s)) return NULL;
    const char* pn = intern_symbol("$b4value");
    const char* names[1];
    CompileType types[1] = { in };
    Expr* body_owned = NULL;
    const Expr* body;
    switch (s.kind) {
        case FN_IDENTITY:
            names[0] = pn; body_owned = expr_new_symbol(pn); body = body_owned; break;
        case FN_HEAD: {
            names[0] = pn;
            Expr* argsym = expr_new_symbol(pn);
            Expr* aargs[1] = { argsym };
            body_owned = expr_new_function(expr_new_symbol(s.head), aargs, 1);   /* head[param] */
            body = body_owned; break;
        }
        case FN_SLOTS:
            names[0] = pn;
            body_owned = subst_slot1(s.body, pn);
            if (!body_owned) return NULL;
            body = body_owned; break;
        case FN_LAMBDA:
            if (s.nparams != 1) return NULL;
            names[0] = s.pname[0]; body = s.body; break;   /* body borrowed from fn */
        default: return NULL;   /* FN_COMPOSE, multi-param, … */
    }
    struct CompiledProgram* cp =
        compile_expr_ex(body, names, types, 1, c->flags & COMPILE_FOLD_GLOBALS);
    if (body_owned) expr_free(body_owned);
    return cp;
}

/* `Slot[k]` (`#`, `#k`) resolved against the live slot frame, or -1.  Shared by
 * both passes so they cannot disagree about which slots are bound. */
int fn_slot_index(const Ctx* c, Expr* const* A, size_t na) {
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
const char* fn_placeholder(int i) {
    static const char* cache[FN_MAX_PARAMS];
    static const char* const names[FN_MAX_PARAMS] = {
        "System`Compile$fn1", "System`Compile$fn2", "System`Compile$fn3", "System`Compile$fn4",
        "System`Compile$fn5", "System`Compile$fn6", "System`Compile$fn7", "System`Compile$fn8"
    };
    if (!cache[i]) cache[i] = intern_symbol(names[i]);
    return cache[i];
}

/* Build `head[$1, ..., $n]` over the reserved placeholder symbols.  Caller frees. */
Expr* fn_head_call(const char* head, int n) {
    Expr* args[FN_MAX_PARAMS];
    for (int i = 0; i < n; i++) {
        args[i] = expr_new_symbol(fn_placeholder(i));
        if (!args[i]) { while (i-- > 0) expr_free(args[i]); return NULL; }
    }
    Expr* hd = expr_new_symbol(head);
    if (!hd) { for (int i = 0; i < n; i++) expr_free(args[i]); return NULL; }
    return expr_new_function(hd, args, (size_t)n);   /* adopts hd and args */
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

bool loop_spec_parse(const Expr* spec, LoopSpec* out) {
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

bool infer_type(Ctx* c, const Expr* e, CompileType* out);

/* Both bounds must be integer-typed: these are integer-counted loops. */
bool loop_spec_int_bounds(Ctx* c, const LoopSpec* s) {
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
bool range_spec(Ctx* c, Expr* const* A, size_t na, LoopSpec* out) {
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
bool emit_loop_bound(Ctx* c, const Expr* b, int dst) {
    Slot z = { 0 };
    if (!b) { Slot one; one.i = 1; ins(c, OP_CONST, (uint32_t)dst, 0, 0, one); return c->ok; }
    Val v; if (!emit(c, b, &v)) return false;
    ins(c, OP_MOVE, (uint32_t)dst, (uint32_t)v.reg, 0, z);
    free_if_tmp(c, v);
    return c->ok;
}

/* Element access opcode for an array's element type.  One place rather than the
 * fourteen `elem == CT_COMPLEX ? _C : _R` conditionals this replaces — adding
 * the integer element type to those by hand is exactly the kind of edit that
 * gets 13 of 14 sites. */
uint16_t a_load_op(CompileType elem) {
    return elem == CT_COMPLEX ? OP_A_LOAD_C
         : elem == CT_INT     ? OP_A_LOAD_I
         : elem == CT_BOOL    ? OP_A_LOAD_B
                              : OP_A_LOAD_R;
}
uint16_t a_store_op(CompileType elem) {
    return elem == CT_COMPLEX ? OP_A_STORE_C
         : elem == CT_INT     ? OP_A_STORE_I
         : elem == CT_BOOL    ? OP_A_STORE_B
                              : OP_A_STORE_R;
}


/* True when any argument is (or infers to) an array — the signal to take an
 * array lowering rather than the scalar one.  This is only a routing hint: a
 * wrong "false" still ends in a clean bail, because every scalar opcode is
 * guarded by scalar_only(). */
bool any_array_arg(Ctx* c, Expr** A, size_t na) {
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
bool try_kernel(Ctx* c, const char* h, Expr** A, size_t na, Val* out) {
    if (na == 1) {
        SymbolDef* d = symtab_lookup(h);
        if (!d || !d->ndarray_unary_kernel) return false;
        const NDUnaryKernel* k = (const NDUnaryKernel*)d->ndarray_unary_kernel;
        /* No real and no complex arm is the degrade sentinel -- unless the
         * kernel is exact-integer-only (to_int), which is a real path over an
         * ARRAY and no path at all over a scalar.  See emit_arr_unary.
         *
         * The test is `to_int`, not `to_int_r`: the integer-only kernels
         * (MoebiusMu, EulerPhi, IntegerLength) have `to_int_i` and NOTHING
         * else, since MoebiusMu of a real is not a machine question -- so
         * keying on the real arm reads exactly the heads this is for as
         * sentinels.  nd_unary_elem picks the arm per element type. */
        bool narrowing_only = (!k->cplx && !k->real && k->to_int);
        if (!k->cplx && !k->real && !narrowing_only) return false;
        Val a; if (!emit(c, A[0], &a)) return false;
        if (CT_IS_ARRAY(a.type)) return emit_arr_unary(c, h, a, out);
        if (narrowing_only) return false;   /* scalar: leave it to a dedicated opcode */
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
        /* No complex arm is the degrade sentinel for a SCALAR body -- but not
         * for an array one, where the exact int64 arms are a real path.  The
         * cheap static test first, so nothing else changes for the kernels
         * that have a complex arm. */
        if (!k->cplx && !k->to_int) return false;
        Val a, b; if (!emit(c, A[0], &a) || !emit(c, A[1], &b)) return false;
        CompileType t = num_common(a.type, b.type); if ((int)t < 0) { c->ok = false; return false; }
        if (CT_IS_ARRAY(t)) {
            /* one array + one broadcast scalar (BesselJ[n,v], ArcTan[v,y], ...);
             * the ND layer has no array-array binary-kernel map. */
            if (CT_IS_ARRAY(a.type) && CT_IS_ARRAY(b.type)) { c->ok = false; return false; }
            CompileType ea = CT_IS_ARRAY(a.type) ? CT_ELEM(a.type) : CT_ELEM(b.type);
            CompileType es = CT_IS_ARRAY(a.type) ? b.type : a.type;
            CompileType er = nd_binary_elem(k, ea, es);
            if ((int)er < 0) { c->ok = false; return false; }
            /* Prepare at the COMMON element type, not the result's and not the
             * array's.  The result's is wrong for a narrowing kernel (CT_INT
             * out of a Real buffer would round the scalar before the kernel
             * ran); the array's is wrong for a complex scalar against a real
             * buffer (`BesselJ[0.5 + I, v]`), where coerce refuses to narrow
             * and the whole body would bail. */
            CompileType ep = CT_ELEM(t);
            arr_prep(c, &a, ep); arr_prep(c, &b, ep);
            Slot z; z.p = k;
            *out = arr_op(c, OP_V_KERN2, a, b, CT_ARRAY(er, CT_RANK(t)), z);
            return c->ok;
        }
        if (!k->cplx) return false;   /* scalar body, sentinel kernel -> bail */
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

EmitMark emit_mark(const Ctx* c) {
    EmitMark m;
    m.n = c->n; m.temp_top = c->temp_top; m.maxreg = c->maxreg;
    m.arr_top = c->arr_top; m.arr_max = c->arr_max; m.nscope = c->nscope; m.ok = c->ok;
    m.tile_top = c->tile_top; m.tile_max = c->tile_max; m.vector_mode = c->vector_mode;
    /* A speculative lowering that fails is not a bail — the caller falls through
     * to another strategy — so its diagnostic must not survive the rollback. */
    m.bail_node = c->bail_node;
    return m;
}
void emit_rollback(Ctx* c, EmitMark m) {
    c->n = m.n; c->temp_top = m.temp_top; c->maxreg = m.maxreg;
    c->arr_top = m.arr_top; c->arr_max = m.arr_max; c->nscope = m.nscope; c->ok = m.ok;
    c->tile_top = m.tile_top; c->tile_max = m.tile_max; c->vector_mode = m.vector_mode;
    c->bail_node = m.bail_node;
}


/* Collect the distinct array-valued symbols in `e`.  Returns false if the tree
 * contains more than FUSE_MAX_LEAVES of them. */
bool fuse_collect(Ctx* c, const Expr* e, FuseLeaves* L) {
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
int fuse_listable(Ctx* c, const Expr* e) {
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

bool emit(Ctx* c, const Expr* e, Val* out);

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
bool emit_apply(Ctx* c, const FnSpec* s, const Val* argv, int n, Val* out) {
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
void emit_max_guard(Ctx* c, int reg, long long lim) {
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
bool emit_iter_count(Ctx* c, const LoopSpec* s, int rlo, int dst) {
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
bool stmt_valued_head(const Expr* e) {
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
uint16_t emit_sameq_op(CompileType t) {
    return t == CT_INT ? OP_EQ_I : t == CT_COMPLEX ? OP_SAMEQ_C : OP_SAMEQ_R;
}

/* `if (R[reg] < 0) fail` — a negative application bound leaves Nest and
 * FixedPoint UNEVALUATED in the interpreter (src/funcprog.c:2153), where a
 * counted loop would silently run zero times and return the seed. */
void emit_nonneg_guard(Ctx* c, int reg) {
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
void emit_nonzero_guard(Ctx* c, int reg) {
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

    /* Association read ops (B1), before fusion and the array lowerings: an
     * association-only head (Lookup/KeyExistsQ/...) is handled or bails here; a
     * head shared with the array path (Length/Values) is handled only when its
     * operand is an association, else try_emit_assoc returns false and the
     * ordinary lowering below takes it. */
    if (try_emit_assoc(c, h, A, na, out)) return c->ok;

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
    { int _r = emit_arith(c, h, e, A, na, out); if (_r) return _r > 0; }
    { int _r = emit_mathfn(c, h, e, A, na, out); if (_r) return _r > 0; }
    { int _r = emit_int(c, h, e, A, na, out); if (_r) return _r > 0; }
    { int _r = emit_array(c, h, e, A, na, out); if (_r) return _r > 0; }
    /* comparisons / Order / And-Or-Xor / Not — compile_emit_logic.c */
    { int _r = emit_logic(c, h, e, A, na, out); if (_r) return _r > 0; }

    { int _r = emit_ctrl(c, h, e, A, na, out); if (_r) return _r > 0; }
    { int _r = emit_iter(c, h, e, A, na, out); if (_r) return _r > 0; }

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
bool emit(Ctx* c, const Expr* e, Val* out) {
    bool r = emit_node(c, e, out);
    if (!r && !c->bail_node) c->bail_node = e;
    return r;
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
int cse_lookup(const Ctx* c, const Expr* e) {
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

static uint32_t patch_reg(uint32_t r, int arr_base, int tile_base, int managed_base) {
    if (r >= (uint32_t)ARR_VREG)  return (uint32_t)arr_base     + (r - (uint32_t)ARR_VREG);
    if (r >= (uint32_t)TILE_VREG) return (uint32_t)tile_base    + (r - (uint32_t)TILE_VREG);
    if (r >= (uint32_t)MGD_VREG)  return (uint32_t)managed_base + (r - (uint32_t)MGD_VREG);
    return r;
}

/* Free the emit-time managed state on the BAIL path (the literal containers have
 * not yet been transferred to a program).  On the success path only mgd_type is
 * freed here; the containers become the program's. */
static void ctx_free_managed(Ctx* c) {
    free(c->mgd_type); c->mgd_type = NULL;
    for (int i = 0; i < c->nmgd_consts; i++) mgd_const_free(&c->mgd_consts[i]);
    free(c->mgd_consts); c->mgd_consts = NULL;
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
    return compile_expr_prec(body, arg_names, arg_types, nargs, flags, 0);
}

CompiledProgram* compile_expr_prec(const Expr* body, const char* const* arg_names,
                                   const CompileType* arg_types, size_t nargs,
                                   unsigned flags, long prec_bits) {
    bail_clear();
    if (!body) { g_bail_reason = "empty body"; return NULL; }
    Ctx c; memset(&c, 0, sizeof(c));
    /* A managed body (arbitrary precision) is lowered by the separate emit_mgd
     * and must skip the optimiser, Expr-CSE, fusion and parallelism, so the
     * optimiser never sees a managed register and the machine path is unchanged. */
    bool managed = (prec_bits > 0) || (flags & COMPILE_BIGINT);
    if (managed) flags |= COMPILE_NO_OPT | COMPILE_NO_CSE | COMPILE_NO_FUSE | COMPILE_NO_PAR;
    c.prec_bits = prec_bits;
    c.allow_bigint = (flags & COMPILE_BIGINT) != 0;
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
    bool ok = (managed ? emit_mgd(&c, body, &res) : emit(&c, body, &res)) && c.ok;
    /* A borrowed argument array — or association (B3) — cannot be the result:
     * the caller owns what it gets back, and freeing an argument would corrupt
     * the caller's value.  A PRODUCED association (KeyDrop/KeyTake) is a tmp and
     * is allowed. */
    if (ok && (CT_IS_ARRAY(res.type) || CT_IS_ASSOC(res.type)) && !res.tmp) {
        ok = false;
        c.bail_node = body;   /* the whole body is a borrowed argument handle */
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
        for (int i = 0; i < c.nassocs; i++) compile_assocspec_free(c.assocs[i]);
        free(c.assocs);
        for (int i = 0; i < c.ncallees; i++) compile_assoccallee_free(c.callees[i]);
        free(c.callees);
        ctx_free_managed(&c);
        free(c.code); free(c.argdep); return NULL;
    }

    /* Four contiguous banks: scalars, then array handles, then strip-mining
     * tiles, then arbitrary-precision managed containers.  A slot therefore has
     * one kind for the whole life of the program, so teardown can never mistake a
     * double for a pointer, and each bank is a range sweep.  A machine program has
     * managed_base == nreg (empty managed bank). */
    int arr_base     = c.maxreg;
    int tile_base    = arr_base + c.arr_max;
    int managed_base = tile_base + c.tile_max;
    int nreg         = managed_base + c.mgd_max;
    for (size_t i = 0; i < c.n; i++) {
        c.code[i].dst = patch_reg(c.code[i].dst, arr_base, tile_base, managed_base);
        c.code[i].a   = patch_reg(c.code[i].a, arr_base, tile_base, managed_base);
        /* `b` is a branch TARGET on the jumping kinds and a register everywhere
         * else.  Asking the kind table rather than listing the opcodes means a
         * new branch opcode cannot have its target silently relocated into the
         * array bank — which is a corruption with no symptom until it jumps. */
        int bk = compile_op_kind[c.code[i].op];
        if (bk != K_JMP && bk != K_JZ && bk != K_LOOP && bk != K_APAR)
            c.code[i].b = patch_reg(c.code[i].b, arr_base, tile_base, managed_base);
    }
    int result_reg = (int)patch_reg((uint32_t)res.reg, arr_base, tile_base, managed_base);

    /* Optimise the emitted bytecode.  Runs after patch_reg so the array bank is
     * already at its final place and `arr_base` means what the optimiser expects.
     * Register numbers are preserved, so `result_reg` stays valid.  A failure
     * here is non-fatal: the unoptimised code is still correct.  Managed programs
     * force COMPILE_NO_OPT, so the optimiser never sees a managed register. */
    if (!(flags & COMPILE_NO_OPT)) compile_optimize(c.code, &c.n, nreg, arr_base, tile_base);

    CompiledProgram* p = calloc(1, sizeof(*p));
    if (!p) {
        for (int i = 0; i < c.nparts; i++) compile_partspec_free(c.parts[i]);
        free(c.parts);
        for (int i = 0; i < c.nassocs; i++) compile_assocspec_free(c.assocs[i]);
        free(c.assocs);
        for (int i = 0; i < c.ncallees; i++) compile_assoccallee_free(c.callees[i]);
        free(c.callees);
        ctx_free_managed(&c);
        free(c.code); free(c.argdep); return NULL;
    }
    /* The general-Part subscript lists become the program's: their literal specs
     * are pointed at from instruction immediates and must outlive the body. */
    p->parts = c.parts; p->nparts = c.nparts;
    p->assocs = c.assocs; p->nassocs = c.nassocs;   /* AssocSpecs, likewise (B1) */
    p->callees = c.callees; p->ncallees = c.ncallees;   /* B4 callees */
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

    /* Managed (arbitrary-precision) wiring.  A managed register needs a container
     * allocated at frame entry: enumerate every one — the managed ARGUMENT
     * registers (which live in the scalar-bank arg range) and every register of
     * the managed bank — with its container type, and hand the program its
     * literal containers.  Zero of this runs for a machine program. */
    p->prec_bits = prec_bits;
    p->managed_base = managed_base;
    p->mgd_consts = c.mgd_consts; p->nmgd_consts = c.nmgd_consts;  /* ownership moves */
    {
        int nmanaged_args = 0;
        for (size_t k = 0; k < nargs; k++) if (CT_IS_MANAGED(arg_types[k])) nmanaged_args++;
        int total = nmanaged_args + c.mgd_max;
        if (total > 0) {
            p->mgd_slots = malloc((size_t)total * sizeof(MgdSlot));
            if (!p->mgd_slots) { free(c.mgd_type); compiled_free(p); return NULL; }
            int j = 0;
            for (size_t k = 0; k < nargs; k++)
                if (CT_IS_MANAGED(arg_types[k])) {
                    p->mgd_slots[j].reg = (int)k;
                    p->mgd_slots[j].type = arg_types[k];
                    j++;
                }
            for (int i = 0; i < c.mgd_max; i++) {
                p->mgd_slots[j].reg = managed_base + i;
                p->mgd_slots[j].type = c.mgd_type[i];
                j++;
            }
            p->nmgd_slots = total;
        }
    }
    free(c.mgd_type);   /* types now live in p->mgd_slots */

    extract_par_loops(p);
    return p;
}

