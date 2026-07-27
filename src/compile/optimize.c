/* Mathilda — Compile[] engine: bytecode optimiser.
 *
 * A bytecode-to-bytecode transform run once, at the end of compilation, between
 * the emitter and the finished program.  The emitter is deliberately naive — it
 * lowers each Expr node independently — so the redundancy it leaves behind is
 * regular and cheap to remove here rather than complicating every lowering.
 *
 * Passes, in order:
 *   1. value numbering (per basic block) — constant folding, CSE, copy
 *      propagation, algebraic identities
 *   2. dead-code elimination, driven by real liveness over the CFG
 *   3. loop-invariant code motion
 *   4. compaction — drop the NOPs left by 1-3 and remap jump targets
 * with 1 and 2 repeated after 3, since hoisting exposes more of both.
 *
 * Invariants every pass must preserve:
 *   - Array registers (>= arr_base) and the impure K_ARR opcodes are never
 *     touched: they allocate, free and can fail, and their ownership is encoded
 *     in instruction flags rather than in the dataflow.
 *   - Instruction *indices* may only change in the compaction/rebuild step,
 *     which remaps every jump target in the same pass.
 *   - Results must be bitwise identical.  Nothing here reassociates floating
 *     point or contracts a multiply-add; folding computes with exactly the
 *     operation the VM would have run.
 */
#include "compile_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* The kind table, generated from the same OPLIST as the opcode enum and the
 * VM's jump table — a new opcode cannot be added to one and missed by another. */
#define X(name, kind) kind,
const unsigned char compile_op_kind[OP__COUNT] = { OPLIST };
#undef X

/* ------------------------------------------------------------------ *
 *  Instruction shape                                                  *
 * ------------------------------------------------------------------ */
/* Which fields of an instruction are register reads, which is the register
 * write, and whether it is safe to remove or reorder. */
typedef struct {
    unsigned char wd;      /* writes R[dst] */
    unsigned char rd;      /* reads  R[dst] */
    unsigned char ra;      /* reads  R[a]   */
    unsigned char rb;      /* reads  R[b]   */
    unsigned char pure;    /* no side effects: removable when the result is dead */
    unsigned char jump;    /* `b` is a branch target, NOT a register */
} OpDesc;

static OpDesc op_desc(unsigned op) {
    OpDesc o;
    memset(&o, 0, sizeof o);
    switch (compile_op_kind[op]) {
        case K_NOP:                                                     break;
        case K_CONST:  o.wd = 1;                            o.pure = 1; break;
        case K_MOVE:   o.wd = 1;             o.ra = 1;      o.pure = 1; break;
        case K_UN:
        case K_KERN1:
        case K_POWI:   o.wd = 1;             o.ra = 1;      o.pure = 1; break;
        case K_BIN:
        case K_KERN2:  o.wd = 1;             o.ra = 1; o.rb = 1;
                                                            o.pure = 1; break;
        case K_INC:    o.wd = 1; o.rd = 1;                  o.pure = 1; break;
        case K_JMP:                                      o.jump = 1;    break;
        case K_JZ:                           o.ra = 1;   o.jump = 1;    break;
        /* increment + test + branch in one: reads and writes the counter,
         * reads the limit, and `b` is a target rather than a register. */
        case K_LOOP:   o.wd = 1; o.rd = 1;   o.ra = 1;   o.jump = 1;    break;
        case K_RET:              o.rd = 1;                              break;
        /* Array ops allocate, free, and can fail.  Treated as reading and
         * writing all three operands and never removable — the scalar passes
         * simply flow around them. */
        case K_ARR:    o.wd = 1; o.rd = 1;   o.ra = 1; o.rb = 1;        break;
        /* Writes array memory (or fails on a shape mismatch): impure, so never
         * removed or hoisted, but it owns nothing, so unlike K_ARR it does not
         * force the surrounding fused loop out of the optimiser's reach. */
        case K_ASTORE:           o.rd = 1;   o.ra = 1; o.rb = 1;        break;
        /* A call reads a RANGE of registers, which this descriptor cannot
         * express, so `ra` stays clear and the range is handled explicitly in
         * liveness and DCE.  Not marked pure: removing or duplicating a call
         * would change which programs fail, and a failing callee is what routes
         * the whole evaluation back to the interpreter. */
        case K_CALL:   o.wd = 1;                                        break;
        /* Strip-mined reduction step: accumulates into dst, so it READS dst as
         * well as writing it — which is also what keeps LICM from hoisting it
         * out of the tile loop. */
        case K_VACC:   o.wd = 1; o.rd = 1;   o.ra = 1;      o.pure = 1; break;
        default:                                                        break;
    }
    return o;
}

/* Immediates are only partly meaningful, and which part depends on the kind, so
 * value numbering compares them per kind rather than memcmp-ing the union (whose
 * unused bytes are indeterminate for everything except CONST, which the emitter
 * zero-fills for exactly this reason). */
static bool imm_eq(unsigned op, Slot x, Slot y) {
    switch (compile_op_kind[op]) {
        case K_CONST: return memcmp(&x, &y, sizeof x) == 0;
        case K_POWI:
        case K_INC:   return x.i == y.i;
        case K_KERN1:
        case K_KERN2: return x.p == y.p;
        default:      return true;             /* immediate carries no value */
    }
}

/* ------------------------------------------------------------------ *
 *  Control-flow graph                                                 *
 * ------------------------------------------------------------------ */
typedef struct {
    size_t   n;            /* instruction count */
    int      nreg;
    int      arr_base;     /* registers in [arr_base, tile_base) are array handles */
    int      tile_base;    /* registers >= this are tile pointers */
    Instr*   code;

    size_t   nblk;
    size_t*  blk_start;    /* [nblk] first instruction of each block */
    size_t*  blk_end;      /* [nblk] one past the last */
    int*     blk_of;       /* [n] block index of each instruction */
    size_t   succ[2];      /* scratch */

    size_t   nw;           /* bitset words per register set */
    uint64_t* live_in;     /* [nblk * nw] */
    uint64_t* live_out;    /* [nblk * nw] */
} Opt;

/* A register whose slot holds a POINTER rather than a number: an array handle or
 * a tile buffer.  The scalar passes must not fold, rewrite or duplicate these —
 * in particular CSE must never turn a tile-producing op into a MOVE, which would
 * copy the pointer and silently alias two tiles onto one buffer. */
static bool is_ptr_reg(const Opt* o, int r) { return r >= o->arr_base; }
static bool is_tile_reg(const Opt* o, int r) { return r >= o->tile_base; }

static void bs_set(uint64_t* s, int r)   { s[r >> 6] |= (uint64_t)1 << (r & 63); }
static void bs_clr(uint64_t* s, int r)   { s[r >> 6] &= ~((uint64_t)1 << (r & 63)); }
static bool bs_get(const uint64_t* s, int r) { return (s[r >> 6] >> (r & 63)) & 1u; }

/* Successors of the block ending at instruction `e`.  Returns the count and
 * fills o->succ. */
static int succ_of(Opt* o, size_t e) {
    unsigned op = o->code[e].op;
    int k = compile_op_kind[op];
    if (k == K_RET) return 0;
    if (k == K_JMP) { o->succ[0] = o->code[e].b; return 1; }
    if (k == K_JZ || k == K_LOOP) { o->succ[0] = o->code[e].b; o->succ[1] = e + 1;
                      return (e + 1 < o->n) ? 2 : 1; }
    if (e + 1 < o->n) { o->succ[0] = e + 1; return 1; }
    return 0;
}

static bool build_cfg(Opt* o) {
    unsigned char* leader = calloc(o->n + 1, 1);
    o->blk_of = malloc(o->n * sizeof(int));
    if (!leader || !o->blk_of) { free(leader); return false; }

    leader[0] = 1;
    for (size_t i = 0; i < o->n; i++) {
        int k = compile_op_kind[o->code[i].op];
        if (k == K_JMP || k == K_JZ || k == K_LOOP) {
            size_t t = o->code[i].b;
            if (t < o->n) leader[t] = 1;
            if (i + 1 < o->n) leader[i + 1] = 1;
        } else if (k == K_RET && i + 1 < o->n) {
            leader[i + 1] = 1;
        }
    }

    o->nblk = 0;
    for (size_t i = 0; i < o->n; i++) if (leader[i]) o->nblk++;
    o->blk_start = malloc(o->nblk * sizeof(size_t));
    o->blk_end   = malloc(o->nblk * sizeof(size_t));
    if (!o->blk_start || !o->blk_end) { free(leader); return false; }

    size_t b = 0;
    for (size_t i = 0; i < o->n; i++) {
        if (leader[i]) { if (b) o->blk_end[b - 1] = i; o->blk_start[b++] = i; }
        o->blk_of[i] = (int)(b - 1);
    }
    o->blk_end[o->nblk - 1] = o->n;
    free(leader);
    return true;
}

/* Backward liveness over the CFG, at register granularity. */
static bool liveness(Opt* o) {
    o->nw = (size_t)((o->nreg + 63) / 64);
    if (o->nw == 0) o->nw = 1;
    size_t words = o->nblk * o->nw;
    uint64_t* use = calloc(words, sizeof(uint64_t));
    uint64_t* def = calloc(words, sizeof(uint64_t));
    o->live_in  = calloc(words, sizeof(uint64_t));
    o->live_out = calloc(words, sizeof(uint64_t));
    if (!use || !def || !o->live_in || !o->live_out) {
        free(use); free(def); return false;
    }

    /* Per block: `use` = read before any write in the block; `def` = written. */
    for (size_t b = 0; b < o->nblk; b++) {
        uint64_t* u = use + b * o->nw;
        uint64_t* d = def + b * o->nw;
        for (size_t i = o->blk_start[b]; i < o->blk_end[b]; i++) {
            const Instr* c = &o->code[i];
            OpDesc s = op_desc(c->op);
            if (s.rd && !bs_get(d, (int)c->dst)) bs_set(u, (int)c->dst);
            if (s.ra && !bs_get(d, (int)c->a))   bs_set(u, (int)c->a);
            if (s.rb && !s.jump && !bs_get(d, (int)c->b)) bs_set(u, (int)c->b);
            if (compile_op_kind[c->op] == K_CALL)
                for (unsigned k = 0; k < c->flags; k++)
                    if (!bs_get(d, (int)c->a + (int)k)) bs_set(u, (int)c->a + (int)k);
            if (s.wd) bs_set(d, (int)c->dst);
        }
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t bi = o->nblk; bi-- > 0; ) {
            uint64_t* out = o->live_out + bi * o->nw;
            int ns = succ_of(o, o->blk_end[bi] - 1);
            for (int s = 0; s < ns; s++) {
                size_t sb = (size_t)o->blk_of[o->succ[s]];
                const uint64_t* si = o->live_in + sb * o->nw;
                for (size_t w = 0; w < o->nw; w++)
                    if ((out[w] | si[w]) != out[w]) { out[w] |= si[w]; changed = true; }
            }
            uint64_t* in = o->live_in + bi * o->nw;
            const uint64_t* u = use + bi * o->nw;
            const uint64_t* d = def + bi * o->nw;
            for (size_t w = 0; w < o->nw; w++) {
                uint64_t v = u[w] | (out[w] & ~d[w]);
                if (v != in[w]) { in[w] = v; changed = true; }
            }
        }
    }
    free(use); free(def);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Pass 1: value numbering (fold / CSE / copy propagation)            *
 * ------------------------------------------------------------------ */
typedef struct { unsigned op; uint32_t a, b; Slot imm; uint32_t dst; bool valid; } VNEntry;

/* Constant-fold `op` over known operand constants.  Computes with exactly the
 * expression the VM would have executed, so the folded value is bit-identical.
 * Returns false for anything not folded (always safe). */
static bool fold_op(unsigned op, Slot a, Slot b, Slot imm, Slot* out) {
    Slot r; memset(&r, 0, sizeof r);
    switch (op) {
        case OP_I2R: r.r = (double)a.i; break;
        case OP_I2C: r.z = (double)a.i; break;
        case OP_R2C: r.z = a.r;         break;
        case OP_ADD_I: r.i = a.i + b.i; break;
        case OP_SUB_I: r.i = a.i - b.i; break;
        case OP_MUL_I: r.i = a.i * b.i; break;
        case OP_ADD_R: r.r = a.r + b.r; break;
        case OP_SUB_R: r.r = a.r - b.r; break;
        case OP_MUL_R: r.r = a.r * b.r; break;
        case OP_DIV_R: r.r = a.r / b.r; break;
        case OP_ADD_C: r.z = a.z + b.z; break;
        case OP_SUB_C: r.z = a.z - b.z; break;
        case OP_MUL_C: r.z = a.z * b.z; break;
        case OP_DIV_C: r.z = a.z / b.z; break;
        case OP_NEG_I: r.i = -a.i;      break;
        case OP_NEG_R: r.r = -a.r;      break;
        case OP_NEG_C: r.z = -a.z;      break;
        case OP_INV_R: r.r = 1.0 / a.r; break;
        case OP_INV_C: r.z = 1.0 / a.z; break;
        case OP_SQRT_R: r.r = sqrt(a.r); break;
        case OP_EXP_R:  r.r = exp(a.r);  break;
        case OP_LOG_R:  r.r = log(a.r);  break;
        case OP_SIN_R:  r.r = sin(a.r);  break;
        case OP_COS_R:  r.r = cos(a.r);  break;
        case OP_ABS_I:  r.i = a.i < 0 ? -a.i : a.i; break;
        case OP_ABS_R:  r.r = fabs(a.r); break;
        case OP_LT_I: r.i = a.i <  b.i; break;
        case OP_LE_I: r.i = a.i <= b.i; break;
        case OP_GT_I: r.i = a.i >  b.i; break;
        case OP_GE_I: r.i = a.i >= b.i; break;
        case OP_LT_R: r.i = a.r <  b.r; break;
        case OP_LE_R: r.i = a.r <= b.r; break;
        case OP_GT_R: r.i = a.r >  b.r; break;
        case OP_GE_R: r.i = a.r >= b.r; break;
        case OP_NOT:  r.i = !a.i;       break;
        case OP_AND:  r.i = a.i && b.i; break;
        case OP_OR:   r.i = a.i || b.i; break;
        case OP_MAX_R: r.r = a.r > b.r ? a.r : b.r; break;
        case OP_MIN_R: r.r = a.r < b.r ? a.r : b.r; break;
        case OP_MAX_I: r.i = a.i > b.i ? a.i : b.i; break;
        case OP_MIN_I: r.i = a.i < b.i ? a.i : b.i; break;
        case OP_POWI_R: {                    /* mirrors the VM's ipow_r exactly */
            long long e = imm.i;
            unsigned long long m = (unsigned long long)(e < 0 ? -e : e);
            double base = a.r, acc = 1.0;
            while (m) { if (m & 1u) acc *= base; base *= base; m >>= 1; }
            r.r = (e < 0) ? 1.0 / acc : acc;
            break;
        }
        default: return false;
    }
    *out = r;
    return true;
}

static bool pass_vn(Opt* o, bool* progress) {
    int nreg = o->nreg;
    unsigned char* kk = calloc((size_t)nreg, 1);      /* register holds a constant */
    Slot*  kv  = calloc((size_t)nreg, sizeof(Slot));
    uint32_t* alias = malloc((size_t)nreg * sizeof(uint32_t));
    VNEntry* vn = calloc(o->n ? o->n : 1, sizeof(VNEntry));
    if (!kk || !kv || !alias || !vn) { free(kk); free(kv); free(alias); free(vn); return false; }

    for (size_t b = 0; b < o->nblk; b++) {
        memset(kk, 0, (size_t)nreg);
        for (int r = 0; r < nreg; r++) alias[r] = (uint32_t)r;
        size_t nvn = 0;

        for (size_t i = o->blk_start[b]; i < o->blk_end[b]; i++) {
            Instr* c = &o->code[i];
            OpDesc s = op_desc(c->op);
            if (compile_op_kind[c->op] == K_ARR) {          /* opaque: reset state */
                memset(kk, 0, (size_t)nreg);
                for (int r = 0; r < nreg; r++) alias[r] = (uint32_t)r;
                nvn = 0;
                continue;
            }

            /* copy propagation: read through any alias established earlier */
            if (s.ra && (int)c->a < nreg && c->a != alias[c->a]) { c->a = alias[c->a]; *progress = true; }
            if (s.rb && !s.jump && (int)c->b < nreg && c->b != alias[c->b]) { c->b = alias[c->b]; *progress = true; }

            /* constant folding */
            if (s.pure && s.wd && !is_ptr_reg(o, (int)c->dst)) {
                bool ka = !s.ra || ((int)c->a < nreg && kk[c->a]);
                bool kb = !s.rb || ((int)c->b < nreg && kk[c->b]);
                if (compile_op_kind[c->op] != K_CONST && s.ra && ka && kb && !s.rd) {
                    Slot va = s.ra ? kv[c->a] : (Slot){ 0 };
                    Slot vb = s.rb ? kv[c->b] : (Slot){ 0 };
                    Slot folded;
                    if (fold_op(c->op, va, vb, c->imm, &folded)) {
                        c->op = OP_CONST; c->a = 0; c->b = 0; c->imm = folded;
                        s = op_desc(c->op);
                        *progress = true;
                    }
                }
            }

            /* common-subexpression elimination.  A redundant computation is
             * rewritten to a MOVE from the register that already holds the
             * value — never deleted outright, so the destination is still
             * written and no later reader can observe a stale slot.  Copy
             * propagation then routes readers past the MOVE and DCE removes it. */
            if (s.pure && s.wd && !s.rd && !is_ptr_reg(o, (int)c->dst)
                && compile_op_kind[c->op] != K_MOVE) {
                bool found = false;
                for (size_t v = 0; v < nvn; v++) {
                    VNEntry* e = &vn[v];
                    if (!e->valid || e->op != c->op) continue;
                    if (s.ra && e->a != c->a) continue;
                    if (s.rb && e->b != c->b) continue;
                    if (!imm_eq(c->op, e->imm, c->imm)) continue;
                    if (e->dst == c->dst) { found = true; break; }  /* already there */
                    c->op = OP_MOVE; c->a = e->dst; c->b = 0;
                    memset(&c->imm, 0, sizeof c->imm);
                    s = op_desc(c->op);
                    *progress = true;
                    found = true;
                    break;
                }
                if (!found && nvn < o->n) {
                    VNEntry* e = &vn[nvn++];
                    e->op = c->op; e->a = s.ra ? c->a : 0; e->b = s.rb ? c->b : 0;
                    e->imm = c->imm; e->dst = c->dst; e->valid = true;
                }
            }

            /* record the write, and invalidate everything it kills */
            if (s.wd && (int)c->dst < nreg) {
                uint32_t d = c->dst;
                for (size_t v = 0; v < nvn; v++)
                    if (vn[v].valid && (vn[v].dst == d || vn[v].a == d || vn[v].b == d))
                        vn[v].valid = false;
                for (int r = 0; r < nreg; r++) if (alias[r] == d) alias[r] = (uint32_t)r;
                /* And d's OWN alias.  Without this, `MOVE d <- a` followed by any
                 * instruction that WRITES d left alias[d] still pointing at a, so
                 * a later read of d was copy-propagated back to a — the value
                 * from before the write.  Every other opcode happened to
                 * re-establish alias[d] just below; a CALL, which writes d
                 * without being a MOVE, is what made it visible. */
                alias[d] = d;
                kk[d] = 0;
                if (compile_op_kind[c->op] == K_CONST && (int)d < o->arr_base) {
                    kk[d] = 1; kv[d] = c->imm;
                } else if (compile_op_kind[c->op] == K_MOVE && (int)c->a < nreg) {
                    alias[d] = c->a;                     /* d is now a copy of a */
                    if (kk[c->a]) { kk[d] = 1; kv[d] = kv[c->a]; }
                }
            }
        }
    }
    free(kk); free(kv); free(alias); free(vn);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Pass 2: dead-code elimination                                      *
 * ------------------------------------------------------------------ */
static bool pass_dce(Opt* o, bool* progress) {
    uint64_t* live = malloc(o->nw * sizeof(uint64_t));
    if (!live) return false;
    for (size_t b = 0; b < o->nblk; b++) {
        memcpy(live, o->live_out + b * o->nw, o->nw * sizeof(uint64_t));
        for (size_t i = o->blk_end[b]; i-- > o->blk_start[b]; ) {
            Instr* c = &o->code[i];
            OpDesc s = op_desc(c->op);
            if (compile_op_kind[c->op] == K_NOP) continue;
            /* A dead tile op IS removable — the tile buffer itself is owned by
             * the frame, not by the instruction, so dropping the write leaks
             * nothing.  Array registers stay untouchable: they own an Expr. */
            if (s.pure && s.wd && (!is_ptr_reg(o, (int)c->dst) || is_tile_reg(o, (int)c->dst))
                && !bs_get(live, (int)c->dst)) {
                c->op = OP_NOP;                      /* result never read */
                *progress = true;
                continue;
            }
            if (s.wd && !s.rd && (int)c->dst < o->nreg) bs_clr(live, (int)c->dst);
            if (s.rd && (int)c->dst < o->nreg) bs_set(live, (int)c->dst);
            if (s.ra && (int)c->a   < o->nreg) bs_set(live, (int)c->a);
            if (s.rb && !s.jump && (int)c->b < o->nreg) bs_set(live, (int)c->b);
            if (compile_op_kind[c->op] == K_CALL)
                for (unsigned k = 0; k < c->flags; k++)
                    if ((int)c->a + (int)k < o->nreg) bs_set(live, (int)c->a + (int)k);
        }
    }
    free(live);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Pass 3: loop-invariant code motion                                 *
 * ------------------------------------------------------------------ */
/* The emitter lowers every loop as a single backward JMP closing a contiguous
 * range, so a natural loop is exactly [target, jmp] for each backward JMP.  An
 * instruction in that range is invariant when it is pure, its operands are not
 * written anywhere inside the range, and its own destination is written only by
 * it.  It is hoisted to the preheader (immediately before the range), which is
 * safe for a zero-trip loop because everything hoisted is side-effect free — but
 * only if the destination is not already live there, or the hoist would clobber
 * a value the loop expects to read. */
static bool pass_licm(Opt* o, size_t* hoist_to, bool* progress) {
    unsigned char* written = malloc((size_t)o->nreg);
    unsigned char* ndef    = malloc((size_t)o->nreg);
    if (!written || !ndef) { free(written); free(ndef); return false; }

    for (size_t j = 0; j < o->n; j++) {
        int jk = compile_op_kind[o->code[j].op];
        if (jk != K_JMP && jk != K_LOOP) continue;
        size_t t = o->code[j].b;
        if (t > j) continue;                                  /* forward jump */

        memset(written, 0, (size_t)o->nreg);
        memset(ndef, 0, (size_t)o->nreg);
        bool has_arr = false;
        for (size_t i = t; i <= j; i++) {
            const Instr* c = &o->code[i];
            OpDesc s = op_desc(c->op);
            if (compile_op_kind[c->op] == K_ARR) has_arr = true;
            if (s.wd && (int)c->dst < o->nreg) {
                written[c->dst] = 1;
                if (ndef[c->dst] < 2) ndef[c->dst]++;
            }
        }
        if (has_arr) continue;               /* array ownership: leave alone */

        const uint64_t* live_at_head = o->live_in + (size_t)o->blk_of[t] * o->nw;

        for (size_t i = t; i <= j; i++) {
            Instr* c = &o->code[i];
            OpDesc s = op_desc(c->op);
            if (!s.pure || !s.wd || s.rd) continue;
            if (compile_op_kind[c->op] == K_NOP) continue;
            /* Tiles may be hoisted (a splat of a loop-invariant scalar is
             * exactly what we want out of the strip loop); array handles may not. */
            if (is_ptr_reg(o, (int)c->dst) && !is_tile_reg(o, (int)c->dst)) continue;
            if (hoist_to[i] != (size_t)-1) continue;             /* already moving */
            if (ndef[c->dst] != 1) continue;      /* destination re-assigned */
            if (bs_get(live_at_head, (int)c->dst)) continue;   /* would clobber */
            if (s.ra && (int)c->a < o->nreg && written[c->a]) continue;
            if (s.rb && !s.jump && (int)c->b < o->nreg && written[c->b]) continue;

            hoist_to[i] = t;
            written[c->dst] = 0;      /* now invariant: lets dependents hoist too */
            *progress = true;
        }
    }
    free(written); free(ndef);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Pass 4: rebuild — apply hoists, drop NOPs, remap jump targets      *
 * ------------------------------------------------------------------ */
static bool pass_rebuild(Opt* o, const size_t* hoist_to) {
    Instr* out = malloc((o->n ? o->n : 1) * sizeof(Instr));
    size_t* map = malloc((o->n + 1) * sizeof(size_t));
    if (!out || !map) { free(out); free(map); return false; }

    size_t m = 0;
    for (size_t i = 0; i < o->n; i++) {
        /* instructions hoisted to this position come first */
        if (hoist_to) {
            for (size_t h = i; h < o->n; h++)
                if (hoist_to[h] == i && compile_op_kind[o->code[h].op] != K_NOP)
                    out[m++] = o->code[h];
        }
        map[i] = m;
        if (compile_op_kind[o->code[i].op] == K_NOP) continue;
        if (hoist_to && hoist_to[i] != (size_t)-1) continue;   /* moved above */
        out[m++] = o->code[i];
    }
    map[o->n] = m;

    /* A jump target that landed on a run of removed instructions maps to the
     * next surviving instruction, which map[] already gives. */
    for (size_t i = 0; i < m; i++) {
        int k = compile_op_kind[out[i].op];
        if (k == K_JMP || k == K_JZ || k == K_LOOP) {
            size_t t = out[i].b;
            out[i].b = (uint32_t)map[t <= o->n ? t : o->n];
        }
    }

    memcpy(o->code, out, m * sizeof(Instr));
    o->n = m;
    free(out); free(map);
    return true;
}

static void opt_free(Opt* o) {
    free(o->blk_start); free(o->blk_end); free(o->blk_of);
    free(o->live_in); free(o->live_out);
    o->blk_start = NULL; o->blk_end = NULL; o->blk_of = NULL;
    o->live_in = NULL; o->live_out = NULL;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */
bool compile_optimize(Instr* code, size_t* np, int nreg, int arr_base, int tile_base) {
    if (!code || !np || *np == 0 || nreg <= 0) return true;

    Opt o;
    memset(&o, 0, sizeof o);
    o.code = code; o.n = *np; o.nreg = nreg;
    o.arr_base = arr_base; o.tile_base = tile_base;

    size_t* hoist = malloc(o.n * sizeof(size_t));
    if (!hoist) return false;

    bool ok = true;
    for (int round = 0; round < 3 && ok; round++) {
        bool progress = false;
        if (!build_cfg(&o) || !liveness(&o)) { ok = false; opt_free(&o); break; }
        ok = pass_vn(&o, &progress) && pass_dce(&o, &progress);
        if (ok && round == 0) {
            for (size_t i = 0; i < o.n; i++) hoist[i] = (size_t)-1;
            ok = pass_licm(&o, hoist, &progress);
            if (ok) ok = pass_rebuild(&o, hoist);
        }
        opt_free(&o);
        if (!progress) break;
    }

    if (ok) {
        if (!build_cfg(&o)) { opt_free(&o); return false; }
        ok = pass_rebuild(&o, NULL);
        opt_free(&o);
    }
    free(hoist);
    if (ok) *np = o.n;
    return ok;
}
