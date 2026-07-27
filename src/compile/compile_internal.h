/* Mathilda — Compile[] engine: internal bytecode representation.
 *
 * Shared between the emitter/VM (compile.c) and the optimiser (optimize.c).
 * Not part of the public API — clients use compile.h.
 *
 * The instruction set is three-address and monomorphic: the opcode carries the
 * operand type, so the VM does no tag dispatch.  Every opcode belongs to exactly
 * one KIND (see OPLIST), which is what tells the optimiser which fields are
 * register reads, which is a register write, and whether the instruction is pure.
 */
#ifndef MATHILDA_COMPILE_INTERNAL_H
#define MATHILDA_COMPILE_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <complex.h>
#include "../expr.h"

/* A register (or an instruction's immediate).  `p` carries a machine-kernel
 * function pointer in an immediate; `arr` carries the OWNED EXPR_NDARRAY handle
 * of an array register (M3a).  The opcode says which member is live. */
typedef union { long long i; double r; double _Complex z; const void* p; Expr* arr; } Slot;

/* scalar kernel signatures exposed by the shared machine-kernel layer */
typedef bool (*kfn_r)(double, double*);
typedef bool (*kfn_c)(double, double, double*, double*);
typedef bool (*kfn_c2)(double, double, double, double, double*, double*);
typedef bool (*kfn_cn)(const double*, const double*, size_t, double*, double*);

typedef struct { uint16_t op, flags; uint32_t dst, a, b; Slot imm; } Instr;

/* ------------------------------------------------------------------ *
 *  Opcode list                                                        *
 * ------------------------------------------------------------------ *
 * ONE list drives three things: the opcode enum, the VM's computed-goto jump
 * table, and the optimiser's instruction-property table.  Adding an opcode
 * anywhere else is a bug that shows up as a VM crash.
 *
 * KIND tells the optimiser the instruction's shape:
 *   K_CONST   dst = imm                       (pure, no register reads)
 *   K_MOVE    dst = a                         (pure)
 *   K_UN      dst = f(a)                      (pure)
 *   K_BIN     dst = f(a, b)                   (pure)
 *   K_POWI    dst = f(a, imm.i)               (pure)
 *   K_KERN1   dst = f(a) via imm.p            (pure)
 *   K_KERN2   dst = f(a, b) via imm.p         (pure)
 *   K_INC     dst = dst + imm.i               (pure, but READS dst)
 *   K_JMP     pc = b                          (b is a TARGET, not a register)
 *   K_JZ      if !a then pc = b               (b is a TARGET, not a register)
 *   K_LOOP    if ++dst < a then pc = b        (b is a TARGET; increment, test and
 *                                              branch in one instruction)
 *   K_APAR    run the strip loop that FOLLOWS this instruction in parallel over
 *             [0, a) and then pc = b; or fall through to run it serially.
 *             dst is the loop index and a the element count, exactly as for the
 *             K_LOOP that closes the range. `b` is a TARGET (one past the loop).
 *             Impure and never moved: it owns the decision to fan out.
 *   K_RET     return dst                      (READS dst)
 *   K_ARR     array op: impure, may allocate/free — never reordered or removed
 *   K_ASTORE  writes array memory or checks a shape: impure and never moved, but
 *             (unlike K_ARR) carries no ownership in its flags, so it does not
 *             stop the optimiser working on the rest of a fused loop
 *   K_VACC    dst += f(a)  (a strip-mined reduction step: READS and WRITES dst)
 *   K_CALL    dst = callee(a .. a+flags-1)   — imm.p is the callee program and
 *             `flags` the argument count, so the operands are a RANGE rather
 *             than the usual fixed fields
 *   K_NARY    dst = kernel(a .. a+flags-1)   — same operand range as K_CALL but
 *             PURE: a machine kernel has no side effects, so it may be CSE'd
 *             and removed when dead
 *   K_NOP     removed by a previous pass; deleted at compaction
 */
enum {
    K_CONST, K_MOVE, K_UN, K_BIN, K_POWI, K_KERN1, K_KERN2,
    K_INC, K_JMP, K_JZ, K_LOOP, K_RET, K_ARR, K_ASTORE, K_VACC, K_CALL,
    K_NARY, K_APAR, K_NOP
};

/* ------------------------------------------------------------------ *
 *  Strip mining                                                       *
 * ------------------------------------------------------------------ *
 * Running the scalar VM once per array element does not pay: at a few ns per
 * bytecode instruction, per-element interpretation costs more than the
 * temporary buffers a fused loop saves.  So a fused elementwise chain is
 * strip-mined instead — each opcode processes a TILE of VBLOCK elements in a
 * tight C loop.  Dispatch is then amortised VBLOCK-fold, the loop bodies are
 * shaped so the C compiler can vectorise them, and every temporary is a tile
 * that stays in L1 instead of a full-length heap buffer that page-faults.
 *
 * VBLOCK is chosen so that a handful of live tiles fit comfortably in L1:
 * 64 elements is 512 B as doubles, 1 KiB as complex.
 *
 * Tile registers live in their own bank (like the array bank) and hold a
 * pointer to their tile storage, fixed for the life of the program.  Arithmetic
 * ops always process the FULL VBLOCK — the tail of a partial tile is padded
 * with 1.0 by the load, so no operation ever sees uninitialised memory and no
 * op needs to know the live length.  Only the load, the store and the reduction
 * consult the live length, and they are pinned inside the loop anyway because
 * they read the loop index. */
#define VBLOCK 64

#define OPLIST \
    X(NOP,   K_NOP)   X(JMP,   K_JMP)  X(JZ,    K_JZ)   X(INC_I, K_INC)   \
    X(LOOP,  K_LOOP)  X(APAR,  K_APAR)                                     \
    X(CONST, K_CONST) X(MOVE,  K_MOVE)                                     \
    X(I2R,   K_UN)    X(I2C,   K_UN)   X(R2C,   K_UN)                      \
    X(ADD_I, K_BIN)   X(ADD_R, K_BIN)  X(ADD_C, K_BIN)                     \
    X(SUB_I, K_BIN)   X(SUB_R, K_BIN)  X(SUB_C, K_BIN)                     \
    X(MUL_I, K_BIN)   X(MUL_R, K_BIN)  X(MUL_C, K_BIN)                     \
    X(DIV_R, K_BIN)   X(DIV_C, K_BIN)  X(MOD_I, K_BIN)  X(QUOT_I, K_BIN)   \
    X(NEG_I, K_UN)    X(NEG_R, K_UN)   X(NEG_C, K_UN)                      \
    X(INV_R, K_UN)    X(INV_C, K_UN)                                       \
    X(POWI_I, K_POWI) X(POWI_R, K_POWI) X(POWI_C, K_POWI)                  \
    X(POW_R, K_BIN)   X(POW_C, K_BIN)                                      \
    X(SQRT_R, K_UN)   X(SQRT_C, K_UN)  X(EXP_R, K_UN)   X(EXP_C, K_UN)     \
    X(LOG_R, K_UN)    X(LOG_C, K_UN)                                       \
    X(SIN_R, K_UN)    X(SIN_C, K_UN)   X(COS_R, K_UN)   X(COS_C, K_UN)     \
    X(TAN_R, K_UN)    X(TAN_C, K_UN)                                       \
    X(SINH_R, K_UN)   X(SINH_C, K_UN)  X(COSH_R, K_UN)  X(COSH_C, K_UN)    \
    X(TANH_R, K_UN)   X(TANH_C, K_UN)                                      \
    X(ASIN_R, K_UN)   X(ASIN_C, K_UN)  X(ACOS_R, K_UN)  X(ACOS_C, K_UN)    \
    X(ATAN_R, K_UN)   X(ATAN_C, K_UN)                                      \
    X(ABS_I, K_UN)    X(ABS_R, K_UN)   X(ABS_C, K_UN)                      \
    X(SIGN_I, K_UN)   X(SIGN_R, K_UN)                                      \
    X(FLOOR_R, K_UN)  X(CEIL_R, K_UN)  X(ROUND_R, K_UN)                    \
    X(RE_C, K_UN)     X(IM_C, K_UN)    X(ARG_C, K_UN)   X(CONJ_C, K_UN)    \
    X(ATAN2_R, K_BIN) X(MAX_I, K_BIN)  X(MAX_R, K_BIN)                     \
    X(MIN_I, K_BIN)   X(MIN_R, K_BIN)                                      \
    X(ERF_R, K_UN)    X(ERFC_R, K_UN)                                      \
    X(KERN_RR, K_KERN1) X(KERN_R2R, K_KERN1) X(KERN_RC, K_KERN1)           \
    X(KERN_CC, K_KERN1) X(KERN_CR, K_KERN1)                                \
    X(KERN2_RR, K_KERN2) X(KERN2_RC, K_KERN2) X(KERN2_CC, K_KERN2)         \
    X(LT_I, K_BIN)    X(LT_R, K_BIN)   X(LE_I, K_BIN)   X(LE_R, K_BIN)     \
    X(GT_I, K_BIN)    X(GT_R, K_BIN)   X(GE_I, K_BIN)   X(GE_R, K_BIN)     \
    X(EQ_I, K_BIN)    X(EQ_R, K_BIN)   X(EQ_C, K_BIN)                      \
    X(NE_I, K_BIN)    X(NE_R, K_BIN)   X(NE_C, K_BIN)                      \
    X(AND, K_BIN)     X(OR, K_BIN)     X(XOR, K_BIN)    X(NOT, K_UN)       \
    X(ARR_FREE, K_ARR) X(V_EW, K_ARR)  X(V_POW, K_ARR)                     \
    X(V_KERN, K_ARR)  X(V_KERN2, K_ARR) X(V_TOTAL, K_ARR) X(V_LEN, K_ARR)  \
    X(A_SIZE, K_UN)   X(A_NEWLIKE, K_ARR) X(A_SHAPECHK, K_ASTORE)          \
    X(A_LOAD_R, K_BIN)   X(A_LOAD_C, K_BIN)                                \
    X(A_STORE_R, K_ASTORE) X(A_STORE_C, K_ASTORE)                          \
    /* strip-mined (tile) forms — one opcode, VBLOCK elements */           \
    X(CALL, K_CALL)   X(KERNN, K_NARY)                                     \
    X(VSETLEN, K_ASTORE)                                                   \
    X(VLOAD_R, K_BIN)   X(VLOAD_C, K_BIN)                                  \
    X(VSTORE_R, K_ASTORE) X(VSTORE_C, K_ASTORE)                            \
    X(VSUM_R, K_VACC)   X(VSUM_C, K_VACC)                                  \
    X(VSPLAT_R, K_UN)   X(VSPLAT_C, K_UN)   X(VR2C, K_UN)                  \
    X(VADD_R, K_BIN)    X(VSUB_R, K_BIN)    X(VMUL_R, K_BIN)               \
    X(VDIV_R, K_BIN)    X(VNEG_R, K_UN)     X(VINV_R, K_UN)                \
    X(VADD_C, K_BIN)    X(VSUB_C, K_BIN)    X(VMUL_C, K_BIN)               \
    X(VDIV_C, K_BIN)    X(VNEG_C, K_UN)     X(VINV_C, K_UN)                \
    X(VPOWI_R, K_POWI)  X(VPOWI_C, K_POWI)                                 \
    X(VPOW_R, K_BIN)    X(VPOW_C, K_BIN)    X(VATAN2_R, K_BIN)             \
    X(VSQRT_R, K_UN)    X(VSQRT_C, K_UN)                                   \
    X(VEXP_R, K_UN)     X(VEXP_C, K_UN)     X(VLOG_R, K_UN)  X(VLOG_C, K_UN) \
    X(VSIN_R, K_UN)     X(VSIN_C, K_UN)     X(VCOS_R, K_UN)  X(VCOS_C, K_UN) \
    X(VTAN_R, K_UN)     X(VTAN_C, K_UN)                                    \
    X(VSINH_R, K_UN)    X(VSINH_C, K_UN)    X(VCOSH_R, K_UN) X(VCOSH_C, K_UN) \
    X(VTANH_R, K_UN)    X(VTANH_C, K_UN)                                   \
    X(VASIN_R, K_UN)    X(VASIN_C, K_UN)    X(VACOS_R, K_UN) X(VACOS_C, K_UN) \
    X(VATAN_R, K_UN)    X(VATAN_C, K_UN)                                   \
    X(VABS_R, K_UN)     X(VABS_C, K_UN)     X(VSIGN_R, K_UN)               \
    X(VRE_C, K_UN)      X(VIM_C, K_UN)      X(VARG_C, K_UN) X(VCONJ_C, K_UN) \
    X(VERF_R, K_UN)     X(VERFC_R, K_UN)                                   \
    X(VKERN_RR, K_KERN1)  X(VKERN_R2R, K_KERN1) X(VKERN_RC, K_KERN1)       \
    X(VKERN_CC, K_KERN1)  X(VKERN_CR, K_KERN1)                             \
    X(VKERN2_RR, K_KERN2) X(VKERN2_RC, K_KERN2) X(VKERN2_CC, K_KERN2)      \
    X(RET, K_RET)

/* The opcode enum, generated from OPLIST so the two cannot drift apart. */
#define X(name, kind) OP_##name,
enum { OPLIST OP__COUNT };
#undef X

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

/* Array and tile temporaries are allocated into virtual ranges and relocated to
 * contiguous banks above the scalar registers at finalize (see patch_reg).  A
 * slot is then always-scalar, always-array or always-tile, so teardown can never
 * mistake a double for a pointer. */
#define ARR_VREG  0x40000000
#define TILE_VREG 0x20000000

/* Kind of each opcode, indexed by opcode.  Defined in optimize.c. */
extern const unsigned char compile_op_kind[OP__COUNT];

/* Optimise `code` (length *n, updated in place) for a program with `nreg`
 * registers whose array bank starts at `arr_base`, and whose result lives in
 * *result_reg (updated if the result instruction moves).  Registers in
 * [arr_base, tile_base) are array handles and are never touched by the scalar
 * passes; registers >= tile_base are strip-mining tiles, which may be hoisted or
 * removed but never copied (a MOVE would alias two tiles onto one buffer).
 *
 * Purely a bytecode-to-bytecode transform: same observable results, fewer
 * instructions.  Returns false only on allocation failure, in which case `code`
 * is left exactly as it was (the caller can still run it). */
bool compile_optimize(Instr* code, size_t* n, int nreg, int arr_base, int tile_base);

#endif /* MATHILDA_COMPILE_INTERNAL_H */
