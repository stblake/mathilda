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
 *   K_RET     return dst                      (READS dst)
 *   K_ARR     array op: impure, may allocate/free — never reordered or removed
 *   K_NOP     removed by a previous pass; deleted at compaction
 */
enum {
    K_CONST, K_MOVE, K_UN, K_BIN, K_POWI, K_KERN1, K_KERN2,
    K_INC, K_JMP, K_JZ, K_RET, K_ARR, K_NOP
};

#define OPLIST \
    X(NOP,   K_NOP)   X(JMP,   K_JMP)  X(JZ,    K_JZ)   X(INC_I, K_INC)   \
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

/* Array temporaries are allocated into a virtual range and relocated to a
 * contiguous bank above the scalar registers at finalize (see patch_reg). */
#define ARR_VREG 0x40000000

/* Kind of each opcode, indexed by opcode.  Defined in optimize.c. */
extern const unsigned char compile_op_kind[OP__COUNT];

/* Optimise `code` (length *n, updated in place) for a program with `nreg`
 * registers whose array bank starts at `arr_base`, and whose result lives in
 * *result_reg (updated if the result instruction moves).  Registers >= arr_base
 * are array handles and are never touched by the scalar passes.
 *
 * Purely a bytecode-to-bytecode transform: same observable results, fewer
 * instructions.  Returns false only on allocation failure, in which case `code`
 * is left exactly as it was (the caller can still run it). */
bool compile_optimize(Instr* code, size_t* n, int nreg, int arr_base);

#endif /* MATHILDA_COMPILE_INTERNAL_H */
