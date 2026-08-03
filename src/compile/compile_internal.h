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
#include <limits.h>
#include <complex.h>
#include "../expr.h"
#include "../checked_int.h"   /* ci_add / ci_sub / ci_mul / ci_neg / ci_abs / ci_powi */
#include "compile.h"      /* CompileType, CompiledProgram — completed below */

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
 * ONE list drives four things: the opcode enum, the VM's computed-goto jump
 * table, the optimiser's instruction-property table, and the disassembler's
 * opcode-name table.  Adding an opcode anywhere else is a bug that shows up as
 * a VM crash.
 *
 * KIND tells the optimiser the instruction's shape:
 *   K_CONST   dst = imm                       (pure, no register reads)
 *   K_MOVE    dst = a                         (pure)
 *   K_UN      dst = f(a)                      (pure)
 *   K_BIN     dst = f(a, b)                   (pure)
 *   K_BINK    dst = f(a, imm)                 (pure) -- a binary op with its
 *             constant operand folded INTO the instruction.  Created by the
 *             optimiser, never by the emitter: materialising a constant into
 *             a register costs a whole instruction, that instruction
 *             re-executes on every call, and in a polynomial body it is a
 *             THIRD of the stream.  Exactly value-preserving (same C operator
 *             on the same two values), so the opt/no-opt bitwise gate covers
 *             it.
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
 *   K_AIDX    dst = dst*dim(b, imm.i) + resolve(a)   — one axis of a subscript,
 *             READS and WRITES dst.  Folding the multiply, the 1-based/negative
 *             resolution and the per-axis range check into a single instruction
 *             is what keeps `u[[i, j]]` at four instructions instead of ten;
 *             the range check is per AXIS, not on the flat index, because
 *             u[[1, n + 5]] is in range flat and reads the wrong row.  Never
 *             removed: dropping it would drop the check.
 *   K_CALL    dst = callee(a .. a+flags-1)   — imm.p is the callee program and
 *             `flags` the argument count, so the operands are a RANGE rather
 *             than the usual fixed fields
 *   K_NARY    dst = kernel(a .. a+flags-1)   — same operand range as K_CALL but
 *             PURE: a machine kernel has no side effects, so it may be CSE'd
 *             and removed when dead
 *   K_ALOAD   dst = a[b], a READ of array memory.  Shaped exactly like K_BIN but
 *             deliberately NOT pure: it was pure while the only writer was
 *             fusion, which stores solely into a fresh result buffer, but once
 *             user code can write a buffer it also reads (`u[[k]] = v` next to
 *             `u[[k]]`) a pure load would be CSE'd across the store, and LICM
 *             would hoist a loop-invariant load out of the loop that mutates it.
 *   K_ARANGE  an array op whose operands are a RANGE: `a` is the first of
 *             `flags` registers, exactly as for K_CALL.  Used where a subscript
 *             list or a dimension list has to reach the VM as values.
 *   K_NOP     removed by a previous pass; deleted at compaction
 */
enum {
    K_CONST, K_MOVE, K_UN, K_BIN, K_BINK, K_POWI, K_KERN1, K_KERN2,
    K_INC, K_JMP, K_JZ, K_LOOP, K_RET, K_ARR, K_ASTORE, K_VACC, K_CALL,
    K_NARY, K_APAR, K_AIDX, K_ARANGE, K_ALOAD, K_NOP
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

/* Application cap for an iteration with no user-supplied bound (FixedPoint,
 * NestWhile).  MUST equal ITER_SAFETY_CAP in src/funcprog.c:2062: on reaching
 * it the interpreter frees its history and leaves the whole expression
 * unevaluated, so a compiled loop that ran on would answer where the
 * interpreter declines.  Reaching it emits OP_FAIL, which drops the call to the
 * interpreter — which then gives up in its turn. */
#define VM_ITER_SAFETY_CAP 1000000LL

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
    /* immediate forms: the optimiser folds a constant operand into the      \
     * instruction and DCE then removes the CONST that materialised it.      \
     * `_RK`/`_IK` = register OP constant; `_KR`/`_KI` = constant OP register \
     * (only needed where the op is not commutative).  A comparison with the  \
     * constant on the LEFT is rewritten by SWAPPING the predicate, which is  \
     * NaN-safe: `k < x` and `x > k` are both false for a NaN x. */          \
    X(ADD_RK, K_BINK) X(SUB_RK, K_BINK) X(SUB_KR, K_BINK)                   \
    X(MUL_RK, K_BINK) X(DIV_RK, K_BINK) X(DIV_KR, K_BINK)                   \
    X(ADD_IK, K_BINK) X(SUB_IK, K_BINK) X(SUB_KI, K_BINK) X(MUL_IK, K_BINK) \
    X(LT_RK, K_BINK)  X(LE_RK, K_BINK)  X(GT_RK, K_BINK)  X(GE_RK, K_BINK)  \
    X(LT_IK, K_BINK)  X(LE_IK, K_BINK)  X(GT_IK, K_BINK)  X(GE_IK, K_BINK)  \
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
    /* Integer-CLOSED forms of heads whose real kernels return a double.     \
     * `Binomial[7,3]` is the Integer 35 and `2^10` the Integer 1024, so     \
     * lowering these through the real kernel would answer 35. and 1024. —   \
     * a different HEAD from the interpreter's, which the engine forbids     \
     * however equal the values print.  Each is overflow-checked like the    \
     * arithmetic opcodes, and each ABANDONS the call outside the domain     \
     * where the interpreter returns an integer at all: a negative exponent  \
     * is a Rational, Factorial of a negative is ComplexInfinity, and the    \
     * interpreter is the only thing that can say so. */                     \
    X(POW_II, K_BIN)  X(FACT_I, K_UN)  X(GAMMA_I, K_UN)                    \
    X(BINOM_I, K_BIN) X(POCH_I, K_BIN) X(FIB_I, K_UN)  X(LUCAS_I, K_UN)    \
    /* Arg of an integer is the Integer 0 at or above zero and the SYMBOL Pi   \
     * below it — a head that depends on the sign of a value only known per     \
     * call, so the negative side abandons the call rather than answering       \
     * 3.14159 where the interpreter says Pi. */                                \
    X(ARG_I, K_UN)                                                          \
    /* Integer-only heads: no real counterpart at all, so these have no double \
     * kernel to fall back on and exist only over CT_INT.  ILEN/IEXP take the  \
     * base as their second operand (the emitter materialises the default 10). \
     * PowerMod is K_NARY because it is ternary and an Instr addresses only two \
     * source registers — same shape as KERNN, a run starting at `a`.          \
     *                                                                          \
     * The Bit* family is deliberately absent: BitAnd/BitOr/BitXor/BitNot/       \
     * BitShiftLeft/BitShiftRight are not implemented in the INTERPRETER (they   \
     * come back unevaluated), and a compiled path that answered them would      \
     * answer where the interpreter declines. */                                 \
    X(GCD_I, K_BIN)   X(LCM_I, K_BIN)  X(ILEN_I, K_BIN) X(IEXP_I, K_BIN)   \
    X(POWMOD_I, K_NARY)                                                     \
    X(FLOOR_R, K_UN)  X(CEIL_R, K_UN)  X(ROUND_R, K_UN)  X(TRUNC_R, K_UN)   \
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
    /* SameQ, which is NOT Equal: expr_eq calls two NaNs the same            \
     * (src/expr.c:622), and that is exactly what stops FixedPoint spinning  \
     * on a NaN orbit.  One opcode rather than a compare/isnan/or chain       \
     * because it sits in the innermost loop and is easy to get subtly wrong. \
     * The complex form applies the rule COMPONENTWISE, mirroring expr_eq's   \
     * recursion into Complex[re, im] — not creal==creal && cimag==cimag,     \
     * which differs on a NaN component. */                                  \
    X(SAMEQ_R, K_BIN) X(SAMEQ_C, K_BIN)                                     \
    /* Abort the call: the result is unusable, so compiled_eval returns false \
     * and the caller interprets.  Used where the interpreter itself gives up \
     * (an unbounded iteration reaching ITER_SAFETY_CAP), so the observable    \
     * behaviour is identical — the compiled path merely burns the iterations  \
     * first.  Impure and never removed. */                                   \
    X(FAIL, K_ASTORE)                                                       \
    X(AND, K_BIN)     X(OR, K_BIN)     X(XOR, K_BIN)    X(NOT, K_UN)       \
    X(ARR_FREE, K_ARR) X(V_EW, K_ARR)  X(V_POW, K_ARR)                     \
    X(V_KERN, K_ARR)  X(V_KERN2, K_ARR) X(V_TOTAL, K_ARR) X(V_LEN, K_ARR)  \
    /* An array -> SCALAR reduction delegated to the interpreter's own entry   \
     * point (ndred_mean, ndred_variance, ...), the reduction counterpart of   \
     * A_NDFN below.  Total has its own opcode because an int64 sum must stay  \
     * exact past 2^53; everything here is real-valued (Mean of an exact       \
     * integer vector is a Rational, which no machine slot holds) and so needs \
     * no such split.  imm.p is the NdRedSpec. */                              \
    X(V_NDRED, K_ARR)                                                       \
    /* Delegated TWO-array heads (COMPILE_MISSING.md §3).  A_NDFN2 reads two    \
     * array registers and produces an array (Dot matrix shapes, LinearSolve,   \
     * Cross, LeastSquares, ListConvolve/Correlate, Join); V_NDFN2 the SCALAR    \
     * form (Dot's vector.vector inner product).  Both rebuild the whole call    \
     * and delegate, like A_NDFN.  imm.p is the NdFn2Spec. */                    \
    X(A_NDFN2, K_ARR) X(V_NDFN2, K_ARR)                                     \
    X(A_SIZE, K_UN)   X(A_NEWLIKE, K_ARR) X(A_SHAPECHK, K_ASTORE)          \
    X(A_LOAD_R, K_ALOAD) X(A_LOAD_C, K_ALOAD) X(A_LOAD_I, K_ALOAD)         \
    X(A_LOAD_B, K_ALOAD)                                                    \
    X(A_STORE_R, K_ASTORE) X(A_STORE_C, K_ASTORE) X(A_STORE_I, K_ASTORE)   \
    X(A_STORE_B, K_ASTORE)                                                  \
    /* Growable append + shrink, for an iteration whose LENGTH is only known    \
     * once it has run (FixedPointList, NestWhileList, and later Select).       \
     * Deliberately separate opcodes rather than a bounds check inside          \
     * A_STORE: that check would then sit in the innermost loop of every fused  \
     * map, forever, for a feature that concerns a handful of lowerings.        \
     * A_PUSH doubles the buffer when the index reaches the end; A_TRUNC sets   \
     * the final length, so the capacity a run happened to reach is never       \
     * observable.  Element kind travels in `flags` (AF_R). */                  \
    X(A_PUSH, K_ASTORE)  X(A_TRUNC, K_ASTORE)                              \
    /* indexed Part (M3c).  A_AXIS resolves ONE subscript against one axis   \
     * and folds it into the running flat index, so the general A_LOAD /      \
     * A_STORE above serve both fusion and user indexing.  A_PART / A_PARTSET \
     * are the escape hatch for every spec that is not a full run of scalar    \
     * subscripts (Span, All, a list of positions, partial indexing): they     \
     * delegate to the interpreter's own ndarray_part / ndarray_part_set, so   \
     * the compiled subset of Part is the interpreted one by construction. */ \
    X(A_AXIS, K_AIDX)    X(A_COPY, K_ARR)    X(A_XFER, K_ARR)             \
    X(A_NEW, K_ARANGE)   X(A_PART, K_ARANGE) X(A_PARTSET, K_ARANGE)        \
    /* A structural head delegated to the interpreter's OWN NDArray entry      \
     * point (ndstruct_reverse, ndred_accumulate, ...), which takes the whole   \
     * call.  Same trick as A_PART, and the same payoff: the compiled subset of \
     * these heads IS the interpreted one by construction.  Their faithful-      \
     * degrade path returns a nested List instead of a packed array, which is    \
     * exactly the signal to hand the call back to the interpreter. */          \
    X(A_NDFN, K_ARANGE)                                                     \
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
/* Operand kinds for the array opcodes. Three bits each (AF_A/AF_B), so there is
 * room; AK_INT was added when int64 buffers stopped being internal to the
 * compiler.
 *
 * AK_INT is load-bearing, not tidiness. Without it a scalar operand of an array
 * op was only ever Real or Complex, so `Compile[{{u, _Integer, 1}}, u * 2]`
 * loaded the literal 2 into a real register and boxed it as a Real -- and the
 * interpreter's elementwise path then correctly widened the int64 buffer to
 * float64, answering {2., 4., 6.} where Range[3] * 2 is {2, 4, 6}. A wrong
 * compiled answer, which is the one thing this engine may not produce. */
enum { AK_ARR = 0, AK_REAL = 1, AK_COMPLEX = 2, AK_INT = 3 };

/* Scalar-integer opcodes: do not abandon the call on overflow, keep the wrapped
 * value (COMPILE_WRAP_INT / "CatchMachineIntegerOverflow" -> False).
 *
 * It shares the `flags` word with the AF_* bits above, which is safe because the
 * integer opcodes never carry any of them — an array op is never one of the ops
 * below.  A top bit rather than a low one so the two groups stay visually
 * distinct in a disassembly.
 *
 * Baking the choice into the instruction is what makes the option free: the VM
 * evaluates this bit only AFTER an overflow has been detected, so the
 * no-overflow path — every iteration of every real program — is byte for byte
 * the path it would take with the option on. */
#define IF_NOCHK      0x8000u

/* The opcodes IF_NOCHK applies to: the ARITHMETIC ones, whose only failure is
 * overflow.  Kept as one predicate so the emitter, the optimiser and the
 * disassembler cannot disagree about which instructions carry the bit.
 *
 * The integer-closed heads (POW_II, FACT_I, BINOM_I, ...) are deliberately NOT
 * here even though their VM bodies also use IOP().  They fail for TWO reasons —
 * overflow, and an argument outside the domain where the interpreter returns an
 * integer at all (`7^-3` is a Rational, `Factorial[-1]` is ComplexInfinity) —
 * and the two are entangled inside one helper.  Wrapping is a coherent answer
 * for `a*b`; there is no coherent wrapped answer for `Factorial[-1]`, so these
 * always defer to the interpreter whatever the option says. */
static inline bool op_is_checked_int(uint16_t op) {
    switch (op) {
        case OP_ADD_I:  case OP_SUB_I:  case OP_MUL_I:
        case OP_ADD_IK: case OP_SUB_IK: case OP_SUB_KI: case OP_MUL_IK:
        case OP_NEG_I:  case OP_ABS_I:  case OP_POWI_I: case OP_INC_I:
            return true;
        default:
            return false;
    }
}

/* The subscript list of a general (non-scalar) Part, carried in A_PART /
 * A_PARTSET's immediate as `imm.p`.
 *
 * A spec like `u[[i, 2 ;; 5, All]]` mixes subscripts known at compile time with
 * ones that are only known per call, so neither a plain Expr nor a plain
 * register list can describe it. `lit[k]` holds the literal spec for axis k, or
 * NULL when that axis's subscript is computed, in which case `reg[k]` names the
 * register holding its integer value and the Expr is built at call time.
 *
 * The literals are DEEP COPIES: the program can outlive the body it was
 * compiled from (a CompiledFunction is assigned to a symbol, the body Expr is
 * freed), and a borrowed pointer here would dangle. The program owns every
 * PartSpec and frees them in compiled_free. */
typedef struct {
    int    n;             /* subscripts given (may be < rank: trailing All) */
    Expr** lit;           /* [n] literal spec, or NULL if computed */
    int*   reg;           /* [n] register with the computed index, else -1 */
    int    rhs_kind;      /* A_PARTSET: AK_ARR / AK_REAL / AK_COMPLEX */
} PartSpec;

void compile_partspec_free(PartSpec* p);

/* ------------------------------------------------------------------ *
 *  The finished program                                               *
 * ------------------------------------------------------------------ *
 * Lives here rather than in compile.c because three modules need to see inside
 * it: the emitter/VM that builds and runs it, and the disassembler (disasm.c)
 * that renders it.  It stays out of the PUBLIC compile.h — clients hold it as
 * an opaque handle and go through the accessors. */

/* A strip-mined MAP loop lifted out as a self-contained program, so a worker
 * thread can run it with the ordinary `vm_run` entry over its own index range.
 *
 * Extracting rather than teaching `vm_run` to stop at an arbitrary pc is what
 * keeps the change out of the hot path entirely: a stop-pc test would cost a
 * comparison on EVERY dispatch, for a feature that concerns one loop. The
 * lifted copy ends in its own OP_RET, so the VM's inner loop is untouched. */
typedef struct {
    Instr*   code;        /* [n], internal jump targets rebased to 0 */
    size_t   n;
    uint32_t ri, rn;      /* loop index / element-count registers, in frame coords */
    /* Enough of the parent's frame shape to build a worker's own frame. Carried
     * here rather than reached through the CompiledProgram so that `vm_run`
     * needs no extra parameter — it is called on every scalar body too. */
    size_t   frame_slots;
    int      nreg, tile_base, ntiles;
} ParLoop;

struct CompiledProgram {
    Instr*      code;
    size_t      n;
    ParLoop*    ploops;       /* [nploops] parallelisable strip loops */
    int         nploops;
    PartSpec**  parts;        /* [nparts] general-Part subscript lists */
    int         nparts;
    int         nreg;
    int         arr_base;     /* array registers are [arr_base, tile_base) */
    int         tile_base;    /* tile registers are [tile_base, nreg) */
    int         result_reg;
    int         ncse;         /* repeated subtrees hoisted by cse_plan */
    CompileType result_type;
    size_t      nargs;
    CompileType* arg_types;   /* [nargs] */
    unsigned char* argdep;    /* [nargs] which args are read */
    size_t      frame_slots;  /* registers + tile storage, in Slots */
    int         ntiles;       /* tile registers; storage is per-CALL, not per-program */
    bool        all_real;     /* every arg + result is CT_REAL, no array temps */
    /* The array result was CONSTRUCTED by the body (ConstantArray, Table, ...)
     * rather than derived from an array argument, so the interpreter running the
     * same body returns a List however the arguments were spelled.  See Val.built
     * in compile.c and the result-kind decision in compiled_function.c. */
    bool        result_built;
};

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

/* The name of the head an A_NDFN / V_NDRED instruction delegates to.
 *
 * The two tables (ND_FNS, ND_REDS) live in compile.c and stay there; the
 * disassembler needs only the NAME out of an entry, so it asks rather than the
 * tables becoming shared state.  Without these an A_NDFN disassembles as
 * "A_NDFN" and a V_NDRED as "V_NDRED(V0, V0)", neither of which says WHICH head
 * was delegated -- and telling Mean from Median in a bytecode dump is the whole
 * point of having a dump.  `imm` is the instruction's Slot.p; NULL-safe. */
const char* nd_fn_head_name(const void* imm);
const char* nd_red_head_name(const void* imm);
/* The head an A_NDFN2 / V_NDFN2 delegates to (ND_FN2S in compile.c). */
const char* nd_fn2_head_name(const void* imm);

#endif /* MATHILDA_COMPILE_INTERNAL_H */
