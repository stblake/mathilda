/* Mathilda — Compile[]: arbitrary-precision special functions.
 *
 * The managed (WorkingPrecision -> n) emitter lowers straight-line arithmetic
 * and the elementary functions onto raw mpfr_* / ncpx_* opcodes.  Every other
 * numeric special function (Erf, BesselK, Hypergeometric1F1, MathieuC, ...)
 * already has a careful arbitrary-precision evaluator behind its builtin — the
 * routine N[f[...], n] reaches.  This module lowers a call to such a head as ONE
 * opcode, OP_MG_SPECIAL, that hands the argument containers to exactly that
 * routine, so compiled == uncompiled by construction and no numerics are
 * duplicated.
 *
 * It is table-driven, mirroring the machine kernel registry: one row per head
 * (name, arity range, whether a complex argument is accepted).  The table is the
 * single source of truth the coverage audit reads
 * (tools/compile_coverage.py --precision).
 *
 * See docs/design/compile-arbitrary-precision.md §19. */
#ifndef MATHILDA_COMPILE_MGD_SPECIAL_H
#define MATHILDA_COMPILE_MGD_SPECIAL_H

#include <stdbool.h>
#include <stddef.h>
#include "../expr.h"

/* Registry flags. */
#define MGS_REAL   0x1u   /* real arguments                                     */
#define MGS_CPLX   0x2u   /* complex arguments (the arb evaluator takes them)   */
#define MGS_COUT   0x4u   /* complex-valued even for real arguments (HankelH1)  */

typedef struct {
    const char*   name;      /* head, as spelled                               */
    unsigned char min_args;  /* arity range accepted                           */
    unsigned char max_args;
    unsigned char flags;     /* MGS_*                                          */
} MgdSpecialEntry;

/* The registry row for head `h` (an INTERNED name) at arity `nargs`, or NULL. */
const MgdSpecialEntry* mgd_special_lookup(const char* h, size_t nargs);

/* The whole table (for the Compile`ArbitraryPrecisionHeads introspection
 * builtin and the audit): *n receives its length. */
const MgdSpecialEntry* mgd_special_table(size_t* n);

/* A lowered call, owned by the compiled program (an MgdConst of kind
 * MGD_CONST_SPECIAL).  Argument i is either the program-owned literal lit[i]
 * (an exact parameter such as the order 2 of BesselJ[2, x], passed as the
 * interpreter would see it) or, when lit[i] is NULL, the next register of the
 * instruction's operand range; cplx[i] says that register holds an ncpx. */
typedef struct {
    const char*    head;     /* interned                                       */
    size_t         nargs;
    Expr**         lit;      /* [nargs]                                        */
    unsigned char* cplx;     /* [nargs]                                        */
    bool           out_cplx; /* the result register is an ncpx                 */
} MgdSpecialCall;

void mgd_special_call_free(MgdSpecialCall* sc);

#ifdef USE_MPFR
/* Run one lowered call: box the operand containers `ops` (pointers to mpfr_t /
 * ncpx, in operand order), evaluate the head through its own evaluator, and
 * write the value into `dst` (an mpfr_t, or an ncpx when dst_cplx).  False
 * when the head does not produce an arbitrary-precision number of the static
 * result type (a complex value where a real was compiled, a machine Real, an
 * unevaluated call, a non-finite value): the VM then fails the call and the
 * caller interprets that point, the same contract as every managed opcode. */
bool mgd_special_invoke(const MgdSpecialCall* sc, void* const* ops,
                        void* dst, bool dst_cplx);
#endif

/* Process-wide counters for tests: lowered calls run, and calls that declined
 * at run time.  Not thread-safe by design (diagnostic only). */
void mgd_special_stats(long* calls, long* declines);

/* Arbitrary-precision compiled_eval calls that declined to the interpreter
 * (compile_vm.c). */
long compiled_managed_eval_failures(void);

/* Registers Internal`CompileArbitraryPrecisionHeads and
 * Internal`CompileArbitraryPrecisionStats (called from compiled_function_init). */
void compile_mgd_special_init(void);

#endif /* MATHILDA_COMPILE_MGD_SPECIAL_H */
