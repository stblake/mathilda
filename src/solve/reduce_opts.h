/*
 * reduce_opts.h
 *
 * Options carried through the `Reduce` engines (REDUCE_PLAN.md, Phase 8 --
 * options).  Reduce reaches its solutions by building internal `Solve[...]`
 * calls and evaluating them, so the "hard" options (Cubics, Quartics,
 * GeneratedParameters) are honored by forwarding them onto those calls
 * (reduce_opts_build_solve), reusing Solve's own solvepoly / solveinv
 * machinery.  Modulus is handled by a top-level pre-pass in builtin_reduce and
 * is deliberately NOT forwarded by reduce_opts_build_solve.  Backsubstitution
 * and WorkingPrecision are consumed directly by the engines that can honor
 * them (reduce_sys / the real sign-decision sites).
 *
 * Keep the defaults (reduce_opts_default) in sync with the Options[Reduce]
 * registry block in options_builtin.c.
 */
#ifndef REDUCE_OPTS_H
#define REDUCE_OPTS_H

#include <stdbool.h>
#include "expr.h"
#include "solvepoly.h"   /* SolvePolyOpts */

typedef struct {
    SolvePolyOpts poly;        /* Cubics / Quartics -> radical vs Root[]      */
    Expr* modulus;             /* borrowed; Modulus -> p (p>=2), else NULL    */
    const char* param_head;    /* GeneratedParameters head; default "C"       */
    bool  backsub;             /* Backsubstitution; default false             */
    Expr* working_precision;   /* borrowed; NULL / Infinity = exact-first     */
    Expr* method;              /* borrowed; NULL / Automatic                  */
} ReduceOpts;

/* Fill *o with the option defaults. */
void reduce_opts_default(ReduceOpts* o);

/* Build Solve[base_args..., <forwarded option rules>] and return it (owned).
 * Consumes the `nbase` owned Expr* in base_args, exactly as expr_new_function
 * does.  Forwards Cubics / Quartics / GeneratedParameters only when they
 * deviate from Solve's own defaults (so a default ReduceOpts reproduces the
 * original 2-/3-arg Solve call byte-for-byte).  Modulus is never forwarded
 * here. */
Expr* reduce_opts_build_solve(Expr** base_args, int nbase, const ReduceOpts* o);

/* If the working precision is a finite number, write its value (decimal
 * digits) to *digits and return true.  Infinity / MachinePrecision / NULL /
 * non-numeric -> false (exact-first: callers keep their default machine
 * tolerance). */
bool reduce_opts_wp_digits(const ReduceOpts* o, double* digits);

#endif /* REDUCE_OPTS_H */
