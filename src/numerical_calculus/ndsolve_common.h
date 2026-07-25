/* Mathilda — NDSolve shared infrastructure.
 *
 * NDSolve is factored into two independent axes ("two-way modularity"):
 *
 *   1. Problem class.  A user problem (ODE IVP; later BVP/DAE/PDE) is compiled
 *      by the front-end in ndsolve.c into a first-order autonomous-or-not system
 *
 *          dY/dt = f(t, Y),   Y in R^d,   Y(t0) = Y0
 *
 *      captured by `NdProblem`.  Higher-order ODEs are reduced to first order;
 *      complex-valued ODEs are realified (Re/Im split, dimension doubled) so the
 *      entire real machinery below is reused.
 *
 *   2. Numerical method.  Each integrator (Euler, Midpoint, RK4, DOPRI5,
 *      BackwardEuler, BDF, Adams, ...) is a self-contained module exposing the
 *      common `NdStepper` vtable and implementing only *single-step* logic.  A
 *      single shared adaptive driver (`nd_integrate`) owns the time loop,
 *      step-size control, rejection, monitors and sample collection.
 *
 * The RHS sampler evaluates the symbolic reduced right-hand sides at numeric
 * points by Block-localizing the reduced-state symbols and the independent
 * variable, following the NIntegrate binding pattern (src/numerical_calculus/
 * nint.c).  Symbolic Jacobians and the dense linear solve reused by the
 * implicit steppers follow the FindRoot pattern (src/numerical_roots/findroot.c).
 */
#ifndef NDSOLVE_COMMON_H
#define NDSOLVE_COMMON_H

#include "../expr.h"
#include "../symtab.h"
#include "../numeric.h"
#include <stdbool.h>
#include <stddef.h>

/* ------------------------------------------------------------------ *
 *  Block-style variable binding (mirrors nint.c's NiBind)             *
 * ------------------------------------------------------------------ */
typedef struct {
    const char* name;      /* interned symbol name (borrowed)         */
    Rule*       saved_own; /* the symbol's OwnValues before binding    */
    uint32_t    saved_attrs;
    bool        valid;
} NdBind;

void nd_bind_snapshot(NdBind* b, const char* name);
void nd_bind_set(NdBind* b, Expr* value);   /* value is copied internally */
void nd_bind_restore(NdBind* b);

/* Coerce a real-valued numeric leaf (Integer/Real/BigInt/Rational/MPFR) to a
 * double.  Returns false, leaving *out untouched, for anything non-real. */
bool nd_to_double(const Expr* e, double* out);

/* ------------------------------------------------------------------ *
 *  Compiled problem: dY/dt = f(t, Y)                                  *
 * ------------------------------------------------------------------ */
typedef struct NdProblem NdProblem;

struct NdProblem {
    size_t   d;             /* reduced first-order dimension                 */

    /* symbolic evaluation context (owned unless noted) */
    NdBind   bind_t;        /* bind for the independent variable t            */
    NdBind*  bind_y;        /* array[d], one per reduced-state symbol         */
    Expr**   f;             /* array[d] reduced RHS bodies f_i(t, Y) (owned)  */
    Expr**   ysym;          /* array[d] reduced-state symbol Exprs (owned)    */
    NumericSpec spec;       /* machine double vs MPFR                          */
    Expr*    eval_monitor;  /* borrowed EvaluationMonitor body, or NULL       */

    /* Jacobian (implicit steppers only; built lazily, owned) */
    Expr***  jac;           /* jac[i][j] = D[f_i, y_j], or NULL entry -> FD    */
    bool     jac_built;

    /* domain + initial state */
    const char* tvar;       /* interned independent-variable name (borrowed)  */
    double   t0, tmin, tmax;
    double*  Y0;            /* array[d] (owned)                               */

    /* resolved tolerances, set by the driver before stepping (for implicit
     * steppers' Newton convergence test). */
    struct { double rtol, atol; } tol;

    /* FSAL (first-same-as-last) state, owned by the adaptive driver during an
     * integration.  For an FSAL stepper (DOPRI5) the last stage of an accepted
     * step equals f at the new point, which is both the next step's first stage
     * and the node slope — so it is computed once, not three times.  `fsal_cur`
     * holds f(t, Y) for the current step (the reusable first stage); the stepper
     * writes f(t+h, Ynew) into `fsal_pending`, which the driver promotes to
     * `fsal_cur` (and the node slope) on acceptance.  fsal==false disables it. */
    bool     fsal;
    double*  fsal_cur;      /* f(t, Y): reusable stage-1 (length d)  */
    double*  fsal_pending;  /* f(t+h, Ynew): stage-7 (length d)      */

    /* result assembly: mapping reduced components -> dependent functions */
    size_t   nfun;          /* number of dependent functions                  */
    const char** fun_names; /* array[nfun] interned names (borrowed)          */
    size_t*  fun_state0;    /* array[nfun] reduced index of each function's y  */
    bool     fun_applied;   /* NDSolve[..,u[x],..] -> return u[x] -> IF form   */
    bool     realified;     /* complex system split into Re/Im (dim doubled)   */
    size_t   d_complex;     /* original complex dimension when realified       */
};

/* Evaluate the reduced RHS at (t, Y): writes out[0..d-1] = f(t, Y).
 * Returns false if any component fails to numericalize to a finite double. */
bool nd_rhs_real(NdProblem* P, double t, const double* Y, double* out);

/* Evaluate the Jacobian J[i*d+j] = df_i/dY_j at (t, Y) into a dense row-major
 * buffer (d*d doubles).  Uses the symbolic Jacobian where available and a
 * central finite-difference fallback otherwise.  Returns false on failure. */
bool nd_jacobian_real(NdProblem* P, double t, const double* Y, double* Jout);

/* ------------------------------------------------------------------ *
 *  Options                                                            *
 * ------------------------------------------------------------------ */
typedef struct {
    const char* method;        /* interned method name, or NULL = Automatic   */
    Expr*       method_subopts;/* borrowed List of Rule suboptions, or NULL   */
    NumericSpec spec;          /* machine vs MPFR                             */
    long        wp_bits;       /* working precision in bits (MPFR)            */
    double      acc_goal;      /* digits; -1 = Automatic; HUGE_VAL = Infinity  */
    double      prec_goal;     /* digits; -1 = Automatic; HUGE_VAL = Infinity  */
    int64_t     max_steps;     /* -1 = Automatic                             */
    double      max_step_size; /* <= 0 => unbounded                          */
    double      max_step_fraction; /* default 1/10                           */
    double      starting_step; /* <= 0 => automatic (Hairer)                 */
    int         interp_order;  /* -1 => Automatic (3, cubic Hermite)          */
    Expr*       step_monitor;  /* borrowed, or NULL                          */
    Expr*       eval_monitor;  /* borrowed, or NULL                          */
    Expr*       norm_function; /* borrowed, or NULL (reserved)               */
} NdOpts;

void nd_opts_default(NdOpts* o);

/* Resolved tolerances from goals + working precision. */
typedef struct { double rtol, atol; } NdTol;
NdTol nd_resolve_tol(const NdOpts* o);

/* ------------------------------------------------------------------ *
 *  Stepper vtable (method axis)                                       *
 * ------------------------------------------------------------------ */
typedef enum {
    ND_ADAPTIVE  = 1u << 0,  /* single_step fills Yerr (embedded estimate)    */
    ND_IMPLICIT  = 1u << 1,  /* needs Jacobian + Newton per step              */
    ND_MULTISTEP = 1u << 2,  /* keeps history; needs startup + reset          */
    ND_FSAL      = 1u << 3   /* first-same-as-last: K[0..d) reusable next step */
} NdStepperFlags;

typedef struct NdStepper NdStepper;

struct NdStepper {
    const char* name;
    unsigned    flags;
    int         order;      /* order p of the propagated solution             */
    int         err_order;  /* order of the embedded/effective error estimate */
    int         stages;     /* s: K scratch length is stages*d                */

    /* Advance one step of size h from (t, Y).  Writes Ynew[0..d).  If
     * ND_ADAPTIVE, also writes Yerr[0..d) (per-component local error).
     * K (length stages*d, may be NULL if stages==0) may cache stage
     * derivatives for FSAL / dense output.  Returns false to signal a failed
     * step (the driver then reduces h). */
    bool (*single_step)(const NdStepper* S, NdProblem* P,
                        double t, const double* Y, double h,
                        double* Ynew, double* Yerr, double* K);

    const char* doc;        /* docstring for NDSolve`Name                     */
};

/* Registry: look up a stepper by (interned or plain) method name.  NULL if
 * unknown.  "Automatic"/NULL resolves to the default non-stiff method. */
const NdStepper* nd_lookup_stepper(const char* name);
const NdStepper* nd_default_stepper(void);
/* Enumerate all steppers (for registration of NDSolve`Name builtins). */
size_t           nd_stepper_count(void);
const NdStepper* nd_stepper_at(size_t i);

/* Steppers defined in the per-method modules. */
extern const NdStepper nd_stepper_explicit_euler;
extern const NdStepper nd_stepper_explicit_midpoint;
extern const NdStepper nd_stepper_rk4;
extern const NdStepper nd_stepper_dopri5;
extern const NdStepper nd_stepper_backward_euler;
extern const NdStepper nd_stepper_implicit_trapezoid;
extern const NdStepper nd_stepper_bdf;
extern const NdStepper nd_stepper_adams;

/* ------------------------------------------------------------------ *
 *  Shared numeric helpers used by steppers                            *
 * ------------------------------------------------------------------ */

/* Weighted RMS norm of e against the per-component scale sc_i = atol+rtol|y_i|,
 * blending the current and proposed states y and ynew for the scale. */
double nd_wrms_norm(size_t d, const double* e, const double* y,
                    const double* ynew, NdTol tol);

/* Dense linear solve A x = b (A is n*n row-major, destroyed; b length n,
 * overwritten with the solution).  Gaussian elimination, partial pivoting.
 * Returns false iff A is singular. */
bool nd_dense_solve(size_t n, double* A, double* b);

/* Newton solve of the implicit stage residual for the theta-method family used
 * by the implicit steppers:  Ynew = Ybase + h*theta*f(t1, Ynew) + rhs_const,
 * i.e. solves  G(Z) = Z - Ybase - h*theta*f(t1, Z) - rhs_const = 0.
 * `rhs_const` may be NULL (treated as 0).  Writes the converged Z into Ynew,
 * returns false on non-convergence.  Reused by BackwardEuler (theta=1,
 * Ybase=Yn), ImplicitTrapezoid (theta=1/2, rhs_const=h/2 f(tn,Yn)) and BDF. */
bool nd_newton_theta(NdProblem* P, double t1, const double* Ybase,
                     double h, double theta, const double* rhs_const,
                     const double* Zguess, double* Ynew, NdTol tol);

/* ------------------------------------------------------------------ *
 *  Solution accumulation + adaptive driver                            *
 * ------------------------------------------------------------------ */
typedef struct {
    size_t  n, cap;
    size_t  d;
    double* ts;   /* n            */
    double* Ys;   /* n*d          */
    double* dYs;  /* n*d (slopes) */
} NdSolution;

void nd_solution_init(NdSolution* s, size_t d);
void nd_solution_free(NdSolution* s);
/* Append an accepted node (t, Y[d], Y'[d]).  Used by the multistep drivers. */
void nd_solution_push(NdSolution* s, double t, const double* Y, const double* dY);

typedef enum {
    ND_OK = 0,
    ND_ERR_MAXSTEPS,     /* MaxSteps exhausted (partial solution returned)     */
    ND_ERR_STEPSIZE,     /* step size collapsed (singularity/stiffness)        */
    ND_ERR_NONCONV,      /* implicit Newton failed to converge                 */
    ND_ERR_SAMPLE        /* RHS could not be evaluated at t0                    */
} NdStatus;

/* Integrate the problem with the given stepper and options, accumulating
 * accepted nodes (t_i, Y_i, Y'_i) into `sol` over [tmin, tmax].  Handles both
 * embedded-adaptive and fixed-order (step-doubling) error control. */
NdStatus nd_integrate(NdProblem* P, const NdStepper* S, const NdOpts* o,
                      NdSolution* sol);

/* Multistep integrators (BDF, Adams) manage their own history and run a
 * dedicated fixed-step loop; nd_integrate dispatches to these for the
 * ND_MULTISTEP steppers.  The initial node must NOT be pre-recorded. */
NdStatus nd_multistep_bdf(NdProblem* P, const NdOpts* o, NdSolution* sol);
NdStatus nd_multistep_adams(NdProblem* P, const NdOpts* o, NdSolution* sol);

/* Shared: choose a fixed step for the multistep methods (Hairer initial step,
 * clamped).  `dir` is +1 (toward tmax) or -1 (toward tmin). */
double nd_fixed_step(NdProblem* P, const NdOpts* o, NdTol tol, double dir);

/* Build the NDSolve result: a rule list {u -> InterpolatingFunction[...], ...}
 * (or {u[x] -> IF[x], ...} when P->fun_applied) from the accumulated solution.
 * Consumes nothing; caller owns the returned Expr*. */
Expr* nd_build_result(NdProblem* P, const NdOpts* o, const NdSolution* sol);

#ifdef USE_MPFR
/* Arbitrary-precision integrator (state, time and step size all MPFR).  Used
 * when WorkingPrecision exceeds MachinePrecision.  Runs the adaptive DOPRI5
 * pair (or fixed RK4 when the requested method is RK4) and returns the assembled
 * MPFR InterpolatingFunction rule list, or NULL on failure. */
Expr* nd_solve_mpfr(NdProblem* P, const NdOpts* o, const NdStepper* S);
#endif

#endif /* NDSOLVE_COMMON_H */
