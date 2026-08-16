/*
 * findmin_internal.h
 *
 * Shared declarations for the split FindMinimum / FindMaximum /
 * NMinimize / NMaximize implementation. findmin.c was one 9.7k-line
 * translation unit; it is now one file per optimizer plus two shared
 * cores (findmin_common.c local, findmin_nm_common.c global). This
 * header carries the shared types, tuning constants, the file-static
 * global registries (as extern), and the prototypes for every symbol
 * used across those translation units. Method-internal helpers stay
 * static in their own file and are NOT declared here.
 *
 * Mirrors the ndsolve_common.h pattern in this directory.
 */
#ifndef FINDMIN_INTERNAL_H
#define FINDMIN_INTERNAL_H

#include "findmin.h"

#include <math.h>
/* M_PI is POSIX, not C99: glibc hides it under -std=c99. See CLAUDE.md §10. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPFR
#  include <mpfr.h>
#endif
#include "arithmetic.h"
#include "attr.h"
#include "compile/compile.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* ================= shared types ================= */


/* ------------------------------------------------------------------ *
 *  Types                                                              *
 * ------------------------------------------------------------------ */

typedef enum {
    FM_METHOD_AUTOMATIC = 0,
    FM_METHOD_BRENT,
    FM_METHOD_QUASINEWTON,   /* BFGS */
    FM_METHOD_CONJGRAD,      /* Polak-Ribière+ */
    FM_METHOD_NEWTON,
    FM_METHOD_LBFGSB,        /* limited-memory BFGS + bounds (Mathilda ext.) */
    FM_METHOD_POWELL,        /* derivative-free conjugate directions */
    FM_METHOD_NELDERMEAD,    /* derivative-free downhill simplex */
    FM_METHOD_TNC,           /* truncated Newton (Hessian-free) + bounds */
    FM_METHOD_SLSQP,         /* sequential least-squares QP (constrained SQP) */
    FM_METHOD_COBYLA,        /* derivative-free linear-approximation (constrained) */
    FM_METHOD_COBYQA,        /* derivative-free quadratic-approximation (constrained) */
    FM_METHOD_NEWTONCG,      /* line-search truncated Newton (Hessian-free CG)       */
    FM_METHOD_DOGLEG,        /* trust-region Powell dogleg (needs SPD Hessian)       */
    FM_METHOD_TRUSTNCG,      /* trust-region Steihaug-Toint CG (Hessian-free)        */
    FM_METHOD_TRUSTEXACT,    /* trust-region Moré-Sorensen (exact subproblem)        */
    FM_METHOD_TRUSTKRYLOV    /* trust-region GLTR/Lanczos (Hessian-free)             */
} FmMethod;

typedef enum {
    FM_PREC_MACHINE = 0
#ifdef USE_MPFR
    , FM_PREC_MPFR
#endif
} FmPrecMode;

typedef enum {
    FM_SPEC_BAD = 0,
    FM_SPEC_VAR_ONLY,   /* {x}                       -> x0 = 1 (MMA-default) */
    FM_SPEC_SINGLE,     /* {x, x0}                                           */
    FM_SPEC_TWO_START,  /* {x, x0, x1}               -> derivative free      */
    FM_SPEC_BRACKET     /* {x, xstart, xmin, xmax}                           */
} FmSpecKind;

typedef struct {
    FmMethod   method;
    FmPrecMode prec_mode;
    long       wp_bits;          /* MPFR bits when prec_mode == MPFR    */
    int64_t    max_iter;         /* default 500                          */
    bool       max_iter_set;     /* true if MaxIterations given explicitly */
    double     acc_goal_digits;  /* filled at WP/2 if Automatic          */
    double     prec_goal_digits;
    Expr*      gradient;         /* borrowed; user-supplied list or NULL */
    Expr*      step_monitor;     /* borrowed; or NULL                    */
    Expr*      eval_monitor;     /* borrowed; or NULL                    */
} FmOpts;

/* Per-variable snapshot for Block-style binding. */
typedef struct {
    const char* name;            /* interned                            */
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} FmVarBind;

/* Per-variable box (used both for 4-elt specs and parsed inequalities).
 * `has_lo`/`has_hi` flag presence. */
typedef struct {
    bool   has_lo, has_hi;
    double lo, hi;
} FmBox;

/* General (non-box) inequality g(x) <= 0 or equality h(x) == 0. The
 * objective during the outer μ loop is f(x) + μ * Σ max(0,g_i)^2
 *                                              + μ * Σ h_j^2.
 * For the inner solver we evaluate each constraint expression directly,
 * and (when present) its symbolic gradient — needed so the augmented
 * objective is differentiated consistently with its value. */
typedef struct {
    Expr*  expr;       /* owned: feasible ≡ (expr <= 0) or (expr == 0)        */
    Expr** grad_exprs; /* owned: ∇expr w.r.t. vars (length n), or NULL → FD   */
    bool   equality;   /* true → equality constraint                          */
} FmGenCon;

/* A disjunctive (Or) constraint: feasible ≡ at least ONE branch is feasible.
 * Stored as its boolean-of-comparisons subtree over the effective variables.
 * The penalty contribution is the MINIMUM branch penalty (fm_bool_penalty), so
 * a point satisfying any one branch scores zero and the (total ≤ NM_FEAS_EPS)
 * feasibility test still means "feasible". Because the min() is non-smooth it is
 * consumed only by the derivative-free global search (NMinimize); the smooth
 * local polish leaves disjunction feasibility to the Deb accept/reject gate, and
 * FindMinimum's gradient penalty method rejects Or outright. Evaluated by the
 * interpreter rather than a compiled program — disjunctive problems are rare. */
typedef struct {
    Expr*  expr;       /* owned: the Or[...] constraint subtree                */
} FmDisjunction;

/* A "one-hot" / assignment group: an equality constraint of the form
 * x_{i1} + x_{i2} + ... + x_{im} == k over DISTINCT binary optimization
 * variables (e.g. each row/column of an assignment matrix summing to 1). The
 * integer local search repairs each group to exactly k ones so the 2-flip swap
 * move can then navigate to full feasibility — the structure single-coordinate
 * descent cannot reach from a random rounded start. `idx` is owned. */
typedef struct {
    size_t* idx;       /* owned: indices (into the effective variable list)     */
    size_t  len;
    int     k;         /* target sum                                            */
} NmOneHot;

/* Limited-memory state: the m most recent correction pairs stored in circular
 * buffers S,Y (each pair occupies one length-n row), rho_i = 1/(y_i·s_i), and
 * the H0 scaling gamma = (s_last·y_last)/(y_last·y_last). `cnt` counts valid
 * pairs (≤ m); `head` is the slot the NEXT pair overwrites (the oldest, once
 * the ring is full). */
typedef struct {
    size_t  n, m, cnt, head;
    double* S;      /* m*n */
    double* Y;      /* m*n */
    double* rho;    /* m   */
    double  gamma;
} FmLbfgsMem;

/* --- Backtracking line search + active-set bound handling -------------- *
 *
 * L-BFGS needs a line search that tries the quasi-Newton unit step (alpha = 1)
 * FIRST — the shared fm_line_search caps the initial step at 1/||d|| for
 * steepest descent, which throttles the well-scaled L-BFGS direction and stalls
 * convergence on ill-conditioned problems. Bounds are handled by an active-set
 * projection:
 * a coordinate resting on a box face with the gradient pushing outward is
 * fixed (masked out of the two-loop and the direction), so the free variables
 * optimise in the reduced subspace; the search is capped at the step that
 * reaches the nearest remaining face, at which point that coordinate joins the
 * active set. This is the practical active-set L-BFGS-B; it reaches the same
 * optima as the reference generalized-Cauchy-point algorithm. */

/* Immutable per-solve context, so the line-search helpers stay tidy. */
typedef struct {
    Expr* f; Expr** g_exprs; FmVarBind* binds; size_t n;
    const FmGenCon* gens; size_t ngens; double mu;
    const FmBox* boxes; const FmOpts* opts; bool augmented;
} FmLbfgsCtx;

/* Quadratic-model operator handed to a subproblem solver: either a dense B
 * (dogleg / trust-exact and the GLTR tridiagonal) or a Hessian-vector context
 * (the CG-based methods). fm_quad_matvec hides which one. */
typedef struct {
    size_t              n;
    const double*       B;        /* dense n·n, or NULL → use the Hv path      */
    const FmLbfgsCtx*   c;        /* Hv context (NULL for a pure dense B)       */
    const double*       xbase;    /* Hv base point                             */
    const double*       gbase;    /* Hv base gradient = ∇f(xbase)              */
    double*             xpert;    /* Hv scratch (n)                            */
    double*             gpert;    /* Hv scratch (n)                            */
    const bool*         active0;  /* all-false mask (n): unconstrained subspace */
} FmQuad;

/* A subproblem solver: given g and Δ, fill the step p and set *hits when the
 * step reaches the trust boundary. The driver recovers pᵀBp itself, so solvers
 * need not return it. Returns false only on an unrecoverable failure. */
typedef bool (*FmSubSolver)(const FmQuad* q, const double* g, double gnorm,
                            double Delta, double* p, bool* hits);

enum { NM_AUTO = 0, NM_DE, NM_NELDERMEAD, NM_RANDOMSEARCH, NM_SA, NM_SHGO,
       NM_DUAL_ANNEALING, NM_DIRECT };

typedef struct {
    int      method;         /* NM_AUTO / NM_DE / ...                        */
    int      search_points;  /* 0 ⇒ auto                                     */
    double   F;              /* DE scaling factor;   <0 ⇒ auto               */
    double   CR;             /* DE crossover prob.;  <0 ⇒ auto               */
    double   reflect_ratio;  /* NelderMead reflection coeff;  <0 ⇒ default 1  */
    double   expand_ratio;   /* NelderMead expansion coeff;   <0 ⇒ default 2  */
    double   contract_ratio; /* NelderMead contraction coeff; <0 ⇒ default .5 */
    double   shrink_ratio;   /* NelderMead shrink coeff;      <0 ⇒ default .5 */
    double   tolerance;      /* simplex convergence tolerance; <0 ⇒ default    */
    int      post_process;   /* -1 auto (on) / 1 on / 0 off (skip polish)    */
    Expr*    init_points;    /* "InitialPoints" -> {{...},...}, borrowed / NULL */
    Expr*    penalty_fn;     /* "PenaltyFunction" -> f, borrowed / NULL ⇒ auto */
    double   perturb_scale;  /* SA "PerturbationScale"; <0 ⇒ default 1.0       */
    Expr*    boltzmann_fn;   /* SA "BoltzmannExponent" -> f, borrowed / NULL ⇒ auto */
    int      level_iterations;/* SA "LevelIterations"; <=0 ⇒ auto (50)         */
    int      shgo_sampling;  /* SHGO "SamplingMethod": 0 Simplicial / 1 Sobol / 2 Halton */
    int      shgo_iters;     /* SHGO "Iterations" (scipy iters); <=0 ⇒ auto (1) */
    double   da_visit;       /* DualAnnealing "VisitingParameter" qv; <=0 ⇒ default 2.62 */
    double   da_accept;      /* DualAnnealing "AcceptanceParameter" qa; 0 ⇒ default -5.0 */
    double   da_init_temp;   /* DualAnnealing "InitialTemperature"; <=0 ⇒ default 5230  */
    double   da_restart_ratio;/* DualAnnealing "RestartTemperatureRatio"; <=0 ⇒ default */
    int      da_local_search;/* DualAnnealing "LocalSearch": -1 auto(on) / 1 on / 0 off  */
    /* DIRECT (DIviding RECTangles), mirroring scipy.optimize.direct. Sentinels:
     * <0 / <=0 ⇒ take the scipy default (see NM_DIRECT_* macros / nm_direct). */
    int      direct_locally_biased; /* "LocallyBiased": -1 auto(on, DIRECT-L) / 1 on / 0 off (original DIRECT) */
    double   direct_eps;     /* "Epsilon": potentially-optimal slack; <0 ⇒ 1e-4    */
    int      direct_max_fun; /* "MaxFunctionEvaluations"; <=0 ⇒ auto 1000*n (capped) */
    int      direct_max_iter;/* "MaxIterations"; <=0 ⇒ default 1000                 */
    double   direct_vol_tol; /* "VolumeTolerance"; <0 ⇒ 1e-16                        */
    double   direct_len_tol; /* "LengthTolerance"; <0 ⇒ 1e-6                         */
    double   direct_fmin;    /* "MinValue": known global minimum; -inf ⇒ inactive   */
    double   direct_fmin_rtol;/* "MinValueTolerance": rel. tol vs MinValue; <0 ⇒ 1e-4 */
    uint64_t seed;
} NmConfig;

/* SplitMix64 — a small, fast, fully deterministic PRNG (mirrors mcint). A
 * fixed default seed makes the stochastic search reproducible so the unit
 * tests are stable; "RandomSeed" overrides it. */
typedef struct { uint64_t s; } NmRng;

/* Block-localise a set of head symbols: snapshot and clear their Own/DownValues
 * so an objective/constraint/variable expression that mentions x[i] evaluates
 * with x free — expanding Table/Sum symbolically without capturing any stray
 * user definition. Restored by nm_heads_restore. */
typedef struct { const char* name; Rule* own; Rule* down; bool valid; } NmHeadSave;

/* ------------------------------------------------------------------ *
 *  Driver context bundle + point evaluation / comparison             *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*        f_raw;          /* objective (borrowed)                     */
    Expr**       vars;           /* variable symbols (borrowed)             */
    size_t       n;
    FmVarBind*   binds;
    Expr**       g_exprs;        /* symbolic ∇f, or NULL → finite diff       */
    FmGenCon*    gens;
    size_t       ngens;
    FmBox*       boxes;
    const FmOpts* opts;
    const bool*  is_int;
    bool         any_int;
    /* True iff at least one CONTINUOUS optimization variable appears in the
     * objective. When false, the continuous variables are constraint-only helpers
     * (e.g. MTZ ordering variables in a TSP), so the expensive continuous-
     * relaxation recovery in nm_local_polish — which exists to optimize/rescue
     * continuous objective variables — is pure overhead and is skipped. */
    bool         cont_in_obj;
    const double* reg_lo;
    const double* reg_hi;
    /* Machine-precision fast path: a bytecode-compiled objective and per-general-
     * constraint programs over the effective variables, or NULL to use the
     * interpreter. Each point evaluation prefers the compiled program and falls
     * back to the interpreter when it is absent or reports a domain/non-finite
     * result — so the compiled path is a pure speedup, never a correctness risk. */
    CompiledProgram*  f_prog;    /* objective, or NULL                          */
    CompiledProgram** g_progs;   /* per-constraint (len ngens), entries may NULL */
    /* "PenaltyFunction" -> f. Borrowed / NULL ⇒ Automatic (built-in squared
     * penalty). When set, the global-search violation of each active constraint
     * is scored as f[violation] rather than violation^2; see nm_eval_pen. */
    Expr*             penalty_fn;
    /* Disjunctive (Or) constraints. Each contributes min-over-branch penalty via
     * the interpreter (fm_bool_penalty); len ndisj, borrowed from the setup. */
    FmDisjunction*    disj;
    size_t            ndisj;
    /* One-hot / assignment groups detected from the equality constraints; owned,
     * freed at driver cleanup. Empty for problems without assignment structure. */
    NmOneHot*         onehots;
    size_t            n_onehots;
} NmDriver;

/* ============================================================================
 *  Dual Annealing — Generalized Simulated Annealing (Tsallis & Stariolo,
 *  Physica A 233 (1996) 395-406; Xiang, Sun, Fan & Gong, Phys. Lett. A 233
 *  (1997) 216-220) with a local search after each Markov chain, mirroring
 *  scipy.optimize.dual_annealing.
 *
 *  Two Tsallis q-generalizations of the classical distributions:
 *    - the VISITING distribution (shape qv) draws a heavy-tailed Cauchy-Lorentz
 *      jump whose scale is set by the visiting temperature T(si). qv -> 1 is
 *      Gaussian (classical SA), qv = 2 is Cauchy (fast SA), qv > 2 has the heavy
 *      tails that let the walk escape deep basins. Jumps beyond a tail limit are
 *      folded back into the box, so a hot chain samples the box near-uniformly.
 *    - the ACCEPTANCE distribution (shape qa) is a generalized Metropolis rule
 *      on the penalized energy, with acceptance temperature T/(step+1).
 *  The temperature follows the GSA cooling schedule and REANNEALS (resets to T0)
 *  once it drops below T0 * restart_ratio, so one chain keeps exploring across
 *  the whole budget. After each chain a local search (Mathilda's shared
 *  nm_local_polish: BFGS / augmented-Lagrangian / integer descent) refines the
 *  incumbent — the second "annealing" that makes the method dual. Deb's
 *  feasibility rules (nm_better) select the reported best, so box + general
 *  constraints + integer domains work with no Dual-Annealing-specific code.
 *
 *  "SearchPoints" -> K runs K independent chains (default 1, matching scipy's
 *  single chain) and keeps the Deb-best. Sub-options: "VisitingParameter" (qv),
 *  "AcceptanceParameter" (qa), "InitialTemperature" (T0),
 *  "RestartTemperatureRatio", "LocalSearch" (on/off), "MaxIterations" (the
 *  per-chain temperature-step budget), "RandomSeed", "PostProcess".
 * ========================================================================== */

/* Visiting-distribution factors that depend only on qv (precomputed once). */
typedef struct {
    double qv;
    double factor4p;   /* sqrt(pi) * factor2 / (factor3 * (3 - qv))           */
    double factor6;    /* pi (1-f5) / sin(pi (1-f5)) / Gamma(2 - f5)          */
} DaFactors;

/* ---- The complex: vertex set + edge list (the graph 1-skeleton) ------------ */

typedef struct {
    size_t  n;          /* dimension                                            */
    double* V;          /* nv * n vertex coordinates (row-major)                */
    double* fv;         /* nv objective values                                  */
    double* pv;         /* nv constraint-violation penalties                    */
    size_t  nv, vcap;
    size_t* E;          /* ne * 2 edge endpoints (duplicates harmless)          */
    size_t  ne, ecap;
} ShgoCx;

/* Variable set: symbols (borrowed), integer-domain mask, and per-variable
 * search-interval hints parsed from {x, lo, hi} / {x, x0, lo, hi} specs. */
typedef struct {
    size_t  n;
    Expr**  vars;
    bool*   is_int;
    double* rlo;
    double* rhi;
    bool*   has_rlo;
    bool*   has_rhi;
    bool    any_int;
} NmVarSet;

/* ============ shared tuning constants ============ */


/* ------------------------------------------------------------------ *
 *  Brent's minimisation (1D, machine precision)                        *
 * ------------------------------------------------------------------ */

#define FM_CGOLD 0.3819660112501051
#define FM_ZEPS  1.0e-12

/* ================================================================== *
 *  Trust-region method family                                         *
 *  (Newton-CG, dogleg, trust-ncg, trust-exact, trust-krylov)          *
 * ================================================================== *
 *
 * Five UNCONSTRAINED smooth minimizers that mirror the corresponding
 * scipy.optimize.minimize methods. Four are genuine trust-region methods
 * sharing one driver (fm_run_trust_region) and differing only in how they
 * solve the trust-region subproblem  min_{||p||<=Δ}  g·p + ½ pᵀB p:
 *
 *   dogleg        Powell dogleg path (needs a positive-definite Hessian B)
 *   trust-ncg     Steihaug-Toint truncated CG (Hessian-vector products only)
 *   trust-exact   Moré-Sorensen near-exact subproblem (indefinite B / hard case)
 *   trust-krylov  GLTR: Lanczos-tridiagonalize B in the Krylov space of g, then
 *                 solve the small tridiagonal subproblem exactly (Moré-Sorensen)
 *
 * Newton-CG is line-search (not trust-region) and gets its own runner.
 *
 * Curvature access splits the family: dogleg and trust-exact form the full
 * Hessian (symbolic via fm_eval_hessian, else finite-differenced column by
 * column from the compiled gradient); the other three never form B and use
 * Hessian-vector products Hv ≈ (∇f(x+hv)-∇f(x))/h — the exact FD kernel
 * fm_tnc_hessvec already used by TNC. scipy's trust-region methods are all
 * unconstrained, so the dispatcher rejects general (non-box) constraints for
 * these five exactly as Powell/NelderMead do; the runners assume ngens == 0.
 *
 * Constants match scipy._trustregion: Δ0 = 1, Δmax = 1000, accept when the
 * actual/predicted-reduction ratio ρ > 0.15, shrink Δ by ¼ when ρ < ¼, expand
 * by 2 when ρ > ¾ AND the step reached the boundary. */

#define FM_TR_DELTA0   1.0
#define FM_TR_DELTAMAX 1000.0
#define FM_TR_ETA      0.15

/* ================================================================== *
 *  NMinimize / NMaximize — numerical GLOBAL optimization             *
 * ================================================================== *
 *
 *  Layered directly on the FindMinimum machinery above. NMinimize reuses
 *  the same Block-style variable binding (fm_bind_*), objective and
 *  constraint evaluation (fm_eval_scalar / fm_eval_penalty), constraint
 *  parsing (fm_collect_constraints → boxes + general FmGenCon[] +
 *  disjunctions), local polishers (fm_run_bfgs / fm_run_penalty, plus
 *  fm_run_bfgs_mpfr for WorkingPrecision), and the result builders. What is new
 *  here is the stochastic *global* search — DifferentialEvolution (Automatic
 *  default), NelderMead, RandomSearch, SimulatedAnnealing — mixed-integer
 *  handling via Element[x, Integers], and the infeasible {Infinity, ...} return.
 *
 *  Constraints are handled with Deb's feasibility rules during the global
 *  search (no penalty-weight tuning): a feasible point always beats an
 *  infeasible one; among feasible points the smaller objective wins; among
 *  infeasible points the smaller total violation wins. The global best is
 *  then polished with the exact local solver and re-checked for feasibility.
 *
 *  Disjunctive (Or) constraints — feasible iff at least one branch holds — are
 *  supported here (but not in FindMinimum's smooth gradient path): each adds its
 *  minimum-branch penalty (fm_bool_penalty), 0 when any branch is satisfied, so
 *  the same Deb gate selects it with no extra tuning. The derivative-free global
 *  search consumes the non-smooth min directly; the local polish optimises only
 *  the conjunctive part and the post-polish feasibility gate keeps the result
 *  inside the disjunctive-feasible region.
 */

#define NM_DEFAULT_SPAN   10.0     /* half-width of the default search box   */
#define NM_BOUND_SPAN     20.0     /* span added when only one bound is known*/
#define NM_FEAS_EPS       1.0e-8   /* penalty ≤ this ⇒ feasible (selection)  */
#define NM_FEAS_FINAL     1.0e-6   /* final feasible-vs-Infinity threshold   */
#define NM_PENALTY_MU     1.0e6    /* fixed penalty weight for NelderMead/SA */
#define NM_SA_TOTAL_CAP   120000   /* SA aggregate iteration cap across chains
                                    * when "SearchPoints" -> K > 1 restarts     */
#define NM_SA_BURN_IN     30       /* trial steps probed per chain to measure
                                    * the objective scale for the temperature   */
#define NM_DEFAULT_SEED   20260814ULL
#define NM_MAX_REGION_EXPAND 4     /* infeasible-region rescue: grow +-SPAN by
                                    * 10^k for k = 1..this (up to +-1e5)       */

/* Dual Annealing (Generalized Simulated Annealing) defaults, matching
 * scipy.optimize.dual_annealing. See nm_dual_annealing. */
#define NM_DA_VISIT          2.62      /* qv: visiting-distribution shape, (1,3] */
#define NM_DA_ACCEPT        (-5.0)     /* qa: acceptance-distribution shape      */
#define NM_DA_INIT_TEMP      5230.0    /* T0: initial (and post-reanneal) temp   */
#define NM_DA_RESTART_RATIO  2.0e-5    /* reanneal when T < T0 * this            */
#define NM_DA_MAXITER        1000      /* temperature-step budget (scipy maxiter)*/
#define NM_DA_MAXFUN         10000000LL/* objective-evaluation safety cap        */
#define NM_DA_STALL_MAX      1000      /* forced local search after this many
                                        * non-improving chains (scipy default)   */
#define NM_DA_TAIL_LIMIT     1.0e8     /* clamp bound on a raw visiting jump      */
#define NM_DA_MIN_VISIT      1.0e-10   /* nudge off the lower bound after wrap    */

/* DIRECT (DIviding RECTangles) defaults, matching scipy.optimize.direct.
 * See nm_direct. */
#define NM_DIRECT_EPS        1.0e-4    /* eps: potentially-optimal slack          */
#define NM_DIRECT_MAXITER    1000      /* maxiter: division-round budget          */
#define NM_DIRECT_MAXFUN_PER_DIM 1000  /* auto maxfun = this * n (scipy rule)     */
#define NM_DIRECT_MAXFUN_CAP 20000000  /* absolute maxfun ceiling (~1 GiB store)  */
#define NM_DIRECT_VOLTOL     1.0e-16   /* stop when incumbent cell volume < this  */
#define NM_DIRECT_LENTOL     1.0e-6    /* stop when incumbent cell half-side < this*/
#define NM_DIRECT_FMINRTOL   1.0e-4    /* rel. tol for the MinValue early stop     */

/* ===== file-static global registries (defined in findmin_common.c) ===== */
extern bool g_fm_quiet;
extern const char* g_fm_name;
extern const double* g_fm_al_lambda;
extern const FmGenCon* g_fm_al_gens;
extern Expr* g_fm_obj_expr;
extern CompiledProgram* g_fm_obj_prog;
extern size_t g_fm_obj_nargs;
extern Expr** g_fm_grad_exprs;
extern CompiledProgram** g_fm_grad_progs;
extern size_t g_fm_grad_n;

/* ============ cross-file function prototypes ============ */

void fm_warn(const char* fn, const char* tag, const char* fmt, ...);
bool fm_expr_to_double_real(Expr* e, double* out);
bool fm_is_option_arg(Expr* e);
bool fm_parse_working_precision(Expr* val,
                                       FmPrecMode* mode, long* bits);
bool fm_parse_goal(Expr* val, double* digits_out);
bool fm_apply_option(Expr* rule, FmOpts* opts);
void fm_bind_snapshot(FmVarBind* b, const char* name);
void fm_bind_clear_temp(FmVarBind* b);
void fm_bind_restore(FmVarBind* b);
void fm_fire_monitor(Expr* monitor);
Expr* fm_eval_with_bindings(Expr* f, FmVarBind* binds,
                                   Expr* const* values, size_t n,
                                   Expr* eval_monitor,
                                   NumericSpec spec);
bool fm_eval_scalar(Expr* f, FmVarBind* binds,
                           const double* x, size_t n,
                           const FmOpts* opts, double* out);
Expr** fm_compute_gradient(Expr* f, Expr** vars, size_t n);
Expr*** fm_compute_hessian(Expr* f, Expr** vars, size_t n);
bool fm_grad_finite_diff(Expr* f, FmVarBind* binds,
                                const double* x, size_t n,
                                const FmOpts* opts, double* g_out);
bool fm_eval_gradient(Expr** g_exprs, FmVarBind* binds,
                             const double* x, size_t n,
                             const FmOpts* opts, double* g_out);
bool fm_eval_hessian(Expr*** H_exprs, FmVarBind* binds,
                            const double* x, size_t n,
                            const FmOpts* opts, double* H_out /* n*n */);
Expr* fm_build_result(double fmin, Expr** vars, const double* vals,
                             size_t n);
Expr* fm_build_result_mpfr(const mpfr_t fmin, Expr** vars,
                                  mpfr_t const* vals, size_t n);
mpfr_t* fm_mpfr_array(size_t count, long bits);
void fm_mpfr_array_free(mpfr_t* arr, size_t count);
bool fm_run_bfgs_mpfr(Expr* f, Expr** vars, size_t n,
                             FmVarBind* binds, Expr** g_exprs,
                             mpfr_t* x, /* in/out */
                             const FmBox* boxes,
                             const FmOpts* opts,
                             mpfr_t fx_out);
FmSpecKind fm_parse_var_spec(Expr* spec, Expr** var_out,
                                    Expr** x0_out, Expr** x1_out,
                                    Expr** xmin_out, Expr** xmax_out);
bool fm_constraint_to_g(Expr* cmp, Expr** expr_out, bool* equality_out);
bool fm_bool_supported(Expr* c);
bool fm_collect_constraints(Expr* cons, Expr** vars, size_t nvars,
                                   FmBox* boxes,
                                   FmGenCon** gens_inout, size_t* ngens_inout,
                                   size_t* gens_cap_inout,
                                   FmDisjunction** disj_inout,
                                   size_t* ndisj_inout, size_t* disj_cap_inout);
void fm_project_box(double* x, size_t n, const FmBox* boxes);
bool fm_eval_penalty(const FmGenCon* gens, size_t ngens,
                            FmVarBind* binds, const double* x, size_t n,
                            const FmOpts* opts, double* pen_out);
bool fm_eval_augmented(Expr* f, FmVarBind* binds,
                              const double* x, size_t n,
                              const FmGenCon* gens, size_t ngens,
                              double mu, const FmOpts* opts, double* out);
bool fm_eval_aug_gradient(Expr* f, Expr** g_exprs,
                                 const FmGenCon* gens, size_t ngens,
                                 double mu,
                                 FmVarBind* binds, const double* x, size_t n,
                                 const FmOpts* opts, double* g_out);
bool fm_line_search(Expr* f, FmVarBind* binds, size_t n,
                           const double* x, const double* d,
                           double f0, double g_dot_d,
                           const FmGenCon* gens, size_t ngens, double mu,
                           const FmBox* boxes,
                           const FmOpts* opts,
                           double* alpha_out, double* f_out, double* x_out);
bool fm_bracket(Expr* f, FmVarBind* binds, const FmOpts* opts,
                       double x0, const FmBox* box1,
                       double* a_out, double* b_out, double* c_out);
bool fm_bracket_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
                            const mpfr_t x0, const FmBox* box1,
                            mpfr_t a, mpfr_t b, mpfr_t c,
                            mpfr_t fa, mpfr_t fb, mpfr_t fc);
bool fm_brent_min_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
                              const mpfr_t a_in, const mpfr_t b_in, const mpfr_t c_in,
                              const FmBox* box1,
                              mpfr_t x_out, mpfr_t fx_out);
bool fm_brent_min(Expr* f, FmVarBind* bind, const FmOpts* opts,
                         double a, double b, double c,
                         const FmBox* box1,
                         double* x_out, double* fx_out);
bool fm_chol_factor(double* H, size_t n, double tau);
void fm_chol_solve(const double* L, size_t n,
                          const double* b, double* x);
bool fm_run_bfgs(Expr* f, Expr** vars, size_t n,
                        FmVarBind* binds, Expr** g_exprs,
                        double* x, /* in/out */
                        const FmGenCon* gens, size_t ngens, double mu,
                        const FmBox* boxes,
                        const FmOpts* opts,
                        double* fx_out);
bool fm_run_cg(Expr* f, Expr** vars, size_t n,
                      FmVarBind* binds, Expr** g_exprs,
                      double* x,
                      const FmGenCon* gens, size_t ngens, double mu,
                      const FmBox* boxes,
                      const FmOpts* opts,
                      double* fx_out);
bool fm_run_newton(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs, Expr*** H_exprs,
                          double* x,
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out);
bool fm_lbfgs_fg(const FmLbfgsCtx* c, const double* p, double* fval, double* g);
bool fm_lbfgs_linesearch(const FmLbfgsCtx* c, const double* x, const double* d,
                                double f0, double dphi0, double alpha_max,
                                double* xa, double* ga, double* f_out, double* a_out);
bool fm_run_lbfgsb(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out);
bool fm_tnc_grad(const FmLbfgsCtx* c, const double* p, double* g);
bool fm_tnc_hessvec(const FmLbfgsCtx* c, const double* x, const double* g_base,
                           const double* v, const bool* active,
                           double* xpert, double* gpert, double* Hv);
void fm_tnc_cg(const FmLbfgsCtx* c, const double* x, const double* g_base,
                      const double* gm, const bool* active,
                      double* p, double* r, double* d, double* Hd,
                      double* xpert, double* gpert);
bool fm_run_tnc(Expr* f, Expr** vars, size_t n,
                       FmVarBind* binds, Expr** g_exprs,
                       double* x, /* in/out */
                       const FmGenCon* gens, size_t ngens, double mu,
                       const FmBox* boxes,
                       const FmOpts* opts,
                       double* fx_out);
double fm_dot(const double* a, const double* b, size_t n);
void fm_matvec(const double* B, size_t n, const double* v, double* out);
bool fm_quad_matvec(const FmQuad* q, const double* v, double* out);
bool fm_tr_dogleg(const FmQuad* q, const double* g, double gnorm,
                         double Delta, double* p, bool* hits);
bool fm_tr_steihaug(const FmQuad* q, const double* g, double gnorm,
                           double Delta, double* p, bool* hits);
bool fm_tr_moresorensen(const FmQuad* q, const double* g, double gnorm,
                               double Delta, double* p, bool* hits);
bool fm_tr_gltr(const FmQuad* q, const double* g, double gnorm,
                       double Delta, double* p, bool* hits);
bool fm_run_trust_region(Expr* f, Expr** vars, size_t n,
                                FmVarBind* binds, Expr** g_exprs, Expr*** H_exprs,
                                double* x, /* in/out */ const FmBox* boxes,
                                const FmOpts* opts, double* fx_out,
                                FmSubSolver solve_sub, bool needs_dense_B);
bool fm_run_newton_cg(Expr* f, Expr** vars, size_t n, FmVarBind* binds,
                             Expr** g_exprs, double* x, /* in/out */
                             const FmGenCon* gens, size_t ngens, double mu,
                             const FmBox* boxes, const FmOpts* opts, double* fx_out);
int fm_slsqp_activeset(const double* L, size_t nv, const double* g,
                              const double* Nrm, const double* b,
                              const int* is_eq, size_t m,
                              double* d, double* mult);
bool fm_slsqp_qp(const double* L, size_t n, const double* g,
                        const FmGenCon* gens, const double* c,
                        const double* Jac, size_t ngens,
                        const FmBox* boxes, const double* x,
                        double* d, double* lam);
bool fm_run_slsqp(Expr* f, Expr** vars, size_t n,
                         FmVarBind* binds, Expr** g_exprs,
                         double* x, /* in/out */
                         const FmGenCon* gens, size_t ngens, double mu,
                         const FmBox* boxes,
                         const FmOpts* opts,
                         double* fx_out);
bool fm_run_powell(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out);
void fm_nm_sort_idx(size_t* idx, const double* fsim, size_t np);
bool fm_run_neldermead(Expr* f, Expr** vars, size_t n,
                              FmVarBind* binds, Expr** g_exprs,
                              double* x, /* in/out */
                              const FmGenCon* gens, size_t ngens, double mu,
                              const FmBox* boxes,
                              const FmOpts* opts,
                              double* fx_out);
bool fm_cobyla_eval(Expr* f, const FmGenCon* gens, size_t ngens,
                           FmVarBind* binds, const FmOpts* opts,
                           const double* p, size_t n, double* fout, double* gv);
bool fm_run_cobyla(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out);
bool fm_run_cobyqa(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs,
                          double* x, /* in/out */
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out);
bool fm_run_penalty(Expr* f, Expr** vars, size_t n,
                           FmVarBind* binds, FmMethod method,
                           Expr** g_exprs, Expr*** H_exprs,
                           double* x, /* in/out */
                           const FmGenCon* gens, size_t ngens,
                           const FmBox* boxes,
                           const FmOpts* opts,
                           double* fx_out);
Expr* builtin_findminimum(Expr* res);
Expr* builtin_findmaximum(Expr* res);
void nm_rng_seed(NmRng* r, uint64_t seed);
uint64_t nm_rng_next(NmRng* r);
double nm_rng_unif(NmRng* r);
double nm_rng_range(NmRng* r, double lo, double hi);
double nm_rng_normal(NmRng* r);
bool nm_is_head(const Expr* e, const char* sym);
bool nm_is_constraint_tree(const Expr* e);
Expr* nm_subst(Expr* e, Expr* const* from, Expr* const* to, size_t n);
Expr* nm_fresh_symbol(void);
void nm_heads_localize(NmHeadSave* sv, const char** names, size_t m);
void nm_heads_restore(NmHeadSave* sv, size_t m);
bool fm_bool_penalty(Expr* c, FmVarBind* binds, const double* x, size_t n,
                            const FmOpts* opts, Expr* penalty_fn, double* out);
bool nm_eval_pen(NmDriver* D, const double* xr, double* out);
void nm_eval(NmDriver* D, const double* x, double* f_out, double* pen_out);
bool nm_better(double fa, double pa, double fb, double pb);
void nm_project(NmDriver* D, double* x);
double nm_phi(NmDriver* D, const double* x);
void nm_detect_onehots(NmDriver* D);
void nm_free_onehots(NmDriver* D);
bool nm_expr_contains_symbol(const Expr* e, const char* sym);
void nm_local_polish(NmDriver* D, double* x, double* f_out, double* pen_out);
void nm_de(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest);
bool nm_simplex_from_points(NmDriver* D, const Expr* pts, size_t n, double* V);
void nm_neldermead(NmDriver* D, const NmConfig* nc, NmRng* rng,
                          double* xbest, double* fbest, double* penbest);
void nm_randomsearch(NmDriver* D, const NmConfig* nc, NmRng* rng,
                            double* xbest, double* fbest, double* penbest);
void nm_sa(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest);
void nm_dual_annealing(NmDriver* D, const NmConfig* nc, NmRng* rng,
                              double* xbest, double* fbest, double* penbest);
void nm_shgo(NmDriver* D, const NmConfig* nc, NmRng* rng,
                    double* xbest, double* fbest, double* penbest);
void nm_direct(NmDriver* D, const NmConfig* nc, NmRng* rng,
                      double* xbest, double* fbest, double* penbest);
bool nm_apply_option(Expr* rule, FmOpts* opts, NmConfig* nc, const char* fn);
void nm_varset_free(NmVarSet* vs);
bool nm_parse_vars(Expr* var_arg, NmVarSet* vs, const char* fn);
Expr* nm_filter_int(Expr* cons, Expr** vars, size_t n, bool* is_int);
Expr* nm_build_result(double fmin, Expr** vars, const double* vals,
                             const bool* is_int, size_t n);
Expr* nm_build_infeasible(Expr** vars, size_t n);
Expr* nm_minimize_driver(Expr* res, const char* fn_name);
Expr* builtin_nminimize(Expr* res);
Expr* builtin_nmaximize(Expr* res);

#endif /* FINDMIN_INTERNAL_H */
