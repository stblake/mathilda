/* Mathilda — NDSolve shared infrastructure (see ndsolve_common.h). */
#include "ndsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../arithmetic.h"
#include "../common.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include "ndsolve_compile.h"    /* nonlinear RHS bytecode compiler */
#ifdef USE_LAPACK
#include "../linalg/lapack.h"   /* dgbtrf_/dgbtrs_/dgetrf_/dgetrs_ */
#endif

/* ------------------------------------------------------------------ *
 *  Block-style variable binding                                       *
 * ------------------------------------------------------------------ */
void nd_bind_snapshot(NdBind* b, const char* name) {
    SymbolDef* def = symtab_get_def(name);
    b->name = name;
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void nd_bind_clear_temp(SymbolDef* def) {
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
}

void nd_bind_set(NdBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    nd_bind_clear_temp(def);
    Expr* sym = expr_new_symbol(b->name);
    /* symtab_add_own_value copies both pattern and replacement, so we own and
     * must free `sym` and `value` afterwards. */
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
    expr_free(value);
}

void nd_bind_restore(NdBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    nd_bind_clear_temp(def);
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Numeric leaf coercion                                              *
 * ------------------------------------------------------------------ */
bool nd_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den; expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        double dd = mpz_get_d(den);
        if (dd == 0.0) { mpz_clears(num, den, NULL); return false; }
        *out = mpz_get_d(num) / dd;
        mpz_clears(num, den, NULL);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Small expr helpers shared by the ODE and MOL/PDE front-ends         *
 * ------------------------------------------------------------------ */
Expr* nd_call1(const char* head, Expr* a) {
    Expr* args[1] = { a }; return expr_new_function(expr_new_symbol(head), args, 1);
}
Expr* nd_call2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b }; return expr_new_function(expr_new_symbol(head), args, 2);
}

bool nd_eval_to_double(Expr* e, NumericSpec spec, double* out) {
    Expr* v = eval_and_free(expr_copy(e));
    if (!v) return false;
    Expr* nv = numericalize(v, spec);
    expr_free(v);
    bool ok = nv && nd_to_double(nv, out) && isfinite(*out);
    expr_free(nv);
    return ok;
}

Expr* nd_replace_all(Expr* body, Expr** lits, Expr** subs, size_t n) {
    Expr** rules = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* ra[2] = { expr_copy(lits[i]), expr_copy(subs[i]) };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* call = nd_call2(SYM_ReplaceAll, body, rl);
    return eval_and_free(call);
}

/* ------------------------------------------------------------------ *
 *  RHS + Jacobian evaluation                                          *
 * ------------------------------------------------------------------ */

/* Bind t and every reduced-state symbol to the current numeric point. */
static void nd_bind_point(NdProblem* P, double t, const double* Y) {
    nd_bind_set(&P->bind_t, expr_new_real(t));
    for (size_t i = 0; i < P->d; i++)
        nd_bind_set(&P->bind_y[i], expr_new_real(Y[i]));
    eval_clock_bump();
}

/* Evaluate a bound expression to a finite double. */
static bool nd_eval_bound(NdProblem* P, Expr* e, double* out) {
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(e));
    arith_warnings_mute_pop();
    if (!raw) return false;
    Expr* num = numericalize(raw, P->spec);
    expr_free(raw);
    bool ok = num && nd_to_double(num, out) && isfinite(*out);
    expr_free(num);
    return ok;
}

bool nd_rhs_real(NdProblem* P, double t, const double* Y, double* out) {
    /* Compiled linear fast path: out = A*Y + s(t), pure arithmetic. */
    if (P->op) {
        NdOperator* op = P->op;
        size_t n = op->n;
        nd_operator_matvec(op, Y, out);       /* out = A·Y (banded/dense BLAS) */
        if (op->time_forcing) {
            nd_bind_set(&P->bind_t, expr_new_real(t));
            for (size_t i = 0; i < n; i++) {
                double s;
                if (!nd_eval_bound(P, op->st[i], &s)) return false;
                out[i] += s;
            }
        } else if (op->s0) {
            for (size_t i = 0; i < n; i++) out[i] += op->s0[i];
        }
        if (P->eval_monitor) {
            nd_bind_point(P, t, Y);
            Expr* m = eval_and_free(expr_copy(P->eval_monitor)); expr_free(m);
        }
        return true;
    }
    /* Compiled nonlinear fast path: run bytecode over the state, no evaluator.
     * Compile lazily on first use; fall back to the symbolic sampler if the RHS
     * uses a construct the compiler does not support, or an EvaluationMonitor is
     * attached (which must run through the evaluator). */
    if (!P->compile_failed && !P->compiled && !P->eval_monitor && P->f) {
        P->compiled = nd_compile_rhs(P);
        if (!P->compiled) P->compile_failed = true;
    }
    if (P->compiled)
        return nd_compiled_eval(P->compiled, t, Y, out);

    nd_bind_point(P, t, Y);
    if (P->eval_monitor) { Expr* m = eval_and_free(expr_copy(P->eval_monitor)); expr_free(m); }
    for (size_t i = 0; i < P->d; i++)
        if (!nd_eval_bound(P, P->f[i], &out[i])) return false;
    return true;
}

/* Build the symbolic Jacobian jac[i][j] = D[f_i, y_j] once (entries may be
 * NULL if D failed).  Returns the d*d table, or NULL on allocation failure. */
static void nd_build_jacobian(NdProblem* P) {
    if (P->jac_built) return;
    size_t d = P->d;
    P->jac = malloc(sizeof(Expr**) * d);
    for (size_t i = 0; i < d; i++) {
        P->jac[i] = malloc(sizeof(Expr*) * d);
        for (size_t j = 0; j < d; j++) {
            Expr* args[2] = { expr_copy(P->f[i]), expr_copy(P->ysym[j]) };
            Expr* call = expr_new_function(expr_new_symbol(SYM_D), args, 2);
            P->jac[i][j] = eval_and_free(call);
        }
    }
    P->jac_built = true;
}

bool nd_jacobian_real(NdProblem* P, double t, const double* Y, double* Jout) {
    size_t d = P->d;
    /* Compiled linear fast path: the Jacobian is exactly the constant A. */
    if (P->op) { memcpy(Jout, P->op->A, sizeof(double) * d * d); return true; }
    /* Compiled nonlinear fast path: colored finite differences over the bytecode
     * RHS (O(bandwidth) evaluations, no evaluator). */
    if (!P->compile_failed && !P->compiled && !P->eval_monitor && P->f) {
        P->compiled = nd_compile_rhs(P);
        if (!P->compiled) P->compile_failed = true;
    }
    if (P->compiled)
        return nd_compiled_jacobian(P->compiled, t, Y, Jout);
    nd_build_jacobian(P);
    /* Try the symbolic Jacobian first. */
    bool ok_all = (P->jac != NULL);
    if (ok_all) {
        nd_bind_point(P, t, Y);
        for (size_t i = 0; i < d && ok_all; i++)
            for (size_t j = 0; j < d && ok_all; j++) {
                if (!P->jac[i][j] || !nd_eval_bound(P, P->jac[i][j], &Jout[i*d + j]))
                    ok_all = false;
            }
    }
    if (ok_all) return true;

    /* Central finite-difference fallback, column by column. */
    double* Yp = malloc(sizeof(double) * d);
    double* fp = malloc(sizeof(double) * d);
    double* fm = malloc(sizeof(double) * d);
    bool ok = true;
    for (size_t j = 0; j < d && ok; j++) {
        double hj = (fabs(Y[j]) + 1.0) * 1.0e-7;
        memcpy(Yp, Y, sizeof(double) * d);
        Yp[j] = Y[j] + hj; if (!nd_rhs_real(P, t, Yp, fp)) { ok = false; break; }
        Yp[j] = Y[j] - hj; if (!nd_rhs_real(P, t, Yp, fm)) { ok = false; break; }
        for (size_t i = 0; i < d; i++) Jout[i*d + j] = (fp[i] - fm[i]) / (2.0 * hj);
    }
    free(Yp); free(fp); free(fm);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Options + tolerances                                               *
 * ------------------------------------------------------------------ */
void nd_opts_default(NdOpts* o) {
    o->method = NULL;
    o->method_subopts = NULL;
    o->spec = numeric_machine_spec();
    o->wp_bits = 53;
    o->acc_goal = -1.0;
    o->prec_goal = -1.0;
    o->max_steps = -1;
    o->max_step_size = 0.0;
    o->max_step_fraction = 0.1;
    o->starting_step = 0.0;
    o->interp_order = -1;
    o->step_monitor = NULL;
    o->eval_monitor = NULL;
    o->norm_function = NULL;
}

NdTol nd_resolve_tol(const NdOpts* o) {
    /* Automatic goals default to WorkingPrecision/2 digits. */
    double wp_digits = o->wp_bits > 0 ? numeric_bits_to_digits(o->wp_bits)
                                      : NUMERIC_MACHINE_PRECISION_DIGITS;
    double pg = o->prec_goal, ag = o->acc_goal;
    if (pg < 0.0) pg = wp_digits * 0.5;
    if (ag < 0.0) ag = wp_digits * 0.5;
    NdTol t;
    t.rtol = (pg >= HUGE_VAL) ? 0.0 : pow(10.0, -pg);
    t.atol = (ag >= HUGE_VAL) ? 0.0 : pow(10.0, -ag);
    /* Keep the scale strictly positive so the WRMS norm is well defined. */
    if (t.rtol <= 0.0 && t.atol <= 0.0) t.atol = 1e-300;
    return t;
}

/* ------------------------------------------------------------------ *
 *  Norm + dense solve                                                 *
 * ------------------------------------------------------------------ */
double nd_wrms_norm(size_t d, const double* e, const double* y,
                    const double* ynew, NdTol tol) {
    double sum = 0.0;
    for (size_t i = 0; i < d; i++) {
        double ay = fabs(y[i]);
        double an = ynew ? fabs(ynew[i]) : ay;
        double sc = tol.atol + tol.rtol * (ay > an ? ay : an);
        if (sc <= 0.0) sc = 1e-300;
        double r = e[i] / sc;
        sum += r * r;
    }
    return sqrt(sum / (double)d);
}

bool nd_dense_solve(size_t n, double* A, double* b) {
    for (size_t col = 0; col < n; col++) {
        size_t piv = col; double maxv = fabs(A[col*n + col]);
        for (size_t r = col + 1; r < n; r++) {
            double v = fabs(A[r*n + col]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        if (maxv == 0.0) return false;
        if (piv != col) {
            for (size_t j = col; j < n; j++) {
                double t = A[col*n + j]; A[col*n + j] = A[piv*n + j]; A[piv*n + j] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        for (size_t r = col + 1; r < n; r++) {
            double m = A[r*n + col] / A[col*n + col];
            for (size_t j = col; j < n; j++) A[r*n + j] -= m * A[col*n + j];
            b[r] -= m * b[col];
        }
    }
    for (size_t i = n; i-- > 0; ) {
        double s = b[i];
        for (size_t j = i + 1; j < n; j++) s -= A[i*n + j] * b[j];
        b[i] = s / A[i*n + i];
    }
    return true;
}

/* Banded no-pivot LU solve on a dense-stored matrix M (nonzeros within
 * bandwidth kl/ku).  Elimination touches only the band, so it costs
 * O(n·kl·ku).  Suitable for the diagonally-dominant iteration matrix
 * I - h·theta·A of a discretized elliptic operator.  Returns false on a small
 * pivot so the caller can retry with the dense solve. */
bool nd_banded_solve(size_t n, int kl, int ku, double* M, double* b) {
    for (size_t k = 0; k < n; k++) {
        double piv = M[k*n + k];
        if (fabs(piv) < 1e-300) return false;
        size_t ilast = (k + (size_t)kl < n - 1) ? k + (size_t)kl : n - 1;
        size_t jlast = (k + (size_t)ku < n - 1) ? k + (size_t)ku : n - 1;
        for (size_t i = k + 1; i <= ilast; i++) {
            double m = M[i*n + k] / piv;
            if (m == 0.0) continue;
            for (size_t j = k + 1; j <= jlast; j++) M[i*n + j] -= m * M[k*n + j];
            b[i] -= m * b[k];
        }
    }
    for (size_t i = n; i-- > 0; ) {
        double s = b[i];
        size_t jlast = (i + (size_t)ku < n - 1) ? i + (size_t)ku : n - 1;
        for (size_t j = i + 1; j <= jlast; j++) s -= M[i*n + j] * b[j];
        b[i] = s / M[i*n + i];
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Iteration-matrix factor + solve  (I - coef*J) x = b                *
 *                                                                     *
 *  Factor ONCE, then back-substitute per Newton iteration.  For the   *
 *  compiled operator the Jacobian J = A is constant so the factor is   *
 *  reused across every iteration of a solve; for a symbolic Jacobian   *
 *  the factor is rebuilt each iteration (J changes).  A banded pivoted *
 *  LAPACK factor is used when the matrix is narrow (the FD structure   *
 *  of a discretized operator), a dense pivoted LAPACK factor for wide  *
 *  bands, and the scalar band/dense LU when built without LAPACK.      *
 * ------------------------------------------------------------------ */
struct NdIterFactor {
    size_t d;
    int    mode;        /* 0 banded-LAPACK, 1 dense-LAPACK, 2 scalar fallback   */
    int    kl, ku, ldab;
    double* ab;         /* factor band (mode 0) or dense col-major LU (mode 1)  */
    int*    ipiv;
    double* M;          /* scalar fallback: I - coef*J, row-major (mode 2)      */
    bool    fb_banded;  /* scalar fallback uses the banded solver               */
    double  coef;
};

/* Detect the half-bandwidths of a dense d x d row-major matrix A. */
static void nd_detect_band(const double* A, size_t d, int* kl, int* ku) {
    int lo = 0, hi = 0;
    for (size_t i = 0; i < d; i++)
        for (size_t j = 0; j < d; j++)
            if (A[i * d + j] != 0.0) {
                if ((long)i - (long)j > lo) lo = (int)((long)i - (long)j);
                if ((long)j - (long)i > hi) hi = (int)((long)j - (long)i);
            }
    *kl = lo; *ku = hi;
}

NdIterFactor* nd_iter_factor(const NdOperator* op, size_t d, double coef,
                             const double* Jdense) {
    /* how to read A(i,j): the operator's packed band, else a dense row-major J */
    const double* AB = (op && op->banded && op->AB) ? op->AB : NULL;
    int mv_ld = AB ? (op->kl + op->ku + 1) : 0;
    const double* Adense = AB ? NULL : (op ? op->A : Jdense);
    if (!AB && !Adense) return NULL;

    int kl, ku;
    if (op) { kl = op->kl; ku = op->ku; }
    else    { nd_detect_band(Adense, d, &kl, &ku); }
    bool banded = ((size_t)(kl + ku + 1) <= d / 2 + 2);

    /* read M(i,j) = (i==j) - coef*A(i,j) */
    #define A_IJ(i, j) (AB ? AB[(size_t)(ku + (long)(i) - (long)(j)) + (size_t)(j) * (size_t)mv_ld] \
                            : Adense[(size_t)(i) * d + (size_t)(j)])

    NdIterFactor* F = calloc(1, sizeof(*F));
    if (!F) return NULL;
    F->d = d; F->kl = kl; F->ku = ku; F->coef = coef;

#ifdef USE_LAPACK
    if (banded) {
        int n = (int)d, ldab = 2 * kl + ku + 1;
        double* ab = calloc((size_t)ldab * d, sizeof(double));
        int* ipiv = malloc(sizeof(int) * d);
        if (!ab || !ipiv) { free(ab); free(ipiv); free(F); return NULL; }
        for (size_t j = 0; j < d; j++) {
            size_t i0 = (j > (size_t)ku) ? j - (size_t)ku : 0;
            size_t i1 = (j + (size_t)kl < d - 1) ? j + (size_t)kl : d - 1;
            for (size_t i = i0; i <= i1; i++)
                ab[(size_t)(kl + ku + (long)i - (long)j) + j * (size_t)ldab] =
                    (i == j ? 1.0 : 0.0) - coef * A_IJ(i, j);
        }
        int info = 0;
        dgbtrf_(&n, &n, &kl, &ku, ab, &ldab, ipiv, &info);
        if (info != 0) { free(ab); free(ipiv); free(F); return NULL; }
        F->mode = 0; F->ab = ab; F->ipiv = ipiv; F->ldab = ldab;
        return F;
    }
    {   /* dense pivoted LU (wide band / full coupling) */
        int n = (int)d;
        double* ab = malloc(sizeof(double) * d * d);   /* col-major */
        int* ipiv = malloc(sizeof(int) * d);
        if (!ab || !ipiv) { free(ab); free(ipiv); free(F); return NULL; }
        for (size_t j = 0; j < d; j++)
            for (size_t i = 0; i < d; i++)
                ab[i + j * d] = (i == j ? 1.0 : 0.0) - coef * A_IJ(i, j);
        int info = 0;
        dgetrf_(&n, &n, ab, &n, ipiv, &info);
        if (info != 0) { free(ab); free(ipiv); free(F); return NULL; }
        F->mode = 1; F->ab = ab; F->ipiv = ipiv;
        return F;
    }
#else
    {   /* no LAPACK: keep M row-major, re-run the hand LU per solve */
        double* M = malloc(sizeof(double) * d * d);
        if (!M) { free(F); return NULL; }
        for (size_t i = 0; i < d; i++)
            for (size_t j = 0; j < d; j++)
                M[i * d + j] = (i == j ? 1.0 : 0.0) - coef * A_IJ(i, j);
        F->mode = 2; F->M = M; F->fb_banded = banded;
        return F;
    }
#endif
    #undef A_IJ
}

bool nd_iter_solve(NdIterFactor* F, double* b) {
    if (!F) return false;
    int n = (int)F->d, nrhs = 1, info = 0;
#ifdef USE_LAPACK
    if (F->mode == 0) {
        dgbtrs_("N", &n, &F->kl, &F->ku, &nrhs, F->ab, &F->ldab, F->ipiv, b, &n, &info);
        return info == 0;
    }
    if (F->mode == 1) {
        dgetrs_("N", &n, &nrhs, F->ab, &n, F->ipiv, b, &n, &info);
        return info == 0;
    }
#endif
    {   /* scalar fallback: solve a fresh copy so F->M survives for reuse */
        double* M = malloc(sizeof(double) * F->d * F->d);
        if (!M) return false;
        memcpy(M, F->M, sizeof(double) * F->d * F->d);
        bool ok = F->fb_banded ? nd_banded_solve(F->d, F->kl, F->ku, M, b) : false;
        if (!ok) { memcpy(M, F->M, sizeof(double) * F->d * F->d); ok = nd_dense_solve(F->d, M, b); }
        free(M);
        return ok;
    }
}

void nd_iter_factor_free(NdIterFactor* F) {
    if (!F) return;
    free(F->ab); free(F->ipiv); free(F->M); free(F);
}

/* ------------------------------------------------------------------ *
 *  Implicit theta-method Newton solve                                 *
 * ------------------------------------------------------------------ */
bool nd_newton_theta(NdProblem* P, double t1, const double* Ybase,
                     double h, double theta, const double* rhs_const,
                     const double* Zguess, double* Ynew, NdTol tol) {
    size_t d = P->d;
    double coef = h * theta;
    bool op_const = (P->op != NULL);   /* constant Jacobian -> factor once */
    double* Z  = malloc(sizeof(double) * d);
    double* f  = malloc(sizeof(double) * d);
    double* G  = malloc(sizeof(double) * d);
    double* J  = op_const ? NULL : malloc(sizeof(double) * d * d);
    double* dZ = malloc(sizeof(double) * d);
    NdIterFactor* F = op_const ? nd_iter_factor(P->op, d, coef, NULL) : NULL;
    if (op_const && !F) { free(Z); free(f); free(G); free(dZ); return false; }
    memcpy(Z, Zguess ? Zguess : Ybase, sizeof(double) * d);
    bool converged = false;
    for (int it = 0; it < 12; it++) {
        if (!nd_rhs_real(P, t1, Z, f)) break;
        for (size_t i = 0; i < d; i++)
            G[i] = Z[i] - Ybase[i] - coef * f[i] - (rhs_const ? rhs_const[i] : 0.0);
        if (!op_const) {
            if (!nd_jacobian_real(P, t1, Z, J)) break;
            nd_iter_factor_free(F);
            F = nd_iter_factor(NULL, d, coef, J);   /* J changes each iteration */
            if (!F) break;
        }
        memcpy(dZ, G, sizeof(double) * d);
        if (!nd_iter_solve(F, dZ)) break;
        for (size_t i = 0; i < d; i++) Z[i] -= dZ[i];
        double nrm = nd_wrms_norm(d, dZ, Z, NULL, tol);
        if (nrm <= 1e-2) { converged = true; break; }
    }
    if (converged) memcpy(Ynew, Z, sizeof(double) * d);
    nd_iter_factor_free(F);
    free(Z); free(f); free(G); free(J); free(dZ);
    return converged;
}

/* ------------------------------------------------------------------ *
 *  Stepper registry                                                   *
 * ------------------------------------------------------------------ */
static const NdStepper* const nd_steppers[] = {
    &nd_stepper_explicit_euler,
    &nd_stepper_explicit_midpoint,
    &nd_stepper_rk4,
    &nd_stepper_dopri5,
    &nd_stepper_backward_euler,
    &nd_stepper_implicit_trapezoid,
    &nd_stepper_bdf,
    &nd_stepper_adams,
};
static const size_t nd_nsteppers = sizeof(nd_steppers) / sizeof(nd_steppers[0]);

size_t nd_stepper_count(void) { return nd_nsteppers; }
const NdStepper* nd_stepper_at(size_t i) { return i < nd_nsteppers ? nd_steppers[i] : NULL; }
const NdStepper* nd_default_stepper(void) { return &nd_stepper_dopri5; }

const NdStepper* nd_lookup_stepper(const char* name) {
    if (!name) return nd_default_stepper();
    if (strcmp(name, "Automatic") == 0) return nd_default_stepper();
    /* "ExplicitRungeKutta" is the adaptive DOPRI5 workhorse. */
    if (strcmp(name, "ExplicitRungeKutta") == 0) return &nd_stepper_dopri5;
    if (strcmp(name, "RungeKutta") == 0) return &nd_stepper_rk4;
    if (strcmp(name, "StiffnessSwitching") == 0) return &nd_stepper_bdf;
    for (size_t i = 0; i < nd_nsteppers; i++)
        if (strcmp(name, nd_steppers[i]->name) == 0) return nd_steppers[i];
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Solution accumulation                                              *
 * ------------------------------------------------------------------ */
void nd_solution_init(NdSolution* s, size_t d) {
    s->n = 0; s->cap = 0; s->d = d; s->ts = NULL; s->Ys = NULL; s->dYs = NULL;
}
void nd_solution_free(NdSolution* s) {
    free(s->ts); free(s->Ys); free(s->dYs);
    s->ts = NULL; s->Ys = NULL; s->dYs = NULL; s->n = s->cap = 0;
}
void nd_solution_push(NdSolution* s, double t, const double* Y, const double* dY) {
    if (s->n == s->cap) {
        size_t nc = s->cap ? s->cap * 2 : 64;
        s->ts  = realloc(s->ts,  sizeof(double) * nc);
        s->Ys  = realloc(s->Ys,  sizeof(double) * nc * s->d);
        s->dYs = realloc(s->dYs, sizeof(double) * nc * s->d);
        s->cap = nc;
    }
    s->ts[s->n] = t;
    memcpy(&s->Ys[s->n * s->d],  Y,  sizeof(double) * s->d);
    memcpy(&s->dYs[s->n * s->d], dY, sizeof(double) * s->d);
    s->n++;
}

/* ------------------------------------------------------------------ *
 *  Adaptive driver                                                    *
 * ------------------------------------------------------------------ */

/* Hairer's starting-step-size heuristic (HNW I.II.4).  Exported so the adaptive
 * multistep drivers (BDF/Adams) can share a well-scaled first step. */
double nd_initial_step(NdProblem* P, const NdOpts* o, NdTol tol,
                       double t0, const double* Y0, const double* f0,
                       int order, double dir) {
    size_t d = P->d;
    double d0 = nd_wrms_norm(d, Y0, Y0, NULL, tol);
    double d1 = nd_wrms_norm(d, f0, Y0, NULL, tol);
    double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * (d0 / d1);
    double* Y1 = calloc(d, sizeof(double));   /* calloc: silences -O3 maybe-uninit */
    double* f1 = malloc(sizeof(double) * d);
    double* df = malloc(sizeof(double) * d);
    double h;
    for (size_t i = 0; i < d; i++) Y1[i] = Y0[i] + dir * h0 * f0[i];
    if (nd_rhs_real(P, t0 + dir * h0, Y1, f1)) {
        for (size_t i = 0; i < d; i++) df[i] = f1[i] - f0[i];
        double d2 = nd_wrms_norm(d, df, Y0, NULL, tol) / h0;
        double dm = (d1 > d2 ? d1 : d2);
        double h1 = (dm <= 1e-15) ? (h0 * 1e-3 > 1e-6 ? h0 * 1e-3 : 1e-6)
                                  : pow(0.01 / dm, 1.0 / (order + 1));
        h = (100.0 * h0 < h1 ? 100.0 * h0 : h1);
    } else {
        h = h0;
    }
    free(Y1); free(f1); free(df);
    /* clamp */
    double span = fabs(P->tmax - P->tmin);
    double cap = o->max_step_fraction * span;
    if (cap > 0.0 && h > cap) h = cap;
    if (o->max_step_size > 0.0 && h > o->max_step_size) h = o->max_step_size;
    if (h <= 0.0) h = 1e-6;
    return dir * h;
}

/* Take one trial step, returning the WRMS error estimate (0 for a pure fixed
 * step) and the proposed Ynew.  `adaptive` requests step-doubling for a
 * non-embedded stepper.  Returns false on a failed sample. */
static bool nd_take_step(const NdStepper* S, NdProblem* P, NdTol tol,
                         double t, const double* Y, double h, bool adaptive,
                         double* Ynew, double* K, double* err_out, int* q_out) {
    size_t d = P->d;
    if (S->flags & ND_ADAPTIVE) {
        double* Yerr = malloc(sizeof(double) * d);
        bool ok = S->single_step(S, P, t, Y, h, Ynew, Yerr, K);
        if (ok) {
            *err_out = nd_wrms_norm(d, Yerr, Y, Ynew, tol);
            int p = S->order, pe = S->err_order;
            *q_out = (p < pe ? p : pe) + 1;
        }
        free(Yerr);
        return ok;
    }
    if (!adaptive) {
        bool ok = S->single_step(S, P, t, Y, h, Ynew, NULL, K);
        *err_out = 0.0; *q_out = S->order + 1;
        return ok;
    }
    /* Step doubling: one step h vs two steps h/2. */
    double* Yb = malloc(sizeof(double) * d);
    double* Ym = malloc(sizeof(double) * d);
    double* e  = malloc(sizeof(double) * d);
    bool ok = S->single_step(S, P, t, Y, h, Yb, NULL, K);
    if (ok) ok = S->single_step(S, P, t, Y, 0.5 * h, Ym, NULL, K);
    if (ok) {
        double* Ym2 = malloc(sizeof(double) * d);
        ok = S->single_step(S, P, t + 0.5 * h, Ym, 0.5 * h, Ym2, NULL, K);
        if (ok) {
            double denom = ldexp(1.0, S->order) - 1.0;   /* 2^p - 1 */
            for (size_t i = 0; i < d; i++) {
                e[i] = (Ym2[i] - Yb[i]) / denom;
                Ynew[i] = Ym2[i] + e[i];                 /* local extrapolation */
            }
            *err_out = nd_wrms_norm(d, e, Y, Ynew, tol);
            *q_out = S->order + 1;
        }
        free(Ym2);
    }
    free(Yb); free(Ym); free(e);
    return ok;
}

/* Integrate in one direction from (start_t, start_Y) toward `target`.
 * The starting node is assumed already recorded by the caller. */
static NdStatus nd_integrate_dir(NdProblem* P, const NdStepper* S, const NdOpts* o,
                                 NdSolution* sol, NdTol tol,
                                 double start_t, const double* start_Y,
                                 const double* start_f, double target,
                                 int64_t max_steps) {
    size_t d = P->d;
    double dir = (target > start_t) ? 1.0 : -1.0;
    double span = fabs(P->tmax - P->tmin);
    if (fabs(target - start_t) <= 16.0 * DBL_EPSILON * (span + 1.0)) return ND_OK;

    bool adaptive = true;   /* honor accuracy goals even for fixed steppers */
    double h;
    if (o->starting_step > 0.0) { h = dir * o->starting_step; }
    else h = nd_initial_step(P, o, tol, start_t, start_Y, start_f, S->order, dir);

    double t = start_t;
    double* Y   = malloc(sizeof(double) * d);
    double* Yn  = malloc(sizeof(double) * d);
    double* fn  = malloc(sizeof(double) * d);
    double* K   = (S->stages > 0) ? malloc(sizeof(double) * S->stages * d) : NULL;
    memcpy(Y, start_Y, sizeof(double) * d);

    /* FSAL setup: the reusable first stage is seeded with f(start_t, start_Y),
     * eliminating one RHS evaluation per step (stage 1) and one per accepted
     * step (the node slope) for DOPRI5 — 6 evals/accepted step, matching ode45. */
    bool use_fsal = (S->flags & ND_FSAL) != 0;
    P->fsal = use_fsal;
    if (use_fsal) {
        P->fsal_cur = malloc(sizeof(double) * d);
        P->fsal_pending = malloc(sizeof(double) * d);
        memcpy(P->fsal_cur, start_f, sizeof(double) * d);
    }

    NdStatus status = ND_OK;
    int64_t steps = 0;
    double h_cap = (o->max_step_size > 0.0) ? o->max_step_size : HUGE_VAL;
    double frac_cap = (o->max_step_fraction > 0.0) ? o->max_step_fraction * span : HUGE_VAL;

    while ((target - t) * dir > 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
        /* clamp step magnitude */
        double hmag = fabs(h);
        if (hmag > h_cap) hmag = h_cap;
        if (hmag > frac_cap) hmag = frac_cap;
        h = dir * hmag;
        /* land exactly on target */
        if ((t + h - target) * dir > 0.0) h = target - t;
        if (fabs(h) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) { status = ND_ERR_STEPSIZE; break; }

        double err; int q;
        bool ok = nd_take_step(S, P, tol, t, Y, h, adaptive, Yn, K, &err, &q);
        steps++;
        if (steps > max_steps) { status = ND_ERR_MAXSTEPS; break; }
        if (!ok) {
            /* Failed step (implicit: Newton diverged; explicit: bad sample).
             * Shrink and retry — do NOT latch an error, since a smaller step
             * usually recovers.  Only a step-size collapse is terminal, and for
             * an implicit method that collapse means "Newton won't converge even
             * at a tiny step", reported as ND_ERR_NONCONV. */
            h *= 0.5;
            if (fabs(h) < 16.0 * DBL_EPSILON * (fabs(t) + 1.0)) {
                status = (S->flags & ND_IMPLICIT) ? ND_ERR_NONCONV : ND_ERR_STEPSIZE;
                break;
            }
            continue;
        }
        if (err <= 1.0) {
            /* accept */
            t += h;
            memcpy(Y, Yn, sizeof(double) * d);
            if (use_fsal) {
                /* the stepper's last stage IS f(t+h, Ynew): reuse it as the node
                 * slope and as the next step's first stage — no extra eval. */
                memcpy(fn, P->fsal_pending, sizeof(double) * d);
                memcpy(P->fsal_cur, P->fsal_pending, sizeof(double) * d);
            } else if (!nd_rhs_real(P, t, Y, fn)) { status = ND_ERR_SAMPLE; break; }
            nd_solution_push(sol, t, Y, fn);
            if (o->step_monitor) { Expr* m = eval_and_free(expr_copy(o->step_monitor)); expr_free(m); }
            double fac = (err > 0.0) ? 0.9 * pow(1.0 / err, 1.0 / q) : 5.0;
            if (fac < 0.2) fac = 0.2;
            if (fac > 5.0) fac = 5.0;
            h *= fac;
        } else {
            /* reject: shrink, never grow */
            double fac = 0.9 * pow(1.0 / err, 1.0 / q);
            if (fac < 0.2) fac = 0.2;
            if (fac > 1.0) fac = 1.0;
            h *= fac;
        }
    }
    free(Y); free(Yn); free(fn); free(K);
    if (use_fsal) { free(P->fsal_cur); free(P->fsal_pending);
                    P->fsal_cur = P->fsal_pending = NULL; P->fsal = false; }
    return status;
}

/* Sort accumulated nodes by ascending t and drop near-duplicates. */
static void nd_solution_sort(NdSolution* s) {
    size_t d = s->d;
    /* simple insertion sort (node lists are already piecewise monotone) */
    for (size_t i = 1; i < s->n; i++) {
        double t = s->ts[i];
        double* yr = malloc(sizeof(double) * d);
        double* dr = malloc(sizeof(double) * d);
        memcpy(yr, &s->Ys[i*d],  sizeof(double) * d);
        memcpy(dr, &s->dYs[i*d], sizeof(double) * d);
        size_t j = i;
        while (j > 0 && s->ts[j-1] > t) {
            s->ts[j] = s->ts[j-1];
            memcpy(&s->Ys[j*d],  &s->Ys[(j-1)*d],  sizeof(double) * d);
            memcpy(&s->dYs[j*d], &s->dYs[(j-1)*d], sizeof(double) * d);
            j--;
        }
        s->ts[j] = t;
        memcpy(&s->Ys[j*d],  yr, sizeof(double) * d);
        memcpy(&s->dYs[j*d], dr, sizeof(double) * d);
        free(yr); free(dr);
    }
    /* drop duplicates in t */
    size_t w = 0;
    for (size_t i = 0; i < s->n; i++) {
        if (w > 0 && fabs(s->ts[i] - s->ts[w-1]) <= 1e-14 * (fabs(s->ts[i]) + 1.0)) continue;
        if (i != w) {
            s->ts[w] = s->ts[i];
            memcpy(&s->Ys[w*d],  &s->Ys[i*d],  sizeof(double) * d);
            memcpy(&s->dYs[w*d], &s->dYs[i*d], sizeof(double) * d);
        }
        w++;
    }
    s->n = w;
}

/* Fixed step for the multistep methods.  Unlike the adaptive one-step driver,
 * these march at a constant step, so the step must both (a) cover the range
 * within the MaxSteps budget and (b) resolve the solution.  We aim for a target
 * number of nodes scaled by the precision goal, clamped by MaxStepSize /
 * MaxStepFraction.  (Hairer's adaptive-first-step heuristic is unsuitable here:
 * when f(t0)=0 it collapses to a tiny step that exhausts the budget.) */
double nd_fixed_step(NdProblem* P, const NdOpts* o, NdTol tol, double dir) {
    double span = fabs(P->tmax - P->tmin);
    if (o->starting_step > 0.0) return dir * o->starting_step;
    if (span <= 0.0) return dir * 1e-3;
    /* target node count: more nodes for tighter tolerances (order-2 methods),
     * but never more than fits the MaxSteps budget (else a fixed march runs out
     * of steps before reaching the endpoint and the result extrapolates). */
    double digits = (tol.rtol > 0.0) ? -log10(tol.rtol) : 8.0;
    double target = 500.0 * pow(10.0, digits / 4.0);   /* grows with accuracy */
    double budget = (o->max_steps > 0) ? (double)o->max_steps : 10000.0;
    double budget_cap = 0.9 * budget;                  /* leave headroom       */
    if (target < 500.0) target = 500.0;
    if (target > budget_cap) target = budget_cap;
    double h = span / target;
    if (o->max_step_size > 0.0 && h > o->max_step_size) h = o->max_step_size;
    double frac = o->max_step_fraction * span;
    if (frac > 0.0 && h > frac) h = frac;
    return dir * h;
}

NdStatus nd_integrate(NdProblem* P, const NdStepper* S, const NdOpts* o, NdSolution* sol) {
    size_t d = P->d;
    NdTol tol = nd_resolve_tol(o);
    P->tol.rtol = tol.rtol; P->tol.atol = tol.atol;
    /* Multistep methods run their own history-managed loop. */
    if (S->flags & ND_MULTISTEP) {
        NdStatus ms = (S == &nd_stepper_adams) ? nd_multistep_adams(P, o, sol)
                                               : nd_multistep_bdf(P, o, sol);
        nd_solution_sort(sol);
        return ms;
    }
    int64_t max_steps = (o->max_steps > 0) ? o->max_steps : 10000;

    double* Y0 = malloc(sizeof(double) * d);
    double* f0 = malloc(sizeof(double) * d);
    memcpy(Y0, P->Y0, sizeof(double) * d);
    if (!nd_rhs_real(P, P->t0, Y0, f0)) { free(Y0); free(f0); return ND_ERR_SAMPLE; }
    nd_solution_push(sol, P->t0, Y0, f0);

    NdStatus st = ND_OK;
    /* forward toward tmax */
    if (P->tmax > P->t0)
        st = nd_integrate_dir(P, S, o, sol, tol, P->t0, Y0, f0, P->tmax, max_steps);
    /* backward toward tmin */
    if (P->tmin < P->t0) {
        NdStatus st2 = nd_integrate_dir(P, S, o, sol, tol, P->t0, Y0, f0, P->tmin, max_steps);
        if (st == ND_OK) st = st2;
    }
    /* also cover the case t0 == tmin (integrate forward only, already done) or
     * t0 == tmax (backward only). */
    if (P->tmax <= P->t0 && P->tmin < P->t0) { /* handled by backward branch */ }

    free(Y0); free(f0);
    nd_solution_sort(sol);
    return st;
}

/* ------------------------------------------------------------------ *
 *  Result assembly (Interpolation -> InterpolatingFunction)           *
 * ------------------------------------------------------------------ */

/* Build the InterpolatingFunction for one scalar component (real path):
 * Hermite triples {{t_i}, y_i, y'_i} passed to the Interpolation builtin. */
static Expr* nd_build_component(const NdSolution* sol, size_t comp) {
    size_t n = sol->n, d = sol->d;
    if (n < 2) return NULL;
    Expr** entries = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* coord_el = expr_new_real(sol->ts[i]);
        Expr* coord = expr_new_function(expr_new_symbol(SYM_List), &coord_el, 1);
        Expr* trip[3] = { coord, expr_new_real(sol->Ys[i*d + comp]),
                          expr_new_real(sol->dYs[i*d + comp]) };
        entries[i] = expr_new_function(expr_new_symbol(SYM_List), trip, 3);
    }
    Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, n);
    free(entries);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1);
    Expr* ifun = eval_and_free(call);
    if (!head_is(ifun, SYM_InterpolatingFunction)) { expr_free(ifun); return NULL; }
    return ifun;
}

Expr* nd_build_result(NdProblem* P, const NdOpts* o, const NdSolution* sol) {
    (void)o;
    if (sol->n < 2) return NULL;
    Expr** rules = malloc(sizeof(Expr*) * P->nfun);
    size_t nr = 0;
    for (size_t k = 0; k < P->nfun; k++) {
        size_t comp = P->fun_state0[k];
        Expr* ifun = nd_build_component(sol, comp);
        if (!ifun) { for (size_t q = 0; q < nr; q++) expr_free(rules[q]); free(rules); return NULL; }
        /* lhs: u  or  u[x] when fun_applied */
        Expr* lhs;
        if (P->fun_applied) {
            Expr* xarg = expr_new_symbol(P->tvar);
            lhs = expr_new_function(expr_new_symbol(P->fun_names[k]), &xarg, 1);
            /* rhs: IF[x] */
            Expr* xarg2 = expr_new_symbol(P->tvar);
            ifun = expr_new_function(ifun, &xarg2, 1);
        } else {
            lhs = expr_new_symbol(P->fun_names[k]);
        }
        Expr* rargs[2] = { lhs, ifun };
        rules[nr++] = expr_new_function(expr_new_symbol(SYM_Rule), rargs, 2);
    }
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), rules, nr);
    free(rules);
    Expr* outer = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
    return outer;
}

Expr* nd_build_result_complex(NdProblem* P, const NdOpts* o, const NdSolution* sol) {
    (void)o;
    if (sol->n < 2) return NULL;
    Expr** rules = malloc(sizeof(Expr*) * P->nfun);
    size_t nr = 0;
    for (size_t k = 0; k < P->nfun; k++) {
        size_t comp = P->fun_state0[k];
        Expr* ifRe = nd_build_component(sol, 2 * comp);
        Expr* ifIm = nd_build_component(sol, 2 * comp + 1);
        if (!ifRe || !ifIm) {
            expr_free(ifRe); expr_free(ifIm);
            for (size_t q = 0; q < nr; q++) expr_free(rules[q]);
            free(rules); return NULL;
        }
        /* body = ifRe[t] + I ifIm[t] */
        Expr* ta = expr_new_symbol(P->tvar);
        Expr* reApp = expr_new_function(ifRe, &ta, 1);
        Expr* tb = expr_new_symbol(P->tvar);
        Expr* imApp = expr_new_function(ifIm, &tb, 1);
        Expr* body = nd_call2("Plus", reApp,
                              nd_call2("Times", expr_new_symbol("I"), imApp));
        Expr* lhs, *rhs;
        if (P->fun_applied) {
            Expr* xa = expr_new_symbol(P->tvar);
            lhs = expr_new_function(expr_new_symbol(P->fun_names[k]), &xa, 1);
            rhs = body;
        } else {
            lhs = expr_new_symbol(P->fun_names[k]);
            Expr* param = expr_new_symbol(P->tvar);
            Expr* plist = expr_new_function(expr_new_symbol(SYM_List), &param, 1);
            Expr* fargs[2] = { plist, body };
            rhs = expr_new_function(expr_new_symbol("Function"), fargs, 2);
        }
        rules[nr++] = nd_call2(SYM_Rule, lhs, rhs);
    }
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), rules, nr);
    free(rules);
    return expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
}
