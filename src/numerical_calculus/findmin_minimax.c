/* findmin_minimax.c -- minimax objectives Max[p_1, ..., p_m] by their epigraph.
 *
 * A best-approximation problem is posed, naturally, as the minimisation of a
 * maximum: NMinimize[Max[Abs[f(x_k) - p(x_k; c)]], c].  The objective is
 * non-smooth exactly where its minimum lies -- at the optimum several pieces
 * are equal (equioscillation) -- so a quasi-Newton polish, whose model assumes
 * a gradient, stalls a few digits short: the SIAM problem 5 discrete minimax
 * (128 complex moduli over 8 coefficients) came back with 7.6 digits.
 *
 * The epigraph form is smooth:
 *
 *     min t   subject to   p_k(x) <= t       (k = 1..m)
 *
 * with an Abs piece split by kind:
 *   - Abs[e], e real for real x (ComplexExpand[Im[e]] == 0):  e - t <= 0 and
 *     -e - t <= 0 (the LP-like case of real Chebyshev approximation);
 *   - Abs[e], e complex:  Re[e]^2 + Im[e]^2 - t^2 <= 0 with t >= 0 (smooth,
 *     convex for e affine in x);
 *   - any other smooth piece p:  p - t <= 0.
 * It is solved by the SQP solver (fm_run_slsqp, every constraint and gradient
 * compiled), which converges super-linearly onto the active set.
 *
 * At WorkingPrecision the active set of the machine solution (the pieces
 * within 1e-8 of the maximum) defines a system of equations -- each active
 * piece equal to t (to t^2 for a complex modulus) -- that is square, n + 1
 * equations in (x, t), under the Haar condition of approximation theory, and
 * over-determined but consistent when more pieces tie.  Newton's method (or
 * Gauss-Newton) in MPFR solves it to the working precision in a handful of
 * steps.  The result is certified rather than trusted: no inactive piece may
 * exceed t at the refined point, and in the square case the KKT multipliers
 * (J^T lambda = -e_t, at machine precision) must be non-negative; otherwise
 * the caller keeps its fallback.  The variables are assumed real (as NMinimize
 * and FindMinimum assume). */
#include "findmin_internal.h"

#define FM_MM_MAX_PIECES 100000

void fm_minimax_free(FmMinimax* mm) {
    if (!mm) return;
    for (size_t k = 0; k < mm->np; k++) {
        expr_free(mm->piece[k]);
        if (mm->grad && mm->grad[k]) {
            for (size_t i = 0; i < mm->n; i++) expr_free(mm->grad[k][i]);
            free(mm->grad[k]);
        }
    }
    free(mm->piece); free(mm->kind); free(mm->grad);
    mm->piece = NULL; mm->kind = NULL; mm->grad = NULL; mm->np = 0;
}

static Expr* mm_fn1(const char* head, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(head), args, 1);
}

static Expr* mm_fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), args, 2);
}

/* Does e contain a non-smooth head (Abs, Max, Min, Sign) or a non-real one? */
static bool mm_nonsmooth(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h->type == EXPR_SYMBOL) {
            const char* s = h->data.symbol.name;
            if (s == SYM_Abs || s == SYM_Max || s == SYM_Min || s == SYM_Sign
                || s == SYM_Re || s == SYM_Im || s == SYM_Complex)
                return true;
        }
        if (mm_nonsmooth(h)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (mm_nonsmooth(e->data.function.args[i])) return true;
    }
    return false;
}

static bool mm_is_zero(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    return false;
}

static bool mm_add_piece(FmMinimax* mm, size_t* cap, Expr* p, int kind) {
    if (mm->np >= *cap) {
        size_t nc = *cap ? 2 * *cap : 16;
        Expr** np_ = (Expr**)realloc(mm->piece, sizeof(Expr*) * nc);
        if (!np_) return false;
        mm->piece = np_;
        int* nk = (int*)realloc(mm->kind, sizeof(int) * nc);
        if (!nk) return false;
        mm->kind = nk;
        *cap = nc;
    }
    mm->piece[mm->np] = p;
    mm->kind[mm->np] = kind;
    mm->np++;
    return true;
}

bool fm_minimax_detect(Expr* f, Expr** vars, size_t n, FmMinimax* mm) {
    memset(mm, 0, sizeof(*mm));
    mm->n = n;
    if (!f || f->type != EXPR_FUNCTION || n == 0) return false;
    const Expr* h = f->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Max) return false;
    size_t m = f->data.function.arg_count;
    if (m < 2 || m > FM_MM_MAX_PIECES) return false;
    size_t cap = 0;
    bool any_var = false;
    for (size_t k = 0; k < m; k++) {
        Expr* a = f->data.function.args[k];
        Expr* piece = NULL;
        int kind = FM_MM_SMOOTH;
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol.name == SYM_Abs
            && a->data.function.arg_count == 1) {
            Expr* e = a->data.function.args[0];
            Expr* re = eval_and_free(mm_fn1(SYM_ComplexExpand, mm_fn1(SYM_Re, expr_copy(e))));
            Expr* im = eval_and_free(mm_fn1(SYM_ComplexExpand, mm_fn1(SYM_Im, expr_copy(e))));
            if (!re || !im) { expr_free(re); expr_free(im); goto fail; }
            if (mm_is_zero(im)) {
                expr_free(im);
                piece = re; kind = FM_MM_REAL_ABS;
            } else {
                Expr* q = mm_fn2(SYM_Plus, mm_fn2(SYM_Power, re, expr_new_integer(2)),
                                           mm_fn2(SYM_Power, im, expr_new_integer(2)));
                piece = eval_and_free(q); kind = FM_MM_COMPLEX_ABS;
            }
        } else {
            piece = expr_copy(a);
        }
        if (!piece || mm_nonsmooth(piece)) { expr_free(piece); goto fail; }
        for (size_t i = 0; i < n && !any_var; i++)
            if (nm_expr_contains_symbol(piece, vars[i]->data.symbol.name)) any_var = true;
        if (!mm_add_piece(mm, &cap, piece, kind)) { expr_free(piece); goto fail; }
    }
    if (!any_var) goto fail;
    mm->grad = (Expr***)calloc(mm->np, sizeof(Expr**));
    if (!mm->grad) goto fail;
    for (size_t k = 0; k < mm->np; k++) {
        mm->grad[k] = fm_compute_gradient(mm->piece[k], vars, n);
        if (!mm->grad[k]) goto fail;
        for (size_t i = 0; i < n; i++)
            if (!mm->grad[k][i] || mm_nonsmooth(mm->grad[k][i])) goto fail;
    }
    return true;
fail:
    fm_minimax_free(mm);
    return false;
}

/* The epigraph SQP.  x (n) in/out; *fx_out = f at the returned point. */
bool fm_minimax_solve(const FmMinimax* mm, Expr* f, Expr** vars, size_t n,
                      FmVarBind* binds, const FmBox* boxes, const FmOpts* opts,
                      double* x, double* fx_out) {
    bool ok = false;
    size_t n1 = n + 1;
    Expr* tsym = nm_fresh_symbol();
    const char* tname = tsym->data.symbol.name;
    FmVarBind* b1 = (FmVarBind*)malloc(sizeof(FmVarBind) * n1);
    Expr** v1 = (Expr**)malloc(sizeof(Expr*) * n1);
    FmBox* bx1 = (FmBox*)calloc(n1, sizeof(FmBox));
    double* x1 = (double*)malloc(sizeof(double) * n1);
    Expr** gobj = (Expr**)malloc(sizeof(Expr*) * n1);
    size_t ngc = 0;
    for (size_t k = 0; k < mm->np; k++) ngc += (mm->kind[k] == FM_MM_REAL_ABS) ? 2 : 1;
    FmGenCon* gens = (FmGenCon*)calloc(ngc ? ngc : 1, sizeof(FmGenCon));
    CompiledProgram** progs = (CompiledProgram**)calloc(ngc * (n1 + 1) + 1, sizeof(CompiledProgram*));
    const char** cn = (const char**)malloc(sizeof(char*) * n1);
    CompileType* ct = (CompileType*)malloc(sizeof(CompileType) * n1);
    size_t mark = fm_progs_mark();
    bool bound_t = false;
    if (!b1 || !v1 || !bx1 || !x1 || !gobj || !gens || !progs || !cn || !ct) goto done;

    for (size_t i = 0; i < n; i++) {
        b1[i] = binds[i]; v1[i] = vars[i];
        if (boxes) bx1[i] = boxes[i];
        x1[i] = x[i];
        gobj[i] = expr_new_integer(0);
    }
    v1[n] = tsym;
    gobj[n] = expr_new_integer(1);
    fm_bind_snapshot(&b1[n], tname);
    bound_t = true;
    bx1[n].has_lo = true; bx1[n].lo = 0.0;      /* every piece here is >= 0 or not: */
    bool all_abs = true;                         /* t >= 0 only when all are moduli  */
    for (size_t k = 0; k < mm->np; k++) if (mm->kind[k] == FM_MM_SMOOTH) all_abs = false;
    if (!all_abs) bx1[n].has_lo = false;

    /* t at the start: the objective there (feasible start) */
    double f0;
    if (!fm_eval_scalar(f, binds, x, n, opts, &f0) || !isfinite(f0)) goto done;
    for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    x1[n] = f0;

    /* the constraints in (x, t) and their gradients */
    {
        size_t q = 0;
        for (size_t k = 0; k < mm->np; k++) {
            Expr* p = mm->piece[k];
            if (mm->kind[k] == FM_MM_COMPLEX_ABS) {
                /* Sqrt[q] - t, not q - t^2: the modulus of an affine e is
                 * convex, so the Lagrangian Hessian stays positive
                 * semidefinite (q - t^2 puts negative curvature on t, which
                 * the damped BFGS cannot model -- the SQP zig-zagged at a
                 * 1e-4 violation) */
                gens[q++].expr = eval_and_free(mm_fn2(SYM_Plus,
                    mm_fn2(SYM_Power, expr_copy(p),
                           mm_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2))),
                    mm_fn2(SYM_Times, expr_new_integer(-1), expr_copy(tsym))));
            } else {
                gens[q++].expr = eval_and_free(mm_fn2(SYM_Plus, expr_copy(p),
                    mm_fn2(SYM_Times, expr_new_integer(-1), expr_copy(tsym))));
                if (mm->kind[k] == FM_MM_REAL_ABS)
                    gens[q++].expr = eval_and_free(mm_fn2(SYM_Plus,
                        mm_fn2(SYM_Times, expr_new_integer(-1), expr_copy(p)),
                        mm_fn2(SYM_Times, expr_new_integer(-1), expr_copy(tsym))));
            }
        }
        for (size_t c = 0; c < ngc; c++) {
            if (!gens[c].expr) goto done;
            gens[c].equality = false;
            gens[c].grad_exprs = fm_compute_gradient(gens[c].expr, v1, n1);
        }
    }
    /* machine-precision programs for every constraint and gradient component */
    if (opts->prec_mode == FM_PREC_MACHINE) {
        for (size_t i = 0; i < n1; i++) { cn[i] = v1[i]->data.symbol.name; ct[i] = CT_REAL; }
        size_t pi = 0;
        for (size_t c = 0; c < ngc; c++) {
            for (size_t i = 0; i <= n1; i++) {
                Expr* e = (i == 0) ? gens[c].expr
                                   : (gens[c].grad_exprs ? gens[c].grad_exprs[i - 1] : NULL);
                if (!e || e->type == EXPR_INTEGER || e->type == EXPR_REAL) continue;
                CompiledProgram* p = compile_expr_ex(e, cn, ct, n1, COMPILE_FOLD_GLOBALS);
                if (p && compiled_result_type(p) != CT_REAL) { compiled_free(p); p = NULL; }
                if (p) { progs[pi++] = p; fm_progs_register(e, p, n1); }
            }
        }
        CompiledProgram* tp = compile_expr_ex(tsym, cn, ct, n1, COMPILE_FOLD_GLOBALS);
        if (tp && compiled_result_type(tp) != CT_REAL) { compiled_free(tp); tp = NULL; }
        if (tp) { progs[pi++] = tp; fm_progs_register(tsym, tp, n1); }
    }
    /* the SQP, to the precision the smooth form now allows */
    {
        FmOpts o = *opts;
        if (o.prec_mode == FM_PREC_MACHINE) {
            if (o.acc_goal_digits  < 14.0) o.acc_goal_digits  = 14.0;
            if (o.prec_goal_digits < 14.0) o.prec_goal_digits = 14.0;
        }
        if (o.max_iter < 1000) o.max_iter = 1000;
        o.step_monitor = NULL;
        bool saved_quiet = g_fm_quiet;
        g_fm_quiet = true;
        double tv = 0.0;
        bool sok = fm_run_slsqp(tsym, v1, n1, b1, gobj, x1, gens, ngc, 0.0, bx1, &o, &tv);
        g_fm_quiet = saved_quiet;
        for (size_t i = 0; i < n1; i++) fm_bind_clear_temp(&b1[i]);
        if (!sok) for (size_t i = 0; i < n; i++) x1[i] = x[i];   /* Newton from the start */
    }
    {
        double fv;
        if (!fm_eval_scalar(f, binds, x1, n, opts, &fv) || !isfinite(fv)) goto done;
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
        /* keep the start if the SQP ended worse (it returns its last iterate) */
        if (fv < f0) {
            for (size_t i = 0; i < n; i++) x[i] = x1[i];
            f0 = fv;
        }
#ifdef USE_MPFR
        /* The SQP converges onto the active set, sometimes slowly (it
         * zig-zags between near-active pieces); Newton on that active set
         * then gives the double answer to rounding.  It runs at 53 bits: the
         * pieces carry machine numbers, so their values are doubles whatever
         * the MPFR precision of the point, and the convergence and
         * certificate tolerances must be those of a double.  Certified like
         * the WorkingPrecision refinement, adopted only when f improves. */
        if (opts->prec_mode == FM_PREC_MACHINE) {
            FmOpts ho = *opts;
            ho.wp_bits = 53;
            mpfr_t* xr = fm_mpfr_array(n, 53);
            mpfr_t fr; mpfr_init2(fr, 53);
            if (fm_minimax_refine_mpfr(mm, f, vars, n, boxes, &ho, x, xr, fr)) {
                for (size_t i = 0; i < n; i++) x1[i] = mpfr_get_d(xr[i], MPFR_RNDN);
                double fh;
                if (fm_eval_scalar(f, binds, x1, n, opts, &fh) && isfinite(fh) && fh <= f0) {
                    for (size_t i = 0; i < n; i++) x[i] = x1[i];
                    f0 = fh;
                }
                for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
            }
            fm_mpfr_array_free(xr, n);
            mpfr_clear(fr);
        }
#endif
        *fx_out = f0;
        ok = true;
    }
done:
    fm_progs_release(mark);
    if (progs) {
        for (size_t i = 0; i < ngc * (n1 + 1) + 1; i++) if (progs[i]) compiled_free(progs[i]);
        free(progs);
    }
    if (gens) {
        for (size_t c = 0; c < ngc; c++) {
            expr_free(gens[c].expr);
            if (gens[c].grad_exprs) {
                for (size_t i = 0; i < n1; i++) expr_free(gens[c].grad_exprs[i]);
                free(gens[c].grad_exprs);
            }
        }
        free(gens);
    }
    if (gobj) { for (size_t i = 0; i < n1; i++) expr_free(gobj[i]); free(gobj); }
    if (bound_t) { fm_bind_clear_temp(&b1[n]); fm_bind_restore(&b1[n]); }
    free(b1); free(v1); free(bx1); free(x1); free(cn); free(ct);
    symtab_remove_symbol(tname);
    expr_free(tsym);
    return ok;
}

#ifdef USE_MPFR
/* value of e at the point (rules), numericalised at bits; false if not real */
static bool mm_eval_mpfr(const Expr* e, Expr* rules, long bits, mpfr_t out) {
    Expr* a[2] = { expr_copy((Expr*)e), expr_copy(rules) };
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll), a, 2));
    NumericSpec sp; sp.mode = NUMERIC_MODE_MPFR; sp.bits = bits;
    Expr* nv = v ? numericalize(v, sp) : NULL;
    mpfr_t im; mpfr_init2(im, bits);
    bool ok = nv && get_approx_mpfr(nv, out, im, NULL) && mpfr_number_p(out)
              && (mpfr_zero_p(im)          /* a rounding-level imaginary part */
                  || (!mpfr_zero_p(out)
                      && mpfr_get_exp(im) < mpfr_get_exp(out) - (mpfr_exp_t)(bits - 8)));
    mpfr_clear(im);
    expr_free(v); expr_free(nv);
    return ok;
}

static Expr* mm_point_rules(Expr** vars, size_t n, mpfr_t* x) {
    Expr** r = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        r[i] = mm_fn2(SYM_Rule, expr_copy(vars[i]), expr_new_mpfr_copy(x[i]));
    Expr* L = expr_new_function(expr_new_symbol(SYM_List), r, n);
    free(r);
    return L;
}

/* Solve the dense m x m system A z = b in place (partial pivoting).  A and b
 * are destroyed; z receives the solution.  False if singular. */
static bool mm_solve_mpfr(mpfr_t* A, mpfr_t* b, size_t m, mpfr_t* z, long bits) {
    mpfr_t t, fct;
    mpfr_inits2(bits, t, fct, (mpfr_ptr)0);
    bool ok = true;
    for (size_t c = 0; c < m && ok; c++) {
        size_t p = c;
        for (size_t r = c + 1; r < m; r++)
            if (mpfr_cmpabs(A[r * m + c], A[p * m + c]) > 0) p = r;
        if (mpfr_zero_p(A[p * m + c])) { ok = false; break; }
        if (p != c) {
            for (size_t j = 0; j < m; j++) mpfr_swap(A[p * m + j], A[c * m + j]);
            mpfr_swap(b[p], b[c]);
        }
        for (size_t r = c + 1; r < m; r++) {
            if (mpfr_zero_p(A[r * m + c])) continue;
            mpfr_div(fct, A[r * m + c], A[c * m + c], MPFR_RNDN);
            for (size_t j = c; j < m; j++) {
                mpfr_mul(t, fct, A[c * m + j], MPFR_RNDN);
                mpfr_sub(A[r * m + j], A[r * m + j], t, MPFR_RNDN);
            }
            mpfr_mul(t, fct, b[c], MPFR_RNDN);
            mpfr_sub(b[r], b[r], t, MPFR_RNDN);
        }
    }
    for (size_t c = m; ok && c-- > 0; ) {
        mpfr_set(t, b[c], MPFR_RNDN);
        for (size_t j = c + 1; j < m; j++) {
            mpfr_mul(fct, A[c * m + j], z[j], MPFR_RNDN);
            mpfr_sub(t, t, fct, MPFR_RNDN);
        }
        mpfr_div(z[c], t, A[c * m + c], MPFR_RNDN);
    }
    mpfr_clears(t, fct, (mpfr_ptr)0);
    return ok;
}

/* the same at machine precision (for the multiplier check) */
static bool mm_solve_d(double* A, double* b, size_t m, double* z) {
    for (size_t c = 0; c < m; c++) {
        size_t p = c;
        for (size_t r = c + 1; r < m; r++) if (fabs(A[r * m + c]) > fabs(A[p * m + c])) p = r;
        if (A[p * m + c] == 0.0) return false;
        if (p != c) {
            for (size_t j = 0; j < m; j++) { double s = A[p * m + j]; A[p * m + j] = A[c * m + j]; A[c * m + j] = s; }
            double s = b[p]; b[p] = b[c]; b[c] = s;
        }
        for (size_t r = c + 1; r < m; r++) {
            double fct = A[r * m + c] / A[c * m + c];
            for (size_t j = c; j < m; j++) A[r * m + j] -= fct * A[c * m + j];
            b[r] -= fct * b[c];
        }
    }
    for (size_t c = m; c-- > 0; ) {
        double t = b[c];
        for (size_t j = c + 1; j < m; j++) t -= A[c * m + j] * z[j];
        z[c] = t / A[c * m + c];
    }
    return true;
}

/* All piece values at x (moduli / absolute values), their maximum in *tmax;
 * sg[k] the sign of a real piece. */
static bool mm_values(const FmMinimax* mm, Expr** vars, size_t n, mpfr_t* x,
                      long bits, mpfr_t* val, int* sg, mpfr_t tmax) {
    Expr* R = mm_point_rules(vars, n, x);
    bool ok = true;
    for (size_t k = 0; k < mm->np; k++) {
        if (!mm_eval_mpfr(mm->piece[k], R, bits, val[k])) { ok = false; break; }
        if (sg) sg[k] = mpfr_sgn(val[k]) < 0 ? -1 : 1;
        if (mm->kind[k] == FM_MM_REAL_ABS) mpfr_abs(val[k], val[k], MPFR_RNDN);
        else if (mm->kind[k] == FM_MM_COMPLEX_ABS) {
            if (mpfr_sgn(val[k]) < 0) mpfr_set_zero(val[k], 1);
            mpfr_sqrt(val[k], val[k], MPFR_RNDN);
        }
        if (k == 0 || mpfr_cmp(val[k], tmax) > 0) mpfr_set(tmax, val[k], MPFR_RNDN);
    }
    expr_free(R);
    return ok;
}

/* Newton / Gauss-Newton on the active set A (na >= n + 1 pieces) from x =
 * (x, t) in place, then the certificate.  True when it converged to a zero
 * residual, no piece exceeds t, and (square case) the multipliers are
 * non-negative; *tmax is then the maximum of all pieces at x. */
/* *why on failure: MM_NEWTON_FAIL, MM_VIOLATED (*which = the piece above t),
 * MM_NEGATIVE (*which = the position in A of the most negative multiplier). */
enum { MM_NEWTON_FAIL = 0, MM_VIOLATED = 1, MM_NEGATIVE = 2, MM_CERTIFIED = 3 };
static bool mm_newton(const FmMinimax* mm, Expr** vars, size_t n, const FmBox* boxes,
                      long bits, const size_t* A, size_t na, const int* sg,
                      mpfr_t* x, mpfr_t* val, mpfr_t tmax, int* why, size_t* which) {
    size_t n1 = n + 1;
    bool ok = false;
    *why = MM_NEWTON_FAIL; *which = 0;
    mpfr_t* J   = fm_mpfr_array(na * n1, bits);
    mpfr_t* F   = fm_mpfr_array(na, bits);
    mpfr_t* N   = fm_mpfr_array(n1 * n1, bits);
    mpfr_t* rhs = fm_mpfr_array(n1, bits);
    mpfr_t* dz  = fm_mpfr_array(n1, bits);
    mpfr_t tmp, tol, lim;
    mpfr_inits2(bits, tmp, tol, lim, (mpfr_ptr)0);
    mpfr_set_ui_2exp(tol, 1, -(mpfr_exp_t)(bits - 12), MPFR_RNDN);
    bool conv = false;
    for (int it = 0; it < 40 && !conv; it++) {
        Expr* R = mm_point_rules(vars, n, x);
        bool eok = true;
        for (size_t a = 0; a < na && eok; a++) {
            size_t k = A[a];
            if (!mm_eval_mpfr(mm->piece[k], R, bits, F[a])) { eok = false; break; }
            for (size_t i = 0; i < n && eok; i++)
                if (!mm_eval_mpfr(mm->grad[k][i], R, bits, J[a * n1 + i])) eok = false;
            if (!eok) break;
            if (mm->kind[k] == FM_MM_COMPLEX_ABS) {           /* q - t^2 */
                mpfr_sqr(tmp, x[n], MPFR_RNDN);
                mpfr_sub(F[a], F[a], tmp, MPFR_RNDN);
                mpfr_mul_si(J[a * n1 + n], x[n], -2, MPFR_RNDN);
            } else {                                          /* s p - t */
                if (mm->kind[k] == FM_MM_REAL_ABS && sg[k] < 0) {
                    mpfr_neg(F[a], F[a], MPFR_RNDN);
                    for (size_t i = 0; i < n; i++) mpfr_neg(J[a * n1 + i], J[a * n1 + i], MPFR_RNDN);
                }
                mpfr_sub(F[a], F[a], x[n], MPFR_RNDN);
                mpfr_set_si(J[a * n1 + n], -1, MPFR_RNDN);
            }
        }
        expr_free(R);
        if (!eok) goto done;
        /* normal equations J^T J dz = -J^T F (Newton itself when square) */
        for (size_t i = 0; i < n1; i++) {
            mpfr_set_zero(rhs[i], 1);
            for (size_t a = 0; a < na; a++) {
                mpfr_mul(tmp, J[a * n1 + i], F[a], MPFR_RNDN);
                mpfr_sub(rhs[i], rhs[i], tmp, MPFR_RNDN);
            }
            for (size_t j = 0; j < n1; j++) {
                mpfr_set_zero(N[i * n1 + j], 1);
                for (size_t a = 0; a < na; a++) {
                    mpfr_mul(tmp, J[a * n1 + i], J[a * n1 + j], MPFR_RNDN);
                    mpfr_add(N[i * n1 + j], N[i * n1 + j], tmp, MPFR_RNDN);
                }
            }
        }
        if (!mm_solve_mpfr(N, rhs, n1, dz, bits)) goto done;
        bool small = true;
        for (size_t i = 0; i < n1; i++) {
            mpfr_add(x[i], x[i], dz[i], MPFR_RNDN);
            mpfr_abs(tmp, x[i], MPFR_RNDN);
            mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
            mpfr_mul(tmp, tmp, tol, MPFR_RNDN);
            if (mpfr_cmpabs(dz[i], tmp) > 0) small = false;
        }
        if (small) conv = true;
    }
    if (!conv) goto done;
    /* a ZERO residual: Gauss-Newton also converges, to a least-squares point,
     * on an over-determined set that holds a piece which is not truly active */
    for (size_t a = 0; a < na; a++) {
        mpfr_abs(lim, x[n], MPFR_RNDN);
        mpfr_add_ui(lim, lim, 1, MPFR_RNDN);
        mpfr_mul_2si(lim, lim, -(long)(bits / 2), MPFR_RNDN);
        if (mpfr_cmpabs(F[a], lim) > 0) goto done;
    }
    /* the refined t is the maximum of ALL pieces */
    if (!mm_values(mm, vars, n, x, bits, val, NULL, tmax)) goto done;
    mpfr_sub(tmp, tmax, x[n], MPFR_RNDN);
    mpfr_abs(lim, x[n], MPFR_RNDN);
    mpfr_add_ui(lim, lim, 1, MPFR_RNDN);
    mpfr_mul_2si(lim, lim, -(long)(bits - 24), MPFR_RNDN);
    if (mpfr_cmpabs(tmp, lim) > 0) {
        *why = MM_VIOLATED;
        for (size_t k = 0; k < mm->np; k++)
            if (mpfr_equal_p(val[k], tmax)) { *which = k; break; }
        goto done;
    }
    /* KKT multipliers (square case): J^T lambda = -e_t, lambda >= 0 */
    if (na == n1) {
        double* Jt = (double*)malloc(sizeof(double) * n1 * n1);
        double* b  = (double*)malloc(sizeof(double) * n1);
        double* lam = (double*)malloc(sizeof(double) * n1);
        bool good = Jt && b && lam;
        for (size_t i = 0; good && i < n1; i++) {
            for (size_t a = 0; a < na; a++) Jt[i * n1 + a] = mpfr_get_d(J[a * n1 + i], MPFR_RNDN);
            b[i] = (i == n) ? -1.0 : 0.0;
        }
        if (good) good = mm_solve_d(Jt, b, n1, lam);
        double lmax = 0.0, lmin = 0.0;
        for (size_t a = 0; good && a < na; a++) if (fabs(lam[a]) > lmax) lmax = fabs(lam[a]);
        for (size_t a = 0; good && a < na; a++)
            if (lam[a] < -1e-8 * lmax) {
                if (good || lam[a] < lmin) { lmin = lam[a]; *which = a; }
                good = false;
            }
        bool negative = Jt && b && lam && !good && lmin < 0.0;
        free(Jt); free(b); free(lam);
        if (!good) { *why = negative ? MM_NEGATIVE : MM_NEWTON_FAIL; goto done; }
    }
    for (size_t i = 0; boxes && i < n; i++) {
        double d = mpfr_get_d(x[i], MPFR_RNDN);
        if ((boxes[i].has_lo && d < boxes[i].lo) || (boxes[i].has_hi && d > boxes[i].hi)) goto done;
    }
    ok = true;
    *why = MM_CERTIFIED;
done:
    fm_mpfr_array_free(J, na * n1);
    fm_mpfr_array_free(F, na);
    fm_mpfr_array_free(N, n1 * n1);
    fm_mpfr_array_free(rhs, n1);
    fm_mpfr_array_free(dz, n1);
    mpfr_clears(tmp, tol, lim, (mpfr_ptr)0);
    return ok;
}

/* KKT Newton for an active set SMALLER than the unknowns (na < n + 1): the
 * optimum of a symmetric problem often ties fewer pieces than there are
 * coefficients (a real function on a conjugate-symmetric point set, fitted
 * with complex coefficients).  The equations piece_a = t alone are then
 * under-determined; with the multipliers they are square:
 *     sum_a lambda_a grad_x F_a = 0,  1 + sum_a lambda_a dF_a/dt = 0,  F_a = 0
 * in (x, t, lambda), F_a = s p_a - t (or q_a - t^2).  Hessians of the active
 * pieces are differentiated on first use (H[k], cached by the caller).  The
 * certificate is the one mm_newton uses, the multipliers here being the
 * Newton unknowns themselves. */
static bool mm_kkt_newton(const FmMinimax* mm, Expr** vars, size_t n, const FmBox* boxes,
                          long bits, const size_t* A, size_t na, const int* sg,
                          Expr*** H, mpfr_t* x, mpfr_t* val, mpfr_t tmax) {
    size_t n1 = n + 1, M = n1 + na;
    bool ok = false;
    for (size_t a = 0; a < na; a++) {
        size_t k = A[a];
        if (H[k]) continue;
        H[k] = (Expr**)calloc(n * n, sizeof(Expr*));
        if (!H[k]) return false;
        for (size_t i = 0; i < n; i++) {
            Expr** row = fm_compute_gradient(mm->grad[k][i], vars, n);
            if (!row) return false;
            for (size_t j = 0; j < n; j++) H[k][i * n + j] = row[j];
            free(row);
        }
    }
    mpfr_t* g   = fm_mpfr_array(na * n, bits);      /* signed grad_x of piece a */
    mpfr_t* Hv  = fm_mpfr_array(na * n * n, bits);  /* signed Hessians          */
    mpfr_t* Fv  = fm_mpfr_array(na, bits);
    mpfr_t* dt  = fm_mpfr_array(na, bits);
    mpfr_t* lam = fm_mpfr_array(na, bits);
    mpfr_t* Jm  = fm_mpfr_array(M * M, bits);
    mpfr_t* r   = fm_mpfr_array(M, bits);
    mpfr_t* dz  = fm_mpfr_array(M, bits);
    mpfr_t tmp, tol, lim;
    mpfr_inits2(bits, tmp, tol, lim, (mpfr_ptr)0);
    mpfr_set_ui_2exp(tol, 1, -(mpfr_exp_t)(bits - 12), MPFR_RNDN);
    bool conv = false;
    for (int it = 0; it < 40 && !conv; it++) {
        Expr* R = mm_point_rules(vars, n, x);
        bool eok = true;
        for (size_t a = 0; a < na && eok; a++) {
            size_t k = A[a];
            bool neg = (mm->kind[k] == FM_MM_REAL_ABS && sg[k] < 0);
            eok = mm_eval_mpfr(mm->piece[k], R, bits, Fv[a]);
            for (size_t i = 0; i < n && eok; i++)
                eok = mm_eval_mpfr(mm->grad[k][i], R, bits, g[a * n + i]);
            for (size_t i = 0; i < n * n && eok; i++)
                eok = mm_eval_mpfr(H[k][i], R, bits, Hv[a * n * n + i]);
            if (!eok) break;
            if (neg) {
                mpfr_neg(Fv[a], Fv[a], MPFR_RNDN);
                for (size_t i = 0; i < n; i++) mpfr_neg(g[a * n + i], g[a * n + i], MPFR_RNDN);
                for (size_t i = 0; i < n * n; i++) mpfr_neg(Hv[a * n * n + i], Hv[a * n * n + i], MPFR_RNDN);
            }
            if (mm->kind[k] == FM_MM_COMPLEX_ABS) {
                mpfr_sqr(tmp, x[n], MPFR_RNDN);
                mpfr_sub(Fv[a], Fv[a], tmp, MPFR_RNDN);
                mpfr_mul_si(dt[a], x[n], -2, MPFR_RNDN);
            } else {
                mpfr_sub(Fv[a], Fv[a], x[n], MPFR_RNDN);
                mpfr_set_si(dt[a], -1, MPFR_RNDN);
            }
        }
        expr_free(R);
        if (!eok) goto done;
        if (it == 0) {
            /* least-squares multipliers: (G G^T) lambda = -dt, G_a = (g_a, dt_a) */
            mpfr_t* GG = fm_mpfr_array(na * na, bits);
            mpfr_t* rh = fm_mpfr_array(na, bits);
            for (size_t a = 0; a < na; a++) {
                mpfr_neg(rh[a], dt[a], MPFR_RNDN);
                for (size_t b = 0; b < na; b++) {
                    mpfr_mul(GG[a * na + b], dt[a], dt[b], MPFR_RNDN);
                    for (size_t i = 0; i < n; i++) {
                        mpfr_mul(tmp, g[a * n + i], g[b * n + i], MPFR_RNDN);
                        mpfr_add(GG[a * na + b], GG[a * na + b], tmp, MPFR_RNDN);
                    }
                }
            }
            bool sok = mm_solve_mpfr(GG, rh, na, lam, bits);
            fm_mpfr_array_free(GG, na * na);
            fm_mpfr_array_free(rh, na);
            if (!sok) goto done;
        }
        /* residual r and Jacobian Jm (row-major M x M); unknowns (x, t, lambda) */
        for (size_t t2 = 0; t2 < M * M; t2++) mpfr_set_zero(Jm[t2], 1);
        for (size_t i = 0; i < n; i++) {
            mpfr_set_zero(r[i], 1);
            for (size_t a = 0; a < na; a++) {
                mpfr_mul(tmp, lam[a], g[a * n + i], MPFR_RNDN);
                mpfr_add(r[i], r[i], tmp, MPFR_RNDN);
                for (size_t j = 0; j < n; j++) {
                    mpfr_mul(tmp, lam[a], Hv[a * n * n + i * n + j], MPFR_RNDN);
                    mpfr_add(Jm[i * M + j], Jm[i * M + j], tmp, MPFR_RNDN);
                }
                mpfr_set(Jm[i * M + n1 + a], g[a * n + i], MPFR_RNDN);
            }
        }
        mpfr_set_ui(r[n], 1, MPFR_RNDN);
        for (size_t a = 0; a < na; a++) {
            mpfr_mul(tmp, lam[a], dt[a], MPFR_RNDN);
            mpfr_add(r[n], r[n], tmp, MPFR_RNDN);
            if (mm->kind[A[a]] == FM_MM_COMPLEX_ABS) {
                mpfr_mul_si(tmp, lam[a], -2, MPFR_RNDN);
                mpfr_add(Jm[n * M + n], Jm[n * M + n], tmp, MPFR_RNDN);
            }
            mpfr_set(Jm[n * M + n1 + a], dt[a], MPFR_RNDN);
            mpfr_set(r[n1 + a], Fv[a], MPFR_RNDN);
            for (size_t j = 0; j < n; j++) mpfr_set(Jm[(n1 + a) * M + j], g[a * n + j], MPFR_RNDN);
            mpfr_set(Jm[(n1 + a) * M + n], dt[a], MPFR_RNDN);
        }
        for (size_t i = 0; i < M; i++) mpfr_neg(r[i], r[i], MPFR_RNDN);
        if (!mm_solve_mpfr(Jm, r, M, dz, bits)) goto done;
        bool small = true;
        for (size_t i = 0; i < M; i++) {
            mpfr_t* zi = (i < n1) ? &x[i] : &lam[i - n1];
            mpfr_add(*zi, *zi, dz[i], MPFR_RNDN);
            if (i >= n1) continue;
            mpfr_abs(tmp, x[i], MPFR_RNDN);
            mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
            mpfr_mul(tmp, tmp, tol, MPFR_RNDN);
            if (mpfr_cmpabs(dz[i], tmp) > 0) small = false;
        }
        if (small && it > 0) conv = true;
    }
    if (!conv) goto done;
    /* certificate: feasibility residual, no piece above t, lambda >= 0 */
    for (size_t a = 0; a < na; a++) {
        mpfr_abs(lim, x[n], MPFR_RNDN);
        mpfr_add_ui(lim, lim, 1, MPFR_RNDN);
        mpfr_mul_2si(lim, lim, -(long)(bits / 2), MPFR_RNDN);
        if (mpfr_cmpabs(Fv[a], lim) > 0) goto done;
    }
    if (!mm_values(mm, vars, n, x, bits, val, NULL, tmax)) goto done;
    mpfr_sub(tmp, tmax, x[n], MPFR_RNDN);
    mpfr_abs(lim, x[n], MPFR_RNDN);
    mpfr_add_ui(lim, lim, 1, MPFR_RNDN);
    mpfr_mul_2si(lim, lim, -(long)(bits - 24), MPFR_RNDN);
    if (mpfr_cmpabs(tmp, lim) > 0) goto done;
    {
        double lmax = 0.0;
        for (size_t a = 0; a < na; a++) {
            double l = mpfr_get_d(lam[a], MPFR_RNDN);
            if (fabs(l) > lmax) lmax = fabs(l);
        }
        for (size_t a = 0; a < na; a++)
            if (mpfr_get_d(lam[a], MPFR_RNDN) < -1e-8 * lmax) goto done;
    }
    for (size_t i = 0; boxes && i < n; i++) {
        double d = mpfr_get_d(x[i], MPFR_RNDN);
        if ((boxes[i].has_lo && d < boxes[i].lo) || (boxes[i].has_hi && d > boxes[i].hi)) goto done;
    }
    ok = true;
done:
    fm_mpfr_array_free(g, na * n);
    fm_mpfr_array_free(Hv, na * n * n);
    fm_mpfr_array_free(Fv, na);
    fm_mpfr_array_free(dt, na);
    fm_mpfr_array_free(lam, na);
    fm_mpfr_array_free(Jm, M * M);
    fm_mpfr_array_free(r, M);
    fm_mpfr_array_free(dz, M);
    mpfr_clears(tmp, tol, lim, (mpfr_ptr)0);
    return ok;
}

bool fm_minimax_refine_mpfr(const FmMinimax* mm, Expr* f, Expr** vars, size_t n,
                            const FmBox* boxes, const FmOpts* opts,
                            const double* xm, mpfr_t* x_out, mpfr_t fx_out) {
    (void)f;
    long bits = opts->wp_bits;
    size_t n1 = n + 1, np = mm->np;
    bool ok = false;
    /* a coordinate on a box bound is outside the equation system's reach */
    for (size_t i = 0; boxes && i < n; i++) {
        if (boxes[i].has_lo && xm[i] <= boxes[i].lo + 1e-9 * (1.0 + fabs(boxes[i].lo))) return false;
        if (boxes[i].has_hi && xm[i] >= boxes[i].hi - 1e-9 * (1.0 + fabs(boxes[i].hi))) return false;
    }
    mpfr_t* x0 = fm_mpfr_array(n1, bits);
    mpfr_t* x = fm_mpfr_array(n1, bits);
    mpfr_t* val = fm_mpfr_array(np ? np : 1, bits);
    size_t* A = (size_t*)malloc(sizeof(size_t) * (np ? np : 1));
    size_t* B = (size_t*)malloc(sizeof(size_t) * (np ? np : 1));
    int* sg = (int*)malloc(sizeof(int) * (np ? np : 1));
    double* vd = (double*)malloc(sizeof(double) * (np ? np : 1));
    mpfr_t tmax;
    mpfr_init2(tmax, bits);
    int why = 0; size_t which = 0;
    Expr*** H = (Expr***)calloc(np ? np : 1, sizeof(Expr**));   /* Hessian cache */
    if (!A || !sg || !vd) goto done;
    for (size_t i = 0; i < n; i++) mpfr_set_d(x0[i], xm[i], MPFR_RNDN);
    if (!mm_values(mm, vars, n, x0, bits, val, sg, tmax)) goto done;
    mpfr_set(x0[n], tmax, MPFR_RNDN);
    for (size_t k = 0; k < np; k++) vd[k] = mpfr_get_d(val[k], MPFR_RNDN);
    /* The active set: the pieces within a relative gap of the maximum.  How
     * close the start is decides how wide that gap must be, so the ladder
     * widens it until a set solves and certifies (a set with a piece that is
     * not active fails the zero-residual test; one missing a piece fails the
     * every-piece-below-t test). */
    {
        static const double thr[] = { 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4 };
        double tm = mpfr_get_d(tmax, MPFR_RNDN);
        size_t last_na = 0;
        for (size_t ti = 0; ti < sizeof thr / sizeof thr[0] && !ok; ti++) {
            size_t na = 0;
            for (size_t k = 0; k < np; k++)
                if (vd[k] >= tm - thr[ti] * fabs(tm)) A[na++] = k;
            if (na == 0 || na == last_na) continue;
            last_na = na;
            for (size_t i = 0; i < n1; i++) mpfr_set(x[i], x0[i], MPFR_RNDN);
            if (na < n1) {          /* degenerate: with the multipliers */
                if (H) ok = mm_kkt_newton(mm, vars, n, boxes, bits, A, na, sg, H, x, val, tmax);
                continue;
            }
            ok = mm_newton(mm, vars, n, boxes, bits, A, na, sg, x, val, tmax, &why, &which);
            /* An over-full set holds a near-active piece that is not active
             * at the optimum (the start sits ~1e-8 from it, the gap to the
             * next piece can be smaller): drop one or two pieces in turn,
             * the exchange step of the Remez algorithm.  The certificate
             * decides, so a wrong drop cannot be returned. */
            size_t extra = na - n1;
            if (!ok && extra >= 1 && extra <= 2 && B) {
                for (size_t d1 = 0; d1 < na && !ok; d1++) {
                    for (size_t d2 = (extra == 2 ? d1 + 1 : na); d2 <= na && !ok; d2++) {
                        if (extra == 2 && d2 == na) break;
                        size_t nb = 0;
                        for (size_t a = 0; a < na; a++)
                            if (a != d1 && a != d2) B[nb++] = A[a];
                        for (size_t i = 0; i < n1; i++) mpfr_set(x[i], x0[i], MPFR_RNDN);
                        ok = mm_newton(mm, vars, n, boxes, bits, B, nb, sg, x, val, tmax, &why, &which);
                        if (extra == 1) break;
                    }
                }
            }
        }
    }
    if (ok) {
        for (size_t i = 0; i < n; i++) mpfr_set(x_out[i], x[i], MPFR_RNDN);
        mpfr_set(fx_out, tmax, MPFR_RNDN);
    }
done:
    fm_mpfr_array_free(x0, n1);
    fm_mpfr_array_free(x, n1);
    fm_mpfr_array_free(val, np ? np : 1);
    mpfr_clear(tmax);
    if (H) {
        for (size_t k = 0; k < np; k++) if (H[k]) {
            for (size_t i = 0; i < n * n; i++) expr_free(H[k][i]);
            free(H[k]);
        }
        free(H);
    }
    free(A); free(B); free(sg); free(vd);
    return ok;
}
#endif /* USE_MPFR */
