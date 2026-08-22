/* findmin_driver.c — FindMinimum / FindMaximum driver + builtins.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Driver                                                              *
 * ------------------------------------------------------------------ */

static Expr* findmin_driver(Expr* res, const char* fn_name) {
    g_fm_name = fn_name;
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn(fn_name, "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && fm_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!fm_is_option_arg(res->data.function.args[i])) {
            fm_warn(fn_name, "badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) {
        fm_warn(fn_name, "argt", "needs exactly 2 positional arguments (got %zu)", pos_end);
        return NULL;
    }
    FmOpts opts;
    opts.method = FM_METHOD_AUTOMATIC;
    opts.prec_mode = FM_PREC_MACHINE;
    opts.wp_bits = 0;
    opts.max_iter = 500;
    opts.max_iter_set = false;
    opts.acc_goal_digits = -1.0;
    opts.prec_goal_digits = -1.0;
    opts.gradient = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;
    for (size_t i = pos_end; i < argc; i++) {
        if (!fm_apply_option(res->data.function.args[i], &opts)) return NULL;
    }
    double wp_digits = (opts.prec_mode == FM_PREC_MACHINE)
        ? NUMERIC_MACHINE_PRECISION_DIGITS
#ifdef USE_MPFR
        : numeric_bits_to_digits(opts.wp_bits);
#else
        : NUMERIC_MACHINE_PRECISION_DIGITS;
#endif
    if (opts.acc_goal_digits  < 0.0) opts.acc_goal_digits  = wp_digits / 2.0;
    if (opts.prec_goal_digits < 0.0) opts.prec_goal_digits = wp_digits / 2.0;

    /* Detect {f, cons} form. */
    Expr* f_arg = res->data.function.args[0];
    Expr* var_arg = res->data.function.args[1];
    Expr* f_raw = f_arg;
    Expr* cons = NULL;
    if (f_arg->type == EXPR_FUNCTION
        && f_arg->data.function.head->type == EXPR_SYMBOL
        && f_arg->data.function.head->data.symbol.name == SYM_List
        && f_arg->data.function.arg_count == 2) {
        f_raw = f_arg->data.function.args[0];
        cons = f_arg->data.function.args[1];
    }

    /* Parse variables. var_arg may be:
     *   {x}  /  {x, x0}  /  {x, x0, x1}  /  {x, xstart, xmin, xmax}      (scalar)
     *   {{x, ...}, {y, ...}, ...}                                        (vector)
     *   {x, y, z}  (each scalar element is a bare symbol, treated as {x_i, 0}) */
    bool is_system = false;
    if (var_arg->type == EXPR_FUNCTION
        && var_arg->data.function.head->type == EXPR_SYMBOL
        && var_arg->data.function.head->data.symbol.name == SYM_List
        && var_arg->data.function.arg_count > 0) {
        size_t na = var_arg->data.function.arg_count;
        bool any_inner = false, all_inner_or_sym = true;
        for (size_t i = 0; i < na; i++) {
            Expr* e = var_arg->data.function.args[i];
            bool is_inner = (e->type == EXPR_FUNCTION
                && e->data.function.head->type == EXPR_SYMBOL
                && e->data.function.head->data.symbol.name == SYM_List);
            bool is_sym = (e->type == EXPR_SYMBOL);
            if (is_inner) any_inner = true;
            if (!is_inner && !is_sym) all_inner_or_sym = false;
        }
        /* {{x,x0},{y,y0}} → system; {x, y, z} (all bare symbols) → system n-D */
        if (any_inner && all_inner_or_sym) is_system = true;
        else if (na >= 2 && var_arg->data.function.args[0]->type == EXPR_SYMBOL) {
            /* Could be either {x, x0} (scalar) or {x, y} (multi-symbol). Check 2nd arg. */
            Expr* a1 = var_arg->data.function.args[1];
            if (a1->type == EXPR_SYMBOL && na > 1) {
                /* {x, y, ...}: all bare symbols → multi-var auto-start. */
                bool all_sym = true;
                for (size_t i = 0; i < na; i++) {
                    if (var_arg->data.function.args[i]->type != EXPR_SYMBOL) {
                        all_sym = false; break;
                    }
                }
                if (all_sym) is_system = true;
            }
        }
    }

    Expr** vars = NULL;
    double* x_vec = NULL;
    FmBox* boxes = NULL;
    FmVarBind* binds = NULL;
    FmGenCon* gens = NULL;
    size_t ngens = 0, gcap = 0;
    Expr** g_exprs = NULL;
    Expr*** H_exprs = NULL;
    Expr* result_out = NULL;
    CompiledProgram*  f_prog = NULL;      /* machine-precision compiled objective */
    CompiledProgram** grad_progs = NULL;  /* per-component compiled exact gradient */
    size_t n = 0;

    if (is_system) {
        n = var_arg->data.function.arg_count;
        vars = (Expr**)calloc(n, sizeof(Expr*));
        x_vec = (double*)calloc(n, sizeof(double));
        boxes = (FmBox*)calloc(n, sizeof(FmBox));
        for (size_t i = 0; i < n; i++) {
            Expr* sub = var_arg->data.function.args[i];
            Expr *u, *x0e = NULL, *x1e = NULL, *xmin = NULL, *xmax = NULL;
            FmSpecKind k = fm_parse_var_spec(sub, &u, &x0e, &x1e, &xmin, &xmax);
            if (k != FM_SPEC_VAR_ONLY && k != FM_SPEC_SINGLE && k != FM_SPEC_BRACKET) {
                fm_warn(fn_name, "ivar", "variable spec %zu malformed", i);
                expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
                goto cleanup;
            }
            vars[i] = u;
            if (!x0e || !fm_expr_to_double_real(x0e, &x_vec[i])) x_vec[i] = 0.0;
            if (k == FM_SPEC_BRACKET) {
                double lo, hi;
                if (fm_expr_to_double_real(xmin, &lo) && fm_expr_to_double_real(xmax, &hi)) {
                    boxes[i].has_lo = true; boxes[i].lo = lo;
                    boxes[i].has_hi = true; boxes[i].hi = hi;
                }
            }
            expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
        }
    } else {
        n = 1;
        vars = (Expr**)calloc(1, sizeof(Expr*));
        x_vec = (double*)calloc(1, sizeof(double));
        boxes = (FmBox*)calloc(1, sizeof(FmBox));
        Expr *u, *x0e = NULL, *x1e = NULL, *xmin = NULL, *xmax = NULL;
        FmSpecKind k = fm_parse_var_spec(var_arg, &u, &x0e, &x1e, &xmin, &xmax);
        if (k == FM_SPEC_BAD) {
            fm_warn(fn_name, "ivar", "variable spec must be {x}, {x, x0}, {x, x0, x1}, or {x, xstart, xmin, xmax}");
            expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
            goto cleanup;
        }
        vars[0] = u;
        if (!x0e || !fm_expr_to_double_real(x0e, &x_vec[0])) x_vec[0] = 0.0;
        if (k == FM_SPEC_BRACKET) {
            double lo, hi;
            if (fm_expr_to_double_real(xmin, &lo) && fm_expr_to_double_real(xmax, &hi)) {
                boxes[0].has_lo = true; boxes[0].lo = lo;
                boxes[0].has_hi = true; boxes[0].hi = hi;
            }
        }
        /* Smuggle TWO_START into method selection via custom path: we
         * encode this by using boxes (lo=x0, hi=x1) if Brent and the
         * caller gave {var, x0, x1}. We'll handle it during method
         * dispatch below by detecting TWO_START separately. */
        if (k == FM_SPEC_TWO_START) {
            double a, b;
            if (fm_expr_to_double_real(x0e, &a) && fm_expr_to_double_real(x1e, &b)) {
                if (a > b) { double t = a; a = b; b = t; }
                boxes[0].has_lo = true; boxes[0].lo = a;
                boxes[0].has_hi = true; boxes[0].hi = b;
                /* For TWO_START with Automatic method we want Brent. */
                if (opts.method == FM_METHOD_AUTOMATIC) opts.method = FM_METHOD_BRENT;
            }
        }
        expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
    }

    /* Now bind variables. */
    binds = (FmVarBind*)calloc(n, sizeof(FmVarBind));
    for (size_t i = 0; i < n; i++) fm_bind_snapshot(&binds[i], vars[i]->data.symbol.name);

    /* Constraints. FindMinimum passes no disjunction sink: its smooth gradient
     * penalty method cannot use the non-smooth min penalty an Or requires. */
    if (cons) {
        if (!fm_collect_constraints(cons, vars, n, boxes, &gens, &ngens, &gcap,
                                    NULL, NULL, NULL))
            goto cleanup;
        /* Best-effort symbolic gradient of each constraint expression. The
         * penalty solver needs ∇(f + μ·Σ penalty) — using a stale ∇f alone
         * gives the inner BFGS/CG an inconsistent value/gradient pair and
         * the penalty term loses all influence on the descent direction. */
        for (size_t k = 0; k < ngens; k++) {
            gens[k].grad_exprs = fm_compute_gradient(gens[k].expr, vars, n);
            /* NULL is fine — fm_eval_aug_gradient will FD that constraint. */
        }
    }

    /* Method selection. */
    FmMethod method = opts.method;
    if (method == FM_METHOD_AUTOMATIC) {
        method = (n == 1) ? FM_METHOD_BRENT : FM_METHOD_QUASINEWTON;
    }

    /* Compute symbolic gradient/Hessian when needed. */
    bool needs_grad = (method == FM_METHOD_QUASINEWTON
                    || method == FM_METHOD_CONJGRAD
                    || method == FM_METHOD_NEWTON
                    || method == FM_METHOD_LBFGSB
                    || method == FM_METHOD_TNC
                    || method == FM_METHOD_SLSQP
                    || method == FM_METHOD_NEWTONCG
                    || method == FM_METHOD_DOGLEG
                    || method == FM_METHOD_TRUSTNCG
                    || method == FM_METHOD_TRUSTEXACT
                    || method == FM_METHOD_TRUSTKRYLOV);
    bool needs_hess = (method == FM_METHOD_NEWTON
                    || method == FM_METHOD_DOGLEG
                    || method == FM_METHOD_TRUSTEXACT);
    if (needs_grad) {
        if (opts.gradient
            && opts.gradient->type == EXPR_FUNCTION
            && opts.gradient->data.function.head->type == EXPR_SYMBOL
            && opts.gradient->data.function.head->data.symbol.name == SYM_List
            && opts.gradient->data.function.arg_count == n) {
            g_exprs = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) g_exprs[i] = expr_copy(opts.gradient->data.function.args[i]);
        } else {
            g_exprs = fm_compute_gradient(f_raw, vars, n);
            if (!g_exprs) {
                /* OK — will fall back to finite differences inside the solver. */
            }
        }
    }
    if (needs_hess) {
        H_exprs = fm_compute_hessian(f_raw, vars, n);
        /* OK if NULL — Newton will fall back to BFGS-style steepest. */
    }

    /* Machine-precision auto-compilation of the objective. The local solvers
     * evaluate f at every trial point (line search, bracketing, function
     * values); lowering it to bytecode once over the variables and running the
     * register machine per point is far cheaper than the interpreter
     * (expr_copy + n OwnValue installs + evaluate + numericalize). Registered in
     * the g_fm_obj_* slots that fm_eval_scalar consults; deregistered and freed
     * at cleanup. Compiled *after* fm_bind_snapshot (which cleared the vars'
     * OwnValues) so they lower as the argument symbols, not folded constants.
     * The symbolic gradient is left exact — only the value path is compiled — so
     * FindMinimum's precision is unchanged. MPFR keeps the exact interpreter
     * path (its solvers never call fm_eval_scalar). A body Compile can't lower
     * stays NULL and the interpreter is used, and every per-point call falls
     * back on a non-finite compiled result, so this is a pure speedup.
     *
     * The exact symbolic gradient `g_exprs` is compiled the same way (each
     * component is a function of all the variables). This is what actually
     * accelerates the QuasiNewton/CG/Newton loop, whose cost is dominated by the
     * per-iteration gradient — the same ∂f/∂x_i, lowered, so the gradient stays
     * exact (no finite differences) and FindMinimum's precision is unchanged.
     * A component Compile can't lower stays NULL and falls back per-component.
     *
     * Skipped when an "EvaluationMonitor" is set: the monitor fires inside
     * fm_eval_with_bindings (per interpreter evaluation of f), which the compiled
     * path bypasses, so compiling would silently stop the monitor from firing.
     * Monitoring is a debugging aid, not a performance path, so falling back to
     * the interpreter there is the right trade. */
    if (opts.prec_mode == FM_PREC_MACHINE && n > 0 && !opts.eval_monitor) {
        const char** cnames = (const char**)malloc(sizeof(char*) * n);
        CompileType* ctypes = (CompileType*)malloc(sizeof(CompileType) * n);
        for (size_t i = 0; i < n; i++) { cnames[i] = vars[i]->data.symbol.name; ctypes[i] = CT_REAL; }
        f_prog = compile_expr_ex(f_raw, cnames, ctypes, n, COMPILE_FOLD_GLOBALS);
        if (f_prog && compiled_result_type(f_prog) != CT_REAL) { compiled_free(f_prog); f_prog = NULL; }
        if (g_exprs) {
            grad_progs = (CompiledProgram**)calloc(n, sizeof(CompiledProgram*));
            for (size_t i = 0; i < n; i++) {
                CompiledProgram* p = compile_expr_ex(g_exprs[i], cnames, ctypes, n,
                                                     COMPILE_FOLD_GLOBALS);
                if (p && compiled_result_type(p) != CT_REAL) { compiled_free(p); p = NULL; }
                grad_progs[i] = p;
            }
        }
        free(cnames); free(ctypes);
    }
    g_fm_obj_expr = f_raw; g_fm_obj_prog = f_prog; g_fm_obj_nargs = n;
    g_fm_grad_exprs = g_exprs; g_fm_grad_progs = grad_progs; g_fm_grad_n = n;

    /* Dispatch. */
    double fx_min = 0.0;
    bool ok = true;
    bool has_general_cons = (ngens > 0);
#ifdef USE_MPFR
    bool mpfr_result = false;
    mpfr_t* x_vec_mpfr = NULL;
    mpfr_t fx_min_mpfr;
    bool use_mpfr = (opts.prec_mode == FM_PREC_MPFR);
    /* Penalty path is not lifted to MPFR yet — fall back to machine
     * precision in that case rather than silently dropping the constraint. */
    if (use_mpfr && has_general_cons) {
        fm_warn(fn_name, "nimpl",
                "general (non-box) constraints at WorkingPrecision > MachinePrecision "
                "are not yet supported; falling back to machine precision");
        use_mpfr = false;
    }
    if (use_mpfr) {
        mpfr_init2(fx_min_mpfr, opts.wp_bits);
        x_vec_mpfr = fm_mpfr_array(n, opts.wp_bits);
        for (size_t i = 0; i < n; i++) mpfr_set_d(x_vec_mpfr[i], x_vec[i], MPFR_RNDN);
        if (method == FM_METHOD_BRENT) {
            if (n != 1) {
                fm_warn(fn_name, "badmeth", "Method \"Brent\" requires a single variable");
                ok = false; goto run_done;
            }
            mpfr_t a_m, b_m, c_m, fa_m, fb_m, fc_m;
            mpfr_init2(a_m, opts.wp_bits); mpfr_init2(b_m, opts.wp_bits);
            mpfr_init2(c_m, opts.wp_bits);
            mpfr_init2(fa_m, opts.wp_bits); mpfr_init2(fb_m, opts.wp_bits);
            mpfr_init2(fc_m, opts.wp_bits);
            bool bracketed = false;
            if (boxes[0].has_lo && boxes[0].has_hi) {
                mpfr_set_d(a_m, boxes[0].lo, MPFR_RNDN);
                mpfr_set_d(c_m, boxes[0].hi, MPFR_RNDN);
                mpfr_add(b_m, a_m, c_m, MPFR_RNDN); mpfr_div_ui(b_m, b_m, 2, MPFR_RNDN);
                /* If user gave a start inside, use it. */
                if (x_vec[0] > boxes[0].lo && x_vec[0] < boxes[0].hi)
                    mpfr_set(b_m, x_vec_mpfr[0], MPFR_RNDN);
                bracketed = true;
            } else {
                bracketed = fm_bracket_mpfr(f_raw, binds, &opts, x_vec_mpfr[0],
                                            &boxes[0], a_m, b_m, c_m, fa_m, fb_m, fc_m);
                if (!bracketed) fm_warn(fn_name, "nlnum", "MPFR bracket-finding failed");
            }
            if (bracketed) {
                mpfr_t xm_m, fmin_m;
                mpfr_init2(xm_m, opts.wp_bits); mpfr_init2(fmin_m, opts.wp_bits);
                ok = fm_brent_min_mpfr(f_raw, binds, &opts, a_m, b_m, c_m,
                                       &boxes[0], xm_m, fmin_m);
                if (ok) {
                    mpfr_set(x_vec_mpfr[0], xm_m, MPFR_RNDN);
                    mpfr_set(fx_min_mpfr, fmin_m, MPFR_RNDN);
                }
                mpfr_clears(xm_m, fmin_m, (mpfr_ptr)0);
            } else {
                ok = false;
            }
            mpfr_clears(a_m, b_m, c_m, fa_m, fb_m, fc_m, (mpfr_ptr)0);
        } else {
            /* n-D path: BFGS handles QuasiNewton; Newton/CG fall back to
             * BFGS at MPFR with a one-shot diagnostic. */
            if (method == FM_METHOD_NEWTON || method == FM_METHOD_CONJGRAD
                || method == FM_METHOD_LBFGSB || method == FM_METHOD_POWELL
                || method == FM_METHOD_NELDERMEAD || method == FM_METHOD_TNC
                || method == FM_METHOD_SLSQP || method == FM_METHOD_COBYLA
                || method == FM_METHOD_COBYQA || method == FM_METHOD_NEWTONCG
                || method == FM_METHOD_DOGLEG || method == FM_METHOD_TRUSTNCG
                || method == FM_METHOD_TRUSTEXACT || method == FM_METHOD_TRUSTKRYLOV) {
                const char* mname = method == FM_METHOD_NEWTON ? "Newton"
                                  : method == FM_METHOD_CONJGRAD ? "ConjugateGradient"
                                  : method == FM_METHOD_LBFGSB ? "LBFGSB"
                                  : method == FM_METHOD_POWELL ? "Powell"
                                  : method == FM_METHOD_NELDERMEAD ? "NelderMead"
                                  : method == FM_METHOD_TNC ? "TNC"
                                  : method == FM_METHOD_SLSQP ? "SLSQP"
                                  : method == FM_METHOD_COBYLA ? "COBYLA"
                                  : method == FM_METHOD_COBYQA ? "COBYQA"
                                  : method == FM_METHOD_NEWTONCG ? "NewtonCG"
                                  : method == FM_METHOD_DOGLEG ? "Dogleg"
                                  : method == FM_METHOD_TRUSTNCG ? "TrustNCG"
                                  : method == FM_METHOD_TRUSTEXACT ? "TrustExact"
                                  : "TrustKrylov";
                fm_warn(fn_name, "nimpl",
                        "Method \"%s\" at WorkingPrecision > MachinePrecision is not yet "
                        "supported; falling back to QuasiNewton", mname);
            }
            ok = fm_run_bfgs_mpfr(f_raw, vars, n, binds, g_exprs,
                                  x_vec_mpfr, boxes, &opts, fx_min_mpfr);
        }
        mpfr_result = ok;
    } else {
#endif
    if (method == FM_METHOD_BRENT) {
        if (n != 1) {
            fm_warn(fn_name, "badmeth", "Method \"Brent\" requires a single variable");
            goto cleanup;
        }
        if (has_general_cons) {
            fm_warn(fn_name, "nimpl", "general constraints with Brent are not supported");
            goto cleanup;
        }
        double a, b, c;
        if (boxes[0].has_lo && boxes[0].has_hi) {
            /* Use the supplied bounds; choose interior start. */
            a = boxes[0].lo; c = boxes[0].hi; b = (a + c) * 0.5;
            /* If user gave a starting point inside, use it. */
            if (x_vec[0] > a && x_vec[0] < c) b = x_vec[0];
        } else {
            if (!fm_bracket(f_raw, binds, &opts, x_vec[0], &boxes[0], &a, &b, &c)) {
                fm_warn(fn_name, "nlnum", "bracket-finding failed");
                ok = false; goto run_done;
            }
        }
        double xm, fm;
        ok = fm_brent_min(f_raw, binds, &opts, a, b, c, &boxes[0], &xm, &fm);
        if (ok) { x_vec[0] = xm; fx_min = fm; }
    } else if (method == FM_METHOD_QUASINEWTON) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_QUASINEWTON,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_bfgs(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_CONJGRAD) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_CONJGRAD,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_cg(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_NEWTON) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_NEWTON,
                                g_exprs, H_exprs, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_newton(f_raw, vars, n, binds, g_exprs, H_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_LBFGSB) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_LBFGSB,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_lbfgsb(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_POWELL) {
        /* Derivative-free; matches scipy's Powell, which supports box bounds
         * but no general nonlinear constraints -- reject those (mirror Brent)
         * so the comparison stays apples-to-apples. */
        if (has_general_cons) {
            fm_warn(fn_name, "nimpl",
                    "general (non-box) constraints are not supported with Method \"Powell\"");
            ok = false; goto cleanup;
        }
        ok = fm_run_powell(f_raw, vars, n, binds, NULL, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
    } else if (method == FM_METHOD_NELDERMEAD) {
        /* Derivative-free downhill simplex; like scipy's Nelder-Mead it supports
         * box bounds but no general nonlinear constraints -- reject those. */
        if (has_general_cons) {
            fm_warn(fn_name, "nimpl",
                    "general (non-box) constraints are not supported with Method \"NelderMead\"");
            ok = false; goto cleanup;
        }
        ok = fm_run_neldermead(f_raw, vars, n, binds, NULL, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
    } else if (method == FM_METHOD_TNC) {
        /* Hessian-free truncated Newton. Gradient-based like LBFGSB, so general
         * (non-box) constraints route through the augmented-Lagrangian wrapper
         * (consistent with the other gradient methods); box-only problems -- the
         * scipy-comparable case -- take the direct path. */
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_TNC,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_tnc(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_SLSQP) {
        /* Sequential least-squares QP.  Unlike the other gradient methods,
         * SLSQP consumes general (non-box) constraints DIRECTLY through its QP
         * subproblem rather than the augmented-Lagrangian penalty wrapper, so it
         * takes a single path for both the constrained and unconstrained cases
         * (fm_run_slsqp handles ngens == 0 as an ordinary empty active set). */
        ok = fm_run_slsqp(f_raw, vars, n, binds, g_exprs, x_vec,
                          gens, ngens, 0.0, boxes, &opts, &fx_min);
    } else if (method == FM_METHOD_COBYLA) {
        /* Derivative-free constrained (Powell's linear-approximation trust
         * region).  Like SLSQP it consumes general constraints DIRECTLY (unlike
         * Powell/NelderMead, which reject them), so a single path serves both the
         * constrained and unconstrained cases. */
        ok = fm_run_cobyla(f_raw, vars, n, binds, NULL, x_vec,
                           gens, ngens, 0.0, boxes, &opts, &fx_min);
    } else if (method == FM_METHOD_COBYQA) {
        /* Derivative-free constrained with quadratic models; direct gens path
         * (native equality + inequality + bounds), like COBYLA/SLSQP. */
        ok = fm_run_cobyqa(f_raw, vars, n, binds, NULL, x_vec,
                           gens, ngens, 0.0, boxes, &opts, &fx_min);
    } else if (method == FM_METHOD_NEWTONCG || method == FM_METHOD_DOGLEG
            || method == FM_METHOD_TRUSTNCG || method == FM_METHOD_TRUSTEXACT
            || method == FM_METHOD_TRUSTKRYLOV) {
        /* Trust-region / truncated-Newton family. Unconstrained in scipy, so
         * reject general (non-box) constraints like Powell/NelderMead; the
         * runners then take a single ngens == 0 path. dogleg/trust-exact form
         * the dense Hessian (H_exprs), the others are Hessian-free. */
        const char* mn = method == FM_METHOD_NEWTONCG ? "NewtonCG"
                       : method == FM_METHOD_DOGLEG   ? "Dogleg"
                       : method == FM_METHOD_TRUSTNCG ? "TrustNCG"
                       : method == FM_METHOD_TRUSTEXACT ? "TrustExact"
                       : "TrustKrylov";
        if (has_general_cons) {
            fm_warn(fn_name, "nimpl",
                    "general (non-box) constraints are not supported with Method \"%s\"", mn);
            ok = false; goto cleanup;
        }
        /* scipy's trust-region methods cannot handle box bounds either — warn
         * and solve unconstrained rather than silently dropping them. */
        bool any_box = false;
        for (size_t bi = 0; bi < n; bi++)
            if (boxes[bi].has_lo || boxes[bi].has_hi) { any_box = true; break; }
        if (any_box)
            fm_warn(fn_name, "nimpl",
                    "bounds are not supported with Method \"%s\"; solving unconstrained", mn);
        if (method == FM_METHOD_NEWTONCG) {
            ok = fm_run_newton_cg(f_raw, vars, n, binds, g_exprs, x_vec,
                                  NULL, 0, 0.0, boxes, &opts, &fx_min);
        } else if (method == FM_METHOD_DOGLEG) {
            ok = fm_run_trust_region(f_raw, vars, n, binds, g_exprs, H_exprs,
                                     x_vec, boxes, &opts, &fx_min, fm_tr_dogleg, true);
        } else if (method == FM_METHOD_TRUSTNCG) {
            ok = fm_run_trust_region(f_raw, vars, n, binds, g_exprs, NULL,
                                     x_vec, boxes, &opts, &fx_min, fm_tr_steihaug, false);
        } else if (method == FM_METHOD_TRUSTEXACT) {
            ok = fm_run_trust_region(f_raw, vars, n, binds, g_exprs, H_exprs,
                                     x_vec, boxes, &opts, &fx_min, fm_tr_moresorensen, true);
        } else {
            ok = fm_run_trust_region(f_raw, vars, n, binds, g_exprs, NULL,
                                     x_vec, boxes, &opts, &fx_min, fm_tr_gltr, false);
        }
    } else {
        fm_warn(fn_name, "nimpl", "method not implemented");
        ok = false;
    }
#ifdef USE_MPFR
    }
#endif
run_done:
    /* Clear temp bindings first so the variable symbol stays free during
     * Rule construction (otherwise `Rule[x, v]` would re-evaluate x to its
     * pre-call value once we restore). */
    if (binds) {
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    }
    if (ok) {
#ifdef USE_MPFR
        if (mpfr_result) {
            result_out = fm_build_result_mpfr(fx_min_mpfr, vars,
                                              (mpfr_t const*)x_vec_mpfr, n);
        } else {
            result_out = fm_build_result(fx_min, vars, x_vec, n);
        }
#else
        result_out = fm_build_result(fx_min, vars, x_vec, n);
#endif
    }
#ifdef USE_MPFR
    if (x_vec_mpfr) fm_mpfr_array_free(x_vec_mpfr, n);
    if (use_mpfr)   mpfr_clear(fx_min_mpfr);
#endif

cleanup:
    g_fm_obj_expr = NULL;      /* deregister objective + gradient before freeing */
    g_fm_obj_prog = NULL;
    g_fm_obj_nargs = 0;
    g_fm_grad_exprs = NULL;
    g_fm_grad_progs = NULL;
    g_fm_grad_n = 0;
    if (f_prog) compiled_free(f_prog);
    if (grad_progs) {
        for (size_t i = 0; i < n; i++) if (grad_progs[i]) compiled_free(grad_progs[i]);
        free(grad_progs);
    }
    if (binds) {
        for (size_t i = 0; i < n; i++) fm_bind_restore(&binds[i]);
        free(binds);
    }
    if (g_exprs) { for (size_t i = 0; i < n; i++) expr_free(g_exprs[i]); free(g_exprs); }
    if (H_exprs) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) expr_free(H_exprs[i][j]);
            free(H_exprs[i]);
        }
        free(H_exprs);
    }
    if (gens) {
        for (size_t k = 0; k < ngens; k++) {
            expr_free(gens[k].expr);
            if (gens[k].grad_exprs) {
                for (size_t i = 0; i < n; i++) expr_free(gens[k].grad_exprs[i]);
                free(gens[k].grad_exprs);
            }
        }
        free(gens);
    }
    free(vars);
    free(x_vec);
    free(boxes);
    return result_out;
}

Expr* builtin_findminimum(Expr* res) {
    return findmin_driver(res, "FindMinimum");
}

/* FindMaximum: build FindMinimum[-f, vars, opts...] internally. */
Expr* builtin_findmaximum(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn("FindMaximum", "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Negate the objective (or the f inside {f, cons}). */
    Expr* f_orig = res->data.function.args[0];
    Expr* neg_f;
    Expr* new_first;
    if (f_orig->type == EXPR_FUNCTION
        && f_orig->data.function.head->type == EXPR_SYMBOL
        && f_orig->data.function.head->data.symbol.name == SYM_List
        && f_orig->data.function.arg_count == 2) {
        /* Wrap inner f only. */
        Expr* inner_f = f_orig->data.function.args[0];
        Expr* cons = f_orig->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(inner_f) };
        neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        Expr* list_args[2] = { neg_f, expr_copy(cons) };
        new_first = expr_new_function(expr_new_symbol(SYM_List), list_args, 2);
    } else {
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(f_orig) };
        neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        new_first = neg_f;
    }
    /* Construct synthetic FindMinimum[new_first, vars, opts...]. */
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * argc);
    new_args[0] = new_first;
    for (size_t i = 1; i < argc; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    Expr* synthetic = expr_new_function(expr_new_symbol(SYM_FindMinimum), new_args, argc);
    free(new_args);
    /* Drive findmin directly so the diagnostic tag is FindMaximum. */
    Expr* min_result = findmin_driver(synthetic, "FindMaximum");
    expr_free(synthetic);
    if (!min_result) return NULL;
    /* min_result is {fmin, {rules}}; negate fmin while preserving its
     * numeric type so a WorkingPrecision -> N run keeps the N-digit MPFR
     * head instead of collapsing back to machine precision. */
    if (min_result->type == EXPR_FUNCTION
        && min_result->data.function.arg_count == 2) {
        Expr* fmin_e = min_result->data.function.args[0];
#ifdef USE_MPFR
        if (fmin_e && fmin_e->type == EXPR_MPFR) {
            long bits = mpfr_get_prec(fmin_e->data.mpfr);
            mpfr_t neg; mpfr_init2(neg, bits);
            mpfr_neg(neg, fmin_e->data.mpfr, MPFR_RNDN);
            expr_free(fmin_e);
            min_result->data.function.args[0] = expr_new_mpfr_copy(neg);
            mpfr_clear(neg);
        } else
#endif
        {
            double fmin;
            if (fm_expr_to_double_real(fmin_e, &fmin)) {
                expr_free(fmin_e);
                min_result->data.function.args[0] = expr_new_real(-fmin);
            }
        }
    }
    return min_result;
}
