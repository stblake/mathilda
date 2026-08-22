/* nm_driver.c — NMinimize / NMaximize driver + builtins.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* nm_minimize_driver(Expr* res, const char* fn_name) {
    g_fm_name = fn_name;
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn(fn_name, "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Peel trailing options (same recogniser FindMinimum uses). */
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
    opts.max_iter = 100;                       /* NMinimize default */
    opts.max_iter_set = false;
    opts.acc_goal_digits = -1.0;
    opts.prec_goal_digits = -1.0;
    opts.gradient = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;

    NmConfig nc;
    nc.method = NM_AUTO;
    nc.search_points = 0;
    nc.F = -1.0;
    nc.CR = -1.0;
    nc.reflect_ratio = -1.0;
    nc.expand_ratio = -1.0;
    nc.contract_ratio = -1.0;
    nc.shrink_ratio = -1.0;
    nc.tolerance = -1.0;
    nc.post_process = -1;
    nc.init_points = NULL;
    nc.penalty_fn = NULL;
    nc.perturb_scale = -1.0;
    nc.boltzmann_fn = NULL;
    nc.level_iterations = 0;
    nc.shgo_sampling = 0;
    nc.shgo_iters = 0;
    nc.da_visit = NM_DA_VISIT;
    nc.da_accept = NM_DA_ACCEPT;
    nc.da_init_temp = NM_DA_INIT_TEMP;
    nc.da_restart_ratio = NM_DA_RESTART_RATIO;
    nc.da_local_search = -1;
    nc.direct_locally_biased = -1;
    nc.direct_eps = -1.0;
    nc.direct_max_fun = 0;
    nc.direct_max_iter = 0;
    nc.direct_vol_tol = -1.0;
    nc.direct_len_tol = -1.0;
    nc.direct_fmin = -HUGE_VAL;
    nc.direct_fmin_rtol = -1.0;
    nc.bh_temp = NM_BH_TEMP;
    nc.bh_step = NM_BH_STEP;
    nc.bh_interval = NM_BH_INTERVAL;
    nc.bh_target_accept = NM_BH_TARGET_ACCEPT;
    nc.bh_step_factor = NM_BH_STEP_FACTOR;
    nc.bh_niter_success = 0;
    nc.seed = NM_DEFAULT_SEED;

    for (size_t i = pos_end; i < argc; i++) {
        if (!nm_apply_option(res->data.function.args[i], &opts, &nc, fn_name))
            return NULL;
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

    /* Parse the variables first — they do not depend on the objective — so a
     * problem passed as a bare symbol (prob = {f, cons}; NMinimize[prob, vars])
     * can be resolved below with those variables protected. */
    Expr* var_arg = res->data.function.args[1];

    /* Parse the variable specification. NMinimize is not HoldAll, so a generator
     * such as Table[x[i], {i, 1, 10}] or Array[x, 10] has already expanded to the
     * variable list before we run, and a symbol bound to a list has already
     * resolved. The eval_and_free calls below are therefore idempotent
     * normalizers (they no-op on an already-evaluated {...}/Element, and still
     * cover a spec reaching the driver unevaluated). A {...} list or Element[...]
     * is used directly. A bare symbol is ambiguous: it may be a single
     * optimization variable (NMinimize[f, x], x unbound) or a symbol that still
     * resolves to a variable list — evaluate it, and if it yields a List/Element
     * use that, otherwise keep the symbol itself as the single variable. */
    NmVarSet vs;
    Expr* var_list_eval = NULL;
    {
        Expr* var_spec = var_arg;
        if (var_arg->type == EXPR_SYMBOL) {
            Expr* ev = eval_and_free(expr_copy(var_arg));
            if (ev && (nm_is_head(ev, SYM_List) || nm_is_head(ev, SYM_Element))) {
                var_list_eval = ev;
                var_spec = ev;
            } else {
                expr_free(ev);          /* unbound / scalar: use the symbol */
            }
        } else if (!nm_is_head(var_arg, SYM_List) && !nm_is_head(var_arg, SYM_Element)) {
            var_list_eval = eval_and_free(expr_copy(var_arg));
            var_spec = var_list_eval;
        }
        if (!var_spec || !nm_parse_vars(var_spec, &vs, fn_name)) {
            expr_free(var_list_eval);
            return NULL;
        }
    }
    size_t n = vs.n;

    /* Effective (scalar-symbol) variables for the solver machinery, and the
     * original variable expressions for the result rules. Indexed variables
     * (x[1], x[2], ...) are rewritten to fresh scalar symbols so the entire
     * symbol-keyed solver applies unchanged; plain symbols pass through. */
    bool indexed = false;
    for (size_t i = 0; i < n; i++)
        if (vs.vars[i]->type != EXPR_SYMBOL) indexed = true;

    Expr**       eff_vars  = (Expr**)calloc(n, sizeof(Expr*));
    Expr**       orig_vars = (Expr**)calloc(n, sizeof(Expr*));
    const char** synth     = (const char**)calloc(n, sizeof(char*));
    const char** heads      = (const char**)malloc(sizeof(char*) * (n ? n : 1));
    size_t       nheads     = 0;
    for (size_t i = 0; i < n; i++) {
        orig_vars[i] = expr_copy(vs.vars[i]);
        if (vs.vars[i]->type == EXPR_SYMBOL) {
            eff_vars[i] = expr_copy(vs.vars[i]);
        } else {
            eff_vars[i] = nm_fresh_symbol();
            synth[i]    = eff_vars[i]->data.symbol.name;
        }
        const char* hn = (vs.vars[i]->type == EXPR_SYMBOL)
            ? vs.vars[i]->data.symbol.name
            : vs.vars[i]->data.function.head->data.symbol.name;
        bool seen = false;
        for (size_t j = 0; j < nheads; j++) if (heads[j] == hn) { seen = true; break; }
        if (!seen) heads[nheads++] = hn;
    }

    /* Resolve and split {f, cons}. NMinimize is not HoldAll, so the problem
     * argument normally arrives already in structural form (an inline {f, cons}
     * list or a scalar objective) and is used directly. The bare-symbol branch
     * below is a defensive resolver: if the argument still reaches the driver as
     * a symbol, it is evaluated once with the variable heads localized — so its
     * structure is exposed without capturing any global variable values — and the
     * {f, cons} list is then split. */
    Expr* f_arg  = res->data.function.args[0];
    Expr* f_eval = NULL;               /* owned resolved objective, or NULL     */
    if (f_arg->type == EXPR_SYMBOL) {
        NmHeadSave* hs = (NmHeadSave*)calloc(nheads ? nheads : 1, sizeof(NmHeadSave));
        nm_heads_localize(hs, heads, nheads);
        f_eval = eval_and_free(expr_copy(f_arg));
        nm_heads_restore(hs, nheads);
        free(hs);
        if (f_eval) f_arg = f_eval;
    }
    Expr* f_raw = f_arg;
    Expr* cons = NULL;         /* borrowed, or points at cons_built           */
    Expr* cons_built = NULL;   /* owned And[...] for a multi-constraint list  */
    if (nm_is_head(f_arg, SYM_List) && f_arg->data.function.arg_count >= 2) {
        /* {f, cons} or {f, c1, c2, ...}: the objective is the first element
         * and every remaining element is a constraint, implicitly And-ed. */
        f_raw = f_arg->data.function.args[0];
        size_t ncons = f_arg->data.function.arg_count - 1;
        if (ncons == 1) {
            cons = f_arg->data.function.args[1];
        } else {
            Expr** cc = (Expr**)malloc(sizeof(Expr*) * ncons);
            for (size_t i = 0; i < ncons; i++)
                cc[i] = expr_copy(f_arg->data.function.args[1 + i]);
            cons_built = expr_new_function(expr_new_symbol(SYM_And), cc, ncons);
            free(cc);
            cons = cons_built;
        }
    }

    /* Expand a held Table/Sum constraint and/or rewrite indexed variables. The
     * objective stays held for the plain path (evaluated per point); it is
     * pre-expanded only when we must rewrite indexed vars inside it. */
    Expr* f_eff = f_raw;   bool f_owned = false;
    Expr* cons_eff = cons; bool cons_owned = false;
    bool infeasible_pre = false;
    bool expand_cons = (cons != NULL) && !nm_is_constraint_tree(cons);
    if (indexed || expand_cons) {
        NmHeadSave* hs = (NmHeadSave*)calloc(nheads ? nheads : 1, sizeof(NmHeadSave));
        nm_heads_localize(hs, heads, nheads);

        if (cons) {
            Expr* ce = eval_and_free(expr_copy(cons));   /* expand, vars free */
            if (ce && ce->type == EXPR_SYMBOL && ce->data.symbol.name == SYM_True) {
                expr_free(ce); ce = NULL;
            } else if (ce && ce->type == EXPR_SYMBOL && ce->data.symbol.name == SYM_False) {
                expr_free(ce); ce = NULL; infeasible_pre = true;
            } else if (nm_is_head(ce, SYM_List)) {
                /* A list of constraints is an implicit And; drop trivially-true
                 * entries, and a False entry makes the whole system infeasible. */
                size_t m = ce->data.function.arg_count;
                Expr** cc = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
                size_t kept = 0;
                for (size_t i = 0; i < m; i++) {
                    Expr* el = ce->data.function.args[i];
                    if (el->type == EXPR_SYMBOL && el->data.symbol.name == SYM_True) continue;
                    if (el->type == EXPR_SYMBOL && el->data.symbol.name == SYM_False) {
                        infeasible_pre = true; continue;
                    }
                    cc[kept++] = expr_copy(el);
                }
                expr_free(ce);
                if (kept == 0)      { ce = NULL; }
                else if (kept == 1) { ce = cc[0]; }
                else                { ce = expr_new_function(expr_new_symbol(SYM_And), cc, kept); }
                free(cc);
            }
            if (ce && indexed) {
                Expr* cs = nm_subst(ce, vs.vars, eff_vars, n);
                expr_free(ce); ce = cs;
            }
            cons_eff = ce; cons_owned = true;
        }

        if (indexed) {
            Expr* fe = eval_and_free(expr_copy(f_raw));
            f_eff = nm_subst(fe, vs.vars, eff_vars, n);
            expr_free(fe);
            f_owned = true;
        }

        nm_heads_restore(hs, nheads);
        free(hs);
    }

    /* Bind variables (Block semantics). */
    FmVarBind* binds = (FmVarBind*)calloc(n, sizeof(FmVarBind));
    for (size_t i = 0; i < n; i++)
        fm_bind_snapshot(&binds[i], eff_vars[i]->data.symbol.name);

    FmBox* boxes = (FmBox*)calloc(n, sizeof(FmBox));
    FmGenCon* gens = NULL;
    size_t ngens = 0, gcap = 0;
    FmDisjunction* disj = NULL;
    size_t ndisj = 0, dcap = 0;
    Expr** g_exprs = NULL;
    Expr* cons2 = NULL;
    double* reg_lo = NULL;
    double* reg_hi = NULL;
    bool*   used_default = NULL;   /* dim used the default +-SPAN (fully unbounded) */
    double* xbest = NULL;
    Expr* result_out = NULL;
    CompiledProgram*  f_prog = NULL;   /* compiled objective (machine prec)     */
    CompiledProgram** g_progs = NULL;  /* compiled general constraints          */

    /* Declared here (before the first `goto cleanup`) so the cleanup path always
     * sees an initialized one-hot list — the remaining fields are filled below on
     * the non-error path. */
    NmDriver D;
    D.onehots = NULL; D.n_onehots = 0;

    /* Extract integer/real domain declarations, then collect the remaining
     * constraints into boxes + general FmGenCon[] + disjunctions. */
    if (cons_eff) {
        cons2 = nm_filter_int(cons_eff, eff_vars, n, vs.is_int);
        for (size_t i = 0; i < n; i++) if (vs.is_int[i]) vs.any_int = true;
        if (cons2 && !fm_collect_constraints(cons2, eff_vars, n, boxes,
                                             &gens, &ngens, &gcap,
                                             &disj, &ndisj, &dcap))
            goto cleanup;
        for (size_t k = 0; k < ngens; k++)
            gens[k].grad_exprs = fm_compute_gradient(gens[k].expr, eff_vars, n);
    }

    /* Resolve the per-dimension search box: box constraints tighten it,
     * else a starting-interval hint, else a default span. Contradictory box
     * bounds (lo > hi, e.g. x > 2 && x < 1) mean an empty feasible set. */
    bool infeasible_box = false;
    reg_lo = (double*)malloc(sizeof(double) * n);
    reg_hi = (double*)malloc(sizeof(double) * n);
    used_default = (bool*)calloc(n ? n : 1, sizeof(bool));
    bool any_default = false;
    for (size_t i = 0; i < n; i++) {
        bool klo = boxes[i].has_lo || vs.has_rlo[i];
        bool khi = boxes[i].has_hi || vs.has_rhi[i];
        double lo = boxes[i].has_lo ? boxes[i].lo : vs.rlo[i];
        double hi = boxes[i].has_hi ? boxes[i].hi : vs.rhi[i];
        if (boxes[i].has_lo && boxes[i].has_hi && boxes[i].lo > boxes[i].hi)
            infeasible_box = true;
        if (klo && khi)      { reg_lo[i] = lo; reg_hi[i] = hi; }
        else if (klo)        { reg_lo[i] = lo; reg_hi[i] = lo + NM_BOUND_SPAN; }
        else if (khi)        { reg_hi[i] = hi; reg_lo[i] = hi - NM_BOUND_SPAN; }
        else                 { reg_lo[i] = -NM_DEFAULT_SPAN; reg_hi[i] = NM_DEFAULT_SPAN;
                               used_default[i] = true; any_default = true; }
        if (reg_hi[i] <= reg_lo[i]) {
            double m = 0.5 * (reg_lo[i] + reg_hi[i]);
            reg_lo[i] = m - 0.5; reg_hi[i] = m + 0.5;
        }
        if (vs.is_int[i]) { reg_lo[i] = floor(reg_lo[i]); reg_hi[i] = ceil(reg_hi[i]); }
    }

    /* Objective gradient (for the continuous local polish; NULL → FD). */
    if (!vs.any_int) g_exprs = fm_compute_gradient(f_eff, eff_vars, n);

    /* Machine-precision auto-compilation: the global search evaluates the
     * objective (and each general constraint) at hundreds–thousands of trial
     * points, so lowering them to bytecode over the effective variables once and
     * running the register machine per point is far cheaper than the interpreter
     * (expr_copy + evaluate + numericalize each call). The variables are already
     * unbound (Block snapshot cleared their OwnValues), so they compile as the
     * argument symbols; COMPILE_FOLD_GLOBALS folds any other machine-valued
     * symbol — safe because these programs live only for this call. A body with
     * a construct Compile can't lower stays NULL and uses the interpreter, and
     * every per-point call falls back to the interpreter on a domain/non-finite
     * result, so this is a pure speedup with no change in answer. MPFR
     * (WorkingPrecision > MachinePrecision) keeps the exact interpreter path. */
    if (opts.prec_mode == FM_PREC_MACHINE) {
        const char** cnames = (const char**)malloc(sizeof(char*) * (n ? n : 1));
        CompileType* ctypes = (CompileType*)malloc(sizeof(CompileType) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            cnames[i] = eff_vars[i]->data.symbol.name;
            ctypes[i] = CT_REAL;
        }
        f_prog = compile_expr_ex(f_eff, cnames, ctypes, n, COMPILE_FOLD_GLOBALS);
        if (f_prog && compiled_result_type(f_prog) != CT_REAL) {
            compiled_free(f_prog); f_prog = NULL;
        }
        if (ngens > 0) {
            g_progs = (CompiledProgram**)calloc(ngens, sizeof(CompiledProgram*));
            for (size_t k = 0; k < ngens; k++) {
                CompiledProgram* p = compile_expr_ex(gens[k].expr, cnames, ctypes,
                                                     n, COMPILE_FOLD_GLOBALS);
                if (p && compiled_result_type(p) != CT_REAL) { compiled_free(p); p = NULL; }
                g_progs[k] = p;
            }
        }
        free(cnames); free(ctypes);
    }

    D.f_raw = f_eff; D.vars = eff_vars; D.n = n; D.binds = binds;
    D.g_exprs = g_exprs; D.gens = gens; D.ngens = ngens; D.boxes = boxes;
    D.opts = &opts; D.is_int = vs.is_int; D.any_int = vs.any_int;
    /* Does any continuous variable appear in the objective? (Gates the heavy
     * continuous-relaxation recovery in nm_local_polish — see that flag's note.) */
    D.cont_in_obj = false;
    for (size_t i = 0; i < n; i++)
        if (!vs.is_int[i] && nm_expr_contains_symbol(f_eff, eff_vars[i]->data.symbol.name)) {
            D.cont_in_obj = true;
            break;
        }
    D.reg_lo = reg_lo; D.reg_hi = reg_hi;
    D.f_prog = f_prog; D.g_progs = g_progs;
    D.penalty_fn = nc.penalty_fn;
    D.disj = disj; D.ndisj = ndisj;
    nm_detect_onehots(&D);   /* assignment groups for the integer repair */

    /* Serve the local polish's objective evaluations from the compiled program
     * instead of the interpreter. RandomSearch runs one local solve per
     * SearchPoint, so this is the difference between 1000 compiled polishes and
     * 1000 interpreter polishes (see g_fm_obj_* and fm_eval_scalar). Cleared to
     * the inactive sentinel below, before f_prog is freed. A plain reset is
     * enough rather than save/restore: the registration is live only while
     * f_prog != NULL, and an objective compilable enough to have a program
     * cannot contain a nested NMinimize call (NMinimize is not a compilable
     * head), so no active registration can ever be clobbered by re-entry. */
    g_fm_obj_expr  = f_eff;
    g_fm_obj_prog  = f_prog;
    g_fm_obj_nargs = n;

    xbest = (double*)malloc(sizeof(double) * n);
    double fbest = 1e300, penbest = 1e300;
    int method = (nc.method == NM_AUTO) ? NM_DE : nc.method;
    bool do_post = (nc.post_process != 0);

    /* Adaptive search-region expansion. If the default +-SPAN sampling region
     * contains no feasible point, the fully-unbounded coordinates are grown by
     * successive powers of ten and the search is retried — so a feasible region
     * whose location is implied by nonlinear constraints rather than stated as
     * variable bounds (e.g. the pressure-vessel MINLP, feasible near x3 ~ 52) is
     * still found instead of reporting {Infinity, ...}. Only fully-unbounded
     * coordinates grow; a coordinate carrying a box bound or a starting-interval
     * hint keeps its resolved region. Attempt 0 is the base region with the base
     * seed, so a problem already feasible there is solved identically to before;
     * expansion triggers only to rescue infeasibility, and stops as soon as a
     * feasible point is found (the smallest region that yields feasibility,
     * which keeps the search from drifting into far, non-physical basins). */
    int max_attempt = (any_default && !infeasible_box && !infeasible_pre)
                    ? NM_MAX_REGION_EXPAND : 0;
    double* xattempt = (double*)malloc(sizeof(double) * n);
    for (int attempt = 0; attempt <= max_attempt; attempt++) {
        if (attempt > 0) {
            double span = NM_DEFAULT_SPAN * pow(10.0, (double)attempt);
            for (size_t i = 0; i < n; i++) if (used_default[i]) {
                reg_lo[i] = vs.is_int[i] ? floor(-span) : -span;
                reg_hi[i] = vs.is_int[i] ? ceil(span)   :  span;
            }
        }
        double fa = 1e300, pa = 1e300;
        NmRng rng;
        nm_rng_seed(&rng, nc.seed + (uint64_t)attempt * 0x100000001B3ULL);
        switch (method) {
            case NM_NELDERMEAD:   nm_neldermead(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_RANDOMSEARCH: nm_randomsearch(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_SA:           nm_sa(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_SHGO:         nm_shgo(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_DUAL_ANNEALING: nm_dual_annealing(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_DIRECT:       nm_direct(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_BASIN_HOPPING: nm_basin_hopping(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_DE:
            default:              nm_de(&D, &nc, &rng, xattempt, &fa, &pa); break;
        }
        /* Polish this attempt's best with the exact local solver, unless the
         * caller disabled it with "PostProcess" -> False. Guard against a
         * penalty/BFGS step that overshoots: if the polished point is worse by
         * Deb's rules than the pre-polish best, keep the latter. */
        if (do_post) {
            double* xsave = (double*)malloc(sizeof(double) * n);
            for (size_t i = 0; i < n; i++) xsave[i] = xattempt[i];
            double fsave = fa, psave = pa;
            nm_local_polish(&D, xattempt, &fa, &pa);
            /* TIGHT: post-polish overshoot guard is a return-path comparison. */
            if (nm_better_return(fsave, psave, fa, pa)) {
                for (size_t i = 0; i < n; i++) xattempt[i] = xsave[i];
                fa = fsave; pa = psave;
            }
            free(xsave);
        }
        /* TIGHT: penbest is what the feasibility gate below judges, so the
         * point that reaches it must be chosen under the return threshold. */
        if (attempt == 0 || nm_better_return(fa, pa, fbest, penbest)) {
            for (size_t i = 0; i < n; i++) xbest[i] = xattempt[i];
            fbest = fa; penbest = pa;
        }
        if (penbest <= NM_FEAS_RETURN) break;  /* genuinely feasible: stop expanding */
    }
    free(xattempt);
    /* The guarantee. A point that cannot meet NM_FEAS_RETURN is NOT handed back
     * as a solution — it falls through to nm_build_infeasible and is reported as
     * {Infinity, x -> Indeterminate}. Returning a constraint-violating point and
     * calling it feasible is the original bug; reporting a worse-than-tolerance
     * result honestly is the fix. */
    bool feasible = !infeasible_box && !infeasible_pre && (penbest <= NM_FEAS_RETURN);

    /* Optional MPFR refinement for WorkingPrecision > MachinePrecision on
     * continuous, general-constraint-free problems (reuses fm_run_bfgs_mpfr). */
#ifdef USE_MPFR
    bool mpfr_built = false;
    bool want_mpfr = (opts.prec_mode == FM_PREC_MPFR);
    bool mpfr_eligible = want_mpfr && do_post && !vs.any_int && ngens == 0;
    mpfr_t* xm = NULL;
    mpfr_t fmv;
    if (want_mpfr && !mpfr_eligible && do_post)
        fm_warn(fn_name, "nimpl",
                "WorkingPrecision > MachinePrecision with general constraints or "
                "integer domains is not supported; using machine precision");
    if (feasible && mpfr_eligible) {
        Expr** gm = fm_compute_gradient(f_eff, eff_vars, n);
        xm = fm_mpfr_array(n, opts.wp_bits);
        for (size_t i = 0; i < n; i++) mpfr_set_d(xm[i], xbest[i], MPFR_RNDN);
        mpfr_init2(fmv, opts.wp_bits);
        bool saved_quiet = g_fm_quiet;
        g_fm_quiet = true;
        bool mok = fm_run_bfgs_mpfr(f_eff, eff_vars, n, binds, gm, xm, boxes, &opts, fmv);
        g_fm_quiet = saved_quiet;
        if (mok)
            mpfr_built = true;
        else { fm_mpfr_array_free(xm, n); xm = NULL; mpfr_clear(fmv); }
        if (gm) { for (size_t i = 0; i < n; i++) expr_free(gm[i]); free(gm); }
    }
#endif

    /* Free the temporary bindings so the variable symbols are unbound while
     * we build the result rules (Rule[x, v] must not re-evaluate x). */
    for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);

#ifdef USE_MPFR
    if (mpfr_built) {
        result_out = fm_build_result_mpfr(fmv, orig_vars, (const mpfr_t*)xm, n);
        fm_mpfr_array_free(xm, n);
        mpfr_clear(fmv);
    } else
#endif
    if (feasible) result_out = nm_build_result(fbest, orig_vars, xbest, vs.is_int, n);
    else          result_out = nm_build_infeasible(orig_vars, n);

cleanup:
    for (size_t i = 0; i < n; i++) fm_bind_restore(&binds[i]);
    free(binds);
    if (g_exprs) { for (size_t i = 0; i < n; i++) expr_free(g_exprs[i]); free(g_exprs); }
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
    if (disj) {
        for (size_t k = 0; k < ndisj; k++) expr_free(disj[k].expr);
        free(disj);
    }
    expr_free(cons2);
    if (cons_owned) expr_free(cons_eff);
    if (f_owned)    expr_free(f_eff);
    expr_free(cons_built);
    expr_free(f_eval);         /* resolved bare-symbol objective; frees f_raw/cons borrows */
    g_fm_obj_expr  = NULL;     /* deregister the objective before f_prog is freed */
    g_fm_obj_prog  = NULL;
    g_fm_obj_nargs = 0;
    if (f_prog) compiled_free(f_prog);
    if (g_progs) {
        for (size_t k = 0; k < ngens; k++) if (g_progs[k]) compiled_free(g_progs[k]);
        free(g_progs);
    }
    nm_free_onehots(&D);
    free(boxes);
    free(reg_lo);
    free(reg_hi);
    free(used_default);
    free(xbest);
    if (eff_vars)  { for (size_t i = 0; i < n; i++) expr_free(eff_vars[i]);  free(eff_vars); }
    if (orig_vars) { for (size_t i = 0; i < n; i++) expr_free(orig_vars[i]); free(orig_vars); }
    if (synth) {
        for (size_t i = 0; i < n; i++) if (synth[i]) symtab_remove_symbol(synth[i]);
        free(synth);
    }
    free(heads);
    expr_free(var_list_eval);
    nm_varset_free(&vs);
    return result_out;
}

Expr* builtin_nminimize(Expr* res) {
    return nm_minimize_driver(res, "NMinimize");
}

/* NMaximize: minimise −f and negate the reported optimum. Mirrors the
 * FindMaximum → FindMinimum wrapper above. */
Expr* builtin_nmaximize(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn("NMaximize", "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    Expr* f_orig = res->data.function.args[0];
    Expr* new_first;
    if (nm_is_head(f_orig, SYM_List) && f_orig->data.function.arg_count == 2) {
        Expr* inner_f = f_orig->data.function.args[0];
        Expr* cons = f_orig->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(inner_f) };
        Expr* neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        Expr* list_args[2] = { neg_f, expr_copy(cons) };
        new_first = expr_new_function(expr_new_symbol(SYM_List), list_args, 2);
    } else {
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(f_orig) };
        new_first = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    }
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * argc);
    new_args[0] = new_first;
    for (size_t i = 1; i < argc; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    Expr* synthetic = expr_new_function(expr_new_symbol(SYM_NMinimize), new_args, argc);
    free(new_args);
    Expr* min_result = nm_minimize_driver(synthetic, "NMaximize");
    expr_free(synthetic);
    if (!min_result) return NULL;
    /* Negate the reported optimum value while preserving its numeric type. */
    if (min_result->type == EXPR_FUNCTION && min_result->data.function.arg_count == 2) {
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
