/*
 * dsolve_nonhomog_vop.c — DSolve`VariationOfParameters (M39).
 *
 * General nonhomogeneous linear-ODE backstop.  For a linear equation
 *     L[y] = c_n(x) y^(n) + ... + c_1(x) y' + c_0(x) y == g(x),   g ≢ 0,
 * whose HOMOGENEOUS part L[y]==0 is solvable by a cascade method that does not
 * itself carry forcing — the Bessel/Airy special forms (`SpecialFunctionForm`)
 * and the change-of-variable family (`ChangeOfVariable`) — this method:
 *   1. recurses `DSolve` on the homogeneous equation to obtain the general
 *      solution C[1] y1 + ... + C[n] yn,
 *   2. extracts and normalises the fundamental set y_k = PowerExpand[Simplify[
 *      D[hom, C[k]]]] (Simplify canonicalises 1+Tan^2 -> Sec^2, ...; PowerExpand
 *      then combines the radicals — e.g. Sqrt[1/(Pi x)]*x -> Sqrt[x]/Sqrt[Pi],
 *      1/Sqrt[Sec^2 x] -> Cos x — so the variation-of-parameters integrals close
 *      instead of staying inert),
 *   3. builds the particular by variation of parameters
 *      (`dsolve_variation_of_parameters`), and returns hom + yp.
 *
 * The constant-coefficient / Euler / undetermined-coefficient forcings are
 * claimed by their own methods earlier in the cascade, so this fires only as a
 * backstop for the variable-coefficient special-function homogeneous set.
 * Placed BEFORE `Kovacic` so the recursive homogeneous solve (which reaches
 * `SpecialFunctionForm` for a Bessel/Airy set) closes it before Kovacic — whose
 * constant-r nonhomogeneous path can churn — is reached on the FULL equation.
 *
 * Bounded exactly as the M12/M14 methods (TimeConstrained sub-solve + wall-clock
 * deadline + re-entry guard + per-top-level decline memo); every returned branch
 * is numerically self-verified (`ds_branch_num_ok`) and rejected if it carries an
 * inert Integrate, since the VoP particular's residual can be undecidable by
 * zero_test (Solve keep-on-undecidable would otherwise pass a non-closed form).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static time_t g_nh_deadline;
static bool nh_expired(void) { return time(NULL) >= g_nh_deadline; }

/* Re-entry guard, ALSO read by dsolve_specialform.c: while this method is solving
 * the homogeneous part (recursively), it needs SpecialFunctionForm's normal-form
 * pre-pass (the Bessel/Airy-with-y'-term recogniser) to fire even at recursion
 * depth > 1.  That pre-pass is otherwise gated to depth 1 because a caller like
 * OperatorFactor composes its mu = Exp[-Int P/2] factor BACK into a lower-order
 * reduction (infinite rewrite); this method does not — it extracts the basis and
 * applies variation of parameters — so the pre-pass is safe here. */
int dsolve_nh_vop_active = 0;
#define g_nh_active dsolve_nh_vop_active

/* decline memo (the evaluator re-invokes a declining builtin ~3x/call) */
#define NH_MEMO_SLOTS 32
static uint64_t nh_epoch = 0;
static int nh_memo_n = 0;
static uint64_t nh_memo[NH_MEMO_SLOTS];
static void nh_memo_sync(uint64_t tid){ if(tid!=nh_epoch){nh_epoch=tid;nh_memo_n=0;} }
static bool nh_memo_seen(uint64_t h){ for(int i=0;i<nh_memo_n;i++) if(nh_memo[i]==h) return true; return false; }
static void nh_memo_add(uint64_t h){ if(nh_memo_n<NH_MEMO_SLOTS && !nh_memo_seen(h)) nh_memo[nh_memo_n++]=h; }

/* PowerExpand[Simplify[e]] — canonicalise (1+Tan^2 -> Sec^2, ...) then reduce the
 * radicals so the VoP integrals close.  e consumed.  NOT wrapped in
 * TimeConstrained: this method is itself reached inside DSolve (and, in the corpus,
 * inside the prelude's TimeConstrained), and a nested TimeConstrained aborts the
 * whole subtree (the documented no-nest hazard); ds_simplify is internally
 * bounded, and the g_nh_deadline check between steps caps the loop. */
static Expr* nh_norm(Expr* e){
    Expr* s = ds_simplify(e);
    return eval_and_free(ds_call1(SYM_PowerExpand, s));
}

/* Solve DSolve[homog, y[x], x] -> body in x (applied form), or NULL.  `homeq`
 * consumed.  Called directly (no TimeConstrained wrapper — see nh_norm): the
 * homogeneous set for the target cases is claimed by SpecialFunctionForm /
 * ChangeOfVariable (each self-bounded) before Kovacic. */
static Expr* nh_solve_homog(Expr* homeq, const char* yname, const char* xv){
    Expr* lhs = ds_call1(yname, expr_new_symbol(xv));                 /* y[x] */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ homeq, lhs, expr_new_symbol(xv) }, 3);
    Expr* r = eval_and_free(call);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count && !body;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==yname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    return body;
}

/* Solve the homogeneous part with the PINNED `DSolve`SpecialFunctionForm` rather than
 * the whole cascade.  Two reasons, both load-bearing for the inert mode:
 *   - it is cheaper (measured 0.05 .. 0.81 s to solve, 0.014 s to self-decline);
 *   - it STRUCTURALLY guarantees a closed-form special-function fundamental set, so a
 *     truncated Frobenius/log series basis — over which an inert integral would be
 *     meaningless — can never reach the variation-of-parameters step.  That is the
 *     gate for the §2.2.34-3393/3394 class, held by construction rather than by a
 *     post-hoc SeriesData check (which is kept below as cheap redundancy).
 * `homeq` consumed. */
static Expr* nh_solve_homog_specialform(Expr* homeq, const char* yname, const char* xv){
    Expr* lhs = ds_call1(yname, expr_new_symbol(xv));
    Expr* call = expr_new_function(expr_new_symbol(intern_symbol("DSolve`SpecialFunctionForm")),
                     (Expr*[]){ homeq, lhs, expr_new_symbol(xv) }, 3);
    Expr* r = eval_and_free(call);
    Expr* body = NULL;
    if (head_is(r, SYM_List) && r->data.function.arg_count >= 1){
        Expr* inner = r->data.function.args[0];
        if (head_is(inner, SYM_List))
            for (size_t k=0;k<inner->data.function.arg_count && !body;k++){
                Expr* rule = inner->data.function.args[k];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count==2){
                    Expr* rl = rule->data.function.args[0];
                    if (rl->type==EXPR_FUNCTION && rl->data.function.head->type==EXPR_SYMBOL
                        && rl->data.function.head->data.symbol.name==yname)
                        body = expr_copy(rule->data.function.args[1]);
                }
            }
    }
    expr_free(r);
    return body;
}

/* Shared core of the two try-fns.  `inert` selects the M61 mode: see
 * dsolve_nonhomog_vop_inert_try below for what it changes and why. */
static Expr** nh_try_core(DSolveProblem* P, size_t* nbranch, bool inert){
    if (P->nfun!=1 || P->neq!=1) return NULL;
    if (P->max_order[0] < 2) return NULL;
    if (g_nh_active) return NULL;                    /* no self-recursion */
    /* An inert particular has no computable value at a point, so its constants cannot
     * be fitted.  Returning an unfitted general solution for an IVP would score as a
     * corpus PASS while being unsolved (the prelude reads the inert residual as UNK),
     * so decline -- the same discipline dsolve_run_parametric applies. */
    if (inert && P->ncond > 0) return NULL;
    const char* xv = P->ind_names[0];
    const char* yname = P->fun_names[0];

    Expr** c = NULL; Expr* g = NULL; int n = 0;
    if (!dsolve_linear_coeffs(P, &c, &g, &n)) return NULL;   /* nonlinear -> decline */
    if (n < 2) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }
    if (inert && n != 2) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }

    /* Homogeneous, or distributional forcing -> owned by other methods. */
    bool skip = ds_is_structural_zero(g)
             || ds_contains(g, intern_symbol("DiracDelta"))
             || ds_has_head(g, intern_symbol("UnitStep"))
             || ds_has_head(g, intern_symbol("Piecewise"));
    if (skip) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }

    if (!inert) {
        /* Fire ONLY for TRANSCENDENTAL (trigonometric) coefficients — the ChangeOfVariable
         * family, whose homogeneous set is found by a method that does not carry forcing.
         * A rational-coefficient nonhomogeneous equation is Kovacic's / Euler's domain
         * (each with its own forcing closure), so restricting to transcendental
         * coefficients keeps the one class this backstop uniquely reaches (e.g. §2.2.19-1822
         * Sin[x] y'' + ... == E^-x) while sparing every rational case a redundant
         * (and, on the large corpus, expensive) recursive homogeneous re-solve. */
        bool transc = false;
        for (int k=0;k<=n && !transc;k++)
            transc = ds_has_head(c[k],SYM_Sin)||ds_has_head(c[k],SYM_Cos)||ds_has_head(c[k],SYM_Tan)||
                     ds_has_head(c[k],SYM_Cot)||ds_has_head(c[k],SYM_Sec)||ds_has_head(c[k],SYM_Csc);
        if (!transc) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }
    } else {
        /* The inert mode's class: a RATIONAL-coefficient second-order equation with a
         * SINGULAR point (leading coefficient vanishing in x) and POLYNOMIAL forcing.
         * That is the regular-singular "series expansion" family whose homogeneous set
         * is Bessel / hypergeometric and whose VoP integrals therefore do not close.
         * Every other forcing shape (arbitrary f[x], transcendental, distributional) is
         * owned elsewhere and is excluded here, so this mode adds nothing to their cost. */
        bool ok = true;
        for (int k=0;k<=n && ok;k++) ok = ds_is_rational_in(c[k], xv);
        if (ok) {
            Expr* pq = eval_and_free(ds_call2(intern_symbol("PolynomialQ"),
                                              expr_copy(g), expr_new_symbol(xv)));
            ok = (pq->type == EXPR_SYMBOL && pq->data.symbol.name == SYM_True);
            expr_free(pq);
        }
        if (ok) ok = !ds_free_of(c[n], xv);          /* a singular point must exist */
        if (!ok) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }
    }

    uint64_t h = expr_hash(P->eq_residuals[0]) ^ (inert ? 0x9e3779b97f4a7c15ULL : 0ULL);
    nh_memo_sync(eval_toplevel_id());
    if (nh_memo_seen(h)) { for (int k=0;k<=n;k++) expr_free(c[k]); free(c); expr_free(g); return NULL; }

    g_nh_deadline = time(NULL) + 12;

    /* homogeneous equation: Sum_k c[k] y^(k)[x] == 0 */
    Expr** terms = malloc((size_t)(n+1)*sizeof(Expr*));
    for (int k=0;k<=n;k++)
        terms[k] = ds_call2(SYM_Times, expr_copy(c[k]), ds_make_funcapp(yname, k, xv));
    Expr* homlhs = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)(n+1)));
    free(terms);
    Expr* homeq = expr_new_function(expr_new_symbol(SYM_Equal),
                      (Expr*[]){ homlhs, expr_new_integer(0) }, 2);

    /* g_nh_active brackets the sub-solve in BOTH modes: it is also read by
     * dsolve_specialform.c to let its normal-form pre-pass (the Bessel/Airy
     * recogniser for an equation carrying a y' term) fire below recursion depth 1. */
    g_nh_active++;
    Expr* hb = inert ? nh_solve_homog_specialform(homeq, yname, xv)
                     : nh_solve_homog(homeq, yname, xv);
    g_nh_active--;

    Expr* body = NULL;
    if (hb && !ds_has_head(hb, intern_symbol("SeriesData")) && !nh_expired()){
        /* fundamental set: b_k = PowerExpand[Simplify[D[hom, C[k]]]] */
        Expr** basis = malloc((size_t)n*sizeof(Expr*));
        for (int k=0;k<n;k++) basis[k] = NULL;
        bool okbasis = true;
        for (int k=1;k<=n && okbasis && !nh_expired();k++){
            Expr* bk = eval_and_free(ds_call2(SYM_D, expr_copy(hb), ds_const(k)));
            bk = nh_norm(bk);
            if (!bk || !ds_free_of(bk, intern_symbol("C"))) { okbasis = false; if (bk) expr_free(bk); }
            else basis[k-1] = bk;
        }
        if (okbasis && !nh_expired()){
            bool was_inert = false;
            Expr* yp = dsolve_variation_of_parameters_mode(basis, (size_t)n, g, c[n], xv,
                           inert ? VP_ALLOW_INERT : VP_ELEMENTARY,
                           inert ? &was_inert : NULL);
            if (!inert) {
                if (yp && !ds_has_head(yp, SYM_Integrate)){
                    Expr* full = eval_and_free(ds_call2(SYM_Plus, expr_copy(hb), yp));
                    if (!ds_has_head(full, SYM_Integrate) && ds_branch_num_ok(P, full))
                        body = full;
                    else expr_free(full);
                } else if (yp) expr_free(yp);
            } else if (yp) {
                Expr* full = eval_and_free(ds_call2(SYM_Plus, expr_copy(hb), yp));
                /* A raw ACTIVE Integrate must never survive: it re-enters the
                 * integration cascade on every later re-evaluation (~1.7 s a pass). */
                if (ds_has_active_integrate(full)) full = ds_inactivate_integrate(full);
                bool has_inert = was_inert || ds_has_inactive_integrate(full);
                bool good = has_inert ? ds_inert_vop_verified(P, full)  /* the only barrier */
                                      : ds_branch_num_ok(P, full);      /* it closed after all */
                if (good) body = full; else expr_free(full);
            }
        }
        for (int k=0;k<n;k++) if (basis[k]) expr_free(basis[k]);
        free(basis);
    }
    if (hb) expr_free(hb);
    for (int k=0;k<=n;k++) expr_free(c[k]);
    free(c); expr_free(g);

    if (!body){ nh_memo_add(h); return NULL; }
    Expr** out = malloc(sizeof(Expr*));
    out[0] = body; *nbranch = 1;
    return out;
}

Expr** dsolve_nonhomog_vop_try(DSolveProblem* P, size_t* nbranch){
    return nh_try_core(P, nbranch, false);
}

/* M61 — the inert variation-of-parameters particular.
 *
 * For a second-order linear ODE with rational coefficients, a singular point and
 * polynomial forcing, the homogeneous fundamental set is a Bessel / hypergeometric
 * pair and the variation-of-parameters integrals do not close in elementary form.
 * Until now the case was simply lost: the homogeneous answer was found (by
 * `SpecialFunctionForm`) and thrown away because the equation had a right-hand side.
 * This returns the particular with its non-closing quadratures held INERT as
 * `Inactive[Integrate]` — Mathematica's own answer shape for this family.
 *
 * It occupies the LAST slot in the scalar cascade, after every other method
 * (including both Frobenius fallbacks) has declined.  That is a deliberate
 * correctness property, not a convenience: the slot can only convert UNEVAL into a
 * solution, so no case that already passes can be lost to it and its latency lands
 * only on equations that were returning nothing.
 *
 * Correctness rests entirely on `ds_inert_vop_verified`: the substrate's own
 * verifiers all KEEP a residual containing an integral, and `PossibleZeroQ` answers
 * True for one, so nothing downstream can catch a wrong inert answer.  See the gate's
 * comment in dsolve_common.h. */
Expr** dsolve_nonhomog_vop_inert_try(DSolveProblem* P, size_t* nbranch){
    return nh_try_core(P, nbranch, true);
}

static Expr* builtin_dsolve_nonhomog_vop(Expr* res){
    /* Pinned entry: try the elementary particular first, then the inert one — the
     * same explicit-then-implicit shape builtin_dsolve_autonomous uses (a declining
     * dsolve_method_builtin leaves `res` intact, so it is reused, not copied). */
    Expr* r = dsolve_method_builtin(res, dsolve_nonhomog_vop_try);
    if (!r) r = dsolve_method_builtin(res, dsolve_nonhomog_vop_inert_try);
    return r;
}

void dsolve_nonhomog_vop_init(void){
    symtab_add_builtin("DSolve`VariationOfParameters", builtin_dsolve_nonhomog_vop);
    symtab_get_def("DSolve`VariationOfParameters")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`VariationOfParameters",
        "DSolve`VariationOfParameters[eqn, y, x] solves a nonhomogeneous linear ODE "
        "L[y] == g(x) whose homogeneous part is solvable by a special-function / "
        "change-of-variable method that does not carry forcing: it solves L[y] == 0, "
        "forms the fundamental set, and builds the particular by variation of "
        "parameters (e.g. 4 x^2 y'' - 4 x y' + (3 - 16 x^2) y == 8 x^(5/2)).  When the "
        "variation-of-parameters quadratures are not elementary -- a second-order "
        "equation with rational coefficients, a singular point and polynomial forcing, "
        "whose homogeneous set is Bessel or hypergeometric -- the particular is returned "
        "with those quadratures held inert as Inactive[Integrate].");
}
