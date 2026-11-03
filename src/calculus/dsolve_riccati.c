/*
 * dsolve_riccati.c — DSolve`Riccati.
 *
 * Solves the Riccati equation  y'[x] == q0(x) + q1(x) y + q2(x) y^2  (q2 != 0)
 * by the classical linearising substitution  y = -u'/(q2 u), which turns it into
 * the homogeneous second-order linear ODE
 *
 *     u'' - (q1 + q2'/q2) u' + (q0 q2) u == 0.
 *
 * That equation is handed back to the scalar DSolve cascade (const-coeff, Euler,
 * Airy/Bessel special forms, Kovacic, or — last — a Frobenius/power series), and
 * the recovered u is mapped back through y = -u'/(q2 u).  The general Riccati
 * solution depends only on the ratio C[1]:C[2] of the two linear constants, so
 * C[2] is collapsed to 1, leaving the single arbitrary constant C[1] (this omits
 * only the lone C[1]->Infinity particular sheet u = u1, the standard
 * one-parameter parametrisation, matching Mathematica).  When the linearised ODE
 * has no closed form the recursion declines and Riccati declines — no
 * implicit/wrong answer is shipped, and the substrate's back-substitution verify
 * is the authoritative backstop (permissive: an undecidable special-function
 * residual is kept, exactly as for Bessel/Kovacic).
 *
 * PARTICULAR-SOLUTION route (tried BEFORE the linearisation).  A
 * Riccati equation with one known solution y_p is solved by y = y_p + 1/v, which
 * leaves the LINEAR first-order equation  v' + (q1 + 2 q2 y_p) v == -q2  (closed
 * by the shared integrating-factor solver).  y_p is searched for among a small,
 * fixed candidate family -- c x^k (k = -2..2), c f[x] for every arbitrary function
 * application f[x] in the equation, and c K[u] for the six trig kernels K of every
 * trig argument u present -- the constant c fixed by Solve at one exact sample
 * point and the candidate kept only when the full symbolic residual is zero.
 * That family is what the corpus' "particular solution by inspection" Riccati
 * equations carry: y = x (x y' == y + (x^2-y^2) f[x]), y = 1/x, y = -1/x
 * (y' == f + x f y + y^2), y = f[x] (y' + f^2 == f' + y^2), y = Cot[x]
 * (2y' + 2Csc^2 == y Csc Sec - y^2 Sec^2) -- §2.2.47/48/49/50.  The final
 * answer is back-substitution verified by the substrate like any other branch.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../common.h"
#include "../core.h"          /* leaf_count_internal */
#include <stdlib.h>

/* base^-1; base consumed. */
static Expr* powneg1(Expr* base) {
    return expr_new_function(expr_new_symbol(SYM_Power),
                             (Expr*[]){ base, expr_new_integer(-1) }, 2);
}

/* True if e mentions a named special function.  Full Simplify can HANG on these
 * (e.g. LegendreP/Q of half-integer order from a Legendre-type u-ODE), and it
 * would not tidy them anyway — the evaluator already folds the -u'/(q2 u) sign
 * and reciprocal.  So the back-mapped body skips Simplify when this holds. */
static bool riccati_has_special(const Expr* e) {
    static const char* names[] = {
        "LegendreP", "LegendreQ", "AiryAi", "AiryBi",
        "BesselJ", "BesselY", "BesselI", "BesselK",
        "Hypergeometric1F1", "Hypergeometric2F1", "HypergeometricPFQ",
        "HypergeometricU", "WhittakerM", "WhittakerW", "MathieuC", "MathieuS", NULL };
    for (int i = 0; names[i]; i++)
        if (ds_has_head(e, intern_symbol(names[i]))) return true;
    return false;
}

/* Extract the RHS of {{u[x] -> expr}} (applied form). */
static Expr* extract_applied(Expr* r, const char* ufun) {
    if (!head_is(r, SYM_List) || r->data.function.arg_count < 1) return NULL;
    Expr* inner = r->data.function.args[0];
    if (!head_is(inner, SYM_List)) return NULL;
    for (size_t k = 0; k < inner->data.function.arg_count; k++) {
        Expr* rule = inner->data.function.args[k];
        if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
            Expr* lhs = rule->data.function.args[0];
            if (lhs->type == EXPR_FUNCTION && lhs->data.function.head->type == EXPR_SYMBOL
                && lhs->data.function.head->data.symbol.name == ufun)
                return expr_copy(rule->data.function.args[1]);
        }
    }
    return NULL;
}

/* Zero test for the Riccati bookkeeping: STRUCTURAL first (Expand of the Together
 * form is 0 -- exact, and the only route when the expression carries an arbitrary
 * function f[x], on which zero_test cannot sample and answers "undecided"), then
 * zero_test for the transcendental identities (Cot' == -Csc^2, ...). */
static bool ric_is_zero(const Expr* e) {
    if (ds_is_structural_zero(e)) return true;
    Expr* t = eval_and_free(ds_call1(SYM_Together, expr_copy((Expr*)e)));
    bool z = ds_is_structural_zero(t);
    expr_free(t);
    return z || ds_is_zero(e);
}

/* Collect the distinct arbitrary-function applications h[x] (a head with no
 * builtin and no DownValues, single argument == x, not the unknown y) of e. */
static void ric_collect_arbfun(const Expr* e, const char* x, const char* yname,
                               Expr** out, int* n, int cap) {
    if (!e || e->type != EXPR_FUNCTION || *n >= cap) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1
        && h->data.symbol.name != yname
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == x) {
        SymbolDef* d = symtab_lookup(h->data.symbol.name);
        if (!d || (!d->builtin_func && !d->down_values)) {
            for (int i = 0; i < *n; i++) if (expr_eq(out[i], (Expr*)e)) return;
            out[(*n)++] = expr_copy((Expr*)e);
            return;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ric_collect_arbfun(e->data.function.args[i], x, yname, out, n, cap);
}

/* Collect the distinct arguments u of the trig kernels Sin..Csc[u] in e. */
static void ric_collect_trigargs(const Expr* e, Expr** out, int* n, int cap) {
    if (!e || e->type != EXPR_FUNCTION || *n >= cap) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Sin || hn == SYM_Cos || hn == SYM_Tan || hn == SYM_Cot
            || hn == SYM_Sec || hn == SYM_Csc) {
            const Expr* u = e->data.function.args[0];
            for (int i = 0; i < *n; i++) if (expr_eq(out[i], (Expr*)u)) return;
            out[(*n)++] = expr_copy((Expr*)u);
            return;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ric_collect_trigargs(e->data.function.args[i], out, n, cap);
}

/* Riccati residual of a candidate y:  y' - (q0 + q1 y + q2 y^2).  y borrowed. */
static Expr* ric_residual(const Expr* y, const Expr* q0, const Expr* q1, const Expr* q2,
                          const char* x) {
    Expr* yp = ds_d(expr_copy((Expr*)y), expr_new_symbol(x));
    Expr* rhs = ds_call2(SYM_Plus, expr_copy((Expr*)q0),
                  ds_call2(SYM_Plus,
                      ds_call2(SYM_Times, expr_copy((Expr*)q1), expr_copy((Expr*)y)),
                      ds_call2(SYM_Times, expr_copy((Expr*)q2),
                          ds_call2(SYM_Power, expr_copy((Expr*)y), expr_new_integer(2)))));
    return eval_and_free(ds_call2(SYM_Subtract, yp, rhs));
}

/* Try y = c*s for a constant c.  The residual is a quadratic r0 + r1 c + r2 c^2 in
 * c with coefficients in x, so c is read off symbolically -- -r0/r1 when r2 == 0,
 * else (-r1 +- PowerExpand[Sqrt[Factor[r1^2 - 4 r0 r2]]])/(2 r2) -- and only an
 * x-FREE root whose full residual is zero is kept.  (Solving at a sample point
 * instead leaves Sqrt[(f'[p] - 2 f[p]^2)^2] unreduced when the equation carries an
 * arbitrary function f, and never exhibits the root c == 1.)  PowerExpand picks a
 * branch of the square root; both signs are tried and each is checked by the
 * residual, so the choice is only a search heuristic.  s borrowed; owned y_p or NULL. */
static Expr* ric_try_scaled(const Expr* s, const Expr* q0, const Expr* q1, const Expr* q2,
                            const char* x) {
    const char* cn = intern_symbol("DSolve`ricC");
    Expr* y = eval_and_free(ds_call2(SYM_Times, expr_new_symbol(cn), expr_copy((Expr*)s)));
    Expr* R = ric_residual(y, q0, q1, q2, x);
    expr_free(y);
    R = eval_and_free(ds_call1(SYM_Expand, R));
    Expr* r[3];
    for (int k = 0; k < 3; k++)
        r[k] = eval_and_free(ds_call1(SYM_Together, expr_new_function(expr_new_symbol(SYM_Coefficient),
                   (Expr*[]){ expr_copy(R), expr_new_symbol(cn), expr_new_integer(k) }, 3)));
    Expr* roots[2] = { NULL, NULL };
    if (ds_is_structural_zero(r[2])) {
        if (!ds_is_structural_zero(r[1]))
            roots[0] = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1),
                           ds_call2(SYM_Times, expr_copy(r[0]), powneg1(expr_copy(r[1])))));
    } else {
        Expr* disc = eval_and_free(ds_call2(SYM_Subtract,
                         ds_call2(SYM_Power, expr_copy(r[1]), expr_new_integer(2)),
                         ds_call2(SYM_Times, expr_new_integer(4),
                             ds_call2(SYM_Times, expr_copy(r[0]), expr_copy(r[2])))));
        Expr* sq = eval_and_free(ds_call1("PowerExpand",
                       ds_call1(SYM_Sqrt, eval_and_free(ds_call1(SYM_Factor, disc)))));
        for (int sg = 0; sg < 2; sg++)
            roots[sg] = eval_and_free(ds_call2(SYM_Times,
                ds_call2(SYM_Plus, ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(r[1])),
                         ds_call2(SYM_Times, expr_new_integer(sg ? -1 : 1), expr_copy(sq))),
                powneg1(ds_call2(SYM_Times, expr_new_integer(2), expr_copy(r[2])))));
        expr_free(sq);
    }
    Expr* found = NULL;
    for (int i = 0; i < 2; i++) {
        if (!roots[i]) continue;
        Expr* c = eval_and_free(ds_call1(SYM_Together, roots[i]));
        if (!found && ds_free_of(c, x) && !ds_is_zero(c)) {
            Expr* Rc = ds_subst(expr_copy(R), expr_new_symbol(cn), expr_copy(c));
            if (ric_is_zero(Rc))
                found = eval_and_free(ds_call2(SYM_Times, expr_copy(c), expr_copy((Expr*)s)));
            expr_free(Rc);
        }
        expr_free(c);
    }
    for (int k = 0; k < 3; k++) expr_free(r[k]);
    expr_free(R);
    return found;
}

/* Search the candidate family (header note) for a particular solution. */
static Expr* ric_particular(const Expr* q0, const Expr* q1, const Expr* q2,
                            const char* x, const char* yname) {
    const Expr* qs[3] = { q0, q1, q2 };      /* the candidates' scan targets */
    /* q0 == 0: y == 0 is the particular solution (the Bernoulli case), giving the
     * standard 1/v family -- prefer it to any nonzero candidate, whose equivalent
     * family parametrises the constant awkwardly (C/(C r^2 - r) for 1/(r(C + r)),
     * which a nested Lie-symmetry sub-solve then cannot carry: §2.2.45-4431). */
    if (ric_is_zero(q0)) return expr_new_integer(0);
    Expr* xs = expr_new_symbol(x);
    Expr* cands[5] = {
        expr_copy(xs),
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(xs), expr_new_integer(-1) }, 2),
        expr_new_integer(1),
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(xs), expr_new_integer(2) }, 2),
        expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(xs), expr_new_integer(-2) }, 2) };
    expr_free(xs);
    Expr* yp = NULL;
    for (int i = 0; i < 5; i++) {
        if (!yp) { Expr* c = eval_and_free(cands[i]); yp = ric_try_scaled(c, q0, q1, q2, x); expr_free(c); }
        else expr_free(cands[i]);
    }
    if (!yp) {
        Expr* af[8]; int naf = 0;
        for (int j = 0; j < 3; j++) ric_collect_arbfun(qs[j], x, yname, af, &naf, 8);
        for (int i = 0; i < naf; i++) {
            if (!yp) yp = ric_try_scaled(af[i], q0, q1, q2, x);
            expr_free(af[i]);
        }
    }
    if (!yp) {
        Expr* ta[4]; int nta = 0;
        for (int j = 0; j < 3; j++) ric_collect_trigargs(qs[j], ta, &nta, 4);
        const char* ks[6] = { SYM_Tan, SYM_Cot, SYM_Sin, SYM_Cos, SYM_Sec, SYM_Csc };
        for (int i = 0; i < nta; i++) {
            for (int k = 0; k < 6 && !yp; k++) {
                Expr* s = eval_and_free(ds_call1(ks[k], expr_copy(ta[i])));
                yp = ric_try_scaled(s, q0, q1, q2, x);
                expr_free(s);
            }
            expr_free(ta[i]);
        }
    }
    return yp;
}

/* y = y_p + 1/v with v' + (q1 + 2 q2 y_p) v == -q2.  All borrowed; owned body or NULL. */
static Expr* ric_from_particular(const Expr* yp, const Expr* q1, const Expr* q2, const char* x,
                                  bool general) {
    Expr* Pc = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)q1),
                   ds_call2(SYM_Times, expr_new_integer(2),
                       ds_call2(SYM_Times, expr_copy((Expr*)q2), expr_copy((Expr*)yp)))));
    Expr* Qc = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)q2)));
    Expr* v = general ? dsolve_linear_factor_solve_inert(Pc, Qc, x)    /* consumes Pc, Qc */
                      : dsolve_linear_factor_solve(Pc, Qc, x);
    if (!v) return NULL;
    Expr* body = eval_and_free(ds_call2(SYM_Plus, expr_copy((Expr*)yp), powneg1(v)));
    return ds_has_head(body, SYM_Integrate) ? body : ds_simplify(body);   /* never Simplify an inert body */
}

Expr** dsolve_riccati_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1) return NULL;
    if (P->max_order[0] != 1) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);        /* y' == F(x, y) */
    if (!F) return NULL;

    const char* Yn = intern_symbol("DSolve`Y");
    Expr* FY = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Yn));

    /* quadratic coefficients:  q0 = FY|_{Y=0},  q1 = FY_Y|_{Y=0},  q2 = FY_YY/2 */
    Expr* dFY = ds_d(expr_copy(FY), expr_new_symbol(Yn));
    Expr* q0  = ds_subst(expr_copy(FY),  expr_new_symbol(Yn), expr_new_integer(0));
    Expr* q1  = ds_subst(expr_copy(dFY), expr_new_symbol(Yn), expr_new_integer(0));
    Expr* q2  = eval_and_free(ds_call2(SYM_Times,
                    ds_d(dFY, expr_new_symbol(Yn)),         /* dFY consumed here */
                    powneg1(expr_new_integer(2))));         /* * 1/2 */

    /* genuine Riccati: coefficients free of Y, q2 != 0, and FY reconstructs
     * exactly as q0 + q1 Y + q2 Y^2 (rejects a cubic-or-higher / rational-in-Y
     * nonlinearity that no first-order quadratic form can carry). */
    bool ok = ds_free_of(q0, Yn) && ds_free_of(q1, Yn) && ds_free_of(q2, Yn)
              && !ds_is_zero(q2);
    if (ok) {
        Expr* Y2 = eval_and_free(ds_call2(SYM_Power, expr_new_symbol(Yn), expr_new_integer(2)));
        Expr* recon = eval_and_free(ds_call2(SYM_Plus,
            expr_copy(q0),
            ds_call2(SYM_Plus,
                ds_call2(SYM_Times, expr_copy(q1), expr_new_symbol(Yn)),
                ds_call2(SYM_Times, expr_copy(q2), Y2))));
        Expr* chk = eval_and_free(ds_call2(SYM_Subtract, expr_copy(FY), recon));
        ok = ric_is_zero(chk);
        expr_free(chk);
    }
    expr_free(FY);
    if (!ok) { expr_free(q0); expr_free(q1); expr_free(q2); return NULL; }

    /* A particular solution by inspection first (header note): when one exists it
     * gives the elementary one-parameter family directly, and spares the recursive
     * second-order solve of the linearised equation -- which can succeed with a
     * special-function form that then fails verify, or spin (2y' + 2Csc[x]^2 ==
     * y Csc Sec - y^2 Sec^2 ran 40 s there; y_p == Cot[x] closes it at once). */
    /* Not when the equation already carries a generated constant C[k] -- a nested
     * sub-solve such as ReducibleFirstIntegral's y' == C[1] y + y^2 (§2.2.45-4431):
     * the particular found there (y_p == -C[1]) and the C[1] the linear solver mints
     * for v are the SAME symbol, collapsing the family to one constant. */
    bool has_C = ds_contains(q0, intern_symbol("C")) || ds_contains(q1, intern_symbol("C"))
                 || ds_contains(q2, intern_symbol("C"));
    if (!has_C) {
        Expr* yp = ric_particular(q0, q1, q2, xvar, yname);
        Expr* pbody = yp ? ric_from_particular(yp, q1, q2, xvar, P->ncond == 0) : NULL;
        if (yp) expr_free(yp);
        if (pbody) {
            expr_free(q0); expr_free(q1); expr_free(q2);
            Expr** out = malloc(sizeof(Expr*));
            out[0] = pbody;
            *nbranch = 1;
            return out;
        }
    }

    /* build the u-ODE:  u'' - (q1 + q2'/q2) u' + (q0 q2) u == 0 */
    const char* ufun = intern_symbol("DSolve`ricU");
    Expr* q2p = ds_d(expr_copy(q2), expr_new_symbol(xvar));                       /* q2' */
    Expr* b1  = eval_and_free(ds_call2(SYM_Plus, expr_copy(q1),
                    ds_call2(SYM_Times, q2p, powneg1(expr_copy(q2)))));           /* q1 + q2'/q2 */
    Expr* b0  = eval_and_free(ds_call2(SYM_Times, expr_copy(q0), expr_copy(q2))); /* q0 q2 */

    Expr* t2 = ds_make_funcapp(ufun, 2, xvar);
    Expr* t1 = ds_call2(SYM_Times,
                   ds_call2(SYM_Times, expr_new_integer(-1), b1),
                   ds_make_funcapp(ufun, 1, xvar));
    Expr* t0 = ds_call2(SYM_Times, b0, ds_make_funcapp(ufun, 0, xvar));
    Expr* lhs = ds_call2(SYM_Plus, t2, ds_call2(SYM_Plus, t1, t0));
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                    (Expr*[]){ lhs, expr_new_integer(0) }, 2);

    /* recurse into the scalar cascade for u (applied form -> {{u[x] -> body}}) */
    Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                     (Expr*[]){ eqn, ds_make_funcapp(ufun, 0, xvar), expr_new_symbol(xvar) }, 3);
    Expr* r = eval_and_free(call);
    Expr* ubody = extract_applied(r, ufun);
    expr_free(r);
    if (!ubody) { expr_free(q0); expr_free(q1); expr_free(q2); return NULL; }

    /* collapse the redundant constant (Riccati has one parameter): C[2] -> 1 */
    ubody = ds_subst(ubody, ds_const(2), expr_new_integer(1));

    /* y = -u'/(q2 u) */
    Expr* up   = ds_d(expr_copy(ubody), expr_new_symbol(xvar));
    Expr* num  = ds_call2(SYM_Times, expr_new_integer(-1), up);
    Expr* den  = eval_and_free(ds_call2(SYM_Times, expr_copy(q2), ubody));   /* consumes ubody */
    Expr* raw  = eval_and_free(ds_call2(SYM_Times, num, powneg1(den)));
    /* Never Simplify across an unevaluated Integrate (as ric_from_particular):
     * Simplify's transforms rewrite inside it and every rewrite re-runs the
     * integration cascade -- y' == 2x - (1 + x^2) y + y^2 (section 2.2.47-4658,
     * u carries Integrate[E^(x + x^3/3)]) spent its whole budget there. */
    Expr* body;
    if (riccati_has_special(raw) || ds_has_head(raw, SYM_Integrate)
        || ds_has_head(raw, SYM_Inactive)) body = raw;
    else if (leaf_count_internal(raw, true) > 600) {
        /* A large logarithmic derivative -- u'/u of a product of powers with radical
         * symbolic exponents, (c x^2 + b x + a)^2 (y' + y^2) + A == 0 (§2.2.50-5000),
         * 1225 leaves -- sends Simplify past any budget (> 20 s); the tidy form is
         * worth at most a bounded attempt, the raw quotient is already the answer. */
        Expr* keep = expr_copy(raw);
        body = eval_and_free(expr_new_function(expr_new_symbol("TimeConstrained"),
                   (Expr*[]){ ds_call1(SYM_Simplify, raw), expr_new_integer(2), keep }, 3));
    } else body = ds_simplify(raw);

    expr_free(q0); expr_free(q1); expr_free(q2);

    Expr** out = malloc(sizeof(Expr*));
    out[0] = body;
    *nbranch = 1;
    return out;
}

static Expr* builtin_dsolve_riccati(Expr* res) {
    return dsolve_method_builtin(res, dsolve_riccati_try);
}

void dsolve_riccati_init(void) {
    symtab_add_builtin("DSolve`Riccati", builtin_dsolve_riccati);
    symtab_get_def("DSolve`Riccati")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Riccati",
        "DSolve`Riccati[eqn, y, x] solves the Riccati equation "
        "y'[x] == q0(x) + q1(x) y + q2(x) y^2 (q2 != 0). It first looks for a particular "
        "solution y_p among c x^k (k = -2..2), c f[x] for arbitrary functions f, and c K[u] "
        "for the trig kernels K of the equation, giving y = y_p + 1/v with v linear; "
        "otherwise the substitution y = -u'/(q2 u) linearises it to a second-order linear "
        "ODE that the scalar cascade solves, and u is mapped back to y.");
}
