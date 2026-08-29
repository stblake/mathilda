/* cherry_driver.c — extended-Liouville dispatch for Cherry's special functions.
 *
 * The C0 seam (CHERRY_PLAN.md §2.3).  A single dispatch point over the
 * special-function form registry (risch_special.c), shared by the outermost
 * Integrate`RischTranscendental dispatch and the Thm 5.4 tower hook inside the
 * field recursion (risch_field_integrate.c).  Because both callers funnel through
 * here, the Cherry engines (ExpIntegralEi / Erf / LogIntegral / PolyLog) fire
 * uniformly on the original integrand AND on every peeled tower-monomial
 * coefficient the field recursion hands back down.
 *
 * The body currently delegates straight to rt_special_case_routed — behaviour is
 * byte-identical to the historical rt_special_case loop when top_mask ==
 * RT_SF_TOP_ANY.  The decision-aggregation seam (a Cherry constant-existence /
 * Sigma-decomposition NON-existence certificate calling rt_dec_nonelem so
 * ElementaryIntegralQ can answer False) lands here as the engines grow their
 * finite generators; today every decline is a soft NULL -> RT_DEC_UNKNOWN, per the
 * A4 hazard-1 rule (never a spurious NONELEMENTARY at a routing decline).
 */

#include "cherry_driver.h"
#include "risch_special.h"
#include "cherry_ei.h"
#include "cherry_li.h"
#include "cherry_dilog.h"
#include "cherry_dilog_exp.h"

#include "expr.h"
#include "symtab.h"
#include "attr.h"

Expr* extended_liouville_solve(Expr* f, Expr* x, unsigned top_mask) {
    /* Dispatch the registered special-function forms narrowed to the applicable
     * top monomials.  A NULL result is a routing decline (UNKNOWN), NOT a Cherry
     * NON-existence verdict — the latter is aggregated here only when a genuine
     * constant-existence / Sigma-decomposition certificate fires (future work). */
    return rt_special_case_routed(f, x, top_mask);
}

/* ==================================================================== */
/* Direct debuggable REPL surfaces for the Cherry engines.               */
/*                                                                       */
/* Each Integrate`Cherry`<Name>[f, x] applies one Cherry engine to [f, x]*/
/* directly, bypassing the integration cascade — a uniform inspection /  */
/* testing / benchmarking hook (the engines are also reached             */
/* automatically from Integrate).  Each returns the engine's fresh,      */
/* diff-back-verified antiderivative or leaves the call unevaluated.      */
/* ==================================================================== */

#define CHERRY_BUILTIN(fn, engine)                                          \
    static Expr* fn(Expr* res) {                                            \
        if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)\
            return NULL;                                                    \
        return engine(res->data.function.args[0],                          \
                      res->data.function.args[1]);                         \
    }
CHERRY_BUILTIN(builtin_cherry_ei,           rt_cherry_ei)
CHERRY_BUILTIN(builtin_cherry_exp_multiterm, rt_cherry_exp_multiterm)
CHERRY_BUILTIN(builtin_cherry_li,           rt_cherry_li)
CHERRY_BUILTIN(builtin_cherry_dilog,        rt_cherry_dilog)
#undef CHERRY_BUILTIN
/* builtin_cherry_dilog_exp is defined in cherry_dilog_exp.c (its engine also
 * carries a top-level depth gate). */

static void reg_cherry(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void cherry_builtins_init(void) {
    reg_cherry("Integrate`Cherry`Ei", builtin_cherry_ei,
        "Integrate`Cherry`Ei[f, x] applies the Cherry base-field exponential-"
        "integral engine (Cherry 1989) to f, bypassing the cascade: for f = g E^h "
        "with g, h rational in x it returns the y E^h + Sum c_i ExpIntegralEi[h+a_i] "
        "antiderivative (also the Erfi Gaussian case), else leaves the call "
        "unevaluated.");
    reg_cherry("Integrate`Cherry`ExpMultiterm", builtin_cherry_exp_multiterm,
        "Integrate`Cherry`ExpMultiterm[f, x] applies the multi-term exponential-"
        "integral engine (Cherry Thm 5.4 case b): a sum Sum_i p_i E^(i w) rational "
        "in a single kernel E^w is integrated term-by-term into ExpIntegralEi, else "
        "the call is left unevaluated.");
    reg_cherry("Integrate`Cherry`Li", builtin_cherry_li,
        "Integrate`Cherry`Li[f, x] applies the Cherry logarithmic-integral engine "
        "(Cherry 1986) to f, bypassing the cascade: a single-log tower theta = "
        "Log[w] gives v(x, theta) + Sum_k d_k LogIntegral[w^k] (with the "
        "transcendental-constant rescaling and Laurent -> ExpIntegralEi cases), "
        "else the call is left unevaluated.");
    reg_cherry("Integrate`Cherry`Dilog", builtin_cherry_dilog,
        "Integrate`Cherry`Dilog[f, x] applies the logarithmic-tower dilogarithm "
        "engine (Cherry degree-2 Sigma-decomposition) to f: R(x) Log[w] -> Log-Log "
        "products + PolyLog[2, g] with g the rational-root interpolants, else the "
        "call is left unevaluated.");
    reg_cherry("Integrate`Cherry`DilogExp", builtin_cherry_dilog_exp,
        "Integrate`Cherry`DilogExp[f, x] applies the exponential-tower dilogarithm "
        "engine (the exp-tower mirror of Integrate`Cherry`Dilog) to f: a rational-"
        "in-E^(c x) x-weighted form (x/(E^x-1)) or an outer-log form (Log[1+E^x]) "
        "-> PolyLog[2, ...], else the call is left unevaluated.");
}
