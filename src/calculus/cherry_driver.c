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

#include "expr.h"

Expr* extended_liouville_solve(Expr* f, Expr* x, unsigned top_mask) {
    /* Dispatch the registered special-function forms narrowed to the applicable
     * top monomials.  A NULL result is a routing decline (UNKNOWN), NOT a Cherry
     * NON-existence verdict — the latter is aggregated here only when a genuine
     * constant-existence / Sigma-decomposition certificate fires (future work). */
    return rt_special_case_routed(f, x, top_mask);
}
