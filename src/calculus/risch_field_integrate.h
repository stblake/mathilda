/* risch_field_integrate.h — the recursive transcendental-tower integrator.
 *
 * The Bronstein field integrator that runs over the differential tower built by
 * risch_tower.c: given F in the tower field K_L, split it into polynomial /
 * special / normal parts and integrate each (Hermite reduction, the
 * Rothstein–Trager / Lazard–Rioboo–Trager log part, the coupled bounded-ansatz
 * for exponential/hypertangent monomials, and the field Risch DE), recursing
 * through the mutually-recursive `rt_field_integrate` entry (declared in
 * integrate_risch_rde.h, since the RDE's antidifferentiation branch calls it).
 *
 * Also hosts the flat multi-kernel "Phase B" tower cases (rt_log_tower_case /
 * rt_exp_tower_case, a SolveAlways ansatz) and the top-level recursive assembly
 * (rt_recursive_tower_case), plus the shared elementary-integrability decision
 * state that the driver in integrate_risch_transcendental.c reads to answer
 * ElementaryIntegralQ.  Defined in risch_field_integrate.c.
 */

#ifndef MATHILDA_RISCH_FIELD_INTEGRATE_H
#define MATHILDA_RISCH_FIELD_INTEGRATE_H

#include <stdbool.h>
#include "expr.h"

/* ---- Elementary-integrability decision state (shared with the driver) ------ */
/* The field integrator is a genuine decision procedure; in decision mode it
 * raises g_rt_decision to RT_DEC_NONELEMENTARY only at authoritative declines
 * (a non-constant residue, a Risch DE with no solution, the hypertangent Dc!=0
 * certificate).  The driver (rt_decide) sets g_rt_decide_mode, runs the same
 * path, and reads the verdict back. */
typedef enum { RT_DEC_UNKNOWN = 0, RT_DEC_ELEMENTARY, RT_DEC_NONELEMENTARY } RtDecision;

extern bool       g_rt_decide_mode;
extern RtDecision g_rt_decision;

/* ---- Integration entry points --------------------------------------------- */
/* (rt_field_integrate itself is declared in integrate_risch_rde.h, where the
 * RDE ↔ integrator mutual recursion is documented.) */

/* Genuine one-extension-at-a-time recursive Risch over a depth->=2 tower. */
Expr* rt_recursive_tower_case(Expr* f, Expr* x);

/* Flat multi-kernel bounded-ansatz fallbacks (Phase B). */
Expr* rt_log_tower_case(Expr* f, Expr* x);
Expr* rt_exp_tower_case(Expr* f, Expr* x);

/* Risch-DE builtin glue: solve D y + f y = g by driving the field integrator. */
Expr* rde_solve_tower(Expr* f, Expr* g, Expr* x);

#endif /* MATHILDA_RISCH_FIELD_INTEGRATE_H */
