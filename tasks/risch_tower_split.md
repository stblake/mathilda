# Split the Risch structure-theorem (tower construction) code out of integrate_risch_transcendental.c

## Context / interpretation

`integrate_risch_transcendental.c` is 5,832 lines. The RDE layer was already
split into `integrate_risch_rde.{c,h}`. The remaining "Risch structure theorem"
code in the big file is the **differential-tower construction** cluster:
`rt_tower_build_min` (the operational structure theorem — collect the log/exp/tan
kernels of the integrand, decide commensurability/independence, and assemble the
`RtTower`) plus the kernel-collection and tower-lifecycle helpers.

NOTE: `risch_structure.c` already exists but is a *different* concern — the
abstract structure-theorem decision builtins (`Risch`RationalSpan/LogReducible/
ExpReducible`). To avoid a naming collision and keep concerns separate, the new
module is named **`risch_tower.{c,h}`**.

The cluster's only upward call is `rt_field_integrate`, already an established
cross-file shared entry point (declared in `integrate_risch_rde.h`), so this
follows the exact pattern of the RDE split — no dependency cycle.

## Functions to MOVE to risch_tower.c
rt_powers_to_exp, rt_collect_exp_exponents, rt_collect_logs, rt_contains,
rt_build_monomial, rt_decode_mono, rt_log_of, rt_expand_logs, rt_tower_free,
rt_collect_tangents, rt_tower_build_min, rt_dt_i, rt_tower_deriv,
rt_build_deriv_rules, rt_has_explog_kernel, rt_expand_exp_sums, rt_subst_kernels

## Types
- Move `RtKind` + `RtTower` typedefs from integrate_risch_rde.h -> risch_tower.h
  (C99 forbids duplicate typedefs; exactly one header defines them).
- integrate_risch_rde.h `#include "risch_tower.h"` for the type + `RdeCtx`.

## Shared helpers that stay DEFINED in transcendental.c (un-static + extern in risch_tower.h)
rt_eval_own, rt_free_of_x, rt_is_rat_const, rt_class_primitive,
rt_find_exp_of_x, rt_find_log_of_x  (rt_field_integrate already shared).

## Steps
- [ ] risch_tower.h (types + moved-fn decls + extern shared-helper decls)
- [ ] risch_tower.c (moved fn bodies + includes)
- [ ] Remove moved bodies from transcendental.c; un-static the 6 shared helpers;
      `#include "risch_tower.h"`; drop now-duplicate externs
- [ ] integrate_risch_rde.h: include risch_tower.h; drop moved typedefs/externs
- [ ] tests/CMakeLists.txt COMMON_SRC += risch_tower.c
- [ ] Build clean (-std=c99 -Wall -Wextra) + tests green
- [ ] Write "further split" modularity suggestions

## Review

DONE. `risch_tower.{c,h}` created (825 + 103 lines); transcendental.c 5832 -> 5052.
17 functions moved verbatim (static stripped; RtKind/RtTower relocated to
risch_tower.h; rde.h includes it). Clean builds of main + all affected tests
(-std=c99 -Wall -Wextra, no warnings on the touched files).

Behavior-preserving, PROVEN:
- transcendental-Risch antiderivatives byte-identical pre/post (checked via main
  binary on the I-laden x Exp[x] Cos[x] case + tower cases).
- integrate_risch_transcendental_tests FAIL set IDENTICAL to HEAD (1 pre-existing
  known-gap soft-fail).
- risch_field / structure / structure_real / hermite / elementaryq / coupled /
  integrals / dispatch / intrischnorman: all green.
- integrate_jeffrey_tests aborts on BOTH HEAD and mine (pre-existing PZQ/Weierstrass
  flakiness, unrelated).

GOTCHA: the incremental tests/build dir carried STALE objects after the rde.h
change (macOS makefile/CMake header-dep hazard) and produced a phantom "8 vs 1"
soft-fail delta; a CLEAN rebuild collapsed it to the true "1 == HEAD". Always
clean-rebuild to confirm after touching a widely-included header.

## Further split (all 5 done)

integrate_risch_transcendental.c 5,832 -> 710-line driver + 6 modules:
  risch_util (362/57), risch_tower (825/97), risch_singleext (1186/31),
  risch_field_integrate (1748/49), risch_trig_frontend (744/26),
  risch_special (443/26).
Shared surface consolidated: rt_ helpers -> risch_util.h; RtKind/RtTower ->
risch_tower.h; RtDecision + g_rt_decide_mode/g_rt_decision -> risch_field_integrate.h.
Each step: verbatim script extraction + per-module header + CMake COMMON_SRC +
CLEAN `make` (c99 -Wall -Wextra, no warnings) + main-binary smoke.

Final verification: risch_field/structure/structure_real/hermite/elementaryq/
coupled/integrals/dispatch/intrischnorman all green; transcendental test FAIL-set
IDENTICAL to HEAD; 13 integrals across all 6 modules diff-back to 0.

PRE-EXISTING (not introduced): integrate_risch_transcendental_tests is built with
asserts ACTIVE and aborts at the first layout-sensitive I-laden Simplify soft-fail
(x Exp[x] Sin[x], inside test_multikernel_case). HEAD aborts identically. Verified
post-abort tower/field tests via the main binary instead.
