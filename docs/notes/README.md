# Development notes

This directory holds Mathilda's internal working documents: implementation
plans, audits, investigation logs and post-mortems written while subsystems
were being built. They record *why* the code is shaped the way it is, and source
comments cite them (for example `docs/notes/REDUCE_PLAN.md, Phase 3`).

They are historical records, not user documentation. A plan may describe work
that has since shipped, changed direction or been abandoned. For current
behaviour, see the [function reference](../../Mathilda_spec.md), the
[architecture overview](../../SPEC.md) and the
[design notes](../design/).

| Area | Documents |
|------|-----------|
| Integration | [`TRAGER_PLAN.md`](TRAGER_PLAN.md), [`ALGEBRAIC_INTEGRAL_PERF.md`](ALGEBRAIC_INTEGRAL_PERF.md), [`MATHILDA_PARALLEL_MIXED_SPECIAL_PLAN.md`](MATHILDA_PARALLEL_MIXED_SPECIAL_PLAN.md), `plans/INTEGRATE_PLAN.md`, `plans/RISCH_*`, `plans/RDE_SPDE_PLAN.md`, `plans/CONTOUR_INTEGRATE_PLAN.md`, `plans/DEFINITE.md` |
| Differential equations | [`DSOLVE_PLAN.md`](DSOLVE_PLAN.md), `plans/GOURSAT_*` |
| Solving and reduction | [`REDUCE_PLAN.md`](REDUCE_PLAN.md), [`RESOLVE_IMPROVEMENTS.md`](RESOLVE_IMPROVEMENTS.md), [`SOLVE_INTEGERS.md`](SOLVE_INTEGERS.md), [`SOLVERAD_IMPROVEMENTS.md`](SOLVERAD_IMPROVEMENTS.md), `plans/SOLVE_PLAN.md` |
| Algebra and simplification | [`SIMPLIFY_IMPROVEMENTS.md`](SIMPLIFY_IMPROVEMENTS.md), `plans/FACTOR_PLAN.md`, `plans/CANCEL_IMPROVEMENT_PLAN.md`, `plans/SIMPLIFY_GAPS.md`, `plans/TRIG_REDUCE_PLAN.md`, `plans/FLINT_*` |
| Series, limits, sums, products | `plans/limit_candidate_spec.md`, `plans/LIMIT_OSCILLATORY_TODO.md`, `plans/SUM_DEVEL_PLAN.md`, `plans/PRODUCT_*`, `plans/NPRODUCT_PLAN.md` |
| Evaluator and symbol table | [`SYMBOL_TABLE_IMPROVEMENTS.md`](SYMBOL_TABLE_IMPROVEMENTS.md), `plans/EVAL_IMPROVEMENTS_PLAN.md`, `plans/MATHILDA_REVIEW.md` |
| Numerics and performance | [`NMINIMIZE_FEASIBILITY_BUG.md`](NMINIMIZE_FEASIBILITY_BUG.md), `plans/HPC_IMPROVEMENT_PLAN.md`, `plans/PACKED_ARRAYS_TODO.md`, `plans/COMPILE_M5_TODO.md`, `plans/FEAST_PLAN.md`, `plans/plan_numerics.md` |
| Correctness logs | [`MATHILDA_DIVERGENCES.md`](MATHILDA_DIVERGENCES.md), [`SEGFAULTS_TO_FIX.md`](SEGFAULTS_TO_FIX.md), [`LEAK.md`](LEAK.md) |
| Documentation and front end | [`BUILTIN_DOCUMENTATION_OVERHAUL.md`](BUILTIN_DOCUMENTATION_OVERHAUL.md), [`frontend-research.md`](frontend-research.md), [`canvas-plan.md`](canvas-plan.md), [`canvas-library-plan.md`](canvas-library-plan.md) |
| Reviews | `reports/association/` (Association audit), `thoughts/` (agent research and plan records), [`KIT-FEEDBACK.md`](KIT-FEEDBACK.md), [`KIT-FEEDBACK-V2.md`](KIT-FEEDBACK-V2.md) |
