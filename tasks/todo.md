# Builtin Documentation Overhaul — task tracker

Full roadmap: `BUILTIN_DOCUMENTATION_OVERHAUL.md` (repo root).
Plan: `~/.claude/plans/i-would-like-all-velvety-bee.md`.

Goal: bring every builtin's docs page to EllipticF-grade (overlay worked
examples + Notes, source-grounded impl note, rich verified examples). Every
example runs through the latest binary.

## Campaign 0 — infrastructure (nearly done)
- [x] Build latest binary (`make`) — v0.266
- [x] Generator change: overlay worked-examples auto-verified (binary supplies
      Out[]) in `site/generate.py` `render_page`. Validated: regen fixed 83 pages
      (fictional N[]-precision outputs, now-included dropped setup lines)
- [x] `site/coverage_report.py` — per-category parity counts (baseline 369/1086)
- [x] `site/verify_docs_examples.py` + `make check-docs-examples`; green on 1086,
      FAILS on planted wrong expectation, ~10s. 7-page EXEMPT (nondeterministic +
      non-transcript-order). CI wiring deferred (needs full-dep job)
- [x] `BUILTIN_DOCUMENTATION_OVERHAUL.md` — living tracker
- [x] Pilot: elementary-functions — impl+overlay for Log10, Log2, UnitStep,
      UnitBox, Ramp; category now 100% graded (26/26). Total 374/1086 (34%)

## Campaigns 1..N — content — COMPLETE (1086/1086, 100%)
All 38 categories at 100% EllipticF-grade, done across 5 subagent waves:
- Wave 1 (math core): elementary-functions, arithmetic, calculus, number-theory,
  linear-algebra, simplification, power-series, solutions-of-equations.
- Wave 2: special-functions, algebra, numerical-calculus, mathematical-constants.
- Wave 3: assignment-and-rules, pattern-matching, scoping-constructs,
  functional-programming, control-flow, structural-manipulation,
  string-operations, lists-and-iteration, time-and-date.
- Wave 4: data-structures, expression-information, file-io,
  random-number-generation, statistics, hypergraphs, graphics, machine-learning,
  fourier-transforms, geometry, packed-arrays, bitwise.
- Wave 5 (bare blocks): graphs, other-advanced, image-processing.

Verification: `make docs` verifies 9150 examples; `make check-docs-examples`
green on 1086 pages (7-page EXEMPT: nondeterministic + rendered-order cases);
strict `mkdocs build` clean after fixing 68 broken same-dir links in file-io
fragments.

### Follow-ups discovered (out of scope — NOT fixed here)
- A few statistics heads (`InterquartileRange`, `Quantile`, `MeanDeviation`,
  `MedianDeviation`) emit `::rectn`/`::q100` via raw `printf`, bypassing the
  `mth_message` funnel — a potential `make check-messages` concern.
- Partial builtins surfaced and documented honestly: `SparseArray` (inert; only
  `Normal` materializes), `UniformDistribution` (PDF/RandomVariate wired;
  Mean/Variance/CDF not), `RatCanonPrototype` (Phase-1), `BesselJZero`
  (symbolic-only). Spec docs now state a false `MaxMemoryUsed >= MemoryInUse`
  claim was removed (they read different OS counters).
- [ ] C2 arithmetic
- [ ] C3 calculus
- [ ] C4 number-theory
- [ ] C5 special-functions (impl gap: 5/47)
- [ ] C6 linear-algebra
- [ ] C7 power-series / simplification / solutions-of-equations / comparisons /
      mathematical-constants
- [ ] non-core partials (expression-information, data-structures, functional-
      programming, structural-manipulation, control-flow, assignment-and-rules,
      scoping, pattern-matching, string-operations, statistics,
      numerical-calculus, flint, file-io, lists-and-iteration, time-and-date, rng)
- [ ] bare blocks (graphs, other-advanced, image-processing, hypergraphs,
      graphics, machine-learning, fourier-transforms, packed-arrays, geometry,
      bitwise)

## Notes
- Docs-only commits: NO `$VersionNumber` bump, NO tag.
- Never hand-edit `site/docs/documentation/**` — regenerate via `make docs`.
- Overlay In[] must be single-line; `(* note *)` first word must not be an
  identifier (it gets capitalised).
