# Task: fix the remaining Mathilda ↔ Mathematica divergences (MATHILDA_DIVERGENCES.md)

The doc is a snapshot at build 0.169; HEAD is v0.209. Verified against the live binary:
A1–A10 and A13 were already fixed by commit `3c12c301`. The genuinely-open items are
A15a, A15b, A11, A14, A12. Fix them one-by-one in priority order; reconcile the doc;
and remove the now-redundant section-D `.m` workarounds (corpus-verified).

Full plan: `/Users/user/.claude/plans/let-s-review-the-deficiencies-staged-prism.md`

## Checklist

- [x] Stage 0 — reconciled `MATHILDA_DIVERGENCES.md` (STATUS banner; committed, no bump).
- [x] Stage 1 — A15a: `{} . {}` segfault → scalar `0`. **v0.210, pushed.**
- [x] Stage 2 — A15b: `Coefficient[…, x, i]` symbolic exponent → NULL. **v0.211, pushed.**
- [x] Stage 3 — A11: capture-avoiding `replace_bindings` for scoping constructs. **v0.212, pushed.** (dsolve_tests SIGALRM confirmed pre-existing via A/B.)
- [x] Stage 4 — A14: **A14a DONE** (v0.213, ToNumberField precision escalation). **A14b DONE**
      (v0.217): `qa_resolve_extension` recognises a `Root[]` object (new `qa_resolve_root_object` +
      `GEN_ROOT` autodetect kind in `autodetect_walk`/`autodetect_build_tower`, src/poly/qafactor.c);
      reused `qaext_from_q_expr`/`expr_qx_is_irreducible`. Far smaller than feared — the tower/GCD
      machinery was already generic over the minpoly. Reducible Root defining polys declined to gcd
      over Q. Verified: target case degree-1 gcd, divisibility over Q(c), numeric root, 9 C suites,
      check-c99, no leak (1x==40x). Tests in test_extension_options.c + test_extension_auto_builtins.c.
- [x] Stage 5 — A12: **DONE** — sub-steps 1+3 (v0.214: pf-resolution + Rothstein–Trager mod-d reduction),
      sub-steps 2+4 (v0.215: general linear denominator + numeric `N[RootSum]`). All five documented
      A12 cases evaluate; each verified against the explicit sum over roots.
- [x] Workaround removal: the mod-p/PadRows workarounds in the uncommitted `.m` ALREADY
      delegate to the fixed builtins (user's own edit, "call sites unchanged"). Nothing to
      remove without churning that deliberate design. Verified no regression: parallelmixedtower_tests
      pass; P1/P2/A4 verify; A6/P1 answers byte-identical with/without A11 (A/B). I did NOT edit `.m`.

Each behaviour-changing stage: unit test → build → audits (if numeric) → version bump
+ tag + changelog (`docs/spec/changelog/2026-09-21.md`) + `docs/spec/builtins/` doc-sync.

## Progress log

(filled in as stages land)

## Review

Reconciled MATHILDA_DIVERGENCES.md (a snapshot from build 0.169) against the live
v0.209 binary: A1–A10 and A13 were already fixed by commit 3c12c301. Fixed the
genuinely-open items one-by-one, each verified and shipped (version bump + tag + push
+ changelog + doc-sync + regression test):

- v0.210 — A15a: `{} . {}` segfault → scalar 0 (`src/linalg/dot.c`; R/S as products
  of non-contracted dims, removing the divide-by-K special case).
- v0.211 — A15b: `Coefficient[…, x, i]` with a non-integer exponent stays unevaluated
  (`src/poly/poly.c`).
- v0.212 — A11: capture-avoiding `replace_bindings` (`src/match.c` + new
  `scoping_capture_avoid`/`expr_is_binding_scope` in `src/modular.c`) — scoping locals
  no longer swallow a free parameter; also fixes the mirror shadow case. Highest impact
  (silent-wrong on the integration hot path). dsolve_tests' pre-existing SIGALRM confirmed
  unchanged by isolated A/B.
- v0.213 — A14a: `ToNumberField[a, theta]` escalates membership precision at every degree
  (`src/poly/flint_qqbar.c`), scoped to the builtin so the dsolve-facing gate is untouched.
- v0.214 — A12 (part 1): RootSum general rational body via Rothstein–Trager reduction +
  Function-valued-symbol resolution (`src/root.c`).
- v0.215 — A12 (part 2): RootSum general linear denominator + numeric `N[RootSum]`
  (`src/root.c`, `src/numeric.c`). All five documented A12 cases evaluate; each checked
  against the explicit sum over roots.

Doc-only (no bump): reconciled MATHILDA_DIVERGENCES.md STATUS.

Workaround removal: the section-D mod-p/PadRows workarounds in the uncommitted `.m` already
delegate to the fixed builtins (deliberate, "call sites unchanged"); nothing to remove
without fighting that design. Certificate engine verified regression-free
(parallelmixedtower_tests pass; A6/P1 answers byte-identical with/without A11).

Remaining open: **A14b** — `PolynomialGCD[…, Extension->Automatic]` with a `Root` generator.
Needs a new algebraic-number generator kind threaded through the radical-oriented
autodetect_walk/qa_resolve_extension/flint_extension_gcd pipeline; deferred to a dedicated
effort (large, subsystem-wide risk, worked around in the `.m`). B-series items are
behavioural/by-design.

`make check-c99` passes. Six version tags pushed (v0.210–v0.215).

### Review addendum (2026-09-27, late): the corpus on the v0.216 binary

Re-ran the 371-integral review corpus on v0.216 (the core fixes above post-date the `.m` port).
Raw: 285 correct but THREE false non-elementary certificates (cube-root `x^-1` binomials) and one
120 s timeout where the port used to crash. Diagnosed to three new core divergences, each with a
repro in MATHILDA_DIVERGENCES.md and a `.m` workaround (no C change, nothing committed):

- A16 `PolynomialGCD[x, e, Extension -> Automatic]` = 1 for an unexpanded algebraic constant e = 0
  -> residue classes silently lost -> unsound certificate. Fixed in the `.m` by expanding before the
  gcd AND by making a lost class an honest failure (soundness guard, also in the research WL,
  Python and Maxima ports).
- A17 `Method -> "OneStepRowReduction"` never finishes on a 20 x 15 parametric matrix
  (`mixed/rowreduce_onestep_repro.m`); the default method is instant. `.m`: default method with parameters.
- A18 `Do` iterators inside a package capture a caller's same-named symbol (`Table` does not) --
  the `{a, Length[monos]}` / `{e, ...}` loops of `AnsatzSystem` substituted loop indices for the
  parameters a, e. `.m`: iterators renamed (`ma`, `ent`, `ne`). Core: make `Do` capture-avoiding
  like `Table`; ~150 single-letter iterators in the `.m` remain exposed until then.

Final: 304 correct / 67 gaps / 0 wrong / 0 false certificates (was 279 / 92 / 0 / 0), no regression,
parameters 22/22. parallelmixedtower_tests pass. Changelog entry in docs/spec/changelog/2026-09-21.md;
`mixed/review_mathilda.json` replaced.
