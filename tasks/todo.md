# Remove the `Locked` attribute (`ATTR_LOCKED`)

Full removal — nothing in open-source Mathilda should be permanently locked.
`Remove`'s self-protection survives via `Protected`.

## A. C core — flag + machinery
- [ ] `src/attr.h` — delete `ATTR_LOCKED` define, leave gap comment
- [ ] `src/attr.c` — set_attributes guard, set/clear_attributes_for_symbol guards, string_to_attribute mapping, builtin_attributes count+emit, scan attr_init
- [ ] `src/core.c` — Remove seeding; Unset/ClearAll/Remove masks; Protect/Unprotect early-returns
- [ ] `src/eval.c` — value-list assign, UpSet, TagSet masks
- [ ] `src/options_builtin.c` — delete SetOptions::locked guard
- [ ] `src/options.h` — fix comment
- [ ] `src/sym_names.{c,h}` — remove SYM_Locked decl/def/intern

## B. C docstrings
- [ ] `src/info.c` — SetOptions, ClearAll, Remove, Protect, Unprotect

## C. Authored docs
- [ ] `SPEC.md` attribute table
- [ ] `docs/spec/builtins/assignment-and-rules.md`
- [ ] `docs/spec/builtins/expression-information.md`
- [ ] `docs/spec/changelog/2026-10-05.md` — add removal entry + reconcile same-week
- [ ] `book/chapters/10-internals.tex` caption

## D. Tests
- [ ] `tests/test_clearall_remove_protect.c` — Attributes[Remove] literal
- [ ] `tests/test_upvalues.c` — convert Locked case to Protected

## E. Generated doc pages (hand-edit)
- [ ] `frontend/public/refpages/` assignment-and-rules + expression-information
- [ ] `site/impl/`, `site/overlays/`, `site/docs/documentation/`

## F. Version + tag
- [ ] `src/version.h` bump to 0.301 (tag on commit, when asked)

## Verification
- [ ] Build clean; check-messages / check-c99
- [ ] Final greps zero
- [ ] Unit tests (incl. -DUP_VALUES)
- [ ] REPL behavior checks
- [ ] Rebuild code-review graph

## Review

**Done — `ATTR_LOCKED` fully removed (v0.301).**

C core: flag deleted from `attr.h` (gap comment left, matching the bit-11
precedent); `SYM_Locked` removed from `sym_names.{c,h}`; the `"Locked"` name
mapping + `Attributes[]` count/emit removed from `attr.c`; the three attribute-
setter guards dropped; `Remove` no longer seeded Locked (`core.c`); all combined
`(ATTR_PROTECTED | ATTR_LOCKED)` guards across `core.c`/`eval.c` reduced to
`ATTR_PROTECTED`; `Protect`/`Unprotect` Locked early-returns deleted; the
`SetOptions::locked` guard + message removed (`options_builtin.c`); comments in
`options.h` fixed. Docstrings updated in `info.c`.

Docs: `SPEC.md` attribute table, `docs/spec/builtins/{assignment-and-rules,
expression-information}.md`, `book/chapters/10-internals.tex` caption, and the
current-week changelog (`2026-10-05.md`) with a removal entry + same-week
reconciliation. 19 generated pages under `frontend/public/refpages/` and `site/`
hand-edited. Historical (pre-this-week) changelogs left as dated records.

Tests: `tests/test_clearall_remove_protect.c` (`Attributes[Remove]` →
`{HoldAll, Protected}`) and `tests/test_upvalues.c` (`test_upset_locked` →
`test_upset_protected`). Both suites green.

Verification: clean `make` build; `make check-messages` and `make check-c99`
pass; `grep ATTR_LOCKED\|SYM_Locked src/` → only the intentional gap comment;
repo-wide `\bLocked\b` audit → only intended + coincidental hits. REPL:
`Attributes[Remove]` = `{HoldAll, Protected}`; `SetAttributes[x, Locked]` is an
ignored no-op (`{}`); a user symbol is Protect→Unprotect→Remove-able;
`$VersionNumber` = 0.301.

Not committed/tagged yet (awaiting user request). On commit: tag `v0.301`.

---

# Radical Simplify hang: generalize radrat to constant algebraic radicals

`Simplify[D[Integrate[(x^2+1)/(x^3 Sqrt[2x^4-2x^2+1]),x],x] - integrand]` hangs.
Root cause: `simp_radical_rational` (radrat.c) — the quotient-ring radical
normal-form decision procedure — DECLINES because it only collects radical
bases containing a symbol, so the constant `Sqrt[2]` is dropped and only the
quartic radical remains (n=1 < 2 -> NULL). Input then falls into the algebraic
`Together` (flint_algebraic_field_together) which blows up. Validated: with
Sqrt[2] included as a generator (s^2=2), radrat's exact algorithm reduces the
numerator to 0 instantly.

## A. Core fix — src/simp/radrat.c
- [x] `rr_collect`: collect constant (symbol-free) radical bases as generators;
      compute `any_symbolic` afterward over the collected bases
- [x] gate (~L326): require `n>=2 && any_symbolic` (pure-numeric radicals stay
      with RootReduce/qqbar Phase 0c)
- [x] relation build (~L402): build `g^q - base` for numeric constant bases too;
      skip only bare EXPR_SYMBOL bases
- [x] `rr_finalize`: short-circuit to Integer 0 when num reduces to literal 0
      (skips the radical back-subst + Factor[den] that re-enter the algebraic path)

## A2. Root-cause fix — src/poly/poly.c (`is_target_power`)
- [x] a LITERAL numeric radical base is a ground constant: in the radical case
      (`A_one`), do NOT match its bare form or integer powers (only `Power[c,p/q]`,
      q>1). Strict literal test (NOT is_number, which includes E/Pi); gated on
      A_one so exponential generators (E^x, 2^x) are unaffected. (Found via the
      `Factor[Exp[2x]+2Exp[x]+1]` regression in radical_polyops_tests.)

## B. Verify
- [x] target: `Simplify[diff] -> 0` instant (0.0s; was infinite hang)
- [x] soundness: `diff + 1/x -> 1/x`, `+ Sqrt[q]/x -> Sqrt[..]/x` (not 0);
      `Sqrt[2]/Sqrt[1+x^2] - Sqrt[3]...` stays `(Sqrt[2]-Sqrt[3])/Sqrt[1+x^2]`
- [x] generality: G1/G2/G3 constant-radical + x-radical zeros -> 0
- [x] regression: 22 targeted suites PASS (simplify_hang, fullsimplify(+corpus),
      zero_test, possiblezeroq, field_together_cancel, ratcanon_*, radical_*,
      polynomialreduce, factorlist, trigexp_zero, trigfactor, integrate_*rad, ...)
- [x] radical_polyops regression (Factor exp) found & fixed; re-PASS
- [x] dsolve_tests in-suite SpecialFunctionForm fail = PRE-EXISTING ordering flake
      (main exits rc=142/alarm before reaching it; standalone returns List; all 4
      dsolve_m{5,17,19,55}_stress_tests PASS)
- [x] leaks (macOS): 64 bytes constant at N=5 and N=55 -> no per-call leak
- [x] make check-c99 clean

## C. Tests / docs / version
- [x] test_simplify_hang.c: reported case in termination battery + exact `-> 0`
      check + generality + soundness + pure-numeric guards
- [x] docs/spec/builtins/simplification.md (multi-generator bullet) +
      docs/spec/changelog/2026-10-05.md (v0.302 entry)
- [x] bump src/version.h -> 0.302 (tag v0.302 at commit time — not committed here)

## D. Optional defense-in-depth (NOT done — follow-up)
- [ ] size/iteration budget in flint_tower_reduce / flint_algebraic_field_normalize
      so a radical case radrat still declines terminates (unsimplified) instead of
      hanging. Not needed for this fix (radrat now covers the class); recommended
      hardening because the algebraic-field Together has no internal deadline.

## Review
Root cause was two-layered. (1) `radrat` — the quotient-ring radical normal-form
decision procedure — already existed but excluded constant radicals from its
generator set, so `function-radical + Sqrt[2]` declined (one generator) and fell
into the algebraic-field `Together`, which is exponential over `Q(x)[R]/(R^2-q)`.
(2) Reusing `poly_subst_radical_to_gen` for a numeric base exposed a latent bug in
`is_target_power`: a numeric base matched its bare form and integer powers, so
every coefficient power of 2 (`2 -> g^2`, the `2x^4` of the radicand included) was
swept into the generator, corrupting the relations; the numerator then never
reduced to 0 and the slow algebraic Cancel recovered 0 only after ~3 s. The fix
makes numeric radical bases contribute only genuine radical powers, so the
numerator reduces cleanly to 0 and `rr_finalize` short-circuits. Minimal, sound
(true relations only + strict score gate), general (any function-radical + any
algebraic-constant radicals), leak-free, no regressions.
