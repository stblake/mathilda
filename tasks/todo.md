# Task: Implement genuine Deléglise–Rivat prime counting (2026-06-27)

Replace the `count_dr → count_lmo` alias in `src/numbertheory/primecount.c`
with a real Deléglise–Rivat implementation (easy/hard/trivial special-leaf
split). Highly efficient, extensive tests, no leaks.

## Plan
- [ ] Rewrite `count_dr`: S2 = trivial + easy (PiTable, O(1)) + hard (segmented sieve)
- [ ] S1 ordinary leaves + P2 shared with LMO shape; free-list mirrors every alloc
- [ ] Tests: exact values, non-power-of-ten, per-integer stress vs Lucy oracle,
      DR-vs-LMO agreement across leaf-class-exercising points
- [ ] Clean -std=c99 -Wall -Wextra build; test_prime passes foreground
- [ ] Benchmark DR vs LMO at 1e11..1e13; wire AUTOMATIC to the winner
- [ ] Valgrind count_dr path vs macOS baseline
- [ ] Docs: header comment, number-theory.md, docstring, changelog 2026-06-22.md

## Review
DONE. `count_dr` is now a genuine Deléglise–Rivat implementation (was
`return count_lmo(x)`).
- **Algorithm**: shared Meissel skeleton + identical S1/P2 as LMO; special leaves
  split trivial/easy (O(1) PiTable via the Lehmer prune) vs hard (Fenwick query),
  branched inline within ONE segmented pass (a separate easy pass would
  re-enumerate ~1e8 leaves and be slower). Strictly ≤ LMO cost.
- **Correct by construction** (easy-leaf identity = the prune phi_rec trusts) and
  verified: exact pi(10^7..10^11) + non-power-of-ten; per-integer windows vs
  Lucy/Sieve oracles (off-by-one detector); all-methods SameQ. test_prime passes.
- **Perf**: DR beats Lucy from ~1e10 (1.6s vs 2.5s @1e12; 9.9s vs 13.5s @1e13).
  Automatic rewired: table ≤1e6, Lucy ≤1e9, DR above (old stale TODO removed).
- **Leaks**: valgrind on the DR path byte-identical to Sin[1.0] baseline; no
  prime-counting frames in any leak stack.
- **Docs**: header comment, number-theory.md, changelog 2026-06-22.md; docstring
  already listed the method. Clean -std=c99 -Wall -Wextra.

---

# Task: Implement Frame -> True for Graphics / Plot

WL spec: `Frame` draws a rectangle along the plot edges (vs. the through-origin
`Axes` cross). Settings: `True`, `False`/`None`, or per-edge
`{{left,right},{bottom,top}}`. With `FrameTicks -> Automatic` (default) ticks
are included whenever a frame edge is drawn. Frame ticks + sub-ticks chosen
carefully.

## Plan
- [x] New symbols SYM_Frame / SYM_FrameStyle / SYM_FrameTicks (sym_names.{c,h}).
- [x] render.c: parse Frame / FrameTicks / FrameStyle into per-edge arrays.
- [x] render.c: draw_frame_lines (world space, inside BeginMode2D) — box edges +
      inward major/minor ticks; draw_frame_labels (screen space) — major labels.
- [x] Two-tier ticks: major on nice_step values; minor subdivisions from the
      step's leading digit (1→5, 2→4, 5→5) via frame_minor_divs().
- [x] plot.c: Frame->True withholds the default Axes->True (Frame->False keeps it).
- [x] Docs: graphics.md option tables + prose; changelog 2026-06-22.md.
- [x] Tests: option passthrough, axes-suppression, frame_minor_divs policy.
- [x] Clean -Wall -Wextra build (graphics on), all graphics_tests pass.

## Review

Done — Frame fully supported across Graphics / Show / Plot:

1. **Settings** — `Frame -> True` (all edges), `False`/`None` (none), and the
   per-edge `{{left,right},{bottom,top}}` form. `FrameTicks -> Automatic/None`
   (and per-edge). `FrameStyle -> RGBColor/GrayLevel` colours box+ticks+labels.
2. **Carefully chosen ticks** — major ticks reuse the axes' `nice_step`; minor
   sub-ticks subdivide each major interval, count read off the step's leading
   digit so minors land on round values (1→5, 2→4, 5→5, every magnitude).
   Minor = half-length, unlabeled; majors labelled on bottom/left (fallback
   top/right), placed just inside the frame (it hugs the window border).
3. **Axes interplay** — in Plot, `Frame -> True` suppresses the injected
   `Axes -> True`; bare `Graphics` already defaults Axes off. Both can coexist
   if `Axes -> True` is given explicitly.

Key factoring: tick-subdivision policy extracted to pure, testable
`frame_minor_divs()` in render.{c,h}; the frame, like the axes, is recomputed
against the visible range each frame so it tracks pan/zoom. No camera/scissor
changes — low-risk, mirrors the existing axes code exactly.

Verified: clean -Wall -Wextra build (USE_GRAPHICS); all graphics_tests pass
(incl. 3 new Frame tests + frame_minor_divs policy test); headless REPL renders
Frame variants (per-edge, FrameStyle, FrameTicks->None) without crashing.

## Revision (2026-06-23, per user feedback)
- [x] Plot must not extend outside the frame → reserved margin + camera fits the
      interior region + BeginScissorMode clips the curve. Unframed = full window
      (no change). Box-zoom fits region too.
- [x] ~5% frame padding (margin) per window dim; bottom/left floored larger to
      hold outside labels and the bottom help line.
- [x] Ticks inside the frame; tick labels (numbers) outside, in the margin.
- [x] Bottom help text now outside the frame (lives in the bottom margin).
- [x] Frame box, ticks and label text at 1.5× the axis hairline; added
      hershey_draw_text_ex(thickness) (rounded joins) for the labels.
- Frame is now a fixed screen rectangle; tick *values* update live as the data
  pans/zooms inside it. Screen-space draw_frame replaced the world-space version.

## Follow-ups (not done)
- FrameLabel (edge axis titles) — deferred; would reuse the now-reserved margin.
- Very long numeric labels (≥6 chars) could exceed the left margin floor (50px)
  on small windows; fine for typical ranges and grows with 5% on larger windows.
- Could not verify the live GUI here (headless mode early-returns; a real window
  blocks the session). Verified via clean build + headless option/policy tests;
  needs an eyeball in an interactive run.

# Task: Implement LegendreP

## Scope (faithful, verifiable core)
- [ ] `src/special_functions/legendre.{c,h}` — new module
- [ ] `LegendreP[n, x]`:
  - exact integer n (any sign, P_{-1-n}=P_n) -> exact monomial polynomial (GMP rationals)
  - `LegendreP[n, 1] = 1` (any n)
  - non-integer n with an inexact arg -> numeric Gauss 2F1(-n,n+1;1;(1-x)/2), real+complex, arbitrary precision
  - else symbolic (NULL)
- [ ] `LegendreP[n, m, x]` (type 1, default): integer n,m>=0 -> (-1)^m (1-x^2)^(m/2) D^m P_n(x); m>n -> 0
- [ ] `LegendreP[n, m, a, x]`: a in {1,2,3}; a=1 -> type 1; a=2/3 -> 2F1Reg core * prefactor
- [ ] Attributes: Listable, NumericFunction, Protected
- [ ] sym_names.{c,h}: SYM_LegendreP
- [ ] core.c: include + legendre_init()
- [ ] info.c: docstring
- [ ] tests/CMakeLists.txt: COMMON_SRC + legendre_tests target
- [ ] tests/test_legendre.c: extensive
- [ ] docs/spec/builtins + changelog
- [ ] valgrind clean

## Deferred (documented): symbolic Series/SeriesCoefficient, D rules, non-integer associated/type forms, |w|>=1 analytic continuation.

## Review
DONE. All checkboxes complete.
- legendre.{c,h} implemented; exact integer polynomials via GMP-rational
  three-term recurrence; numeric non-integer order via Gauss 2F1 on ncpx
  (real/complex, machine + MPFR, verified to 300 digits vs WL); associated
  type-1 via Rodrigues derivative; types 2/3 via regularized-2F1 core.
- All WL example outputs reproduced exactly (poly forms, 9.58312, complex
  5.20466+0.299479 I, N[..,50], 300-digit timing, 10/2/x, types).
- Attributes {Listable, NumericFunction, Protected}; docstring; SYM_LegendreP.
- tests/test_legendre.c (13 groups) passes, exit 0, no FAIL.
- valgrind: identical to Sin[1.0] baseline (12,800B/400 — known macOS noise);
  zero LegendreP frames.
- Strict C99 -Wall -Wextra clean.
- Deferred (left symbolic, documented): non-integer associated/type forms,
  negative m, Series/SeriesCoefficient, D[] rules, |(1-x)/2|>=1 continuation.
