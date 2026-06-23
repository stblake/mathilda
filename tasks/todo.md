# Task: Implement AspectRatio option for Plot[]

WL spec settings to support:
- `Automatic` — ratio from actual coordinate values (true geometry)
- `Full` — stretch graphics to fill the enclosing region
- `a` — explicit height-to-width ratio (number, incl. symbolic like 1/GoldenRatio)
- default for Plot: `1/GoldenRatio`

## Current state (already committed)
- plot.c injects `AspectRatio -> 1/GoldenRatio` real default. OK
- render.c honours numeric `a > 0`; `<= 0` sentinel = Automatic. OK
- Gaps: no `Full`; `Automatic`/`Full` symbols not parsed; symbolic numeric
  values (1/GoldenRatio) not resolved on the Plot passthrough; no docstring.

## Plan
- [ ] Add `SYM_Full` to sym_names.{c,h} + intern in sym_names_init.
- [ ] plot.c split_options: handle AspectRatio in the loop — pass
      Automatic/Full verbatim, numericalize anything else to a real.
- [ ] render.c gfx_options_parse: add `aspect_full`; parse Automatic/Full/number;
      resolve Full -> height/width after the option loop.
- [ ] info.c: terse AspectRatio docstring (no examples) + ATTR_PROTECTED.
- [ ] docs/spec/builtins/graphics.md + changelog 2026-06-22.md.
- [ ] Build clean (graphics on), review.

## Review

Done — all four asks addressed:
1. AspectRatio settings: numeric `a` (incl. symbolic 1/GoldenRatio via N[]),
   `Automatic` (data geometry), `Full` (fill ImageSize box). Default 1/GoldenRatio.
2. Window now reshapes to the ratio (was: only y-scale, letterboxed). Default
   plot 800x494; Automatic = data ratio; Full keeps box; ImageSize->{w,h} pins.
3. Docstrings for AspectRatio AND ImageSize in info.c (terse, no examples).
4. Unit tests in tests/test_graphics.c for every setting: option resolution at
   the Graphics level + the pure gfx_window_height() policy (USE_GRAPHICS).

Key factoring: window-height policy extracted to pure, testable
`gfx_window_height()` in render.{c,h}. plot.c numericalizes AspectRatio via N[]
(renderer has no evaluator). New SYM_Full.

Verified: clean -Wall -Wextra build; all graphics_tests + graphics_sampling_tests
pass; window dims confirmed via temporary instrumentation (since removed);
valgrind byte-identical to Sin[1.0] baseline (no new leaks, no src frames).
