# ComplexPlot fixes (2026-08-24)

User report:
1. `ComplexPlot[(z^2-1)/(z^2+1), {z,-2-2I,2+2I}]` renders "far too dark".
2. `ComplexPlot[(z^3-3)/z, {z,-2-2I,2+2I}, ColorFunction -> (Hue[#8+0.5]&)]`
   — ColorFunction "does not work"; prints `Power::infy: 1/0`.

Root causes (all in `src/graphics/complexplot.c`):

- [ ] **A. Too dark.** `cp_ramp_color` multiplies the saturated hue by
  `bright = |w|/(1+|w|)`. At |w|=1 that is 0.5, so a bounded function
  (|f|~1 over most of the plane) renders at ~50% brightness everywhere.
  Fix: treat `L = |w|/(1+|w|)` as an HSL lightness — fade toward BLACK for
  L<0.5 (zeros) and toward WHITE for L>0.5 (poles). |w|=1 → full hue.
  Keeps default == `ColorFunction->"Cyclic"` (both go through cp_ramp_color).

- [ ] **B. Spurious `Power::infy`.** Held body evaluated at pole grid points
  (z=0, z=±i) prints `Power::infy` to stderr. Wrap the grid eval in
  `arith_warnings_mute_push/pop()` (already honored by every infy site).

- [ ] **C. ColorFunction gets wrong args.** Custom fn is called `f[re,im]`
  (Re/Im of f only), so `#8` is an unfilled Slot and the fn "does not work".
  Mathematica supplies 8 SCALED args in order:
  #1 Re[z] #2 Im[z] #3 Abs[z] #4 Arg[z] #5 Re[f] #6 Im[f] #7 Abs[f] #8 Arg[f].
  Fix: pass all 8, per-arg scaled to [0,1] when ColorFunctionScaling->True.

- [ ] **D. Compiled poles marked valid.** `cp_eval`'s autocompiled path skips
  the `isfinite` guard, so a compiled pole returns inf and pollutes the
  color-scaling ranges. Add the guard (return false → invalid cell), matching
  the interpreter path.

- [ ] **E. is_color_head too loose.** Accepts `Hue[0.5 + #8]` (symbolic).
  Require all color components be real numbers.

Verify: no `Power::infy`; the Hue ColorFunction plot emits varied concrete
Hue colors; issue-1 plot mean luminance materially higher; tests build+pass.

## Review (done 2026-08-24)

All five fixes landed in `src/graphics/complexplot.c`; docstrings
(`graphics_init.c`), `docs/spec/builtins/graphics.md`, and the weekly changelog
updated. Verified:

- **A. Brightness.** `cp_ramp_color` now treats `L = |w|/(1+|w|)` as HSL
  lightness (black below 1/2, white above). Issue-1 mean cell luminance
  `0.265 → 0.523` (measured by rebuilding the old binary and comparing).
  Default stays identical to `ColorFunction -> "Cyclic"` (`=== True`).
- **B. Power::infy.** Muted via `arith_warnings_mute_push/pop` around the
  sampling loop — both original cases and the pole sweeps print nothing to
  stderr.
- **C. 8-arg ColorFunction.** `(Hue[#8+0.5]&)` renders 159 996 distinct Hue
  colours on the 400² grid; unscaled `#8` is `Arg[f]` in radians, scaled maps
  to [0,1]. Per-arg scaling via new `CFRange`/`grid_cfrange`/`cf_eight`.
- **D.** compiled-pole `isfinite` guard added (parity with interpreter path).
- **E.** `is_color_head` now requires real components.

Tests: `autocompile_tests`, `graphics_tests`, `graphics_sampling_tests` pass.
Memory: differential `valgrind` — new build's "definitely lost" is byte-for-byte
equal to the prior build (14,120 bytes / 433 blocks; all pre-existing baseline).
Clean `-std=c99 -Wall -Wextra` build.
