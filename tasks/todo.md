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
