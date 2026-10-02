# Notebook UX: 70s green theme + evaluate-and-advance (v0.259)

Four requests, all frontend (`frontend/`). No CAS-engine code.

## 1. 1970s green-phosphor syntax scheme (adaptive, bold builtins)
- [ ] `App.svelte`: add `--cm-grn-*` ramp — bright CRT greens in `:root` (dark),
      deep greens in `html.light` (light). Adaptive, mirrors the `--cm-gray-*` ramp.
- [ ] `schemes.ts`: add `seventies` scheme (label "70s (Green Phosphor)"),
      `builtinBold: true`, all `cm` values `var(--cm-grn-*)`. Placed before `eighties`.
- [ ] `menuCommands.ts`: add `scheme-seventies` to `MENU_IDS` + a dispatcher `case`.
- [ ] `src-tauri/src/lib.rs`: add `sc_seventies` CheckMenuItem, submenu entry, vec entry.

## 2. Evaluate-then-advance to next cell (Shift-Enter walks the notebook)
- [ ] `NotebookCard.svelte`: `handleRun` focuses the next cell after evaluating —
      the appended cell when on the last row (existing), else the next existing cell
      via `seekFocusable`. Add `advanceToNextCell`; `continueBelowIfLast` returns the id.

## 3. Shift-Enter dismisses the autocomplete dropdown and evaluates
- [ ] `CellShell.svelte`: `closeCompletion(view)` in the Shift-Enter and Mod-Enter
      key handlers (import from `@codemirror/autocomplete`).

## 4. Release bookkeeping
- [ ] `src/version.h`: 0.258 -> 0.259 (number + string).
- [ ] `docs/spec/changelog/2026-09-28.md`: changelog entry.

## Verify
- [ ] `python3 tools/check_menu_ids.py` (menu wiring gate)
- [ ] `cd frontend && npm run check` (svelte-check + tsc)
- [ ] `cd frontend/src-tauri && cargo check` (Rust compiles)

## Review

All four done and verified (v0.259).

- **70s theme** — `schemes.ts` `seventies` scheme (`builtinBold: true`), adaptive
  `--cm-grn-*` ramp in `App.svelte` (bright on dark, deep on light), wired into
  `menuCommands.ts` (`MENU_IDS` + case) and native menu in `lib.rs`. Placed before
  the `eighties` scheme so the menu reads 70s → 80s. Error stays red (a
  monochrome-green scheme can't signal "wrong" in green).
- **Evaluate-and-advance** — `NotebookCard.handleRun` focuses the next cell after
  eval: appended cell on the last row, else next existing cell via `seekFocusable`.
  Kernel call not awaited, so the caret moves at once. Mod-Enter invariant (one
  cell, no doubles) preserved.
- **Shift-Enter dismisses autocomplete** — `closeCompletion(view)` in the
  Shift-Enter and Mod-Enter handlers (`CellShell.svelte`); no-op when none open.
- **Bookkeeping** — `src/version.h` 0.258→0.259; changelog entry in
  `docs/spec/changelog/2026-09-28.md`.

### Verification
- `python3 tools/check_menu_ids.py` → 64/64/64, wiring complete.
- `npm run check` (svelte-check + tsc) → 0 errors (7 warnings, all pre-existing).
- `cargo check` → compiles (1 pre-existing warning in kernel.rs).
- `npm run build` (vite) → built.

Not committed/tagged — awaiting go-ahead. When ready: commit + `git tag v0.259`.
