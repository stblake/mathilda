---
schema_version: 2
title: Wolfram-style notebook view — grouped toolbar, properties sidebar, split panes
slug: notebook-toolbar
status: in-review
source: direct-user-request
owner: Michael Sollami
issue: pending
pull_request: https://github.com/stblake/mathilda/pull/57
started: 2026-08-13
last_updated: 2026-08-13
blocked_by: none
goal_lock:
  # Released 2026-08-13 22:10: implementation is complete and in review as
  # stblake/mathilda#57. The lock exists to keep an in-flight implementation
  # inside its agreed plan boundary; leaving it active once the work has shipped
  # would gate every unrelated feature in the repo. The follow-up work listed
  # under Follow-Up will be locked by its own tracking file.
  status: released
  stamped: 2026-08-13 17:33
  scope:
    - "frontend/src/**"
    - "frontend/package.json"
  success_criteria:
    - "AC-1 — The focused view shows labelled, vertically-ruled toolbar groups instead of seven unlabeled glyphs."
    - "AC-2 — A [← Canvas] control is always visible and returns to canvas with a mouse alone."
    - "AC-3 — The toolbar knows which cell the caret is in, and clicking a toolbar button does not move or lose the caret."
    - "AC-4 — The Text-or-Code group swaps contents by the active cell's type."
    - "AC-5 — Focused mode tiles 1–4 notebooks horizontally, vertically, or 2×2, with draggable dividers and independently scrolling panes."
    - "AC-6 — Toolbar buttons act on the active pane only."
    - "AC-7 — Abort leaves the kernel alive, and Cmd+. does the same as the menu item."
    - "AC-8 — The single-pane focused view is visually unchanged apart from the toolbar, with no overflow scrollbar."
    - "AC-9 — No control in the shipped toolbar is dead."
    - "AC-10 — Text cells render Markdown, and B/I/U round-trip through a save."
---

# Wolfram-style notebook view — grouped toolbar, properties sidebar, split panes

## Schema

- Schema version: `2`
- ID conventions: `AC-N` acceptance criteria, `NFR-N` non-functional requirements, `RISK-N` risk-register rows.

## Tracking Metadata

- Title: `Wolfram-style notebook view — grouped toolbar, properties sidebar, split panes`
- Slug: `notebook-toolbar`
- Status: `in-review`
- Source: `direct-user-request`
- Owner: `Michael Sollami`
- Issue / ticket: `pending`
- Pull request: `stblake/mathilda#57 — FindClusters in n dimensions, runnable reference pages, and a Wolfram-style notebook toolbar`
- Started: `2026-08-13`
- Last updated: `2026-08-13`
- Blocked by: `none`

## Feature Definition

- One-line goal: `Bring the single-notebook (focused) view up to Wolfram 14's grouped-toolbar legibility, add a properties sidebar, and let 1–4 notebooks tile side by side for simultaneous editing and execution.`
- Problem: `The focused view is a 34px bar of seven unlabeled Unicode glyphs (App.svelte:158-197) with no grouping, no labels, no formatting, no evaluation control beyond "run all", and no way to see two notebooks at once. Cell type is buried behind a 12px badge in the cell gutter (CellShell.svelte:255).`
- Requested by: `direct user request, with three Wolfram 14 toolbar screenshots as the reference`
- Related links:
  - `/Users/67840/.claude/plans/floofy-percolating-walrus.md` (approved plan)

## Working Description

Replace the focused-mode app bar with a row of labelled, vertically-ruled groups
matching Wolfram 14's arrangement, one of which is context-sensitive on the active
cell's type. Add a left-hand Notebook Properties panel. Convert focused mode from
"exactly one notebook fills the window" to "1–4 notebooks tiled h / v / 2×2 with
draggable dividers and independent scrolling", where the toolbar drives whichever
pane is active.

Four confirmed decisions shaped the design:

1. Text cells become **Markdown in `source`** — B/I/U wrap the selection, rendered
   with `marked` (already a dependency), raw on focus and rendered on blur. Chosen
   over storing HTML because it round-trips through `serialize()` unchanged and
   gains headings/lists/links for free.
2. Evaluation stays **frontend-only** — no Rust changes and no faked parallelism.
3. The toolbar **replaces** the app-bar contents in focused mode (one bar, not two).
4. Panes come from a toolbar control and **share one kernel**.

## Product Requirements

### Business requirement

A user editing Mathilda notebooks can see and act on every notebook-level command
without hunting through unlabeled glyphs, and can work two notebooks side by side
in one window.

### User stories

- As a notebook author, I want toolbar controls grouped under labels so that I can
  find a command without hovering every icon to read its tooltip.
- As a notebook author, I want the toolbar to reflect the cell I am editing so that
  formatting and cell commands apply to what I have selected, with zero misfires.
- As a notebook author, I want two notebooks side by side so that I can copy a
  result from one into the other without leaving the focused view.
- As a notebook author, I want to abort a runaway evaluation so that I regain
  control, and I want to be told plainly if that costs me my kernel state.
- As a notebook author, I want bold and italic in text cells so that prose survives
  a save with its emphasis intact.

### Actors

`Toolbar` (App.svelte + Toolbar.svelte), `Focused view` (Canvas.svelte),
`NotebookCard`, `CellShell`, `active-cell store` (active.ts), `canvas store`
(canvas.ts), `Tauri kernel bridge` (ipc.ts).

### Functional requirements — basic

| ID | Requirement |
|---|---|
| REQ-001 | The toolbar shall render labelled, vertically-ruled groups in focused mode. |
| REQ-002 | The toolbar shall always render a canvas-exit control in focused mode. |
| REQ-003 | The system shall record which cell holds the caret, and which notebook owns it. |
| REQ-004 | The context-sensitive group shall render Text controls for prose cells and Code controls for code cells. |
| REQ-005 | The focused view shall tile 1–4 notebooks in horizontal, vertical, or 2×2 arrangements. |
| REQ-006 | Each pane shall scroll independently of the others. |
| REQ-007 | Toolbar commands shall act on the active pane only. |
| REQ-008 | Abort shall stop the running evaluation and leave a live kernel. |
| REQ-009 | Text cells shall store and render Markdown. |

### Functional requirements — failure paths

- `REQ-001-R: When the viewport is too narrow for every group, the Toolbar shall move lower-priority controls into the … overflow menu rather than clipping or wrapping.`
- `REQ-002-R: When any other toolbar control is unavailable, the canvas-exit control shall remain enabled — focused mode must never become inescapable.`
- `REQ-003-R: When the recorded active cell no longer exists (deleted, notebook closed, library reloaded), the Toolbar shall resolve it as absent on read and disable dependent controls, without erroring.`
- `REQ-004-R: When there is no active cell, the context-sensitive group shall be hidden rather than rendered with disabled controls.`
- `REQ-005-R: When a notebook in a pane is closed elsewhere, the focused view shall drop that pane and re-equalize the survivors; when the last pane goes, it shall return to canvas.`
- `REQ-006-R: When a pane is narrower than its content's intrinsic width, the pane shall scroll rather than push the grid past its container.`
- `REQ-007-R: When the active pane's action registry is momentarily absent (between activation and the card's next flush), toolbar clicks shall no-op rather than throw.`
- `REQ-008-R: When the abort kills the kernel process, the system shall respawn it automatically and label the control so the cost is stated before the click.`
- `REQ-009-R: When a text cell holds Markdown that fails to parse, the cell shall render its raw source rather than blanking.`

### Non-functional requirements

Graduated to `NFR-N` under Success Criteria below.

### Open points

- `none` — the three forks (rich-text model, kernel depth, split mechanism) were resolved with the user before planning.

## Current State Study

- Relevant existing files:
  - `frontend/src/App.svelte` — the app bar (`:158-197`), both root palettes (`:209-241`), the `menu:interrupt` no-op (`:64`)
  - `frontend/src/lib/Canvas.svelte` — focused-mode render (`:561-572`), `.focused-view` CSS (`:737-755`)
  - `frontend/src/lib/canvas.ts` — `canvasState` (`:84-95`), `FocusedActions` (`:70-82`), `openRefpage` (`:325-433`)
  - `frontend/src/lib/NotebookCard.svelte` — the actions publish block (`:599-615`), `runAll`/`runCell` (`:446-477`)
  - `frontend/src/lib/CellShell.svelte` — CodeMirror setup (`:83-138`), type picker (`:254-268`), prose cells (`:304-347`)
  - `frontend/src/lib/notebook.ts` — `selectedCells` (`:54`), `setCellType` (`:183-188`), `serialize` (`:247-251`)
  - `frontend/src/lib/refpages.ts` — the Cmd+click/F1 symbol lookup (`:217-348`)
  - `frontend/src-tauri/src/kernel.rs` — single kernel, `interrupt` (`:195-199`)
- Adjacent modules touched:
  - `frontend/src/app.css` (height and typographic vars), `frontend/src/lib/platform.ts` (`isTouchDevice`)
- Existing behavior to preserve:
  - Canvas mode in full — pan, zoom, rubber-band select, group drag, minimap
  - `serializeLibrary` / `loadLibraryData` round-trip
  - `openRefpage` and `addQueryNotebook` card creation
  - The single-pane focused view's appearance
  - `setFocused(id | null)`'s signature, called from three pane-unaware sites
- Constraints from the current codebase:
  - Svelte 4 idiom throughout (`export let`, `$:`) — not runes
  - `selectedCells` and `kernelStatus` are module-global singletons shared by all nine notebook stores
  - `.nb-card` has `backdrop-filter`, which creates a containing block for fixed-position descendants
  - `App.svelte:130` sets root `font-size` for Cmd+= zoom, so rem sizing inside a fixed-height bar overflows
- Existing tests / validation paths:
  - `npx svelte-check` — no frontend unit-test suite exists; validation is svelte-check plus manual gates

## Implementation Spec

- New files to create:
  - `frontend/src/lib/active.ts` — active-cell store + cell-handle registry
  - `frontend/src/lib/Toolbar.svelte` — composes groups, owns menu-open state
  - `frontend/src/lib/ToolbarGroup.svelte` — one labelled, ruled group
  - `frontend/src/lib/Menu.svelte` — the single popover primitive
  - `frontend/src/lib/Icon.svelte` — name → inline SVG
  - `frontend/src/lib/PropertiesPanel.svelte` — the left sidebar panel
  - `frontend/src/lib/cellCommands.ts` — pure split / merge / duplicate / delete / wrap / indent / comment
- Existing files to update:
  - `frontend/src/lib/canvas.ts`, `frontend/src/lib/Canvas.svelte`, `frontend/src/App.svelte`,
    `frontend/src/lib/NotebookCard.svelte`, `frontend/src/lib/CellShell.svelte`,
    `frontend/src/lib/refpages.ts`, `frontend/src/app.css`
- Files to delete:
  - `frontend/src/lib/CodeCell.svelte` — 318 lines, zero importers, superseded by CellShell
- Data models / contracts:
  - `ActiveCell { notebookId, cellId, cellType, focused }` — global, single record
  - `CellHandle { view, el, focus }` — kept in a plain `Map`, never a store
  - `canvasState` focus fields: `focusedIds[]`, `focusedLayout`, `focusedActiveId`, `focusedSizes[]`, `focusedGrid{x,y}` replacing `focusedId`
  - `PaneActions` = `FocusedActions` + `notebookId`, `store`, `runCell`, `runRange`, `focusCell`; held in `panes: Map<string, PaneActions>` with a derived `activeActions`
- Import directions:
  - `active.ts` imports only a type from `notebook.ts` (leaf module, no cycle)
  - `Toolbar.svelte` imports `active.ts`, `canvas.ts`, `cellCommands.ts`, `Menu.svelte`, `Icon.svelte`, `ToolbarGroup.svelte`
  - `CellShell.svelte` imports `active.ts`; `NotebookCard.svelte` passes `notebookId` down
  - `Menu.svelte` renders only from `Toolbar.svelte`, never from inside a card
- Execution flow:
  - Caret enters a cell → `updateListener`'s `focusChanged` (code) or `on:focus` (prose) → `setActiveCell`
  - Toolbar reads → validate-on-read re-resolves the cell from the active pane's store
  - Toolbar click → `pointerdown|preventDefault` keeps the caret → handler reads `activeActions` / `activeHandle()`
  - Pane click or focusin → `setFocusedActive` → `activeActions` re-derives → toolbar retargets

## Scope

### In Scope

- Grouped, labelled toolbar replacing the app-bar contents in focused mode
- Active-cell tracking, with the blur problem solved
- Groups: canvas-exit, Sidebar, Evaluation, Cell Style, Cells, Text-or-Code, Insert, Notebook, overflow
- Notebook Properties sidebar panel
- 1–4 pane tiling (h / v / 2×2), draggable dividers, independent scrolling, active-pane routing
- Markdown text cells with inline formatting commands
- Evaluation control: run cell / notebook / from here / above, restart, abort-with-respawn
- Notebook-wide search; math template palette; inline TeX; hyperlink
- Fixing `--surface-2`, the `menu:interrupt` no-op, the two unlistened menu items, the
  `min-height:100vh` overflow, the light-mode type picker, and the missing `commentTokens`
- Deleting `CodeCell.svelte` and de-duplicating symbol-under-cursor into `symbolAtSelection()`

### Out Of Scope

- Real parallel kernels — needs `MathildaKernel` restructured into a keyed map plus a
  kernel id threaded through `evaluate_cell` and `ipc.ts` (`kernel.rs:14-24`)
- Cooperative abort in the C evaluator (`src/eval.c`)
- `Out[n]` rendering, and therefore its toggle — no such markup exists
- Sketch / freehand insert — `serialize()` is `{type, source}`; a drawing could not persist
- Publish / cloud export — no target or backend
- Per-cell styles and themes — `Cell` has no style field
- Fixing the two documented global-singleton hazards (`selectedCells` shared across
  stores; the two conflicting root palettes) — documented, not repaired
- Canvas-mode changes beyond what the focus refactor requires
- Dragging cards together on canvas to form a split (rejected in favour of a toolbar control)

## Success Criteria

- `AC-1` — In focused mode the toolbar renders labelled, vertically-ruled groups; every control has a visible group caption above it.
- `AC-2` — A `[← Canvas]` control is visible at all times in focused mode and returns to canvas using a mouse alone, with no trackpad gesture and no keyboard.
- `AC-3` — Clicking into a code cell makes Cell Style read "Code"; clicking into a text cell makes it read "Text"; and clicking any toolbar button leaves the caret where it was and applies the command to that cell.
- `AC-4` — The Text-or-Code group renders B/I/U for prose cells, indent/outdent/comment for code cells, and is hidden entirely when no cell is active.
- `AC-5` — Focused mode tiles 1–4 notebooks in horizontal, vertical, and 2×2 arrangements; each pane scrolls independently; dividers drag without scrolling the page or selecting text.
- `AC-6` — With two panes open, the Run control runs the active pane's notebook and not the other's.
- `AC-7` — Abort during a running evaluation stops it and leaves a live kernel that evaluates the next cell successfully; `Cmd+.` behaves identically to the menu item.
- `AC-8` — With one pane, the focused view is visually unchanged from before this work apart from the toolbar itself, and shows no overflow scrollbar.
- `AC-9` — Every control in the shipped toolbar performs a real action; nothing is rendered that cannot work.
- `AC-10` — Text cells render Markdown; bold and italic applied via the toolbar survive a save-quit-reopen cycle.

### Non-functional

- `NFR-1` — The toolbar shall not overflow or clip its controls at UI scales from 0.5× to 2.0× (`App.svelte:130`'s Cmd+= range).
- `NFR-2` — Dragging a pane divider shall write the canvas store at most once per drag (on pointerup), not per pointermove, so a 60 Hz drag causes 1 store write rather than ~60.
- `NFR-3` — A keystroke in a cell shall not invalidate the whole toolbar: the per-pane action registry's stable methods are registered once, with volatile flags in a separate store.
- `NFR-4` — `npx svelte-check --threshold error` shall report 0 errors at every phase boundary.
- `NFR-5` — Both light and dark themes shall render every toolbar caption, rule, and menu legibly — no undefined custom property left without a fallback.

## Tests

- Must-have tests:
  - `npx svelte-check --threshold error` clean at every phase boundary
  - Manual gate 1 — single-pane no-regression, no scrollbar
  - Manual gate 2 — escapability with a mouse alone, re-run every commit
  - Manual gate 3 — active cell tracks the caret; toolbar clicks do not steal it
  - Manual gate 4 — split: independent scroll, divider drag, per-pane run routing
  - Manual gate 5 — abort leaves a live kernel; two-pane runs serialize without error
  - Manual gate 6 — save / quit / reopen round-trip, including Markdown formatting
  - Manual gate 7 — light mode legibility
  - Manual gate 8 — Cmd+= to 2.0× without toolbar overflow
  - Manual gate 9 — retyping a code cell with output states its output-clearing behaviour
- Areas requiring full coverage:
  - Active-cell lifecycle: focus, blur, cell deleted, notebook closed, library reloaded
  - Focus-state invariants: `normalizeFocus` under pane removal, duplicate ids, layout demotion
- Out of scope for testing:
  - Automated unit tests for Svelte components — the repo has no frontend test runner, and
    standing one up (vitest + @testing-library/svelte) is its own piece of work. `identifierAt`
    is extracted as a pure function so it is unit-testable when a runner exists.
  - Real parallel evaluation — out of scope as a feature
- Test framework / runner:
  - `svelte-check` for static verification; manual gates for behaviour. No frontend unit runner exists.

## Task List

- [ ] `P0 — foundations: palette, --toolbar-h, active.ts, Icon, Menu, CellShell focus hooks + handle registry + notebookId + commentTokens, delete CodeCell, extract symbolAtSelection, fix menu:interrupt + two unlistened items, light-mode type picker` | `independent` | `pending`
- [ ] `P1 — focus-state refactor: five canvasState fields, normalizeFocus, initialFocus, setFocused shim, panes Map + activeActions, volatile-flag split, all read sites` | `depends on: P0` | `pending`
- [ ] `P2 — Toolbar + ToolbarGroup replace the app bar; [← Canvas], Cell Style combo, Notebook group, overflow; height contract and min-height fixes` | `depends on: P1` | `pending`
- [ ] `P3 — two-pane horizontal split: grid container, keyed each, pane scroller CSS, pane header, activation, divider drag` | `depends on: P2` | `pending`
- [ ] `P4 — Evaluation group: run + caret menu, kernel chip + menu, abort-with-respawn` | `depends on: P3` | `pending`
- [ ] `P5 — Cells group + cellCommands.ts: split, merge, duplicate, delete` | `depends on: P4` | `pending`
- [ ] `P6 — layouts h/v/grid, add-pane picker, openRefpage and addQueryNotebook in split mode` | `depends on: P3` | `pending`
- [ ] `P7 — PropertiesPanel + Sidebar group` | `depends on: P2` | `pending`
- [ ] `P8 — Code sub-group: context-sensitivity mechanism, indent, outdent, comment` | `depends on: P2` | `pending`
- [ ] `P9 — Markdown text cells + Text sub-group (one phase, together)` | `depends on: P8` | `pending`
- [ ] `P10 — Insert group + notebook search` | `depends on: P9` | `pending`

## Test Results

- Command: `pending`
- Outcome: `pending`
- Summary: `pending`
- Failures fixed:
  - `pending`
- Known exceptions:
  - `pending`

## Checkpoints

- [x] start | completed: `2026-08-13 17:33`
- [x] spec / plan created | completed: `2026-08-13 17:33`
- [ ] threat-model stamped | completed: `pending`
- [x] implementation started | completed: `2026-08-13 17:40`
- [x] implementation complete | completed: `2026-08-13 21:14`
- [x] critic pass | completed: `2026-08-13 21:32`
- [x] risk-register reviewed | completed: `2026-08-13 21:32`
- [ ] feature validated | completed: `pending`
- [x] PR created | completed: `2026-08-13 21:34`
- [ ] closeout complete | completed: `pending`

## PR Updates

- `2026-08-13 21:34` Opened stblake/mathilda#57 from msollami:feat/findclusters-ndim against main. 15 commits: 7 pre-existing (FindClusters n-dim, reference pages), 8 from this feature. Scoped as one PR at the user's decision — the toolbar commits build on the reference-page ones and touch the same eight frontend files, so splitting would have meant conflict resolution across all of them.

## Decisions

- `Text cells store Markdown in source rather than HTML — it round-trips through the existing serialize() ({type, source}) with no format change, reuses marked which is already a dependency for RefPage, and gains headings/lists/links for free. HTML would have needed sanitizing and would put markup in .lb files.`
- `Evaluation stays frontend-only — real parallelism needs MathildaKernel restructured into a keyed kernel map with an id threaded through evaluate_cell and ipc.ts. Deferred rather than faked, because a Parallel group that does not parallelize is worse than its absence.`
- `Abort is interruptKernel() then restartKernel(), labelled "Abort (restarts kernel)" — interrupt_kernel is literally self.kill() with no respawn (kernel.rs:197-199), so an unlabelled Abort would silently destroy all definitions. The honest fix here is the label, not the cut.`
- `No frontend queue for multi-pane evaluation — kernel.rs:103-105 already holds a tokio Mutex across the whole request/response loop, so concurrent evaluate_cell calls serialize FIFO in Rust. A frontend queue would only add latency in front of the existing one.`
- `Three enumerated pane layouts (h / v / grid) instead of a binary split tree — with n ≤ 4, grid is the only non-linear arrangement, and the three collapse to two CSS grid templates plus a spanning rule. A tree would buy arbitrary nesting no requirement asks for at the cost of recursive components, size arithmetic, delete-collapse, and serialization.`
- `Two height vars (--appbar-h 34px for canvas, --toolbar-h 46px for focused) instead of making --appbar-h reactive — a reactive var would put a JS style write in the middle of the {#if} branch swap between .canvas-stage and .focused-view, and it keeps .canvas-stage out of the change set entirely.`
- `setFocused(id | null) is kept as a shim over focusedIds rather than removed — it is called from three sites that should not know panes exist (NotebookCard.svelte:685, App.svelte:179, Canvas.svelte:241), making it the largest available back-compat lever.`
- `A notebook may appear in at most one pane — not for the Map key, but because selectedCells and kernelStatus are module-global singletons, so two cards over one store would share cell selection, double-register cellFocusFns, and open two refpages from one Cmd+click.`
- `The blur problem is solved primarily with pointerdown|preventDefault on toolbar buttons, not by keeping a sticky record alone — suppressing the default keeps the editor from ever blurring, so the live text selection survives, which is what B/I/U and comment-toggle actually need. The sticky record plus validate-on-read is the second layer, for controls that must take focus.`
- `Stale active-cell ids are handled by validate-on-read rather than cleanup hooks — re-resolving the cell from the pane's store self-heals across cell deletion, notebook close, and loadLibraryData replacing every store, with zero teardown code.`
- `Icons are inline SVG with stroke="currentColor" rather than Unicode glyphs — Unicode has no acceptable glyph for split-cell or comment-toggle, and glyph metrics vary across the four WebViews this ships to, which is why the current glyph row is visually ragged. Typographic marks (B, I, U, f[], …) stay as text.`
- `The In-label toggle moves from the Cells group to the properties panel — it is a display preference, not a cell action, and it needs a per-notebook flag threaded to every CellShell.`
- `Five controls cut as unbuildable rather than shipped dead: the Out-label toggle (no Out[n] markup exists), sketch and publish (serialize() is {type, source}; no export target), the cell-style swatch (Cell has no style field), and the text size caret (Markdown has no font size, and the honest mapping is the Cell Style combo).`

## Risks And Unknowns

- `The [← Canvas] control is load-bearing: setFocused(null) at App.svelte:178-179 is currently the only mouse-reachable exit from focused mode, the alternative being a pinch-out gesture. If it is ever moved into the overflow or disabled, focused mode becomes a trap. Re-verify at every commit.`
- `canvas.ts:587 (loadLibraryData) is a full-object canvasState.set — omitting any one of the five new focus fields makes it undefined and throws on the next render. Mitigated with a shared initialFocus() literal used by both it and the module-level writable.`
- `The Svelte 4 reactive block at NotebookCard.svelte:599-615 derives its dependencies from the identifiers it syntactically references. Extracting the published object into a helper function would silently stop the block re-running and freeze the toolbar's icon states.`
- `46px × several labelled groups will not fit a 390px viewport. Phase 2 uses isTouchDevice plus a hardcoded priority order rather than a measured ResizeObserver collapse — deterministic but approximate.`
- `Two conflicting root palettes (app.css:1-34 vs App.svelte:209-241) disagree on --text/--bg/--border/--accent, App.svelte winning by load order. New vars are routed to the block that owns their category, but this remains a landmine for anyone adding a var later.`
- `setCellType (notebook.ts:186) clears output and execIdx. The Cell Style combo puts this one click away instead of behind a small badge, so the output-clearing behaviour needs deliberate handling.`
- `Whether Canvas.svelte:241's pinch-out-to-canvas handler still receives wheel events in focused mode is unverified — .focused-view is position:fixed above .canvas-stage, where on:wheel is bound. If it is already broken, the split-pane work must not be blamed for it.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Denial of service (self-inflicted, STRIDE-D) | med | high | Abort hard-kills the kernel with no respawn (`kernel.rs:197-199`), destroying all user definitions. Mitigated by pairing it with `restartKernel()` and labelling the control "Abort (restarts kernel)" so the cost is stated before the click — see `AC-7`. | User still loses in-session definitions on abort; a cooperative abort in the C evaluator is the real fix and is out of scope. | Michael Sollami |
| RISK-2 | Availability / usability lockout | med | high | Losing the only mouse-reachable exit from focused mode traps the user (`App.svelte:178-179`). Mitigated by making `[← Canvas]` always visible, never in the overflow, never disabled — `AC-2`, re-verified at every commit. | None if `AC-2` holds. | Michael Sollami |
| RISK-3 | Tampering (data integrity) | low | med | Markdown text cells render through `marked` with `{@html}`, as `RefPage.svelte:99` already does. Notebook content is user-authored and local, but a `.lb` library file received from elsewhere becomes an injection vector. Mitigation: treat loaded Markdown as untrusted and restrict the rendered surface (no raw HTML pass-through beyond the `<u>` tag the U control emits). | Residual risk equals the pre-existing `RefPage` surface; not made worse. | Michael Sollami |
| RISK-4 | Information disclosure | low | low | The properties panel surfaces the library's absolute filesystem path (`libraryPath`). Local-only desktop app, shown to the user who opened the file. Accepted. | Accepted — no remote surface. | Michael Sollami |

## Dependencies / Blockers

- Dependencies:
  - `marked` (present), `katex` (present), `@codemirror/commands` (present)
  - `@codemirror/autocomplete` and `@codemirror/search` — currently transitive in `node_modules`; promote to `package.json` if P10 uses them, since relying on a hoisting accident will break
- Blockers:
  - `none`

## Proofs Of Completion

- Completion timestamp: `pending`
- Verification commands and results:
  - `pending`
- Verification exceptions:
  - `pending`
- Artifacts:
  - `pending`
- Review / verification links:
  - `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `pending`
- Existing tech debt noticed:
  - `CodeCell.svelte — 318 lines, zero importers (being deleted in P0)`
  - `KernelStatus.svelte — imported at App.svelte:13 and never rendered (being used in P4)`
  - `.out-label CSS at CellShell.svelte:561-567 — orphaned, no matching markup`
  - `--surface-2 referenced with no fallback at App.svelte:302 and defined nowhere; --ok never defined`
  - `Two conflicting root palettes (app.css vs App.svelte)`
  - `selectedCells is a module-global singleton shared by all nine notebook stores`
  - `Canvas.svelte:605's on:focusNotebook handler is dead — NotebookCard never dispatches it`
  - `menu:run-all and menu:add-cell registered in lib.rs:50,67 with no listener; menu:interrupt listened but a no-op`
  - `@codemirror/lang-javascript in package.json, imported nowhere`
  - `min-height:100vh inside an --appbar-h-inset .focused-view guarantees an overflow scrollbar`
- Mitigations taken:
  - `pending`
- Follow-up needed:
  - `pending`

## Activity Log

- `2026-08-13 17:33` Tracking file created. Three Wolfram 14 screenshots supplied as the reference for the target UI.
- `2026-08-13 17:33` Research complete: three parallel Explore agents mapped the kernel IPC surface, the canvas/focus architecture, and cell/output rendering. Three forks surfaced and were resolved with the user — Markdown for rich text, frontend-only evaluation, toolbar-driven split with a shared kernel.
- `2026-08-13 17:33` Two parallel Plan agents designed the split-pane focus restructuring and the toolbar plus active-cell tracking. Both returned corrections: a flat pane-size array cannot express a 2×2; `toggleComment` exists but is a silent no-op without `commentTokens` language data; pane queueing needs no work because Rust already serializes on a mutex; and `[← Canvas]` is a hard requirement because the current exit disappears when the app bar is replaced.
- `2026-08-13 17:33` Plan approved by the user. Stamped `start` and `spec / plan created`; `goal_lock` set active over `frontend/src/**` and `frontend/package.json`. 11 phases created as tasks.
- `2026-08-13 17:40` Implementation started. P0 foundations and P1 focus-state refactor landed together in `7f6fdfa7`, verified by 34 assertions over the focus logic (bundled for node with vite; harness kept out of the repo since a frontend test runner is out of scope).
- `2026-08-13 18:02` P2 toolbar shell landed in `cda8e262`. Found and fixed a pre-existing sizing bug: the focused card asked for `min-height: 100vh` inside a view inset by the bar height, guaranteeing a scrollbar. A percentage min-height cannot replace it — percentages resolve against the parent's *height*, which stays `auto` however many min-heights are stacked — so the view became a flex column. Confirmed by sampling rendered pixels.
- `2026-08-13 20:15` **Scope added at the user's request** after seeing the toolbar in use: an optional status bar (kernel state, last-op timing, session totals, active cell, notebook shape) and a fix for the info button, which was disabled whenever no documented symbol sat under the caret and so read as broken. Required new evaluation-timing instrumentation — nothing measured duration before. Landed with the Evaluation group in `d7fcdc78`.
- `2026-08-13 21:05` **Scope added at the user's request**: make splitting easy from the canvas as well as the toolbar, and open documentation as a side-by-side pane on the right. Landed with the tiled pane view in `d261ea8d`. `openRefpage` was never broken — it positioned its card on the canvas, which is off-screen in focused mode, so Cmd+click appeared to do nothing.
- `2026-08-13 21:14` P5 Cells and P8 Code groups landed in `0c876515`, together with a fix for `setCellType` silently discarding a cell's output on retype — tolerable behind a 12px gutter badge, not with the cell-style control one click away.
- `2026-08-13 21:32` Critic pass: a high-effort code review over the whole diff returned six findings, all real and all fixed in `8714c825` — one HIGH (retyping a cell to Code through the toolbar produced a dead, uneditable cell, because only the gutter badge rebuilt the editor), three MEDIUM (`focusCell` fired before the cell existed so Split/Duplicate/Insert never moved the caret; the pane's first flags publish was discarded so a restored notebook's toolbar had no flags; documentation lookup forced up to five synchronous layouts per keystroke), and two LOW (group rules broke while a menu was open; the caret probe could strand an undocumented name and suppress its fallback). The HIGH fix was verified at runtime.
- `2026-08-13 21:19` Implementation complete for this PR's scope. P7 properties panel, P9 Markdown text cells with the Text sub-group, and P10 Insert plus notebook search are deliberately deferred to a follow-up; see Follow-Up.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- Deferred from this PR, in the order they should land: `P7` Notebook Properties
  panel and its Sidebar group; `P9` Markdown text cells plus the Text half of the
  context-sensitive group (one phase — until cells render Markdown, a bold button
  inserts literal `**` and leaves the asterisks on screen); `P10` Insert group
  (math template palette via `@codemirror/autocomplete`'s `snippet()`, inline TeX,
  hyperlink) and notebook-wide search.
- Interactive verification of the toolbar was not possible in the authoring
  environment: synthetic mouse events do not reach the WKWebView without
  input-synthesis permission, and `view.focus()` on a background window produces no
  focus event. Rendering, disabled states and group visibility were verified from
  screenshots; the enabled paths need a human click-through. Worth confirming
  explicitly: that clicking a toolbar button leaves the caret in the cell, and that
  `[← Canvas]` returns to the canvas with a mouse alone.
- `--surface-2` and `--ok` were referenced but undefined and are now defined, but the
  two conflicting root palettes (`app.css:1-34` vs `App.svelte`'s `:global(:root)`,
  disagreeing on `--text`/`--bg`/`--border`/`--accent`, App.svelte winning by load
  order) remain. It shows as a tonal step between a cell row and the empty card
  space below it. Worth its own cleanup.
- `selectedCells` in `notebook.ts` is a module-global singleton shared by every
  notebook store. Currently masked because a notebook may occupy only one pane —
  that constraint is doing load-bearing work, and should be removed properly if
  duplicate panes are ever wanted.
- `KernelStatus.svelte` is still unused and hardcodes light-mode colours. Either
  make it theme-aware and use it or delete it; `StatusBar.svelte` renders its own
  kernel indicator.

## Team Addendum

- `none`
