# Mathilda Notebook — Desktop Frontend

Next-generation CAS notebook UI for [Mathilda](../README.md), built with
**Tauri 2 + Svelte 5 + TypeScript**. Communicates with the Mathilda C binary
over stdio using a line-delimited NDJSON protocol.

## Stack

| Layer | Technology |
|---|---|
| Desktop shell | Tauri 2 (Rust) |
| Frontend | Svelte 5 + TypeScript + Vite |
| Code editor | CodeMirror 6 |
| Math rendering | KaTeX |
| Plot rendering | Plotly.js |
| IPC | stdio pipes, NDJSON |

## Prerequisites

- Rust toolchain (`rustup`)
- Node.js 18+
- Tauri CLI: `cargo install tauri-cli`
- Mathilda build deps: GMP, GNU Readline (`brew install gmp readline`)

## Dev Setup

```bash
# 1. Build the Mathilda C binary and install as Tauri sidecar
./build-sidecar.sh

# 2. Install JS dependencies
npm install

# 3. Launch the app in dev mode (hot-reload)
cargo tauri dev
```

## Production Build

```bash
./build-sidecar.sh
cargo tauri build
```

The bundled `.app` (macOS) / `.deb` / `.msi` will be in `src-tauri/target/release/bundle/`.

## File Format (.mathilda)

Notebooks are plain-text files with no stored outputs. Each cell is a stanza:

```
(* cell: code *)
Integrate[x^2, {x, 0, 1}]

(* cell: code *)
Factor[x^4 - 1]
```

This format is fully Git-diffable and can also be loaded directly into the
Mathilda terminal REPL (`./Mathilda < notebook.mathilda`).

## Keyboard Shortcuts

| Shortcut | Action |
|---|---|
| Shift+Enter | Run current cell |
| Ctrl/Cmd+Enter | Run current cell and insert new cell below |

## Architecture

```
Svelte UI (src/)
    +  @tauri-apps/api invoke + Channel
Tauri Rust layer (src-tauri/src/)
    +  stdio pipes (NDJSON)
Mathilda C binary (../Mathilda)
```

See docs/frontend-research.md for the full design rationale.

### Focused-mode surfaces

The window has two modes. On the canvas, the top strip is a 34px name-and-theme
bar (`--appbar-h`). Focusing one or more notebooks swaps it for the 46px notebook
toolbar (`--toolbar-h`), and three surfaces then belong to that mode only:

| File | What it owns |
|------|--------------|
| `lib/Toolbar.svelte` | the labelled, ruled control groups; one `Menu.svelte` instance is shared and its `items` swapped |
| `lib/StatusBar.svelte` | the optional 22px bottom strip: kernel state, last evaluation time, session totals |
| `lib/PropertiesPanel.svelte` | the sidebar that slides in from the left: notebook name and size, cell counts, kernel status with Restart/Abort, display preferences, pane layout |

The panel reports only what the model actually holds. A canvas notebook has a
title and **no file** — `saveNotebook` takes a path from a dialog and nothing
writes it back — so Location says the notebook is unsaved rather than inventing a
path, and Size counts characters of source rather than quoting a file size for a
file that does not exist.

Two things follow from the toolbar being *verbs*: a preference that lasts the rest
of the session belongs in the panel instead, which is where the `In[n]` label
toggle lives (`lib/properties.ts`), and the Sidebar group is **one** button —
Mathematica's equivalent group carries a chat panel as its second, and a button
that opened nothing would be worse than the asymmetry.

### Markdown text cells

A `text` cell shows **rendered Markdown** when it is not being edited and its raw
source while it is — the two states Jupyter has, and what the cell-type picker has
described as "Prose / markdown" since it was written, back when nothing rendered
any. Rendering is `lib/prose.ts` over `marked`, with `breaks: true` so a single
newline is a line break: a cell is typed like prose, and needing two trailing
spaces to end a line would read as the cell ignoring Return.

One element switches between the two states rather than two elements swapping, so
the existing handlers, handle registration and arrow-key navigation are untouched
— only what is painted into it and whether it is `contenteditable` change. An
empty cell starts in edit mode, since a rendered empty cell is a zero-height
click target.

The toolbar's **Text** group wraps the selection in `**`, `*`, `` ` `` or a link,
via `document.execCommand('insertText')` — deprecated, and still the only API that
edits a contenteditable while keeping the browser's native undo stack, and it
fires `input` so the cell's existing handler saves the new source with no extra
plumbing. The group appears for `text` only: a section or subsection is an
`<h1>`/`<h2>` that is not Markdown-rendered, so `**bold**` there would display its
own asterisks.

There is no bullet-list button and no caret-at-click-point, both for the same
reason — each needs to map a position through a contenteditable whose line boxes
may be text nodes, `<div>`s or `<br>`s, and a bullet landing mid-word or a caret
landing at the wrong offset is worse than the button not being there. Rendered
prose reaches the DOM through `{@html}`: a notebook is an executable document
whose code cells already evaluate arbitrary Mathilda, so HTML in its prose is not
a new capability, and a regex pass would look like sanitisation without being it.

`npm run check:prose` covers the pure half — the renderer and the marker
constants, read out of `prose.ts` so the checks cannot drift from it. The DOM half
needs a real pointer, since synthetic events do not reach the WKWebView.

Both panel stores are plain writables with no persistence, matching `darkMode` in
`theme.ts`. Nothing in the app persists UI state yet, and making one preference
the only setting that survives a restart would be a surprise rather than a
feature.
