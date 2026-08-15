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

Both panel stores are plain writables with no persistence, matching `darkMode` in
`theme.ts`. Nothing in the app persists UI state yet, and making one preference
the only setting that survives a restart would be a surprise rather than a
feature.
