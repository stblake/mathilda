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
