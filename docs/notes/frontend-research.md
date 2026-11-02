# Mathilda Desktop Notebook Frontend: Research Report

*Research for the next-generation CAS notebook frontend — Rust + Tauri 2.*
*Use this document as the sole context for planning and implementation.*

---

## 0. Mathilda Background (for frontend context)

- Mathilda is a C99 Mathematica-like CAS (~22 kLoC) with a GNU Readline REPL
- Evaluation pipeline: `String → Parser → Evaluator → Printer → String`
- The REPL reads one expression per line, evaluates to a fixed point, prints `Out[n]=`, frees memory
- Output is currently plain text (e.g. `x^2 + 1`, `{1, 2, 3}`)
- Graphics subsystem: `Plot[]`/`Show[]`/`Graphics[]` via Raylib (optional)
- Build: `make USE_ECM=0` → `./Mathilda` binary
- The binary is the sidecar the frontend will spawn and own

---

## 1. Why Existing Tools Fail: Wolfram & Jupyter Pain Points

### 1.1 Hidden State / Out-of-Order Execution (Most Critical)

Both systems maintain a live kernel accumulating state independently of visual cell order. Running cells out of order leaves variables in states no visual inspection can reconstruct.

**Mathematica (documented by Wolfram):** "the kernel remembers all operations in the order they are sent to it, and this order may have nothing to do with the order in which commands are displayed in the front end." A rule `h[m,e] = 10` stored before `m` has a numeric value silently produces different results on each subsequent evaluation.

**Joel Grus (JupyterCon 2018):** the format simultaneously encourages iterative development and requires cells to run in exact order — a direct contradiction.

### 1.2 Reproducibility: The Quantitative Record

| Study | Sample | Failed execution |
|---|---|---|
| Pimentel et al. (MSR 2019) | 1.4M GitHub notebooks | ~76% (95.97% failed full reproducibility) |
| Samuel & Mietchen (GigaScience 2023) | 27,271 biomedical notebooks | 92.4% (accompanying peer-reviewed publications) |
| Chattopadhyay et al. CHI 2020 | 20 interviews + 156-person survey | 9 named structural pain points documented |

### 1.3 Version Control: The JSON Format Problem

`.ipynb` is JSON with interleaved source, base64-encoded outputs, execution counts, and kernel metadata.

- One matplotlib plot embeds as hundreds of lines of base64
- Re-running any cell (without code change) increments `execution_count` → spurious diff
- Two people opening the same notebook generate kernel metadata conflicts without touching code
- Git conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`) break JSON validity → Jupyter cannot open the file
- `mathematica-notebook-filter` README: "Wolfram Research unfortunately does not appear to offer any specification to their language or file formats" — all tooling is reverse-engineered

**nbstripout (4,000+ stars) and nbdime exist specifically to paper over these problems.** Their popularity is evidence of ecosystem-wide severity.

### 1.4 Performance (Measured)

**Jupyter:**
- 1-second lag per keystroke in larger notebooks
- 5-second sidebar expansion lag
- 10-second tab move lag
- Root causes: MathJax LaTeX rendering; CSS custom properties injected into `:root` for every cell

**Mathematica:**
- 3 seconds per character typed in large cells (29-page, 67K-character cell)
- Problem caused by the syntax-coloring renderer — persisted from v7 through v12.3 (decade without resolution)
- "Formatting Notebook Contents" startup: 40 seconds on v10.1

### 1.5 Vendor Lock-In (Mathematica)

- ~$1,500/person/year industrial licensing
- French physics institute (IPhT Saclay): went from 200 licenses at €25K/year to 10 licenses at same cost after policy change
- Cannot be freely embedded in Docker, CI/CD, or cloud deployments
- Paul Romer (Nobel laureate): PDF output "typography so bad that someone must have worked at making it bad"; sharing required a 1.3 GB CDF Player download

### 1.6 What Each Does Well (to replicate, not discard)

**Mathematica strengths:**
- First-class symbolic math with typeset output
- Dynamic interactivity (`Manipulate`, `DynamicModule`)
- Pattern matching, rule rewriting, infinite evaluation semantics
- Publication-quality math layout natively

**Jupyter strengths:**
- Open ecosystem: custom kernels, extensions
- Python library access
- VS Code integration

### 1.7 The Opportunity

Take Mathematica's evaluation semantics + symbolic output quality. Deliver it without: vendor lock-in, opaque format, Python/JS ecosystem isolation, and the state-corruption failure modes of both tools. The reproducibility bar is low — 96% of Jupyter notebooks fail it.

---

## 2. Tauri 2 + Rust Capabilities

### 2.1 What Changed from Tauri 1 → Tauri 2 (Stable: October 2024)

**IPC — Complete Rewrite**
- v1 forced serialization of all messages to strings (slow)
- v2 uses custom protocols resembling HTTP with binary payload support
- `tauri::ipc::Channel<T>`: new type for ordered, high-throughput streaming from Rust → frontend. Frontend creates `Channel<T>`, passes to Rust via `invoke()`, Rust calls `.send()` repeatedly. **This is the correct primitive for streaming CAS computation steps.**
- Raw payloads: v2 supports `Vec<u8>`, BSON, Protobuf — no JSON encoding overhead for binary data

**Security Model**
- v1 `allowlist` → v2 ACL capability files in `src-tauri/capabilities/*.json`
- Every command denied by default; each window must be granted explicit permission

**Multi-Window / Multi-Webview**
- Multiple webviews per window available behind `unstable` feature flag
- New types: `tauri::Webview`, `tauri::WebviewWindow`

**Breaking Changes to Know**
- `tauri::api::process` removed → replaced by `tauri-plugin-shell`
- Windows assets served on `http://tauri.localhost` (not `https://`) — resets IndexedDB/LocalStorage on upgrade
- Linux minimum: webkit2gtk 4.1 (Ubuntu 22.04+)
- Sidecar binaries must have target-triple suffix: e.g. `mathilda-aarch64-apple-darwin`
- `CommandEvent::Stdout` delivers `Vec<u8>`, not `String` — decode with `String::from_utf8_lossy`

### 2.2 Spawning and Communicating with the Mathilda C Process

Tauri's `tauri-plugin-shell` handles sidecar spawning:

```rust
use tauri_plugin_shell::ShellExt;
use tauri_plugin_shell::process::CommandEvent;

let (mut rx, mut child) = app.shell()
    .sidecar("mathilda")
    .unwrap()
    .spawn()
    .expect("Failed to spawn Mathilda");

tauri::async_runtime::spawn(async move {
    while let Some(event) = rx.recv().await {
        match event {
            CommandEvent::Stdout(bytes) => {
                let line = String::from_utf8_lossy(&bytes);
                // parse NDJSON response, emit to frontend via Channel
            }
            CommandEvent::Terminated(_) => {
                // restart the process with exponential backoff
                break;
            }
            _ => {}
        }
    }
});

// Write: append \n since C reads line-by-line
child.write(b"{\"id\": 1, \"expr\": \"Integrate[x^2, x]\"}\n").unwrap();
```

**Critical C-side requirement — the most common gotcha:** when stdout is piped (not a terminal), libc buffers it with a 4KB buffer. Responses sit in libc's buffer indefinitely and Tauri never receives them. Fix in Mathilda's `main()`:

```c
setvbuf(stdout, NULL, _IONBF, 0);  // unbuffered stdout — mandatory
```

Also call `fflush(stdout)` after each response as belt-and-suspenders.

For lower-level control, use `tokio::process::Command` directly with `stdin(Stdio::piped())` and `stdout(Stdio::piped())`. Spawn stdin writer and stdout reader as **separate Tokio tasks** — never write-all-then-read-all sequentially or the pipe buffer fills and deadlocks. Use `kill_on_drop(true)` and always call `wait()` after `kill()` to reap zombies.

**processkit crate**: higher-level Tokio subprocess manager with no-orphan guarantees (Linux cgroup v2, Windows Job Object, POSIX process group), built-in readiness probes, supervision/restart with backoff. Worth evaluating over raw `tokio::process`.

### 2.3 Code Editor: CodeMirror 6 vs Monaco

| | CodeMirror 6 | Monaco (VS Code) |
|---|---|---|
| Bundle size | ~300 KB | ~10 MB |
| CSP compatibility | Clean — no workers, no `blob:`, no `unsafe-inline` | Requires `'unsafe-inline'` and `blob:` relaxations |
| Tauri CSP | Works with default strict policy | Requires policy holes |
| Extensibility | Modular, composable | Monolithic |

**CodeMirror 6 is the only viable choice for Tauri.** Monaco's CSP requirements are incompatible with Tauri's security model without deliberate weakening.

### 2.4 Syntax Highlighting

Run `tree-sitter-highlight` v0.25.x server-side in Rust. A custom tree-sitter grammar for Mathilda syntax enables highlights computed in Rust and injected into the webview as HTML spans.

### 2.5 Rich Text / Prose Cells

**TipTap** (ProseMirror-based) + `@tiptap/extension-mathematics`: supports KaTeX inline (`$...$`) and block (`$$...$$`) math in prose cells — handles the "mixed prose + math" cell type Mathematica's Text cells provide.

### 2.6 Asset Serving

`wry`'s `register_uri_scheme_protocol` serves dynamic cell outputs (SVG, plot data) from Rust into the webview without filesystem writes.

---

## 3. Math Rendering

### 3.1 KaTeX vs MathJax 3

| | KaTeX v0.16 | MathJax 3.x |
|---|---|---|
| Complex page render | ~137 ms | ~137–300 ms |
| Bundle (min+gzip) | ~100–110 kB | ~59–150 kB (on-demand) |
| Rendering mode | **Synchronous** | Async by default |
| `align` environment | No (use `aligned`) | Yes |
| Equation numbering (`\eqref`) | No | Yes |
| SVG output | No | Yes |
| `katex-rs` Rust crate | **Yes** — server-side pre-rendering | No |

**Mathematica-style notation in KaTeX:** `\frac`, `\dfrac`, `\cfrac`, definite integrals with limits, sums/products with sub/superscript, all matrix bracket styles, `\begin{cases}` piecewise, `\lim_{x\to 0}` — all supported. Gap: `align` environments (use `aligned`) and equation numbering.

**Key Tauri advantage of KaTeX:** `katex-rs` Rust crate enables server-side pre-rendering. Rust generates the LaTeX string, calls `katex_rs::render()`, sends finished HTML to the frontend — zero client-side rendering latency, no flash of unstyled content.

**Recommendation:** KaTeX for the common case. Add MathJax 3 only if `align`, SVG output, deep accessibility, or equation cross-referencing become required. Bundle KaTeX CSS/fonts locally via Vite — never load from CDN for offline use.

### 3.2 Plot Rendering

- **Canvas 2D**: practical up to ~50K data points at 60 FPS
- **WebGL**: needed for dense 3D surfaces, large parametric plots, real-time manipulation
- **SVG**: appropriate for static diagrams; degrades above ~2,000 DOM elements

**Recommendation: Plotly.js** with `scattergl` type (Canvas/WebGL). Interactive 2D/3D, handles 1M+ points, has `surface` for 3D. Use `Plotly.react()` for incremental updates without full re-render. Add D3.js v7 only for custom mathematical geometry that Plotly cannot express.

---

## 4. Cell Execution Model

### 4.1 Reactive vs Linear: The Critical CAS Finding

Observable (JS reactive) and Marimo (Python reactive) both use static AST analysis to build a dependency DAG between cells and auto-re-execute dependents on change.

**This model fundamentally cannot work for a CAS without major caveats:**

- If cell A adds `f[x_] := x^2` to the symbol table and cell B calls `f[3]`, the dependency exists but is **invisible to AST analysis** — no shared variable names
- Non-deterministic algorithms (Monte Carlo, randomized primality) silently produce non-reproducible results on reactive re-runs
- Non-idempotent side effects (`Export[]`, file writes) fire unexpectedly on auto-re-run
- Reactive cascades fire on every keystroke — expensive CAS operations (`NIntegrate`, `Factor`) run constantly
- The 2025 paper "When Are Reactive Notebooks Not Reactive?" quantified a ~1.44x spurious re-run ratio

**Conclusion:** Do not make reactive mode the default. Offer it as opt-in with lazy mode (mark stale, user confirms before executing).

### 4.2 IPC Latency Benchmarks

| Mechanism | p50 RTT | p90 RTT |
|---|---|---|
| AF_UNIX domain socket | 1.4 µs | 1.7 µs |
| Anonymous pipe | 4.3 µs | 5.2 µs |
| TCP loopback | 7.3 µs | 7.9 µs |
| Shared memory | ~0.2 µs | — |

At UI timescales all three are imperceptible vs. CAS computation time. **Use stdio pipes** — zero setup, cross-platform (works on Windows without AF_UNIX complications), and `tauri-plugin-shell`'s `CommandEvent::Stdout` is already designed for them.

Do not implement the Jupyter kernel protocol (5-socket ZMQ, HMAC-SHA256 signing, multipart wire framing) — overengineered for a single-user local desktop tool, ~8ms minimum latency, substantial complexity for zero benefit.

### 4.3 Wire Protocol: NDJSON with Correlation IDs

**Request** (Tauri Rust → Mathilda stdin, one line):
```json
{"id": 1, "expr": "Integrate[x^2, {x, 0, 1}]"}
```

**Response** (Mathilda stdout → Tauri, one object per line):
```json
{"id": 1, "type": "stream", "text": "Evaluating..."}
{"id": 1, "type": "expr",   "payload": ["Rational", 1, 3]}
{"id": 1, "type": "plot",   "payload": {"data": [...], "layout": {...}}}
{"id": 1, "type": "done"}
{"id": 1, "type": "error",  "message": "...", "trace": "..."}
```

Output types:
- `"expr"` — **MathJSON** (CortexJS S-expression format: `["Add", 1, "x"]`). Mathilda's `Expr*` trees map almost 1:1. Enables: streaming, frontend symbolic manipulation (click subexpression → send subtree back to CAS), semantic preservation without re-parsing LaTeX.
- `"plot"` — Plotly JSON spec; frontend calls `Plotly.react(div, data, layout)`
- `"html"` — prose, tables, error messages with markup
- `"stream"` — partial/incremental text during long computations
- `"done"` — signals end of response for this ID

**Important:** a partial LaTeX string is syntactically invalid — cannot render `\frac{x^2` without closing `}`. Always buffer until a complete expression is delimited before calling KaTeX. MathJSON sub-expressions can be streamed token-by-token via an incremental JSON parser.

### 4.4 Tauri Streaming: Channel API

```rust
#[tauri::command]
async fn evaluate_cell(
    expr: String,
    channel: tauri::ipc::Channel<serde_json::Value>,
    mathilda: tauri::State<'_, MathildaProcess>,
) -> Result<(), String> {
    let process = mathilda.inner();
    process.send_expr(&expr).await?;
    while let Some(msg) = process.recv_line().await {
        channel.send(msg).unwrap();
        if msg["type"] == "done" { break; }
    }
    Ok(())
}
```

Frontend:
```typescript
const channel = new Channel<OutputMessage>();
channel.onmessage = (msg) => {
    if (msg.type === "expr")   renderMathJSON(cellId, msg.payload);
    else if (msg.type === "stream") appendStream(cellId, msg.text);
    else if (msg.type === "plot")   renderPlotly(cellId, msg.payload);
};
await invoke('evaluate_cell', { expr: cellSource, channel });
```

### 4.5 Recommended Execution Model

**Primary: Linear document with explicit execution**
- Cells execute in document order
- Running a cell evaluates it in the current kernel state (accumulated symbol table) — matches Mathematica's model and user expectations
- "Run All": restarts the Mathilda process, executes all cells top-to-bottom in a fresh context — guarantees reproducibility
- "Run Cell": evaluates in the live kernel

**Optional: Reactive mode (opt-in, lazy)**
- Cells declare explicit `deps` in metadata
- Do not auto-detect symbol table dependencies via AST — the mutable symbol table is invisible to static analysis
- Lazy mode: cells marked stale, user confirms before expensive CAS operations fire

**What to explicitly not implement:**
- Automatic hidden-dependency tracking through the symbol table (unsolvable without heavy runtime instrumentation)
- The Jupyter kernel protocol
- Shared memory IPC

### 4.6 Process Lifecycle

```
App startup:
  → spawn Mathilda via tauri-plugin-shell sidecar
  → store CommandChild in State<Mutex<Option<CommandChild>>>
  → spawn stdout reader task
  → send {"type": "ping"} and await {"type": "pong"} before accepting cell evals

Cell evaluation:
  → serialize to NDJSON request line → write to child stdin
  → stdout reader task dispatches responses via Channel to frontend

Process death:
  → CommandEvent::Terminated fires
  → restart with exponential backoff (max 3 attempts)
  → emit "kernel_died" / "kernel_restarting" status event to frontend

Notebook close:
  → send {"type": "quit"} to Mathilda stdin
  → wait for graceful exit (2s timeout) → SIGKILL → wait() to reap zombie
```

---

## 5. Recommended Technology Stack

| Layer | Choice | Rationale |
|---|---|---|
| Desktop framework | **Tauri 2** (stable Oct 2024) | Binary size, security model, Channel API for streaming |
| Frontend language | **TypeScript** | Type safety across IPC boundary |
| Frontend framework | **Svelte 5** or React 18 | Both integrate cleanly with Tauri invoke/Channel |
| Code cells | **CodeMirror 6** | CSP-safe, ~300KB, modular, custom Mathilda grammar |
| Math rendering | **KaTeX** + `katex-rs` (pre-render in Rust) | Synchronous, all Mathematica-style notation, offline-safe |
| Prose cells | **TipTap** + `@tiptap/extension-mathematics` | Prose + inline math in one component |
| Syntax highlighting | **tree-sitter-highlight** (Rust, server-side) | Custom Mathilda grammar, pre-computed |
| Plot rendering | **Plotly.js** (`scattergl`, `surface`) | Interactive 2D/3D, 1M+ points |
| IPC: Tauri ↔ Mathilda | **stdio pipes** via `tauri-plugin-shell` | Zero setup, cross-platform, 4µs RTT |
| IPC protocol | **NDJSON** with correlation IDs | Line-delimited, debuggable, maps to `Expr*` trees |
| Output format | **Discriminated JSON envelope** (MathJSON / Plotly / HTML) | Semantic, streaming-safe, frontend-interactable |
| Streaming to webview | **Tauri v2 `Channel` API** | Ordered, high-throughput |
| Subprocess management | **processkit** or `tokio::process` | No-orphan guarantees, supervision/restart |
| Execution model | **Linear document + explicit execution** | Predictable, matches CAS semantics |
| Notebook file format | **Plain text `.mathilda`** | Git-diffable, no embedded outputs |

---

## 6. The Single Biggest Win Over Existing Tools

**Store no output in the notebook file.**

Outputs are ephemeral, computed on demand. The notebook file is source code only — pure text, clean diffs, no JSON blob problems, valid if it parses.

This eliminates:
- The entire version-control / merge-conflict class of problems (Jupyter `.ipynb`)
- File bloat from embedded base64 graphics (Mathematica `.nb`)
- The need for nbstripout, nbdime, or mathematica-notebook-filter

The file format is `input cell source, one per stanza` — think `.m` files Mathilda already understands.

---

## 7. What Needs to Change in the Mathilda C Binary

To support the frontend, the Mathilda REPL needs a small, backward-compatible protocol mode:

1. **`setvbuf(stdout, NULL, _IONBF, 0)`** at startup when stdout is a pipe (detect via `isatty(fileno(stdout))`)
2. **NDJSON input mode**: detect `{"id": N, "expr": "..."}` as input → evaluate → emit NDJSON responses
3. **Structured output**: instead of `printf("Out[%d]= %s\n", ...)`, emit `{"id": N, "type": "expr", "payload": "..."}` — initially as a LaTeX or plain-text string; later as MathJSON
4. **Plot output**: when `Plot[]`/`Show[]` is evaluated, serialize the plot data to Plotly JSON instead of opening a Raylib window
5. **Ping/pong**: respond to `{"type": "ping"}` with `{"type": "pong"}` for readiness detection
6. **Graceful quit**: respond to `{"type": "quit"}` by flushing and exiting cleanly

The existing Readline REPL mode is preserved when `isatty(stdin)` is true — no breakage for terminal users.

---

## Sources

- [Tauri 2.0 Stable Release](https://v2.tauri.app/blog/tauri-20/)
- [Tauri v2 Sidecar Documentation](https://v2.tauri.app/develop/sidecar/)
- [Tauri v2 Channel API](https://v2.tauri.app/develop/calling-frontend/)
- [Upgrade from Tauri 1.0 Guide](https://v2.tauri.app/start/migrate/from-tauri-1/)
- [Jupyter Kernel Messaging Protocol](https://jupyter-client.readthedocs.io/en/latest/messaging.html)
- [Pimentel et al. MSR 2019 — 1.4M notebooks](https://leomurta.github.io/papers/pimentel2019a.pdf)
- [Samuel & Mietchen GigaScience 2023](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad113/7516267)
- [Chattopadhyay et al. CHI 2020](https://dl.acm.org/doi/fullHtml/10.1145/3313831.3376729)
- [Joel Grus — I Don't Like Notebooks (JupyterCon 2018)](https://www.youtube.com/watch?v=7jiPeIFXb6U)
- [Paul Romer — Jupyter vs. Mathematica](https://paulromer.net/jupyter-mathematica-and-the-future-of-the-research-paper/)
- [Marimo Dataflow Architecture](https://marimo.io/blog/dataflow)
- [Observable Runtime](https://github.com/observablehq/runtime)
- [arxiv: When Are Reactive Notebooks Not Reactive? (2025)](https://arxiv.org/html/2511.21994)
- [IPC Latency Benchmarks](https://kamalmarhubi.com/blog/2015/06/10/some-early-linux-ipc-latency-data/)
- [processkit crate](https://crates.io/crates/processkit)
- [mathematica-notebook-filter README](https://github.com/JP-Ellis/mathematica-notebook-filter)
