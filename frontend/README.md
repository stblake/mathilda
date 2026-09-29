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
npm run build:app      # kernel sidecar + JS deps + release .app, one command
```

That is `build-app.sh`: it runs `build-sidecar.sh`, installs the JS dependencies
if they are missing, then `tauri build --bundles app` (the `.app` only — the
`.dmg` is for shipping to someone else and costs minutes more). The bundle lands
in `src-tauri/target/release/bundle/`. For every target, including `.dmg` /
`.deb` / `.msi`, run `npm run tauri build` with no flags.

One side effect worth knowing: `build-sidecar.sh` builds the kernel with
`USE_ECM=0` (the bundled app must not need a `libecm` the user has not installed),
and the makefile writes one output path — so it **replaces the repo's
`./Mathilda`** with that degraded build. The script now says so on the way past;
restore the default with `make -j` at the repo root.

### The clickable launcher

**`MathildaNotebook.command`**, at the repository root beside the makefile —
double-click it in Finder to open the notebook. It is the one-click entry point for someone who does not want a
terminal:

- if the release bundle is missing, it runs `build-app.sh` (the same path
  `npm run build:app` takes, not a second copy of it) and reports progress while
  it builds, then opens the app;
- if the bundle is there but any `.c`/`.h`/`.rs`/`.ts`/`.svelte` source is newer
  than it, it names one of them, prints the refresh command, and opens the app
  anyway — a click means "open it", but a notebook run from a stale bundle
  answers with a stale kernel, and every result then looks like a bug in code that
  has since changed;
- it launches the **release** bundle, which carries its own kernel sidecar. For
  hot reload use `npm run tauri dev` instead.

The `.command` suffix is what makes it clickable at all: Finder runs a `.command`
file in Terminal, where an extensionless executable would open in a text editor.
Finder hides the suffix unless you have asked it to show all extensions, so it
reads as `MathildaNotebook` in the folder. `build-app.sh` re-applies the
executable bit on every build, because losing that bit is exactly what turns a
double-click into "open in TextEdit".

## Opening the app

It opens on an empty, untitled notebook — the same thing File > New gives you.
The two guided canvases still exist as opt-in calls, `loadStartupContent()` (the
tour) and `loadDemoContent()` (the demo) in `src/lib/canvas.ts`; neither runs on
mount. Each builds the cards it needs rather than assuming the nine the canvas
once shipped with (`ensureStarterCards`), so seeding one notebook did not break
them.

## File Formats

**File > Open** reads either format, told apart by extension; **File > Save As**
writes either.

| Extension | Holds | Open | Save As |
|---|---|---|---|
| `.lb` | the whole canvas (every notebook, positions, sources) as JSON | replaces the canvas | the library; Cmd+S then saves back to it |
| `.mnb` | **one** notebook ("Mathilda notebook"), plain text | adds it to the canvas and focuses it | exports the current notebook (active pane, else the top card) |

Neither stores outputs: a result is what the kernel computed this session, and a
file that carried one would claim a computation the reader's kernel never did.

`.mathilda` is the extension a notebook was written under before v0.235. Open
still accepts it, so no file on disk stops working; Save As offers only `.mnb`.

### The `.mnb` stanza format

Each cell is a stanza, one blank line between them:

```
(* cell: code *)
Integrate[x^2, {x, 0, 1}]

(* cell: text *)
Some **Markdown** prose.
```

Parsing and writing live in Rust (`src-tauri/src/notebook_format.rs`, behind the
`load_notebook` / `save_notebook` commands) with unit tests. A marker must be the
whole line, so a comment that merely begins `(* cell:` stays in the cell; an
unknown type reads as `code`; text before the first marker becomes a leading code
cell; a file with no markers at all is one code cell, so an ordinary `.m` script
opens as a notebook. Side-by-side cells are written as consecutive stanzas.

The format is Git-diffable, and because the marker is a Mathilda comment, a
notebook whose cells are all code is also a script: `./Mathilda -file
notebook.mnb`. (Piping it to stdin does not work: a non-tty stdin selects the
NDJSON protocol below.)

## Keyboard Shortcuts

| Shortcut | Action |
|---|---|
| Shift+Enter | Run current cell |
| Ctrl/Cmd+Enter | Run current cell and insert new cell below |
| Ctrl/Cmd+L | Copy the input from above into the caret (Mathematica's Copy Input from Above) |

Cmd+L inserts the nearest **non-empty code cell** above, in document order, at the
insertion point — replacing the selection, not the cell, so a cell you have
started typing in keeps what you wrote. Empty cells are skipped: "the input
above" means the last thing you actually typed, and the fresh cell that
evaluating the last cell leaves behind would otherwise always be the answer. It
is also Edit > Copy Input from Above in the native menu. The lookup is
`previousInputSource` in `cellCommands.ts`, split out so it is testable without
an editor.

### Cell styles

Eight styles, offered by the toolbar's style control and by Cell > Convert to…:
`code`, `text`, and the heading ladder `title`, `subtitle`, `chapter`, `section`,
`subsection`, `subsubsection`. `CELL_STYLES` in `src/lib/notebook.ts` is the
single source of truth — the toolbar, the menu and the Rust `.mnb` writer's
`KNOWN_TYPES` are all checked against it, because `normalise_type` runs on
**serialize** as well as parse, so a style the Rust side does not know would be
written out as `code` and the heading destroyed.

The HTML tag and the visual size are deliberately decoupled: `section` stays `h1`
and `subsection` stays `h2` because every generated reference page is built from
them and their anchors are linked. The size ladder is carried by per-style CSS
(`.heading-title` … `.heading-subsubsection`) instead.

### What happens after you evaluate

- **A fresh cell below, when there was none.** Evaluating the last cell of a
  notebook leaves an empty code cell under it with the caret in it, as Mathematica
  and Jupyter both do. Only the single-cell paths do this; "Evaluate Notebook"
  does not, or it would end by appending a stray cell every time. Refused for
  anything but a non-empty code cell — the same condition that decides whether to
  evaluate at all — so holding Shift+Enter cannot fill the notebook. Rows, not
  cells: a sibling *beside* the evaluated cell is not below it.
- **`In[n]` is the kernel's line number.** `%`, `%%`, `%n`, `In[n]` and `Out[n]`
  all resolve against the kernel's `$Line`, and one kernel serves every notebook
  on the canvas — so the label cannot be a counter the front end keeps, or a cell
  would read `In[2]` while `%2` addressed a line typed in another pane. The kernel
  reports the line it used (a `line` protocol message per statement) and the cell
  is labelled with the first one it reports. The local counter survives only as
  the placeholder that appears the instant you press Shift+Enter.

## Architecture

```
Svelte UI (src/)
    +  @tauri-apps/api invoke + Channel
Tauri Rust layer (src-tauri/src/)
    +  stdio pipes (NDJSON)
Mathilda C binary (../Mathilda)
```

See docs/frontend-research.md for the full design rationale.

### Evaluating a cell

The Rust side sends each cell as `{"id": N, "expr": "...", "cell": true}`; the
protocol is documented at the top of `pipe_mode_loop`'s section in `src/repl.c`.
The `cell` flag gives the request a notebook input cell's semantics:

- **Several statements.** One per line, or separated by `;`, each evaluated in
  turn and each non-Null result shown, as in a Mathematica input cell. A statement
  ending in `;` shows no result. The whole cell is syntax-checked first, so an
  error anywhere evaluates nothing.
- **Print output and messages reach the cell.** The kernel captures each
  statement's `Print` text and its messages and sends them as `stream` and
  `message` lines before that statement's result. They render as plain text and
  as an amber warning block respectively.
- **No length limit.** Requests are read at any length (the kernel used to read
  them into a 10 KB buffer and cut a larger cell off).
- **A numbered session history.** One `$Line` per *statement*, with `In[n]` and
  `Out[n]` recorded, reported to the front end as a `line` message before that
  statement's output. This is what makes `%` work: the parser turns `%` into
  `Out[-1]` and `builtin_out` resolves it against `$Line`, so before this the
  notebook printed a literal `Out[-1]`. Per statement rather than per cell, so `%`
  means "the previous result" whether it came from the cell above or the line
  above inside the same cell; and recorded even for a result suppressed by `;` or
  equal to `Null`, because `x = 5;` then `%` is 5.

Plain requests without the flag behave exactly as before, which is what the site
generator and the audit tools depend on (`make check-pipe-protocol` pins both
modes). As a fallback for an older kernel, the Rust side also forwards any
non-JSON stdout line as `stream` text and stderr as a `message`, instead of
dropping and logging them. `kernel.rs`'s routing and the request line are unit
tested, including one round trip through the real binary when it is built.

### Output rendering

Output wraps rather than scrolling sideways: a long result reflows onto the next
line where it used to run off to the right with a horizontal scrollbar. Plain
text, errors and messages get `pre-wrap` plus `overflow-wrap: anywhere`, so even
an unbroken 400-character symbol name wraps.

Typeset math is the case that needs more than CSS. KaTeX lays an expression out as
inline-block boxes and **cannot** line-break, so a formula wider than the cell
would still overflow. `Output.svelte` therefore measures the rendered element
against its container (a `ResizeObserver` action) and, when the math genuinely
does not fit, re-renders that item as wrapped plain text — the re-parseable form
the kernel also sends. The decision is sticky and is reconsidered only when the
container's measured width actually *increases*: the swap changes the element's
height, so a symmetric rule would oscillate. Measured, never guessed from a
character count.

Tables, diagrams and code blocks keep their own `overflow-x: auto` container; the
page body never scrolls horizontally.

### Stacking order

Every surface positioned against the window takes its `z-index` from one scale of
`--z-*` tokens at the top of `app.css` (canvas, focused view, status dock,
overlays, minimap, app bar, popovers, banner, menus, lowest first). The properties
panel and the find bar were once 40 and 45 against the focused view's 50, so they
opened *behind* the notebook they act on; on the shared scale they sit above it.
Numbers inside one component's own stacking context stay local.

`npm run check:notebook` pins that order, the Cell menu's Convert to items (which
now convert the cell through the same `convertCell` the toolbar's cell-style
control uses, instead of only relabelling the toolbar), and the find bar's mark.

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

### Notebook search

`Cmd+F` in focused mode opens a find bar that searches **every cell** of the
notebook in the active pane. It is deliberately not `@codemirror/search`: that
package's `openSearchPanel` searches one editor — whichever cell holds focus — and
a bar that silently ignores the other forty cells while calling itself notebook
search is worse than no bar, because you would believe its "No matches". So
`lib/search.ts` matches over the notebook's own model via `store.allCells()`.

`Cmd+F` is free because `@codemirror/search` is *not* installed, so no editor
claims the binding. Enter and Shift+Enter walk the matches, wrapping both ways
(the first Enter lands on the match the count shows); typing only updates the
count, since jumping per keystroke would scroll the notebook out from under
someone still typing.

**The find field keeps focus.** Each jump *marks* the current match and scrolls it
into view without focusing its cell (`lib/searchHighlight.ts`): a code cell gets a
CodeMirror mark decoration held in its own state, a prose cell shows its source
and the range is painted with the CSS Custom Highlight API. It used to select the
match and focus the cell's editor, so the second Enter went to the cell and
replaced the match with a newline. Escape (or the close button) is what puts the
caret on the current match, as every editor's find does.

It still does not highlight every match at once, only the current one; the count
says how many there are.

`npm run check:search` compiles `lib/search.ts` with the project's own `tsc` and
imports the result, so its 24 checks exercise the shipped functions rather than a
paraphrase. Two properties carry it: matches are **non-overlapping** (`"aa"` in
`"aaaa"` is two matches, not three, which is what separates the loop from the
naive `from = at + 1`), and stepping wraps in **both** directions — in JavaScript
`-1 % 3` is `-1`, so the naive form sends Shift+Enter at the first match to a
negative index and the bar reads "0 of 3".

### The Insert palette

The **Insert** group offers CodeMirror snippet templates — Table, Matrix, Sum,
Integrate, Solve, Plot, Module, a definition, and so on — inserted at the caret
with `${field}` placeholders that Tab walks. Code cells only: every template is a
Mathilda expression, and offering `Table[]` for a prose cell would insert text that
never evaluates. The menu's hint is the *expansion*, so it shows what will land in
the cell rather than only what it is called.

`@codemirror/autocomplete` is now an explicit dependency. It was already in
`node_modules` transitively through the `codemirror` meta-package and importing it
on that basis works right up until the day the meta-package reorganises.

The templates are the risk, not the machinery: one with an unbalanced bracket
produces a broken cell every time the button is pressed, and reading
`{{${a}, ${b}}, {${c}, ${d}}}` is a poor way to notice. So
`npm run check:snippets` expands each template and feeds it to **the real Mathilda
binary** inside `Hold[...]`, which parses without evaluating — `Plot[]` must not
open a window and a definition must not define anything. A template whose expansion
does not parse fails the check, which is a stronger guarantee than counting
brackets: it is the language's own parser agreeing that what the button inserts is
valid. Bracket counts are checked too, because they say *which* kind is unbalanced
where the parser only says the line is wrong. If the C binary is not built the
parser half reports SKIP rather than failing.

### Interface scale

The scale rows in the properties panel drive the same store the `Cmd+=` / `Cmd+-`
/ `Cmd+0` bindings drive (root `font-size`, via `uiScale` in `lib/properties.ts`).
It was a local in `App.svelte`; lifting it to a store is what lets the panel show
the number the keyboard changes, rather than becoming a second scale that drifts
from the first. The panel's steps are exactly representable in binary (0.75, 1,
1.25, 1.5) so a button reads as selected on an exact equality, and only on one —
the keyboard still moves in 0.1, and highlighting the nearest step while sitting
between two would misreport the state. `Cmd+0` had been documented in that
handler's own comment without being implemented; the store made it one line.

Inline TeX renders through KaTeX: `$…$` inline, `$$…$$` display. The math is
pulled OUT before Markdown runs and put back after, and that order is the whole
point — handed to Markdown first, `$a_1 * b_2$` loses `_1 * b_2` to emphasis and
nothing downstream can recover it. Extraction is a hand-written scan rather than a
regex because everything that must *not* be math is context: a `$` inside a code
span or fence is a literal dollar, and `\$` is one anywhere. Two standard
heuristics (pandoc and markdown-it use the same pair) keep prose about money from
becoming equations — an opening `$` must be followed by a non-space and a closing
`$` preceded by one, so "it cost $5 and $6" stays prose. Unclosed delimiters are
left exactly as typed rather than swallowing the rest of the cell. The Text group's
fifth button writes `$…$`, which is how the feature gets discovered at all; the
alternative is knowing to type it.

`npm run check:prose` compiles `lib/prose.ts` with the project's own `tsc` and
imports the result, so its 39 checks exercise the shipped `renderProse`,
`extractMath` and marker constants rather than a paraphrase — a paraphrase can
agree with its test and disagree with the app. It keeps one direct `marked` import
for a single negative control: that without `breaks: true` there is no `<br>` at
all, so the option is load-bearing rather than decorative. The DOM half needs a
real pointer, since synthetic events do not reach the WKWebView.

Both panel stores are plain writables with no persistence, matching `darkMode` in
`theme.ts`. Nothing in the app persists UI state yet, and making one preference
the only setting that survives a restart would be a surprise rather than a
feature.
