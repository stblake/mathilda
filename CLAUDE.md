## Workflow Orchestration

### 1. Plan Node Default
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- Plan mode is the up-front gate: the plan check-in *is* the approval step. Once
  the plan is approved, execute it continuously without pausing for feedback (see
  Project Specifics) — re-plan only if something goes sideways
- If something goes sideways, STOP and re-plan immediately - don't keep pushing
- Use plan mode for verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Subagent Strategy
- Use subagents liberally to keep main context window clean
- Offload research, exploration, and parallel analysis to subagents
- For complex problems, throw more compute at it via subagents
- One tack per subagent for focused execution

### 3. Self-Improvement Loop
- After ANY correction from the user: update `tasks/lessons.md` with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project
- Two lesson stores coexist: `tasks/lessons.md` (this repo's long-form,
  investigation-specific log) and the harness file-based memory (`MEMORY.md`
  index + `memory/`, auto-loaded each session for durable, generalizable
  patterns). Put durable one-fact patterns in memory; keep verbose task notes in
  `tasks/`

### 4. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between main and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness

### 5. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution"
- Skip this for simple, obvious fixes - don't over-engineer
- Challenge your own work before presenting it

### 6. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding
- Point at logs, errors, failing tests - then resolve them
- Zero context switching required from the user
- Go fix failing CI tests without being told how

## Task Management

1. **Plan First**: Write plan to `tasks/todo.md` with checkable items
2. **Verify Plan**: Check in before starting implementation
3. **Track Progress**: Mark items complete as you go
4. **Explain Changes**: High-level summary at each step
5. **Document Results**: Add review section to `tasks/todo.md`
6. **Capture Lessons**: Update `tasks/lessons.md` after corrections

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Changes should only touch what's necessary. Avoid introducing bugs.

## Project Specifics

The goal of this project is to use Claude to create a highly capable and efficient computer algebra system (CAS), called Mathilda. The CAS should be a faithful recreation of the core architecture (parser, pattern matcher, symbol table, evaluator) of Mathematica (or the Wolfram Language) and a recreation of the core simple mathematical functions of Mathematica (Plus, Times, Power, Divide, etc.)

-- When working on a task, do not pause to ask for feedback. Proceed through the plan continuously. If tests pass, move to the next phase immediately.

-- No code changes should be made to any libraries in @src/external/

-- Before any coding takes place the document @SPEC.md should be read to get an understanding of the system. 

-- The code should be well documented, with performance and scalability in mind.

-- Every time a builtin function is implemented, we should add it to the symbol table so it's accessible in the repl. 

-- Internal symbols should be defined in `src/sym_names.c` (the `SYM_*` interned-name pointers are declared in `src/sym_names.h` and defined in `src/sym_names.c`). 

-- Every time a builtin function is implemented or modified, we should update `Mathilda_spec.md`. The spec is now an overview file that points into `docs/spec/`; edit the relevant per-category file under `docs/spec/builtins/` (and the matching weekly `docs/spec/changelog/<YYYY-MM-DD>.md` — where `<YYYY-MM-DD>` is the Monday of the current ISO week — for a change summary), then update the overview only if a new top-level section or category was added.

-- Every time a builtin function is implemented we should also assign the appropriate Attributes to that function. 

-- Every new builtin that operates on numbers MUST, wherever the operation is numeric, support the packed/NDArray fast paths AND bytecode compilation. Treat these as correctness surfaces, not "optimizations to add later": a visible `NDArray[...]` argument left unevaluated, or a packed buffer silently materialised/truncated, is a WRONG answer (or a 10–40× regression that spreads to every consumer), not merely a slow one. The three surfaces, each with a concrete wiring point and an enforcing audit:
    - **NDArray[] and packed arrays.** A head that maps element-wise over numbers (or reads a numeric buffer / reduces one) needs an ND kernel — `src/ndkernels.c` for elementary/special functions, `src/ndinteger.c` for exact-integer kernels, `src/ndreduce.c` / `src/ndstruct.c` for reductions/structural ops — and an entry on the `AWARE` list (and `INT64_OK`, when the int64 answer is exact) in `src/pack.c`. Registering a kernel (`REG_U`/`REG_B`, via `symtab_set_ndarray_*_kernel`) also sets `packed_aware`. Without this, a visible `NDArray` is left unevaluated and a packed argument is un-packed into boxed `Expr`s before the head runs. Mirror an existing sibling: `UnitStep`/`Sign`/`Floor` for narrowing integer-result heads, `Sin`/`Exp` for real kernels.
    - **Compile[] and auto-compilation.** The same head must lower in `Compile[]` at both scalar and rank-1 array shapes — type inference in `src/compile/compile_infer.c`, codegen in the matching `src/compile/compile_emit_*.c`. Auto-compilation (`Table`/`Sum`/`NIntegrate`/… via `src/compile/autocompile.c`) reuses that lowering, so it follows for free — there is no separate registry. The compilable subset is a CLIFF: one unlowered head sends the WHOLE body to the interpreter, so a numeric head fast at the REPL but unlowerable in `Compile[]` is a contradiction. Verify with `CompileDiagnostics[{{v, _Real, 1}}, H[v]]` (array) and `CompileDiagnostics[{{x, _Real}}, H[x]]` (scalar) — both must report `Compiled -> True`.
    - **Enforcement & escape hatch.** `make check-packed-aware`, `check-nd-surfaces`, `check-array-exactness`, `check-fastpath-sweep`, and `check-compile-coverage` gate these (see §9 of `SPEC.md` and `docs/design/packed_arrays.md`). "Wherever possible" is real: a purely symbolic/structural head, or one returning a non-machine object (e.g. `InterpolatingFunction`), genuinely cannot support a surface — but that is then a DELIBERATE, DOCUMENTED decision, recorded in the audit tool's `EXEMPT`/`BASELINE`/`OFF_BUFFER` list with a one-line reason, never a silent omission.

-- Efficient and careful memory management is important. The system should track memory usage and leaks with valgrind. 

After any change or improvement to the system is made, a summary of the features should be given under an appropriate heading in the current week's `docs/spec/changelog/<YYYY-MM-DD>.md` file, where `<YYYY-MM-DD>` is the Monday of the current ISO week (Mon – Sun). Create the file with a short `# Changelog: week of <Mon> (Mon) – <Sun> (Sun)` heading if it does not yet exist. The top-level `Mathilda_spec.md` itself stays terse — it is the navigational entry point only; add a row to its changelog table when you create a new weekly file.

-- Every builtin function should have an Information string that gives a concise, but complete description of the function (via symtab_set_docstring)

-- The book in `book/` carries a back-of-book **Index**. Regenerate it after any major book update (a new chapter/section, substantial new prose, or a large example rework). Two mechanisms: **builtin names index themselves** — the `\B{...}` macro (`book/mathilda.sty`) auto-emits an `\index{}` entry in code font (`\mcode`), and `\usagebox` marks the definition page in bold, so referencing a builtin is all that is needed; **concept terms are indexed by hand** with `\index{topic!subtopic}` at their defining mention (follow the `Mathilda!goals` / `interval arithmetic!outward rounding` style). So "regenerate the Index" means: add/refresh `\index{}` concept entries for the new material, then rebuild (`cd book && make pdf`, which reruns `makeindex`). The index convention is documented in `book/CONTEXT.md`.

-- Code must compile cleanly under strict C99 on Linux (`gcc -std=c99 -Wall -Wextra`). glibc hides the entire POSIX surface under `-std=c99` while macOS exposes it implicitly, so an unguarded use compiles clean here and breaks only for Linux users (this has shipped as issues #36, #37, #40). `SPEC.md` §10 holds the full rationale and paste-ready guard blocks; the rules, terse:
    - **`<math.h>` constants** (`M_PI`, `M_E`, `M_PI_2`, `M_LN2`, …) are POSIX, not C99: add an `#ifndef` fallback right after `#include <math.h>` (see `src/trig.c`, `src/numeric.c`). Never `#define _USE_MATH_DEFINES` — that is MSVC/Windows-only and does nothing on glibc.
    - **POSIX functions** (`jn`, `yn`, `strdup`, `strndup`, `asprintf`, `getline`, `fileno`, `popen`, `clock_gettime`, anything from `<unistd.h>`) need a feature-test macro (`_XOPEN_SOURCE` / `_POSIX_C_SOURCE 200809L` / `_GNU_SOURCE`) placed **before any `#include`** (see `src/ndkernels.c`, `src/core.c`, `src/repl.c`, `src/loadmodule.c`). Below the first `#include` it does nothing — the header is already parsed with the wrong namespace.
    - **`int64_t` is not `long long`** — same width, distinct types (`long long` on macOS, `long` on glibc). `src/checked_int.h` ships two non-interchangeable families: `ci_mul`/`ci_powi`/… take `long long` (`src/compile/` only), `ci_mul_i64`/`ci_powi_i64`/… take `int64_t` (everywhere else). Mixing them compiles clean on macOS and is a GCC 14 error on Linux.
    - Avoid `ssize_t` (loop `n` down to `1`, index `i - 1`) and GNU/BSD extensions (nested functions, statement expressions, unguarded `__attribute__`).
    - `make check-c99` (`tools/check_c99_portability.py`) catches all of these and prints the exact fix; the Linux CI job in `.github/workflows/build.yml` compiles the whole tree against glibc on every push/PR as backstop. A new POSIX symbol the checker doesn't know is one line in its `FUNCTIONS` table.

<!-- code-review-graph MCP tools -->
## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

The short names below (`detect_changes`, `query_graph`, …) are shorthand;
the actual invocable tools carry the server prefix and a `_tool` suffix —
e.g. `mcp__code-review-graph__detect_changes_tool`,
`mcp__code-review-graph__query_graph_tool`.

### When to use graph tools FIRST

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with callers_of/callees_of/imports_of/tests_for
- **Architecture questions**: `get_architecture_overview` + `list_communities`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.

### Key Tools

| Tool | Use when |
|------|----------|
| `detect_changes` | Reviewing code changes — gives risk-scored analysis |
| `get_review_context` | Need source snippets for review — token-efficient |
| `get_impact_radius` | Understanding blast radius of a change |
| `get_affected_flows` | Finding which execution paths are impacted |
| `query_graph` | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes` | Finding functions/classes by name or keyword |
| `get_architecture_overview` | Understanding high-level codebase structure |
| `refactor_tool` | Planning renames, finding dead code |

### Workflow

1. The graph is refreshed by a `SessionStart` hook that asks Claude to run
   `build_or_update_graph_tool` (an incremental diff since `HEAD~1`) — it is
   **not** automatic on every edit, so rebuild after a batch of changes or
   after pulling.
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.
