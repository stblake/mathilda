
The goal of this project is to use Gemini cli to create a small (pico) computer algebra system (CAS), called picocas. The CAS should be a faithful recreation of the core architecture (parser, pattern matcher, symbol table, evaluator) of Mathematica (or the Wolfram Language) and a recreation of the core simple mathematical functions of Mathematica (Plus, Times, Power, Divide, etc.)

-- No code changes should be made to any libraries in @src/external/

-- Before any coding takes place the document @SPEC.md should be read to get an understanding of the system. 

-- The code should be well documented, with performance and scalability in mind.

-- Every time a builtin function is implemented, we should add it to the symbol table so its accessible in the repl. 

-- Internal symbols should be defined in sym_names.c 

-- Every time a builtin function is implemented or modified, we should update the picocas_spec.md file. 

-- Every time a builtin function is implemented we should also assign the appropriate Attributes to that function. 

-- Efficient and careful memory management is important. The system should track memory usage and leaks with valgrind. 

After any change or improvement to the system is made, a summary of the features should be given in picocas_spec.md under an appropriate heading/section etc. 

-- Every builtin function should have an Information string that gives a concise, but complete description of the function (via symtab_set_docstring)

-- Code must compile cleanly under strict C99 on Linux (gcc -std=c99 -Wall -Wextra). Do NOT use POSIX-only or non-C99 types/functions without a portable fallback. In particular:
    - Avoid `ssize_t` (POSIX, not C99). For reverse iteration over a `size_t` count, loop from `n` down to `1` and index with `i - 1`, or use a fixed-width signed type like `int64_t` from `<stdint.h>`.
    - Avoid `strdup`, `strndup`, `asprintf`, `getline`, `fileno`, `popen`, etc. without guarded fallbacks — these are not in C99.
    - Avoid GNU/BSD extensions (e.g., nested functions, statement expressions, `__attribute__` without a portability guard).
    - Darwin (macOS) headers often expose POSIX symbols implicitly that Linux glibc hides under `-std=c99`. Test any new system-header usage on both platforms, or include the correct feature-test macros / headers explicitly.
    - `<math.h>` constants `M_PI`, `M_E`, `M_PI_2`, `M_LN2`, etc. are POSIX, NOT C99. glibc hides them under `-std=c99` (macOS exposes them anyway, masking the bug). When a new file needs one, add a guarded fallback right after `#include <math.h>`, matching the pattern already used in `src/trig.c` and `src/numeric.c`:
        ```c
        #ifndef M_PI
        #define M_PI 3.14159265358979323846
        #endif
        ```
      Do NOT use `#define _USE_MATH_DEFINES` — that is an MSVC/Windows mechanism and has no effect on glibc.

<!-- code-review-graph MCP tools -->
## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

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

1. The graph auto-updates on file changes (via hooks).
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.
