---
source: src/context.c
---
Not a builtin: `$Context` is an ordinary symbol whose OwnValue mirrors the internal context state held in `src/context.c` (`g_current`, a heap string always ending in a backtick). Whenever `Begin`/`BeginPackage`/`End`/`EndPackage` change state, `republish_state` calls `publish_own("$Context", context_as_string())`, which clears the symbol's OwnValues and reinstalls a single string-valued rule. The symbol is marked `ATTR_PROTECTED` in `context_init` so the user cannot reassign it directly — mutation only flows through the context-stack builtins. `context_resolve_name` reads `g_current` to qualify bare identifiers at parse/lookup time.
