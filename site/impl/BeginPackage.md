---
source: src/context.c
---
`builtin_begin_package` (`src/context.c`) accepts `"ctx`"` and an optional `Needs` list of context strings, dispatching to `context_begin_package`. That pushes a `FRAME_PACKAGE` frame (snapshotting context and path), sets `g_current` to the absolute package context (relative `` `... `` contexts are rejected), clears the search path and rebuilds it as `{ctx`, System`}` plus any valid, non-duplicate `Needs` entries. `republish_state` then republishes `$Context`/`$ContextPath`, and the new context is returned. `EndPackage[]` (`context_end_package`) pops the frame and prepends the closed context to `$ContextPath`. Invalid specs emit `BeginPackage::cxt`.
