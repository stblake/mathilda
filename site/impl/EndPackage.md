---
source: src/context.c
---
`builtin_end_package` (0-arg) calls `context_end_package`, which pops the top `CtxFrame` from the context stack via `frame_pop` (restoring `$Context` and the saved `$ContextPath`), then *prepends* the just-closed package context to `$ContextPath` with `path_prepend` (unless already present) so symbols defined inside the package remain visible after the package closes. State is republished. Unlike `End[]`, the builtin returns the symbol `Null`. On an empty stack it emits `EndPackage::noctx` and returns NULL. The frame stack is shared by `Begin`/`BeginPackage`/`End`/`EndPackage`; the `FrameKind` tag distinguishes package frames but pop logic is common.
