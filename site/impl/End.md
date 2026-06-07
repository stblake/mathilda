---
source: src/context.c
---
`builtin_end` (0-arg) calls `context_end`, which pops the top frame of the context stack (`g_stack`, a linked list of `CtxFrame` pushed by `Begin`/`BeginPackage`). `frame_pop` restores `$Context` (`g_current`) and the `$ContextPath` snapshot (`saved_path`) that the matching `Begin[]` saved, frees the popped frame, and republishes the live state. The builtin returns the *closed* context string (the one that was current just before the pop). If the stack is empty it emits `End::noctx` and returns NULL (unevaluated). Context names are owned as plain `char*` via the file's `ctx_strdup` (a C99-safe `strdup` replacement).
