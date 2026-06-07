---
source: src/context.c
---
`builtin_begin` (`src/context.c`) takes a single string `"ctx`"` and calls `context_begin`, which pushes a `FRAME_BEGIN` onto the internal `g_stack` (snapshotting the current context and search path), then sets `g_current` to the new context. A leading-backtick argument like `"`Private`"` is interpreted relative to the current context (`MyPkg``` + `Private``). `republish_state` refreshes the `$Context`/`$ContextPath` OwnValues, and the new current context string is returned. Invalid (non-backtick-terminated) specs emit `Begin::cxt`. The matching `End[]` (`builtin_end`) pops the frame to restore the saved context.
