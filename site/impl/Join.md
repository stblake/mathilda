---
source: src/list.c
---
`builtin_join` (in `src/list.c`) concatenates its arguments via the helper `join_at_level`. A trailing integer argument is interpreted as a level specification (default 1): at level 1 the arguments' top-level elements are spliced into a single result sharing the first list's head; deeper levels splice element-wise at the corresponding depth. Returns `NULL` if no lists remain or the level is below 1.
