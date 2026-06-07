---
source: src/list.c
---
`UpTo[n]` is an inert specification object with no builtin handler — it is interpreted by the consumers that accept a count or position. List extractors (`Take`/`Drop` in `src/list.c`, `Part` ranges in `src/part.c`) detect `UpTo[n]` with an integer `n` and clamp the request to whatever is available: if at least `n` elements/positions exist all `n` are used, otherwise only those present, without raising the out-of-range error a bare `n` would. (Note: `SVD`'s `"UpTo"` target clamps to matrix rank, a separate use of the name.)
