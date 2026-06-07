---
source: src/arithmetic.c
---
`builtin_subtract` is a thin two-argument rewrite: `a - b` becomes `Plus[a, Times[-1, b]]`. It does no arithmetic itself — the returned `Plus`/`Times` tree is canonicalised and folded by the evaluator's `Plus`/`Times` machinery. Non-binary calls return `NULL`.
