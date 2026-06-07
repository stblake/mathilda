---
source: src/core.c
---
`builtin_leafcount` (`src/core.c`) returns `leaf_count_internal`, which counts 1 per non-`EXPR_FUNCTION` (atomic) node and recurses into function arguments. By default heads are counted too; the option `Heads -> False` suppresses head counting.
