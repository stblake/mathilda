---
source: src/core.c
---
`builtin_bytecount` (`src/core.c`) recursively sums `byte_count_internal` over the tree: `sizeof(Expr)` per node, plus `strlen+1` for symbol/string payloads and `sizeof(Expr*) * arg_count` for each function's argument array, descending into the head and all arguments. It returns an integer; the count is a structural estimate and does not account for GMP/MPFR limb storage of bigints/reals.
