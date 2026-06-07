---
source: src/replace.c
---
**Algorithm.** `builtin_replace_repeated` (`//.`) is a two-argument wrapper over `apply_replace_repeated_nested`. If the rule argument is a list of *lists* of rules (no top-level `Rule`), it threads over them, producing a list of independently-repeated results. Otherwise it iterates a fixed-point loop: each pass applies `apply_replace_all_nested` (the same top-down `ReplaceAll` traversal used by `/.`), then `eval_and_free`s the result, and compares against the previous expression with `expr_eq`. It returns when the expression stops changing, with a hard cap of 65536 iterations as a divergence guard. The input is deep-copied up front so the original is not mutated.
