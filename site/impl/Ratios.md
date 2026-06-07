---
source: src/list.c
---
`builtin_ratios` returns successive ratios of list elements: `Ratios[{a, b, c, ...}]` gives `{b/a, c/b, ...}`. `Ratios[list, n]` applies the ratio operation n times (`ratio_n`); `Ratios[list, {n1, n2, ...}]` takes ratios along the given levels (`ratio_levels`). Each ratio is formed as a division, so the usual numeric/symbolic simplification follows. Level/count specs must be non-negative integers; bad shapes return `NULL`.
