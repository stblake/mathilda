---
source: src/funcprog.c
---
`builtin_freeq` (`src/funcprog.c`) returns `True` iff no subexpression of the first argument matches the pattern `form`. `freeq_at_level` walks the tree depth-first, calling the pattern matcher `match` at each level permitted by the level spec (default `{0, Infinity}`), and short-circuits to `False` on the first match. `Rational`/`Complex` are treated as atomic; `Heads -> False` excludes function heads from the search.
