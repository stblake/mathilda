---
source: src/assoc.c
---
**Algorithm.** `builtin_keyvaluemap` builds `List[f[k1, v1], f[k2, v2], …]` — one
`f` applied to the key and value of each entry, in order. The applications are
*not* forced inside the builtin; the enclosing evaluator reduces each `f[ki, vi]` as
usual, so an operator like `#1 -> #2^2 &` produces ordinary rules.

**Data structures.** One `List` of freshly constructed `f[k, v]` nodes; no hash index
(it is a straight ordered walk).

**Complexity / limits.** O(n) to construct, plus whatever the evaluator spends
reducing each application. Unlike `KeyMap` (transforms keys, result is an
association) and `Map` over an association (transforms values), `KeyValueMap` sees
both parts and returns a plain `List`.
