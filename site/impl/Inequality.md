---
source: src/comparisons.c
---
**Algorithm.** `builtin_inequality` handles the variadic chained-comparison head the parser emits for mixed chains like `a < b <= c == d`. The argument layout is strict: `2k+1` args, even slots values, odd slots bare operator symbols (`Less`/`LessEqual`/`Greater`/`GreaterEqual`/`Equal`); a malformed layout returns NULL. Each adjacent pair `(v_i, op_i, v_{i+1})` is resolved by `decide_pair`, which routes through `compare_numeric`/`expr_eq` and returns 1 (True), 0 (False), or -1 (undecidable). Any False pair short-circuits the whole chain to `False`. If all pairs are True → `True`. Otherwise the proven-True pairs are dropped and a residual `Inequality[...]` is rebuilt from the undecided pairs (collapsing to a single binary `op[a,b]` head when exactly one survives, and de-duplicating a shared boundary value so the residual stays a well-formed `2k+1` chain).

**Data structures.** A `bool* undecided` scratch array (one slot per pair) flags which pairs to keep; survivors are deep-copied into a fresh argument array that becomes the residual chain.
