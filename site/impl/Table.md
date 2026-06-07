---
source: src/list.c
---
**Algorithm.** `builtin_table` (in `src/list.c`) implements iterator-driven list construction. `Table` carries `HoldAll | Protected` (set in `src/attr.c`), so the body and iterator specs arrive unevaluated. Multi-iterator forms `Table[expr, spec1, ..., specN]` are handled by *rewriting*: the handler peels off the outermost iterator, wraps the remaining iterators in a fresh `Table[Table[expr, ..., specN], spec1]`, and recurses through the evaluator — so nesting depth equals iterator count and the rightmost spec varies fastest.

For a single iterator the spec is parsed by the shared `iter_spec_parse` (`src/iter.c`) into an `IterSpec` discriminated by `kind`: `ITER_KIND_COUNT` (`{n}` or bare `n` — repeat the body), `ITER_KIND_LIST` (`{i, {v1, v2, ...}}` — iterate over explicit values), or `ITER_KIND_RANGE` (`{i, imax}`, `{i, imin, imax}`, `{i, imin, imax, di}`). Range bounds are resolved to doubles via `iter_spec_resolve_numeric` (with `allow_inf=false`, since `Table` never iterates to `Infinity`).

**Localization.** The iteration variable is dynamically scoped, not renamed. `iter_spec_shadow(var)` saves and clears the symbol's `own_values` chain; each step calls `symtab_add_own_value` to bind the current index, `evaluate`s the body, then the next index is computed symbolically (`evaluate(Plus[curr, di])`) so exact integers/rationals are preserved while a parallel `double` accumulator drives loop termination (with a `1e-14` tolerance and a 1,000,000-step safety cap). `iter_spec_restore` reinstalls the saved `own_values` afterward, so the variable's outer value is untouched.

**Data structures.** Results accumulate in a geometrically grown `Expr**` buffer, finally wrapped as `List[...]`.
