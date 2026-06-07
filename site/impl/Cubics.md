---
source: src/solve.c
---
`Cubics` is not a function — it is an option symbol for `Solve`/`Roots` (and forwarded by the eigenvalue solver). It has no builtin handler; it is interned in `sym_names.c` and registered only with a docstring. Option parsing in `solve.c` (`is_*_option` / the option setter) reads `Cubics -> True|False` into `opts->poly.cubics_radical`, which controls whether degree-3 factors are solved by radicals (Cardano) or returned as held `Root[]` objects. The actual radical-vs-`Root` decision lives in the polynomial solver (`src/poly/solvepoly.c`); the default is `Cubics -> False`.
