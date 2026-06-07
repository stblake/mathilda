---
source: src/solve.c
---
`Quartics` is not a builtin function — it is an option symbol for `Solve`, registered (with its docstring) in `solve_init`. It is recognised by `is_known_option_name` and consumed in `apply_option`, which sets `opts.poly.quartics_radical = (rhs === True)`. With the default `Quartics -> False`, quartic equations are returned as held `Root[]` objects; `Quartics -> True` requests explicit radical (Ferrari) formulas from the polynomial specialist (`src/poly/solvepoly.c`). The same option is forwarded by `Eigenvalues`/`Eigensystem` (`src/linalg/`) so the characteristic-polynomial roots can likewise be returned as radicals or held `Root[]`s.
