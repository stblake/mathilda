---
source: src/linalg/hilbertmat.c
---
**Algorithm.** `builtin_hilbertmatrix` constructs the `m×n` Hilbert matrix with entry `(i,j) = 1/(i+j-1)`. The dimension spec (`hm_parse_dims`) is a positive integer `n` (square) or a pair `{m, n}` of positive integers; bad specs emit `HilbertMatrix::dims`, zero arguments emit `HilbertMatrix::argx`. The only recognised option is `WorkingPrecision` (`hm_parse_working_precision`, last-valid-setting-wins): `Infinity` (default) yields exact `Rational` entries via `make_rational`; `MachinePrecision` (or a digit count at/below machine precision) yields machine-precision `Real`s; a larger digit count yields MPFR reals (`mpfr_div_ui`) when built with `USE_MPFR`, degrading to machine reals otherwise (`HilbertMatrix::wprec`). Any trailing non-option argument triggers `HilbertMatrix::nonopt`.

**Data structures.** A `List` of `List`s built row by row; each entry is created by `hm_entry` according to the selected `hm_prec_mode` (`EXACT`/`MACHINE`/`MPFR`). Complexity is `O(mn)` entry constructions.
