---
source: src/core.c
---
**Algorithm.** `builtin_machineintegerq` (`src/core.c`) is a one-line type test: it takes a
single argument and returns `True` exactly when that argument's tagged-union type is
`EXPR_INTEGER` — Mathilda's 64-bit machine-word integer — and `False` otherwise. It does no
arithmetic and never looks at the value's magnitude, because the representation already
encodes the distinction: a value too large for a machine word has been promoted to
`EXPR_BIGINT` (a GMP `mpz_t`) by the arithmetic heads, so it is no longer `EXPR_INTEGER`.

**Data structures.** None beyond the input `Expr`. The function inspects `res->data.function
.args[0]->type` and constructs a fresh `True`/`False` symbol. It is registered `Protected`
in `core_init`; a non-unary call returns `NULL` and is left unevaluated.

**Complexity / limits.** `O(1)`. This is the sharp difference from `IntegerQ`:
`MachineIntegerQ[2^100]` is `False` (the value is a `BigInt`), whereas `IntegerQ[2^100]` is
`True`. It is likewise `False` for a real such as `3.0`, and `False` for `2^63`, which
overflows a signed 64-bit word and so is stored as a `BigInt`. Use it to decide whether a
value will travel on the packed/fixnum fast paths rather than the arbitrary-precision path.
