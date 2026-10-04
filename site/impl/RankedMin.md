---
source: src/sort.c
---
**Algorithm.** `builtin_ranked_min` selects the `n`-th smallest element of a list
— `RankedMin[list, n]`, with `RankedMin[list, -n]` the `n`-th largest, so
`RankedMin[list, 1]` is `Min[list]` and `RankedMin[list, -1]` is `Max[list]`. It
requires an integer, nonzero `n`, converts it to a 1-based ascending rank `r`
(`r = n` if positive, `m + n + 1` if negative; an out-of-range `r` or empty list
leaves the call unevaluated), and runs `ranked_select`. One scan chooses the
comparison path: if every element is a real numeric atom it compares exactly with
`expr_compare` (BigInt/Rational-safe), otherwise it builds a `double` key per
element with `ranked_numeric_key` (`±HUGE_VAL` for `±Infinity`, the value of a
real atom, else the machine-precision numericalisation of a symbolic real such as
`Pi`, `E`, `Sqrt[2]`) — and a non-real element (free symbol, non-real complex)
makes it decline, since a definite result needs every element to be a real
number. `ranked_select_idx` then quickselects the `r`-th position with a Hoare
partition and middle pivot (`O(m)` average), the original index breaking value
ties stably, and the exact element at that position is returned.

**Data structures.** A `size_t* idx` index array quickselected in place; an
optional `double* keys` for the approximate path; a `RankedCtx {Expr** elem;
const double* keys;}` read by `ranked_cmp`. Only the `r`-th slot is guaranteed
settled on return. The result is an `expr_copy` of the selected element in its
exact form.

**Complexity / limits.** `O(m)` average quickselect. An `NDArray` argument takes
the buffer order-statistic path `ndred_ranked_min` (`int64` exactly, reals via an
`O(m)` quickselect) — `RankedMin` is on `pack.c`'s `AWARE` and `INT64_OK` lists.
It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, RankedMin[v, 2]]` reports `Compiled -> True`.
`Protected`.
