# LucasL

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LucasL[n]
    gives the nth Lucas number L_n.
LucasL[n, x]
    gives the nth Lucas polynomial L_n(x).
Exact integer orders are computed via GMP fast doubling (numbers, using
L_m = 2 F_{m+1} - F_m) or the recurrence L_k = x L_{k-1} + L_{k-2} with
L_0 = 2, L_1 = x (polynomials); negative orders use L_{-n} = (-1)^n L_n.
Inexact or complex orders evaluate the generalized closed form
L_n = phi^n + Cos[Pi n] phi^-n (phi = GoldenRatio) numerically.
Listable; symbolic orders stay unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_lucasl` mirrors `Fibonacci`. For exact integer order it fast-doubles the Fibonacci pair `(F_m, F_{m+1})` in GMP and derives `L_m = 2 F_{m+1} - F_m` in `O(log n)`, with negative orders via `L_{-m} = (-1)^m L_m`. `LucasL[n, x]` (Lucas polynomial) iterates `L_k = x L_{k-1} + L_{k-2}` from `L_0 = 2, L_1 = x`, Expand-ing each step. For inexact/non-integer order it builds the closed form `phi^n + Cos[Pi n] phi^-n` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2`) and hands it to `numericalize`. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` integer pair; `Expr` trees through `eval_and_free` for the polynomial and closed-form branches.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/lucas.c`](https://github.com/stblake/mathilda/blob/main/src/lucas.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
