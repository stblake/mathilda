# LucasL

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LucasL[n]`**

gives the nth Lucas number L\_n.

**`LucasL[n, x]`**

gives the nth Lucas polynomial L\_n(x).

<details>
<summary>Notes</summary>

Exact integer orders are computed via GMP fast doubling (numbers, using L\_m = 2 F\_{m+1} - F\_m) or the recurrence L\_k = x L\_{k-1} + L\_{k-2} with L\_0 = 2, L\_1 = x (polynomials); negative orders use L\_{-n} = (-1)^n L\_n. Inexact or complex orders evaluate the generalized closed form L\_n = phi^n + Cos\[Pi n\] phi^-n (phi = GoldenRatio) numerically. Listable; symbolic orders stay unevaluated.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Table[LucasL[n], {n, 10}]
Out[1]= {1, 3, 4, 7, 11, 18, 29, 47, 76, 123}

In[2]:= LucasL[7, x]
Out[2]= 7 x + 14 x^3 + 7 x^5 + x^7

In[3]:= LucasL[-11.]
Out[3]= -199.0

In[4]:= N[LucasL[11/3], 50]
Out[4]= 5.92396265296195541013569786219401262875198554223617
```

### Applications (5)

```mathematica
In[5]:= LucasL[10]
Out[5]= 123

In[6]:= LucasL[5, x]
Out[6]= 5 x + 5 x^3 + x^5

In[7]:= LucasL[100]
Out[7]= 792070839848372253127

In[8]:= LucasL[-7]
Out[8]= -29

In[9]:= N[LucasL[20]/LucasL[19], 20]
Out[9]= 1.6180339887498948482
```

## Algorithm

Mathilda -- Lucas numbers and Lucas polynomials.

Surface forms -------------

```text
  LucasL[n]      L_n, the nth Lucas number.
  LucasL[n, x]   L_n(x), the nth Lucas polynomial.
```

Evaluation strategy ------------------- The builtin follows the same two-tier philosophy as Fibonacci (and the rest of the system): exact arithmetic is computed directly in C, while inexact / symbolic-constant reductions are expressed as expression trees and handed back to the evaluator (and, for numeric requests, to

```text
`numericalize`).

  * Exact integer order n:
      - LucasL[n]      : GMP fast-doubling of the Fibonacci pair
                         (F_m, F_{m+1}), then L_m = 2 F_{m+1} - F_m,
                         O(log n) big-integer math. Negative orders via
                         L_{-m} = (-1)^m L_m.
      - LucasL[n, x]   : the recurrence L_k = x L_{k-1} + L_{k-2} with
                         L_0 = 2, L_1 = x, evaluated (Expand-ed) at each
                         step so the partial result stays a canonical
                         (and, for numeric x, fully reduced) expression.

  * Inexact / non-integer order (Real, MPFR, or Complex with an inexact
    part) -- the generalized closed forms, built symbolically and then
    numericalized at the precision carried by the inputs:

        LucasL[n]    = phi^n + Cos[Pi n] phi^-n,  phi = GoldenRatio.
        LucasL[n, x] = beta^n + Cos[Pi n] beta^-n,
                       beta = (x + Sqrt[x^2 + 4]) / 2.

    Because GoldenRatio / Pi / Cos / Sqrt / Power already carry numeric
    paths (machine and MPFR), `numericalize` drives the whole reduction;
    this also yields complex results for complex arguments for free.

  * Anything else (purely symbolic order, exact non-integer order with no
    N applied) -> NULL, leaving the call unevaluated.
```

Memory: the builtin honours the ownership contract -- it never frees `res`, returns a fresh Expr* on success or NULL otherwise, and clears every GMP / temporary tree it allocates.

## Implementation notes

**Algorithm.** `builtin_lucasl` mirrors `Fibonacci`. For exact integer order it fast-doubles the Fibonacci pair `(F_m, F_{m+1})` in GMP and derives `L_m = 2 F_{m+1} - F_m` in `O(log n)`, with negative orders via `L_{-m} = (-1)^m L_m`. `LucasL[n, x]` (Lucas polynomial) iterates `L_k = x L_{k-1} + L_{k-2}` from `L_0 = 2, L_1 = x`, Expand-ing each step. For inexact/non-integer order it builds the closed form `phi^n + Cos[Pi n] phi^-n` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2`) and hands it to `numericalize`. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` integer pair; `Expr` trees through `eval_and_free` for the polynomial and closed-form branches.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Complex](../../arithmetic/Complex/), [N](../../arithmetic/N/)

- Source: [`src/lucas.c`](https://github.com/stblake/mathilda/blob/main/src/lucas.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_lucas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lucas.c)

## Notes & additional examples

### Notes

`LucasL[n]` is the `n`th Lucas number `L_n` (`L_0 = 2`, `L_1 = 1`, `L_k = L_{k-1} + L_{k-2}`); `LucasL[n, x]` is the Lucas polynomial `L_n(x)`. Integer orders use GMP fast doubling for arbitrary size (so `LucasL[100]` is exact), and negative orders follow `L_{-n} = (-1)^n L_n`. Consecutive ratios `L_{n+1}/L_n` converge to `GoldenRatio` (last example). Inexact or complex orders use the closed form `phi^n + Cos[Pi n] phi^-n`. Listable; symbolic orders stay unevaluated.
