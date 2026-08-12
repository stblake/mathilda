# Fibonacci

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Fibonacci[n]`**

gives the nth Fibonacci number F\_n.

**`Fibonacci[n, x]`**

gives the nth Fibonacci polynomial F\_n(x).

<details>
<summary>Notes</summary>

Exact integer orders are computed via GMP fast doubling (numbers) or the recurrence F\_k = x F\_{k-1} + F\_{k-2} (polynomials); negative orders use F\_{-n} = (-1)^(n+1) F\_n. Inexact or complex orders evaluate the generalized closed form numerically. Listable; symbolic orders stay unevaluated.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Table[Fibonacci[n], {n, 10}]
Out[1]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55}

In[2]:= Fibonacci[7, x]
Out[2]= 1 + 6 x^2 + 5 x^4 + x^6

In[3]:= Fibonacci[5.8, 3]
Out[3]= 283.483

In[4]:= N[Fibonacci[15/17], 50]
Out[4]= 0.956519913924311225085822634276922986486069690120617

In[5]:= Fibonacci[1/2, 0]
Out[5]= 1/2

In[6]:= Fibonacci[1/2, 3.2]
Out[6]= 0.494833
```

### Applications (5)

```mathematica
In[1]:= Fibonacci[10]
Out[1]= 55
```

```mathematica
In[1]:= Fibonacci[100]
Out[1]= 354224848179261915075
```

```mathematica
In[1]:= Fibonacci[10, x]
Out[1]= 5 x + 20 x^3 + 21 x^5 + 8 x^7 + x^9
```

```mathematica
In[1]:= Fibonacci[200]/Fibonacci[199] // N
Out[1]= 1.61803
```

```mathematica
In[1]:= Sum[Fibonacci[k], {k, 1, 10}] == Fibonacci[12] - 1
Out[1]= True
```

## Algorithm

Mathilda -- Fibonacci numbers and Fibonacci polynomials.

Surface forms -------------

```text
  Fibonacci[n]      F_n, the nth Fibonacci number.
  Fibonacci[n, x]   F_n(x), the nth Fibonacci polynomial.
```

Evaluation strategy ------------------- The builtin follows the same two-tier philosophy as the rest of the system: exact arithmetic is computed directly in C, while inexact / symbolic-constant reductions are expressed as expression trees and handed back to the evaluator (and, for numeric requests, to `numericalize`).

```text
  * Exact integer order n:
      - Fibonacci[n]      : GMP fast-doubling, O(log n) big-integer math.
                            Negative orders via F_{-m} = (-1)^{m+1} F_m.
      - Fibonacci[n, x]   : the recurrence F_k = x F_{k-1} + F_{k-2},
                            evaluated at each step so the partial result
                            stays a canonical (and, for numeric x, fully
                            reduced) expression.

  * Inexact / non-integer order (Real, MPFR, or Complex with an inexact
    part) -- the generalized closed forms, built symbolically and then
    numericalized at the precision carried by the inputs:

        Fibonacci[n]    = (phi^n - Cos[Pi n] phi^-n) / Sqrt[5],
                          phi = GoldenRatio.
        Fibonacci[n, x] = (beta^n - Cos[Pi n] beta^-n) / Sqrt[x^2 + 4],
                          beta = (x + Sqrt[x^2 + 4]) / 2.

    Because GoldenRatio / Pi / Cos / Sqrt / Power already carry numeric
    paths (machine and MPFR), `numericalize` drives the whole reduction;
    this also yields complex results for complex arguments for free.

  * Fibonacci[n, x] with an exact non-integer order n (e.g. a Rational)
    and an exact numeric x evaluates the SAME closed form exactly, but
    only keeps the result when it collapses to a number -- so
    Fibonacci[1/2, 0] -> 1/2 while Fibonacci[1/2, x] (symbolic x) and a
    non-collapsing exact x stay unevaluated, matching the one-argument
    convention that exact non-integer orders stay symbolic until they
    reduce to a value.

  * Anything else (purely symbolic order, one-argument exact non-integer
    order with no N applied) -> NULL, leaving the call unevaluated.
```

Memory: the builtin honours the ownership contract -- it never frees `res`, returns a fresh Expr* on success or NULL otherwise, and clears every GMP / temporary tree it allocates.

## Implementation notes

**Algorithm.** `builtin_fibonacci` follows the system's two-tier split. For exact integer order `n`, `Fibonacci[n]` uses GMP **fast doubling** (`fib_mpz`: `F(2k)=F(k)(2F(k+1)-F(k))`, `F(2k+1)=F(k+1)^2+F(k)^2`) in `O(log n)` big-integer operations, with negative orders via `F(-m)=(-1)^{m+1}F(m)`. `Fibonacci[n, x]` (Fibonacci polynomial) iterates the recurrence `F_k = x F_{k-1} + F_{k-2}`, evaluating each step so the partial result stays canonical. For inexact or non-integer order, it builds the generalized closed form `(phi^n - Cos[Pi n] phi^-n)/Sqrt[5]` (`phi = GoldenRatio`, or `beta = (x+Sqrt[x^2+4])/2` for the polynomial) as an expression tree and lets `numericalize` drive the reduction at the inputs' precision — which also yields complex results for complex order for free. Purely symbolic order returns `NULL`.

**Data structures.** GMP `mpz_t` for the integer fast-doubling pair; `Expr` trees (built with `mk_fn1`/`mk_fn2` helpers, reduced via `eval_and_free`) for polynomial and closed-form paths.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[Complex](../../arithmetic/Complex/), [N](../../arithmetic/N/)

## References

- Source: [`src/fibonacci.c`](https://github.com/stblake/mathilda/blob/main/src/fibonacci.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_fibonacci.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fibonacci.c)

## Notes & additional examples

### Notes

`Fibonacci[n]` uses GMP fast-doubling, so even `Fibonacci[100]` (21 digits) is returned exactly and instantly. The two-argument form gives the Fibonacci polynomial `F_n(x)` from the recurrence `F_k = x F_{k-1} + F_{k-2}`; `F_10(x)` factors the integer Fibonacci numbers (`F_10(1) = 55`). The ratio of consecutive terms converges to the golden ratio `φ = 1.61803...`, and the telescoping identity `∑_{k=1}^{n} F_k = F_{n+2} - 1` holds exactly.
