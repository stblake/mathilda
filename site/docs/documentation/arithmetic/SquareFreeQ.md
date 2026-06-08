# SquareFreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SquareFreeQ[expr]
    gives True if expr is a square-free polynomial or number, and False otherwise.
SquareFreeQ[expr, vars]
    gives True if expr is square-free with respect to the variables vars.
Option GaussianIntegers -> True | False | Automatic switches to Gaussian integers.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SquareFreeQ[10]
Out[1]= True

In[2]:= SquareFreeQ[4]
Out[2]= False

In[3]:= SquareFreeQ[20]
Out[3]= False

In[4]:= SquareFreeQ[3 + 2 I]
Out[4]= True

In[5]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[5]= False

In[6]:= SquareFreeQ[2/3]
Out[6]= True

In[7]:= SquareFreeQ[6 + 6 x + x^2]
Out[7]= True

In[8]:= SquareFreeQ[x^3 - x^2 y]
Out[8]= False
```

## Implementation notes

**Algorithm.** `builtin_squarefreeq` is a `*Q` predicate (always `True`/`False` on a valid call; wrong arg count emits `SquareFreeQ::argb`, malformed options `SquareFreeQ::nonopt`, both returning NULL). `sqfree_dispatch` routes by argument kind: an integer `n` is factored with `FactorInteger` and is square-free iff every prime exponent ≤ 1 (with `0 -> False`, `±1 -> True`); a rational `p/q` iff both numerator and denominator are; a Gaussian integer `a + b I` (with `GaussianIntegers -> True`, or auto-detected from a `Complex[Integer, Integer]`) by factoring the norm `a^2 + b^2` over Z and dispatching on each rational prime's residue mod 4 (`sqfree_gaussian`); a polynomial in `vars` by, for each variable `x_i` of positive degree, computing `PolynomialGCD(p, ∂p/∂x_i)` and requiring it to be degree 0 in `x_i` (independent of the variable). Everything else (Real, symbolic) is `False`.

**Data structures.** `Expr*`; integer/Gaussian work on GMP `mpz_t` via `FactorInteger`; the polynomial path uses the expand + `PolynomialGCD` machinery. `Modulus` is parsed but only `Modulus -> 0` is honoured (non-zero emits `SquareFreeQ::modnotimpl` and leaves the call unevaluated).

- `Protected`. Not `Listable` -- passing a list of inputs treats the list as

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/poly/squarefreeq.c`](https://github.com/stblake/mathilda/blob/main/src/poly/squarefreeq.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
