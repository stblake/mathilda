# FromContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FromContinuedFraction[{a1, a2, ..., an}]
    reconstructs a1 + 1/(a2 + 1/(a3 + ... + 1/an)). The terms may be
    symbolic; the result is the convergent in nested (un-expanded) form.
FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]
    gives the exact quadratic irrational whose continued-fraction terms
    begin with the ai then cycle through the bi forever; all ai and bi
    must be integers. FromContinuedFraction[{}] is 0. It is the inverse
    of ContinuedFraction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FromContinuedFraction[{2, 1, 3, 4}]
Out[1]= 47/17

In[2]:= FromContinuedFraction[{a, b, c, d}]
Out[2]= (1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)

In[3]:= FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]
Out[3]= Sqrt[71]

In[4]:= FromContinuedFraction[{{1, 2, 3, 4}}]
Out[4]= 1/15 (9 + 2 Sqrt[39])

In[5]:= FromContinuedFraction[ContinuedFraction[Pi, 3]]
Out[5]= 333/106
```

## Implementation notes

**Algorithm.** `builtin_from_continued_fraction` reconstructs the value from a `List` of terms. A non-periodic list uses the standard convergent recurrence `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}` (`fcf_simple`), evaluating each step so numeric terms collapse and symbolic terms stay in convergent form; the result is `h_{n-1}/k_{n-1}`. A trailing sub-list marks the cyclic (period) block, requiring all integer terms: `fcf_periodic` builds the period's convergents in GMP, forms the quadratic `A x^2 + B x + C = 0` for the purely-periodic tail (with `A=k_{k-1}`, `B=k_{k-2}-h_{k-1}`, `C=-h_{k-2}`), solves via the discriminant with largest-square extraction (`fcf_extract_square`, trial division up to `10^6`), then applies the leading terms as a Möbius transform `(Hx+H')/(Kx+K')` and rationalises to `(P + Q Sqrt[R])/S` (`fcf_qirr_to_expr`).

**Data structures.** Convergents are GMP `mpz_t` registers; symbolic non-periodic convergents are `Expr` trees via `eval_and_free`. Output is an Integer/Rational or a rationalised quadratic-irrational expression.

**Complexity / limits.** Linear in the number of terms. The periodic path requires exact integer terms; a residual `p^2 q` with both `p, q > 10^6` is left un-reduced (astronomically unlikely for reconstructed CF data).

- `Protected` (not `Listable` — the argument is the whole term list).
- The `ai` of the finite form may be **symbolic**; the result is the convergent

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/contfrac.c`](https://github.com/stblake/mathilda/blob/main/src/contfrac.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
