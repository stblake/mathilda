---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Moebius function and Moebius inversion (Chapter 16)."
---

### Worked examples

```mathematica
In[1]:= MoebiusMu[30]  (* A square-free integer with three prime factors gives (-1)^3 *)
```

```mathematica
In[1]:= MoebiusMu[12]  (* Any squared prime factor sends the value to 0 *)
```

```mathematica
In[1]:= MoebiusMu[1]  (* The empty product is 1 *)
```

```mathematica
In[1]:= MoebiusMu[Range[12]]  (* Threaded element-wise over a range by Listable *)
```

```mathematica
In[1]:= Table[Total[MoebiusMu[Divisors[n]]], {n, 1, 8}]  (* Summing mu over the divisors of n is 1 only at n = 1 -- the identity behind Moebius inversion *)
```

```mathematica
In[1]:= MoebiusMu[1 + I]  (* A non-real argument is factored over the Gaussian integers *)
```

### Notes

`MoebiusMu[n]` is `0` when `n` has a squared prime factor, and otherwise `(-1)^k` for `k`
distinct primes, so it detects square-freeness (`MoebiusMu[n] != 0` iff `n` is square-free,
the test [`SquareFreeQ`](SquareFreeQ.md) makes directly). The sign of `n` is ignored. The
summatory identity `sum_{d | n} mu(d) = [n == 1]` is the kernel of Moebius inversion, which
recovers an arithmetic function from its divisor sums.

`MoebiusMu` is the only member of its factor-counting family — alongside
[`PrimeNu`](PrimeNu.md), [`PrimeOmega`](PrimeOmega.md) and
[`LiouvilleLambda`](LiouvilleLambda.md) — with a machine-integer array kernel, so a packed
or `NDArray` `int64` argument is factored element-wise on the buffer, and the head also
lowers inside `Compile[]` at a rank-1 integer-array shape. A non-real Gaussian-integer
argument, or `GaussianIntegers -> True`, switches the factorisation to `Z[i]`. `MoebiusMu[0]`
is left unevaluated.
