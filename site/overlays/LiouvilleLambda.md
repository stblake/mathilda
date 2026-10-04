---
references:
  - "G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the Liouville function lambda and its summatory divisor identity (Chapter 17, 22)."
---

### Worked examples

```mathematica
In[1]:= LiouvilleLambda[12]  (* lambda(n) = (-1)^Omega(n); 12 = 2^2 * 3 has Omega = 3 *)
```

```mathematica
In[1]:= LiouvilleLambda[210]  (* With 210 = 2 * 3 * 5 * 7, Omega = 4, so lambda = 1 *)
```

```mathematica
In[1]:= LiouvilleLambda[1]  (* The empty product is 1 *)
```

```mathematica
In[1]:= LiouvilleLambda[Range[10]]  (* Threaded element-wise over a range by Listable *)
```

```mathematica
In[1]:= LiouvilleLambda[6 * 35] == LiouvilleLambda[6] LiouvilleLambda[35]  (* Complete multiplicativity, with no coprimality needed *)
```

```mathematica
In[1]:= Table[Total[LiouvilleLambda[Divisors[n]]], {n, 1, 10}]  (* Summing lambda over the divisors of n is 1 exactly when n is a perfect square *)
```

### Notes

The Liouville function `LiouvilleLambda[n] = (-1)^Omega(n)` is `+1` when `n` has an even
number of prime factors counted with multiplicity and `-1` when odd, so it equals
`(-1)^PrimeOmega[n]` (see [`PrimeOmega`](PrimeOmega.md)). It is *completely* multiplicative —
`lambda(ab) = lambda(a) lambda(b)` for all `a, b`, unlike [`MoebiusMu`](MoebiusMu.md), which
is only multiplicative over coprime arguments — and the sign of `n` is ignored.

Its summatory divisor identity `sum_{d | n} lambda(d) = [n is a perfect square]` is the
analogue of the Moebius identity and detects squares. A non-real Gaussian-integer argument,
or `GaussianIntegers -> True`, counts Gaussian prime factors over `Z[i]`.
`LiouvilleLambda[0]` is left unevaluated.
