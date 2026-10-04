### Worked examples

```mathematica
In[1]:= HurwitzZeta[2, 1]  (* reduces to Zeta[2] when the second argument is 1 *)
```

```mathematica
In[1]:= HurwitzZeta[s, 1/2]  (* the half-integer shift: (2^s - 1) Zeta[s] *)
```

```mathematica
In[1]:= HurwitzZeta[0, 0]  (* read off the Bernoulli polynomial, -BernoulliB[1, 0] *)
```

```mathematica
In[1]:= HurwitzZeta[4, 5]  (* positive integer second argument: Zeta[4] minus a finite power sum *)
```

```mathematica
In[1]:= N[HurwitzZeta[3, 1/4], 20]  (* arbitrary precision via Euler-Maclaurin *)
```

```mathematica
In[1]:= HurwitzZeta[2, 1.5]  (* inexact argument triggers the numeric kernel *)
```

### Notes

`HurwitzZeta[s, a] = Sum_{k>=0} (k+a)^{-s}` (continued to all `s != 1`). Unlike the
two-argument `Zeta`, it sums the *principal-branch* powers `(k+a)^{-s}`, so it
retains poles at `a = 0, -1, -2, …` and disagrees with `Zeta[s, a]` for
non-positive real `a`.

Exact reductions cover `s == 1` (`ComplexInfinity`), `a == 1` (`Zeta[s]`),
`a == 1/2` (`(2^s - 1) Zeta[s]`), positive integer `a`
(`Zeta[s] - Sum_{k=1}^{a-1} k^{-s}`), and non-positive integer `s` at integer `a`
(the entire value `-BernoulliB[1-s, a]/(1-s)`). Any inexact operand routes to an
Euler–Maclaurin complex-MPFR kernel (MPFR has no native Hurwitz zeta).
`HurwitzZeta` carries a real `NDArray` kernel and lowers under `Compile[]` at both
scalar and rank-1 array shapes.
