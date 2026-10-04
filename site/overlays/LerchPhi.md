### Worked examples

```mathematica
In[1]:= LerchPhi[1, 2, 1]  (* at z = 1 the Lerch transcendent is a Hurwitz zeta *)
```

```mathematica
In[1]:= LerchPhi[0, s, a]  (* only the k = 0 term survives *)
```

```mathematica
In[1]:= LerchPhi[z, 0, a]  (* a geometric sum, independent of a *)
```

```mathematica
In[1]:= LerchPhi[-1, 1, 1]  (* the alternating harmonic series *)
```

```mathematica
In[1]:= LerchPhi[z, -1, 1]  (* a negative integer s gives a rational function *)
```

```mathematica
In[1]:= LerchPhi[2, 3, 1]  (* reduces to a polylogarithm *)
```

```mathematica
In[1]:= N[LerchPhi[1/2, 2, 1], 20]
```

### Notes

`LerchPhi[z, s, a]` is the Lerch transcendent `Sum_{k>=0} z^k/(k+a)^s`, the
common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`:
`LerchPhi[1, s, a] = Zeta[s, a]` and `z LerchPhi[z, s, 1] = PolyLog[s, z]`.

Many argument shapes reduce in closed form — `z = 0` gives `a^-s`, `s = 0` gives
`1/(1-z)`, a negative integer `s` gives a rational function of `z`, and a
positive integer `a` shifts onto `PolyLog`. Otherwise a numeric value is summed
from the defining series for `|z| < 1`; `|z| > 1` uses an analytic continuation
within its domain and is otherwise left symbolic. The options
`DoublyInfinite -> True` and `IncludeSingularTerm -> True` select the two-sided
sum and the singular-term convention.
