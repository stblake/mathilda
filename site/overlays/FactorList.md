### Worked examples

```mathematica
In[1]:= FactorList[x^2 - 1]  (* each factor with its multiplicity *)
```

```mathematica
In[2]:= FactorList[2 x^2 + 4 x + 2]  (* the overall numerical factor leads *)
```

```mathematica
In[3]:= FactorList[x^4 - 2, Extension -> Sqrt[2]]  (* options are forwarded to Factor *)
```

### Notes

`FactorList[poly]` gives the irreducible factors of `poly` as `{factor,
exponent}` pairs. It is a thin wrapper over `Factor` — every option
(`GaussianIntegers -> True`, `Extension -> {a1, ...}`) is forwarded verbatim. The
first pair is always the overall numerical factor `{c, 1}` (it is `{1, 1}` when
there is none), and the denominator factors of a rational function appear with
negative exponents. An irreducible radical factor such as `Sqrt[base]` carries
exponent `1`, since only an integer `Power` exponent is read as a multiplicity.
