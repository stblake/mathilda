### Worked examples

```mathematica
In[1]:= FromCoefficientRules[{{2, 0} -> 1, {1, 1} -> 2, {0, 2} -> 1}, {x, y}]  (* rebuild the polynomial *)
```

```mathematica
In[2]:= FromCoefficientRules[{{2} -> 3, {0} -> 1}, x]  (* a single variable *)
```

```mathematica
In[3]:= FromCoefficientRules[CoefficientRules[1 + x^3 + 7 x y^2, {x, y}], {x, y}]  (* round-trips CoefficientRules *)
```

### Notes

`FromCoefficientRules[{expvec -> coeff, ...}, {x1, ..., xk}]` reconstructs the
polynomial `Sum[coeff x1^e1 ... xk^ek]` — the inverse of `CoefficientRules`. Each
exponent vector must have exactly `k` integer components (matching the variable
list); the variable list may also be a single bare variable. An exponent of `0`
drops the variable and an exponent of `1` drops the `Power` wrapper, and an empty
rule list reconstructs `0`. A malformed rule (wrong length or a non-integer
exponent) leaves the call unevaluated.
