---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational function arithmetic and common denominators."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on polynomial GCDs used in cancellation."
---
### Worked examples

```mathematica
In[1]:= Together[1/x + 1/y]
Out[1]= (x + y)/(x y)
```

```mathematica
In[1]:= Together[1/(x-1) + 1/(x+1)]
Out[1]= (2 x)/(-1 + x^2)
```

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d)/(b d)
```

### Notes

Together combines a sum of fractions over a single common denominator, the inverse
operation of `Apart`. It computes the least common denominator and cancels common
polynomial factors via GCD, so `1/(x-1) + 1/(x+1)` collapses to `(2 x)/(x^2-1)`
with the cross terms cleared. Distinct symbolic denominators are simply multiplied,
giving `(b c + a d)/(b d)` for `a/b + c/d`. Unlike a raw expansion, Together keeps
the denominator factored where possible and does not multiply out `Power[Plus, n]`
factors unnecessarily.
