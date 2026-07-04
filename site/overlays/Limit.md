---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 3."
---
### Worked examples

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1
```

```mathematica
In[1]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[1]= 2
```

```mathematica
In[1]:= Limit[(1 + a/x)^x, x -> Infinity]
Out[1]= E^a
```

```mathematica
In[1]:= Limit[(Sin[x] - x + x^3/6)/x^5, x -> 0]
Out[1]= 1/120
```

```mathematica
In[1]:= Limit[(x^x - x)/(1 - x + Log[x]), x -> 1]
Out[1]= -2
```

```mathematica
In[1]:= Limit[x - Sqrt[x^2 + x], x -> Infinity]
Out[1]= -1/2
```

```mathematica
In[1]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[1]= 5
```

```mathematica
In[1]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[1]= Infinity

In[2]:= Limit[1/x, x -> 0, Direction -> "FromBelow"]
Out[2]= -Infinity
```

### Choosing a method

`Limit` evaluates by running a cascade of strategy layers in turn; each either
resolves the limit or hands the problem to the next. `Method -> m` restricts the
top-level call to a single strategy group. The default `Method -> Automatic`
runs the whole cascade in the order below. A named method computes *only* that
group; if it does not apply to the given expression, the `Limit` is left
unevaluated (an unrecognised method name is reported and also left unevaluated).

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[1]= 1

In[2]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[2]= 2

In[3]:= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
Out[3]= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
```

| `Method` | Strategy | Typical use |
|----------|----------|-------------|
| `Automatic` | run every strategy below, in order | default — best all-rounder |
| `"Substitution"` | continuity / direct substitution (via `Together`), `Abs` kink resolution, atom-substitution and one-sided probes | removable singularities, `Abs`, essential-singularity ratios |
| `"RationalFunction"` | leading-degree comparison for `P(x)/Q(x)` | rational functions at a point or at `Infinity` |
| `"Series"` | Taylor / Laurent / Puiseux expansion, reading the leading term | the workhorse — most `0/0` and `∞/∞` forms |
| `"LHospital"` | L'Hospital's rule with growth guardrails | `0/0`, `∞/∞` where `Series` cannot expand |
| `"Asymptotic"` | dominant-term / `Log` / exponential reductions at infinity, including `f^g` via `Exp[g Log f]` | limits at `Infinity`, `(1 + a/x)^x`, `Log`-of-sum |
| `"Bounded"` | squeeze / bounded-envelope to 0 and bounded-oscillation `Interval` returns | `Sin[x^2]/x`, oscillatory numerators |

The method restriction applies only to the outermost call: recursive
sub-limits — one-sided probes, L'Hospital iterations, `Abs` splitting — always
run the full cascade, so e.g. `Method -> "Series"` still resolves a two-sided
pole by falling back to its one-sided branches.

### Notes

`Limit[f, x -> a]` resolves the standard removable-singularity and indeterminate forms, including the classic `(1 + 1/x)^x -> E` and `0/0` cancellations such as `(x^2 - 1)/(x - 1)`. The `Direction` option selects one-sided (`"FromAbove"`/`"FromBelow"`) or complex approaches; the default is two-sided. The `Method` option (see above) selects a specific internal strategy, defaulting to `Automatic`. Results may be a finite value, `Infinity`, `ComplexInfinity`, `Indeterminate`, an `Interval`, or the original expression unevaluated when the limit cannot be determined. Iterated and joint multivariate limits are supported through the list forms of the second argument.
