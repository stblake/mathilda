# Integrate

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

```text
Integrate[f, x] gives the indefinite integral of f with respect to x.
Integrate[f, x, Method -> "<name>"] dispatches directly to a single
subroutine, bypassing the default cascade.  Accepted method names:
  "Automatic"          — try BronsteinRational, then RischNorman, then CRCTable (default)
  "BronsteinRational"  — Integrate`BronsteinRational (polynomial / rational)
  "DerivativeDivides"  — Integrate`DerivativeDivides (substitution u(x); direct + Eliminate/Solve)
  "LinearRadicals"     — Integrate`LinearRadicals (rationalise radicals of a x + b)
  "QuadraticRadicals"  — Integrate`QuadraticRadicals (Euler substitution for Sqrt[a x^2 + b x + c])
  "LinearRatioRadicals" — Integrate`LinearRatioRadicals (rationalise radicals of (a x + b)/(c x + d))
  "RischNorman"        — Integrate`RischNorman (Bronstein pmint heuristic)
  "CRCTable"           — Integrate`CRCTable (lazy-loaded CRC integral table)
  "Undefined"          — Integrate`Undefined (unknown functions u[x], u'[x]; Roach §1.7)
Named methods are strict: failure returns unevaluated, with no fallback.
The CRCTable rules are loaded from disk on first use only.
An applied 1-D InterpolatingFunction integrates to its antiderivative
InterpolatingFunction (mirroring D).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Integrate[3 + 5 x + 2 x^2, x]
Out[1]= 3 x + 5/2 x^2 + 2/3 x^3

In[2]:= Integrate[2 x/(x^2 + 1), x]
Out[2]= Log[1 + x^2]

In[3]:= Integrate[1/(x - a)^2, x]
Out[3]= -1/(-a + x)

In[4]:= Integrate[(2x+3)/(x^2+3x+5)^2, x]
Out[4]= -1/(5 + 3 x + x^2)

In[5]:= Integrate[1/((x-1)(x-2)(x-3)), x]          (* Phase 2 LRT closes this *)
Out[5]= -Log[-2 + x] + 1/2 Log[3 - 4 x + x^2]

In[6]:= Integrate[1/(x^2 + 1), x]                  (* Phase 4 LogToReal *)
Out[6]= ArcTan[x]

In[7]:= Integrate[1/(x^4 + x^2 + 1), x]            (* two quadratic factors *)
Out[7]= 1/4 Log[1 + x + x^2] + 1/2 ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] + 1/2 ArcTan[(1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/4 Log[1 - x + x^2]

In[8]:= Integrate[Sin[x], x, Method -> "RischNorman"]  (* strict, no fallback *)
Out[8]= -(1 + Cos[x])
```

## Implementation notes

- `Protected`, `Listable`.
- Eight-stage dispatch cascade (`DerivativeDivides`, `LinearRadicals`,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Partial** — implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## References

- Bronstein, "Symbolic Integration I: Transcendental Functions", 2nd ed. (Springer, 2005).
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 11–12.
- Source: [`src/calculus/integrate.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/integrate.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Integrate[1/(1 + x^2), x]
Out[1]= ArcTan[x]
```

```mathematica
In[1]:= Integrate[1/x, x]
Out[1]= Log[x]
```

```mathematica
In[1]:= Integrate[Cos[x], x]
Out[1]= Sin[x]
```

```mathematica
In[1]:= Integrate[x^3 + x, x]
Out[1]= 1/2 x^2 + 1/4 x^4
```

### Notes

`Integrate[f, x]` computes the indefinite integral via a cascade: Bronstein's rational-function algorithm, then the Risch–Norman (`pmint`) heuristic, then the lazy-loaded CRC integral tables; `Method -> "<name>"` pins a single subroutine. Antiderivatives are returned without an integration constant and are not always simplified — for example `Integrate[Sin[x], x]` returns `-(1 + Cos[x])` rather than `-Cos[x]`. Definite integration is **not** supported in this build: `Integrate[x^2, {x, 0, 1}]` threads `Integrate` over the bound list and returns the garbage form `{1/3 x^3, Integrate[x^2, 0], Integrate[x^2, 1]}` instead of `1/3`. Restrict use to the indefinite, single-variable form.
