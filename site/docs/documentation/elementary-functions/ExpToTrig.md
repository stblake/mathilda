# ExpToTrig

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpToTrig[expr]
    rewrites exponentials and logarithms in expr in terms of circular and
    hyperbolic trigonometric functions when possible.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]

In[2]:= ExpToTrig[Log[1 + I x] - Log[1 - I x]]
Out[2]= (2*I) ArcTan[x]

In[3]:= ExpToTrig[Exp[I x] == -1]
Out[3]= Cos[x] + I Sin[x] == -1
```

## Implementation notes

**Algorithm.** `builtin_exptotrig` is the partial inverse of TrigToExp. It runs
a four-stage pipeline on the argument: (1) `ReplaceRepeated` with
`exp_to_trig_rules` — `E^(I x) :> Cos[x] + I Sin[x]`, `E^(-x) :> Cosh[x] -
Sinh[x]`, `E^x :> Cosh[x] + Sinh[x]`, plus `Log`-combination patterns that fold
back to `ArcTan`/`ArcTanh`/`ArcSin`/`ArcCsch`/... ; (2) `Together` to combine
over a common denominator; (3) `Cancel` to reduce the rational; (4) a final
`ReplaceRepeated` with `exp_to_trig_simp` that recovers reciprocal heads
(`Sin/Cos :> Tan`, `Cos^-1 :> Sec`, etc.). The trig canonicalizer is suppressed
across the whole pipeline (`trig_canon_suppress_inc`/`dec`) so the intermediate
`Sin/Cos` forms survive long enough for `Together`/`Cancel` to act.

**Data structures.** Three static rule lists (`exp_to_trig_rules`,
`exp_to_trig_simp`, plus the shared `trig_factor_*` machinery) parsed once in
`trigsimp_init`. Results are routed through the active `FactorMemo` by the public
`builtin_exptotrig` wrapper path via `trig_memo_call` semantics shared with the
other trig builtins.

- `Listable`, `Protected`.
- Tries when possible to give results that do not involve explicit complex numbers.
- ExpToTrig is natively the inverse of `TrigToExp`.
- Automatically threads over lists, equations, inequalities, and logic functions.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ExpToTrig[Exp[x]]
Out[1]= Cosh[x] + Sinh[x]
```

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]
```

```mathematica
In[1]:= ExpToTrig[(E^x - E^(-x))/2]
Out[1]= Sinh[x]
```

```mathematica
In[1]:= ExpToTrig[E^(I x) + E^(-I x)]
Out[1]= 2 Cos[x]
```

### Notes

`ExpToTrig` rewrites exponentials in terms of circular and hyperbolic
functions: a real exponent yields `Cosh + Sinh`, an imaginary one yields
`Cos + I Sin`. It recognises the classical combinations, folding the
half-difference `(E^x - E^(-x))/2` back to `Sinh[x]` and the sum
`E^(I x) + E^(-I x)` to `2 Cos[x]`.
