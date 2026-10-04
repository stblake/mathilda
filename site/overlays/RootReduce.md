### Worked examples

```mathematica
In[1]:= RootReduce[Sqrt[2] + Sqrt[3]]  (* a degree-4 algebraic number: a Root object *)
```

```mathematica
In[1]:= RootReduce[Sqrt[8]]  (* pulled apart into 2 Sqrt[2] *)
```

```mathematica
In[1]:= RootReduce[Sqrt[2] Sqrt[3]]  (* combined into Sqrt[6] *)
```

```mathematica
In[1]:= RootReduce[Sqrt[3 + 2 Sqrt[2]]]  (* a nested radical, denested *)
```

```mathematica
In[1]:= RootReduce[(1 + Sqrt[5])/2 - GoldenRatio]  (* decided to be exactly zero *)
```

```mathematica
In[1]:= RootReduce[Sqrt[2] + Sqrt[3] == Sqrt[5 + 2 Sqrt[6]]]  (* an equality, decided exactly *)
```

### Notes

`RootReduce[expr]` canonicalises an algebraic expression to a single
representative. A constant algebraic number becomes a rational, a quadratic
radical `(a + b Sqrt[c])/q`, or a `Root[poly &, k]` object for degree three and
up; the representative is unique, so two spellings of the same number reduce to
the same thing and their difference reduces to `0`.

Reduction runs on FLINT's exact `qqbar` engine — minimal polynomial plus an
isolating enclosure, with **no numeric zero oracle** — so `RootReduce` also
*decides* equations and inequalities between constant algebraic numbers exactly,
and threads over lists, `Solve`-result rules, and logical combinations. It
leaves anything with no algebraic content unchanged. Option
`Method -> "Automatic" | "Recursive" | "NumberField"`.
