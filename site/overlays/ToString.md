### Worked examples

```mathematica
In[1]:= ToString[x^2 + y^3]  (* InputForm by default *)
```

```mathematica
In[1]:= ToString[x^2 + y^3, FullForm]
```

```mathematica
In[1]:= ToString[x^2 + y^3, TeXForm]
```

```mathematica
In[1]:= Head[ToString[42]]  (* the result is always a String *)
```

```mathematica
In[1]:= ToString[C[10], TeXForm]  (* a generated constant becomes a subscript *)
```

### Notes

`ToString[expr]` returns a `String` holding the printed form of `expr` in
`InputForm`; `ToString[expr, form]` selects the form, with `FullForm` and
`TeXForm` supported and `StandardForm` / `OutputForm` accepted as aliases of
`InputForm`. All of the formatting is shared with the standard printer, so the
string matches what the REPL would display.

An unsupported form leaves the call unevaluated (e.g. `ToString[x, FooForm]`),
making a typo visible at the call site rather than silently downgrading it. Under
`TeXForm` the generated constants `C[k]` from `DSolve` / `Reduce` / `Integrate`
render as the subscripted `c_k` (single-character subscripts bare, longer ones
braced), matching Mathematica and the notebook LaTeX renderer.
