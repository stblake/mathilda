### Worked examples

```mathematica
In[1]:= Function[x, If[x > 0, Return[pos], Return[neg]]][3]  (* Return exits the enclosing Function body *)
```

```mathematica
In[1]:= h[n_] := Module[{s = 0}, Do[s += i; If[s > 10, Return[i]], {i, 1, n}]]; h[10]  (* breaks out of a loop inside a Module *)
```

```mathematica
In[1]:= Module[{}, Do[Return[5, Block], {3}]]  (* the two-argument form targets a named boundary; with no Block here it survives *)
```

### Notes

`Return[expr]` yields `expr` from the innermost enclosing scope or loop boundary —
`Function`, `Module`, `Block`, `With`, `Do`, `For` or `While`. `Return[]` is
shorthand for `Return[Null]`.

`CompoundExpression` and the `Hold`-free heads (`If`, `Which`, `Switch`, …) let
the marker bubble through unchanged so it can reach the enclosing boundary. The
two-argument `Return[expr, h]` skips past intervening boundaries to the nearest
one whose head is `h`; if none matches, the marker survives at top level as a
literal expression (as in the third example, where there is no enclosing `Block`).
