### Worked examples

```mathematica
In[1]:= $Assumptions  (* the default: no assumptions in force *)
```

```mathematica
In[1]:= Assuming[x > 0, Refine[Sqrt[x^2]]]  (* Assuming extends $Assumptions for the body *)
```

```mathematica
In[1]:= Refine[Sqrt[x^2], x > 0]  (* a positional assumption is conjoined with $Assumptions *)
```

```mathematica
In[1]:= Assuming[Element[n, Integers], Simplify[Sin[n Pi]]]  (* Simplify reads the extended $Assumptions *)
```

```mathematica
In[1]:= Block[{$Assumptions = x > 0}, Refine[Abs[x]]]  (* a direct rebinding, which is what Assuming desugars to *)
```

```mathematica
In[1]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]  (* conjunction of facts, consumed per symbol *)
```

### Notes

`$Assumptions` is a global symbol — the default setting for the `Assumptions` option used by
`Simplify`, `Refine`, `Element` and the other functions that take assumptions. It is not a
builtin; it carries an `OwnValue` that defaults to `True` (no assumptions). The consumers
read its value directly rather than evaluating it, so a bound fact such as
`Element[x, Reals]` does not re-fire while being read.

`Assuming[fact, body]` temporarily extends it — effectively
`Block[{$Assumptions = $Assumptions && fact}, body]` — so nested `Assuming` calls compose and
the rebinding is restored on exit. The accumulated value is flattened into a per-variable fact
set when a consumer needs it. The precedence a consumer applies is uniform: a positional
assumption is conjoined with `$Assumptions`, whereas an explicit `Assumptions -> X` option
replaces it.
