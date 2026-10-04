### Worked examples

```mathematica
In[1]:= s = 0; Do[If[i > 3, Break[]]; s = s + i, {i, 10}]; s  (* the loop stops the moment the guard fires *)
```

```mathematica
In[1]:= For[i = 1, i <= 10, i++, If[i == 4, Break[]]]; i  (* Break escapes the loop, leaving i at its exit value *)
```

```mathematica
In[1]:= n = 1; While[True, If[n > 5, Break[]]; n++]; n  (* the standard way out of a While[True] loop *)
```

```mathematica
In[1]:= Break[]  (* outside any loop it is inert: reported and wrapped in Hold *)
```

### Notes

`Break[]` takes effect as soon as it is evaluated — even from inside an `If`
within the body — and escapes only the *innermost* enclosing `Do`, `For` or
`While`, which then yields `Null`.

It is recognised by its head at the loop boundary, not short-circuited the way a
`Throw` is, so `Print[Break[]]` does **not** escape the loop. `Table` deliberately
does not honour `Break`. A `Break[]` that reaches top level with no loop to
consume it prints `Break::nofwd` and comes back as the inert `Hold[Break[]]`, so
feeding that result back in does not re-trigger anything.
