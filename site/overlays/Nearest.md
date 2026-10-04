### Worked examples

```mathematica
In[1]:= Nearest[{1, 2, 3, 4, 5}, 2.3]  (* the single closest element *)
```

```mathematica
In[1]:= Nearest[{10, 20, 30, 40}, 25]  (* a tie returns both, in input order *)
```

```mathematica
In[1]:= Nearest[{3 + 4 I, 1}, 0]  (* complex elements compare by modulus *)
```

### Notes

`Nearest[list, x]` returns the element (or elements) of `list` closest to the
target `x`, using `Abs[element - x]` as the distance. **Every** element tied at
the minimum distance is returned, in the order it appears in the input, so the
result is always a list. `Nearest[{}, x]` is `{}`.

Every distance must be a real number, or the whole call is left unevaluated — a
symbolic element, a symbolic target, a non-real complex, or even a symbolic real
such as `Pi` all decline, rather than silently dropping out of the result.
Complex elements are handled through their modulus.
