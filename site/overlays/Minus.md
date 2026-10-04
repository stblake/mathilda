### Worked examples

```mathematica
In[1]:= Minus[5]  (* unary negation, -x *)
```

```mathematica
In[1]:= Minus[a + b]  (* -1 distributes through the sum *)
```

```mathematica
In[1]:= Minus[2 + 3 I]  (* negates a complex number *)
```

```mathematica
In[1]:= Minus[{1, -2, 3}]  (* Listable: threads over the list *)
```

```mathematica
In[1]:= SortBy[{3, -1, 2, -5}, Minus]  (* as a sort key, orders descending *)
```

```mathematica
In[1]:= FullForm[Minus[x]]  (* the single rewrite Minus[x] -> Times[-1, x] *)
```

### Notes

`Minus[x]` is the functional form of unary negation, equivalent to `-x`. The
parser already turns the `-` operator into `Times[-1, x]`, so `Minus` exists
mainly as a head you can pass by name: `SortBy[list, Minus]` and
`KeySortBy[assoc, Minus]` sort by the negated value (i.e. descending), and
`Map[Minus, list]` negates each element. The rewrite to `Times[-1, x]` means
`Minus` inherits everything `Times` can do — distributing over a `Plus`,
negating a `Complex`, and operating on packed or `NDArray` arguments directly on
the buffer. It also lowers inside `Compile[]`.

`Minus` is strictly unary. `Minus[x, y]` (or any count other than one) emits the
message `Minus::argx` and is left unevaluated; use `Subtract[x, y]` (or `x - y`)
for a binary difference.
