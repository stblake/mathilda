### Worked examples

```mathematica
In[1]:= PowerExpand[Sqrt[a^2 b^2]]
Out[1]= a b
```

`PowerExpand` turns the logarithm of a product into a sum of logarithms:

```mathematica
In[1]:= PowerExpand[Log[a b c]]
Out[1]= Log[a] + Log[b] + Log[c]
```

It pulls exponents out of `Log` of a power:

```mathematica
In[1]:= PowerExpand[Log[x^n]]
Out[1]= n Log[x]
```

Restricting to a variable list expands only with respect to those variables (here `Sqrt[x^2] -> x`):

```mathematica
In[1]:= PowerExpand[Sqrt[x^2], {x}]
Out[1]= x
```

With `Assumptions -> True` the result is universally correct: the omitted branch-cut term reappears as an explicit `Floor[...]` correction:

```mathematica
In[1]:= PowerExpand[(a b)^n, Assumptions -> True]
Out[1]= a^n b^n E^((2*I) Pi Floor[1/2 - 1/2 (Arg[a] + Arg[b])/Pi] n)
```

### Notes

`PowerExpand[expr]` rewrites `(a b)^c` as `a^c b^c`, `(a^b)^c` as `a^(b c)`, and
expands `Log` and `Arg` of products and powers. These transformations are valid
in general only when `c` is an integer or `a, b` are positive reals; by default
`PowerExpand` disregards branch cuts. Use the variable-list form
`PowerExpand[expr, {x1, ...}]` to expand selectively, or
`Assumptions -> True` to obtain a branch-cut-correct result in which the
discarded `E^(2 Pi I Floor[...])` factor is made explicit.
