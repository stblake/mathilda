### Worked examples

The generalized outer product — every pairing of elements fed to `f`:

```mathematica
In[1]:= Outer[Times, {1, 2, 3}, {a, b, c}]
Out[1]= {{a, b, c}, {2 a, 2 b, 2 c}, {3 a, 3 b, 3 c}}
```

With a symbolic `f` you get the full table of combinations:

```mathematica
In[1]:= Outer[f, {1, 2}, {x, y}]
Out[1]= {{f[1, x], f[1, y]}, {f[2, x], f[2, y]}}
```

`Outer[Times, v, v]` builds the rank-one matrix `v v^T`. Taking the basis
`{1, x, x^2}` produces the symmetric Hankel-style table of monomial products —
the kind of structure that underlies Vandermonde and Gram constructions:

```mathematica
In[1]:= Outer[Times, {1, x, x^2}, {1, x, x^2}]
Out[1]= {{1, x, x^2}, {x, x^2, x^3}, {x^2, x^3, x^4}}
```

Every such outer product of two vectors is rank one, so its determinant must
vanish — a fact the symbolic linear algebra confirms exactly:

```mathematica
In[1]:= Det[Outer[Times, {a, b, c}, {1, 1, 1}]]
Out[1]= 0
```

Outer products of three or more lists nest one level deeper per list:

```mathematica
In[1]:= Outer[List, {1, 2}, {a, b}, {X, Y}]
Out[1]= {{{{1, a, X}, {1, a, Y}}, {{1, b, X}, {1, b, Y}}}, {{{2, a, X}, {2, a, Y}}, {{2, b, X}, {2, b, Y}}}}
```

### Notes

`Outer[f, l1, l2, ...]` forms every combination of lowest-level elements, one
from each list, and applies `f` to it. The result has nesting depth equal to the
combined depth of the inputs, so `Outer` of two vectors is a matrix, of three is
a rank-3 array, and so on. With `f = Times` it is the tensor (outer) product;
with a symbolic head it is a complete combination table. The optional level
arguments `Outer[f, l1, l2, ..., n]` (or per-list `n1, n2, ...`) control which
sublists are treated as the separate elements to combine.
