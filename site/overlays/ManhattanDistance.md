### Worked examples

```mathematica
In[1]:= ManhattanDistance[{0, 0}, {3, 4}]  (* sum of coordinate gaps: 3 + 4 = 7 *)
```

```mathematica
In[1]:= ManhattanDistance[{a}, {b}]  (* a symbolic pair stays symbolic *)
```

### Notes

`ManhattanDistance[u, v]` is `Sum Abs[u_i - v_i]`, the taxicab (city-block)
distance — the distance travelled along axis-aligned moves rather than in a
straight line. Both arguments must be scalars or lists of equal length.

Exact input stays exact, complex components contribute their modulus, and a
symbolic pair returns `Abs[a - b]` rather than an error, exactly as Mathematica
answers. Note that a bare rational *scalar* such as `ManhattanDistance[1/2, 3/2]`
is left unevaluated, because a `Rational` is a compound expression rather than an
atomic scalar — wrap the values in lists (`ManhattanDistance[{1/2}, {3/2}]`) to
measure them.
