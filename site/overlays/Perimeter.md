### Worked examples

```mathematica
In[1]:= Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]  (* an exact perimeter keeps the radical *)
```

```mathematica
In[1]:= Perimeter[Polygon[{{0, 0}, {1/2, 0}, {0, 1}}]]  (* rational coordinates, a canonicalised Sqrt sum *)
```

```mathematica
In[1]:= Perimeter[Polygon[{{0, 0}, {3., 0}, {3., 4.}}]]  (* a real coordinate gives a machine length *)
```

### Notes

`Perimeter[Polygon[{{x1, y1}, …}]]` is the sum of the edge lengths, including the closing edge
from the last vertex back to the first.

On exact (`Integer`/`Rational`) coordinates each edge length is a `Sqrt` of an exact squared
distance, and the sum is handed to the evaluator, so radical canonicalisation produces the
Wolfram-shaped exact answer (`2 + Sqrt[2]`, `3/2 + 1/2 Sqrt[5]`, …). Any `Real` coordinate
makes the whole result a machine `Real`, summed with `hypot`.

Like `Area`, consecutive and wrap-around duplicate vertices are collapsed first, and a polygon
with fewer than three distinct vertices is `Undefined`. Simple polygons only; `Protected`, not
`Listable`.
