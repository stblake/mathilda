### Worked examples

```mathematica
In[1]:= Cross[{1,0,0},{0,1,0}]
Out[1]= {0, 0, 1}
```

```mathematica
In[1]:= Cross[{a1,a2,a3},{b1,b2,b3}]
Out[1]= {-a3 b2 + a2 b3, -(-a3 b1 + a1 b3), -a2 b1 + a1 b2}
```

```mathematica
In[1]:= Cross[{2,1,-1},{1,-1,2}]
Out[1]= {1, -5, -3}
```

```mathematica
In[1]:= Cross[{1,2}]
Out[1]= {-2, 1}
```

```mathematica
In[1]:= Cross[{1,2,3,4},{5,6,7,8},{9,10,11,13}]
Out[1]= {4, -8, 4, 0}
```

### Notes

`Cross[a, b]` is the usual 3-vector cross product, returned in fully symbolic cofactor form so it works for indeterminate components. The general `Cross[a1, ..., a(n-1)]` form gives the unique vector in *n* dimensions orthogonal to all *n*-1 inputs, computed as the signed cofactor minors of the matrix whose remaining row is the basis vectors — so `Cross[{1,2}]` in the plane returns the 90-degree rotation `{-2, 1}`, and three 4-vectors yield a fourth vector perpendicular to all of them.
