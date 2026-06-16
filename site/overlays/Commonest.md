### Worked examples

```mathematica
In[1]:= Commonest[{1, 2, 2, 3, 3, 3, 4}]
Out[1]= {3}
```

```mathematica
In[1]:= Commonest[{a, b, a, c, b, a, d, b}, 2]
Out[1]= {a, b}
```

```mathematica
In[1]:= Commonest[Table[Mod[k^2, 7], {k, 0, 20}]]
Out[1]= {1, 4, 2}
```

### Notes

`Commonest[list]` returns every element tied for the highest frequency; `Commonest[list, n]` returns the `n` most common. The quadratic-residue example above shows the three nonzero residues mod 7 each occurring equally often.
