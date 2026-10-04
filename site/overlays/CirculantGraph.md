### Worked examples

```mathematica
In[1]:= EdgeList[CirculantGraph[6, 1]]  (* a single offset of 1 gives the 6-cycle *)
```

```mathematica
In[1]:= EdgeList[CirculantGraph[6, 3]]  (* the offset n/2 joins antipodes, each edge counted once *)
```

```mathematica
In[1]:= EdgeList[CirculantGraph[5, -1]]  (* a negative offset is reduced mod n, so this is the 5-cycle again *)
```

```mathematica
In[1]:= EdgeList[CirculantGraph[8, {1, 2}]]  (* a list of offsets gives the square of the 8-cycle *)
```

```mathematica
In[1]:= VertexDegree[CirculantGraph[10, {1, 3}]]  (* each distinct offset contributes two to every degree *)
```

```mathematica
In[1]:= EdgeCount[CirculantGraph[12, {1, 2, 3}]]  (* the count is n times the number of offsets when none is n/2 *)
```

### Notes

`CirculantGraph[n, j]` joins `i` to `i ± j (mod n)` on vertices `1..n`; the offset may be one integer or a list. Offsets are reduced modulo `n`, so `j` and `n - j` give the same graph, and an offset that is a multiple of `n` adds nothing.

The result is an undirected graph whose edge list is sorted lexicographically. Non-integer arguments, or a non-positive `n`, leave the call unevaluated.
