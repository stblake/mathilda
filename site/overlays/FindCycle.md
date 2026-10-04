### Worked examples

```mathematica
In[1]:= FindCycle[CycleGraph[5]]  (* the only cycle of a pentagon *)
```

```mathematica
In[1]:= FindCycle[PathGraph[{1, 2, 3}]]  (* a tree has no cycle, so the result is empty *)
```

```mathematica
In[1]:= FindCycle[CompleteGraph[4], {4}]  (* a cycle of exactly four edges *)
```

```mathematica
In[1]:= FindCycle[CompleteGraph[4], {3}, All]  (* every triangle, each reported once *)
```

```mathematica
In[1]:= FindCycle[{CycleGraph[5], 2}]  (* a cycle through a given vertex *)
```

```mathematica
In[1]:= FindCycle[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* directed cycle *)
```

```mathematica
In[1]:= FindCycle[CompleteGraph[4], Infinity, 2]  (* up to two cycles of any length *)
```

### Notes

A cycle is returned as a list of edges, and the result is a list of cycles: `{}` when there is none. A length specification counts edges and may be `k` (at most), `{k}` (exactly), `{kmin, kmax}` or `Infinity`; a third argument asks for up to `n` cycles or `All`.

The plain form is linear time. Exact-length and enumerating forms are exponential in the worst case, so they give up under a step budget and then stay unevaluated. For those forms, which cycles are reported first is Mathilda's own deterministic choice.
