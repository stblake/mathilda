---
references:
  - "G. E. Andrews, *The Theory of Partitions*, Cambridge University Press, 1998 — the standard reference on partitions and their generating functions."
---

### Partitions of an integer

A *partition* of `n` is a way of writing it as a sum of positive integers, order disregarded
— pictured as a *Young diagram* of left-justified rows. `IntegerPartitions[n]` lists them in
reverse-lexicographic order; the second and later arguments restrict the number of parts and
the allowed parts. The count of unrestricted partitions is [`PartitionsP`](PartitionsP.md),
so `Length[IntegerPartitions[n]] == PartitionsP[n]`.

### Worked examples

```mathematica
In[1]:= IntegerPartitions[4]
Out[1]= {{4}, {3, 1}, {2, 2}, {2, 1, 1}, {1, 1, 1, 1}}
```

```mathematica
In[1]:= Length[IntegerPartitions[10]]
Out[1]= 42
```

There are five partitions of `4`; the ten has `42` partitions, matching `PartitionsP[10]`.
