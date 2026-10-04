### Worked examples

```mathematica
In[1]:= Counts[{a, b, a, c, b, a}]
```

```mathematica
In[1]:= Counts[{1, 1, 2, 3, 3, 3}]
```

```mathematica
In[1]:= Counts[Characters["mississippi"]]  (* letter frequencies *)
```

### Notes

`Counts[list]` returns `<|element -> count, ...|>` in order of first appearance —
the association form of `Tally`, which gives the same information as a list of
`{element, count}` pairs. It is the standard histogram primitive: counting
characters, residues, or category labels. A packed integer or real buffer is
counted on its machine words (through `Tally`'s direct-indexed or hashed count)
and relabelled as rules, so large numeric data stay on the buffer. Use `CountsBy`
to count by a function of each element.
