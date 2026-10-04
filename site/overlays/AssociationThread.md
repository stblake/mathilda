### Worked examples

```mathematica
In[1]:= AssociationThread[{a, b, c}, {1, 2, 3}]
```

```mathematica
In[1]:= AssociationThread[{a, b} -> {1, 2}]  (* the keys -> values rule form *)
```

### Notes

`AssociationThread[keys, values]` zips two equal-length lists into an
association, pairing `keys[[i]]` with `values[[i]]`. The single-argument rule form
`AssociationThread[keys -> values]` is equivalent. It is the inverse of taking
`Keys` and `Values` apart, and the quickest way to build an association from two
parallel lists. Duplicate keys collapse with the usual *first position, last
value* rule; the two lists must have the same length.
