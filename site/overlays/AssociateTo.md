### Worked examples

```mathematica
In[1]:= a = <|x -> 1, y -> 2|>;
```

```mathematica
In[1]:= AssociateTo[a, z -> 3]  (* a new key is appended at the end *)
```

```mathematica
In[1]:= AssociateTo[a, x -> 9]  (* an existing key is updated in place *)
```

```mathematica
In[1]:= a
```

### Notes

`AssociateTo[s, key -> val]` is the in-place analogue of `AppendTo` for
associations: it reads the current association stored in `s`, inserts or updates
the entry, assigns the result back to `s`, and returns it. The second argument
may also be a list of rules. Because the key set is de-duplicated with the usual
*first position, last value* rule, associating an existing key overwrites its
value without moving it, while a genuinely new key is added at the end. `s` must
already hold an association.
