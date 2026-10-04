### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomSample[{a, b, c, d, e}]  (* a random permutation: no element repeats *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomSample[Range[10], 5]  (* five distinct draws from 1..10 *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomSample[{a, b, c, d, e}, UpTo[3]]  (* UpTo clamps to the length *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomSample[{1, 1, 1, 10} -> {a, b, c, d}, 2]  (* weighted, without replacement *)
```

### Notes

`RandomSample` selects **without replacement**: no element is chosen twice, so the
count may not exceed the list length (use `UpTo[n]`, which clamps). With no count
it returns a full random permutation of the list.

The uniform form runs a partial Fisher–Yates shuffle — each of the first `n` slots
is swapped with a uniformly chosen later slot — so it is `O(n)` even when drawing a
few elements from a very long list. The weighted form `{w1, ...} -> {e1, ...}`
draws sequentially by inverse-CDF, zeroing each chosen weight so it cannot recur.

All draws come from the stream `SeedRandom` controls, so the samples above are
reproducible across runs.
