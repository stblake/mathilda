# RandomSample

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomSample[{e1, e2, ...}, n]`**

gives a pseudorandom sample of n of the ei, without replacement.

**`RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]`**

gives a weighted pseudorandom sample of n of the ei.

**`RandomSample[{e1, e2, ...}]`**

gives a pseudorandom permutation of the ei.

**`RandomSample[list, UpTo[n]]`**

gives a sample of n of the ei, or as many as are available.

<details>
<summary>Notes</summary>

RandomSample never samples any element more than once. Use SeedRandom to seed the pseudorandom generator for reproducible results.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (10)

```mathematica
In[1]:= SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]
Out[1]= {e, c, a}

In[2]:= Sort[RandomSample[{1, 2, 3, 4, 5}, 5]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= Length[RandomSample[{a, b, c, d, e}]]
Out[3]= 5

In[4]:= RandomSample[{a, b, c}, 0]
Out[4]= {}

In[5]:= Length[RandomSample[{a, b, c, d, e}, UpTo[10]]]
Out[5]= 5

In[6]:= RandomSample[{1, 0, 0} -> {a, b, c}, 1]
Out[6]= {a}

In[7]:= Sort[RandomSample[{1, 1, 0} -> {a, b, c}, 2]]
Out[7]= {a, b}

In[8]:= Sort[RandomSample[{1, 2, 3} -> {a, b, c}]]
Out[8]= {a, b, c}

In[9]:= RandomSample[{a, b}, 5]
Out[9]= RandomSample[{a, b}, 5]

In[10]:= RandomSample[x]
Out[10]= RandomSample[x]
```

### Applications (4)

A random permutation: no element repeats

```mathematica
In[11]:= SeedRandom[1]; RandomSample[{a, b, c, d, e}]
Out[11]= {e, c, b, d, a}
```

Five distinct draws from 1..10

```mathematica
In[12]:= SeedRandom[1]; RandomSample[Range[10], 5]
Out[12]= {9, 8, 3, 1, 6}
```

UpTo clamps to the length

```mathematica
In[13]:= SeedRandom[1]; RandomSample[{a, b, c, d, e}, UpTo[3]]
Out[13]= {e, c, b}
```

Weighted, without replacement

```mathematica
In[14]:= SeedRandom[1]; RandomSample[{1, 1, 1, 10} -> {a, b, c, d}, 2]
Out[14]= {d, c}
```

## Options & behaviour

> **Packed arrays.** `RandomSample` and `RandomChoice` gather from a packed
> list or an `NDArray` directly, drawing from the **same generator sequence**
> the ordinary path uses — so `SeedRandom[n]` gives the same answer whether
> the argument is packed or not.

## Implementation notes

**Algorithm.** `builtin_randomsample` (in `src/random.c`) selects elements *without replacement*. The uniform form uses `fisher_yates_sample(total, n)`: a partial Fisher–Yates shuffle of an index array `[0..total)` that performs only the first `n` swaps (each swap picks `j` uniformly from the remaining suffix via `random_index`), returning the first `n` shuffled indices. With no count it returns a full random permutation. The size argument may be an integer or `UpTo[n]` (clamped to the list length via `is_upto`).

The weighted form `RandomSample[{w1,...}->{e1,...}, n]` uses `weighted_sample_without_replacement`, which repeatedly draws by inverse-CDF over the live weights (`u = U(0,1)*total`, linear scan accumulating cumulative weight) and zeroes the chosen weight so it cannot be picked again — i.e. sequential weighted sampling without replacement, O(n·count). Both paths draw from the shared Mersenne Twister state and deep-copy the selected elements into a `List[...]`. Requesting more than the available count (without `UpTo`) returns `NULL` (unevaluated).

- `Protected`.
- `RandomSample[{e1, e2, ...}, n]` never samples any of the ei more than once.
- `RandomSample[{e1, e2, ...}, n]` samples each of the ei with equal probability.
- `RandomSample[{e1, e2, ...}, UpTo[n]]` gives a sample of n of the ei, or as many as are available.
- RandomSample gives a different sequence of pseudorandom choices whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Requesting n greater than the list length (without UpTo) returns unevaluated.
- Uses the Fisher-Yates shuffle for uniform sampling without replacement.
- Weighted sampling removes selected elements and renormalizes weights.

**Attributes:** `Protected`.

## References

**See also:** [RandomChoice](../../random-number-generation/RandomChoice/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)
- Tests: [`tests/test_random.c`](https://github.com/stblake/mathilda/blob/main/tests/test_random.c)

## Notes & additional examples

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
