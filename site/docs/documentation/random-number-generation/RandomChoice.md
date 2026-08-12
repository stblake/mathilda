# RandomChoice

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomChoice[{e1, e2, ...}]`**

gives a pseudorandom choice of one of the ei.

**`RandomChoice[list, n]`**

gives a list of n pseudorandom choices.

**`RandomChoice[list, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of pseudorandom choices.

**`RandomChoice[{w1, w2, ...} -> {e1, e2, ...}]`**

gives a pseudorandom choice weighted by the wi.

**`RandomChoice[wlist -> elist, n]`**

gives a list of n weighted choices.

**`RandomChoice[wlist -> elist, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of weighted choices.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= SeedRandom[42]; RandomChoice[{a, b, c, d, e}]
Out[1]= e

In[2]:= SeedRandom[42]; RandomChoice[{a, b, c}, 5]
Out[2]= {c, a, c, c, c}

In[3]:= SeedRandom[42]; Dimensions[RandomChoice[{1, 2, 3}, {3, 4}]]
Out[3]= {3, 4}

In[4]:= RandomChoice[{1, 0, 0} -> {a, b, c}]
Out[4]= a

In[5]:= RandomChoice[{1, 0} -> {x, y}, 5]
Out[5]= {x, x, x, x, x}

In[6]:= RandomChoice[x]
Out[6]= RandomChoice[x]
```

## Implementation notes

**Algorithm.** `builtin_randomchoice` (in `src/random.c`) selects elements *with replacement*. For the uniform form `RandomChoice[{e1,...,ek}]`, `random_index(k)` draws an index with `mpz_urandomm` over the shared Mersenne Twister state and the chosen element is deep-copied. For the weighted form `RandomChoice[{w1,...}->{e1,...}]`, the handler builds a cumulative-weight array and `weighted_random_index` performs inverse-CDF sampling: it draws `u = U(0,1) * total` and **binary-searches** for the first index whose cumulative weight exceeds `u` (O(log k) per draw). The `RandomChoice[spec, n]` and `RandomChoice[spec, {n1,...}]` forms recurse over the dimension spec via `random_choice_array` / `weighted_choice_array`, drawing one element per leaf. Selection is always with replacement; for sampling without replacement see `RandomSample`.

- `Protected`.
- `RandomChoice[{e1, e2, ...}]` chooses with equal probability between all of the ei.
- RandomChoice gives a different sequence of pseudorandom choices whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Weighted selection uses cumulative weight binary search for efficient O(log n) per choice.
- Weights must be non-negative real numbers with a positive total.

**Attributes:** `Protected`.

## References

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_random.c`](https://github.com/stblake/mathilda/blob/main/tests/test_random.c)
