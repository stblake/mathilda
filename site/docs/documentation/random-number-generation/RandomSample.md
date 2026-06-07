# RandomSample

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RandomSample[{e1, e2, ...}, n]
    gives a pseudorandom sample of n of the ei, without replacement.
RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]
    gives a weighted pseudorandom sample of n of the ei.
RandomSample[{e1, e2, ...}]
    gives a pseudorandom permutation of the ei.
RandomSample[list, UpTo[n]]
    gives a sample of n of the ei, or as many as are available.
RandomSample never samples any element more than once.
Use SeedRandom to seed the pseudorandom generator for reproducible results.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]
Out[1]= {b, c, d}

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
```

## Implementation notes

- `Protected`.
- `RandomSample[{e1, e2, ...}, n]` never samples any of the ei more than once.
- `RandomSample[{e1, e2, ...}, n]` samples each of the ei with equal probability.
- `RandomSample[{e1, e2, ...}, UpTo[n]]` gives a sample of n of the ei, or as many as are available.
- RandomSample gives a different sequence of pseudorandom choices whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Requesting n greater than the list length (without UpTo) returns unevaluated.
- Uses the Fisher-Yates shuffle for uniform sampling without replacement.
- Weighted sampling removes selected elements and renormalizes weights.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
