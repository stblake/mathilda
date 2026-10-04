# Lookup

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Lookup[assoc, key]`**

Gives the value for key, or Missing\["KeyAbsent", key\].

**`Lookup[assoc, key, default]`**

Uses default when key is absent.

**`Lookup[assoc, {k1, k2, ...}]`**

Looks up several keys at once (O(n+m)).

**`Lookup[{assoc1, assoc2, ...}, key, ...]`**

Threads over a list of associations and/or lists of rules.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Lookup[<|"a" -> 1, "b" -> 2|>, "b"]
Out[1]= 2

In[2]:= Lookup[<|"a" -> 1|>, "z", 0]
Out[2]= 0

In[3]:= Lookup[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3|>}, "a", 0]
Out[3]= {1, 3}

In[4]:= Lookup[<|"a" -> 1|>, "a", Print["never printed"]; 0]
Out[4]= 1

In[5]:= Lookup["a"][<|"a" -> 7|>]
Out[5]= 7
```

### Applications (4)

```mathematica
In[6]:= Lookup[<|a -> 1, b -> 2|>, a]
Out[6]= 1
```

A default for the absent key

```mathematica
In[7]:= Lookup[<|a -> 1, b -> 2|>, c, 0]
Out[7]= 0
```

Several keys at once

```mathematica
In[8]:= Lookup[<|a -> 1, b -> 2, c -> 3|>, {c, a}]
Out[8]= {3, 1}
```

No default: a Missing object

```mathematica
In[9]:= Lookup[<|a -> 1|>, c]
Out[9]= Missing["KeyAbsent", c]
```

## Implementation notes

**Algorithm.** `Lookup` is `HoldAll` so that the default stays lazy: `builtin_lookup`
evaluates only the association and the key(s), holding the default and evaluating it
once *per absent key*. `lookup_core` then dispatches on the key shape:

1. a single key → `assoc_lookup_value` (one hash probe);
2. `Key[k]` → the one literal key `k`, even when `k` is itself a list;
3. a `List` of keys → one `KeyIndex` is built over the association, then each key is
   an O(1) probe (O(n + m) overall);
4. a `List` of associations → the lookup threads, one `Lookup` per element.

A hit copies the stored value; a miss yields the evaluated default, or
`Missing["KeyAbsent", key]` when no default was given.

**Data structures.** The association's persistent `AssocIndex`, built lazily and
cached on the first single-key read (eager attachment does not survive the
fixed-point evaluation step, so it is deferred to the reader); the compiled
evaluator pre-builds it at its marshalling boundary so no worker thread mutates a
shared node. The multi-key path uses a transient open-addressing `KeyIndex`.

**Complexity / limits.** O(1) amortised per single key, O(n + m) for a list of `m`
keys. When the arguments do not form a valid lookup the call is returned with its
first two arguments evaluated.

- `HoldAll`, as in Mathematica: the association and key are evaluated, the
  default only when a key is actually absent (once per absent key), so a
  default with side effects runs only when needed.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Keys](../../data-structures/Keys/), [Values](../../data-structures/Values/), [HoldAll](../../expression-information/HoldAll/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

`Lookup` is `HoldAll`, so the default is evaluated only when the key is actually
absent (and once per absent key); `Lookup[<|a -> 1|>, a, Print["x"]]` prints
nothing. A list of keys shares a single index build, so it costs O(n + m) rather than
m separate scans. `Key[k]` looks up one literal key even when `k` is a list.
