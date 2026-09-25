# Part

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

expr\[\[i\]\] or Part\[expr, i\] gives the i-th part of expr. expr\[\[-i\]\] counts from the end. expr\[\[0\]\] gives the head of expr. expr\[\[i, j, ...\]\] or Part\[expr, i, j, ...\] is equivalent to expr\[\[i\]\]\[\[j\]\]..., descending into nested parts. expr\[\[{i1, i2, ...}\]\] gives a list of the parts i1, i2, ... of expr (wrapped in the head of expr). expr\[\[m;;n\]\] / expr\[\[m;;n;;s\]\] gives the span of parts m through n (with optional step s); ;; alone or All means all parts. Part is treated as atomic on Integer, Real, String, Symbol, Rational\[n, d\], and Complex\[re, im\]; Part\[atom, i\] for i != 0 stays unevaluated. Indices are 1-based and may be negative; out-of-range indices leave the expression unevaluated.

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>[[2 ;; 3]]
Out[1]= <|"b" -> 2, "c" -> 3|>

In[2]:= <|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>[[{-1, 1}]]
Out[2]= <|"d" -> 4, "a" -> 1|>

In[3]:= <|"a" -> 1, "b" -> 2, "c" -> 3, "d" -> 4|>[[{"c", "z"}]]
Out[3]= <|"c" -> 3, "z" -> Missing["KeyAbsent", "z"]|>

In[4]:= <|"x" -> {1, 2}, "y" -> {3, 4}|>[[All, 1]]
Out[4]= <|"x" -> 1, "y" -> 3|>

In[5]:= <|"a" -> 1|>[Key["a"]]
Out[5]= Missing["KeyAbsent", Key["a"]]
```

### Applications (9)

```mathematica
In[6]:= {a, b, c, d}[[2]]
Out[6]= b

In[7]:= {a, b, c, d}[[-1]]
Out[7]= d

In[8]:= {a, b, c, d}[[0]]
Out[8]= List

In[9]:= {{1, 2}, {3, 4}}[[2, 1]]
Out[9]= 3

In[10]:= {a, b, c, d}[[{1, 3}]]
Out[10]= {a, c}

In[11]:= {a, b, c, d, e, f}[[1 ;; 6 ;; 2]]
Out[11]= {a, c, e}

In[12]:= m[[All, 2]]
Out[12]= {2, 5, 8}

In[13]:= Tr[m[[{1, 3}, {1, 3}]]]
Out[13]= 10

In[14]:= (a + b + c)[[2]]
Out[14]= b
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Join two 2x10^6 | 0.899 s | 0.6 s | 0.397 s |
| RotateLeft 4x10^6 by 1000 | 0.855 s | 0.307 s | 0.456 s |

## Implementation notes

**Algorithm.** `builtin_part` extracts elements by index path, delegating to the recursive
`expr_part(expr, indices, nindices)`. Each index level may be: a positive or negative integer
(`-k` resolves to `len + k + 1`); `0`, which extracts the head and is allowed even on atoms; a
`List` of indices (extract several, returning a list); `All`; or a `Span` (`i;;j;;k`) built by
the parser from `;;` syntax, which `expr_part` resolves into an explicit element range with the
given start/end/step (negative endpoints wrap, `UpTo`/`All` endpoints clamp). Index paths apply
left to right, descending one structural level per index. Out-of-range or non-integer indices on
atoms yield `NULL` (unevaluated).

- An absent key gives `Missing["KeyAbsent", spec]` with the spec as written:
  `<|"a" -> 1|>[[Key["b"]]]` is `Missing["KeyAbsent", Key["b"]]`; inside a key
  list it becomes the entry `k -> Missing[...]`.
- As in Mathematica 15, the Part stays unevaluated for a position out of range,
  a span running off either end (an empty span just inside, such as `3 ;; 2`,
  gives `<||>`), positions mixed with keys, a bare symbol or real as the spec,
  or a deeper index that fails on any selected value. (Mathematica also prints
  `Part::partw`, `Part::take`, `Part::pmix`, `Part::pkspec1` or `Part::partd`;
  Mathilda does not print these messages.)
- Sub-associations share their entries with the source, so `RuleDelayed`
  values stay delayed.

**Attributes:** `NHoldRest`, `Protected`.

## References

**See also:** [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`expr[[i]]` (= `Part[expr, i]`) is 1-based; `expr[[-i]]` counts from the end and
`expr[[0]]` returns the head. Multi-index `expr[[i, j, ...]]` descends through
nested parts, a list of indices `expr[[{i1, i2, ...}]]` gathers several (rewrapped
in the original head), and spans `m ;; n ;; s` (with `All` or `;;` meaning "every
part") slice ranges. Part operates on any head, not just `List`, so it indexes
sums, products, and arbitrary symbolic structure uniformly. It is atomic on
Integer, Real, String, Symbol, Rational, and Complex; out-of-range or
unsupported indices leave the expression unevaluated rather than erroring.
