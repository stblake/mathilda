# Implies

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Implies[p, q]`**

The material implication p \\[Implies\] q, equivalent to !p || q. Implies\[False, q\] and Implies\[p, True\] are True, Implies\[True, q\] is q, and Implies\[p, False\] is !p.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Implies[True, q]
Out[1]= q

In[2]:= Implies[p, False]
Out[2]= Not[p]

In[3]:= LogicalExpand[Implies[p, q]]
Out[3]= Not[p] || q
```

### Applications (5)

A true premise reduces to the conclusion

```mathematica
In[4]:= Implies[True, q]
Out[4]= q
```

A false premise implies anything

```mathematica
In[5]:= Implies[False, q]
Out[5]= True
```

P implies False is the negation of p

```mathematica
In[6]:= Implies[p, False]
Out[6]= Not[p]
```

A statement implies itself

```mathematica
In[7]:= Implies[p, p]
Out[7]= True
```

The definition, !p || q

```mathematica
In[8]:= LogicalExpand[Implies[p, q]]
Out[8]= Not[p] || q
```

## Implementation notes

**Algorithm.** `Implies` is `Protected` (binary, no `Flat`/`Orderless`, since
material implication is neither associative nor commutative). `builtin_implies`
requires exactly two arguments and tests a short ladder of literal and structural
cases: `Implies[False, q]` and `Implies[p, True]` are `True`, `Implies[True, q]`
is `q`, `Implies[p, False]` is `Not[p]`, and `Implies[p, p]` (by `expr_eq`) is
`True`. Anything else returns `NULL` and stays symbolic.

**Data structures.** None; the two argument pointers are read in place and only
the chosen result is constructed.

**Complexity / limits.** `O(1)` — a fixed set of comparisons. No further
reasoning is done at this layer; `LogicalExpand` and `Reduce` rewrite
`Implies[p, q]` to its definition `!p || q` when a boolean normal form is wanted.

**Attributes:** `Protected`.

## References

**See also:** [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [Reduce](../../solutions-of-equations/Reduce/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Implies[p, q]` is material implication, `p ⟹ q`, logically `!p || q`. It takes
exactly two arguments and is neither associative nor commutative, so (unlike `And`
/`Or`/`Xor`/`Equivalent`) it is not `Flat` or `Orderless`.

It simplifies the literal and structural cases — `Implies[False, q]` and
`Implies[p, True]` are `True`, `Implies[True, q]` is `q`, `Implies[p, False]` is
`Not[p]`, and `Implies[p, p]` is `True` — and otherwise stays symbolic.
`LogicalExpand` and `Reduce` rewrite it to `!p || q` when a boolean normal form is
wanted.
