# AssociationMap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AssociationMap[f, {k1, k2, ...}]`**

Gives \<|k1 -\> f\[k1\], k2 -\> f\[k2\], ...|\>.

**`AssociationMap[f, assoc]`**

Applies f to each rule k -\> v of assoc; the results (rules, lists of rules or associations) form the new association.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= AssociationMap[#^2 &, {1, 2, 3, 4}]
Out[1]= <|1 -> 1, 2 -> 4, 3 -> 9, 4 -> 16|>

In[2]:= AssociationMap[Reverse, <|"a" -> 1, "b" -> 2|>]
Out[2]= <|1 -> "a", 2 -> "b"|>
```

### Applications (3)

```mathematica
In[3]:= AssociationMap[f, {a, b, c}]
Out[3]= <|a -> f[a], b -> f[b], c -> f[c]|>
```

Keys are the inputs, values the results

```mathematica
In[4]:= AssociationMap[#^2 &, {1, 2, 3}]
Out[4]= <|1 -> 1, 2 -> 4, 3 -> 9|>
```

Over an association, f acts on each key -> value rule

```mathematica
In[5]:= AssociationMap[Reverse, <|a -> 1, b -> 2|>]
Out[5]= <|1 -> a, 2 -> b|>
```

## Implementation notes

**Algorithm.** `AssociationMap` has two forms. The key form
`AssociationMap[f, {k1, ...}]` (`builtin_associationmap`, `assoc.c`) builds
`<|k1 -> f[k1], ...|>`, each `f[k]` left for the evaluator. The association form
`AssociationMap[f, assoc]` (`ops_associationmap`, `assoc_ops.c`) applies `f` to
each entry *as a rule* `k -> v` and splices the result: `f` may return a rule, a
list of rules, an association, or `Nothing` (which contributes no entry). Both
forms canonicalise the collected rules through `assoc_from_rules`.

**Data structures.** The rule array feeds `assoc_from_rules`, whose transient
`KeyIndex` open-addressing hash set de-duplicates keys (first position, last
value). The association form keeps the raw applications in an unevaluated
`Association` if `f` ever returns something that is not a rule / rule-list /
association.

**Complexity / limits.** `O(n)` applications of `f` plus `O(n)` hashing. The key
form requires a `List` of keys; an invalid application in the association form
raises `AssociationMap::invrlf` through the message funnel.

**Attributes:** `Protected`.

## References

**See also:** [Nothing](../../lists-and-iteration/Nothing/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`AssociationMap[f, {k1, ...}]` builds `<|k1 -> f[k1], ...|>` — the keys are the
list elements and the values are `f` applied to them, which is the natural way to
tabulate a function over a set of inputs. Given an association instead of a key
list, `f` is applied to each entry *as a rule* `k -> v`, and its result (a rule,
a list of rules, an association, or `Nothing`) is spliced into the output — so
`AssociationMap[Reverse, assoc]` swaps keys and values.
