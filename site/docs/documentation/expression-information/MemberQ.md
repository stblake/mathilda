# MemberQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MemberQ[list, form] returns True if an element of list matches form, and False otherwise.`**

**`MemberQ[list, form, levelspec] tests all parts of list specified by levelspec.`**

**`MemberQ[form] represents an operator form of MemberQ that can be applied to an expression.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= MemberQ[{1, 3, 4, 1, 2}, 2]
Out[1]= True

In[2]:= MemberQ[{x^2, y^2, x^3}, x^_]
Out[2]= True

In[3]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0, 2]
Out[3]= True

In[4]:= MemberQ[{{1, 1, 3, 0}, {2, 1, 2, 2}}, 0]
Out[4]= False
```

### Applications (5)

```mathematica
In[5]:= MemberQ[{1, 2, 3}, 2]
Out[5]= True

In[6]:= MemberQ[{1, 2, 3}, _Integer]
Out[6]= True

In[7]:= MemberQ[{x^2, y^3, z}, _^_]
Out[7]= True

In[8]:= MemberQ[{{1, 2}, {3, 4}}, 3, {2}]
Out[8]= True

In[9]:= MemberQ[#, 0] & /@ {{1, 2}, {0, 3}}
Out[9]= {False, True}
```

## Implementation notes

`builtin_memberq` (`src/patterns.c`) returns `True`/`False` for whether any part of the first argument matches the second (a pattern), tested by `do_member_at_level` which applies the pattern matcher `match` at each position within the level spec (default level 1 — immediate elements). It supports integer/`{min,max}`/`All`/`Infinity` level specs and the `Heads` option; the one-argument form returns an operator `Function[MemberQ[#1, patt]]`.

- `Protected`.
- Default option: `Heads -> False`.
- `form` can be a structural pattern.
- The first argument of `MemberQ` can have any head, not necessarily `List`.
- Returns immediately upon finding the first match.
- Standard level specifications are supported. The default value for `levelspec` in `MemberQ` is `{1}`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/)

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_array_flatten.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_flatten.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_complement.c`](https://github.com/stblake/mathilda/blob/main/tests/test_complement.c)

## Notes & additional examples

### Notes

`MemberQ[list, form]` tests whether any element of `list` matches `form`. The
second argument is a full pattern, not just a literal value — `MemberQ[list,
_Integer]` checks for any integer and `MemberQ[{x^2, y^3, z}, _^_]` finds the
first element stored as a `Power`. The optional level specification controls how
deep to look: `MemberQ[{{1, 2}, {3, 4}}, 3, {2}]` searches at level 2 (inside the
sub-lists) and so detects the `3`, whereas the default level 1 would not. Because
`MemberQ` is an ordinary predicate, it maps cleanly over collections, as in the
last example which flags which rows contain a zero.
