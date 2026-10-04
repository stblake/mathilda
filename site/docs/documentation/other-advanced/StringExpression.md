# StringExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringExpression[p1, p2, ...] or p1 ~~ p2 ~~ ...`**

Represents a sequence of string patterns to be matched consecutively.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

~~ builds a StringExpression, flattened by Flat

```mathematica
In[1]:= FullForm["a" ~~ "b" ~~ "c"]
Out[1]= StringExpression["a", "b", "c"]
```

Letters then digits

```mathematica
In[2]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]
Out[2]= True
```

A two-character run

```mathematica
In[3]:= StringReplace["hello", "l" ~~ "l" -> "L"]
Out[3]= "heLo"
```

## Implementation notes

**Definition.** `StringExpression[p1, p2, ...]` represents a sequence of string
patterns to be matched consecutively — the long form of the `~~` operator
(`p1 ~~ p2 ~~ ...`). It is an **inert head**: `regex_init` gives it no builtin, only
attributes and a docstring. It is the glue of the string-pattern language, carrying a
run of literals and pattern atoms (`__`, `DigitCharacter`, `LetterCharacter`,
`NumberString`, `StartOfString`, ...) that the string engine then compiles and matches.

**Representation.** A bare `EXPR_FUNCTION` with head `StringExpression`. Its attributes
are `Flat | OneIdentity | Protected`: `Flat` means nested `~~` sequences collapse into
one flat argument list when evaluated (the parser builds `~~` right-associatively, so
`a ~~ b ~~ c` parses as `StringExpression[StringExpression[a, b], c]` and flattens to
`StringExpression[a, b, c]` on evaluation), and `OneIdentity` means a one-element
sequence is identified with its single element for pattern-matching purposes. The
string-pattern heads (`StringMatchQ`, `StringCases`, `StringReplace`, `StringSplit`,
...) consume a `StringExpression` by walking its arguments left to right through the
shared string-pattern engine.

**Usage & limits.** `Protected`. It is consumed, not evaluated to a value, so on its own
a `StringExpression[...]` stays as data; its meaning is supplied by whichever
string-pattern head receives it. A rule `lhs -> rhs` built inside a `~~` sequence is the
replacement spec used by `StringReplace`.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## References

- Source: [`src/strings/regex/regex_init.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regex_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringposition.c)

## Notes & additional examples

### Notes

`StringExpression[p1, p2, ...]`, written `p1 ~~ p2 ~~ ...`, represents a run of string
patterns matched consecutively. It is the backbone of the string-pattern language,
threading literals together with pattern atoms such as `__`, `DigitCharacter`,
`LetterCharacter` and `NumberString`.

It carries the attributes `Flat`, `OneIdentity` and `Protected`: `Flat` collapses
nested `~~` into one flat argument list on evaluation (so `"a" ~~ "b" ~~ "c"` becomes
`StringExpression["a", "b", "c"]`), and `OneIdentity` identifies a one-element sequence
with its single element for matching. It is consumed, not evaluated to a value — its
meaning is supplied by whichever string head (`StringMatchQ`, `StringCases`,
`StringReplace`, `StringSplit`, ...) receives it.
