---
source: src/strings/regex/regex_init.c
---
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
