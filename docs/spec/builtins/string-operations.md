# String Operations

## StringLength

Gives the number of characters in a string.

- `StringLength["string"]`: Returns the number of characters as an integer.
- Returns unevaluated for non-string arguments.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= StringLength["tiger"]
Out[1]= 5

In[2]:= StringLength[""]
Out[2]= 0

In[3]:= StringLength["hello world"]
Out[3]= 11

In[4]:= StringLength[{\"ABC\", \"DE\", \"F\"}]
Out[4]= {3, 2, 1}

In[5]:= StringLength[x]
Out[5]= StringLength[x]
```

## Characters

Gives a list of the characters in a string.

- `Characters["string"]`: Returns a `List` of single-character strings.
- Each character is given as a length-1 string.
- Returns unevaluated for non-string arguments.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= Characters["ABC"]
Out[1]= {"A", "B", "C"}

In[2]:= Characters["A string."]
Out[2]= {"A", " ", "s", "t", "r", "i", "n", "g", "."}

In[3]:= Characters[""]
Out[3]= {}

In[4]:= Characters[{"ABC", "DEF", "XYZ"}]
Out[4]= {{"A", "B", "C"}, {"D", "E", "F"}, {"X", "Y", "Z"}}

In[5]:= Characters[x]
Out[5]= Characters[x]
```

## StringJoin

Concatenates strings together.

- `StringJoin["s1", "s2", ...]`: Joins all string arguments into a single string.
- `StringJoin[{"s1", "s2", ...}]`: Flattens all lists recursively and joins enclosed strings.
- `StringJoin[]`: Returns the empty string `""`.
- The infix form is `"s1" <> "s2" <> ...`.
- Returns unevaluated if any non-string, non-list leaf is encountered.
- **Attributes**: `Flat`, `OneIdentity`, `Protected`.

```mathematica
In[1]:= StringJoin["abcd", "ABCD", "xyz"]
Out[1]= "abcdABCDxyz"

In[2]:= "abcd" <> "ABCD" <> "xyz"
Out[2]= "abcdABCDxyz"

In[3]:= StringJoin[{{"AB", "CD"}, "XY"}]
Out[3]= "ABCDXY"

In[4]:= StringJoin[]
Out[4]= ""

In[5]:= StringJoin["a", x]
Out[5]= StringJoin["a", x]

In[6]:= StringJoin[Characters["hello"]]
Out[6]= "hello"
```

## StringPart

Extracts characters from a string by position.

- `StringPart["string", n]`: Gives the nth character in `"string"` as a single-character string.
- `StringPart["string", -n]`: Counts from the end to give the nth-to-last character.
- `StringPart["string", {n1, n2, ...}]`: Gives a list of the specified characters.
- `StringPart["string", m;;n]`: Gives a list of characters from position m through n.
- `StringPart["string", m;;n;;s]`: Gives characters from m through n in steps of s.
- `StringPart[{s1, s2, ...}, spec]`: Gives the list of results for each of the si.
- Negative indices count from the end of the string.
- In `StringPart["string", m;;n;;s]`, m, n, and/or s can be negative.
- Returns unevaluated for invalid arguments or out-of-bounds indices.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringPart["abcdefghijklm", 6]
Out[1]= "f"

In[2]:= StringPart["abcdefghijklm", {1, 3, 5}]
Out[2]= {"a", "c", "e"}

In[3]:= StringPart["abcdefghijklm", -4]
Out[3]= "j"

In[4]:= StringPart["abcdefghijklm", 1;;6]
Out[4]= {"a", "b", "c", "d", "e", "f"}

In[5]:= StringPart["abcdefghijklm", 1;;-1;;2]
Out[5]= {"a", "c", "e", "g", "i", "k", "m"}

In[6]:= StringPart["abcdefghijklm", -1;;1;;-2]
Out[6]= {"m", "k", "i", "g", "e", "c", "a"}

In[7]:= StringPart[{"abcd", "efgh", "ijklm"}, 1]
Out[7]= {"a", "e", "i"}

In[8]:= StringPart[{"abcd", "efgh", "ijklm"}, {1, -1}]
Out[8]= {{"a", "d"}, {"e", "h"}, {"i", "m"}}

In[9]:= StringPart["abcdef", -3;;-1]
Out[9]= {"d", "e", "f"}

In[10]:= StringPart["abcde", 5;;1;;-1]
Out[10]= {"e", "d", "c", "b", "a"}

In[11]:= StringPart[x, 1]
Out[11]= StringPart[x, 1]
```

## StringTake

Gives a substring of a string.

- `StringTake["string", n]`: Gives a string containing the first n characters.
- `StringTake["string", -n]`: Gives the last n characters.
- `StringTake["string", {n}]`: Gives the nth character.
- `StringTake["string", {m, n}]`: Gives characters m through n.
- `StringTake["string", {m, n, s}]`: Gives characters m through n in steps of s.
- `StringTake["string", UpTo[n]]`: Gives n characters, or as many as are available.
- `StringTake[{s1, s2, ...}, spec]`: Gives the list of results for each si.
- `StringTake["string", 0]`: Returns the empty string `""`.
- Negative indices count from the end of the string.
- In `StringTake["string", {m, n, s}]`, m, n, and/or s can be negative.
- Returns unevaluated for invalid arguments or out-of-bounds counts.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringTake["abcdefghijklm", 6]
Out[1]= "abcdef"

In[2]:= StringTake["abcdefghijklm", -4]
Out[2]= "jklm"

In[3]:= StringTake["abcdefghijklm", {5, 10}]
Out[3]= "efghij"

In[4]:= StringTake["abcdefghijklm", {6}]
Out[4]= "f"

In[5]:= StringTake["abcdefghijklm", {1, -1, 2}]
Out[5]= "acegikm"

In[6]:= StringTake[{"abcdef", "stuv", "xyzw"}, -2]
Out[6]= {"ef", "uv", "zw"}

In[7]:= StringTake["abc", UpTo[4]]
Out[7]= "abc"

In[8]:= StringTake["abcdef", {-3, -1}]
Out[8]= "def"

In[9]:= StringTake["abcde", {5, 1, -1}]
Out[9]= "edcba"

In[10]:= StringTake["abc", 0]
Out[10]= ""

In[11]:= StringTake[x, 1]
Out[11]= StringTake[x, 1]
```

## StringDrop

Gives a string with a specified subset of its characters removed. `StringDrop`
is the complement of `StringTake`: the same standard Wolfram Language sequence
specification selects the characters to *remove*, and the remaining characters
are concatenated in order.

- `StringDrop["string", n]`: Gives `"string"` with its first n characters dropped.
- `StringDrop["string", -n]`: Gives `"string"` with its last n characters dropped.
- `StringDrop["string", {n}]`: Gives `"string"` with its nth character dropped.
- `StringDrop["string", {m, n}]`: Gives `"string"` with characters m through n dropped.
- `StringDrop["string", {m, n, s}]`: Drops characters m through n in steps of s.
- `StringDrop["string", UpTo[n]]`: Drops n characters, or as many as are available.
- `StringDrop[{s1, s2, ...}, spec]`: Gives the list of results for each of the si.
- In `StringDrop["string", {m, n, s}]`, m, n, and/or s can be negative.
- A decreasing range (m > n) is empty and drops nothing.
- A call whose argument count is not two emits `StringDrop::argrx` and is left unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringDrop["abcdefghijklm", 4]
Out[1]= "efghijklm"

In[2]:= StringDrop["abcdefghijklm", -4]
Out[2]= "abcdefghi"

In[3]:= StringDrop["abcdefghijklm", {5, 10}]
Out[3]= "abcdklm"

In[4]:= StringDrop["abcdefghijklm", {3}]
Out[4]= "abdefghijklm"

In[5]:= StringDrop["abcdefghijklm", {1, -1, 2}]
Out[5]= "bdfhjl"

In[6]:= StringDrop[{"abcdef", "xyzw", "stuv"}, -2]
Out[6]= {"abcd", "xy", "st"}

In[7]:= StringDrop["abc", UpTo[4]]
Out[7]= ""

In[8]:= StringDrop["abcdefghijklm", {5, -4}]
Out[8]= "abcdklm"

In[9]:= StringDrop[]
StringDrop::argrx: StringDrop called with 0 arguments; 2 arguments are expected.
Out[9]= StringDrop[]
```

## StringReverse

Reverses the order of the characters in a string.

- `StringReverse["string"]`: Reverses the order of the characters in `"string"`.
- `StringReverse[{s1, s2, ...}]`: Gives the list of results for each of the si (via `Listable`).
- Returns unevaluated for non-string arguments.
- A call with any number of arguments other than one emits `StringReverse::argx` and is left unevaluated.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= StringReverse["abcdef"]
Out[1]= "fedcba"

In[2]:= StringReverse[{"cat", "dog", "fish", "coelenterate"}]
Out[2]= {"tac", "god", "hsif", "etaretneleoc"}

In[3]:= StringReverse[""]
Out[3]= ""

In[4]:= StringReverse[x]
Out[4]= StringReverse[x]

In[5]:= StringReverse[]
StringReverse::argx: StringReverse called with 0 arguments; 1 argument is expected.
Out[5]= StringReverse[]
```

## StringInsert

Inserts one string into another at one or more positions.

- `StringInsert["string", "snew", n]`: Makes the first character of `"snew"` the nth character of the result (i.e. inserts before the original nth character).
- `StringInsert["string", "snew", -n]`: Makes the last character of `"snew"` the nth character from the end of the result.
- `StringInsert["string", "snew", {n1, n2, ...}]`: Inserts a copy of `"snew"` at each of the positions ni.
- `StringInsert[{s1, s2, ...}, "snew", spec]`: Gives the list of results for each of the si.
- The ni all refer to positions in `"string"` *before* any insertion is done. Negative positions count from the end.
- A position n is valid for `1 <= n <= StringLength["string"] + 1` (or the corresponding negative range); an out-of-range or non-integer position leaves the call unevaluated.
- A call whose argument count is not three emits `StringInsert::argrx` and is left unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringInsert["abcdefghijklm", "XYZ", 4]
Out[1]= "abcXYZdefghijklm"

In[2]:= StringInsert["abcdefghijklm", "XYZ", -4]
Out[2]= "abcdefghijXYZklm"

In[3]:= StringInsert["abcdefghijklm", "XYZ", {2, 3, 7}]
Out[3]= "aXYZbXYZcdefXYZghijklm"

In[4]:= StringInsert["1234567890123456", ".", Range[4, 16, 3]]
Out[4]= "123.456.789.012.345.6"

In[5]:= StringInsert["1234567890123456", ".", Range[-16, -4, 3]]
Out[5]= "1.234.567.890.123.456"

In[6]:= StringInsert[{"abc", "de"}, "X", 2]
Out[6]= {"aXbc", "dXe"}

In[7]:= StringInsert[]
StringInsert::argrx: StringInsert called with 0 arguments; 3 arguments are expected.
Out[7]= StringInsert[]
```

## StringReplacePart

Replaces one or more ranges of characters in a string with new strings.

- `StringReplacePart["string", "snew", {m, n}]`: Replaces the characters at positions m through n by `"snew"`.
- `StringReplacePart["string", "snew", {{m1, n1}, {m2, n2}, ...}]`: Inserts a copy of `"snew"` at each of the given ranges.
- `StringReplacePart["string", {"snew1", "snew2", ...}, {{m1, n1}, ...}]`: Replaces the characters at each range by the corresponding new string. The list of new strings must be the same length as the list of positions.
- `StringReplacePart[{s1, s2, ...}, snew, part]`: Gives the list of results for each of the si.
- `StringReplacePart[new, part][old]`: Operator form, equivalent to `StringReplacePart[old, new, part]`.
- Position specifications use the form returned by `StringPosition`: each is a pair `{m, n}` of first/last character positions. Negative positions count from the end. All positions refer to `"string"` *before* any replacement is done.
- Positions are not allowed to overlap: a range that overlaps a previously accepted one emits `StringReplacePart::ovlp` and its new string is not inserted.
- `StringReplacePart["string", "", ...]` deletes the selected characters.
- A malformed or out-of-range position, or a new-string/position length mismatch, leaves the call unevaluated.
- A call whose argument count is not two (operator form) or three emits `StringReplacePart::argt` and is left unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {2, 5}]
Out[1]= "aABCDEFGHfghijk"

In[2]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {{1, 1}, {3, 5}, {-3, -1}}]
Out[2]= "ABCDEFGHbABCDEFGHfghABCDEFGH"

In[3]:= StringReplacePart["abcdefghijk", "ABCDEFGH", {-3, -2}]
Out[3]= "abcdefghABCDEFGHk"

In[4]:= StringReplacePart["abcdefghijk", {"XYZ", "ABCD"}, {{2, 3}, {-2, -2}}]
Out[4]= "aXYZdefghiABCDk"

In[5]:= StringReplacePart["ABCDEFGH", {2, 5}]["abcdefghijk"]
Out[5]= "aABCDEFGHfghijk"

In[6]:= StringReplacePart["abcde", "", {2, 4}]
Out[6]= "ae"

In[7]:= StringReplacePart["abcde", "XYZ", {{1, 3}, {3, 5}}]
StringReplacePart::ovlp: Position {3,5} overlaps previous positions; new string XYZ will not be inserted.
Out[7]= "XYZde"

In[8]:= StringReplacePart[]
StringReplacePart::argt: StringReplacePart called with 0 arguments; 2 or 3 arguments are expected.
Out[8]= StringReplacePart[]
```

## RegularExpression

Represents a class of strings given by a PCRE regular expression, for use in
`StringMatchQ`, `StringCases`, `StringReplace`, and `StringSplit`. Backed by
PCRE2 (the same engine the Wolfram Language uses).

- `RegularExpression["regex"]`: an inert head that evaluates to itself and
  carries the pattern.
- Supported syntax: `.` `[c1c2]` `[c1-c2]` `[^...]`, quantifiers `p*` `p+` `p?`
  `p{m,n}` and their non-greedy forms `*?` `+?` `??`, groups `(...)`, and
  alternation `p1|p2`; classes `\d \D \s \S \w \W` and `[[:name:]]` (alnum,
  alpha, ascii, blank, cntrl, digit, graph, lower, print, punct, space, upper,
  word, xdigit); anchors `^ $ \b \B`; and inline options `(?i)` `(?m)` `(?s)`.
- In a replacement right-hand side, `$n` stands for the n-th captured group and
  `$0` for the whole match; `$$` is a literal `$`.
- Requires PCRE2 at build time; without it these builtins warn and stay
  unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}
```

## StringMatchQ

Tests whether a whole string matches a pattern.

- `StringMatchQ["string", patt]`: `True` if all of `"string"` matches `patt`
  (the pattern is anchored to the whole string), else `False`.
- `StringMatchQ[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `patt` may be `RegularExpression["re"]`, a literal string (exact match), or a
  list of alternatives. A non-string subject leaves the call unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}
```

## StringCases

Extracts the substrings of a string that match a pattern.

- `StringCases["string", patt]`: the list of non-overlapping substrings that
  match `patt`, from left to right.
- `StringCases["string", patt -> rhs]`: the `rhs` for each match, with `$n`
  replaced by the n-th captured group and `$0` by the whole match.
- `StringCases[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `patt` may be `RegularExpression["re"]`, a literal string, or a list of
  patterns/rules (leftmost match wins, ties broken by order).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringCases["a13b12c17a32", RegularExpression["[^a1]"]]
Out[1]= {"3", "b", "2", "c", "7", "3", "2"}

In[2]:= StringCases["AaBBccDDeefG", RegularExpression["[[:upper:]]+"]]
Out[2]= {"A", "BB", "DD", "G"}
```

## StringReplace

Replaces matches of a pattern in a string.

- `StringReplace["string", patt -> rep]`: replaces each non-overlapping match
  of `patt` by `rep`, with `$n`/`$0` expanded; unmatched text is copied.
- `StringReplace["string", {patt1 -> rep1, ...}]`: applies a list of rules; at
  each position the leftmost match wins, ties broken by rule order.
- `StringReplace[{s1, s2, ...}, rules]`: gives the list of results for each
  `si`.
- Zero-width matches (e.g. `\b`) insert `rep` at each boundary without dropping
  characters.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"
```

## StringSplit

Splits a string at matches of a delimiter pattern.

- `StringSplit["string"]`: splits at runs of whitespace (equivalent to
  `StringSplit["string", Whitespace]`), dropping leading/trailing empties.
- `StringSplit["string", patt]`: the list of substrings between non-overlapping
  matches of the delimiter `patt`.
- `StringSplit["string", {p1, p2, ...}]`: splits at any of the `pi`.
- `StringSplit["string", patt -> val]`: inserts `val` at the position of each
  delimiter (`patt :> val` evaluates `val` per match; a named `x:patt :> f[x]`
  binds the matched text). Rejoining the pieces reproduces `StringReplace`.
- `StringSplit["string", patt, n]`: splits into at most `n` substrings.
- `StringSplit["string", patt, All]`: keeps the leading/trailing empty
  substrings that are otherwise dropped.
- `StringSplit[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- The empty-string delimiter `""` splits at every character. Zero-length
  substrings between two adjacent interior delimiters are kept.
- Option `IgnoreCase -> True` matches delimiters case-insensitively.
- `patt` accepts the full shared string-pattern vocabulary (see below), plus
  `RegularExpression["re"]` and literal strings. Zero-width delimiters (e.g.
  `(?m)^`) split at positions.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringSplit["a bbb  cccc aa   d"]
Out[1]= {"a", "bbb", "cccc", "aa", "d"}

In[2]:= StringSplit["a-b:c-d:e-f-g", {":", "-"}]
Out[2]= {"a", "b", "c", "d", "e", "f", "g"}

In[3]:= StringSplit["a b::c d::e f g", "::" -> "--"]
Out[3]= {"a b", "--", "c d", "--", "e f g"}

In[4]:= StringSplit["This is a sentence, which goes on.", Except[WordCharacter] ..]
Out[4]= {"This", "is", "a", "sentence", "which", "goes", "on"}
```

## String patterns

Beyond `RegularExpression["re"]` and literal strings, `StringMatchQ`,
`StringCases`, `StringReplace`, and `StringSplit` share a translator
(`string_pattern.c`) that turns symbolic Wolfram string patterns into PCRE:

| Pattern | Matches |
|---------|---------|
| `Whitespace` | a run of whitespace (`\s+`) |
| `WhitespaceCharacter` | one whitespace character |
| `LetterCharacter` / `DigitCharacter` / `WordCharacter` | one letter / digit / letter-or-digit |
| `NumberString` | a signed integer or decimal number |
| `p1 ~~ p2 ~~ ...` (`StringExpression`) | `p1`, then `p2`, ... in sequence |
| `p1 \| p2` (`Alternatives`) | any of the alternatives |
| `p ..` / `p ...` (`Repeated` / `RepeatedNull`) | one-or-more / zero-or-more of `p` |
| `Except[p]` | one character that does not begin a `p`-match |
| `_` (`Blank[]`) | any single character |
| `x : p` (`Pattern`) | `p`, capturing the matched text as `x` |
| `_?LetterQ` (and `DigitQ`/`UpperCaseQ`/`LowerCaseQ`) | one character satisfying the predicate |

An unsupported pattern leaves the call unevaluated. `~~` (`StringExpression`)
is a `Flat` operator with precedence just below `Alternatives`.

