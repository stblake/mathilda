# String Operations

Strings are UTF-8, and every String function counts, indexes, slices and
reverses in **characters** (Unicode code points), as the reference does:
`StringLength["h\:00e9llo"]` is `5`, `StringReverse` keeps each character whole,
`StringPosition` reports character positions, and the string-pattern classes
are Unicode-aware (`LetterCharacter` matches `"\:00fc"`, `DigitCharacter`
matches any decimal digit, `IgnoreCase` folds `"\:00dc"`/`"\:00fc"`).
Pure-ASCII strings take a byte-speed path. A byte that does not begin a valid
UTF-8 sequence counts as one character, so malformed input never crashes or
splits.

`Style[expr, ...]` is presentation only. In OutputForm (`Print`,
`ToString`) it shows as `expr`, because there is no front end to style it.
`FullForm` and `InputForm` keep the wrapper.

String literals accept the reference's character escapes: `\[Name]` (the 1100+
named characters, e.g. `"\[Alpha]"`), `\:hhhh`, `\|hhhhhh`, `\.hh` and `\ooo`.
A backslash at the end of a line is a line continuation, and both characters
are dropped: `"ab\` followed by a newline and `cd"` is `"abcd"`. Leading
whitespace on the next line is kept. Before a CRLF only the CR is dropped, so
the LF stays. An escaped backslash `\\` before a newline gives a backslash
followed by a real newline.

Invalid arguments emit the reference's diagnostics through the message system
(`StringLength::string`, `StringPadLeft::intnm`, `StringRepeat::intp`,
`StringJoin::string`, `...::strse`, ...) and leave the call unevaluated.

## StringQ

- `StringQ[expr]`: `True` if `expr` is a string, `False` otherwise.
- Not `Listable`: a list is tested as one object. Map `StringQ` for an
  element-wise test.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringQ[""]
Out[1]= True

In[2]:= StringQ[{"a", "b"}]
Out[2]= False

In[3]:= StringQ /@ {"a", 1, x, "b"}
Out[3]= {True, False, False, True}
```

## StringLength

Gives the number of characters in a string.

- `StringLength["string"]`: Returns the number of characters (Unicode code
  points) as an integer.
- A non-string argument emits `StringLength::string` and stays unevaluated.
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

- `Characters["string"]`: Returns a `List` of single-character strings, one per
  Unicode code point.
- Each character is given as a length-1 string (single ASCII characters are
  shared, refcounted strings, so this allocates nothing per character).
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
- A non-string, non-list leaf emits `StringJoin::string` and leaves the call
  unevaluated (`StringJoin["a", 1]`).
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
- `StringPart["string", All]`, and spans containing `All`. Positions are
  characters.
- Out-of-range positions emit `StringPart::partw` (or `::take` for a span);
  an invalid spec emits `::pkspec1`. The call stays unevaluated.
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
- `StringTake["string", All]` / `None`, and a list of specs
  `StringTake["string", {{1, 2}, {4, 5}}]` (one result per spec).
- `StringTake[{s1, s2, ...}, spec]`: Gives the list of results for each si.
- `StringTake["string", 0]`: Returns the empty string `""`.
- Negative indices count from the end of the string.
- In `StringTake["string", {m, n, s}]`, m, n, and/or s can be negative.
- Positions are characters (code points). The range rules are the reference's:
  an empty result needs adjacent bounds (`{4, 3}`), a literal `0` start means
  "past the end".
- Out-of-range positions emit `StringTake::take`, a malformed spec
  `StringTake::mseqs`, a subject that is not a string or list of strings
  `::strse`; the call stays unevaluated. Over a list, an out-of-range element
  stays as an unevaluated `StringTake[si, spec]`.
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
is the complement of `StringTake`: the same standard reference sequence
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
- `StringDrop["string", All]` / `None`. Positions are characters, with
  `StringTake`'s range rules: `{4, 3}` (adjacent) drops nothing, while
  `{5, 4}` on `"abc"` or `{4, 2}` is out of range and emits `StringDrop::drop`
  (an invalid spec emits `::seqs`); the call stays unevaluated.
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

- `StringReverse["string"]`: Reverses the order of the characters in `"string"`
  (each multi-byte UTF-8 character stays whole).
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
- A position n is valid for `1 <= n <= StringLength["string"] + 1` (or the corresponding negative range), counting characters; an out-of-range position emits `StringInsert::ins`, a non-integer one `StringInsert::psl`, and the call stays unevaluated.
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
- `StringReplacePart["string", "snew", n]` means `{1, n}` (`-n`: `{-n, -1}`), and `{m, m-1}` inserts before character m.
- It does not thread over a list of strings (the reference emits `StringReplacePart::string`).
- `StringReplacePart[new, part][old]`: Operator form, equivalent to `StringReplacePart[old, new, part]`.
- Position specifications use the form returned by `StringPosition`: each is a pair `{m, n}` of first/last character positions. Negative positions count from the end. All positions refer to `"string"` *before* any replacement is done.
- Positions are not allowed to overlap: a range that overlaps a previously accepted one emits `StringReplacePart::ovlp` and its new string is not inserted.
- `StringReplacePart["string", "", ...]` deletes the selected characters.
- Positions are characters. Out-of-range positions emit `::repart`, a malformed one `::spos` / `::psl`, a new-string/position length mismatch `::plen`; the call stays unevaluated.
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
PCRE2 (the same engine the reference uses).

- `RegularExpression["regex"]`: an inert head that evaluates to itself and
  carries the pattern.
- Supported syntax: `.` `[c1c2]` `[c1-c2]` `[^...]`, quantifiers `p*` `p+` `p?`
  `p{m,n}` and their non-greedy forms `*?` `+?` `??`, groups `(...)`, and
  alternation `p1|p2`; classes `\d \D \s \S \w \W` and `[[:name:]]` (alnum,
  alpha, ascii, blank, cntrl, digit, graph, lower, print, punct, space, upper,
  word, xdigit); anchors `^ $ \b \B`; and inline options `(?i)` `(?m)` `(?s)`.
- In a replacement right-hand side, `$n` stands for the n-th captured group and
  `$0` for the whole match; `$$` is a literal `$`. The expansion applies to every
  string in the right-hand side, before it is evaluated: `:> {"$1", "$2"}`,
  `:> f["$1", {"$2"}]` and `:> ToUpperCase["$1"]` all see the captures.
- Requires PCRE2 at build time; without it these builtins warn and stay
  unevaluated.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}

In[3]:= StringReplace["John Smith", RegularExpression["(\\w+) (\\w+)"] -> "$2 $1"]
Out[3]= "Smith John"

In[4]:= StringCases["caf\:00e9 na\:00efve", RegularExpression["\\w+"]]
Out[4]= {"café", "naïve"}
```

When the matcher's backtracking limit is reached (`(a+)+$` on a long run of `a`s ending
in `b`), `RegularExpression::limit` says so and the call stays unevaluated; a
`StringContainsQ` no longer answers `False` for an undecided match.

## StringMatchQ

Tests whether a whole string matches a pattern.

- `StringMatchQ["string", patt]`: `True` if all of `"string"` matches `patt`
  (the pattern is anchored to the whole string), else `False`.
- `StringMatchQ[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `patt` may be any string pattern (see [String patterns](#string-patterns)),
  `RegularExpression["re"]`, a literal string, or a list of alternatives. A
  non-string subject leaves the call unevaluated.
- **Wildcards**: in `StringMatchQ` (only), `*` in a literal string matches any
  run of characters and `@` one or more characters that are not uppercase
  letters; `\\*` and `\\@` are literal. `StringMatchQ["abc", "a*"]` is `True`.
- **Options**: `IgnoreCase -> True` folds case.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}

In[3]:= StringMatchQ["apple.txt", "*.txt"]
Out[3]= True

In[4]:= StringMatchQ["aXb", x_ ~~ "X" ~~ x_]
Out[4]= False
```

`SpellingCorrection -> False` is accepted (and in `Options[StringMatchQ]`); no
spelling correction is ever applied.

## StringContainsQ

Tests whether *any substring* of a string matches a pattern. Where
`StringMatchQ` anchors to the whole string, `StringContainsQ` searches inside
it.

- `StringContainsQ["string", patt]`: `True` if some substring matches `patt`,
  else `False`.
- `StringContainsQ["string", {p1, p2, ...}]`: `True` if some substring matches
  any of the `pi` — equivalent to `StringContainsQ["string", p1 | p2 | ...]`.
  An empty list matches nothing.
- `StringContainsQ[{s1, s2, ...}, patt]`: gives the list of results for each
  `si`. A non-string subject — or a subject list holding one — leaves the call
  unevaluated rather than passing the element through.
- `StringContainsQ[patt]`: the operator form, applicable to a string or a list
  of strings. A curried `IgnoreCase` option is carried through.
- Equivalent to `!StringFreeQ["string", patt]` and to
  `StringMatchQ["string", ___ ~~ patt ~~ ___]`.
- **Options** (`Options[StringContainsQ]` → `{IgnoreCase -> False}`;
  `SetOptions` changes the default): `IgnoreCase -> True` folds case. An
  explicit option in the call overrides the registered default.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringContainsQ["bcde", "b" ~~ __ ~~ "e"]
Out[1]= True

In[2]:= StringContainsQ[{"a", "b", "ab", "abcd", "bcde"}, "a"]
Out[2]= {True, False, True, True, False}

In[3]:= StringContainsQ["bac 123", RegularExpression["a.*"] ~~ DigitCharacter ..]
Out[3]= True

In[4]:= StringContainsQ["abcd", "BC", IgnoreCase -> True]
Out[4]= True

In[5]:= Select[{"abc", "xyz", "bat"}, StringContainsQ["a"]]
Out[5]= {"abc", "bat"}
```

## StringFreeQ

Tests whether *no substring* of a string matches a pattern — the negation of
`StringContainsQ`, sharing its implementation.

- `StringFreeQ["string", patt]`: `True` if no substring matches `patt`, else
  `False`.
- `StringFreeQ["string", {p1, p2, ...}]`: `True` if no substring matches any of
  the `pi`.
- `StringFreeQ[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `StringFreeQ[patt]`: the operator form.
- **Options**: `IgnoreCase -> True` folds case (default `False`).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringFreeQ["abcd", "a"]
Out[1]= False

In[2]:= StringFreeQ["abcade", x_ ~~ x_]
Out[2]= True

In[3]:= StringFreeQ[{"ability", "listable", "argument"}, "a" ~~ __ ~~ "t" ~~ ___]
Out[3]= {False, True, False}

In[4]:= StringFreeQ["ac", IgnoreCase -> True]["BACCD"]
Out[4]= False
```

## StringStartsQ

Tests whether a *prefix* of a string matches a pattern.

- `StringStartsQ["string", patt]`: `True` if `"string"` begins with a substring
  matching `patt`.
- `StringStartsQ["string", {p1, p2, ...}]`: `True` if it begins with a match of
  any of the `pi`.
- `StringStartsQ[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `StringStartsQ[patt]`: the operator form.
- Equivalent to `StringContainsQ["string", StartOfString ~~ patt]`.
- **Options**: `IgnoreCase -> True` folds case (default `False`).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringStartsQ["commit", "co"]
Out[1]= True

In[2]:= StringStartsQ["commit", "om"]
Out[2]= False

In[3]:= StringStartsQ[{"apple", "banana"}, "a"]
Out[3]= {True, False}

In[4]:= StringStartsQ["a123", LetterCharacter ~~ DigitCharacter ..]
Out[4]= True
```

## StringEndsQ

Tests whether a *suffix* of a string matches a pattern.

- `StringEndsQ["string", patt]`: `True` if `"string"` ends with a substring
  matching `patt`.
- `StringEndsQ["string", {p1, p2, ...}]`: `True` if it ends with a match of any
  of the `pi`.
- `StringEndsQ[{s1, s2, ...}, patt]`: gives the list of results for each `si`.
- `StringEndsQ[patt]`: the operator form.
- Equivalent to `StringContainsQ["string", patt ~~ EndOfString]`.
- **Options**: `IgnoreCase -> True` folds case (default `False`).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringEndsQ["commit", "it"]
Out[1]= True

In[2]:= StringEndsQ["commit", "mi"]
Out[2]= False

In[3]:= StringEndsQ[{"apple", "banana"}, "a"]
Out[3]= {False, True}

In[4]:= StringEndsQ["a123", DigitCharacter ..]
Out[4]= True
```

## StringCases

Extracts the substrings of a string that match a pattern.

- `StringCases["string", patt]`: the list of non-overlapping substrings that
  match `patt`, from left to right.
- `StringCases["string", patt -> rhs]` / `patt :> rhs`: `rhs` for each match,
  with every pattern name in `patt` bound to the text it matched:
  `StringCases["abc123", x__ ~~ "123" :> x]` is `{"abc"}`, and `rhs` may be any
  expression (`x:DigitCharacter.. :> ToExpression[x]` gives integers). A name in
  an alternative that did not take part in the match is bound to `""`.
- For a `RegularExpression` pattern every string in `rhs` is a template: `$n` is
  the n-th capture group, `$0` the whole match, `$$` a literal `$` --
  `RegularExpression["([a-z])=(\\d+)"] :> {"$1", "$2"}` gives
  `{{"x", "1"}, {"y", "22"}}` on `"x=1, y=22"`. With any other pattern `"$1"` is
  literal text, as in the reference.
- `StringCases[{s1, s2, ...}, patt]`: gives the list of results for each `si`;
  non-string elements pass through unchanged.
- `patt` may be `RegularExpression["re"]`, a literal string, or a list of
  patterns/rules (leftmost match wins, ties broken by order).
- **Options** (`Options[StringCases]` → `{IgnoreCase -> False, Overlaps ->
  False}`; `SetOptions` changes the defaults): `Overlaps -> False` (default)
  does not treat overlapping substrings as separate; `Overlaps -> True` counts
  overlaps but keeps only the first matching substring at each start;
  `Overlaps -> All` gives every matching substring at every start.
  `IgnoreCase -> True` folds case. An explicit option in the call overrides the
  registered default.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringCases["a13b12c17a32", RegularExpression["[^a1]"]]
Out[1]= {"3", "b", "2", "c", "7", "3", "2"}

In[2]:= StringCases["AaBBccDDeefG", RegularExpression["[[:upper:]]+"]]
Out[2]= {"A", "BB", "DD", "G"}

In[3]:= StringCases["AAAA", "AA", Overlaps -> True]
Out[3]= {"AA", "AA", "AA"}

In[4]:= StringCases["AAA", x__, Overlaps -> All]
Out[4]= {"AAA", "AA", "A", "AA", "A", "A"}

In[5]:= StringCases["aAbB", "a", IgnoreCase -> True]
Out[5]= {"a", "A"}

In[6]:= StringCases["abc123", (x__ ~~ "123") :> x]
Out[6]= {"abc"}

In[7]:= StringCases["a1b22c333", x:DigitCharacter.. /; StringLength[x] > 1 :> ToExpression[x]]
Out[7]= {22, 333}
```

## StringCount

Counts the substrings of a string that match a pattern — the counting-only
companion of `StringCases`, running the same match enumeration but without
building the matched substrings.

- `StringCount["string", "sub"]`: the number of times `"sub"` appears as a
  substring of `"string"`.
- `StringCount["string", patt]`: the number of substrings of `"string"` that
  match the general string expression `patt`.
- `StringCount["string", {p1, p2, ...}]`: counts occurrences of any of the `pi`.
- `StringCount[{s1, s2, ...}, patt]`: gives the list of results for each `si`;
  non-string elements pass through unchanged.
- `patt` may also be a `Rule` or `RuleDelayed` — only the left-hand side matters,
  since a replacement cannot change a count.
- **Options** (`Options[StringCount]` → `{IgnoreCase -> False, Overlaps ->
  False}`; `SetOptions` changes the defaults): identical in meaning to the
  `StringCases` options above.
- Equal to `Length[StringCases[...]]` for every pattern and option setting: both
  builtins, and `StringPosition`, share one scanner (`regex_scan` in
  `src/strings/regex/regex_common.c`).
- A non-string, non-list subject leaves the call unevaluated; an arity other
  than two positional arguments emits `StringCount::argrx`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringCount["the cat sat on the mat", "at"]
Out[1]= 3

In[2]:= StringCount["AAAA", "AA"]
Out[2]= 2

In[3]:= StringCount["AAAA", "AA", Overlaps -> True]
Out[3]= 3

In[4]:= StringCount["AAAA", x__, Overlaps -> All]
Out[4]= 10

In[5]:= StringCount[{"a1", "b22", "ccc"}, DigitCharacter]
Out[5]= {1, 2, 0}

In[6]:= StringCount["aAbB", "a", IgnoreCase -> True]
Out[6]= 2

In[7]:= StringCount["x=1.5, y=-2", NumberString]
Out[7]= 2
```

## StringReplace

Replaces matches of a pattern in a string.

- `StringReplace["string", patt -> rep]`: replaces each non-overlapping match
  of `patt` by `rep`, with the pattern names of `patt` bound to the text they
  matched (as in `StringCases`); unmatched text is copied. For a
  `RegularExpression` pattern a string `rep` may use `$n`/`$0`.
- A replacement that is not a string cannot be spliced into text, so the result
  becomes a `StringExpression` of the text runs and those values, as in
  the reference: `StringReplace["abc123", x:DigitCharacter :> ToExpression[x]]`
  is `"abc"~~1~~2~~3`. `Sequence[...]` replacements splice; `Sequence[]` deletes.
- `StringReplace["string", {patt1 -> rep1, ...}]`: applies a list of rules; at
  each position the leftmost match wins, ties broken by rule order.
- `StringReplace[{s1, s2, ...}, rules]`: gives the list of results for each
  `si`.
- Zero-width matches (e.g. `\b`) insert `rep` at each boundary without dropping
  characters.
- `StringReplace["string", {}]` returns the string unchanged.
- **Options**: `IgnoreCase -> False`, `MetaCharacters -> None`.
  - `MetaCharacters -> None` (the default) treats a literal string pattern
    literally.
  - `MetaCharacters -> Automatic` reads it as a regular expression, so
    `StringReplace["a.b", "." -> "x", MetaCharacters -> Automatic]` is `"xxx"`.
  - A custom set (a 3-character string or a list of three strings) is not
    supported and stays unevaluated.
  - Any other value emits `StringReplace::meta`.
- A backtracking regular expression on a long subject (e.g. `(a|b)*` over a
  few kB) would overflow the PCRE2 JIT's fixed stack. The match is rerun in
  the interpreter rather than failing.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"

In[3]:= StringReplace["abc123", a__ ~~ "123" :> a <> "!"]
Out[3]= "abc!"

In[4]:= StringReplace["abc123", x:DigitCharacter :> ToExpression[x]]
Out[4]= "abc"~~1~~2~~3

In[5]:= StringReplace["aaaa", "a" -> "b", 2]
Out[5]= "bbaa"
```

## StringSplit

Splits a string at matches of a delimiter pattern.

- `StringSplit["string"]`: splits at runs of whitespace (equivalent to
  `StringSplit["string", Whitespace]`), dropping leading/trailing empties.
- `StringSplit["string", patt]`: the list of substrings between non-overlapping
  matches of the delimiter `patt`.
- `StringSplit["string", {p1, p2, ...}]`: splits at any of the `pi`.
- `StringSplit["string", patt -> val]`: inserts `val` at the position of each
  delimiter (`patt :> val` evaluates `val` per match, with every pattern name
  bound to the text it matched, e.g. `"-" ~~ x:DigitCharacter :> x`). Rejoining
  the pieces reproduces `StringReplace`.
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

`StringSplit[s, "", n]` gives at most n pieces, the last taking the rest.

## StringTrim

Trims substrings matching a pattern from the beginning and end of a string.

- `StringTrim["string"]`: trims runs of whitespace (spaces, tabs, newlines)
  from both ends.
- `StringTrim["string", patt]`: removes one substring matching the string
  pattern `patt` from the start and one from the end, as in the reference:
  `StringTrim["xxabxx", "x"]` is `"xabx"`; use `"x"..` to remove runs. The end
  match must begin at or after the end of the start match.
- Option `IgnoreCase -> True` matches case-insensitively.
- `StringTrim[{s1, s2, ...}, ...]`: gives the list of results for each `si`;
  non-string list elements pass through unchanged.
- `patt` accepts the full shared string-pattern vocabulary (see below),
  `RegularExpression["re"]`, and literal strings. A pattern that anchors itself
  to the start (e.g. `RegularExpression["^ *"]`) therefore only trims the front.
- A non-string subject leaves the call unevaluated; an argument count other than
  1 or 2 emits `StringTrim::argt`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringTrim["   aaa bbb ccc   "]
Out[1]= "aaa bbb ccc"

In[2]:= StringTrim["++++aaa bbb ccc----", ("+" | "-") ..]
Out[2]= "aaa bbb ccc"

In[3]:= StringTrim["   aaa bbb ccc   ", RegularExpression["^ *"]]
Out[3]= "aaa bbb ccc   "

In[4]:= StringTrim["007bond007", DigitCharacter ..]
Out[4]= "bond"
```

## StringExtract

Splits a string into blocks and selects blocks by position. Built on top of
`StringSplit`, so the split behaviour (whitespace-run collapsing, empty-end
trimming, the full string-pattern vocabulary) is identical.

- `StringExtract["string", spec]`: splits at runs of whitespace and selects with
  `spec`, one of `n` (nth block), `-n` (nth from the end), `{n1, n2, ...}` (a
  collection), `n1;;n2` (a span), or `All` (every block).
- `StringExtract["string", sep -> pos]`: delimits blocks with the string pattern
  `sep`. `sep -> All` is exactly `StringSplit["string", sep]`.
- `StringExtract["string", pos1, pos2, ...]`: extracts across levels, delimiting
  with whitespace at the lowest level, a single `"\n"` at the next, `"\n\n"` at
  the next, and so on for higher levels.
- `StringExtract["string", sep1 -> pos1, sep2 -> pos2, ...]`: uses `sepi` as the
  separator for successive levels, giving a nested list of blocks.
- Out-of-range blocks yield `Missing["PartAbsent", pos]`. `{pos...}` is
  equivalent to `Part[StringSplit["string"], {pos...}]` for in-range positions.
- A list of strings threads. A non-string, non-list subject leaves the call
  unevaluated; fewer than 2 arguments emits `StringExtract::argm`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringExtract["a bbb  cccc aa   d", 2]
Out[1]= "bbb"

In[2]:= StringExtract["a bbb  cccc aa   d", 2 ;; 4]
Out[2]= {"bbb", "cccc", "aa"}

In[3]:= StringExtract["a--bbb--ccc--dddd", "--" -> 3]
Out[3]= "ccc"

In[4]:= StringExtract["a 1\nb 2\nc 3 x", All, 3]
Out[4]= {Missing["PartAbsent", 3], Missing["PartAbsent", 3], "x"}
```

An integer span is the list of its positions: past the last block each is
`Missing["PartAbsent", k]` (`StringExtract["a b c", 2;;5]`).

## StringPosition

Gives the `{start, end}` character positions at which substrings matching a
string pattern occur, in the 1-based inclusive form used by `StringTake`,
`StringDrop`, and `StringReplacePart`. Positions count characters, and the
matching is the same PCRE2 pattern engine as `StringCases`.

- `StringPosition["string", patt]`: positions of every substring matching the
  string pattern `patt` (a literal, `RegularExpression[...]`, a general string
  expression such as `x_ ~~ x_`, or a character class).
- `StringPosition["string", patt, n]`: only the first `n` occurrences (with
  `Overlaps -> False` the scan stops once it has them).
- `StringPosition[patt]`: operator form.
- `StringPosition["string", {p1, p2, ...}]`: matches of each `pi`, merged and
  ordered by start position (a list of patterns is not the same as
  `Alternatives`).
- `StringPosition[{s1, s2, ...}, patt]`: threads, giving one result list per
  subject; non-string elements pass through unchanged.
- **Options** (`Options[StringPosition]` → `{IgnoreCase -> False, Overlaps ->
  True}`; `SetOptions` changes the defaults): `Overlaps -> True` (default)
  includes overlapping substrings but only the first match starting at each
  position; `Overlaps -> False` excludes overlaps (greedy, global); `Overlaps ->
  All` includes every matching substring at every start. `IgnoreCase -> True`
  folds case. An explicit option in the call overrides the registered default.
- A non-string, non-list subject leaves the call unevaluated; zero arguments
  emits `StringPosition::argb`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringPosition["abXYZaaabXYZaaaaXYZXYZ", "XYZ"]
Out[1]= {{3, 5}, {10, 12}, {17, 19}, {20, 22}}

In[2]:= StringPosition["AABBBAABABBCCCBAAA", x_ ~~ x_]
Out[2]= {{1, 2}, {3, 4}, {4, 5}, {6, 7}, {10, 11}, {12, 13}, {13, 14}, {16, 17}, {17, 18}}

In[3]:= StringPosition["AAAAA", "AA", Overlaps -> False]
Out[3]= {{1, 2}, {3, 4}}

In[4]:= StringPosition["abAB", "a", IgnoreCase -> True]
Out[4]= {{1, 1}, {3, 3}}
```

An occurrence count of 0 gives `{}`.

## StringPartition

Partitions a string into fixed-length substrings. Lengths count characters.

- `StringPartition["string", n]`: non-overlapping substrings of length `n`
  (offset `d = n`). Trailing characters that cannot fill a full block are
  dropped.
- `StringPartition["string", n, d]`: length-`n` substrings starting every `d`
  characters. All substrings are length `n`; if `d > n` characters in the middle
  are skipped, and if `d < n` the substrings overlap.
- `StringPartition["string", UpTo[n]]`: substrings of length up to `n` (offset
  `d = n`), allowing the final substring to be shorter so every character
  appears. Combinable with an explicit offset as `UpTo[n], d`.
- With an explicit offset `d < n` and `UpTo[n]`, the trailing blocks shorten as
  the reference's do (`StringPartition["abcdefg", UpTo[3], 2]` is
  `{"abc", "cde", "efg", "g"}`).
- Does not thread over a list (the reference emits `StringPartition::string`).
  A bad length or offset emits `::intpm` / `::intnm` / `::intp`; other arities
  emit `::argtu` / `::argt`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringPartition["123456789123456789", 9]
Out[1]= {"123456789", "123456789"}

In[2]:= StringPartition["123456789", 2, 1]
Out[2]= {"12", "23", "34", "45", "56", "67", "78", "89"}

In[3]:= StringPartition["123456789", UpTo[2]]
Out[3]= {"12", "34", "56", "78", "9"}

In[4]:= StringPartition["ababababab", 3]
Out[4]= {"aba", "bab", "aba"}
```

## StringRepeat

Builds a string by repeating a base string. Lengths count characters.

- `StringRepeat["string", n]`: `"string"` concatenated `n` times.
- `StringRepeat["string", n, max]`: up to `n` copies of `"string"`, truncated so
  the total length is at most `max` characters (a partial final copy is allowed).
- `n == 0`, an empty base string, or `max == 0` gives `""`. A non-string
  subject emits `StringRepeat::string`, a non-integer or negative `n`/`max`
  `StringRepeat::intp` (both leave the call unevaluated); other arities emit
  `StringRepeat::argt` (2 or 3 arguments expected).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringRepeat["a", 50]
Out[1]= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

In[2]:= StringRepeat["abc", 10]
Out[2]= "abcabcabcabcabcabcabcabcabcabc"

In[3]:= StringRepeat["ab", 10, 19]
Out[3]= "abababababababababa"
```

## StringPadLeft / StringPadRight

Make a string a given length by padding it (with spaces or a repeated pad
string) or truncating it. Lengths count characters.
`StringPadLeft` pads on the left and truncates keeping the last `n` characters;
`StringPadRight` pads on the right and truncates keeping the first `n`.

- `StringPad*["string", n]`: length `n`, padded with spaces or truncated.
- `StringPad*["string", n, "padding"]`: pads with `"padding"` repeated
  cyclically and aligned to the result, as in the reference: on the right,
  result character `i` is `padding[[i mod len]]`; on the left the cycle is
  anchored at the result's end. `StringPadRight["ab", 6, "xyz"]` is
  `"abzxyz"`, `StringPadLeft["ab", 7, "xyz"]` is `"zxyzxab"`.
- `StringPad*["string"]`: a single string needs a length -- emits
  `::strlist` and stays unevaluated (the reference).
- `StringPad*[{s1, s2, ...}]`: pads each string to the length of the longest,
  so all become the same length.
- `StringPad*[{s1, s2, ...}, n, ...]`: pads or truncates each string to length `n`.
- Invalid arguments emit the reference's messages and leave the call
  unevaluated: a subject that is not a string or a list of strings
  (`::strlist`), a non-integer or negative `n` (`::intnm`), an empty or
  non-string `padding` (`::stringnz`). `StringPadLeft[{}, n]` is `""`. Other
  arities emit `::argb` (1 to 3 arguments expected).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringPadLeft["abcde", 10]
Out[1]= "     abcde"

In[2]:= StringPadRight["abcde", 10, "."]
Out[2]= "abcde....."

In[3]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}]
Out[3]= {"    a", "   ab", "  abc", " abcd", "abcde"}

In[4]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}, 3]
Out[4]= {"  a", " ab", "abc", "bcd", "cde"}
```

## StringRiffle

Assembles a string from a (possibly nested) list by inserting separators between
elements — the inverse of `StringSplit`. Non-string leaves are converted with
`ToString` (OutputForm, so `27` becomes `"27"` and `f["a"]` becomes `"f[a]"`);
string leaves are used verbatim. The result is assembled in one buffer.

- `StringRiffle[list]`: joins with the default scheme — a single space at the
  innermost level, one extra newline per level above it. A 2-D array joins rows
  with `"\n"` and cells with `" "`; a 3-D array joins blocks with `"\n\n"`, rows
  with `"\n"`, cells with `" "`.
- `StringRiffle[list, sep]`: inserts the string `sep` between the top-level
  elements; deeper levels fall back to the default scheme.
- `StringRiffle[list, {"left", "sep", "right"}]`: a 3-string list is a delimiter
  triple — joins with `sep` and wraps the result in `left`…`right`.
- `StringRiffle[list, sep1, sep2, ...]`: inserts separator `sep_i` (a string or a
  delimiter triple) between elements at level `i` (level 1 is outermost); levels
  deeper than the supplied separators use the default scheme.
- An empty list gives `""` (or just the delimiters for a triple). A first
  argument that is neither a list nor a string emits `StringRiffle::list`, and a
  separator that is neither a string nor a 3-string list `StringRiffle::string`;
  both leave the call unevaluated. No arguments emits
  `StringRiffle::argm` (1 or more arguments expected).
- Any nesting depth is riffled (the traversal is iterative): `StringRiffle[Nest[List, "a", 100000]]` is `"a"` (audit StringRiffle#1).
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringRiffle[{"a", "b", "c", "d", "e"}]
Out[1]= "a b c d e"

In[2]:= StringRiffle[{"a", "b", "c", "d", "e"}, ", "]
Out[2]= "a, b, c, d, e"

In[3]:= StringRiffle[{"a", "b", "c", "d", "e"}, {"(", " ", ")"}]
Out[3]= "(a b c d e)"

In[4]:= StringRiffle[{{"a", 27}, {"b", 28}, {"c", 29}}, {"{", ", ", "}"}, ": "]
Out[4]= "{a: 27, b: 28, c: 29}"
```

## String patterns

Beyond `RegularExpression["re"]` and literal strings, every pattern-taking
string function (`StringMatchQ`, `StringCases`, `StringCount`,
`StringPosition`, `StringReplace`, `StringSplit`, `StringTrim`, and the
substring predicates `StringContainsQ` and friends) shares one translator
(`string_pattern.c`) that turns symbolic reference string patterns into PCRE2:

| Pattern | Matches |
|---------|---------|
| `Whitespace` | a run of whitespace |
| `WhitespaceCharacter` | one whitespace character |
| `LetterCharacter` / `DigitCharacter` / `WordCharacter` | one letter / digit / letter-or-digit |
| `HexadecimalCharacter` / `PunctuationCharacter` | one hex digit / punctuation character |
| `NumberString` | a signed integer or decimal number (`"+5"`, `"-.5"`, `"1."`) |
| `CharacterRange["a", "z"]` | one character in the range |
| `StartOfString` / `EndOfString` | zero-width; the start / end of the whole string |
| `StartOfLine` / `EndOfLine` / `WordBoundary` | zero-width; a line start / line end / word boundary |
| `p1 ~~ p2 ~~ ...` (`StringExpression`) | `p1`, then `p2`, ... in sequence |
| `p1 \| p2` (`Alternatives`) | any of the alternatives |
| `p ..` / `p ...` (`Repeated` / `RepeatedNull`) | one-or-more / zero-or-more of `p` |
| `Repeated[p, n]`, `Repeated[p, {n}]`, `Repeated[p, {m, n}]` | up to `n` / exactly `n` / `m` to `n` repetitions (`RepeatedNull` likewise) |
| `Shortest[p]` / `Longest[p]` | `p` with the shortest / longest repetitions |
| `Except[p]` / `Except[p, q]` | one character / a `q`-match that does not begin a `p`-match |
| `_` / `__` / `___` | one character / one or more / zero or more (newlines included) |
| `x : p`, `x_`, `x__` (`Pattern`) | `p`, naming the matched text `x` |
| `p ? f` (`PatternTest`) | `p`, if `f[text]` is `True` (`_?LetterQ`-style tests are a character class) |
| `p /; cond` (`Condition`) | `p`, if `cond` is `True` with the names bound so far |

Pattern names follow the reference: a name's first occurrence captures, and
every later occurrence must match the same text -- `x_ ~~ x_` matches two equal
characters, and in a repeated pattern `(x:LetterCharacter)..` every repetition
equals the first (`"aaa"`, not `"abc"`). A failed `Condition` or `PatternTest`
backtracks, so `StringCases["aabc", (x_ ~~ y_) /; x != y]` finds `"ab"`. The
blanks match newlines; `.` inside a `RegularExpression` keeps PCRE's meaning.

An unsupported pattern leaves the call unevaluated. `~~` (`StringExpression`)
is a `Flat` operator with precedence just below `Alternatives`; `x:p` binds
tighter than `~~` (`x__ ~~ y:DigitCharacter..` is
`StringExpression[x__, Pattern[y, Repeated[DigitCharacter]]]`).

`StringExpression` joins adjacent literal strings (`"a" ~~ "b"` is `"ab"`) and a lone
string is itself; matching is unchanged.

## StringDelete

- `StringDelete["string", patt]`: deletes every substring matching `patt`
  (any string pattern, or a list of them) -- `StringReplace[s, patt -> ""]`.
- `StringDelete[{s1, s2, ...}, patt]`: one result per string.
- `StringDelete[patt]`: operator form.
- **Options**: `IgnoreCase`. **Attributes**: `Protected`.

```mathematica
In[1]:= StringDelete["a1b22c", DigitCharacter..]
Out[1]= "abc"

In[2]:= StringDelete["Hello World", "l", IgnoreCase -> True]
Out[2]= "Heo Word"

In[3]:= StringDelete[{"x1", "y22"}, DigitCharacter]
Out[3]= {"x", "y"}
```

## StringReplaceList

- `StringReplaceList["string", rules]`: the list of strings obtained by
  replacing a single match at a time -- one result for every match start of
  every rule (overlapping), in order of position, then rule.
- `StringReplaceList["string", rules, n]`: at most `n` results; `n = Infinity`
  is no limit (it gave `{}`, audit StringReplaceList#1). `StringReplace[s, rules, 0]`
  replaces nothing.
- **Options**: `IgnoreCase`. **Attributes**: `Protected`.

```mathematica
In[1]:= StringReplaceList["aaa", "a" -> "b"]
Out[1]= {"baa", "aba", "aab"}

In[2]:= StringReplaceList["aaa", "aa" -> "X"]
Out[2]= {"Xa", "aX"}

In[3]:= StringReplaceList["abc", {"a" -> "1", "b" -> "2"}]
Out[3]= {"1bc", "a2c"}
```

## StringTakeDrop

- `StringTakeDrop["string", spec]`: `{StringTake["string", spec],
  StringDrop["string", spec]}` for any `StringTake` spec.
- `StringTakeDrop[{s1, s2, ...}, spec]`: one pair per string.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringTakeDrop["abcdef", {2, 4}]
Out[1]= {"bcd", "aef"}

In[2]:= StringTakeDrop["abcdef", 2]
Out[2]= {"ab", "cdef"}

In[3]:= StringTakeDrop[{"abcd", "wxyz"}, -1]
Out[3]= {{"d", "abc"}, {"z", "wxy"}}
```

## StringRotateLeft / StringRotateRight

- `StringRotateLeft["string", n]`: cycles the characters `n` places to the
  left (`n` may be negative or exceed the length); `StringRotateLeft["string"]`
  rotates by one. `StringRotateRight` rotates the other way.
- **Attributes**: `Protected`.

```mathematica
In[1]:= StringRotateLeft["abcdef", 2]
Out[1]= "cdefab"

In[2]:= StringRotateRight["abcdef", 2]
Out[2]= "efabcd"

In[3]:= StringRotateLeft["abcdef", -8]
Out[3]= "efabcd"
```

## ToUpperCase / ToLowerCase

- `ToUpperCase["string"]`, `ToLowerCase["string"]`: full Unicode case mapping
  (`ToUpperCase["\:00df"]` is `"SS"`, `ToLowerCase["\[CapitalSigma]"]` is
  `"\[Sigma]"`). Non-string arguments stay unevaluated.
- **Attributes**: `Listable`, `Protected`.

```mathematica
In[1]:= ToUpperCase["straße"]
Out[1]= "STRASSE"

In[2]:= ToLowerCase[{"ABC", "ΣΟΦΙΑ"}]
Out[2]= {"abc", "σοφια"}

In[3]:= ToUpperCase[x]
Out[3]= ToUpperCase[x]
```

## Capitalize / Decapitalize

- `Capitalize["string"]`: the first character in uppercase;
  `Decapitalize["string"]` in lowercase. A list of strings maps; anything else
  emits `::strse`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= Capitalize["hello world"]
Out[1]= "Hello world"

In[2]:= Capitalize["élan"]
Out[2]= "Élan"

In[3]:= Decapitalize[{"Alpha", "Beta"}]
Out[3]= {"alpha", "beta"}
```

`Capitalize[s, "AllWords"]` capitalizes the first letter of every word;
`"TitleCase"` all but short function words not at either end
(`"The Lord of the Rings"`). `Decapitalize[s, "AllWords"]` likewise.

## LetterQ / DigitQ / UpperCaseQ / LowerCaseQ

- `LetterQ["string"]`: every character is a letter of any script.
- `DigitQ["string"]`: every character is an ASCII digit `0`-`9` (as in
  the reference, even though the pattern `DigitCharacter` matches any Unicode
  decimal digit).
- `UpperCaseQ["string"]` / `LowerCaseQ["string"]`: every character is an
  uppercase / lowercase letter.
- All are `True` for `""`. A non-string emits `::string`.
- In a string pattern, `_?LetterQ`, `_?DigitQ`, `_?UpperCaseQ` and
  `_?LowerCaseQ` compile to a character class that accepts exactly what the
  predicate does: `StringMatchQ["\:0661", _?DigitQ]` is `False` like
  `DigitQ["\:0661"]` (it used to be `True`), while
  `StringMatchQ["\:0661", DigitCharacter]` stays `True`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= LetterQ["αβγ"]
Out[1]= True

In[2]:= DigitQ["\:0662"]
Out[2]= False

In[3]:= UpperCaseQ["ABC"]
Out[3]= True

In[4]:= LowerCaseQ["abc1"]
Out[4]= False
```

## ToCharacterCode / FromCharacterCode

- `ToCharacterCode["string"]`: the Unicode code points, as a packed integer
  array (above the packing threshold). A list of strings gives a list of lists.
- Code point 0 is an ordinary character: `FromCharacterCode[{65, 0, 66}]` has
  `StringLength` 3 and round-trips through `ToCharacterCode` (held internally
  as the bytes C0 80, its UTF-8 byte is 0; audit FromCharacterCode#1).
- `FromCharacterCode[n]`, `FromCharacterCode[{n1, ...}]`,
  `FromCharacterCode[{{...}, ...}]`: the inverse. It reads a packed or `NDArray`
  integer buffer directly. A code outside `0`..`1114111` emits
  `FromCharacterCode::notunicode` (`::intnm` for a negative or non-integer).
- `ToCharacterCode[s, "enc"]` / `FromCharacterCode[codes, "enc"]`: the same in
  an encoding. `"UTF8"` (or `"UTF-8"`) is the byte sequence of the UTF-8 form
  (a packed integer array); decoding reports a byte that starts no valid
  sequence as `$CharacterEncoding::utf8` and keeps it as the code point of the
  same value. `"ISOLatin1"` and `"ASCII"` are single-byte: `ToCharacterCode`
  gives `None` for a character the encoding cannot represent. `"Unicode"` is the
  default. Any other encoding name leaves the call unevaluated.
- No `Compile[]` lowering: one side of each is a `String`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= ToCharacterCode["h\:00e9"]
Out[1]= {104, 233}

In[2]:= FromCharacterCode[{72, 105}]
Out[2]= "Hi"

In[3]:= FromCharacterCode[{{72, 105}, {955, 8364}}]
Out[3]= {"Hi", "λ€"}

In[4]:= ToCharacterCode["a\[Alpha]", "UTF8"]
Out[4]= {97, 206, 177}

In[5]:= ToCharacterCode["a\[Alpha]", "ASCII"]
Out[5]= {97, None}
```

## CharacterRange / Alphabet

- `CharacterRange["c1", "c2"]`: the characters with codes from `c1`'s to
  `c2`'s; `CharacterRange[n1, n2]` by code. `{}` when the range is empty.
- `Alphabet[]`: `{"a", ..., "z"}`.
- **Attributes**: `Protected`.

```mathematica
In[1]:= CharacterRange["a", "e"]
Out[1]= {"a", "b", "c", "d", "e"}

In[2]:= CharacterRange[65, 70]
Out[2]= {"A", "B", "C", "D", "E", "F"}

In[3]:= Length[Alphabet[]]
Out[3]= 26
```

`Alphabet["English"]`, `Alphabet["Greek"]` (24 letters, no final sigma) and
`Alphabet["Russian"]` (33) give the lowercase letters.

## Operator forms and occurrence limits

- Operator forms: `StringCases[patt][s]`, `StringReplace[rules][s]`,
  `StringMatchQ[patt][s]`, `StringPosition[patt][s]`, `StringDelete[patt][s]`
  (in addition to the `StringContainsQ` family). `StringSplit`, `StringCount`
  and `StringTrim` have none, as in the reference.
- Occurrence limits: `StringCases[s, patt, n]`, `StringReplace[s, rules, n]`,
  `StringPosition[s, patt, n]`, `StringReplaceList[s, rules, n]` (`n` an
  integer >= 0 or `Infinity`). With non-overlapping matching the scan stops as
  soon as it has `n` matches.
- `IgnoreCase` is an option of every pattern function (`Options[StringReplace]`
  etc.; `SetOptions` changes the default).

## Performance notes

- Compiled patterns are cached (64 entries, keyed by the pattern expression and
  flags), so applying one pattern to many strings -- `Map[StringMatchQ[#, p] &,
  list]`, `Select[list, StringContainsQ[p]]` -- compiles it once.
  `RegularExpression["re"]` validation shares the cache.
- A rule set is scanned with a leftmost-match cursor that remembers each rule's
  next match, so `StringReplace`/`StringSplit`/`StringCases` with several rules
  cost O(n · rules) rather than re-searching every rule at every position.
- `StringTrim` matches its front `PCRE2_ANCHORED` and its back
  `PCRE2_ENDANCHORED`; `StringJoin`/`StringRiffle` assemble the result in one
  buffer; results are handed to the String constructor without a second copy;
  single ASCII characters are shared, refcounted strings.

