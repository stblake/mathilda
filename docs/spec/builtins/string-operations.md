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

