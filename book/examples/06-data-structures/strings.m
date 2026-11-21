# 11.4 Strings: atoms, and the basic operations
s = "Hello, world";
{AtomQ[s], Head[s], StringLength[s]}
StringTake[s, 5]
StringTake[s, -5]
"Hello" <> ", " <> "there"
StringSplit["the quick  brown fox"]
StringSplit["a,b,,c", ","]
StringReplace[s, "world" -> "there"]
ToUpperCase[s]
Characters["stressed"]
StringReverse["stressed"]
StringRiffle[{"red", "green", "blue"}, ", "]
ToString[x^2 + 1]
ToExpression["x^2 + 1"]
ToExpression["Sum[k, {k, 1, 100}]"]
# Canonical order of strings: punctuation is ignored at first, lower case before upper (v0.530)
Sort[{"b", "a-b", "B", "ab", ".hid", "a.m", "A"}]
