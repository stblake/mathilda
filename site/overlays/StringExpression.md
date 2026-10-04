### Worked examples

```mathematica
In[1]:= FullForm["a" ~~ "b" ~~ "c"]  (* ~~ builds a StringExpression, flattened by Flat *)
```

```mathematica
In[1]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]  (* letters then digits *)
```

```mathematica
In[1]:= StringReplace["hello", "l" ~~ "l" -> "L"]  (* a two-character run *)
```

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
