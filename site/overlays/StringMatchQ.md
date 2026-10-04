### Worked examples

```mathematica
In[1]:= StringMatchQ["hello", "h" ~~ ___]  (* ~~ builds a string expression; ___ is zero or more characters *)
```

```mathematica
In[1]:= StringMatchQ["abc123", LetterCharacter .. ~~ DigitCharacter ..]  (* letters followed by digits *)
```

```mathematica
In[1]:= StringMatchQ[{"cat", "dog"}, "c" ~~ __]  (* threads over the subjects *)
```

### Notes

`StringMatchQ` tests whether the *whole* string matches: the pattern is compiled
with a `\A(?:...)\z` wrap so both ends are anchored. The pattern may be a literal
string, `RegularExpression[...]`, a general string expression, or a list of
alternatives (a match if any one matches).

A list of subjects threads, giving a list of `True`/`False`; a non-string subject
leaves the call unevaluated.
