### Worked examples

```mathematica
In[1]:= StringSplit["a b c"]  (* default: split at runs of whitespace *)
```

```mathematica
In[1]:= StringSplit["a,b,c", ","]  (* a literal delimiter *)
```

```mathematica
In[1]:= StringSplit["a1b2c", DigitCharacter]  (* a character-class delimiter *)
```

### Notes

`StringSplit` returns the substrings between non-overlapping matches of the
delimiter, which is translated to PCRE by the shared string-pattern engine
(literals, `RegularExpression`, character classes, `~~`, `|`, `..`, `Except`, …).

Zero-length substrings between two adjacent interior delimiters are kept; leading
and trailing empties are dropped unless `All` is given as a third argument. The
empty-string delimiter `""` splits at every character; `IgnoreCase -> True` folds
case.
