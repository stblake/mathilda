### Worked examples

```mathematica
In[1]:= StringReplace["hello world", "o" -> "0"]  (* every match is replaced *)
```

```mathematica
In[1]:= StringReplace["a1b2c3", DigitCharacter -> "#"]  (* a string-pattern rule *)
```

```mathematica
In[1]:= StringReplace["2024-01", RegularExpression["(\\d+)-(\\d+)"] -> "$2.$1"]  (* $n reorders captured groups *)
```

### Notes

`StringReplace` scans left to right, replacing each non-overlapping match by its
rule's right-hand side (a string, with `$0`/`$n` expanded) and copying the
unmatched text. With several rules, the leftmost match wins at each position,
ties broken by rule order.

A zero-width match (e.g. `\b`) inserts the replacement without consuming a
character. A list of subjects threads.
