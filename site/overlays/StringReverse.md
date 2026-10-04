### Worked examples

```mathematica
In[1]:= StringReverse["abcdef"]  (* reverse the characters *)
```

```mathematica
In[1]:= StringReverse["racecar"]  (* a palindrome is its own reverse *)
```

```mathematica
In[1]:= StringReverse[{"ab", "cd"}]  (* Listable, so a list threads *)
```

### Notes

`StringReverse` reverses the bytes of a string into a fresh buffer. It is
`Listable`, so the evaluator threads it over a list before the builtin runs and
each call sees a single string.

A call with a number of arguments other than one emits `StringReverse::argx`; a
non-string argument leaves the call unevaluated.
