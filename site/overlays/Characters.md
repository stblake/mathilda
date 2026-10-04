### Worked examples

```mathematica
In[1]:= Characters["cat"]  (* a string becomes a list of length-1 strings *)
```

```mathematica
In[1]:= StringJoin[Characters["abcd"]]  (* Characters and StringJoin are inverse *)
```

```mathematica
In[1]:= Length[Characters["hello"]]  (* one element per byte, as in StringLength *)
```

### Notes

`Characters` is byte-oriented — one element per `char`, with no UTF-8 codepoint
decoding — so its result length always equals `StringLength`. It is `Listable`,
so a list of strings gives a list of character lists, and the empty string gives
`{}`.
