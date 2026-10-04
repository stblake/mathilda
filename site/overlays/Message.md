### Worked examples

```mathematica
In[1]:= Quiet[Check[Message[g::x]; 7, flagged]]  (* Message fires, so the enclosing Check takes its failure branch *)
```

```mathematica
In[1]:= Message[myfun::warn]  (* on its own a Message just returns Null *)
```

### Notes

`Message[sym::tag, e1, …]` notes that a diagnostic fired — so an enclosing `Check`
registers it — and, unless messages are suppressed, prints the text defined for
`sym::tag`. It is `HoldFirst` (the message name is held, so that evaluating it
yields its template string) and always returns `Null`.

The *firing* is the point: it is the hook that `Check` detects and `Quiet`
silences, both through the same message funnel. The first example shows a bare
`Message` inside a `Quiet[Check[…]]`, where it flips the result to the failure
branch without printing anything.
