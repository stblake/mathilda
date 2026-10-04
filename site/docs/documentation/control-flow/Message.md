# Message

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Message[sym::tag, args] prints the named message unless messages are suppressed.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Quiet[Check[Message[g::x]; 7, flagged]]
Out[1]= flagged
```

Fires a diagnostic; returns Null

```mathematica
In[2]:= Message[myfun::warn]
```

### Applications (2)

Message fires, so the enclosing Check takes its failure branch

```mathematica
In[3]:= Quiet[Check[Message[g::x]; 7, flagged]]
Out[3]= flagged
```

On its own a Message just returns Null

```mathematica
In[4]:= Message[myfun::warn]
```

## Implementation notes

**Algorithm.** `Message` is `HoldFirst, Protected`: the message name `sym::tag`
(a `MessageName`) is held, so that evaluating it yields the template string the
user defined for it, while the trailing arguments evaluate normally.
`builtin_message` first calls `mth_msg_note_fired()` unconditionally — so an
enclosing `Check` registers the diagnostic whether or not it is displayed — then,
unless messages are suppressed, evaluates the first argument and, when it resolves
to a string, prints it to stderr. It always returns `Null`.

**Data structures.** It shares the message subsystem's single fired-counter and
suppression depth with `Check` and `Quiet`; there is no per-message table beyond
the `MessageName` own-value that holds a template string.

**Complexity / limits.** Argument substitution into the template (the `` `1` ``
slots) is **not** performed: the port uses `Message` only in error branches that
are immediately followed by `Throw`, so the printed text is diagnostic rather than
load-bearing. The firing is what matters — it is the hook that `Check` detects and
`Quiet` silences, both through the same funnel.

- `HoldFirst`.
- Returns `Null`.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [Check](../../control-flow/Check/), [HoldFirst](../../other-advanced/HoldFirst/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)

## Notes & additional examples

### Notes

`Message[sym::tag, e1, …]` notes that a diagnostic fired — so an enclosing `Check`
registers it — and, unless messages are suppressed, prints the text defined for
`sym::tag`. It is `HoldFirst` (the message name is held, so that evaluating it
yields its template string) and always returns `Null`.

The *firing* is the point: it is the hook that `Check` detects and `Quiet`
silences, both through the same message funnel. The first example shows a bare
`Message` inside a `Quiet[Check[…]]`, where it flips the result to the failure
branch without printing anything.
