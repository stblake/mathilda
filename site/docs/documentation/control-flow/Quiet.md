# Quiet

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Quiet[expr] evaluates expr with messages suppressed and returns its value. Quiet[expr, spec] suppresses all messages regardless of spec.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Quiet[1/0]
Out[1]= ComplexInfinity

In[2]:= Quiet[Log[0]]
Out[2]= -Infinity
```

### Applications (2)

The Power::infy message is silenced; the value is unchanged

```mathematica
In[3]:= Quiet[1/0]
Out[3]= ComplexInfinity
```

The value comes back with no printed warning

```mathematica
In[4]:= Quiet[Log[0]]
Out[4]= -Infinity
```

## Implementation notes

**Algorithm.** `Quiet` is `HoldAll, Protected`, so the argument is evaluated
*under* the suppression rather than before it. `builtin_quiet` brackets the
evaluation between `mth_msg_suppress_push()` and `mth_msg_suppress_pop()` and
returns the value (an in-flight `Throw` sentinel propagates unchanged). The
optional `spec` of the two-argument form — a message name or list — is accepted
and ignored: all messages are suppressed, a harmless superset of any requested
set.

**Data structures.** A single static nesting depth `g_msg_suppress_depth`
(`> 0` ⇒ suppress). The depth is a *count*, so `Quiet` nests correctly. Inside the
funnel `mth_message_v` still calls `mth_msg_note_fired()` before testing
`mth_msg_suppressed()`, so a message under `Quiet` still *fires* (an enclosing
`Check` sees it) and only its printing is silenced. A paired
`mth_msg_suppress_depth_save`/`_load` lets `TimeConstrained`'s `siglongjmp` restore
the depth if a timeout unwinds out of a `Quiet` region, which would otherwise
leave messages silenced for the rest of the session.

**Complexity / limits.** `O(1)` push/pop around the inner evaluation. `Quiet`
affects only the *display* of diagnostics; it changes neither the value computed
nor whether a message is considered to have fired.

- `HoldAll`, so the argument is evaluated under the suppression, not before it.
- A message still *fires* while suppressed (an enclosing `Check` sees it); only the printing is silenced.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Check](../../control-flow/Check/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
- Tests: [`tests/test_parallelmixedspecial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedspecial.c)

## Notes & additional examples

### Notes

`Quiet[expr]` evaluates `expr` and returns its value with any messages suppressed.
It is `HoldAll`, so the suppression is in force *during* the evaluation. The
optional second argument (a message name or list) is accepted and ignored — all
messages are suppressed, which is a harmless superset of any requested set.

Suppression silences only the *printing*. A message still fires while quieted, so
an enclosing `Check` still sees it — which is exactly what the common
`Quiet[Check[expr, failexpr]]` pattern relies on.
