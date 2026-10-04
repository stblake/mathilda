# Check

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Check[expr, failexpr] returns failexpr if a message is generated while evaluating expr, otherwise the value of expr.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Check[2 + 2, failed]
Out[1]= 4

In[2]:= Check[1/0, failed]
Out[2]= failed

In[3]:= Quiet[Check[1/0, caughtit]]
Out[3]= caughtit
```

### Applications (3)

No message fires, so the value passes through

```mathematica
In[4]:= Check[2 + 2, failed]
Out[4]= 4
```

The division emits Power::infy, so failexpr is returned

```mathematica
In[5]:= Check[1/0, failed]
Out[5]= failed
```

Detect the failure without printing its message

```mathematica
In[6]:= Quiet[Check[1/0, caughtit]]
Out[6]= caughtit
```

## Implementation notes

**Algorithm.** `Check` is `HoldAll, Protected`, so its argument arrives
unevaluated and is evaluated under `Check`'s watch. `builtin_check` snapshots the
global message-fired counter (`mth_msg_fired_count`), evaluates `expr`, and
compares the counter afterward: if any diagnostic fired during the evaluation it
frees `expr`'s value and returns the evaluated `failexpr`, otherwise it returns
`expr`'s value. The optional third `spec` argument is accepted and ignored — any
message counts.

**Data structures.** A single `unsigned long` counter incremented by
`mth_msg_note_fired()` from inside the message funnel (`mth_message_v`). Crucially
the funnel notes a firing *even while `Quiet[]` is suppressing the print*, which
is what makes the `Quiet[Check[expr, failexpr]]` idiom — detect a failure without
showing its message — work. This is why every user-facing diagnostic must route
through the funnel: a raw `fprintf(stderr, …)` would be invisible to this counter
and `Check` would silently take the wrong branch.

**Complexity / limits.** `O(1)` bookkeeping around the inner evaluation. A `Throw`
inside `expr` is not a message: `builtin_check` tests `eval_is_inflight_throw` on
the value first and lets the sentinel propagate unchanged, so a throw escapes
`Check` rather than being reported as a failure.

- `HoldAll`. A `Throw` inside `expr` propagates (it is not a message).
- Typically wrapped as `Quiet[Check[expr, failexpr]]` to detect a failure without printing its message.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Throw](../../control-flow/Throw/)

- Source: [`src/message.c`](https://github.com/stblake/mathilda/blob/main/src/message.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)
- Tests: [`tests/test_findmin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin.c)
- Tests: [`tests/test_findmin_cobyla.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_cobyla.c)
- Tests: [`tests/test_findmin_slsqp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_slsqp.c)

## Notes & additional examples

### Notes

`Check[expr, failexpr]` returns `failexpr` if *any* message is generated while
`expr` is evaluated, otherwise the value of `expr`. It is `HoldAll`, so `expr` is
evaluated under `Check`'s watch rather than before it.

`Check` keys off whether a message *fired*, which is independent of whether it was
*printed* — so `Quiet[Check[expr, failexpr]]` is the idiom for detecting a failure
silently. A `Throw` inside `expr` is not a message: it propagates straight out of
`Check` rather than tripping the failure branch.
