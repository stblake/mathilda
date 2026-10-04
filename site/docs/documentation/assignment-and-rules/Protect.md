# Protect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Protect[s1, s2, ...]`**

sets the attribute Protected for the named symbols and returns the list of their names. Protect\[{s1, s2, ...}\] accepts a list of specs.

<details>
<summary>Notes</summary>

Protect has attribute HoldAll; Locked symbols are not affected.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= f[x_] := x^2; Protect[f]
Out[1]= {"f"}

In[2]:= Unprotect[f]
Out[2]= {"f"}
```

### Applications (2)

Returns the names whose Protected attribute was newly set

```mathematica
In[3]:= Protect[foo, bar]
Out[3]= {"foo", "bar"}
```

```mathematica
In[4]:= MemberQ[Attributes[foo], Protected]
Out[4]= True
```

## Implementation notes

**Algorithm.** `builtin_protect` (`src/core.c`) calls the shared driver
`core_protect_unprotect(res, protecting = true)`. That driver walks the argument
list — each spec being a symbol, a string, or a flat `List` of them — and applies
`core_protect_one` to every name. `core_protect_one` leaves a `Locked` symbol
untouched and does nothing if the symbol is already `Protected`; otherwise it
sets the `ATTR_PROTECTED` bit and bumps the rule epoch (invalidating the
evaluation cache), returning `true` only when the bit was *newly* set.

The driver collects the names whose state actually changed into a growable
`Expr**` buffer and returns them as a `List` of strings — matching the Wolfram
convention where `Protect` reports exactly what it altered, so re-protecting an
already-protected symbol yields `{}`. The `Protected` attribute is what the
evaluator and `Set` consult to refuse redefinition of a symbol.

**Attributes & limits.** `Protect` carries `HoldAll | Protected`, so its
arguments reach the handler as unevaluated symbol names.

- Both have attributes `{HoldAll, Protected}` and hold their arguments.
- Neither affects symbols with the attribute `Locked`.
- The typical sequence for adding rules to an existing symbol is
  `Unprotect[f]; definition; Protect[f]`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Unprotect](../../assignment-and-rules/Unprotect/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)

## Notes & additional examples

### Notes

`Protect[s1, s2, ...]` sets the `Protected` attribute on each symbol, which is
what the evaluator and `Set` consult to refuse redefinition. It returns a list of
the names (as strings) whose state actually *changed*, so re-protecting an
already-protected symbol gives `{}`.

Arguments may be symbols, strings, or a flat list of them. A `Locked` symbol is
left untouched. `Protected` is the attribute every built-in carries; `Unprotect`
is its inverse and the usual first step before extending a built-in's behaviour.
