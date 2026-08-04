# Compile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Compile[{x, ...}, expr] or Compile[{{x, _Real}, ...}, expr] builds a CompiledFunction that evaluates expr over machine numbers (types _Real, _Integer, _Complex; default _Real), falling back to the interpreter for symbolic arguments or non-compilable bodies. With RuntimeAttributes -> Listable the object threads over List arguments; the default is RuntimeAttributes -> {}. RuntimeOptions -> {"CatchMachineIntegerOverflow" -> False} (or the shorthand RuntimeOptions -> "Speed") lets machine-integer arithmetic wrap instead of falling back to the interpreter, which is faster and gives a different answer from the interpreter once a result leaves the machine-integer range; the default True never does.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= f = Compile[{{x, _Real}}, x^2 + 1]
Out[1]= CompiledFunction[{x}, x^2 + 1]

In[2]:= f[3.0]
Out[2]= 10.0

In[3]:= f[a]                 (* symbolic argument -> interpreter fallback *)
Out[3]= 1 + a^2

In[4]:= Compile[{{z, _Complex}}, z^2][1.0 + 2.0 I]
Out[4]= -3.0 + 4.0*I

In[5]:= Compile[{{m, _Real, 2}}, m[[All, 1]]][{{1., 2.}, {3., 4.}}]
Out[5]= {1.0, 3.0}
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/compile/compiled_function.c`](https://github.com/stblake/mathilda/blob/main/src/compile/compiled_function.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
