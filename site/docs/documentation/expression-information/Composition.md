# Composition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Composition[f1, f2, f3, ...]
    represents a composition of the functions f1, f2, f3, ....

Composition allows you to build up compositions of functions which can
later be applied to specific arguments. Applied to arguments, the
composition acts innermost-first:
    Composition[f, g, h][x, y]  ->  f[g[h[x, y]]].

Composition has the attributes Flat and OneIdentity.
Composition can be entered in the form f1 @* f2 @* ....

Composition objects containing Identity or InverseFunction[f] are
automatically simplified when possible:
    Composition[]                       ->  Identity
    Composition[f]                      ->  f
    Composition[f, Identity, g]         ->  Composition[f, g]
    Composition[f, InverseFunction[f]]  ->  Identity.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Composition[f, g, h][x, y]
Out[1]= f[g[h[x, y]]]

In[2]:= f @* g @* h @ x
Out[2]= f[g[h[x]]]

In[3]:= Composition[f, Identity, g]
Out[3]= Composition[f, g]

In[4]:= Composition[f, InverseFunction[f]][x]
Out[4]= x

In[5]:= Composition[f, g] @* Composition[a, b]
Out[5]= Composition[f, g, a, b]
```

## Implementation notes

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
