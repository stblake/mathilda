# Cubics

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Cubics is an option for Solve that controls whether cubic
    equations are solved via explicit radical formulas
    (Cubics -> True) or returned as held Root[] objects
    (default Cubics -> False).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Cubics` is not a function — it is an option symbol for `Solve`/`Roots` (and forwarded by the eigenvalue solver). It has no builtin handler; it is interned in `sym_names.c` and registered only with a docstring. Option parsing in `solve.c` (`is_*_option` / the option setter) reads `Cubics -> True|False` into `opts->poly.cubics_radical`, which controls whether degree-3 factors are solved by radicals (Cardano) or returned as held `Root[]` objects. The actual radical-vs-`Root` decision lives in the polynomial solver (`src/poly/solvepoly.c`); the default is `Cubics -> False`.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
