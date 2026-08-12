# List

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

{e1, e2, ...} or List\[e1, e2, ...\] represents an ordered list of the elements ei. List is the fundamental container head: vectors are lists, matrices are lists of lists, and the structural operators (Part, Map, Take, Drop, Length, ...) act on List. Elements are evaluated normally and kept in the given order (List has no Orderless attribute). The parser writes the {...} syntax to List, and the printer renders List\[...\] back as {...}.

## Examples

_No verified examples yet for this function._

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Dot 6x6 x 6x6 x 10000 | 338 s | 6.36 s | 4.01 s |
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Inverse 3x3 x 5000 | 2.7 s | 4.75 s | 7.67 s |

## Implementation notes

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_array_flatten.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_flatten.c)
- Tests: [`tests/test_backtrack.c`](https://github.com/stblake/mathilda/blob/main/tests/test_backtrack.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_constant_array.c`](https://github.com/stblake/mathilda/blob/main/tests/test_constant_array.c)
