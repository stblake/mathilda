---
source: src/linalg/construct.c
---
`builtin_identitymatrix` accepts either an integer `n` (square `n×n`) or a pair `{m, n}` of integers, and constructs a `List` of `List`s with `Integer` `1` on the main diagonal (`i == j`) and `0` elsewhere. Non-integer or malformed dimension specs are returned unevaluated (`expr_copy(res)`). The output is exact integer entries; no numeric or symbolic processing is involved.
