---
source: src/list.c
---
`builtin_matrixq` (`src/list.c`) returns `True` when the argument is a non-empty `List` of equal-length `List` rows whose entries are non-lists; with an optional `test` argument each entry must instead satisfy `test[entry]` (evaluated to `True`). Any ragged row, empty outer list, or failing entry yields `False`.
