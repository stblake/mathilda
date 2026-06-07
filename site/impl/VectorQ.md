---
source: src/list.c
---
`builtin_vectorq` (`src/list.c`) returns `True` when the argument is a `List` (`is_listq`) none of whose elements is itself a `List`. With a second argument `test`, each element is instead required to satisfy `test[elem]` (evaluated, must yield `True`). Returns `False` as soon as a check fails.
