---
source: src/list.c
---
**Algorithm.** `builtin_rotateright` cyclically shifts elements toward the back by `n` (default
1). It negates the shift spec (an integer or a per-level `List` of integers) and delegates to
the same `rotate_rec` worker used by `RotateLeft`, so the offset wrapping and per-level nested
behaviour are identical with the opposite sign.
