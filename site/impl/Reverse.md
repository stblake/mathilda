---
source: src/list.c
---
**Algorithm.** `builtin_reverse` recurses through the expression with `reverse_rec`, reversing
the argument order at the levels selected by an optional level spec (`should_reverse_at_level`
matches an integer level, or any level listed in a `List` spec; default is level 1). At each
visited node it builds a new function with the same head, drawing children either forward or
reversed depending on whether the current level is selected, and recursing into each child.
