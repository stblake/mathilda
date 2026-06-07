---
source: src/core.c
---
**Algorithm.** `builtin_prepend` builds a new function node with the same head whose first
element is a copy of the new element followed by copies of all original arguments. Two-argument
form only; returns `NULL` if the first argument is an atom. (`PrependTo` is the mutating
variant that writes the result back to a symbol's OwnValue.)
