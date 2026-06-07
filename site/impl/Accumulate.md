---
source: src/list.c
---
`builtin_accumulate` returns the list of cumulative sums (prefix sums), keeping the input expression's head. The default path folds `running = Plus[running, elem]` left to right via the evaluator, so it works on any addable elements (integers, rationals, symbolics, matrix rows). When the optional `Method -> "CompensatedSummation"` is supplied and every element is a machine number, it instead runs **Kahan compensated summation** in `double` precision (tracking a running correction term `c`), emitting `Real` partial sums. An empty list returns a copy unchanged.
