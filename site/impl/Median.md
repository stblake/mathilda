---
source: src/stats.c
---
**Algorithm.** `builtin_median` requires a `List`. If the first element is itself a `List` it treats the input as a matrix/tensor and reduces column-wise through `apply_columnwise` (`Map[Median, Transpose[...]]`). For a 1-D vector it first verifies every element is a real numeric via the helper `is_real_numeric` (which checks `NumericQ` and `FreeQ[#, I]`); non-real data prints `Median::rectn` and leaves the call unevaluated. It then evaluates `Sort[data]`: for odd `n` it returns the middle element (`sorted[n/2]`); for even `n` it returns `(sorted[n/2-1] + sorted[n/2]) / 2`, built as `Plus` then `Divide` and re-evaluated so the result stays exact (rational) when the inputs are exact. `ATTR_PROTECTED`.
