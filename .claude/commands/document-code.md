---
description: Add or update docstrings and module-level documentation for a file or directory. Pass the target path as an argument.
argument-hint: <file-or-directory-path>
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

Add or update documentation for the code at `$ARGUMENTS`. The goal is
*why* and *contract*, not *what* — the code already shows what.

## Step 1 — Read and assess

```bash
cat "$ARGUMENTS"
```

Identify:
- Public functions/classes with missing or outdated docstrings.
- Non-obvious module-level context (what does this module do? what are its dependencies?).
- Complex logic that needs a "why" comment (algorithm choice, workaround, invariant).

Skip:
- Private helpers where the name is self-documenting.
- Trivial getters/setters.
- Anything the type signatures already express.

## Step 2 — Write docstrings

Match the project's existing docstring style. Default to Google-style for Python:

```python
def segment_customers(df: pd.DataFrame, *, thresholds: list[float] | None = None) -> pd.DataFrame:
    """Assign each customer to a revenue segment.

    Uses log-linear quantile boundaries by default. Pass custom ``thresholds``
    (sorted ascending, in the same currency as ``revenue``) to override.

    Args:
        df: DataFrame with at minimum ``customer`` (str) and ``revenue`` (float) columns.
        thresholds: Optional explicit segment boundaries. Defaults to [10_000, 100_000, 1_000_000].

    Returns:
        Input DataFrame with a new ``segment`` column:
        ``"small"`` / ``"mid"`` / ``"large"`` / ``"strategic"``.

    Raises:
        ValueError: If ``revenue`` column is missing or contains negative values.
    """
```

TypeScript / JSDoc:

```typescript
/**
 * Assigns each customer to a revenue segment.
 *
 * Uses log-linear quantile boundaries by default.
 *
 * @param customers - Array of customer records with `revenue` field.
 * @param thresholds - Optional sorted ascending breakpoints (same currency as revenue).
 * @returns The input array with a `segment` field added to each record.
 * @throws {Error} If any revenue value is negative.
 */
```

## Step 3 — Add module-level docstring

For Python modules that lack one:

```python
"""Customer segmentation utilities.

Provides revenue-based segmentation logic used by the analytics pipeline.
All monetary values are in USD unless the caller specifies otherwise.

Depends on: pandas >= 2.0, numpy >= 1.26.
"""
```

## Step 4 — Add "why" comments for non-obvious logic

```python
# Clamp to [0, 1] before log to avoid -inf; values < 1e-9 treated as zero
log_revenue = np.log1p(np.clip(df["revenue"], 0, None))
```

## Step 5 — Verify no regressions

```bash
# Python: check for syntax errors introduced
python -m py_compile "$ARGUMENTS"

# TypeScript: type-check
npx tsc --noEmit 2>/dev/null || true
```

## Output

Report:
- N docstrings added / updated.
- Any complex logic that was left uncommented with a reason (e.g., "explanation
  would duplicate the type annotation").
- Any functions that looked complex enough to warrant a refactor instead of
  a comment (flag only — don't refactor without being asked).

## References

- [`../shared/guidelines/pr-and-review.md`](../shared/guidelines/pr-and-review.md):
  "Document the why, not the what."
