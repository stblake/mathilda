---
description: Generate a test suite for a specific file, module, or function. Pass the target path as an argument.
argument-hint: <file-or-module-path>
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

Generate comprehensive tests for the code at `$ARGUMENTS`.

## Step 1 — Read the target

Read the file at `$ARGUMENTS`. Identify:
- All public functions, methods, and classes.
- Their signatures, return types, and documented/inferred side effects.
- Any existing tests (check `tests/`, `__tests__/`, or `*.test.*` alongside the file).

```bash
cat "$ARGUMENTS"
# find existing tests
find . -name "test_$(basename $ARGUMENTS .py).py" -o -name "*.test.ts" | head -20
```

## Step 2 — Map the behaviour surface

For each public unit, list:
1. **Happy path** — normal input → expected output.
2. **Boundary values** — empty, zero, None/null, very large, exact limits.
3. **Error paths** — invalid input, missing dependencies, upstream failures.
4. **State effects** — DB writes, file mutations, external calls.
5. **Concurrency** (if async) — what happens under concurrent calls?

## Step 3 — Write the tests

Use the [`test-generation`](../skills/test-generation/SKILL.md) skill:
- One test per scenario, named `test_<unit>_<scenario>_<expected>`.
- Mock only external dependencies (DB, HTTP, filesystem) — not your own code.
- Each test must be independently runnable (no shared mutable state).
- Use `pytest.mark.parametrize` for value tables.

Write tests to a file alongside the target or in the project's `tests/` directory,
matching the existing convention.

## Step 4 — Run the suite

```bash
pytest path/to/test_file.py -v
# or for TypeScript:
# npx vitest run path/to/test_file.test.ts
```

Fix any failures before reporting done.

## Step 5 — Check coverage

```bash
pytest path/to/test_file.py --cov=$(dirname $ARGUMENTS) --cov-report=term-missing
```

Report coverage. If below 80% on the target module, identify what's missing
and add it.

## Output

Report:
- List of test cases written (unit → scenario).
- Coverage percentage on the target module.
- Any behaviours that could not be tested without additional fixtures or refactoring.

## Related

- Skill: [`test-generation`](../skills/test-generation/SKILL.md)
- Command: [`review-pr`](review-pr.md) — verifies test quality in code review.
