---
name: test-generation
description: Use when the user asks to add tests, improve test coverage, generate a test suite for a module, or write tests for a specific function or class.
allowed-tools: Read, Grep, Glob, Bash, Edit, Write
---

# test-generation

Full-loop test writing: read the code under test, identify untested
behaviours and edge cases, write tests that assert observable outcomes
(not implementation details), and leave the suite green.

## When to invoke

- User says "add tests for X", "write a test suite", "increase coverage".
- A function or module is new and has no tests yet.
- Existing tests are brittle (break when internals change) or incomplete.
- User asks to cover a specific edge case or bug that was recently fixed.

## Method

1. **Read the code.** Understand the public interface: inputs, return values,
   side effects, error cases. Read existing tests to match style.
2. **Map the behaviour surface.** List: happy paths, boundary values (empty,
   zero, max), error paths (invalid input, dependency failure), state changes,
   and any concurrency concerns.
3. **Pick the right test type** for each behaviour:
   - **Unit** — pure functions, no I/O; mock external dependencies.
   - **Integration** — real DB / filesystem / queue; use fixtures.
   - **End-to-end** — full request/response cycle; use `TestClient` (FastAPI) or `httpx`.
4. **Write the tests.** One assertion per test concept. Name the test so it reads
   as a specification: `test_<unit>_<scenario>_<expected>`.
5. **Run the suite.** Fix any failures before reporting done.
6. **Check coverage.** `pytest --cov` — aim for 80%+ on the changed module.

## Skeletons

### Python / pytest

```python
# test_<module>.py
import pytest
from myapp.<module> import <function_or_class>


# --- Happy path ---

def test_<function>_returns_expected_result():
    result = <function>(valid_input)
    assert result == expected_output


# --- Boundary ---

@pytest.mark.parametrize("value,expected", [
    (0,   ...),
    (-1,  ...),
    (None, ...),
])
def test_<function>_boundary_values(value, expected):
    assert <function>(value) == expected


# --- Error path ---

def test_<function>_raises_on_invalid_input():
    with pytest.raises(ValueError, match="descriptive message"):
        <function>(invalid_input)


# --- Side effects / state ---

def test_<function>_writes_to_db(db_session):
    <function>(db_session, payload)
    record = db_session.query(Model).one()
    assert record.field == expected_value


# --- Async ---

@pytest.mark.asyncio
async def test_<async_function>_happy_path():
    result = await <async_function>(input)
    assert result.status == "ok"
```

### TypeScript / Vitest

```typescript
import { describe, it, expect, vi } from 'vitest';
import { myFunction } from '../src/myModule';

describe('myFunction', () => {
  it('returns expected result for valid input', () => {
    expect(myFunction('valid')).toBe('expected');
  });

  it('throws for invalid input', () => {
    expect(() => myFunction(null)).toThrow('descriptive message');
  });

  it('calls dependency once with correct args', () => {
    const spy = vi.spyOn(dependency, 'method');
    myFunction('input');
    expect(spy).toHaveBeenCalledWith('input');
    expect(spy).toHaveBeenCalledOnce();
  });
});
```

### React Testing Library

```tsx
import { render, screen, userEvent } from '@testing-library/react';
import { MyComponent } from '../src/MyComponent';

it('shows error message when form submitted empty', async () => {
  render(<MyComponent />);
  await userEvent.click(screen.getByRole('button', { name: /submit/i }));
  expect(screen.getByRole('alert')).toHaveTextContent('Required');
});
```

## Guardrails

- Follow [`../../shared/guidelines/pr-and-review.md`](../../shared/guidelines/pr-and-review.md):
  tests assert observable behaviour, not internal names.
- Never use `.skip` or `xfail` without a linked issue.
- Don't mock your own module's internals — only mock external dependencies.
- Integration tests must use isolated state (transaction rollback, temp DB, tmp dir).
- Don't share mutable state between tests — each test must be independently runnable.

## Anti-patterns

- Tests that only verify the call happened (`assert mock.called`) without
  checking what was called with or what was returned.
- Snapshot tests for logic-heavy components — they catch everything and explain nothing.
- Testing private methods directly — refactor if the private logic is complex enough to need it.
- `time.sleep` in tests — use `freeze_time` or mock the clock.
- Testing third-party libraries instead of your code.

## Extending

- Mutation testing: `mutmut run` to verify tests actually catch regressions.
- Property-based testing: `hypothesis` for functions with large input spaces.
- Contract testing: `pact` for API consumer/provider boundary validation.
- Coverage enforcement: add `--cov-fail-under=80` to CI so coverage can't regress.
