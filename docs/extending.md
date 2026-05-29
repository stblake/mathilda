# Extending Mathilda

How to add new built-in functions, modules, pattern constructs, internal `.m`
rules, and parser operators. The architectural background lives in
[`SPEC.md`](../SPEC.md); the contributor workflow (testing, valgrind, coding
standards) lives in [`CLAUDE.md`](../CLAUDE.md).

---

## 1. Adding a New Built-in Function

The most common extension. Six steps.

### Step 1 — Write the C implementation

Built-ins have signature:

```c
Expr* builtin_myfunc(Expr* res);
```

Ownership contract (mirrors [`SPEC.md` §4](../SPEC.md#4-memory-management)):

- Return `NULL` when the function cannot evaluate (wrong arity, symbolic args,
  etc.). The evaluator retains ownership of `res`; **do not** free it.
- Return a new `Expr*` on success. The evaluator frees `res` for you. Do **not**
  call `expr_free(res)` yourself — see `feedback_builtin_res_ownership` in
  auto-memory if in doubt.

Template:

```c
Expr* builtin_myfunc(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_INTEGER) return NULL;

    int64_t v = arg->data.integer;
    return expr_new_integer(v * v);
}
```

Handle `EXPR_BIGINT` alongside `EXPR_INTEGER` when the function should work for
large integers, and `EXPR_REAL` when it should handle floats. For
`Rational[n, d]` or `Complex[re, im]` operands check `EXPR_FUNCTION` with the
appropriate head.

### Step 2 — Declare it

Add the prototype to the module's `.h` file:

```c
Expr* builtin_myfunc(Expr* res);
```

### Step 3 — Register in the module init

```c
void mymodule_init(void) {
    symtab_add_builtin("MyFunc", builtin_myfunc);
    symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("MyFunc",
        "MyFunc[x]\n\tComputes the square of x.");
}
```

Docstrings must be **terse** — no inline examples; put those in the relevant
`docs/spec/builtins/*.md` file and the changelog entry.

If this is a new module, add `mymodule_init()` to `core_init()` in `core.c`.

### Step 4 — Assign attributes

| Attribute | When to use |
|-----------|-------------|
| `ATTR_PROTECTED` | Almost always — prevents user redefinition. |
| `ATTR_LISTABLE` | Threads over lists: `f[{a, b}]` → `{f[a], f[b]}`. |
| `ATTR_NUMERICFUNCTION` | Operates on numbers. |
| `ATTR_FLAT` | Associative: `f[f[a, b], c]` → `f[a, b, c]`. |
| `ATTR_ORDERLESS` | Commutative: `f[b, a]` → `f[a, b]`. |
| `ATTR_HOLDALL` / `ATTR_HOLDFIRST` / `ATTR_HOLDREST` | Suppress evaluation of some/all arguments. |

`*Q` predicates must return `True` or `False`, never symbolic and never `NULL`
(see `feedback_q_predicates_return_bool` in auto-memory).

### Step 5 — Write tests

Add to an existing `tests/test_*.c` or create a new one and register it in
`tests/CMakeLists.txt`.

```c
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

void test_myfunc_basic(void) {
    assert_eval_eq("MyFunc[5]", "25", 0);
}

void test_myfunc_symbolic(void) {
    assert_eval_eq("MyFunc[x]", "MyFunc[x]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_myfunc_basic);
    TEST(test_myfunc_symbolic);
    return 0;
}
```

Helpers in `tests/test_utils.h`: `TEST`, `ASSERT`, `ASSERT_STR_EQ`,
`ASSERT_MSG`, `assert_eval_eq(input, expected, is_fullform)`.

### Step 6 — Update documentation

Edit the matching `docs/spec/builtins/*.md` and add a changelog entry under
the current week's `docs/spec/changelog/<YYYY-MM-DD>.md`, where
`<YYYY-MM-DD>` is the Monday of the ISO week (Mon – Sun) the change lands
in. Create the file with a `# Changelog: week of <Mon> (Mon) – <Sun> (Sun)`
header (and a corresponding row in `Mathilda_spec.md`'s changelog table) if
it does not yet exist. The top-level `Mathilda_spec.md` only needs touching
if you introduce a new top-level category or a new weekly changelog file.

---

## 2. Adding a New Module

1. Create `src/mymodule.c` and `src/mymodule.h` (or a subdirectory analogous to
   `src/poly/`, `src/linalg/`, `src/simp/`, `src/calculus/`).
2. Define and declare `void mymodule_init(void)`.
3. `#include "mymodule.h"` in `core.c` and call `mymodule_init()` from
   `core_init()`.
4. New `.c` files in `src/` (and subdirectories already wired into the
   makefile) are picked up automatically by `SRC = $(wildcard $(SRC_DIR)/*.c)`.

---

## 3. Adding a New Pattern Construct

1. Extend `match.c` to recognize the new pattern during tree unification.
2. If the construct needs new syntax, extend `parse.c` — define a precedence,
   add token recognition, and emit a standard `Expr*` function call.
3. Add tests in `tests/` covering both match and non-match scenarios.
4. If the evaluator would prematurely evaluate the new pattern expression,
   guard the surrounding builtin with the appropriate Hold attribute.

---

## 4. Adding Rules via Internal `.m` Files

For mathematical identities that are naturally rewrite rules:

1. Create `src/internal/yourfile.m`.
2. Add `Get["src/internal/yourfile.m"]` to `src/internal/init.m`.
3. Define rules using normal Mathematica syntax:

```mathematica
MyFunc[0] := 0;
MyFunc[x_ + y_] := MyFunc[x] + MyFunc[y];
MyFunc[n_Integer x_] := n MyFunc[x];
```

This approach is ideal for derivative rules, identities, and simplification
rules where the pattern-matching system does the heavy lifting (see
`src/internal/deriv.m`).

---

## 5. Adding New Parser Operators

1. In `parse.c`, define a new precedence constant (consult the Mathematica
   documentation for the standard value).
2. Add token recognition in the lexer section of the parser.
3. Implement the infix / prefix / postfix logic in `parse_expression_prec()`.
4. The operator should produce a normal function-call expression (e.g., `a <> b`
   should produce `StringJoin[a, b]`).
