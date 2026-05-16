# Mathilda Technical Specifications & Developer Guide

## 1. Introduction

**Mathilda** is a tiny, symbolic computer algebra system (CAS) written in C, heavily inspired by the core architecture and evaluation semantics of Mathematica (the Wolfram Language). 

This document serves as a comprehensive technical overview for new developers (or AI assistants) joining the project. It outlines the core architecture, data structures, evaluation semantics, and conventions used throughout the codebase.

## 2. Core Architecture

The system is composed of several interdependent subsystems:
1.  **Expression Representation** (`expr.c`, `expr.h`)
2.  **Parser** (`parse.c`, `parse.h`)
3.  **Symbol Table** (`symtab.c`, `symtab.h`)
4.  **Evaluator** (`eval.c`, `eval.h`)
5.  **Pattern Matcher** (`match.c`, `match.h`)
6.  **Rule Engine** (`replace.c`, `replace.h`)
7.  **Built-in Function Libraries** (e.g., `core.c`, `list.c`, `arithmetic.c`, `part.c`)

### 2.1. Expression Representation (`Expr`)

Everything in Mathilda is an expression, represented by the `Expr` struct. It is implemented as a tagged union.

*   **Types (`ExprType`)**:
    *   `EXPR_INTEGER`: A 64-bit signed integer (`int64_t`).
    *   `EXPR_REAL`: A 64-bit floating-point number (`double`).
    *   `EXPR_SYMBOL`: A named identifier (e.g., `x`, `Plus`, `True`).
    *   `EXPR_STRING`: A string literal.
    *   `EXPR_FUNCTION`: A compound expression consisting of a `head` (which is itself an `Expr*`) and an array of `args` (`Expr**`).

All data structures are immutable-by-convention during evaluation. Functions that manipulate expressions generally return new allocated trees or unmodified references if no changes were made. 

**Memory Management:** Explicit manual memory management is required. `expr_new_*` allocates nodes. `expr_copy` performs a deep copy. `expr_free` performs a deep deallocation. **Crucial Rule:** Built-in functions take ownership of the expression passed to them. If a built-in returns a new expression, it *must* free the input expression (or reuse its nodes). If it returns `NULL` (meaning it remains unevaluated), the evaluator retains ownership.

### 2.2. Symbol Table (`symtab.c`)

The symbol table (`SymbolDef`) stores the global state associated with symbols. Each symbol can have:
*   **Attributes**: Bitflags (e.g., `ATTR_FLAT`, `ATTR_ORDERLESS`, `ATTR_LISTABLE`, `ATTR_PROTECTED`, `ATTR_HOLDALL`, `ATTR_NUMERICFUNCTION`).
*   **OwnValues**: Immediate assignments (e.g., `x = 5`).
*   **DownValues**: Delayed assignments with pattern matching (e.g., `f[x_] := x^2`).
*   **Builtin C Function**: A pointer to a C function (`BuiltinFunc`) that implements native evaluation logic (e.g., `builtin_plus`).

### 2.3. Parser (`parse.c`)

The parser translates raw string input into `Expr` trees.
*   **Lexical Analysis**: Handled inline within the parsing routines by advancing a `ParserState` pointer and skipping whitespace.
*   **Pratt Parsing**: Operator precedence (Infix, Prefix, Postfix) is implemented using a Pratt parser (`parse_expression_prec`). Operator precedence values exactly mirror Mathematica's standard precedences (e.g., `OP_CALL` is 1000, `OP_TIMES` is 400, `OP_PLUS` is 310, `OP_PREFIX` (`@`) is 620, `OP_POSTFIX` (`//`) is 70).

### 2.4. Evaluator (`eval.c`)

The `evaluate(Expr* e)` function is the heart of Mathilda. It follows Mathematica's infinite evaluation semantics: expressions are repeatedly evaluated until a fixed point is reached (the expression no longer changes).

**Evaluation Sequence for Functions (`f[arg1, arg2]`):**
1.  **Evaluate Head**: The head `f` is evaluated first.
2.  **Check Attributes**: The attributes of the evaluated head are retrieved.
3.  **Evaluate Arguments**: Arguments are evaluated standardly, *unless* the head possesses Hold attributes (`ATTR_HOLDFIRST`, `ATTR_HOLDREST`, `ATTR_HOLDALL`, `ATTR_HOLDALLCOMPLETE`).
4.  **Apply Listable**: If the head is `ATTR_LISTABLE` and any evaluated argument is a `List`, the function automatically threads over the lists (e.g., `Mod[{1, 2}, 3] -> {Mod[1, 3], Mod[2, 3]}`).
5.  **Apply Flat / Orderless**: If `ATTR_FLAT` (associative), nested calls to the same head are flattened. If `ATTR_ORDERLESS` (commutative), arguments are lexically sorted.
6.  **Apply Built-ins**: If a C-level built-in function is registered, it is called.
7.  **Apply User Rules**: If no built-in handles it (or returns `NULL`), the evaluator checks the `DownValues` of the head and applies the first matching pattern replacement.

### 2.5. Pattern Matcher (`match.c`)

Implements structural pattern matching (`MatchQ`).
*   `Blank[]` (`_`): Matches any single expression.
*   `BlankSequence[]` (`__`): Matches 1 or more expressions in a sequence.
*   `BlankNullSequence[]` (`___`): Matches 0 or more expressions in a sequence.
*   `Pattern[name, pattern]` (`name_`): Binds the matched sub-expression to `name` inside a `MatchEnv`.
*   The matcher recursively unifies trees and handles sequence segmenting through backtracking.

### 2.6. Rule Engine (`replace.c`)

Implements expression transformations via rules (`Rule` `->` and `RuleDelayed` `:>`).
*   `ReplaceAll` (`/.`): Traverses the tree top-down, applying rules to sub-expressions.
*   Uses `match.c` to determine if a rule's LHS matches the current expression, and if so, substitutes bindings into the RHS.

## 3. Built-in Modules

Built-in functionality is heavily modularized:
*   `core.c`: Basic types and utilities (`AtomQ`, `Mod`, `Quotient`, `Re`, `Im`).
*   `list.c`: List generation and manipulation (`Table`, `Range`, `Array`, `Take`, `Drop`).
*   `part.c`: Structural extraction and modification (`Part` `[[ ]]`, `Head`, `First`, `Insert`, `Delete`).
*   `arithmetic.c` / `plus.c` / `times.c` / `power.c`: Core algebra, simplification, and complex number mechanics.
*   `boolean.c` / `comparisons.c`: Logical operators (`And`, `Or`, `Not`) and inequalities (`Equal`, `Less`, `SameQ`).

## 4. Development Workflow & Adding Features

### Adding a new Built-in Function
1.  **C Implementation**: Write the evaluation logic in an appropriate `.c` file (e.g., `core.c`). The signature must be `Expr* builtin_myfunc(Expr* res)`. 
    *   *Return `NULL`* if the function cannot evaluate the given arguments.
    *   *Return a new `Expr*`* if successful, and remember to `expr_free(res)` (or reuse its components) if you consume it.
2.  **Registration**: Register the function in the corresponding init function (e.g., `core_init()`) using `symtab_add_builtin("MyFunc", builtin_myfunc);`.
3.  **Attributes**: Assign necessary attributes in the init function via `symtab_get_def("MyFunc")->attributes |= ...`.
    *   If it threads over lists, add `ATTR_LISTABLE`.
    *   If it operates strictly on numbers, add `ATTR_NUMERICFUNCTION`.
    *   Almost all built-ins should be `ATTR_PROTECTED`.
    *   Don't forget to also add it to `builtin_attrs` array in `eval.c` so the attributes apply dynamically during recursive evaluation.
4.  **Testing**: Add test cases to the appropriate suite in `tests/` using the `TEST(test_my_feature)` macro block. Compile and run with `make -C tests && for t in tests/*_tests; do ./$t; done`.
5.  **Documentation**: Update `Mathilda_spec.md` with the new function, describing its behavior and providing sample inputs/outputs.

### Testing Conventions
All tests use a centralized macro suite defined in `tests/test_utils.h`. 
*   Use `ASSERT(condition)`.
*   Use `ASSERT_STR_EQ(a, b)`.
*   Register tests inside `main()` using `TEST(test_function_name);`.

## 5. Coding Standards
*   **C Standard**: C99.
*   **Memory Safety**: Always be vigilant about memory ownership. Unfreed intermediate evaluations represent leaks. Avoid double-frees by carefully tracing `expr_copy()` usage.
*   **Documentation**: Keep `Mathilda_spec.md` perfectly in sync with the current system capabilities.
