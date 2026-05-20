# Mathilda System Architecture Specification

## 1. Overview

Mathilda is a symbolic computer algebra system (CAS) written in C99, closely modeled on the core architecture and evaluation semantics of Mathematica (the Wolfram Language). It implements a recursive expression model, structural pattern matching with backtracking, rewriting rules, and an extensive library of built-in mathematical functions.

The name "Mathilda" pays homage to David Stoutemyer's seminal PICOMATH-80 tiny computer algebra system.

**Key characteristics:**
- Written in ~22,000 lines of C99 across 44 source modules
- Arbitrary-precision integer arithmetic via GMP
- Interactive REPL with GNU Readline support
- Mathematica-compatible syntax and evaluation semantics
- Licensed under GPLv3

**External dependencies:**
- **GMP** (GNU Multiple Precision Arithmetic Library) -- arbitrary-precision integers
- **GNU Readline** -- interactive command-line editing and history
- **GMP-ECM** (Elliptic Curve Method) -- advanced integer factorization (bundled in `src/external/ecm/`)

---

## 2. Project Structure

```
Mathilda/
├── src/                        # All source code (44 .c files, 44 .h files)
│   ├── expr.c / expr.h         # Expression representation (tagged union AST)
│   ├── parse.c / parse.h       # Pratt parser
│   ├── eval.c / eval.h         # Core evaluator (infinite evaluation loop)
│   ├── symtab.c / symtab.h     # Symbol table (values, rules, builtins)
│   ├── match.c / match.h       # Pattern matcher with backtracking
│   ├── replace.c / replace.h   # Rule engine (ReplaceAll, Replace)
│   ├── attr.c / attr.h         # Attribute system (bitflags)
│   ├── repl.c                  # REPL and main()
│   ├── print.c / print.h       # Expression output formatting
│   ├── core.c / core.h         # Core builtins + module initialization hub
│   ├── plus.c / times.c / power.c  # Specialized arithmetic heads
│   ├── arithmetic.c            # Basic numeric operations
│   ├── list.c                  # List generation and manipulation
│   ├── part.c                  # Part, Head, First, Insert, Delete
│   ├── poly/                   # Polynomials (see "Polynomials" in §7)
│   │   ├── poly.{c,h}          # PolynomialQ, GCD/LCM, Quotient/Remainder,
│   │   │                       #   Coefficient(List), Collect, Decompose,
│   │   │                       #   HornerForm, Resultant, Discriminant
│   │   ├── facpoly.{c,h}       # Factor, FactorSquareFree, FactorTerms[List]
│   │   ├── qafactor.{c,h}      # Trager algebraic-number factorisation
│   │   ├── qa.{c,h}            # Q(α) field elements
│   │   ├── qaupoly.{c,h}       # Q(α)[x] univariate polynomials
│   │   ├── mvfactor.{c,h}      # Multivariate factor orchestration
│   │   ├── mvfactor3.{c,h}     # Multifactor Hensel lift (3+ factors)
│   │   ├── bpoly.{c,h}         # Bivariate Z[x, y] polynomials
│   │   ├── mpoly.{c,h}         # Multivariate Z[x1,...,xn] polynomials
│   │   ├── zupoly.{c,h}        # Univariate Z[x] with GMP coefficients
│   │   └── solvepoly.{c,h}     # Solve`SolvePolynomialEquality
│   ├── facint.c                # Integer factorization (multiple algorithms)
│   ├── trig.c                  # Trigonometric functions
│   ├── hyperbolic.c            # Hyperbolic functions
│   ├── logexp.c                # Log, Exp
│   ├── complex.c               # Complex number arithmetic
│   ├── linalg/                 # Linear algebra (see "Linear Algebra" in §7)
│   │   ├── linalg.{c,h}        # Module init (registers all builtins)
│   │   ├── util.c              # Shared tensor helpers (dims, flatten, exact div)
│   │   ├── dot.c, det.c, cross.c, norm.c, tr.c   # Per-builtin
│   │   ├── construct.c, matpow.c                 # IdentityMatrix/DiagonalMatrix, MatrixPower
│   │   ├── inv.{c,h}           # Inverse, PseudoInverse
│   │   ├── linsolve.{c,h}      # LinearSolve, RowReduce
│   │   ├── lstsq.{c,h}         # LeastSquares
│   │   ├── eigen.{c,h}         # Eigenvalues/Eigenvectors entry + dispatch
│   │   ├── eigen_internal.h    # Shared MatD/MatM/EigenOpts + dispatcher decls
│   │   ├── eigen_common.c      # Shared eigen helpers, MPFR array ops
│   │   ├── eigen_direct.c      # Direct kernel (Hessenberg+QR), machine + MPFR
│   │   ├── eigen_arnoldi.c     # Arnoldi (Krylov) kernel, machine + MPFR
│   │   ├── eigen_banded.c      # Banded kernel, machine + MPFR
│   │   └── eigen_feast.c       # FEAST (contour integral) kernel, machine + MPFR
│   ├── funcprog.c              # Functional programming (Map, Apply, Select, etc.)
│   ├── purefunc.c              # Pure functions (Function, Slot, SlotSequence)
│   ├── boolean.c               # And, Or, Not
│   ├── comparisons.c           # Equal, Less, SameQ, etc.
│   ├── cond.c                  # If, Which
│   ├── iter.c                  # Do, For, While loops
│   ├── patterns.c              # Cases, Position, Count, MatchQ
│   ├── stats.c                 # Mean, Variance, Median, etc.
│   ├── sort.c                  # Sort, OrderedQ
│   ├── rat.c                   # Rational function operations
│   ├── parfrac.c               # Partial fractions (Apart)
│   ├── expand.c                # Expand, Collect
│   ├── simp/                   # Simplification (see "Simplification" in §7)
│   │   ├── simp.{c,h}          # Simplify, SimplifyCount, Assuming, Element
│   │   ├── trigsimp.{c,h}      # TrigToExp, ExpToTrig, TrigExpand, TrigFactor, TrigReduce
│   │   └── trigrat.{c,h}       # Trig rationalisation helpers
│   ├── calculus/               # Calculus (see "Calculus" in §7)
│   │   ├── deriv.{c,h}         # D, Dt, Derivative
│   │   ├── series.{c,h}        # Series, Normal
│   │   ├── limit.{c,h}         # Limit
│   │   ├── integrate.{c,h}     # Integrate (dispatcher)
│   │   ├── intrat.{c,h}        # Rational integration (Bronstein pipeline)
│   │   ├── intrat_internal.h   # Shared private surface between intrat and intsimp
│   │   ├── intrischnorman.{c,h}# Risch-Norman algorithm
│   │   └── intsimp.{c,h}       # Integration simplification helpers
│   ├── modular.c               # PowerMod, modular arithmetic
│   ├── piecewise.c             # Piecewise conditional expressions
│   ├── load.c                  # File I/O (Get)
│   ├── info.c                  # Information, ? queries
│   ├── datetime.c              # Date/time functions
│   ├── default_helper.c        # Default value handling
│   ├── internal.c              # Internal system functions
│   ├── internal/               # Mathematica-syntax initialization scripts
│   │   ├── init.m              # Bootstrap (loads other .m files)
│   │   ├── deriv.m             # Derivative rules (D[...])
│   │   └── CRCMathTablesIntegrals.m  # Reference integral tables
│   └── external/
│       └── ecm/                # Bundled GMP-ECM library (DO NOT MODIFY)
├── tests/                      # Unit test suite
│   ├── CMakeLists.txt          # CMake build for tests
│   ├── test_utils.h            # Shared test macros (TEST, ASSERT, ASSERT_STR_EQ)
│   └── test_*.c                # ~30+ test source files
├── makefile                    # Primary build system
├── Mathilda_spec.md             # Comprehensive function reference
├── CLAUDE.md                   # AI development guidelines
└── README.md                   # User guide
```

---

## 3. Core Architecture

Mathilda is built from seven interdependent subsystems that together implement a Mathematica-like evaluation pipeline. The flow is:

```
Input String
    │
    ▼
┌─────────┐
│  Parser  │  parse.c -- Pratt parser, produces Expr* AST
└────┬─────┘
     │
     ▼
┌──────────┐
│ Evaluator│  eval.c -- Infinite evaluation loop to fixed point
└────┬─────┘
     │  Uses:
     ├──► Symbol Table (symtab.c) -- OwnValues, DownValues, builtins
     ├──► Attribute System (attr.c) -- Controls evaluation behavior
     ├──► Pattern Matcher (match.c) -- Structural matching + backtracking
     ├──► Rule Engine (replace.c) -- ReplaceAll, user-defined rules
     └──► Built-in Libraries -- C functions for mathematical operations
            │
            ▼
        Output Expr*
            │
            ▼
┌──────────┐
│ Printer  │  print.c -- Formats Expr* for display
└──────────┘
```

### 3.1 Expression Representation (`expr.c`, `expr.h`)

Everything in Mathilda is an `Expr`, implemented as a tagged union:

```c
typedef enum {
    EXPR_INTEGER,   // int64_t
    EXPR_REAL,      // double
    EXPR_SYMBOL,    // char* (e.g., "x", "Plus", "Pi")
    EXPR_STRING,    // char*
    EXPR_FUNCTION,  // head (Expr*) + args (Expr**) + arg_count
    EXPR_BIGINT     // mpz_t (GMP arbitrary-precision integer)
} ExprType;

typedef struct Expr {
    ExprType type;
    union {
        int64_t integer;
        double real;
        char* symbol;
        char* string;
        struct {
            struct Expr* head;
            struct Expr** args;
            size_t arg_count;
        } function;
        mpz_t bigint;
    } data;
} Expr;
```

**Design principles:**
- Expressions are **immutable by convention** during evaluation. Functions that transform expressions return new allocated trees.
- The `EXPR_FUNCTION` type is the universal compound expression. A call like `f[x, y]` is represented as a function node with head=`Symbol("f")` and args=`[Symbol("x"), Symbol("y")]`. Lists `{a, b}` are `List[a, b]` internally.
- `EXPR_BIGINT` automatically promotes from `EXPR_INTEGER` when arithmetic overflows 64-bit limits. The helper `expr_bigint_normalize()` demotes back to `EXPR_INTEGER` when the value fits.

**Key API:**
| Function | Purpose |
|----------|---------|
| `expr_new_integer(val)` | Allocate an integer node |
| `expr_new_real(val)` | Allocate a real node |
| `expr_new_symbol(name)` | Allocate a symbol node (copies the string) |
| `expr_new_string(str)` | Allocate a string node |
| `expr_new_function(head, args, count)` | Allocate a compound expression |
| `expr_new_bigint_from_mpz(val)` | Allocate a bigint node |
| `expr_copy(e)` | Deep copy an expression tree |
| `expr_free(e)` | Deep deallocation of an expression tree |
| `expr_eq(a, b)` | Structural equality test |
| `expr_compare(a, b)` | Canonical ordering comparison |
| `expr_hash(e)` | Hash for use in hash tables |
| `expr_to_mpz(e, out)` | Convert integer-like Expr to GMP mpz_t |
| `expr_is_integer_like(e)` | True for EXPR_INTEGER or EXPR_BIGINT |

### 3.2 Parser (`parse.c`, `parse.h`)

The parser translates input strings into `Expr*` AST trees using a **Pratt parser** (top-down operator precedence).

**Key characteristics:**
- Lexical analysis is handled inline by advancing a `ParserState` pointer and skipping whitespace. There is no separate tokenizer pass.
- Operator precedences exactly mirror Mathematica's standard precedences.
- Supports infix, prefix, and postfix operators.
- Supports Mathematica pattern syntax: `_`, `__`, `___`, `x_`, `x_h`, `/;` (Condition).

**Operator precedence table (selected):**

| Precedence | Operators |
|------------|-----------|
| 1000 | `f[x]` (function call) |
| 730 | `_`, `__`, `___` (blanks) |
| 680 | `?` (PatternTest) |
| 620 | `@` (Prefix), `@@` (Apply), `/@` (Map) |
| 590 | `^` (Power, right-associative) |
| 470 | `/` (Divide) |
| 400 | `*` (Times) |
| 310 | `+`, `-` (Plus) |
| 290 | `==`, `===`, `<`, `>`, etc. |
| 100 | `[[]]` (Part) |
| 90 | `&` (Function) |
| 70 | `//` (Postfix) |
| 40 | `=`, `:=`, `->`, `:>` |
| 10 | `;` (CompoundExpression) |

**Public API:**
- `parse_expression(input)` -- Parse a complete string, return `Expr*` or `NULL` on error.
- `parse_next_expression(input_ptr)` -- Parse and advance pointer, for multi-expression input.

### 3.3 Symbol Table (`symtab.c`, `symtab.h`)

The symbol table stores global state for all symbols. Each symbol is represented by a `SymbolDef`:

```c
typedef struct {
    char* symbol_name;
    Rule* own_values;       // OwnValues: immediate assignments (x = 5)
    Rule* down_values;      // DownValues: pattern rules (f[x_] := x^2)
    BuiltinFunc builtin_func; // C function pointer, or NULL
    uint32_t attributes;    // Bitflag attributes
    char* docstring;        // Help text for ?symbol
} SymbolDef;

typedef Expr* (*BuiltinFunc)(Expr* res);

typedef struct Rule {
    Expr* pattern;
    Expr* replacement;
    struct Rule* next;      // Linked list of rules
} Rule;
```

**How values are stored:**
- **OwnValues** are immediate symbol assignments. `x = 5` stores a rule matching the symbol `x` to the value `5`.
- **DownValues** are pattern-based rules. `f[x_] := x^2` stores a rule on `f` that pattern-matches any single-argument call.
- Rules are stored as linked lists (newest first for DownValues, which gives "most specific rule first" semantics when rules are added in order).

**Key API:**
| Function | Purpose |
|----------|---------|
| `symtab_init()` | Initialize the global symbol table |
| `symtab_clear()` | Free the entire symbol table |
| `symtab_get_def(name)` | Get or create a SymbolDef for a symbol |
| `symtab_add_builtin(name, func)` | Register a C built-in function |
| `symtab_set_docstring(name, doc)` | Set the `?name` help text |
| `symtab_add_own_value(name, pattern, replacement)` | Add an OwnValue rule |
| `symtab_add_down_value(name, pattern, replacement)` | Add a DownValue rule |
| `apply_own_values(expr)` | Try OwnValue rules on a symbol; return new Expr* or NULL |
| `apply_down_values(expr)` | Try DownValue rules on a function; return new Expr* or NULL |

### 3.4 Attribute System (`attr.c`, `attr.h`)

Attributes are bitflags on symbols that control how the evaluator processes expressions with that symbol as a head:

```c
#define ATTR_HOLDFIRST       (1 << 0)   // Don't evaluate the first argument
#define ATTR_HOLDREST        (1 << 1)   // Don't evaluate arguments after the first
#define ATTR_HOLDALL         (ATTR_HOLDFIRST | ATTR_HOLDREST)
#define ATTR_HOLDALLCOMPLETE (1 << 2)   // Don't evaluate any args, don't flatten Sequence
#define ATTR_FLAT            (1 << 3)   // Associative: flatten nested same-head calls
#define ATTR_ORDERLESS       (1 << 4)   // Commutative: sort arguments canonically
#define ATTR_LISTABLE        (1 << 5)   // Auto-thread over List arguments
#define ATTR_NUMERICFUNCTION (1 << 6)   // Numeric operation
#define ATTR_PROTECTED       (1 << 7)   // Cannot be redefined by user
#define ATTR_ONEIDENTITY     (1 << 8)   // f[x] equivalent to x for pattern matching
#define ATTR_NHOLDREST       (1 << 9)   // Numeric hold on rest arguments
#define ATTR_LOCKED          (1 << 10)  // Attributes themselves cannot be changed
#define ATTR_READPROTECTED   (1 << 11)  // Definition not readable
#define ATTR_TEMPORARY       (1 << 12)  // Created by Module, auto-cleaned
#define ATTR_SEQUENCEHOLD    (1 << 13)  // Don't flatten Sequence in arguments
```

**Common attribute combinations:**
- `Plus`: `ATTR_FLAT | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS`
- `Times`: `ATTR_FLAT | ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_ONEIDENTITY | ATTR_ORDERLESS`
- `Table`: `ATTR_HOLDALL | ATTR_PROTECTED`
- `If`: `ATTR_HOLDREST | ATTR_PROTECTED`
- `Sin`: `ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED`

### 3.5 Evaluator (`eval.c`, `eval.h`)

The evaluator is the heart of Mathilda. It implements Mathematica's **infinite evaluation semantics**: expressions are repeatedly evaluated until a fixed point is reached (the expression no longer changes).

**`evaluate(Expr* e)` -- the main evaluation loop:**

The outer loop calls `evaluate_step()` repeatedly, comparing each result to the previous one via `expr_eq()`. It terminates when the expression stabilizes or a recursion limit is hit.

**`evaluate_step(Expr* e)` -- single evaluation pass:**

For a function expression `f[arg1, arg2, ...]`:

1. **Evaluate Head** -- Evaluate `f` first (it might itself be a computed expression).
2. **Check Attributes** -- Retrieve the attributes of the evaluated head.
3. **Evaluate Arguments** -- Evaluate each argument, *unless* the head has Hold attributes:
   - `ATTR_HOLDFIRST`: Don't evaluate arg 1.
   - `ATTR_HOLDREST`: Don't evaluate args 2..N.
   - `ATTR_HOLDALL`: Don't evaluate any args.
   - `ATTR_HOLDALLCOMPLETE`: Don't evaluate any args and don't flatten `Sequence`.
4. **Apply Listable** -- If `ATTR_LISTABLE` and any argument is a `List`, thread the function element-wise. E.g., `Sin[{1, 2, 3}]` becomes `{Sin[1], Sin[2], Sin[3]}`.
5. **Apply Flat** -- If `ATTR_FLAT`, flatten nested calls: `Plus[a, Plus[b, c]]` becomes `Plus[a, b, c]`.
6. **Apply Orderless** -- If `ATTR_ORDERLESS`, sort arguments into canonical order.
7. **Apply Built-in** -- If a C function is registered for this head, call it. The built-in returns a new `Expr*` (success) or `NULL` (cannot evaluate, remain unevaluated).
8. **Apply User Rules** -- If no built-in handled it, try DownValues pattern matching. The first matching rule's replacement is substituted.

For a symbol expression `x`:
- Try OwnValues. If a rule matches, return the replacement.

**Key API:**
| Function | Purpose |
|----------|---------|
| `evaluate(e)` | Full evaluation to fixed point |
| `evaluate_step(e)` | Single evaluation step |
| `eval_flatten_args(e, head_name)` | Flatten nested same-head calls |
| `eval_sort_args(e)` | Sort arguments into canonical order |

### 3.6 Pattern Matcher (`match.c`, `match.h`)

The pattern matcher implements structural tree unification with sequence-matching and backtracking.

**Pattern forms:**
| Pattern | Internal Form | Matches |
|---------|--------------|---------|
| `_` | `Blank[]` | Any single expression |
| `_h` | `Blank[h]` | Any expression with head `h` |
| `__` | `BlankSequence[]` | 1 or more expressions |
| `___` | `BlankNullSequence[]` | 0 or more expressions |
| `x_` | `Pattern[x, Blank[]]` | Any single expression, bound to `x` |
| `x__` | `Pattern[x, BlankSequence[]]` | 1+ expressions, bound to `x` |
| `p ? test` | `PatternTest[p, test]` | `p` if `test[matched]` is `True` |
| `p /; cond` | `Condition[p, cond]` | `p` if `cond` evaluates to `True` |
| `p..` | `Repeated[p]` | Sequence of 1+ expressions each matching `p` |
| `p...` | `RepeatedNull[p]` | Sequence of 0+ expressions each matching `p` |
| `p:def` | `Optional[p, def]` | `p`, or `def` if argument is absent |

**Match environment (`MatchEnv`):**

Pattern bindings are stored in a `MatchEnv` structure:

```c
typedef struct MatchEnv {
    char** symbols;     // Bound variable names
    Expr** values;      // Corresponding matched expressions
    size_t count;
    size_t capacity;
    bool (*callback)(struct MatchEnv*, void*);  // Optional match callback
    void* callback_data;
} MatchEnv;
```

**Key API:**
| Function | Purpose |
|----------|---------|
| `match(expr, pattern, env)` | Match expression against pattern, populate bindings |
| `replace_bindings(expr, env)` | Substitute bindings into an expression template |
| `env_new()` | Create empty match environment |
| `env_free(env)` | Free match environment |
| `env_set(env, name, value)` | Add a binding |
| `env_get(env, name)` | Retrieve a binding |

**Backtracking:** The matcher handles sequence patterns (`__`, `___`) through backtracking. When matching function arguments against a pattern list that contains sequence blanks, the matcher tries different partitions of the argument list, backtracking when a partition fails to match subsequent pattern elements.

### 3.7 Rule Engine (`replace.c`, `replace.h`)

The rule engine applies transformation rules to expressions:

- **`Rule` (`->`):** Immediate rules. The RHS is evaluated at rule creation time.
- **`RuleDelayed` (`:>`):** Delayed rules. The RHS is evaluated only after pattern variables are substituted.
- **`ReplaceAll` (`/.`):** Top-down traversal applying rules to each subpart. Once a rule matches a subpart, its children are not further traversed.
- **`ReplaceRepeated` (`//.`):** Applies `ReplaceAll` repeatedly until the expression stops changing.
- **`Replace`:** Applies rules only at specified levels.

---

## 4. Memory Management

Mathilda uses **explicit manual memory management**. This is the most critical aspect of the codebase for contributors to understand.

### 4.1 Ownership Rules

**The Fundamental Contract for Built-in Functions:**

```
Expr* builtin_myfunc(Expr* res)
```

- The built-in function **takes ownership** of `res`.
- **If it returns a new `Expr*`:** The built-in must `expr_free(res)` (or structurally reuse its nodes). The caller (evaluator) takes ownership of the returned value.
- **If it returns `NULL`:** The built-in could not evaluate the expression. The evaluator retains ownership of `res` and it remains unevaluated. The built-in must NOT free `res` in this case.

### 4.2 Common Patterns

**Returning early when arguments don't match:**
```c
Expr* builtin_myfunc(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;  // Don't free!
    Expr* arg0 = res->data.function.args[0];
    if (arg0->type != EXPR_INTEGER) return NULL;         // Don't free!
    // ... process ...
    Expr* result = expr_new_integer(42);
    expr_free(res);  // Free input since we're returning a new result
    return result;
}
```

**Reusing parts of the input tree:**
```c
// Steal arg0 from res before freeing res
Expr* arg0 = res->data.function.args[0];
res->data.function.args[0] = NULL;  // Prevent double-free
expr_free(res);
// arg0 is now owned by us, use it in the result
return arg0;
```

**Deep copying when you need to keep a reference:**
```c
Expr* copy = expr_copy(res->data.function.args[0]);
// ... use copy independently of res ...
```

### 4.3 Memory Safety Checklist

- Never access an `Expr*` after calling `expr_free()` on it.
- Never free the same `Expr*` twice. Use the NULL-out-before-free pattern above.
- Always trace ownership: if you `expr_copy()` something, you must eventually `expr_free()` the copy.
- Use `valgrind` to detect leaks: `valgrind --leak-check=full ./Mathilda`

---

## 5. Module Initialization Chain

At startup, `main()` in `repl.c` calls:

```c
symtab_init();   // Initialize the global symbol table
core_init();     // Register all built-in functions
```

`core_init()` in `core.c` serves as the initialization hub, registering its own builtins and then calling every other module's init function:

```
core_init()
├── [registers core builtins: AtomQ, Mod, Quotient, Re, Im, Abs, ...]
├── parfrac_init()
├── modular_init()
├── facint_init()
├── comparisons_init()
├── boolean_init()
├── list_init()
├── replace_init()
├── patterns_init()
├── cond_init()
├── iter_init()
├── complex_init()
├── trig_init()
├── simp_init()
├── hyperbolic_init()
├── logexp_init()
├── piecewise_init()
├── attr_init()
├── purefunc_init()
├── stats_init()
├── poly_init()
├── facpoly_init()
├── rat_init()
├── expand_init()
├── info_init()
├── datetime_init()
├── linalg_init()
└── load_init()
```

After C initialization, `main()` loads `src/internal/init.m`, which bootstraps additional rules defined in Mathematica syntax (e.g., derivative rules in `deriv.m`).

---

## 6. Internal Definition Files

The `src/internal/` directory contains Mathematica-syntax definition files that are loaded at startup:

- **`init.m`** -- Bootstrap script, loads other `.m` files via `Get[]`.
- **`deriv.m`** -- Derivative rules for `D[expr, x]`. This is a pure pattern-matching implementation: rules like `D[Sin[f_], x_] := Cos[f] * D[f, x]` are defined as DownValues using Mathilda's own syntax, demonstrating the power of the rule system. Covers all elementary functions, chain rule, product rule, and higher-order derivatives.
- **`CRCMathTablesIntegrals.m`** -- Reference integral tables from the CRC Mathematical Tables.

This architecture demonstrates a key design pattern: core operations are implemented in C for performance, while higher-level mathematical rules can be defined in the system's own language via DownValues.

---

## 7. Built-in Module Reference

Each module follows the same pattern: a `.c`/`.h` pair, with a `*_init()` function that registers builtins and assigns attributes.

| Module | File | Key Functions |
|--------|------|---------------|
| **Core** | `core.c` | `AtomQ`, `NumberQ`, `IntegerQ`, `Mod`, `Quotient`, `GCD`, `LCM`, `Re`, `Im`, `Abs`, `Conjugate`, `Arg`, `Factorial`, `Binomial`, `PrimeQ`, `PrimePi`, `NextPrime`, `EulerPhi`, `Numerator`, `Denominator`, `Floor`, `Ceiling`, `Round` |
| **Arithmetic** | `arithmetic.c` | Low-level numeric operations, overflow detection, bigint promotion |
| **Plus** | `plus.c` | Symbolic addition: term collection, constant folding, `Overflow[]` handling |
| **Times** | `times.c` | Symbolic multiplication: coefficient merging, base grouping |
| **Power** | `power.c` | Exponentiation: radical simplification, integer power distribution |
| **Lists** | `list.c` | `Table`, `Range`, `Array`, `Take`, `Drop`, `Flatten`, `Partition`, `Split`, `Union`, `DeleteDuplicates`, `Tally`, `Commonest`, `Append`, `Prepend`, `Reverse`, `RotateLeft`, `RotateRight`, `Total`, `Min`, `Max`, `Tuples`, `Permutations` |
| **Parts** | `part.c` | `Part` (`[[]]`), `Head`, `Length`, `Dimensions`, `First`, `Last`, `Most`, `Rest`, `Insert`, `Delete`, `Depth`, `Level`, `Extract`, `Span`, `LeafCount`, `ByteCount` |
| **Pattern Ops** | `patterns.c` | `Cases`, `Position`, `Count`, `MemberQ`, `FreeQ` |
| **Matching** | `match.c` | `MatchQ`, `Blank`, `BlankSequence`, `BlankNullSequence`, `Pattern`, `PatternTest`, `Condition`, `Repeated`, `RepeatedNull`, `Optional`, `Default` |
| **Comparisons** | `comparisons.c` | `Equal`, `Unequal`, `Less`, `Greater`, `LessEqual`, `GreaterEqual`, `SameQ`, `UnsameQ`, `OrderedQ` |
| **Boolean** | `boolean.c` | `And`, `Or`, `Not`, `Xor`, `TrueQ`, `BooleanQ` |
| **Conditions** | `cond.c` | `If`, `Which`, `Switch`, `Piecewise` |
| **Iteration** | `iter.c` | `Do`, `For`, `While`, `NestWhile` |
| **Functional** | `funcprog.c` | `Map`, `MapAll`, `Apply`, `Select`, `Through`, `Distribute` |
| **Pure Functions** | `purefunc.c` | `Function`, `Slot` (`#`), `SlotSequence` (`##`), `CompoundExpression` |
| **Complex** | `complex.c` | Complex arithmetic, `Complex[re, im]` construction |
| **Trig** | `trig.c` | `Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`, `ArcSin`, `ArcCos`, `ArcTan` (including 2-arg form), etc. Exact values for rational multiples of Pi |
| **Hyperbolic** | `hyperbolic.c` | `Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` and inverses |
| **Log/Exp** | `logexp.c` | `Log`, `Log[b, z]`, `Exp`. Euler's formula for `Exp[I*q*Pi]` |
| **Simplification** | `simp/` | `Simplify`, `SimplifyCount`, `Assuming`, `Element`, `TrigToExp`, `ExpToTrig`, `TrigExpand`, `TrigFactor`, `TrigReduce`. Simplify pipeline (`simp/simp.c`) drives candidate-set search over a transform table; trig-specific rewrites live in `simp/trigsimp.c` and `simp/trigrat.c`. |
| **Expansion** | `expand.c` | `Expand`, `Collect` |
| **Polynomials** | `poly/` | `Coefficient`, `CoefficientList`, `PolynomialGCD`, `PolynomialExtendedGCD`, `PolynomialLCM`, `PolynomialQuotient`, `PolynomialRemainder`, `PolynomialQ`, `PolynomialMod`, `Resultant`, `Discriminant`, `Variables`, `HornerForm`, `Collect`, `Decompose`, `SubresultantPolynomialRemainders`. See `src/poly/poly.{c,h}` for the public surface; helper substrates (`zupoly`, `bpoly`, `mpoly`, `qa`, `qaupoly`) live alongside. |
| **Poly Factoring** | `poly/facpoly.c` | `Factor` (over the integers), `FactorSquareFree`, `FactorTerms`, `FactorTermsList`. Bivariate / multivariate Hensel pipeline in `poly/{bpoly,mvfactor,mvfactor3}.c`; algebraic-number Trager pipeline in `poly/qafactor.c`. |
| **Poly Solving** | `poly/solvepoly.c` | `` Solve`SolvePolynomialEquality `` (single-variable polynomial-equality solver backing the public `Solve` router). |
| **Int Factoring** | `facint.c` | `FactorInteger` with methods: Trial Division, Pollard's Rho, Pollard P-1, Williams P+1, Fermat, CFRAC, SQUFOF, ECM |
| **Rational Funcs** | `rat.c` | `Cancel`, `Together` |
| **Partial Fractions** | `parfrac.c` | `Apart` |
| **Modular** | `modular.c` | `PowerMod` (exponentiation, inverse, roots) |
| **Calculus** | `calculus/` | `D`, `Dt`, `Derivative` (in `calculus/deriv.c`); `Series`, `Normal` (in `calculus/series.c`); `Limit` (in `calculus/limit.c`); `Integrate` dispatcher (in `calculus/integrate.c`) routing to the rational-integration Bronstein pipeline (`calculus/intrat.c`) and the Risch-Norman algorithm (`calculus/intrischnorman.c`); integration simplification helpers in `calculus/intsimp.c`. |
| **Linear Algebra** | `linalg/` | `Dot`, `Inner`, `Outer`, `Det`, `Inverse`, `PseudoInverse`, `MatrixPower`, `Cross`, `Norm`, `Tr`, `RowReduce`, `LinearSolve`, `LeastSquares`, `IdentityMatrix`, `DiagonalMatrix`, `Transpose`, `Eigenvalues`, `Eigenvectors`. Split into per-builtin files; `eigen.c` dispatches the Direct / Arnoldi / Banded / FEAST kernels (machine + MPFR) in `eigen_*.c`. |
| **Statistics** | `stats.c` | `Mean`, `Median`, `Variance`, `StandardDeviation`, `RootMeanSquare`, `Quartiles` |
| **Sorting** | `sort.c` | `Sort`, `OrderedQ` |
| **Rules/Replace** | `replace.c` | `Replace`, `ReplaceAll`, `ReplaceRepeated`, `ReplacePart`, `ReplaceList` |
| **Piecewise** | `piecewise.c` | `Piecewise` |
| **Information** | `info.c` | `Information`, `?` shortcut |
| **File I/O** | `load.c` | `Get` (loads and evaluates .m files) |
| **DateTime** | `datetime.c` | `Timing`, `RepeatedTiming`, date functions |
| **Scoping** | (in `eval.c`/`purefunc.c`) | `Module`, `Block`, `With` |
| **Printing** | `print.c` | `Print`, `FullForm`, `InputForm`, `HoldForm` |

---

## 8. How to Extend Mathilda

### 8.1 Adding a New Built-in Function

This is the most common extension. Follow these steps:

#### Step 1: Write the C Implementation

Choose the appropriate module file (or create a new one if the function belongs to a new mathematical domain). The function signature must be:

```c
Expr* builtin_myfunc(Expr* res)
```

**Template:**

```c
Expr* builtin_myfunc(Expr* res) {
    // 1. Validate arguments
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return NULL;  // Return NULL = can't evaluate

    Expr* arg = res->data.function.args[0];

    // 2. Check argument types
    if (arg->type != EXPR_INTEGER) return NULL;

    // 3. Compute the result
    int64_t val = arg->data.integer;
    int64_t result_val = val * val;  // Example: squaring

    // 4. Build the result expression
    Expr* result = expr_new_integer(result_val);

    // 5. Free the input (we took ownership)
    expr_free(res);

    return result;
}
```

**Important guidelines:**
- Return `NULL` if the function cannot evaluate (wrong number of args, symbolic args to a numeric-only function, etc.). Do NOT free `res` in this case.
- Return a new `Expr*` on success, and `expr_free(res)`.
- Handle `EXPR_BIGINT` alongside `EXPR_INTEGER` if your function should work on large integers.
- Handle `EXPR_REAL` if your function should work on floating-point numbers.
- For functions that should handle `Rational[n, d]` or `Complex[re, im]`, check for `EXPR_FUNCTION` with the appropriate head.

#### Step 2: Declare the Function

Add the declaration to the appropriate header file:

```c
// In myfunc_module.h
Expr* builtin_myfunc(Expr* res);
```

#### Step 3: Register in the Init Function

In the module's `*_init()` function:

```c
void mymodule_init(void) {
    // Register the function
    symtab_add_builtin("MyFunc", builtin_myfunc);

    // Set attributes
    symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;

    // Set docstring (accessible via ?MyFunc in the REPL)
    symtab_set_docstring("MyFunc",
        "MyFunc[x]\n\tComputes the square of x.");
}
```

If this is a new module, add the `mymodule_init()` call to `core_init()` in `core.c`.

#### Step 4: Assign Appropriate Attributes

Choose attributes based on the function's mathematical properties:

| Attribute | When to Use |
|-----------|-------------|
| `ATTR_PROTECTED` | Almost always. Prevents user from redefining the built-in. |
| `ATTR_LISTABLE` | If the function should auto-thread over lists: `f[{a, b}]` → `{f[a], f[b]}` |
| `ATTR_NUMERICFUNCTION` | If the function operates on numbers |
| `ATTR_FLAT` | If the function is associative: `f[f[a, b], c]` → `f[a, b, c]` |
| `ATTR_ORDERLESS` | If the function is commutative: `f[b, a]` → `f[a, b]` |
| `ATTR_HOLDALL` | If arguments should not be evaluated before the function runs |
| `ATTR_HOLDFIRST` | If only the first argument should be held unevaluated |
| `ATTR_HOLDREST` | If arguments after the first should be held |

#### Step 5: Write Tests

Add tests to the appropriate test file in `tests/`, or create a new one:

```c
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

void test_myfunc_basic() {
    // Method 1: Using the assert_eval_eq helper (preferred for simple cases)
    assert_eval_eq("MyFunc[5]", "25", 0);

    // Method 2: Manual construction (for detailed assertions)
    Expr* input = parse_expression("MyFunc[3]");
    Expr* result = evaluate(input);
    ASSERT(result->type == EXPR_INTEGER);
    ASSERT(result->data.integer == 9);
    expr_free(input);
    expr_free(result);
}

void test_myfunc_symbolic() {
    // Should remain unevaluated for symbolic input
    assert_eval_eq("MyFunc[x]", "MyFunc[x]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_myfunc_basic);
    TEST(test_myfunc_symbolic);

    printf("All tests passed!\n");
    return 0;
}
```

**Test macros available in `test_utils.h`:**
- `TEST(func_name)` -- Run a test and print its name
- `ASSERT(cond)` -- Assert a condition
- `ASSERT_STR_EQ(a, b)` -- Assert string equality
- `ASSERT_MSG(cond, fmt, ...)` -- Assert with formatted error message
- `assert_eval_eq(input, expected, is_fullform)` -- Parse, evaluate, and compare string output

If you created a new test file, add it to `tests/CMakeLists.txt`.

#### Step 6: Update Documentation

Update `Mathilda_spec.md` with the new function: describe its behavior, list its attributes, and provide sample inputs/outputs following the existing format.

### 8.2 Adding a New Module

If your functions represent a new mathematical domain:

1. Create `src/mymodule.c` and `src/mymodule.h`.
2. Define a `void mymodule_init(void)` function in the `.c` file.
3. Declare `void mymodule_init(void)` in the `.h` file.
4. `#include "mymodule.h"` in `core.c` and add `mymodule_init()` to `core_init()`.
5. The makefile automatically picks up new `.c` files in `src/` via the wildcard `SRC = $(wildcard $(SRC_DIR)/*.c)`.

### 8.3 Adding New Pattern Constructs

To add a new pattern form:

1. Modify `match.c` to recognize and handle the new pattern during tree unification.
2. If the pattern requires new syntax (e.g., a new operator), modify `parse.c` to parse it and produce the correct `Expr*` tree.
3. Add test cases in `tests/` covering match and non-match scenarios.
4. Ensure the evaluator in `eval.c` doesn't accidentally evaluate the pattern expression before the matcher sees it (you may need Hold attributes).

### 8.4 Adding Rules via Internal .m Files

For mathematical identities that are naturally expressed as rewrite rules, you can define them in Mathematica syntax:

1. Create a new `.m` file in `src/internal/`.
2. Add a `Get["src/internal/yourfile.m"]` line to `src/internal/init.m`.
3. Define rules using standard syntax:

```mathematica
(* Example: rules for a new function *)
MyFunc[0] := 0;
MyFunc[x_ + y_] := MyFunc[x] + MyFunc[y];
MyFunc[n_Integer x_] := n MyFunc[x];
```

This approach is ideal for derivative rules, identities, and simplification rules where the pattern-matching system does the heavy lifting.

### 8.5 Adding New Operators to the Parser

To add a new syntactic operator:

1. In `parse.c`, define a new precedence constant (consult the Mathematica documentation for the correct value).
2. Add the operator's token recognition in the lexer section of the parser.
3. Implement the infix/prefix/postfix parsing logic in `parse_expression_prec()`.
4. The operator should produce a standard `Expr*` function call (e.g., `a <> b` should produce `StringJoin[a, b]`).

---

## 9. Build System

### Building Mathilda

```bash
make -j$(nproc)     # Build the Mathilda executable
./Mathilda            # Run the REPL
```

The makefile uses:
- `gcc` with `-std=c99 -O3`
- Automatic source file discovery: `SRC = $(wildcard src/*.c)`
- Links against `-lreadline -lgmp -lm`
- Optionally builds and links the ECM library (`USE_ECM=1` by default)

### Building and Running Tests

```bash
cd tests
mkdir -p build && cd build
cmake ..
make -j$(nproc)

# Run all test binaries
for t in *_tests; do ./$t; done
```

### Checking for Memory Leaks

```bash
valgrind --leak-check=full ./Mathilda
```

---

## 10. Coding Standards

- **C Standard:** C99 strictly. No C11+ features.
- **Memory Safety:** Always trace ownership of `Expr*` pointers. Every `expr_new_*` or `expr_copy` must have a corresponding `expr_free`. Use valgrind regularly.
- **No modifications to `src/external/`:** The ECM library is a vendored dependency and must not be modified.
- **Docstrings:** Every built-in function must have a docstring set via `symtab_set_docstring()` so it is accessible via `?FunctionName` in the REPL.
- **Attributes:** Every built-in function must have appropriate attributes assigned in its module init function.
- **Documentation:** Keep `Mathilda_spec.md` in sync with the current system capabilities. Every new or modified function must be documented there.

---

## 11. REPL Architecture

The REPL (`repl.c`) implements a standard Read-Eval-Print Loop:

1. **Read:** Uses GNU Readline for line input with history. Supports multiline input via trailing `\`.
2. **Parse:** `parse_expression(input)` converts the string to an `Expr*` AST.
3. **Evaluate:** `evaluate(parsed)` runs the expression through the evaluation pipeline.
4. **Store:** Input is stored as `In[n]` and output as `Out[n]` via DownValues, accessible in later expressions.
5. **Print:** `expr_print(evaluated)` formats and displays the result.
6. **Clean up:** Both the parsed input and evaluated result are freed.

The REPL exits on `Quit[]` or EOF (Ctrl-D).

---

## 12. Expression Printing (`print.c`)

The printer converts `Expr*` trees back into human-readable Mathematica-like syntax. Key features:

- **Infix notation:** `Plus[a, b]` prints as `a + b`, `Times[a, b]` as `a b` or `a*b`.
- **Precedence-aware parenthesization:** Inserts parentheses only when needed for correct precedence.
- **Special forms:** `Rational[1, 2]` prints as `1/2`, `Complex[a, b]` prints as `a + b I`, `List[a, b]` prints as `{a, b}`.
- **FullForm:** `FullForm[expr]` prints the raw tree structure for debugging.
- **InputForm:** `InputForm[expr]` prints in a form suitable for re-parsing.

---

## 13. Key Design Patterns

### Everything is an Expression
There is no separate type system for lists, matrices, rules, or patterns. A matrix is simply a `List` of `List`s. A rule `a -> b` is `Rule[a, b]`. A pattern `x_` is `Pattern[x, Blank[]]`. This uniformity means any expression can be manipulated with the same tools (`Part`, `Map`, `ReplaceAll`, etc.).

### Attributes Drive Evaluation
Rather than special-casing each function in the evaluator, behavior is controlled by attributes. Adding `ATTR_LISTABLE` to a function is all it takes to make it thread over lists. Adding `ATTR_FLAT | ATTR_ORDERLESS` makes it associative and commutative. This keeps the evaluator generic and the per-function code focused on mathematical logic.

### C for Performance, Rules for Mathematics
The system uses a two-tier approach: core operations (parsing, evaluation, pattern matching, arithmetic) are implemented in C for speed, while higher-level mathematical rules (derivatives, identities) are defined in the system's own Mathematica-like syntax. This makes the system self-hosting for rule-based mathematics.

### Fixed-Point Evaluation
The infinite evaluation loop means rules compose naturally. If `f[x_] := g[x]` and `g[x_] := x^2`, then `f[3]` automatically evaluates to `9` in two evaluation passes. This eliminates the need for explicit rule-chaining logic.

### Built-in Returns NULL for "I Can't Evaluate This"
The convention of returning `NULL` when a built-in can't handle its arguments is elegant: it lets symbolic arguments flow through naturally. `Sin[x]` remains as `Sin[x]` because `builtin_sin` returns `NULL` for a symbolic argument, and the evaluator simply keeps the expression as-is.
