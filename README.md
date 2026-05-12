# PicoCAS

PicoCAS is a tiny, AI Agent-generated, symbolic computer algebra system (CAS) heavily inspired by the core architecture and evaluation semantics of Mathematica. Written entirely in C99, PicoCAS implements a recursive expression model, structural pattern matching, rewriting rules, and a surprisingly extensive library of built-in functions.

The name "PicoCAS" pays homage to David Stoutemyer's seminal PICOMATH-80 tiny computer algebra system. PicoCAS has the goal of implementing the core of the Mathematica programming language, and the algebraic and calculus capabilities of the Derive computer algebra system.

## 🌟 Key Features

* **Infinite Evaluation Semantics:** Expressions are repeatedly evaluated top-down until a fixed point is reached.
* **First-Class Pattern Matching:** Supports `Blank` (`_`), `BlankSequence` (`__`), `BlankNullSequence` (`___`), conditional matching (`/;`), and named bindings (`x_`).
* **Rule Engine:** Built-in support for transformation rules (`->`, `:>`) and repeated replacements (`/.`, `//.`).
* **Arbitrary Precision Arithmetic:** Powered natively by the GNU Multiple Precision Arithmetic Library (GMP) for handling operations on exceedingly large integers.
* **Advanced Integer Factorization:** Includes a unified, automatic factorization pipeline, alongside explicit execution of algorithms such as:
    * Pollard's Rho, Pollard's $P-1$, Williams' $P+1$, Fermat, CFRAC, Dixon's Method, and the Elliptic Curve Method (ECM).
* **Extensive Standard Library:** Built-in functions covering linear algebra, calculus, polynomial manipulation, combinatorics, statistics, list structures, and symbolic tensor operations.

---

## 🚀 Getting Started

### Prerequisites

To build and run PicoCAS, you need the following development libraries installed on your system:
* A C99-compliant compiler (`gcc` or `clang`)
* **GMP** (`libgmp` / `gmp-dev`)
* **Readline** (`libreadline` / `readline-dev`)
* **CMake** (Optional, only required if you intend to build the test suite)

*Note: The GMP-ECM library is utilized for factorization methods and is included as an external source dependency. It will be built automatically during the compilation process.*

### Building PicoCAS

The `makefile` automatically handles the configuration and compilation of internal dependencies before linking the main executable.

1. Clone the repository (including the bundled GMP-ECM submodule):
   ```bash
   git clone --recurse-submodules https://github.com/stblake/picocas.git
   cd picocas
   ```
   If you already cloned without `--recurse-submodules`, run
   `git submodule update --init --recursive` from inside the repo
   before building. (The build script will also attempt this
   automatically if it detects an uninitialized submodule.)
2. Build the project using `make`:
   ```bash
   make -j$(nproc)
   ```
3. Start the interactive REPL:
   ```bash
   ./picocas
   ```

### Running the Test Suite

PicoCAS boasts a comprehensive C-based unit test suite.

```bash
cd tests
mkdir -p build && cd build
cmake ..
make -j$(nproc)

# Run all test binaries
for t in *_tests; do ./$t; done
```

---

## 🛠️ Developer Guide & Architecture

Everything in PicoCAS is represented by an immutable-by-convention `Expr` AST node, implemented as a tagged union (`Integer`, `Real`, `BigInt`, `Rational`, `Complex`, `Symbol`, `String`, `Function`).

The system is modularized into several independent subsystems:
1. **Parser (`parse.c`)**: A robust Pratt parser mirroring Mathematica's exact operator precedences.
2. **Evaluator (`eval.c`)**: Manages the infinite evaluation loop. It processes attributes like `HoldAll`, `Flat`, `Orderless`, and `Listable` before recursively evaluating function arguments.
3. **Symbol Table (`symtab.c`)**: Stores global variable definitions, `DownValues` (user-defined rules), and native C built-in function pointers.
4. **Pattern Matcher & Rule Engine (`match.c`, `replace.c`)**: Implements structural tree unification and sequence segmenting through backtracking.

### Extending PicoCAS: Adding a New Built-in Function

Adding new functionality to PicoCAS is straightforward:

1. **Write the C Implementation:**
   Create your evaluation logic in the appropriate `.c` module (e.g., `core.c`). Your function signature must be `Expr* builtin_myfunc(Expr* res)`.
   * **Memory Rule:** Built-in functions take *ownership* of the expression passed to them. If you successfully evaluate the input and return a new `Expr*`, you MUST free the original expression (or structurally reuse its nodes). If the function cannot evaluate the arguments (e.g., symbolic inputs to a purely numeric function), return `NULL`, and the evaluator will automatically retain ownership.

   ```c
   Expr* builtin_myfunc(Expr* res) {
       if (res->data.function.arg_count != 1) return NULL;
       // ... mathematical logic ...
       expr_free(res); // Free input if returning a new expression
       return expr_new_integer(42);
   }
   ```

2. **Register the Function:**
   In the module's initialization routine (e.g., `core_init()`), register the function and assign its documentation string so it is available via `?MyFunc`:
   ```c
   symtab_add_builtin("MyFunc", builtin_myfunc);
   symtab_set_docstring("MyFunc", "MyFunc[x]\n\tComputes the ultimate answer.");
   ```

3. **Assign Attributes:**
   If your function automatically threads over lists, operates symmetrically, or holds its arguments unevaluated, assign the corresponding attributes during initialization:
   ```c
   symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
   ```

4. **Test and Document:**
   * Add test cases to the appropriate suite in the `tests/` directory using the `TEST(...)` macro.
   * Update `picocas_spec.md` with your new function, describing its behavior and providing sample inputs/outputs.

---

## 📜 Open Source & License

PicoCAS is open-source software licensed under the **GNU General Public License v3.0 (GPLv3)**. 

You are heavily encouraged to explore the codebase, submit pull requests, report issues, and expand the CAS with new mathematical algorithms! Please see the `LICENSE` file for more details.