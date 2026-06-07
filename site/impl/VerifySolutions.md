---
source: src/solve.c
---
`VerifySolutions` is an option *symbol* for `Solve`, not a callable function. In `src/solve.c` it is listed in `is_known_option_name`, so `VerifySolutions -> _` is peeled off the argument list as a valid option (`is_option_arg`) rather than treated as a variable. However `apply_option` currently does nothing with it — the value is parsed and accepted but not yet wired into the polynomial solver (the docstring notes `Default: Automatic. Reserved.`). It exists so user code can pass the option without error.
