---
source: src/core.c
---
A Protected OwnValue registered in `system_constants_init` (`src/core.c`). Under `USE_MPFR` it is the smallest positive value at machine precision (`mpfr_set_zero` then `mpfr_nextabove`, stored via `expr_new_mpfr_move`); without MPFR it collapses to `expr_new_real(DBL_MIN)`.
