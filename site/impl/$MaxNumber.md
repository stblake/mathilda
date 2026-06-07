---
source: src/core.c
---
A Protected OwnValue registered in `system_constants_init` (`src/core.c`). In a `USE_MPFR` build it is the largest finite value at machine precision (`DBL_MANT_DIG` bits), computed by `mpfr_set_inf` then `mpfr_nextbelow` and stored via `expr_new_mpfr_move`; without MPFR it collapses to `expr_new_real(DBL_MAX)`.
