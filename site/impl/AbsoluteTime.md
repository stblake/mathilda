---
source: src/datetime.c
---
**Algorithm.** `builtin_absolute_time` returns seconds since the reference epoch 1900-01-01. With no argument it reads `time(NULL)` + `localtime`, splits the broken-down time into y/m/d/h/m/s, and combines `days_since_1900(...)*86400 + h*3600 + m*60 + s` (returned as `EXPR_REAL`). A numeric argument passes through as an existing absolute-time spec. A `{y,m,d,h,m,s}` `List` (1–6 elements, defaults `{y,1,1,0,0,0}`) is converted the same way; year and month must be integral, while day/hour/minute/second may be fractional (the day's fractional part is folded into seconds). When all inputs are integral and the total is representable, it returns an exact `EXPR_INTEGER`, otherwise an `EXPR_REAL`.

**Data structures.** Day numbers come from `days_since_1900`, which applies the Fliegel–Van Flandern Julian-Day-Number formula (offset by the JDN of 1900-01-01) and normalises out-of-range months explicitly. `ATTR_PROTECTED`.
