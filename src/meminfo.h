/* meminfo.h -- process memory introspection: MemoryInUse, MaxMemoryUsed. */
#ifndef MATHILDA_MEMINFO_H
#define MATHILDA_MEMINFO_H

#include <stdbool.h>
#include <stdint.h>

void meminfo_init(void);

/* Current resident set size in bytes. False means the platform gives no answer, which
 * callers must treat as "no answer" rather than as zero. */
bool meminfo_current(uint64_t* bytes);

/* Peak resident set size in bytes, over the life of the process. */
bool meminfo_peak(uint64_t* bytes);

#endif /* MATHILDA_MEMINFO_H */
