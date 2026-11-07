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

/* Total physical memory of the machine, in bytes.  False when the platform
 * offers no way to ask. */
bool meminfo_physical(uint64_t* bytes);

/* Number of logical processors online (backs $ProcessorCount).  Always >= 1:
 * a platform that cannot answer reports 1, the safe serial value. */
int meminfo_processor_count(void);

#endif /* MATHILDA_MEMINFO_H */
