/*
 * checked_alloc.h -- allocations whose size a user chose.
 *
 * A builtin that sizes a buffer from its arguments (IdentityMatrix[n],
 * Array[f, n], PadLeft[l, n], CompleteGraph[n], RandomImage[.., {w, h}],
 * Series[f, {x, 0, n}], ...) must not multiply those sizes unchecked: the
 * product wraps (an undersized buffer, then a heap overrun), or it is honest
 * and enormous (calloc fails and the NULL is dereferenced, or the OS
 * overcommits and the kernel thrashes for minutes).  The 2026-11-14 audit
 * found ~25 crashes of that one shape.
 *
 * The contract here is the reference one: a request that cannot be met
 * reports General::nomem and aborts the computation (the builtin returns
 * mth_nomem_abort(), the Abort[] sentinel, which the top level prints as
 * $Aborted and CheckAbort can catch).  "Cannot be met" means any of:
 *   - a negative count, or count * size overflowing 64 bits;
 *   - more bytes than the machine's physical memory (meminfo_physical) --
 *     a buffer that large can only be paged, never computed with;
 *   - an active MemoryConstrained frame refusing it (mth_mem_request).
 *
 * Every check is a few integer operations; the physical size is read once.
 */
#ifndef MATHILDA_CHECKED_ALLOC_H
#define MATHILDA_CHECKED_ALLOC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct Expr;

/* a * b for non-negative counts.  False on a negative operand or overflow. */
bool mth_count_mul(int64_t a, int64_t b, int64_t* out);

/* Whether `count` elements of `elsize` bytes can be allocated: no overflow,
 * within physical memory, and allowed by MemoryConstrained.  On success the
 * byte count is stored in *bytes (when non-NULL).  No message. */
bool mth_alloc_ok(int64_t count, size_t elsize, size_t* bytes);

/* malloc / calloc of `count` elements of `elsize` bytes, checked as above.
 * NULL (after General::nomem) when the request cannot be met or the
 * allocator fails.  A zero count allocates one byte, so NULL always means
 * failure.  The memory is released with free(). */
void* mth_malloc_n(int64_t count, size_t elsize);
void* mth_calloc_n(int64_t count, size_t elsize);

/* Emit General::nomem (routed through mth_message, so Quiet and Check see
 * it) and return a fresh Abort[] sentinel for the builtin to return. */
struct Expr* mth_nomem_abort(void);

/* Emit General::nomem only. */
void mth_nomem_message(void);

#endif /* MATHILDA_CHECKED_ALLOC_H */
