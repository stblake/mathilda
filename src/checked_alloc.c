/*
 * checked_alloc.c -- allocations whose size a user chose.  See checked_alloc.h.
 */
#include "checked_alloc.h"

#include <stdlib.h>

#include "checked_int.h"
#include "expr.h"
#include "memacct.h"
#include "meminfo.h"
#include "message.h"
#include "sym_names.h"

/* Physical memory in bytes, read once; 0 until read, UINT64_MAX when the
 * platform cannot say (then only overflow and MemoryConstrained limit). */
static uint64_t g_phys_bytes = 0;

static uint64_t phys_bytes(void) {
    if (g_phys_bytes == 0) {
        uint64_t b = 0;
        g_phys_bytes = (meminfo_physical(&b) && b > 0) ? b : UINT64_MAX;
    }
    return g_phys_bytes;
}

bool mth_count_mul(int64_t a, int64_t b, int64_t* out) {
    int64_t r = 0;
    if (a < 0 || b < 0) return false;
    if (ci_mul_i64(a, b, &r)) return false;     /* true = overflow */
    if (out) *out = r;
    return true;
}

bool mth_alloc_ok(int64_t count, size_t elsize, size_t* bytes) {
    int64_t n = 0;
    if (count < 0 || elsize == 0) return false;
    if ((uint64_t)elsize > (uint64_t)INT64_MAX) return false;
    if (!mth_count_mul(count, (int64_t)elsize, &n)) return false;
    if ((uint64_t)n > phys_bytes()) return false;
    if ((uint64_t)n > (uint64_t)SIZE_MAX) return false;
    if (!mth_mem_request((size_t)n)) return false;
    if (bytes) *bytes = (size_t)n;
    return true;
}

void mth_nomem_message(void) {
    mth_message(SYM_General, "nomem",
                "The current computation was aborted because there was "
                "insufficient memory available to complete the computation.");
}

Expr* mth_nomem_abort(void) {
    mth_nomem_message();
    return expr_new_function(expr_new_symbol(SYM_Abort), NULL, 0);
}

static void* alloc_n(int64_t count, size_t elsize, bool zero) {
    size_t bytes = 0;
    void* p;
    if (!mth_alloc_ok(count, elsize, &bytes)) {
        mth_nomem_message();
        return NULL;
    }
    if (bytes == 0) bytes = 1;
    p = zero ? calloc(1, bytes) : malloc(bytes);
    if (!p) mth_nomem_message();
    return p;
}

void* mth_malloc_n(int64_t count, size_t elsize) { return alloc_n(count, elsize, false); }
void* mth_calloc_n(int64_t count, size_t elsize) { return alloc_n(count, elsize, true); }
