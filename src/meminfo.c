/* meminfo.c -- MemoryInUse, MaxMemoryUsed, $SystemMemory and $ProcessorCount.
 *
 * WHAT THESE MEASURE, AND HOW IT DIFFERS FROM MATHEMATICA. Wolfram's MemoryInUse[] reports
 * the bytes used to store the data of the current session -- expressions, definitions,
 * caches -- and nothing else. Mathilda reports the process's RESIDENT SET SIZE, which also
 * includes the binary, the shared libraries (GMP, MPFR, LAPACK, Readline, and on a graphics
 * build Raylib), the stacks, and whatever the allocator is holding but has not returned to
 * the OS.
 *
 * That is a real difference and it is documented rather than papered over, because the
 * numbers are not interchangeable: on a freshly started kernel Wolfram's figure is small
 * while this one is tens of megabytes of mapped libraries. Reporting session-data bytes
 * exactly would mean routing every allocation in the tree through a counting wrapper --
 * about 500 modules, all of them calling malloc directly today -- and the resulting number
 * would still miss what the allocator retains.
 *
 * RSS is also the more USEFUL quantity for the thing that prompted this: a notebook status
 * bar. It is the number Activity Monitor and top show, so a user comparing the two sees them
 * agree, where a smaller session-only figure would look like an under-report.
 *
 * PORTABILITY, which has two traps in it.
 *
 * `getrusage` is POSIX, not C99, so the feature-test macro has to come BEFORE any include --
 * below the first one the header has already been parsed with the wrong namespace and the
 * macro does nothing. glibc hides the symbol under -std=c99 while Darwin exposes it anyway,
 * so getting this wrong compiles clean here and breaks only on Linux. Same pattern as
 * src/ndkernels.c, src/core.c, src/repl.c.
 *
 * `ru_maxrss` IS NOT IN THE SAME UNIT ON BOTH PLATFORMS. Darwin reports BYTES; Linux
 * reports KILOBYTES. A single unconverted use is wrong by a factor of 1024 on one of the two
 * and silently plausible on both -- exactly the class of bug the int64_t/long long note in
 * CLAUDE.md is about. Hence the explicit per-platform scaling below.
 *
 * AND THE FEATURE-TEST MACRO RUNS THE OPPOSITE WAY FROM EVERY OTHER FILE IN THE TREE, which
 * is why the guard below is conditional rather than the usual unconditional _XOPEN_SOURCE.
 * `ru_maxrss` is a BSD extension, not POSIX -- POSIX only requires ru_utime and ru_stime. So
 * on Darwin, defining _XOPEN_SOURCE selects the strict POSIX subset and HIDES the field:
 * `no member named 'ru_maxrss' in 'struct rusage'`, from the very macro that every other
 * module needs in order to see POSIX at all. Meanwhile glibc needs a feature macro to expose
 * `getrusage` itself. Each platform therefore gets the macro that widens ITS namespace --
 * _DARWIN_C_SOURCE for the BSD fields, _DEFAULT_SOURCE plus _XOPEN_SOURCE for glibc -- and
 * both still come before any include, where a feature-test macro has to be to do anything.
 */
#if defined(__APPLE__)
#  ifndef _DARWIN_C_SOURCE
#  define _DARWIN_C_SOURCE
#  endif
#elif !defined(_WIN32)
#  ifndef _DEFAULT_SOURCE
#  define _DEFAULT_SOURCE
#  endif
#  ifndef _XOPEN_SOURCE
#  define _XOPEN_SOURCE 600
#  endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined(_WIN32)
#include <sys/resource.h>   /* getrusage -- POSIX, absent on native Windows */
#endif

#if defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>     /* sysctlbyname("hw.memsize") */
#endif
#if defined(__linux__)
#include <unistd.h>
#endif
#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>          /* GetProcessMemoryInfo; link -lpsapi on Windows */
#endif

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "meminfo.h"

bool meminfo_current(uint64_t* bytes) {
    if (!bytes) return false;

#if defined(__APPLE__)
    /* mach_task_basic_info gives resident_size directly, in bytes. */
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  (task_info_t)&info, &count) != KERN_SUCCESS)
        return false;
    *bytes = (uint64_t)info.resident_size;
    return true;

#elif defined(__linux__)
    /* /proc/self/statm: total pages, then RESIDENT pages. In pages, so scale by the page
     * size -- assuming 4096 would be wrong on a 16K-page arm64 kernel. */
    FILE* f = fopen("/proc/self/statm", "r");
    if (!f) return false;
    unsigned long total = 0, resident = 0;
    int got = fscanf(f, "%lu %lu", &total, &resident);
    fclose(f);
    if (got != 2) return false;
    long page = sysconf(_SC_PAGESIZE);
    if (page <= 0) return false;
    *bytes = (uint64_t)resident * (uint64_t)page;
    return true;

#elif defined(_WIN32)
    /* Current resident set == working set, in bytes. */
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return false;
    *bytes = (uint64_t)pmc.WorkingSetSize;
    return true;

#else
    /* No portable way to ask. Declining is the honest answer; a zero would read as "no
     * memory in use", which is worse than no answer at all. */
    (void)bytes;
    return false;
#endif
}

bool meminfo_peak(uint64_t* bytes) {
    if (!bytes) return false;

#if defined(_WIN32)
    /* No getrusage on native Windows; the peak working set is the direct
     * equivalent of ru_maxrss, already in bytes. */
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return false;
    *bytes = (uint64_t)pmc.PeakWorkingSetSize;
    return true;
#else
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return false;

#if defined(__APPLE__)
    *bytes = (uint64_t)ru.ru_maxrss;                /* Darwin: already bytes */
#else
    *bytes = (uint64_t)ru.ru_maxrss * (uint64_t)1024; /* Linux: kilobytes */
#endif
    /* Linux keeps the high-water mark lazily (hiwater_rss is folded in from the
     * split per-CPU RSS counters), so ru_maxrss can read a few pages BELOW the
     * live /proc/self/statm figure.  A peak is never below the current value. */
    uint64_t cur = 0;
    if (meminfo_current(&cur) && cur > *bytes) *bytes = cur;
    return true;
#endif
}

bool meminfo_physical(uint64_t* bytes) {
    if (!bytes) return false;

#if defined(__APPLE__)
    uint64_t mem = 0;
    size_t len = sizeof mem;
    if (sysctlbyname("hw.memsize", &mem, &len, NULL, 0) != 0 || mem == 0) return false;
    *bytes = mem;
    return true;

#elif defined(__linux__)
    /* Pages times page size; _SC_PHYS_PAGES is exposed by _DEFAULT_SOURCE above. */
    long pages = sysconf(_SC_PHYS_PAGES);
    long page = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || page <= 0) return false;
    *bytes = (uint64_t)pages * (uint64_t)page;
    return true;

#elif defined(_WIN32)
    MEMORYSTATUSEX ms;
    ms.dwLength = sizeof ms;
    if (!GlobalMemoryStatusEx(&ms)) return false;
    *bytes = (uint64_t)ms.ullTotalPhys;
    return true;

#else
    (void)bytes;
    return false;
#endif
}

/* meminfo_processor_count -- logical processors online.
 *
 * Logical (hardware threads), not physical cores: it is the count the
 * in-tree parallel kernels (src/linalg/mp_parallel.c, fi_par.c, gmet_core.c)
 * already size their thread pools by via sysconf(_SC_NPROCESSORS_ONLN), so
 * $ProcessorCount and the work those kernels actually schedule agree.  It is
 * read once at startup; the result is clamped to >= 1. */
int meminfo_processor_count(void) {
    long n = 0;
#if defined(__APPLE__)
    int v = 0;
    size_t len = sizeof v;
    if (sysctlbyname("hw.logicalcpu", &v, &len, NULL, 0) == 0) n = v;
#elif defined(__linux__)
    n = sysconf(_SC_NPROCESSORS_ONLN);   /* exposed by _DEFAULT_SOURCE above */
#elif defined(_WIN32)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    n = (long)si.dwNumberOfProcessors;
#endif
    if (n < 1) n = 1;
    if (n > 1048576) n = 1048576;
    return (int)n;
}

/* MemoryInUse[] -- current resident bytes.
 *
 * Takes no arguments. Mathematica's one-argument form reports a subkernel's usage, and there
 * are no subkernels here, so anything other than MemoryInUse[] returns unevaluated rather
 * than quietly ignoring what was passed. */
static Expr* builtin_memoryinuse(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    uint64_t b = 0;
    if (!meminfo_current(&b)) return NULL;
    return expr_new_integer((int64_t)b);
}

/* MaxMemoryUsed[] -- peak resident bytes over the life of the process.
 *
 * A genuine high-water mark from the OS, not the largest value some previous call to
 * MemoryInUse happened to observe. That distinction matters: a polled maximum would miss a
 * spike between two polls, and a status bar polling once a second would miss almost every
 * spike worth knowing about. */
static Expr* builtin_maxmemoryused(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    uint64_t b = 0;
    if (!meminfo_peak(&b)) return NULL;
    return expr_new_integer((int64_t)b);
}

void meminfo_init(void) {
    symtab_add_builtin("MemoryInUse", builtin_memoryinuse);
    symtab_get_def("MemoryInUse")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MemoryInUse",
        "MemoryInUse[] gives the number of bytes of memory currently resident for the "
        "Mathilda process. This is the process resident set size, so unlike "
        "Mathematica's MemoryInUse it also counts the binary, the shared libraries and "
        "whatever the allocator holds without returning it to the system -- it is the "
        "figure Activity Monitor and top report, not a count of session data alone. "
        "Returns unevaluated on a platform that offers no way to ask, rather than "
        "reporting zero.");

    /* $SystemMemory -- total physical memory, fixed for the session.  Left
     * undefined (rather than 0) on a platform that offers no way to ask. */
    uint64_t phys = 0;
    if (meminfo_physical(&phys)) {
        Expr* sym = expr_new_symbol("$SystemMemory");
        Expr* val = expr_new_integer((int64_t)phys);
        symtab_add_own_value("$SystemMemory", sym, val);
        expr_free(sym);
        expr_free(val);
    }
    symtab_get_def("$SystemMemory")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("$SystemMemory",
        "$SystemMemory\n\tgives the total number of bytes of physical memory on the computer.");

    /* $ProcessorCount -- logical processors online, fixed for the session. */
    {
        Expr* sym = expr_new_symbol("$ProcessorCount");
        Expr* val = expr_new_integer((int64_t)meminfo_processor_count());
        symtab_add_own_value("$ProcessorCount", sym, val);
        expr_free(sym);
        expr_free(val);
    }
    symtab_get_def("$ProcessorCount")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("$ProcessorCount",
        "$ProcessorCount\n\tgives the number of processor cores available on the computer system on "
        "which Mathilda is being run (logical processors online, as the parallel kernels count them). "
        "It is used by Parallelize and other parallel functions.");

    symtab_add_builtin("MaxMemoryUsed", builtin_maxmemoryused);
    symtab_get_def("MaxMemoryUsed")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("MaxMemoryUsed",
        "MaxMemoryUsed[] gives the peak number of bytes resident for the Mathilda process "
        "over its lifetime. The high-water mark comes from the operating system, so it "
        "catches spikes that occurred between two calls to MemoryInUse rather than only "
        "the largest value previously observed.");
}
