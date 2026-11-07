/* test_meminfo.c -- MemoryInUse, MaxMemoryUsed, $SystemMemory, and the session constants.
 *
 * Memory figures are the awkward opposite of the usual test: the VALUE is machine-dependent
 * and cannot be pinned at all, so every assertion here has to be a relation that must hold
 * whatever the number is. The four that do the work:
 *
 *   - it is a positive integer, which catches a failed platform query returning zero;
 *   - the peak is never below the current, which is the definition of a high-water mark;
 *   - it RESPONDS to allocation, which is the one that distinguishes a real measurement from
 *     a plausible constant -- a hard-coded 16 MB would satisfy every other row here;
 *   - the argument form declines, rather than ignoring what it was passed.
 */
#include <stdio.h>
#include <limits.h>
#include "test_utils.h"
#include "meminfo.h"

extern void symtab_init(void);
extern void core_init(void);

static void test_memory_in_use_is_a_positive_integer(void) {
    /* A failed platform query must return unevaluated, never zero -- a zero would read as
     * "no memory in use", which is both false and the kind of wrong that looks fine in a
     * status bar. So this row also pins the decline-rather-than-lie contract. */
    assert_eval_eq("IntegerQ[MemoryInUse[]] && MemoryInUse[] > 0", "True", 0);
    assert_eval_eq("IntegerQ[MaxMemoryUsed[]] && MaxMemoryUsed[] > 0", "True", 0);
    /* Sanity on the ORDER OF MAGNITUDE, which is what catches a unit error. ru_maxrss is
     * bytes on Darwin and kilobytes on Linux, so a missing conversion is wrong by 1024x --
     * and this kernel links GMP, MPFR, LAPACK and Readline, so it cannot plausibly be under
     * a megabyte or over sixteen gigabytes. Wide enough never to be flaky, narrow enough
     * that a factor of 1024 in either direction fails it. */
    assert_eval_eq("1024^2 < MemoryInUse[] < 16*1024^3", "True", 0);
    assert_eval_eq("1024^2 < MaxMemoryUsed[] < 16*1024^3", "True", 0);
}

static void test_peak_is_never_below_current(void) {
    /* The defining property of a high-water mark. It is also what fails if the two are ever
     * wired to different sources with different units. */
    assert_eval_eq("MaxMemoryUsed[] >= MemoryInUse[]", "True", 0);
    /* And the peak cannot go DOWN between two calls. */
    assert_eval_eq("Module[{a, b}, a = MaxMemoryUsed[]; b = MaxMemoryUsed[]; b >= a]",
                   "True", 0);
}

static void test_it_responds_to_a_real_allocation(void) {
    /* THE row that proves this is a measurement rather than a constant. A hard-coded value
     * passes every other assertion in this file; only this one can tell the difference.
     *
     * Three million machine integers is 24 MB at 8 bytes each, which is far too much for the
     * allocator to satisfy from free space it already held, so resident memory has to grow.
     * The bound is deliberately loose -- 4 MB against an expected 24 -- because the exact
     * growth depends on page granularity and on what the allocator had in hand, and a tight
     * bound here would be flaky for no gain. The claim being tested is "it moves, and in the
     * right direction, by an amount of the right order", not any particular figure. */
    assert_eval_eq("Module[{before, after, big},"
                   " before = MemoryInUse[];"
                   " big = Range[3000000];"
                   " after = MemoryInUse[];"
                   " after - before > 4*1024^2]", "True", 0);
    /* The peak must have absorbed that spike too. */
    assert_eval_eq("Module[{big, peak},"
                   " big = Range[3000000];"
                   " peak = MaxMemoryUsed[];"
                   " peak >= MemoryInUse[]]", "True", 0);
}

static void test_argument_forms_decline(void) {
    /* Mathematica's one-argument MemoryInUse reports a subkernel's usage. There are no
     * subkernels here, so an argument is refused rather than silently ignored -- accepting
     * and dropping it would hide a real mistake in a caller that thought it meant
     * something. */
    assert_eval_eq("Head[MemoryInUse[1]]", "MemoryInUse", 0);
    assert_eval_eq("Head[MemoryInUse[a, b]]", "MemoryInUse", 0);
    assert_eval_eq("Head[MaxMemoryUsed[1]]", "MaxMemoryUsed", 0);
}

static void test_both_are_protected(void) {
    /* Every builtin gets its attributes; these are read-only system quantities, so a user
     * redefining them would make a status bar quietly report fiction. */
    assert_eval_eq("MemberQ[Attributes[MemoryInUse], Protected]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[MaxMemoryUsed], Protected]", "True", 0);
}

static void test_system_memory(void) {
    /* Total physical memory: a positive integer at least as large as what this process
     * holds resident, and fixed for the session. */
    assert_eval_eq("IntegerQ[$SystemMemory] && $SystemMemory > 0", "True", 0);
    assert_eval_eq("$SystemMemory >= MaxMemoryUsed[]", "True", 0);
    assert_eval_eq("$SystemMemory === $SystemMemory", "True", 0);
    /* A package's memory budget derived from it is an integer. */
    assert_eval_eq("IntegerQ[Floor[$SystemMemory/10^6/4]]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[$SystemMemory], Protected]", "True", 0);
}

static void test_session_constants(void) {
    /* A standalone kernel: no notebook, no front end, no cloud, no subkernels. */
    assert_eval_eq("{$Notebooks, $FrontEnd, $CloudEvaluation, $KernelID, $KernelCount}",
                   "{False, Null, False, 0, 0}", 0);
    assert_eval_eq("$PathnameSeparator", "\"/\"", 0);
    assert_eval_eq("{Attributes[$Notebooks], Attributes[$FrontEnd], Attributes[$CloudEvaluation], "
                   "Attributes[$KernelID], Attributes[$KernelCount], Attributes[$PathnameSeparator]}",
                   "{{Protected}, {Protected}, {Protected}, {Protected}, {Protected}, {Protected}}", 0);
    /* Protected: an assignment is refused and the value survives. */
    assert_eval_eq("Quiet[$KernelCount = 4]; $KernelCount", "0", 0);
}

static void test_system_identity(void) {
    assert_eval_eq("{StringQ[$System], StringQ[$SystemID]}", "{True, True}", 0);
    /* $SystemID: alphanumerics and dashes only. */
    assert_eval_eq("StringMatchQ[$SystemID, RegularExpression[\"[A-Za-z0-9-]+\"]]", "True", 0);
    assert_eval_eq("{Attributes[$System], Attributes[$SystemID]}", "{{Protected}, {Protected}}", 0);
    assert_eval_eq("Quiet[$System = \"x\"]; $System === \"x\"", "False", 0);
    assert_eval_eq("Quiet[$SystemID = \"x\"]; $SystemID === \"x\"", "False", 0);
    /* $SystemWordLength: pointer width in bits. */
    {
        char want[8];
        snprintf(want, sizeof want, "%d", (int)(CHAR_BIT * sizeof(void*)));
        assert_eval_eq("$SystemWordLength", want, 0);
    }
    assert_eval_eq("MemberQ[{32, 64}, $SystemWordLength]", "True", 0);
    /* $ProcessorCount: a positive integer, the count meminfo reports. */
    assert_eval_eq("IntegerQ[$ProcessorCount] && $ProcessorCount >= 1", "True", 0);
    {
        char want[16];
        snprintf(want, sizeof want, "%d", meminfo_processor_count());
        assert_eval_eq("$ProcessorCount", want, 0);
    }
    assert_eval_eq("{Attributes[$SystemWordLength], Attributes[$ProcessorCount]}",
                   "{{Protected}, {Protected}}", 0);
    assert_eval_eq("Quiet[$ProcessorCount = 0]; $ProcessorCount >= 1", "True", 0);
#if defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
    assert_eval_eq("{$System, $SystemID}", "{\"Mac OS X ARM (64-bit)\", \"MacOSX-ARM64\"}", 0);
#elif defined(__APPLE__) && defined(__x86_64__)
    assert_eval_eq("{$System, $SystemID}", "{\"Mac OS X x86 (64-bit)\", \"MacOSX-x86-64\"}", 0);
#elif defined(__linux__) && defined(__x86_64__)
    assert_eval_eq("{$System, $SystemID}", "{\"Linux x86 (64-bit)\", \"Linux-x86-64\"}", 0);
#endif
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_memory_in_use_is_a_positive_integer);
    TEST(test_peak_is_never_below_current);
    TEST(test_it_responds_to_a_real_allocation);
    TEST(test_argument_forms_decline);
    TEST(test_both_are_protected);
    TEST(test_system_memory);
    TEST(test_session_constants);
    TEST(test_system_identity);

    printf("All meminfo tests passed.\n");
    return 0;
}
