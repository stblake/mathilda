#!/usr/bin/env python3
"""Experiment 30 -- String operations and regular expressions (Python column).

Same eight kernels as ``string_ops.m``, same order and sizes.

NOT COVERED BY THE EXISTING CORPUS AT ALL -- every one of the twenty existing
experiments is numeric, and Mathilda's string subsystem (28 files) has never been
timed against anything.

Python's ``re`` is a mature C engine and ``str`` methods are heavily optimised, so
unlike the sympy and networkx columns this IS a real compiled baseline.  Mathilda
routes through PCRE2 when USE_REGEX is compiled in, so the regex rows compare two
C engines.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import re

from harness import bench, check, require

require(["re:sub", "re:findall"])

unit = "the quick brown fox jumps over the lazy dog "
big = unit * 5000

VOWELS = re.compile("[aeiou]+")
VOWEL = re.compile("[aeiou]")

bench("StringLength of 200k chars", lambda: len(big))
check("StringLength of 200k chars", len(unit))

bench("StringReplace literal, 200k chars",
      lambda: big.replace("quick", "slow"))
check("StringReplace literal, 200k chars", len(unit.replace("quick", "slow")))

bench("StringSplit on space, 200k chars", lambda: big.split())
check("StringSplit on space, 200k chars", len(unit.split()))

bench("StringCount substring, 200k chars", lambda: big.count("the"))
check("StringCount substring, 200k chars", unit.count("the"))

bench("StringReverse 200k chars", lambda: big[::-1])
check("StringReverse 200k chars", len(unit[::-1]))

bench("StringCases regex, 200k chars", lambda: VOWELS.findall(big), reps=1)
check("StringCases regex, 200k chars", len(VOWELS.findall(unit)))

bench("StringReplace regex, 200k chars",
      lambda: VOWEL.sub("*", big), reps=1)
check("StringReplace regex, 200k chars", len(VOWEL.sub("*", unit)))

bench("Characters of 200k chars", lambda: list(big), reps=1)
check("Characters of 200k chars", len(list(unit)))
