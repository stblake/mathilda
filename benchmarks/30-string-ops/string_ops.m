(* Experiment 30 -- String operations and regular expressions.
   MEASURES src/strings/ (17 files) and src/strings/regex/ (11 files, PCRE2-backed
   when USE_REGEX is compiled in).

   NOT COVERED BY THE EXISTING CORPUS AT ALL.  Every one of the twenty existing
   experiments is numeric; the string subsystem has never been timed against
   anything.  Python's `re` is a mature C engine and its str methods are heavily
   optimised, so this is a real execution baseline rather than a pure-Python one.

   Checks are LENGTHS and COUNTS on small deterministic inputs. *)

Get["../harness.m"];
Get["../data.m"];

require[{"StringJoin", "StringSplit", "StringReplace", "StringLength",
         "StringTake", "StringPosition", "StringCases", "StringCount",
         "StringReverse", "StringRiffle", "RegularExpression", "StringMatchQ",
         "ToUpperCase", "Characters"}];

(* A deterministic ~200000-character subject built by repetition. *)
unit = "the quick brown fox jumps over the lazy dog ";
big = StringJoin[Table[unit, {5000}]];

bench["StringLength of 200k chars", StringLength[big];];
check["StringLength of 200k chars", StringLength[unit]];

bench["StringReplace literal, 200k chars",
  StringReplace[big, "quick" -> "slow"];];
check["StringReplace literal, 200k chars",
  StringLength[StringReplace[unit, "quick" -> "slow"]]];

bench["StringSplit on space, 200k chars", StringSplit[big];];
check["StringSplit on space, 200k chars", Length[StringSplit[unit]]];

bench["StringCount substring, 200k chars", StringCount[big, "the"];];
check["StringCount substring, 200k chars", StringCount[unit, "the"]];

bench["StringReverse 200k chars", StringReverse[big];];
check["StringReverse 200k chars", StringLength[StringReverse[unit]]];

benchIf["StringCases regex, 200k chars", "RegularExpression",
  StringCases[big, RegularExpression["[aeiou]+"]], 1];
checkIf["StringCases regex, 200k chars", "RegularExpression",
  Length[StringCases[unit, RegularExpression["[aeiou]+"]]]];

benchIf["StringReplace regex, 200k chars", "RegularExpression",
  StringReplace[big, RegularExpression["[aeiou]"] -> "*"], 1];
checkIf["StringReplace regex, 200k chars", "RegularExpression",
  StringLength[StringReplace[unit, RegularExpression["[aeiou]"] -> "*"]]];

bench["Characters of 200k chars", Characters[big], 1];
check["Characters of 200k chars", Length[Characters[unit]]];
