# Module locals are reclaimed when nothing refers to them any more (v0.6056-58)
before = Names["Global`*"];
Do[Module[{a, b}, a = b = i; a + b], {10000}]
Complement[Names["Global`*"], before]
kept = Module[{c}, c = 3; Hold[c]]
Complement[Names["Global`*"], before]
Do[Module[{a}, a = Range[100]], {1000}]; m0 = MemoryInUse[];
Do[Module[{a}, a = Range[100]], {10000}]
MemoryInUse[] - m0
