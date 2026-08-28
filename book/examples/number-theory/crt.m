# Section 4.6.5 (worked example): the Chinese Remainder Theorem.
# ChineseRemainder solves the system of congruences directly; the by-hand
# construction that follows shows the classical recipe it implements.
mods = {3, 5, 7}
res = {2, 3, 2}
x = ChineseRemainder[res, mods]
Mod[x, mods]
n = Times @@ mods
cofactor = n / mods
inverse = MapThread[PowerMod[#1, -1, #2] &, {cofactor, mods}]
Mod[Total[res cofactor inverse], n]
ChineseRemainder[res, mods, 100]
ChineseRemainder[{1, 2}, {6, 10}]
