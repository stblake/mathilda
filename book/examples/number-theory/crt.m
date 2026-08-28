# Section 4.6.5 (worked example): the Chinese Remainder Theorem, composed by hand.
mods = {3, 5, 7}
res = {2, 3, 2}
n = Times @@ mods
cofactor = n / mods
inverse = MapThread[PowerMod[#1, -1, #2] &, {cofactor, mods}]
x = Mod[Total[res cofactor inverse], n]
Mod[x, mods]
