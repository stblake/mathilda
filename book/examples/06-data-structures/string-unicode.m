# 11.4 Strings are sequences of characters, not bytes (outputs kept ASCII for the listings)
g = "\[Alpha]\[Beta]\[Gamma]";
StringLength[g]
ToCharacterCode[g]
ToCharacterCode[StringTake[g, {2}]]
ToCharacterCode[ToUpperCase[g]]
"\:03b1" === "\[Alpha]" === FromCharacterCode[945]
StringLength["caf\[EAcute]"]
StringMatchQ[g, LetterCharacter ..]
ToCharacterCode[g, "UTF8"]
# Named characters are read in source code too, not only inside strings (v0.587).
\[Alpha] + 1 /. \[Alpha] -> 2
ToCharacterCode[SymbolName[\[Mu]]]
{\[Pi] === Pi, FullForm[a \[Rule] b], FullForm[x \[Element] Reals]}
# The linear box syntax \!\( ... \) reads boxes as the expressions they display (v0.6001).
\!\(\*FractionBox[\(1\), \(2\)]\) + \!\(\*SqrtBox[\(16\)]\)
FullForm[\!\(\*SubscriptBox[\(x\), \(1\)]\)]
