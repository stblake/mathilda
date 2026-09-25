# 9.3 Separators: reading delimited data such as CSV.
# WordSeparators controls what breaks one token from the next.
c = OpenWrite["/tmp/mathilda_book_data.csv"];
WriteString[c, "1,2,3\n4,5,6\n"];
Close[c];
ReadList["/tmp/mathilda_book_data.csv", Number, WordSeparators -> {","}]
ReadList["/tmp/mathilda_book_data.csv", {Number, Number, Number}, WordSeparators -> {","}]
# TokenWords force a string to be read as its own word even with no surrounding
# separators, so an operator can be split from its operands.
t = OpenWrite["/tmp/mathilda_book_tok.txt"];
WriteString[t, "a+b\n"];
Close[t];
ReadList["/tmp/mathilda_book_tok.txt", Word]
ReadList["/tmp/mathilda_book_tok.txt", Word, TokenWords -> {"+"}]
