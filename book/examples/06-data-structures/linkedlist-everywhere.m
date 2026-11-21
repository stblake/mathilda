# Every head that works on a List works on a LinkedList
ll = LinkedList[{5, 3, 8, 1, 9, 2}];
Sort[ll]
Normal[Sort[ll]]
{Max[ll], Mean[ll], Total[ll]}
Position[ll, 8]
Normal[Select[ll, EvenQ]]
Normal /@ Normal[Partition[ll, 2]]
Normal[ll^2 + 1]
# Only Normal converts to a List; the counter proves it
c0 = Internal`LinkedListConversionCount[];
{Sort[ll], Total[ll], Sin[N[ll]], Partition[ll, 3]}; Internal`LinkedListConversionCount[] - c0
Normal[ll]; Internal`LinkedListConversionCount[] - c0
# Iterators see a snapshot taken when the loop begins
Do[AppendTo[ll, -i], {i, ll}]; Normal[ll]
