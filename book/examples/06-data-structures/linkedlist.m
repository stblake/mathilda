# LinkedList: the object
ll = LinkedList[]
Do[AppendTo[ll, i^2], {i, 5}]; ll
{AtomQ[ll], ListQ[ll], Head[ll], Length[ll]}
Normal[ll]
{ll[[2]], ll[[-1]], ll[[2 ;; 4]]}
# Reference semantics: assignment shares the one object
ll2 = ll;
PrependTo[ll2, 0];
Normal[ll]
ll === ll2
copy = LinkedList[ll]; Push[copy, 99]; {Length[ll], Length[copy]}
LinkedList[{1, 2}] === LinkedList[{1, 2}]
Normal[LinkedList[{1, 2}]] === Normal[LinkedList[{1, 2}]]
