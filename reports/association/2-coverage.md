## 2. Feature coverage vs Mathematica 15

**Method.** I wrote 333 probe expressions into two scripts (`/tmp/assocreport/cases.m`, `/tmp/assocreport/cases2.m`). Each one runs through `ToString[Quiet[Check[ToExpression[s], $Failed]], InputForm]` in both `./mathilda` and `wolframscript` (Mathematica 15), and a script diffs the two outputs line by line. **219 of the 333 probes (66%) give identical output.** To check whether a head is defined, I used Mathilda's `Names[...]` and grepped `symtab_add_builtin` in `src/`. Association builtins are registered in `src/assoc.c` (25 heads). Heads are grouped into 81 inventory rows.

A row is **Full** when every probed form matches Mathematica. It is **Partial** when the head exists but at least one common form is missing or gives a different answer. It is **Missing** when the head is undefined or never evaluates on associations. It is **N/A** when Mathematica itself rejects the form, or the difference is not about associations.

**Summary: Full 31 · Partial 30 · Missing 18 · N/A 2 (81 rows).**

A pattern cuts across many rows: **no curried operator form evaluates.** None of `Lookup[k][a]`, `KeyTake[ks][a]`, `KeyDrop[k][a]`, `KeySelect[p][a]`, `KeyMap[f][a]`, `KeyExistsQ[k][a]`, `KeyValueMap[f][a]`, `Select[p][a]`, `GroupBy[f][l]`, `CountsBy[f][l]`, `AssociationMap[f][l]`, `Merge[f][l]` or `Apply[f][x]` evaluates. Each is left unevaluated. The rows below count this against each head as "op-form missing".

| Head | Status | Notes |
|---|---|---|
| Association / `<\|..\|>` | Partial | Construction, duplicate keys (last wins), flattening of nested `<\|a->1,<\|b->2\|>\|>`, `:>` values and `Association[{rules}]` all match. `<\|a->1, 5\|>` stays `Association[a->1, 5]`, and `AssociationQ` of it is `True` (Mathematica: `False`). Structural differences are in the "atomicity" row. |
| AssociationQ | Partial | `AssociationQ[<\|a->1, 5\|>]` gives MT `True`, WL `False`. The other probes match. |
| Keys | Partial | `Keys[{<\|a->1\|>,<\|b->2\|>}]` gives WL `{{a},{b}}`, MT unevaluated. `Keys[a, f]` gives WL `{f[a],f[b]}`, MT unevaluated. `Keys[{a->1,b->2}]` works. |
| Values | Partial | `Values[a, f]` is unevaluated (WL `{f[1],f[2]}`). The rule-list form works. |
| KeyValueMap | Partial | The 2-argument form matches. Op-form missing: `KeyValueMap[List][<\|a->1\|>]` gives WL `{{a,1}}`. |
| Lookup | Partial | Single key, key list, default, list of associations and rule lists all match. Op-form missing: `Lookup[a][<\|a->5\|>]` gives WL `5`. **The default is evaluated eagerly:** `Lookup[<\|a->1\|>, a, Print["side"]; 0]` prints `side` in MT but not in WL. |
| KeyExistsQ | Partial | The 2-argument form matches. Op-form `KeyExistsQ[a][assoc]` is unevaluated (WL `True`). |
| MemberQ | Full | Tests values, not keys, in both systems. |
| KeyMemberQ | Partial | `KeyMemberQ[<\|a->1,b->2\|>, _]` gives MT `False`, WL `True`: the key is not treated as a pattern. |
| KeyFreeQ | Full | |
| KeyTake | Partial | Result order follows the association, not the key list: `KeyTake[<\|a->1,b->2,c->3\|>, {c,a}]` gives WL `<\|c->3, a->1\|>`, MT `<\|a->1, c->3\|>`. Op-form missing. |
| KeyDrop | Partial | Single key and key list match. Op-form `KeyDrop[b][assoc]` is unevaluated. |
| KeySelect | Partial | The 2-argument form matches. Op-form missing. |
| KeySort | Partial | The 1-argument form matches. `KeySort[a, Order]` and `KeySort[a, -Order[#1,#2]&]` are unevaluated (WL `<\|a->2,b->3,c->1\|>` and `<\|c->1,b->3,a->2\|>`). |
| KeySortBy | Full | Correct with `-#&`. `KeySortBy[a, Minus]` gives the wrong order only because `Minus[3]` does not evaluate in Mathilda (a general gap, not an association one). |
| KeyMap | Partial | `KeyMap[f, a]` and `KeyMap[<\|a->x\|>, a]` match. Op-form missing. |
| KeyUnion | Partial | The 1-argument form matches. `KeyUnion[{..}, 0&]` is unevaluated (WL fills in missing keys with `0`). |
| KeyIntersection | Missing | Undefined. WL: `KeyIntersection[{<\|a->1,b->2\|>,<\|b->3,c->4\|>}]` gives `{<\|b->2\|>,<\|b->3\|>}`. |
| KeyComplement | Missing | Undefined. WL: `KeyComplement[{<\|a->1,b->2\|>,<\|b->3\|>}]` gives `<\|a->1\|>`. |
| AssociateTo | Full | Single rule and rule list, including overwriting a key. |
| KeyDropFrom | Missing | Only interned in `sym_names.c`, never registered. `KeyDropFrom[x, a]` returns itself unevaluated and leaves `x` unchanged. |
| AssociationMap | Partial | The list input matches. `AssociationMap[Reverse, <\|a->1,b->2\|>]` is unevaluated (WL `<\|1->a,2->b\|>`). Op-form missing. |
| AssociationThread | Full | `{k..},{v..}`, `{k..}->{v..}` and duplicate keys match. With a length mismatch, MT leaves it unevaluated with no message, where WL emits a message. |
| Merge | Partial | `Merge[{assocs}, f]` matches for `Total`, `Identity`, `Join`, `First`, `Last` and `Max`. `Merge[{a->1, a->2}, Total]` (rule-list input) is unevaluated (WL `<\|a->3\|>`). Op-form missing. `Merge[.., Apply[Join]]` gives `<\|a -> Apply[Join][{{1},{2}}]\|>` because `Apply`'s op-form is missing. |
| GroupBy | Partial | `f`, `f->g`, `f->g, reducer`, association input and `Key[..]->Key[..]` all match. **The multi-level form is wrong, not unevaluated:** `GroupBy[{1,..,6}, {OddQ, #>3&}]` gives MT `<\|{OddQ,#1>3&}[1]->{1}, ...\|>`, WL `<\|True-><\|False->{1,3},True->{5}\|>, ...\|>`. Op-form missing. |
| Counts | Full | |
| CountsBy | Partial | The 2-argument form matches. Op-form missing. |
| PositionIndex | Partial | List input matches. Association input is unevaluated (WL `<\|a->{x,z}, b->{y}\|>`). |
| Normal | Partial | Top-level and nested `<\|a-><\|b->1\|>\|>` match. It does not reach associations inside other expressions: `Normal[{<\|a->1\|>, f[<\|b->2\|>]}]` gives WL `{{a->1}, f[{b->2}]}`, MT the input unchanged. |
| Query | Missing | Undefined. WL: `Query["a"][<\|"a"->1\|>]` gives `1`, `Query[All, Key[a]][{..}]` gives `{1,2}`. |
| Dataset | Missing | Undefined. `Dataset[..][All, a]` stays inert. |
| Part (read) | Partial | `[[1]]`, `[[-1]]`, `[["a"]]`, `[[Key[b]]]`, nested `[[Key[a], Key[x]]]`, `a[k1, k2]` and `list[[All, Key[a]]]` match. **Several forms give wrong answers.** `[[{1,2}]]` gives MT `{1,2}`, WL `<\|a->1,b->2\|>`. `[[{Key[a]}]]` gives MT `{1}`, WL `<\|a->1\|>`. `[[{-1,1}]]` gives MT `{3,1}`, WL `<\|c->3,a->1\|>`. `[[All]]` gives MT `Missing["KeyAbsent", All]`. `[[2;;]]` gives MT `Missing["KeyAbsent", Span[2, All]]` (WL `<\|b->2\|>`). `[[Key[b]]]` on an absent key gives MT `Missing["KeyAbsent", b]`, WL `...Key[b]]`. |
| Assignment `x[k]=v`, `x[[Key[k]]]=v`, `AppendTo`/`PrependTo`, `AppendTo[x[k], v]` | Full | Includes `Do[x[i]=i^2, ..]` and nested `x[a, b] = 5`. |
| Compound mutation (`x[k]+=`, `-=`, `*=`, `++`, `--`, `++x[k]`, `x[k]=.`, `x[a][b]=v`, `x[[k, i]]=v`, `x[k][[i]]=v`) | Missing | Every form fails. The operators print `AddTo::rvalue` / `Increment::rvalue` / ... and the association is left unchanged. The counter idiom `If[KeyExistsQ[c,#], c[#]+=1, c[#]=1]` over `{x,y,x}` gives MT `<\|x->1,y->1\|>`, WL `<\|x->2,y->1\|>`, silently. |
| Extract | Full | `Key[b]`, `{Key[b]}` and a positional index. |
| ReplacePart | Full | `Key[a]->v` and `2->v`. MT also accepts a bare `a->v`, which WL rejects. |
| Map / MapIndexed / MapAt | Full | `MapAt` accepts both `Key[a]` and positional indices. |
| Apply | Partial | `f @@ a` matches. `Apply[f, <\|a->{1,2}\|>, {1}]` gives MT `Association[f[a, {1,2}]]` (malformed), WL `<\|a->f[1,2]\|>`. |
| Select | Partial | The 2- and 3-argument forms match. Op-form missing. |
| Take / Drop / Most / Rest | Full | |
| First / Last | Full | `First[<\|\|>, d]` matches. |
| Sort | Partial | `Sort[a]` matches. **The ordering function is silently ignored:** `Sort[<\|a->3,b->1,c->2\|>, Greater]` and `Sort[a, #1>#2&]` give MT `<\|a->3,b->1,c->2\|>`, WL `<\|a->3,c->2,b->1\|>`. `Sort[{3,1,2}, Greater]` is correct. |
| SortBy | Full | Works with `-#&`, `Key["a"]` and `#["a"]&`. `SortBy[a, Minus]` differs only because of the general `Minus[3]` gap. |
| ReverseSort / ReverseSortBy | Full | |
| Reverse | Full | |
| Length / Total / Max / Mean (one association) | Full | |
| Join | Full | Later keys win and keep their first position. |
| Append / Prepend | Full | Rule, existing key, and association argument. |
| Insert | Partial | Positional form matches. `Insert[a, c->3, Key[b]]` is a no-op in MT (WL `<\|a->1,c->3,b->2\|>`). |
| Delete | Partial | Positional forms match. **`Delete[a, Key[a]]` is a silent no-op** in MT (WL `<\|b->2\|>`). |
| KeyValuePattern | Full | Single rule, rule list, `/;` condition, `Replace` with bindings, and matching rule lists. |
| Missing | Full | Not in `Names[]` but works (`<\|a->1\|>[b]` gives `Missing["KeyAbsent", b]`). |
| MissingQ | Missing | Undefined. `MissingQ[Missing[]]` is unevaluated. |
| DeleteMissing | Full | On associations and lists. |
| Key | Full | `Key[a][assoc]` and absent-key `Missing`. Not listed in `Names[]`. |
| TakeLargest / TakeSmallest / TakeLargestBy / TakeSmallestBy | Full | |
| MaximalBy / MinimalBy | Full | |
| Catenate | Partial | **Wrong answer:** `Catenate[{<\|a->1\|>,<\|b->2\|>}]` gives MT `<\|a->1,b->2\|>`, WL `{1,2}`. An association of associations is unevaluated. |
| JoinAcross | Missing | Undefined. |
| Transpose (associations) | Missing | `Transpose[<\|a-><\|x->1,..\|>, ..\|>]` is unevaluated (WL swaps the key levels to `<\|x-><\|a->1,b->3\|>, ..\|>`). A list of associations is also unevaluated. |
| DeleteDuplicates / DeleteDuplicatesBy | Full | |
| ApplyTo | Missing | Undefined. `ApplyTo[x, f]` leaves `x` unchanged. |
| Nothing | Full | Both as a literal value and returned from `Map`. |
| Named slots `#name` | Missing | **This is a parser gap:** `#name` is read as `#1*name`. `#name&[<\|"name"->1\|>]` gives MT `<\|"name"->1\|> name`, WL `1`. `Query[Select[#a>1&]]` is also mis-parsed. |
| `#["key"]&` and `Slot["x"]` calls | Full | `#["a"]&` inside `Map`, `Select` and `SortBy` matches. Only the printed form of an unapplied `Slot["x"]` differs. |
| Arithmetic threading over associations | Missing | `<\|a->1,b->2\|> + 1` gives MT `1 + <\|..\|>`, WL `<\|a->2,b->3\|>`. `a*2`, `a + a2`, `Sin[a]`, `Total[{a1, a2}]` and `Mean[{a1, a2}]` all stay symbolic. |
| Equal / Unequal | Partial | Identical associations give `True`. Different ones (`<\|a->1,b->2\|> == <\|b->2,a->1\|>`, `Unequal[<\|a->1\|>,<\|a->2\|>]`) stay unevaluated (WL `False` / `True`). `SameQ` matches. |
| SelectFirst / AllTrue / AnyTrue / NoneTrue | Full | |
| Discard | Missing | Undefined. |
| CountDistinct | Missing | Undefined. |
| DeleteCases / Cases / Position / Count / FreeQ | Full | |
| Pick | Partial | `Pick[a, {True, False}]` matches. `Pick[a, <\|a->True, b->False\|>]` builds a malformed `Association[Rule[1], Rule[]]`. |
| MapThread | Full | Over a list of associations. |
| Union / Intersection / Complement | Full | |
| SubsetQ | Missing | `SubsetQ[<\|a->1,b->2\|>, <\|a->1\|>]` is unevaluated (WL `True`). |
| Atomicity (AtomQ / Depth / Level / ReplaceAll into keys) | Partial | Mathilda treats an association as an ordinary compound expression. `AtomQ` gives MT `False` (WL `True`). `Depth` gives MT `3` (WL `2`). `Level[a,{1}]` gives MT `{a->1, b->{2}}` (WL `{1,{2}}`). **`<\|a->1\|> /. a->z` gives MT `<\|z->1\|>`**, while WL leaves keys alone. |
| ExportString / ImportString "RawJSON" | Missing | Undefined, so there is no JSON round trip for associations. |
| StringTemplate / TemplateApply | Missing | Undefined. These take associations as their slot source. |
| Splice | Missing | Undefined. |
| Differences / Accumulate / Partition / GatherBy on an association | N/A | Mathilda extensions (`Accumulate[<\|a->1,b->4\|>]` gives `<\|a->1,b->5\|>`). WL raises an error on all four. |
| Through | N/A | Both systems route `<\|..\|>[x]` through key lookup, giving `Missing`. |

### Top 10 missing or partial items, by real-world usage

1. **Compound mutation of association entries** (`c[k] += 1`, `c[k]++`, `c[k] =.`). This is the core counter/accumulator idiom. It fails silently: the association is left unchanged and the code keeps going with wrong data.
2. **Named slots `#name`.** This is the most common way to read record fields in `Select`, `SortBy`, `Map` and `Query`. It is mis-parsed as multiplication, so it gives wrong answers rather than failing loudly. `#["name"]` works.
3. **Operator (curried) forms** of `Lookup`, `Select`, `KeyTake`, `KeyDrop`, `KeySelect`, `KeyMap`, `GroupBy`, `CountsBy`, `KeyValueMap`, `AssociationMap`, `Merge` and `Apply`. These are used constantly in `/@`, `//` and `Query` pipelines.
4. **Part with lists and spans** (`a[[{1,2}]]`, `a[[2;;]]`, `a[[All]]`). These return a List or `Missing` where Mathematica returns an Association.
5. **Arithmetic threading over associations** (`assoc + 1`, `2 assoc`, `Total[{assoc1, assoc2}]`, `Mean` of a list of associations). This is common in statistics and aggregation code.
6. **`Query` / `Dataset`.** These are the main tabular-data interface. They are entirely absent.
7. **`Delete[a, Key[k]]` / `KeyDropFrom`.** The first is a silent no-op and the second is undefined. Both are the standard ways to remove a key in place.
8. **JSON (`ImportString`/`ExportString` "RawJSON").** This is the usual way associations enter and leave a program.
9. **`Sort[a, p]` ignoring `p`, `KeyTake` ordering, and multi-level `GroupBy`.** All three return plausible but wrong results.
10. **`Keys`/`Values` over lists of associations and with a function argument, plus `KeySort[a, p]`.** Simple variants that come up often in reporting code.

### Key findings

- [strength] The core association surface is solid, and the order and last-wins semantics match Mathematica. This covers construction, `Keys`/`Values`, `Lookup` (defaults, key lists, lists of associations), `KeyExistsQ`, `Counts`, the 2-argument `GroupBy` family, `Merge`, `AssociationThread`, `Join`/`Append`/`Prepend`, `Map`/`Select`/`Take`/`Drop`, `TakeLargest`/`MaximalBy`, `KeyValuePattern` and simple `x[k] = v` assignment. 219 of the 333 probes are byte-identical.
- [risk] Silent wrong answers matter more than missing heads. Five cases return plausible data with no error: failed `+=`/`++` on entries, `#name` parsed as `#1*name`, `Delete[a, Key[k]]` as a no-op, `Sort[a, p]` ignoring `p`, and `Catenate` of associations returning an association. Multi-level `GroupBy` also returns garbage keys, and `Pick`/`Apply` at level 1 build malformed `Association[...]` expressions.
- [gap] No curried operator forms evaluate. About 12 heads are affected, and so is `Apply[f]` inside combiners such as `Merge[.., Apply[Join]]`. This single cross-cutting gap breaks most pipeline-style code.
- [gap] Missing heads: `Query`, `Dataset`, `KeyDropFrom`, `KeyIntersection`, `KeyComplement`, `JoinAcross`, `ApplyTo`, `MissingQ`, `Discard`, `CountDistinct`, `Splice`, `Transpose` of associations, RawJSON import/export, and `StringTemplate`.
- [risk] Associations are not atomic. `AtomQ` is `False`, `Depth`/`Level` see the internal rules, and `ReplaceAll` rewrites keys (`<\|a->1\|> /. a->z` gives `<\|z->1\|>`). Code that runs rules over data containing associations will corrupt keys where Mathematica would not.
- [gap] `Part` diverges from Mathematica's association-preserving semantics for lists of positions and keys, `All`, and spans. `AssociationQ` accepts invalid `<\|a->1, 5\|>`. `Lookup` evaluates its default even when the key is present.
- [gap] `Listable`/`NumericFunction` heads (`Plus`, `Times`, `Sin`, `Total` over lists) do not thread through association values, and `==` on unequal associations stays unevaluated.
- [strength] Mathilda goes beyond Mathematica in a few places: `Accumulate`, `Differences`, `Partition` and `GatherBy` on associations, and a bare-key `ReplacePart`. These are harmless, but they are non-standard and code written against them will not port back to Mathematica.
