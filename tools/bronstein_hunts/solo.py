import sys,collections
cur=None; solo=collections.Counter(); combos=collections.Counter(); n=0
for l in sys.stdin:
    l=l.rstrip()
    if l=='CASE': cur=set(); continue
    if l.startswith('VERDICT '):
        if l[8:]=='undec' and cur is not None:
            n+=1; combos[tuple(sorted(cur))]+=1
        cur=None; continue
    if cur is not None and l.startswith('[rde-cert] '): cur.add(l[11:].replace(' (no scope)',''))
print('undecided',n)
for k,v in combos.most_common(12): print(v, ' + '.join(k) if k else '(NO MARK: silent decline)')
