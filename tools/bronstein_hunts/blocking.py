import sys,collections
cnt=collections.Counter(); und=0; tot=0; cur=None
for l in sys.stdin:
    l=l.rstrip()
    if l=='CASE': cur=set(); continue
    if l.startswith('VERDICT '):
        tot+=1
        if l[8:]=='undec' and cur is not None:
            und+=1
            for m in cur: cnt[m]+=1
        cur=None; continue
    if cur is not None and l.startswith('[rde-cert] '):
        cur.add(l[11:].replace(' (no scope)',''))
print('cases',tot,'undecided',und)
for m,c in cnt.most_common(): print(c, m)
