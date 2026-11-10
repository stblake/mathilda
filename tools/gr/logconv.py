"""
logconv.py -- Mathilda's logarithm convention for G&R antiderivatives.

Mathilda's Integrate[1/x, x] is Log[x], not Log[Abs[x]]: an antiderivative is
analytic, and on every interval Log[u] and Log[Abs[u]] differ by a constant
(0 or I Pi).  G&R prints the real form ln|u|.  Every INDEFINITE record is
therefore rewritten Log[Abs[u]] -> Log[u], in its integrand and its value
alike, before it reaches the corpus (build_corpus.py) or the table
(build_table.py).  The rewritten formula is still exact: Log[u] has the same
derivative u'/u as Log[Abs[u]], and a Log[Abs[u]] in the integrand becomes
Log[u] on both sides consistently.  The table self-checks every answer by
differentiation at run time, so a record where the rewrite were not an
identity could not be served.

Definite values keep Abs: there Log[Abs[a/b]] is part of a real number, not a
constant of integration.
"""


def _close(s, i):
    """Index of the bracket closing the one opened at s[i] == '['."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == '[':
            depth += 1
        elif s[j] == ']':
            depth -= 1
            if depth == 0:
                return j
    return -1


def log_abs_to_log(s):
    """Rewrite every Log[Abs[u]] (the Log's whole argument an Abs) as Log[u]."""
    if not s or 'Log[Abs[' not in s:
        return s
    out, i = [], 0
    while True:
        k = s.find('Log[Abs[', i)
        if k < 0:
            out.append(s[i:])
            return ''.join(out)
        a = k + len('Log[Abs')           # the '[' of Abs
        e = _close(s, a)
        if e < 0 or e + 1 >= len(s) or s[e + 1] != ']':
            # Log[Abs[u]^2], Log[Abs[u]/2]: not the plain form, left as printed
            out.append(s[i:k + 4])
            i = k + 4
            continue
        out.append(s[i:k] + 'Log[' + log_abs_to_log(s[a + 1:e]) + ']')
        i = e + 2


def convert_record(v):
    """A copy of an oracle record with the convention applied (indefinite only)."""
    if v.get('kind') != 'indefinite':
        return v
    v = dict(v)
    for k in ('integrand', 'value', 'inst_integrand'):
        if v.get(k):
            v[k] = log_abs_to_log(v[k])
    if v.get('probes'):
        v['probes'] = [dict(p, f=log_abs_to_log(p.get('f', ''))) for p in v['probes']]
    return v


if __name__ == '__main__':
    for t in ['Log[Abs[x]]', 'x*Log[Abs[x^2 - a^2]] - 2*x', 'Log[Abs[(a + x)/(Abs[b] - x)]]',
              'Log[Abs[x]^2]', 'ArcTan[Abs[Sinh[x]]]', 'Log[Abs[Log[Abs[x]]]]']:
        print(t, '->', log_abs_to_log(t))
