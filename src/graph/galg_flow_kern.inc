/* galg_flow_kern.inc -- the capacity-typed kernels of galg_flow.c: Dinic's
 * blocking flow on a residual network, the residual reachability walks, and
 * the Nagamochi-Ibaraki global minimum cut.
 *
 * Included twice by galg_flow.c: once with GFK_T = int64_t (the common case,
 * every capacity sum below 2^61) and, where the compiler has a 128-bit
 * integer, once with GFK_T = gf_int for exact Integer capacities beyond that
 * (audit FindMinimumCut#3-4, FindMaximumFlow#3, EdgeConnectivity#2).  The
 * wide instance costs ~25-40% on the same graph, so the narrow one is used
 * whenever the network fits (GfNet.wide; gf_ni_mincut's dispatch).
 *
 *   GFK_T       capacity / flow type
 *   GFK_MAX     the "no limit" sentinel, above any flow
 *   GFK(name)   the instance's name for kernel `name`
 */

#define GFK_RES(N) ((GFK_T*)(N)->res)   /* the network's residuals, typed */

/* BFS levels from s over positive residual arcs; stops once t's level is
 * settled. Returns 1 iff t is reachable. */
static int GFK(gf_bfs)(GfNet* N, int s, int t) {
    int n = N->n;
    for (int v = 0; v < n; v++) N->level[v] = -1;
    int head = 0, tail = 0;
    N->level[s] = 0; N->queue[tail++] = s;
    while (head < tail) {
        int v = N->queue[head++];
        if (N->level[t] >= 0 && N->level[v] >= N->level[t]) break;
        int lv = N->level[v] + 1;
        for (long e = N->first[v]; e < N->first[v + 1]; e++) {
            int w = N->to[e];
            if (GFK_RES(N)[e] > 0 && N->level[w] < 0) {
                N->level[w] = lv;
                N->queue[tail++] = w;
            }
        }
    }
    return N->level[t] >= 0;
}

/* Dinic from s to t, stopping once the flow reaches `limit`. Returns the flow
 * pushed by this call (the residuals keep it). */
static GFK_T GFK(gf_dinic)(GfNet* N, int s, int t, GFK_T limit) {
    GFK_T flow = 0;
    if (s == t) return 0;
    long polls = 0;
    while (flow < limit && GFK(gf_bfs)(N, s, t)) {
        if ((++polls & 63) == 0) galg_poll();
        for (int v = 0; v < N->n; v++) N->it[v] = N->first[v];
        int top = 0, v = s;
        for (;;) {
            if (v == t) {
                GFK_T f = limit - flow;
                for (int i = 0; i < top; i++) if (GFK_RES(N)[N->pstack[i]] < f) f = GFK_RES(N)[N->pstack[i]];
                int cut = -1;
                for (int i = 0; i < top; i++) {
                    long e = N->pstack[i];
                    GFK_RES(N)[e] -= f; GFK_RES(N)[N->rev[e]] += f;
                    if (cut < 0 && GFK_RES(N)[e] == 0) cut = i;
                }
                flow += f;
                if (flow >= limit) return flow;
                top = cut;                           /* retreat to the bottleneck */
                v = top == 0 ? s : N->to[N->pstack[top - 1]];
                continue;
            }
            long e = N->it[v], end = N->first[v + 1];
            int lv = N->level[v] + 1;
            while (e < end && !(GFK_RES(N)[e] > 0 && N->level[N->to[e]] == lv)) e++;
            N->it[v] = e;
            if (e == end) {                          /* dead end: prune v */
                N->level[v] = -1;
                if (top == 0) break;
                top--;
                v = top == 0 ? s : N->to[N->pstack[top - 1]];
                N->it[v]++;
                continue;
            }
            N->pstack[top++] = e;
            v = N->to[e];
        }
    }
    return flow;
}

/* mark[v] = 1 for v reachable from s in the residual network. */
static void GFK(gf_reach_from)(const GfNet* N, int s, char* mark) {
    memset(mark, 0, (size_t)N->n);
    int head = 0, tail = 0;
    mark[s] = 1; N->queue[tail++] = s;
    while (head < tail) {
        int v = N->queue[head++];
        for (long e = N->first[v]; e < N->first[v + 1]; e++) {
            int w = N->to[e];
            if (GFK_RES(N)[e] > 0 && !mark[w]) { mark[w] = 1; N->queue[tail++] = w; }
        }
    }
}

/* mark[v] = 1 for v that can reach t in the residual network. */
static void GFK(gf_reach_to)(const GfNet* N, int t, char* mark) {
    memset(mark, 0, (size_t)N->n);
    int head = 0, tail = 0;
    mark[t] = 1; N->queue[tail++] = t;
    while (head < tail) {
        int w = N->queue[head++];
        for (long e = N->first[w]; e < N->first[w + 1]; e++) {
            int v = N->to[e];                       /* arc w->v; v->w is rev */
            if (GFK_RES(N)[N->rev[e]] > 0 && !mark[v]) { mark[v] = 1; N->queue[tail++] = v; }
        }
    }
}


/* Undirected global min cut by Nagamochi-Ibaraki. Weighted edges
 * (eu[k], ev[k], w[k]) over n >= 2 vertices; parallel edges are summed. Writes
 * side[v] = 1 for one shore of a minimum cut and returns its value; -1 on
 * allocation failure. */
static GFK_T GFK(gf_ni_mincut)(int n, long m, const int* eu, const int* ev, const GFK_T* w,
                            char* side) {
    GFK_T result = -1;
    size_t nn = (size_t)n + 1, mm = (size_t)(2 * m) + 1;
    /* current contracted graph (CSR) */
    int* off = malloc(nn * sizeof(int));
    int* adj = malloc(mm * sizeof(int));
    GFK_T* wt = malloc(mm * sizeof(GFK_T));
    int* noff = malloc(nn * sizeof(int));
    int* nadj = malloc(mm * sizeof(int));
    GFK_T* nwt = malloc(mm * sizeof(GFK_T));
    int* owner = malloc(nn * sizeof(int));     /* original v -> current node */
    int* uf = malloc(nn * sizeof(int));
    int* newid = malloc(nn * sizeof(int));
    GFK_T* key = malloc(nn * sizeof(GFK_T));
    char* done = malloc(nn);
    int* order = malloc(nn * sizeof(int));
    GFK_T* acc = malloc(nn * sizeof(GFK_T));
    int* stamp = malloc(nn * sizeof(int));
    int* heapv = malloc((mm + nn) * sizeof(int));
    GFK_T* heapk = malloc((mm + nn) * sizeof(GFK_T));
    if (!off || !adj || !wt || !noff || !nadj || !nwt || !owner || !uf || !newid || !key
        || !done || !order || !acc || !stamp || !heapv || !heapk) goto out;

    /* initial graph: merge parallel edges by stamping */
    {
        int* cnt = calloc(nn, sizeof(int));
        int* tu = malloc(mm * sizeof(int));
        GFK_T* tw = malloc(mm * sizeof(GFK_T));
        if (!cnt || !tu || !tw) { free(cnt); free(tu); free(tw); goto out; }
        for (long k = 0; k < m; k++) if (eu[k] != ev[k]) { cnt[eu[k]]++; cnt[ev[k]]++; }
        int s = 0;
        for (int v = 0; v < n; v++) { off[v] = s; s += cnt[v]; }
        off[n] = s;
        for (int v = 0; v < n; v++) cnt[v] = off[v];
        for (long k = 0; k < m; k++) {
            int a = eu[k], b = ev[k];
            if (a == b) continue;
            tu[cnt[a]] = b; tw[cnt[a]++] = w[k];
            tu[cnt[b]] = a; tw[cnt[b]++] = w[k];
        }
        for (int v = 0; v < n; v++) stamp[v] = -1;
        int wpos = 0;
        for (int v = 0; v < n; v++) {
            int lo = off[v], hi = off[v + 1], start = wpos;
            off[v] = start;
            for (int j = lo; j < hi; j++) {
                int x = tu[j];
                if (stamp[x] >= start) wt[stamp[x]] += tw[j];
                else { stamp[x] = wpos; adj[wpos] = x; wt[wpos] = tw[j]; wpos++; }
            }
        }
        off[n] = wpos;
        free(cnt); free(tu); free(tw);
    }
    for (int v = 0; v < n; v++) owner[v] = v;

    int cn = n;                 /* current node count */
    GFK_T lambda = GFK_MAX;
    long polls = 0;
    while (cn > 1) {
        if ((++polls & 7) == 0) galg_poll();
        /* bound by weighted degrees */
        int bestx = -1;
        for (int x = 0; x < cn; x++) {
            GFK_T d = 0;
            for (int j = off[x]; j < off[x + 1]; j++) d += wt[j];
            if (d < lambda) { lambda = d; bestx = x; }
        }
        if (bestx >= 0) for (int v = 0; v < n; v++) side[v] = owner[v] == bestx;
        if (lambda == 0) break;

        /* maximum-adjacency order with a lazy max-heap; contract edges whose
         * attachment reaches lambda. */
        for (int x = 0; x < cn; x++) { key[x] = 0; done[x] = 0; uf[x] = x; }
        int hn = 0, cnt_order = 0;
        for (int start = 0; start < cn; start++) {
            if (done[start]) continue;
            /* push start (a new component of the current graph) */
            heapv[hn] = start; heapk[hn] = 0; hn++;
            while (hn > 0) {
                /* pop max */
                int x = heapv[0]; GFK_T kx = heapk[0];
                hn--;
                if (hn > 0) {
                    int lv = heapv[hn]; GFK_T lk = heapk[hn]; int i = 0;
                    for (;;) {
                        int c = 2 * i + 1;
                        if (c >= hn) break;
                        if (c + 1 < hn && heapk[c + 1] > heapk[c]) c++;
                        if (heapk[c] <= lk) break;
                        heapv[i] = heapv[c]; heapk[i] = heapk[c]; i = c;
                    }
                    heapv[i] = lv; heapk[i] = lk;
                }
                if (done[x] || kx != key[x]) continue;   /* stale */
                done[x] = 1;
                order[cnt_order++] = x;
                for (int j = off[x]; j < off[x + 1]; j++) {
                    int y = adj[j];
                    if (done[y]) continue;
                    key[y] += wt[j];
                    if (key[y] >= lambda) {
                        /* q(x,y) >= lambda: contract */
                        int a = x, b = y;
                        while (uf[a] != a) a = uf[a] = uf[uf[a]];
                        while (uf[b] != b) b = uf[b] = uf[uf[b]];
                        if (a != b) uf[a] = b;
                    }
                    /* push (y, key[y]) */
                    int i = hn++;
                    while (i > 0) {
                        int p = (i - 1) / 2;
                        if (heapk[p] >= key[y]) break;
                        heapv[i] = heapv[p]; heapk[i] = heapk[p]; i = p;
                    }
                    heapv[i] = y; heapk[i] = key[y];
                }
            }
        }
        /* the last two of the order: lambda(s, t) = its attachment = deg(t)
         * >= lambda, so they may always be merged (guarantees progress) */
        {
            int a = order[cn - 2], b = order[cn - 1];
            while (uf[a] != a) a = uf[a] = uf[uf[a]];
            while (uf[b] != b) b = uf[b] = uf[uf[b]];
            if (a != b) uf[a] = b;
        }
        /* relabel and rebuild */
        int nc = 0;
        for (int x = 0; x < cn; x++) {
            int r = x;
            while (uf[r] != r) r = uf[r];
            uf[x] = r;
        }
        for (int x = 0; x < cn; x++) newid[x] = -1;
        for (int x = 0; x < cn; x++) if (uf[x] == x) newid[x] = nc++;
        for (int x = 0; x < cn; x++) newid[x] = newid[uf[x]];
        /* group old nodes by new id */
        int* grp = order;           /* reuse: old nodes sorted by new id */
        {
            int* c2 = stamp;        /* reuse as counts */
            for (int y = 0; y <= nc; y++) c2[y] = 0;
            for (int x = 0; x < cn; x++) c2[newid[x] + 1]++;
            for (int y = 0; y < nc; y++) c2[y + 1] += c2[y];
            for (int x = 0; x < cn; x++) grp[c2[newid[x]]++] = x;
            /* c2[y] is now the end of group y; starts are c2[y-1] */
        }
        for (int y = 0; y < nc; y++) stamp[y] = -1;
        int wpos = 0, gi = 0;
        for (int y = 0; y < nc; y++) {
            int start = wpos;
            noff[y] = start;
            while (gi < cn && newid[grp[gi]] == y) {
                int x = grp[gi++];
                for (int j = off[x]; j < off[x + 1]; j++) {
                    int z = newid[adj[j]];
                    if (z == y) continue;
                    if (stamp[z] >= start) nwt[stamp[z]] += wt[j];
                    else { stamp[z] = wpos; nadj[wpos] = z; nwt[wpos] = wt[j]; wpos++; }
                }
            }
        }
        noff[nc] = wpos;
        { int* t = off; off = noff; noff = t; }
        { int* t = adj; adj = nadj; nadj = t; }
        { GFK_T* t = wt; wt = nwt; nwt = t; }
        for (int v = 0; v < n; v++) owner[v] = newid[owner[v]];
        cn = nc;
    }
    result = lambda == GFK_MAX ? 0 : lambda;
out:
    free(off); free(adj); free(wt); free(noff); free(nadj); free(nwt); free(owner);
    free(uf); free(newid); free(key); free(done); free(order); free(acc); free(stamp);
    free(heapv); free(heapk);
    return result;
}

#undef GFK_RES
