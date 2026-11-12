/*
 * term_status.c -- see term_status.h.
 *
 * Threads and data races (C99 has no memory model; this is the contract):
 *   - The control state (state, the per-evaluation parameters) is guarded by
 *     one mutex + condition variable; begin/end hand over under it.
 *   - The status thread reads the progress stack, the dispatch head, the
 *     worker count and the interrupt counter WITHOUT locking. Each is a
 *     word-sized object written by one thread; a stale value only changes
 *     one frame of the drawing (term_progress.h has the details).
 *   - The status thread never touches an Expr, the symbol table or stdio;
 *     it writes with write(2) to its own descriptor.
 * Signals: the thread is created with every signal blocked, so SIGINT (the
 * REPL's handler), SIGPROF (TimeConstrained) and the fatal signals are always
 * delivered to the evaluating thread, exactly as before this module existed.
 */
#if defined(__APPLE__)
#  ifndef _DARWIN_C_SOURCE
#  define _DARWIN_C_SOURCE        /* TIOCGWINSZ, alongside POSIX */
#  endif
#elif !defined(_WIN32)
#  ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 200809L
#  endif
#  ifndef _DEFAULT_SOURCE
#  define _DEFAULT_SOURCE         /* TIOCGWINSZ on glibc */
#  endif
#endif

#include "term_status.h"
#include "term_spinner.h"
#include "term_progress.h"
#include "term_color.h"
#include "term_interrupt.h"

#include "expr.h"
#include "eval.h"
#include "meminfo.h"
#include "mtime.h"
#include "context.h"
#include "sym_names.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#if defined(MATHILDA_THREADS) && !defined(_WIN32)
#define TS_LIVE 1
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <pthread.h>
#include <signal.h>
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#else
#define TS_LIVE 0
#endif

/* ---- pure decisions -------------------------------------------------------- */

bool term_status_decide(const TermStatusEnv* e) {
    if (!e || !e->interactive || !e->threads || !e->ttys) return false;
    if (!e->term || !*e->term || strcmp(e->term, "dumb") == 0) return false;
    if (e->spinner && strcmp(e->spinner, "0") == 0) return false;
    return true;
}

static bool has_utf8(const char* s) {
    char low[64];
    size_t n = 0;
    for (; s[n] && n + 1 < sizeof low; n++) low[n] = (char)tolower((unsigned char)s[n]);
    low[n] = '\0';
    return strstr(low, "utf-8") || strstr(low, "utf8");
}

bool term_status_locale_utf8(const char* lc_all, const char* lc_ctype, const char* lang) {
    const char* v = (lc_all && *lc_all) ? lc_all
                  : (lc_ctype && *lc_ctype) ? lc_ctype
                  : (lang && *lang) ? lang : NULL;
    return v && has_utf8(v);
}

const char* term_status_label_of(const Expr* e) {
    for (int guard = 0; e && guard < 8; guard++) {
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return NULL;
        const char* h = e->data.function.head->data.symbol.name;
        size_t n = e->data.function.arg_count;
        if ((h == SYM_Set || h == SYM_SetDelayed) && n == 2) { e = e->data.function.args[1]; continue; }
        /* Monitor[expr, mon] shows mon in the line; the label is expr's */
        if (h == SYM_Monitor && n == 2) { e = e->data.function.args[0]; continue; }
        if (h == SYM_CompoundExpression && n > 0) {
            const Expr* pick = NULL;
            for (size_t i = n; i > 0; i--) {
                const Expr* a = e->data.function.args[i - 1];
                if (a->type == EXPR_FUNCTION) { pick = a; break; }
            }
            e = pick;
            continue;
        }
        return h;
    }
    return NULL;
}

/* ---- shared state ---------------------------------------------------------- */

static bool g_statusline = false;
static bool g_unicode = false;
static bool g_stats = false;         /* begin/end measure (interactive session)  */
static int  g_nest = 0;              /* begin/end nesting                         */

/* per evaluation, main thread */
static double   g_t0 = 0, g_c0 = 0;
static uint64_t g_peak0 = 0;
static bool     g_have_peak0 = false;

static double cpu_now(void) {
    clock_t c = clock();
    return c == (clock_t)-1 ? -1.0 : (double)c / (double)CLOCKS_PER_SEC;
}

bool term_status_unicode(void) { return g_unicode; }
bool term_status_statusline(void) { return g_statusline; }
void term_status_set_statusline(bool on) { g_statusline = on; }

void term_status_print_statusline(double last) {
    if (!g_statusline) return;
    uint64_t mem = 0;
    if (!meminfo_current(&mem)) mem = 0;
    char buf[512];
    TermSpinnerStyle st = { term_color_stdout(), g_unicode };
    term_status_line(buf, sizeof buf, last, mem, context_current(), st);
    fputs(buf, stdout);
    fflush(stdout);
}

#if TS_LIVE && !defined(_WIN32)   /* TS_LIVE implies it; spelled out for check-c99 */

/* ---- the forwarder --------------------------------------------------------- */

static bool g_color = false;
static bool g_enabled = false;       /* decided (or opened) and not switched off */

enum { ST_IDLE, ST_RUN, ST_STOP, ST_QUIT };

static struct {
    bool open;
    int tty;                    /* where the thread writes (a dup of the terminal) */
    int route[2], n_route;      /* descriptors pointed at the pipe                 */
    int saved[2];               /* their real targets while redirected             */
    int echo_fd;                /* termios target, -1 = leave alone                */
    int pr, pw;                 /* data pipe                                       */
    int cr, cw;                 /* control pipe (wakes poll)                       */
    int delay_ms, period_ms, cols_fallback;
    bool thread_started;
    int state;
    bool redirected;            /* fds currently point at the pipe                 */
    bool suspended;
    /* parameters of the current run (written under mu before ST_RUN) */
    const char* label;
    double run_t0, run_c0;
    /* result of the current run (written by the thread before ST_IDLE) */
    uint64_t sampled_peak;
} S = { false, -1, { -1, -1 }, 0, { -1, -1 }, -1, -1, -1, -1, -1, 200, 100, 80,
        false, ST_IDLE, false, false, NULL, 0, 0, 0 };
static pthread_t g_th;
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t g_cv = PTHREAD_COND_INITIALIZER;

/* Signal-handler view: plain sig_atomic_t flags and a copy of the termios. */
static volatile sig_atomic_t g_sig_active = 0;    /* the line may be on screen */
static volatile sig_atomic_t g_sig_echo_off = 0;  /* termios must be restored  */
static volatile sig_atomic_t g_sig_tty = -1;
static volatile sig_atomic_t g_sig_color = 0;
static volatile sig_atomic_t g_sig_exiting = 0;   /* an exit path owns the tty now */
static volatile sig_atomic_t g_sig_width = 0;     /* cells of the visible line     */
static volatile sig_atomic_t g_sig_stopped = 0;   /* Ctrl-Z: stopped, do not draw  */
static const char SPACES[256 + 1] =
    "                                                                "
    "                                                                "
    "                                                                "
    "                                                                ";
static struct termios g_tio_saved;

static void write_all(int fd, const char* p, size_t n) {
    while (n > 0) {
        long k = (long)write(fd, p, n);
        if (k < 0) {
            if (errno == EINTR) continue;
            return;
        }
        p += k;
        n -= (size_t)k;
    }
}

static int term_cols(void) {
    struct winsize ws;
    if (ioctl(S.tty, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0) return (int)ws.ws_col;
    return S.cols_fallback;
}

/* The erase sequence for a line `w` cells wide: EL with colour, spaces without
 * (NO_COLOR promises no escape sequences at all). */
static void erase_line(int w) {
    if (g_color) { write_all(S.tty, "\r\033[2K", 5); return; }
    char sp[512];
    if (w > (int)sizeof sp - 2) w = (int)sizeof sp - 2;
    sp[0] = '\r';
    memset(sp + 1, ' ', (size_t)w);
    sp[w + 1] = '\r';
    write_all(S.tty, sp, (size_t)w + 2);
}

typedef struct {
    bool visible;
    int width;
    bool col0;
    double last_draw;
    unsigned frame;
    long int_seen;
    TermEta eta;
    uint64_t budget;
    double quiet_until;   /* after a Ctrl-C: let the abort unwind undrawn */
} Fwd;

static void draw(Fwd* f, double now, const char* label, double t0, double c0) {
    if (g_sig_exiting || g_sig_stopped) return;
    TermSpinnerSnap s;
    memset(&s, 0, sizeof s);
    s.frame = f->frame++;
    s.label = label;
    s.live = eval_dispatch_head();
    s.elapsed = now - t0;
    double c = cpu_now();
    s.cpu = (c >= 0 && c0 >= 0) ? c - c0 : -1;
    MthProgressView o, in;
    int depth = term_progress_sample(&o, &in);
    s.frac = -1;
    s.eta = -1;
    if (depth > 0 && o.head) {
        s.outer.head = o.head; s.outer.done = o.done; s.outer.total = o.total;
        if (depth >= 2 && in.head) { s.inner.head = in.head; s.inner.done = in.done; s.inner.total = in.total; }
        s.frac = term_spinner_fraction(&s.outer, &s.inner, depth);
        s.eta = term_eta_update(&f->eta, o.gen, s.frac, now);
    } else {
        f->eta.have = false;
    }
    uint64_t m = 0;
    if (meminfo_current(&m)) {
        s.mem = m;
        if (m > S.sampled_peak) S.sampled_peak = m;
    }
    s.peak = S.sampled_peak;
    uint64_t lim = mth_progress_memory_limit();
    s.budget = lim ? lim : f->budget;
    s.workers = mth_progress_workers();
    char mon[MTH_MONITOR_TEXT];
    if (mth_monitor_active()) {
        g_mth_monitor_due = 1;               /* the main thread refreshes it */
        if (mth_monitor_text(mon, sizeof mon)) s.monitor = mon;
    }

    char line[1024];
    TermSpinnerStyle st = { g_color, g_unicode };
    int w = term_spinner_render(&s, term_cols(), st, line, sizeof line);
    char out[2200];
    size_t n = 0;
    if (g_color) {
        n = (size_t)snprintf(out, sizeof out, "\r\033[2K%s", line);
    } else {
        /* no escapes: overwrite, blank the rest of a longer previous line,
         * then write the line again so the cursor sits at its end */
        int pad = f->visible && f->width > w ? f->width - w : 0;
        if (pad) n = (size_t)snprintf(out, sizeof out, "\r%s%*s\r%s", line, pad, "", line);
        else     n = (size_t)snprintf(out, sizeof out, "\r%s", line);
    }
    if (n >= sizeof out) n = sizeof out - 1;
    write_all(S.tty, out, n);
    f->visible = true;
    f->width = w;
    g_sig_width = w;
    f->last_draw = now;
    g_sig_active = 1;
}

static void hide(Fwd* f) {
    if (!f->visible) return;
    erase_line(f->width);
    f->visible = false;
    g_sig_active = 0;
}

/* Copy one chunk to the terminal, with the "^C" the disabled echo would have
 * shown, before the handler's newline. */
static void forward(Fwd* f, const char* buf, size_t n) {
    if (g_sig_exiting) return;
    hide(f);
    long ic = term_interrupt_eval_count();
    if (ic != f->int_seen) {
        f->int_seen = ic;
        if (S.echo_fd >= 0) write_all(S.tty, "^C", 2);
        /* an abort usually unwinds in milliseconds; a frame drawn meanwhile
         * would only flash between "^C" and the result */
        f->quiet_until = mth_mono_seconds() + 0.5;
    }
    write_all(S.tty, buf, n);
    f->col0 = buf[n - 1] == '\n';
}

static void drain_nonblocking(Fwd* f) {
    char buf[8192];
    for (;;) {
        long k = (long)read(S.pr, buf, sizeof buf);
        if (k > 0) { forward(f, buf, (size_t)k); continue; }
        if (k < 0 && errno == EINTR) continue;
        break;   /* EAGAIN: empty */
    }
}

static void run_one(void) {
    pthread_mutex_lock(&g_mu);
    const char* label = S.label;
    double t0 = S.run_t0, c0 = S.run_c0;
    pthread_mutex_unlock(&g_mu);

    Fwd f;
    memset(&f, 0, sizeof f);
    f.col0 = true;
    f.int_seen = term_interrupt_eval_count();
    uint64_t phys = 0;
    f.budget = meminfo_physical(&phys) ? phys : 0;
    const double delay = S.delay_ms / 1000.0, period = S.period_ms / 1000.0;

    const double min_gap = 0.05;   /* redraw at most 20 times a second after output */
    for (;;) {
        double now = mth_mono_seconds();
        double wait;
        if (now < t0 + delay)  wait = t0 + delay - now;
        else if (!f.col0)      wait = period;            /* a partial line: never overwrite */
        else if (!f.visible)   wait = f.last_draw + min_gap - now;
        else                   wait = f.last_draw + period - now;
        int timeout = (int)(wait * 1000.0 + 0.5);
        if (timeout < 0) timeout = 0;
        if (timeout > S.period_ms) timeout = S.period_ms;
        struct pollfd p[2];
        p[0].fd = S.pr; p[0].events = POLLIN; p[0].revents = 0;
        p[1].fd = S.cr; p[1].events = POLLIN; p[1].revents = 0;
        int r = poll(p, 2, timeout);
        if (r < 0 && errno != EINTR) break;
        if (r > 0 && (p[0].revents & (POLLIN | POLLHUP))) {
            char buf[8192];
            long k = (long)read(S.pr, buf, sizeof buf);
            if (k > 0) forward(&f, buf, (size_t)k);
        }
        if (r > 0 && (p[1].revents & POLLIN)) {
            char c;
            while (read(S.cr, &c, 1) == 1) { }
            pthread_mutex_lock(&g_mu);
            int st = S.state;
            pthread_mutex_unlock(&g_mu);
            if (st != ST_RUN) {
                /* the main thread has restored fd 1/2: everything it wrote is
                 * already in the pipe */
                drain_nonblocking(&f);
                hide(&f);
                long ic = term_interrupt_eval_count();
                if (ic != f.int_seen && S.echo_fd >= 0) write_all(S.tty, "^C\n", 3);
                break;
            }
        }
        now = mth_mono_seconds();
        if (now - t0 >= delay && f.col0 && now >= f.quiet_until &&
            (f.visible ? now - f.last_draw >= period * 0.9 : now - f.last_draw >= min_gap))
            draw(&f, now, label, t0, c0);
    }
    g_sig_active = 0;
}

static void* thread_main(void* arg) {
    (void)arg;
    pthread_mutex_lock(&g_mu);
    for (;;) {
        while (S.state == ST_IDLE) pthread_cond_wait(&g_cv, &g_mu);
        if (S.state == ST_QUIT) break;
        S.sampled_peak = 0;
        pthread_mutex_unlock(&g_mu);
        run_one();
        pthread_mutex_lock(&g_mu);
        if (S.state == ST_QUIT) break;
        S.state = ST_IDLE;
        pthread_cond_broadcast(&g_cv);
    }
    S.state = ST_IDLE;
    pthread_cond_broadcast(&g_cv);
    pthread_mutex_unlock(&g_mu);
    return NULL;
}

static bool start_thread(void) {
    if (S.thread_started) return true;
    sigset_t all, old;
    sigfillset(&all);
    pthread_sigmask(SIG_SETMASK, &all, &old);
    int rc = pthread_create(&g_th, NULL, thread_main, NULL);
    pthread_sigmask(SIG_SETMASK, &old, NULL);
    S.thread_started = rc == 0;
    return S.thread_started;
}

/* ---- terminal state the cleanup paths restore ------------------------------ */

static void echo_off(void) {
    if (S.echo_fd < 0) return;
    struct termios t;
    if (tcgetattr(S.echo_fd, &g_tio_saved) != 0) return;
    t = g_tio_saved;
    if (!(t.c_lflag & ECHO)) return;          /* already off: leave it */
    t.c_lflag &= ~(tcflag_t)ECHO;
    if (tcsetattr(S.echo_fd, TCSANOW, &t) == 0) g_sig_echo_off = 1;
}

static void echo_restore(void) {
    if (!g_sig_echo_off) return;
    tcsetattr(S.echo_fd, TCSANOW, &g_tio_saved);
    g_sig_echo_off = 0;
}

/* Async-signal-safe: termios, the line, and whatever output is still queued. */
static void sig_cleanup(void) {
    int tty = (int)g_sig_tty;
    if (g_sig_echo_off) { tcsetattr(S.echo_fd, TCSANOW, &g_tio_saved); g_sig_echo_off = 0; }
    if (tty < 0) return;
    if (g_sig_active) {
        if (g_sig_color) { if (write(tty, "\r\033[2K", 5) < 0) { /* best effort */ } }
        else {
            int w = (int)g_sig_width;
            if (w > 256) w = 256;
            if (write(tty, "\r", 1) < 0 || write(tty, SPACES, (size_t)w) < 0 ||
                write(tty, "\r", 1) < 0) { /* best effort */ }
        }
        g_sig_active = 0;
    }
}

/* The double Ctrl-C exit (term_interrupt_set_exit_fn). */
static void status_exit_fn(int status) {
    g_sig_exiting = 1;
    if (S.redirected && g_sig_tty >= 0) {
        sig_cleanup();
        if (write((int)g_sig_tty, "^C\n", 3) < 0) { /* exiting */ }
    }
    _Exit(status);
}

/* Fatal signals: restore the terminal, flush what is queued, re-raise. */
static const int FATAL[] = { SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT, SIGTERM, SIGHUP, SIGQUIT };
#define N_FATAL ((int)(sizeof FATAL / sizeof FATAL[0]))
static struct sigaction g_fatal_old[N_FATAL];
static bool g_fatal_installed[N_FATAL];

static void on_fatal(int sig) {
    g_sig_exiting = 1;
    if (S.redirected) {
        sig_cleanup();
        /* forward what the pipe still holds (read/write are signal-safe) */
        int fl = fcntl(S.pr, F_GETFL);
        if (fl >= 0) fcntl(S.pr, F_SETFL, fl | O_NONBLOCK);
        char buf[4096];
        for (int i = 0; i < 64; i++) {
            long k = (long)read(S.pr, buf, sizeof buf);
            if (k <= 0) break;
            if (write((int)g_sig_tty, buf, (size_t)k) < 0) break;
        }
        for (int i = 0; i < S.n_route; i++)
            if (S.saved[i] >= 0) dup2(S.saved[i], S.route[i]);
    } else {
        sig_cleanup();
    }
    for (int i = 0; i < N_FATAL; i++)
        if (FATAL[i] == sig && g_fatal_installed[i]) sigaction(sig, &g_fatal_old[i], NULL);
    raise(sig);
}

static void install_fatal(void) {
    for (int i = 0; i < N_FATAL; i++) {
        struct sigaction cur;
        if (sigaction(FATAL[i], NULL, &cur) != 0) continue;
        if (cur.sa_handler != SIG_DFL || (cur.sa_flags & SA_SIGINFO)) continue;  /* someone else's */
        struct sigaction sa;
        memset(&sa, 0, sizeof sa);
        sa.sa_handler = on_fatal;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = SA_RESETHAND;
        if (sigaction(FATAL[i], &sa, &g_fatal_old[i]) == 0) g_fatal_installed[i] = true;
    }
}

/* Ctrl-Z during an evaluation (ISIG is on, so the tty sends SIGTSTP to the
 * main thread): erase the line and give the shell back its echoing terminal,
 * stop with the default action, and on SIGCONT put ECHO-off back; the
 * forwarder redraws on its next tick. Outside an evaluation (the editor
 * raises SIGTSTP itself after restoring the terminal) it is a plain stop. */
static struct sigaction g_tstp_old;
static bool g_tstp_installed = false;

static void on_tstp(int sig) {
    int saved_errno = errno;
    bool echo_was_off = g_sig_echo_off != 0;
    g_sig_stopped = 1;
    sig_cleanup();
    struct sigaction dfl, ours;
    memset(&dfl, 0, sizeof dfl);
    dfl.sa_handler = SIG_DFL;
    sigemptyset(&dfl.sa_mask);
    sigaction(SIGTSTP, &dfl, &ours);
    sigset_t m;
    sigemptyset(&m);
    sigaddset(&m, SIGTSTP);
    raise(sig);                                  /* pending: blocked in here */
    pthread_sigmask(SIG_UNBLOCK, &m, NULL);      /* delivered: stopped       */
    pthread_sigmask(SIG_BLOCK, &m, NULL);        /* continued                */
    sigaction(SIGTSTP, &ours, NULL);
    if (echo_was_off && S.echo_fd >= 0) {
        struct termios t = g_tio_saved;
        t.c_lflag &= ~(tcflag_t)ECHO;
        if (tcsetattr(S.echo_fd, TCSANOW, &t) == 0) g_sig_echo_off = 1;
    }
    g_sig_stopped = 0;
    errno = saved_errno;
}

static void install_tstp(void) {
    struct sigaction cur;
    if (sigaction(SIGTSTP, NULL, &cur) != 0) return;
    if (cur.sa_handler != SIG_DFL || (cur.sa_flags & SA_SIGINFO)) return;   /* someone else's */
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_tstp;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    if (sigaction(SIGTSTP, &sa, &g_tstp_old) == 0) g_tstp_installed = true;
}

static void at_exit_cleanup(void) {
    if (g_nest > 0 && S.redirected) {
        TermStatusStats st;
        g_nest = 1;
        term_status_end(&st);
    }
    term_status_shutdown();
}

/* ---- opening ------------------------------------------------------------- */

static int cloexec_pipe(int fds[2]) {
    if (pipe(fds) != 0) return -1;
    for (int i = 0; i < 2; i++) {
        int fl = fcntl(fds[i], F_GETFD);
        if (fl >= 0) fcntl(fds[i], F_SETFD, fl | FD_CLOEXEC);
    }
    return 0;
}

bool term_status_open(const TermStatusConfig* c) {
    if (S.open || !c || c->tty_fd < 0) return S.open;
    int dp[2], cp[2];
    if (cloexec_pipe(dp) != 0) return false;
    if (cloexec_pipe(cp) != 0) { close(dp[0]); close(dp[1]); return false; }
    /* the write end must NOT be close-on-exec once dup2'd onto 1/2 -- dup2
     * clears FD_CLOEXEC on the target, so children of Run[] still write
     * into the pipe; the originals stay private */
    int fl = fcntl(dp[0], F_GETFL);
    if (fl >= 0) fcntl(dp[0], F_SETFL, fl | O_NONBLOCK);
    fl = fcntl(cp[0], F_GETFL);
    if (fl >= 0) fcntl(cp[0], F_SETFL, fl | O_NONBLOCK);
    S.pr = dp[0]; S.pw = dp[1];
    S.cr = cp[0]; S.cw = cp[1];
    S.tty = c->tty_fd;
    S.n_route = c->n_route < 0 ? 0 : (c->n_route > 2 ? 2 : c->n_route);
    for (int i = 0; i < S.n_route; i++) S.route[i] = c->route[i];
    S.echo_fd = c->echo_fd;
    S.delay_ms = c->delay_ms < 0 ? 0 : c->delay_ms;
    S.period_ms = c->period_ms <= 0 ? 100 : c->period_ms;
    S.cols_fallback = c->cols > 0 ? c->cols : 80;
    g_color = c->color;
    g_unicode = c->unicode;
    g_sig_tty = S.tty;
    g_sig_color = c->color ? 1 : 0;
    S.open = true;
    g_enabled = true;
    g_stats = true;
    return true;
}

static bool same_terminal(int a, int b) {
    struct stat sa, sb;
    if (fstat(a, &sa) != 0 || fstat(b, &sb) != 0) return false;
    return sa.st_rdev == sb.st_rdev && sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
}

void term_status_init(bool interactive) {
    g_stats = interactive;
    g_unicode = term_status_locale_utf8(getenv("LC_ALL"), getenv("LC_CTYPE"), getenv("LANG"));
    const char* sl = getenv("MATHILDA_STATUSLINE");
    g_statusline = interactive && sl && strcmp(sl, "1") == 0;
    if (interactive) term_progress_init();
    TermStatusEnv e;
    e.interactive = interactive;
    e.threads = true;
    e.ttys = isatty(0) && isatty(1) && isatty(2) && same_terminal(1, 2);
    e.term = getenv("TERM");
    e.spinner = getenv("MATHILDA_SPINNER");
    if (!term_status_decide(&e)) return;

    TermStatusConfig c;
    memset(&c, 0, sizeof c);
    c.tty_fd = dup(1);
    if (c.tty_fd < 0) return;
    int fl = fcntl(c.tty_fd, F_GETFD);
    if (fl >= 0) fcntl(c.tty_fd, F_SETFD, fl | FD_CLOEXEC);
    c.route[0] = 1; c.route[1] = 2; c.n_route = 2;
    c.echo_fd = 0;
    c.color = term_color_stdout();
    c.unicode = g_unicode;
    c.delay_ms = 200;
    const char* d = getenv("MATHILDA_SPINNER_DELAY_MS");
    if (d && *d) {
        long v = strtol(d, NULL, 10);
        if (v >= 0 && v <= 10000) c.delay_ms = (int)v;
    }
    c.period_ms = 100;
    const char* col = getenv("COLUMNS");
    c.cols = col ? atoi(col) : 0;
    if (!term_status_open(&c)) { close(c.tty_fd); return; }
    term_interrupt_set_exit_fn(status_exit_fn);
    install_fatal();
    install_tstp();
    atexit(at_exit_cleanup);
}

/* ---- begin / end ---------------------------------------------------------- */

static void redirect_on(const char* label, double t0, double c0) {
    if (!start_thread()) { g_enabled = false; return; }
    fflush(stdout);
    fflush(stderr);
    for (int i = 0; i < S.n_route; i++) {
        S.saved[i] = dup(S.route[i]);
        if (S.saved[i] >= 0) {
            int fl = fcntl(S.saved[i], F_GETFD);
            if (fl >= 0) fcntl(S.saved[i], F_SETFD, fl | FD_CLOEXEC);
        }
        dup2(S.pw, S.route[i]);
    }
    S.redirected = true;
    echo_off();
    pthread_mutex_lock(&g_mu);
    S.label = label;
    S.run_t0 = t0;
    S.run_c0 = c0;
    S.state = ST_RUN;
    pthread_cond_broadcast(&g_cv);
    pthread_mutex_unlock(&g_mu);
}

static uint64_t redirect_off(void) {
    if (!S.redirected) return 0;
    echo_restore();
    fflush(stdout);
    fflush(stderr);
    for (int i = 0; i < S.n_route; i++) {
        if (S.saved[i] >= 0) { dup2(S.saved[i], S.route[i]); close(S.saved[i]); S.saved[i] = -1; }
    }
    S.redirected = false;
    pthread_mutex_lock(&g_mu);
    S.state = ST_STOP;
    if (write(S.cw, "s", 1) < 0) { /* the pipe has room: one byte per stop */ }
    while (S.state != ST_IDLE) pthread_cond_wait(&g_cv, &g_mu);
    uint64_t peak = S.sampled_peak;
    pthread_mutex_unlock(&g_mu);
    return peak;
}

void term_status_begin(const char* label) {
    if (g_nest++ > 0) return;
    g_t0 = mth_mono_seconds();
    g_c0 = g_stats ? cpu_now() : -1;
    g_have_peak0 = g_stats && meminfo_peak(&g_peak0);
    if (g_enabled && S.open && !S.suspended) redirect_on(label, g_t0, g_c0);
}

void term_status_end(TermStatusStats* out) {
    if (g_nest <= 0) { g_nest = 0; if (out) memset(out, 0, sizeof *out); return; }
    if (--g_nest > 0) return;
    uint64_t sampled = redirect_off();
    S.suspended = false;
    if (!out) return;
    out->elapsed = mth_mono_seconds() - g_t0;
    double c = g_c0 >= 0 ? cpu_now() : -1;
    out->cpu = (c >= 0 && g_c0 >= 0) ? c - g_c0 : -1;
    uint64_t pk = 0;
    out->peak = sampled;
    if (g_have_peak0 && meminfo_peak(&pk) && pk > g_peak0 && pk > out->peak) out->peak = pk;
}

void term_status_suspend(void) {
    if (g_nest <= 0 || !S.redirected) return;
    redirect_off();
    S.suspended = true;
}

void term_status_resume(void) {
    if (g_nest <= 0 || !S.suspended) return;
    S.suspended = false;
    if (g_enabled && S.open) redirect_on("Dialog", mth_mono_seconds(), cpu_now());
}

bool term_status_enabled(void) { return g_enabled && S.open; }

void term_status_set_enabled(bool on) {
    if (g_nest > 0) return;              /* never mid-evaluation */
    g_enabled = on && S.open;
}

void term_status_shutdown(void) {
    if (!S.open) return;
    if (S.redirected) redirect_off();
    if (S.thread_started) {
        pthread_mutex_lock(&g_mu);
        S.state = ST_QUIT;
        pthread_cond_broadcast(&g_cv);
        pthread_mutex_unlock(&g_mu);
        pthread_join(g_th, NULL);
        S.thread_started = false;
    }
    close(S.pr); close(S.pw); close(S.cr); close(S.cw);
    S.pr = S.pw = S.cr = S.cw = -1;
    if (S.tty > 2) close(S.tty);
    S.tty = -1;
    g_sig_tty = -1;
    S.open = false;
    g_enabled = false;
    term_interrupt_set_exit_fn(NULL);
    if (g_tstp_installed) { sigaction(SIGTSTP, &g_tstp_old, NULL); g_tstp_installed = false; }
}

#else  /* !TS_LIVE: no thread -- begin/end only measure */

void term_status_init(bool interactive) {
    g_stats = interactive;
    g_unicode = term_status_locale_utf8(getenv("LC_ALL"), getenv("LC_CTYPE"), getenv("LANG"));
    const char* sl = getenv("MATHILDA_STATUSLINE");
    g_statusline = interactive && sl && strcmp(sl, "1") == 0;
}

bool term_status_open(const TermStatusConfig* cfg) { (void)cfg; return false; }

void term_status_begin(const char* label) {
    (void)label;
    if (g_nest++ > 0) return;
    g_t0 = mth_mono_seconds();
    g_c0 = g_stats ? cpu_now() : -1;
    g_have_peak0 = g_stats && meminfo_peak(&g_peak0);
}

void term_status_end(TermStatusStats* out) {
    if (g_nest <= 0) { g_nest = 0; if (out) memset(out, 0, sizeof *out); return; }
    if (--g_nest > 0) return;
    if (!out) return;
    out->elapsed = mth_mono_seconds() - g_t0;
    double c = g_c0 >= 0 ? cpu_now() : -1;
    out->cpu = (c >= 0 && g_c0 >= 0) ? c - g_c0 : -1;
    uint64_t pk = 0;
    out->peak = (g_have_peak0 && meminfo_peak(&pk) && pk > g_peak0) ? pk : 0;
}

void term_status_suspend(void) { }
void term_status_resume(void) { }
bool term_status_enabled(void) { return false; }
void term_status_set_enabled(bool on) { (void)on; }
void term_status_shutdown(void) { }

#endif /* TS_LIVE */
