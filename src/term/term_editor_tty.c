/*
 * term_editor_tty.c -- see term_editor_tty.h.
 */
/* termios, poll, sigaction and TIOCGWINSZ are POSIX / BSD, hidden by glibc
 * under -std=c99 and partly by Darwin under a strict _XOPEN_SOURCE (the
 * src/meminfo.c pattern: each platform gets the macro that widens its own
 * namespace, before any include). */
#if defined(__APPLE__)
#  ifndef _DARWIN_C_SOURCE
#  define _DARWIN_C_SOURCE
#  endif
#elif !defined(_WIN32)
#  ifndef _DEFAULT_SOURCE
#  define _DEFAULT_SOURCE
#  endif
#  ifndef _XOPEN_SOURCE
#  define _XOPEN_SOURCE 700
#  endif
#endif

#include "term_editor_tty.h"
#include "term_editor.h"
#include "term_history.h"
#include "term_highlight.h"

#include <stdlib.h>
#include <string.h>

static bool g_session_active = false;

bool term_editor_session_active(void) { return g_session_active; }

bool term_editor_select(const char* env_editor, const char* env_term,
                        bool stdin_tty, bool stdout_tty) {
    if (!stdin_tty || !stdout_tty) return false;
    if (!env_term || !*env_term || strcmp(env_term, "dumb") == 0) return false;
    if (env_editor && *env_editor) {
        static const char* const OFF[] = { "readline", "libedit", "plain", "none", "off", "0" };
        for (size_t i = 0; i < sizeof(OFF) / sizeof(OFF[0]); i++)
            if (strcmp(env_editor, OFF[i]) == 0) return false;
    }
    return true;
}

#ifdef _WIN32

TermReader* term_reader_open_editor(const TermEditorOptions* o) {
    (void)o;
    return NULL;
}

#else /* POSIX */

#include <errno.h>
#include <poll.h>
#include <signal.h>
#include <stdio.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

/* ---- terminal state shared with the signal and exit paths ---------------- */

static struct termios        g_saved;          /* the termios of the read's start */
static volatile sig_atomic_t g_raw = 0;        /* raw mode is on                  */
static volatile sig_atomic_t g_winch = 0;      /* SIGWINCH arrived                */
static bool                  g_atexit = false;

static const int CRASH_SIGS[] = { SIGTERM, SIGHUP, SIGQUIT, SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGABRT };
#define N_CRASH ((int)(sizeof(CRASH_SIGS) / sizeof(CRASH_SIGS[0])))
static struct sigaction g_old_crash[N_CRASH];
static struct sigaction g_old_winch;
static bool             g_handlers = false;

static const char PASTE_ON[]  = "\033[?2004h";
static const char PASTE_OFF[] = "\033[?2004l";

/* Async-signal-safe: write() and tcsetattr() only. */
static void restore_terminal(void) {
    if (!g_raw) return;
    g_raw = 0;
    if (write(STDOUT_FILENO, PASTE_OFF, sizeof(PASTE_OFF) - 1) < 0) { /* nothing to do */ }
    tcsetattr(STDIN_FILENO, TCSANOW, &g_saved);
}

static void at_exit_restore(void) { restore_terminal(); }

static void on_crash(int sig) {
    restore_terminal();
    for (int i = 0; i < N_CRASH; i++) {
        if (CRASH_SIGS[i] == sig) { sigaction(sig, &g_old_crash[i], NULL); break; }
    }
    raise(sig);   /* delivered with the previous disposition once we return */
}

static void on_winch(int sig) {
    (void)sig;
    g_winch = 1;
}

static void install_handlers(void) {
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sigemptyset(&sa.sa_mask);
    sa.sa_handler = on_crash;
    for (int i = 0; i < N_CRASH; i++) sigaction(CRASH_SIGS[i], &sa, &g_old_crash[i]);
    sa.sa_handler = on_winch;
    sa.sa_flags = 0;                  /* no SA_RESTART: poll() returns EINTR */
    sigaction(SIGWINCH, &sa, &g_old_winch);
    g_handlers = true;
}

static void remove_handlers(void) {
    if (!g_handlers) return;
    for (int i = 0; i < N_CRASH; i++) sigaction(CRASH_SIGS[i], &g_old_crash[i], NULL);
    sigaction(SIGWINCH, &g_old_winch, NULL);
    g_handlers = false;
}

static bool enter_raw(void) {
    if (tcgetattr(STDIN_FILENO, &g_saved) != 0) return false;
    struct termios t = g_saved;
    t.c_iflag &= ~(tcflag_t)(BRKINT | ICRNL | INLCR | IGNCR | INPCK | ISTRIP | IXON);
    t.c_lflag &= ~(tcflag_t)(ECHO | ICANON | IEXTEN | ISIG);
    t.c_cflag |= CS8;
    t.c_cc[VMIN] = 1;
    t.c_cc[VTIME] = 0;
    /* TCSANOW, not TCSAFLUSH: keys typed during the evaluation are kept. */
    if (tcsetattr(STDIN_FILENO, TCSANOW, &t) != 0) return false;
    g_raw = 1;
    return true;
}

static void leave_raw(void) {
    restore_terminal();
}

static void write_all(const char* s, size_t n) {
    while (n > 0) {
        long w = (long)write(STDOUT_FILENO, s, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            return;
        }
        s += w;
        n -= (size_t)w;
    }
}

static void tty_write(void* ctx, const char* s, size_t n) {
    (void)ctx;
    write_all(s, n);
}

static void term_size(int* cols, int* rows) {
    struct winsize ws;
    *cols = 80;
    *rows = 24;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0) {
        if (ws.ws_col > 0) *cols = ws.ws_col;
        if (ws.ws_row > 0) *rows = ws.ws_row;
    }
}

/* ---- the reader ---------------------------------------------------------- */

typedef struct {
    TermEditorOptions opt;
    TermEditor*       ed;
    char*             pre;          /* bytes read while waiting for a DSR reply */
    size_t            pre_len, pre_cap;
    int               dsr_failures; /* after two, stop asking                  */
} EdImpl;

static void pre_add(EdImpl* im, const char* s, size_t n) {
    if (im->pre_len + n > im->pre_cap) {
        size_t nc = im->pre_cap ? im->pre_cap * 2 : 256;
        while (nc < im->pre_len + n) nc *= 2;
        char* np = realloc(im->pre, nc);
        if (!np) return;
        im->pre = np;
        im->pre_cap = nc;
    }
    memcpy(im->pre + im->pre_len, s, n);
    im->pre_len += n;
}

static long now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long)ts.tv_sec * 1000L + ts.tv_nsec / 1000000L;
}

/* Ask the terminal where the cursor is (ESC [ 6 n). Bytes that are not the
 * reply (typeahead) are kept for the editor. Returns false on no reply. */
static bool query_cursor(EdImpl* im, int* row, int* col) {
    *row = -1;
    *col = 0;
    if (im->dsr_failures >= 2) return false;
    write_all("\033[6n", 4);
    char buf[256];
    size_t n = 0;
    long deadline = now_ms() + 80;
    for (;;) {
        long left = deadline - now_ms();
        if (left <= 0) break;
        struct pollfd pfd = { STDIN_FILENO, POLLIN, 0 };
        int pr = poll(&pfd, 1, (int)left);
        if (pr < 0 && errno == EINTR) continue;
        if (pr <= 0) break;
        long r = (long)read(STDIN_FILENO, buf + n, sizeof(buf) - 1 - n);
        if (r <= 0) break;
        n += (size_t)r;
        /* Look for ESC [ rows ; cols R. */
        for (size_t i = 0; i + 1 < n; i++) {
            if (buf[i] != '\033' || buf[i + 1] != '[') continue;
            size_t j = i + 2;
            int a = 0, b = 0;
            bool da = false, db = false;
            while (j < n && buf[j] >= '0' && buf[j] <= '9') { a = a * 10 + (buf[j++] - '0'); da = true; }
            if (j >= n || buf[j] != ';') continue;
            j++;
            while (j < n && buf[j] >= '0' && buf[j] <= '9') { b = b * 10 + (buf[j++] - '0'); db = true; }
            if (j >= n || buf[j] != 'R' || !da || !db) continue;
            pre_add(im, buf, i);
            pre_add(im, buf + j + 1, n - j - 1);
            *row = a - 1;
            *col = b - 1;
            im->dsr_failures = 0;
            return true;
        }
        if (n >= sizeof(buf) - 1) break;
    }
    pre_add(im, buf, n);
    im->dsr_failures++;
    return false;
}

static TermReadStatus ed_read(TermReader* rd, const char* prompt, char** out) {
    EdImpl* im = (EdImpl*)rd->impl;
    fflush(stdout);
    fflush(stderr);
    if (!enter_raw()) return TERM_READ_EOF;
    install_handlers();
    write_all(PASTE_ON, sizeof(PASTE_ON) - 1);

    int cols, rows;
    term_size(&cols, &rows);
    term_editor_set_size(im->ed, cols, rows);
    int row, col;
    im->pre_len = 0;
    (void)query_cursor(im, &row, &col);
    term_editor_begin(im->ed, prompt, row, col);
    TermEdStatus st = term_editor_feed(im->ed, im->pre, im->pre_len);
    im->pre_len = 0;

    char buf[4096];
    bool hangup = false;
    while (st == TERM_ED_CONTINUE) {
        if (g_winch) {
            g_winch = 0;
            term_size(&cols, &rows);
            term_editor_resize(im->ed, cols, rows);
        }
        struct pollfd pfd = { STDIN_FILENO, POLLIN, 0 };
        int pr = poll(&pfd, 1, term_editor_pending(im->ed) ? 50 : -1);
        if (pr < 0) {
            if (errno == EINTR) continue;
            hangup = true;
            break;
        }
        if (pr == 0) { st = term_editor_timeout(im->ed); continue; }
        long r = (long)read(STDIN_FILENO, buf, sizeof(buf));
        if (r < 0 && (errno == EINTR || errno == EAGAIN)) continue;
        if (r <= 0) { hangup = true; break; }
        st = term_editor_feed(im->ed, buf, (size_t)r);
        if (st == TERM_ED_SUSPEND) {
            leave_raw();
            remove_handlers();
            raise(SIGTSTP);                 /* stops here until fg */
            if (!enter_raw()) { hangup = true; break; }
            install_handlers();
            write_all(PASTE_ON, sizeof(PASTE_ON) - 1);
            term_size(&cols, &rows);
            term_editor_set_size(im->ed, cols, rows);
            term_editor_redraw(im->ed, false);
            st = TERM_ED_CONTINUE;
        }
    }
    leave_raw();                    /* also switches bracketed paste off */
    remove_handlers();

    if (hangup) return TERM_READ_EOF;
    if (st == TERM_ED_EOF) return TERM_READ_EOF;
    if (st == TERM_ED_INTERRUPT) return TERM_READ_INTERRUPTED;
    size_t n;
    const char* text = term_editor_text(im->ed, &n);
    char* copy = malloc(n + 1);
    if (!copy) return TERM_READ_EOF;
    memcpy(copy, text, n + 1);
    *out = copy;
    return TERM_READ_OK;
}

static void ed_history_add(TermReader* rd, const char* entry) {
    EdImpl* im = (EdImpl*)rd->impl;
    if (im->opt.history) (void)term_history_add(im->opt.history, entry);
}

static void ed_close(TermReader* rd) {
    EdImpl* im = (EdImpl*)rd->impl;
    restore_terminal();
    remove_handlers();
    if (im) {
        term_editor_free(im->ed);
        free(im->pre);
        free(im);
    }
    g_session_active = false;
}

static const TermReaderOps ED_OPS = { ed_read, ed_history_add, ed_close };

TermReader* term_reader_open_editor(const TermEditorOptions* o) {
    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO)) return NULL;
    struct termios probe;
    if (tcgetattr(STDIN_FILENO, &probe) != 0) return NULL;
    TermReader* rd = malloc(sizeof(*rd));
    EdImpl* im = calloc(1, sizeof(*im));
    if (!rd || !im) { free(rd); free(im); return NULL; }
    if (o) im->opt = *o;
    int cols, rows;
    term_size(&cols, &rows);
    TermEditorConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.cols = cols;
    cfg.rows = rows;
    cfg.color = im->opt.color;
    cfg.history = im->opt.history;
    cfg.classify = term_hl_classify_symtab;
    cfg.write = tty_write;
    im->ed = term_editor_new(&cfg);
    if (!im->ed) { free(rd); free(im); return NULL; }
    if (!g_atexit) { atexit(at_exit_restore); g_atexit = true; }
    rd->ops = &ED_OPS;
    rd->kind = TERM_BACKEND_EDITOR;
    rd->impl = im;
    g_session_active = true;
    return rd;
}

#endif /* POSIX */
