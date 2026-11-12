/*
 * repl_cmds.h -- the interactive REPL's slash commands, themes and mods.
 *
 * A line that starts with '/' is a command (term_slash.h), handled here and
 * never evaluated: it gets no In[n]/Out[n] and leaves $Line alone. Built-in
 * commands:
 *
 *   /help                    the commands and the keys
 *   /theme [name]            the theme picker (live preview; Return applies
 *                            and saves, Esc reverts), or set one by name
 *   /history [n]             the last n (20) inputs
 *   /timing [on|off|always]  the "took" note: for slow results (on, the
 *                            default), never, or after every result
 *   /statusline [on|off]     the status line after each result
 *   /clear                   clear the screen
 *   /mods                    the mods loaded, their commands and hooks
 *   /quit                    end the session, as Quit[] does
 *
 * Mods: at start-up the session evaluates <dir>/init.m and then every
 * <dir>/mods/NAME.m in name order (<dir> from term_config.h; MATHILDA_MODS=0
 * skips them). A mod is ordinary Mathilda code: it can define functions, set
 * the session hooks $PreRead / $Pre / $Post / $PrePrint, and add commands
 * with SlashCommand["name", f, "summary"]. A mod that fails -- a message, an
 * abort, a Quit[] -- is reported once and the session carries on; a failing
 * command prints its messages and returns to the prompt.
 *
 * Only the interactive session (editor or readline) uses any of this; pipe
 * mode, -file and the stdio loop of a NO_READLINE build never call it.
 */
#ifndef MATHILDA_REPL_CMDS_H
#define MATHILDA_REPL_CMDS_H

#include <stdbool.h>

#include "term/term_input.h"

struct TermHistory;

/* Register the built-in commands and remember the session's reader, history
 * and per-user directory (term_config_dir; NULL for none, copied). */
void repl_cmds_init(TermReader* rd, struct TermHistory* hist, const char* config_dir);
bool repl_cmds_enabled(void);

/* Load init.m and the mods (after the banner). */
void repl_cmds_load_mods(void);

/* Run the command line `line` (line[0] == '/'). */
void repl_cmds_run(const char* line);

/* Should the "took" note follow a result that took `seconds`? */
bool repl_cmds_show_took(double seconds);

void repl_cmds_shutdown(void);

#endif /* MATHILDA_REPL_CMDS_H */
