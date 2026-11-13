/*
 * term_config.h -- the interactive REPL's per-user directory and config file.
 *
 * The directory is $MATHILDA_CONFIG_DIR, else $HOME/.mathilda; an empty
 * MATHILDA_CONFIG_DIR turns the whole mechanism off (no config, no init.m,
 * no mods), which is what tests and scripted sessions use. It holds:
 *
 *   config      "key = value" lines; '#' starts a comment. Keys today: theme.
 *   init.m      Mathilda code evaluated at the start of an interactive session
 *   mods/NAME.m the same, one file per mod, in name order, after init.m
 *
 * Only the interactive session reads any of it; pipe mode and -file never do.
 * Pure C99 except term_config_set, which creates the directory (mkdir, POSIX;
 * on Windows _mkdir) and replaces the file through a temporary + rename.
 */
#ifndef MATHILDA_TERM_CONFIG_H
#define MATHILDA_TERM_CONFIG_H

#include <stdbool.h>
#include <stddef.h>

/* The directory (malloc'd, no trailing slash), or NULL when disabled or
 * HOME is unset. Pure: the caller passes the environment. */
char* term_config_dir(const char* env_dir, const char* env_home);

/* dir + "/" + name, malloc'd; NULL when dir is NULL. */
char* term_config_join(const char* dir, const char* name);

/* The value of `key` in the file at `path` into buf (trimmed). False when the
 * file or the key is missing. */
bool  term_config_get(const char* path, const char* key, char* buf, size_t cap);

/* Set `key` to `value` in the file at `path`: the first "key = ..." line is
 * replaced (later duplicates dropped), other lines and comments are kept, a
 * missing key is appended. Creates the file and its directory. */
bool  term_config_set(const char* path, const char* key, const char* value);

/* The same edit on text in memory (malloc'd result), for the tests. */
char* term_config_edit(const char* text, const char* key, const char* value);

#endif /* MATHILDA_TERM_CONFIG_H */
