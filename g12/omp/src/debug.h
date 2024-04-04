// Macros to help us debug stuff
// Since we can only print the final result to stdout, we'll be using the stderr for debug messages
#pragma once

#if MSG_LEVEL >= 3
#define debug(...) fprintf(stderr, "[DEBG] " __VA_ARGS__)
#else
#define debug(...)
#endif

#if MSG_LEVEL >= 2
#define info(...) fprintf(stderr, "[INFO] " __VA_ARGS__)
#else
#define info(...)
#endif

#if MSG_LEVEL >= 1
#define warn(...) fprintf(stderr, "[WARN] " __VA_ARGS__)
#else
#define warn(...)
#endif

#define error(...) fprintf(stderr, "[ERROR] " __VA_ARGS__)
#define syserr(...) perror("[ERROR] " __VA_ARGS__)
