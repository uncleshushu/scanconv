#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>
#include <assert.h>

#ifndef NDEBUG
#define DBG_PRINT(FORMAT, ...) fprintf(stderr,"[DEBUG] %s(%d), `%s`: "FORMAT, __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define DBG_PRINT(FORMAT, ...)
#endif /* #ifndef NDEBUG */

#endif /* #ifndef DEBUG_H */