#ifndef UTIL_H
#define UTIL_H

#include <stddef.h>

void memxchgshift(void *buf1, void *buf2, size_t count, int shift);

void memshift(void *buf, size_t size, int shift);

#endif // UTIL_H
