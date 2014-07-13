#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"

/* Implementation without temporary buffer when shift != 0 */
void memxchgshift(void *buf1, void *buf2, size_t count, int shift)
{
	assert(shift);

	if(shift > 0) {
		buf1+= count;
		buf2+= count;

		while(count--) {
			*((char*) (buf1 + shift)) = *((char*) buf2);
			*((char*) (buf2 + shift)) = *((char*) buf1);
			buf1--; buf2--;
		}
	} else {
		shift = -shift;
		while(count--) {
			*((char*) buf1) = *((char*) (buf2 + shift));
			*((char*) buf2) = *((char*) (buf1 + shift));
			buf1++; buf2++;
		}
	}

	memset(buf1, 0x0, shift);
	memset(buf2, 0x0, shift);
}

void memshift(void *buf, size_t size, int shift)
{
	if(shift > 0) {
		memmove(buf + shift, buf, size - shift);
		memset(buf, 0x0, shift);
	} else {
		shift = abs(shift);
		memmove(buf, buf + shift, size - shift);
		memset(buf + size - shift, 0x0, shift);
	}
}

