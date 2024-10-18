#ifndef _XMALLOC_C
#define _XMALLOC_C

#include <stdlib.h>

#include "error.c"

static void *xmalloc(size_t size)
{
	if (size == 0)
		error("xmalloc: zero size");
	void *new = malloc(size);
	if (!new)
	{
		double sm = size / (0x100000 * 1.0);
		error("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return new;
}

//static void *xrealloc(void *p, size_t s)
//{
//	void *r = realloc(p, s);
//	if (!r) error("realloc failed");
//	return r;
//}

static void xfree(void *p)
{
	if (!p)
		error("thou shalt not free a null pointer!");
	free(p);
}

#endif//_XMALLOC_C
