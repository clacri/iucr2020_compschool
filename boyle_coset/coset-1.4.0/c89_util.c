/* contains implementation for C89 compatibility functions
 * for compilers which don't yet support ANSI C99.
 *
 * Copyright (C) 2008, 2012, 2013 Paul D. Boyle
 *
 * Written by:
 *      Paul D. Boyle
 *      Department of Chemistry
 *      University of Western Ontario
 *      London, Ontario CANADA N6A 5B7
 *
 *      February 2008
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version
 * 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#undef _ISOC99_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <errno.h>

#ifdef NEED_C89_COMPATIBILITY
#include "c89_util.h"


/* TMPBUF_LEN should be long enough for must uses, but
 * be prepared to change it if there are truncation problems.
 */
#define TMPBUF_LEN 256 

/* C89 compatible version of isblank() */
int isblank( int c )
{
    int ret = 0;

    if( (c == ' ') || (c == '\t') )
        ret = 1;

    return ret;
}

/* my implementation of the C99 cbrt() "cube root" function */

double cbrt( double x )
{
   return pow( x, (1.0/3.0) );
}


/* my own version of vsnprintf() which is just a wrapper around
 * vsprintf().  Uses a dynamically allocated temporary buffer
 * (tmpbuf) with a relatively large memory allocation to catch all
 * the conversions (hopefully). Then the function strncpy()'s 
 * the contents of the buffer to the 'str' parameter up to 'size'-1
 * bytes (leaving room for the terminating NUL character.  The dynamic
 * memory allocation is completely internal to the function, so the caller
 * does not have to free() the memory; the function itself cleans up 
 * after itself.  The problem with the this function is if the formatted
 * output exceeds TMPBUF_LEN, the data will be truncated.
 *
 */

int vsnprintf(char *str, size_t size, const char *format, va_list ap)
{
    int c = 0,
        ret = 0;
    char *tmpbuf;

    errno = 0;
    tmpbuf = malloc( (size_t)TMPBUF_LEN );
    if( NULL == tmpbuf )
        return -1;

    memset( tmpbuf, 0, (size_t)TMPBUF_LEN );
    memset( str, 0, size );
    c = vsprintf( tmpbuf, format, ap );

    strncpy( str, tmpbuf, size );
    str[size-1] = '\0';

    if( c >= (int)size ) {
        ret = c;
    }
    else {
        ret = (int)strlen( str );
    }

    free( tmpbuf );
    return ret;
}

/* my own version of snprintf() which is just a wrapper around
 * my version of vsnprintf().
 */

int snprintf( char *str, size_t size, const char *format, ... )
{
    int ret = 0;
    va_list ap;

    memset( str, 0, size );
    va_start( ap, format );
    ret = vsnprintf( str, size, format, ap );
    va_end( ap );

    return ret;
}


#undef TMPBUF_LEN
#endif
