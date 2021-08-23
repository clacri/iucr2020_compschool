/* contains public interface for C89 compatibility typedefs
 * and functions for compilers which don't yet support ANSI C99.
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

#ifndef C89_UTIL_H
#define C89_UTIL_H
#undef _ISOC99_SOURCE

#include <stdarg.h>
/* typedef a bool so C99's <stdbool.h> doesn't need to be included. */
typedef unsigned int bool;

/* prototypes for C99 functions */

int isblank( int c );
int snprintf(char *str, size_t size, const char *format, ...);
int vsnprintf(char *str, size_t size, const char *format, va_list ap);
double cbrt( double x );

#endif


