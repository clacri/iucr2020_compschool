/* dupstr.c: contains the implementation of string duplication
 * utilities.
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
#include <stdlib.h>
#include <string.h>

#include "dupstr.h"

char *dupnstr( const char *s, size_t len )
{
    char *ret = NULL;

    /* start with some sanity checks */
    if( (NULL == s) || (0 == len)  )
        return NULL;

    ret = malloc( len + 1 );
    if( !ret )
        return NULL;

    strncpy( ret, s, len);
    ret[len] = '\0';
    return ret;
}  /* end dupnstr() */

char *dupstr( const char *s )
{
    size_t len = strlen( s );

    return dupnstr( s, len );
}  /* end dupstr() */


char *dupsubstr( const char *s, size_t offset, size_t len )
{
    const char *start = NULL;
    size_t s_len = strlen( s );

    if( (NULL == s ) ||
        ( offset > s_len) ||
        ((offset + len) > s_len) ) 
        return NULL;

     start = s + offset;
     return dupnstr( start, len );
} /* end dupsubstr() */

char *duppsubstr( const char *s, char *start, size_t len )
{
   size_t offset;
   if( NULL == start )
       return NULL;
   
   offset = start - s;
   return dupsubstr( s, offset, len );
} /* end duppsubstr() */
