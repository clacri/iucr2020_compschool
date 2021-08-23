/* implementation file for the interfaces defined in double_util.h 
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

#define _ISOC99_SOURCE
#include <math.h>

#include "float_util.h"

double relative_difference( double a, double b )
{
    double c = ABS(a);
    double d = ABS(b);

    d = MAX( c, d );
    return d == 0.0 ? 0.0 : ABS( a - b) / d;
} 

/* returns a 1 if 'a' is essentially zero and 0 if it is
 * non-zero.  This is weird.
 */
int is_zero( double a )
{
   return relative_difference( a, 0.0 ) <= TOLERANCE ? 1 : 0;
} 

/* returns a 1 if 'a' and 'b' are equal within TOLERANCE
 * and zero if they are non-equal.
 */

int is_equal( double a, double b )
{
    double x;
    x = a - b;
    return is_zero( x );
}

/* takes a double and returns an int which has been rounded to the nearest integer */

int round_to_nearest_int( double f )
{
    int ret;

    if ( f > 0.0 )
         ret = (int)( f + 0.500 );
    else
         ret = (int)( f - 0.500 );

    return ret;
}

/* the next 2 functions were added by Paul Boyle in November 2017. */

/* signof():  returns -1 for negative numbers and a +1 for positive numbers */
int signof( double a )
{
    return signbit( a ) ? -1 : 1;  /*signbit() is C99 macro */
}

/* match_signs(): returns 1 if signs of two doubles match or 0 if they do not */
int match_signs( double a, double b )
{
    return signof(a) == signof(b) ? 1 : 0;  
}
