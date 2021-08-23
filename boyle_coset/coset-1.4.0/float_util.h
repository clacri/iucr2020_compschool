/* some useful macros and functions for dealing with double and floating point
 * and other sorts of "mathy" operations in C.
 * 
 * See the comp.lang.c FAQ question 14.5 for information regarding
 * the ABS(), MAX() macros as well as the relative_difference() 
 * function.
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
#ifndef FLOAT_UTIL_H
#define FLOAT_UTIL_H


#define ABS(x)     ((x < (double)0.0 ? -(x) : (x) ))
#define MAX(a, b)  ((a) > (b) ? (a) : (b))


#if ! defined  _XOPEN_SOURCE && ! defined  __GSL_MATH_H__

/* value of PI (denoted as M_PI in <math.h>)
 * from http://www.chemie.fu-berlin.de/chemnet/use/info/libc/libc_13.html
 */

#define M_PI  3.14159265358979323846264338327
#endif


#define TOLERANCE 0.0000001  /* good enough for most crystallographic applications */
      

/* convert deg -> rad and back with these macros */
#define DEG2RAD(deg) ( (M_PI/180.) * (deg) )
#define RAD2DEG(rad) ( (180./M_PI) * (rad) )

/* here is a function like macro for seeing of a value is between two others */
#define is_between(v,low,high) ((v) > (low) ? ((v) < (high) ? 1 : 0 ): 0 )


double relative_difference( double a, double b );
/* is_zero() returns 1 if tested quantity is 0 and 0 if it is not zero (huh?) :-) */
int is_zero( double a );

/* is_equal() returns 1 if 'a' and 'b' are equal and 0 if they are not */
int is_equal( double a, double b );

/* rounds a double to the nearest integer. */
int round_to_nearest_int( double f );

/* signof():  returns -1 for negative numbers and a +1 for positive numbers */
int signof( double a );

/* match_signs(): returns 1 if signs of two doubles match or 0 if they do not */
int match_signs( double a, double b );
#endif
