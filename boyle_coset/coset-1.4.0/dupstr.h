/* dupstr.h: contains the public interface for some string duplication
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
 * Note:
 * these functions all call malloc(), and the caller must free the
 * memory allocated by the functions.  These functions return NULL
 * on any error.
 */
#ifndef DUPSTR_H
#define DUPSTR_H

/* dupstr() -- duplicates a string */
char *dupstr( const char *s );

/* dupnstr() -- duplicates the first 'len' bytes of a string */
char *dupnstr( const char *s, size_t len );

/* dupsubstr() -- duplicates a substring of 'len' bytes of an input
 * string starting from an offset from the start of the buffer.
 */

char *dupsubstr( const char *s, size_t off_set, size_t len );

/* duppsubstr() -- duplicates a substring of 'len' bytes of an
 * input string starting at a specified pointer.
 */

char *duppsubstr( const char *s, char *start, size_t len );
#endif
