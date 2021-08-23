/* dynamic_sll.c: contains the implementation of singly linked list
 * wrapper functions that use dynamic memory allocations for the 
 * linked list pointers.
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

#include <stdlib.h>

#include "sll.h"
#include "dynamic_sll.h"


SLinkedList *alloc_list_init( SLinkedList *sll, void (*free_func)(void *) )
{
    sll = malloc( sizeof(*sll) );
    if( NULL != sll ) {
        sll_init( sll, free_func );
    }
    return sll;
}

void dealloc_list( SLinkedList *s )
{
    if( NULL != s ) {
        sll_destroy( s );
        free( s );
    }
    return;
}

