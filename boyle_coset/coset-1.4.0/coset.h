/* contains interface used for the Flack left coset decomposition
 * algorithms as outlined in Acta Cryst. (1987), A43, 564-568, by H. D.
 * Flack.
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

#ifndef COSET_H
#define COSET_H
#include <stdio.h>

#include "symm_mat.h"

/* prototypes */
void coset_decomposition_A( struct symm_op *G, struct symm_op *H );
void coset_decomposition_B( struct symm_op *G, struct symm_op *H );
void output_cosets( FILE *out, struct symm_op *s, int start );

#endif
