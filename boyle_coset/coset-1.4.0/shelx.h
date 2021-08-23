/* public interface for functions used for reading, editing, and writing
 * SHELX .ins files for the Flack left coset decomposition program
 * which uses alogorithms outlined in Acta Cryst. (1987), A43,
 * 564-568, by H. D. Flack.
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

#ifndef SHELX_OUTPUT_H
#define SHELX_OUTPUT_H

#define _ISOC99_SOURCE

#define TWIN_INS_BUF_LEN 160
#define BASF_SCALE_FORMAT " %4.2f"
#define BASF_BUF_LEN      6
#define TWIN_INS_FORMAT "TWIN %7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f %d\n"
#define SHELX_LINE_LEN  80

#include "symm_mat.h"
#include "sll.h"
#include "dynamic_sll.h"

/* prototypes */
struct symm_op *pick_twin_laws( struct symm_op *s );
SLinkedList *twin_ins_list( struct symm_op *s );
SLinkedList *read_shelx_ins_file( char *ins_file_name );
char *get_basename(char *filename, int delim_char);
SLinkedList *write_new_ins_files( char *base_name, SLinkedList *twin_laws, SLinkedList *ins );

#endif
