/* contains the 'task' public interface for the Flack left coset
 * decomposition program which uses alogorithms outlined in
 * Acta Cryst. (1987), A43, 564-568, by H. D. Flack.
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

#define _ISOC99_SOURCE
#ifndef TASK_H
#define TASK_H

#include <stdio.h>

#include "symm_mat.h"


#define GROUP_NAME_LEN 6

struct task {
       void (*coset_decomp)(struct symm_op *, struct symm_op * );
       char *title;
       struct symm_op *super;
       struct symm_op *sub;
       char super_name[GROUP_NAME_LEN];
       char sub_name[GROUP_NAME_LEN];
       char algorithm_name;
       int n_subgroup_mats;
       double trans_mat[3][3];
       unsigned int trans_mat_bcm;
       char *outfile;
       char *shelx_ins_file;
       char *new_base_name;
       char *shelx_executable;
       };

void init_task( struct task *t );
void dealloc_task( void *task );
void process_task( struct task *t );
#endif

