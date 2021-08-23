/* read input interface for Flack 'coset' program using a Finite
 *  State Machine (FSM).
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

#ifndef INPUT_H
#define INPUT_H

#include "queue.h"
#include "task.h"

#define INPUT_LINE_LEN 80
#define ERR_MSG_LEN    72
#define NIBBLE_LEN 4

typedef struct fsm {
        unsigned int flags;
        int line_num;
        struct fsm *(*next)(struct fsm * );
        char *input_filename;
        FILE *inp;
        char line[INPUT_LINE_LEN];
        int rmats_read;
        int last_err;
        char err_msg[ERR_MSG_LEN];
        struct task *tsk;
        Queue *task_queue;
        } fsm;


struct keyword_table {
       char *keywd;
       fsm *(*func)( struct fsm *f );
       };


/* function prototypes */
Queue *read_input_file( char *fname );
void fsm_init( struct fsm *f );
Queue *fsm_pass_task_queue( struct fsm *f );

#endif
