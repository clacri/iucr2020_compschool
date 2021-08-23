/* contains the 'main' module for the Flack left coset decomposition
 * program which uses alogorithms outlined in Acta Cryst. (1987), A43,
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

#define _ISOC99_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "input.h"
#include "queue.h"
#include "task.h"

#ifdef PYTHON_EXTENSION_MODULE
#include <Python.h>
#define MSG_BUF_SZ 256

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif
#endif


void usage( void );


#ifdef PYTHON_EXTENSION_MODULE
static PyObject *decomp( PyObject *self, PyObject *args )
#else  /* just go with an ANSI declaration for main() */
int main( int argc, char **argv )
#endif
{
    Queue *task_queue;
    void *task_data;
    struct task *t;
    char *filename;
    int n_tasks = 0;
#ifdef PYTHON_EXTENSION_MODULE
    char msg[MSG_BUF_SZ] = {0};

    if( ! PyArg_Parse( args, "(s)", &filename ) ) {
        return NULL;
    }
#else
    if( argc == 2 ) {
        filename = argv[argc-1];
    }
    else {
        usage();
        exit( EXIT_FAILURE );
    }
#endif
    
    errno = 0;
    task_queue = read_input_file( filename );
    if( NULL == task_queue ) {
        fprintf( stderr, "task_queue pointer is NULL, most recent error: %s",
                 errno != 0 ? strerror(errno) : "was not recorded by errno" );
        exit(EXIT_FAILURE );
    }

    n_tasks = queue_size( task_queue );

    while( 0 < queue_size(task_queue) ) {
           queue_dequeue( task_queue, &task_data );
           t = (struct task *)task_data;
           process_task( t );
           dealloc_task( t );
    }

    queue_destroy( task_queue );
    free( task_queue );

#ifdef PYTHON_EXTENSION_MODULE
    snprintf( msg, MSG_BUF_SZ, "Program processed %d tasks input from file %s", n_tasks, filename ); 
    return Py_BuildValue( "s", msg );
#else
    fprintf( stdout, "Program processed %d tasks input from file %s\n", n_tasks, filename ); 
    fflush( stdout );
    exit( EXIT_SUCCESS );
#endif
}
#ifdef PYTHON_EXTENSION_MODULE
struct PyMethodDef flack_coset_methods[] = {
                   { "decomp", decomp, METH_VARARGS, "Takes input file as only argument" },
                   { NULL, NULL, 0, NULL}
                   };


PyMODINIT_FUNC initFlackCoset( void )
{
   Py_InitModule( "FlackCoset", flack_coset_methods );
}
#endif
