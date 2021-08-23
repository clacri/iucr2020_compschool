/* contains the 'task' implementation module for the Flack left coset
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#ifdef USE_SYSTEM_FUNCTION
#define USE_NONSTANDARD_FOPEN  /* for nonconforming C implementations of fopen() */
#endif

#include "matrix.h"
#include "symm_mat.h"
#include "coset.h"
#include "task.h"
#include "version.h"


#include "sll.h"
#include "dynamic_sll.h"
#include "shelx.h"
#include "shelx_exec.h"

#define TIME_FORMAT "%d %b %Y at %H:%M:%S"
#define TIME_BUFLEN 24  /* this should be long enough to take TIME_FORMAT and a NUL */

static void get_time( char *outbuf, size_t len )
{
    time_t now;
    size_t sz;
    struct tm *bdt; /* broken down time */

    now = time( NULL );
    bdt = localtime( &now ); 

    sz = strftime( outbuf, len, TIME_FORMAT, bdt );

    if( 0 == sz ) {
        fprintf( stderr, "%s:%d: strftime() returned zero bytes.", __FILE__, __LINE__ );
    }
    return;
}


static void print_task_header( FILE *out, struct task *t )
{
   
    char timebuf[TIME_BUFLEN] = { 0 };
    size_t sz = TIME_BUFLEN;

    get_time( timebuf, sz);
    fprintf( out, "COSET Decomposition Program (version %s - %s eigen code) run on: %s\n\n", VERSION, EIGEN_CODE, timebuf );
    fprintf( out, "Task Description: %s\n", t->title );
    fprintf( out, "Metrically Available Supergroup's Symmetry: %s\n", t->super_name );
    fprintf( out, "Crystal's Pointgroup (Subgroup): %s (%s)\n", t->sub_name,
                   is_centric(t->sub) ? "centric" : "acentric" );
    fprintf( out, "Flack Algorithm: %c\n", t->algorithm_name );
    fprintf( out, "Matrix which transforms Subgroup's Lattice to Supergroup's Lattice:\n" );
    print_matrix( out, t->trans_mat, "%8.4f%8.4f%8.4f\n" );
    fputc( '\n', out );
    if( NULL != t->shelx_ins_file ) 
        fprintf( out, "Original SHELX .ins file: %s\n", t->shelx_ins_file );
    if( NULL != t->new_base_name )
        fprintf( out, "New SHELX .ins files to be created with this basename: %s\n",
                       t->new_base_name );

    if( NULL != t->shelx_executable )
        fprintf( out, "SHELXL refinements on trial twin laws will done\nwith the executable file: %s\n",
                       t->shelx_executable );

    fputc( '\n', out );
    return;
}


void init_task( struct task *t )
{
    int i, j;

    t->coset_decomp = NULL;
    t->title = NULL;
    t->super = NULL;
    t->sub = NULL;
    t->n_subgroup_mats = 0;

/* as a default the indentity matrix is used for the transformation matrix */
    for( i = 0; i < 3; i++ ) {
         for( j = 0; j < 3; j++ ) {
              if( i == j ) {
                  t->trans_mat[i][j] = 1.0;
              }
              else {
                 t->trans_mat[i][j] = 0.0;
              }
         }
    }

    t->trans_mat_bcm = encode_matrix( t->trans_mat );
    t->outfile = NULL;
    t->shelx_ins_file = NULL;
    t->new_base_name = NULL;
    t->shelx_executable = NULL;

    return;
}

void dealloc_task( void *task )
{
    struct task *t;

    if( NULL == task )
        return;

    t = (struct task *)task;

    if( NULL != t->title )
        free( t->title );

    if( NULL != t->super )
        free( t->super );

    if( NULL != t->sub )
        free( t->sub );

    if( NULL != t->outfile )
        free( t->outfile );

    if( NULL != t->shelx_ins_file)
        free( t->shelx_ins_file);

    if( NULL != t->new_base_name )
        free( t->new_base_name );

    if( NULL != t->shelx_executable )
        free( t->shelx_executable );

    free( t );
    return;
}


void process_task( struct task *t )
{
    FILE *coset_out;
    int jobs_run = 0;
    SLinkedList *orig_ins_file = NULL, 
                *twin_shelx_instr = NULL,
                *new_ins_file_list = NULL, 
                *job_list = NULL;

    struct symm_op *duped = NULL;
    double det;
    double inverted_trans_mat[3][3] = { { 0.0 } };
#ifdef USE_NONSTANDARD_FOPEN
    const char *mode = "at";
#else  /* use only ANSI C Standard flags for mode */
    const char *mode = "a";
#endif

/* give the user a little information that his/her program is
 * actually running.
 */
   fprintf( stdout, "Processing Task: %s ...\n", t->title );

/* if an OUTFILE is not specified, write the results to stdout */
    if( NULL == t->outfile ) {
        coset_out = stdout; 
    }
    else {
       errno = 0;
       coset_out = fopen( t->outfile, mode );
       if( NULL == coset_out ) {
           fprintf( stderr, "%s: %s\n", t->outfile,
                     errno != 0 ? strerror(errno) : "couldn't open file." );
           fputs( "Writing results to stdout.\n", stderr );
           coset_out = stdout;
       }
    }

/* invert the transformation matrix to prepare for transforming the 
 * system of representatives back to the crystal's lattice setting.
 */
    det = determinant( t->trans_mat );
    invert_matrix( det, t->trans_mat, inverted_trans_mat );

/* print out the input */
    print_task_header( coset_out, t );

/* save a copy of the original subgroup  and then transform the subgroup 
 * the subgroups' truefalse value is initialized to 'False', so set it to
 * 'True' so that print_2_symm_ops() will work.
 */

    set_truth_value( t->sub, True, 0 );
    duped = duplicate_ops( t->sub );
    if( NULL == duped ) {  /* symm_op duplication didn't work */
        fprintf( stderr, "%s:%d: %s\n", __FILE__, __LINE__, strerror(errno) );
    }
    transform_group( t->sub, t->trans_mat );
    if( NULL != duped ) {
        print_2_symm_ops( coset_out, "Subgroup Symmetry Matricies",
                          "Subgroup Symmetry Matrices Transformed to Supergroup's Lattice", duped, t->sub );
    }

    free( duped ); /* don't need it anymore */ 

/* do the actual coset decomposition here */
    if( NULL != t->coset_decomp ) {
        t->coset_decomp( t->super, t->sub );
    }
    else {  /* no algorithm  selected */
        fputs( "### End of COSET Output ###\n", coset_out );
        fclose( coset_out );
        return;
    }

/* duplicate the supergroup to prepare for printing out the untransformed 
 * and transformed potential twin laws.
 */
    duped = duplicate_ops( t->super );
    if( NULL == duped ) {  /* symm_op duplication didn't work */
        fprintf( stderr, "%s:%d: %s\n", __FILE__, __LINE__, strerror(errno) );
    }

    transform_group( t->super, inverted_trans_mat );
 
/* determine types of symmetry operators which are potential twin laws */
    analyze_symm_group( t->super );

/* print out potential twin laws */
    fputs("\n*** Potential Twin Laws for this Subgroup-Supergroup Relationship ***\n", coset_out );
    fputs( "Use matricies in right hand column for creating SHELX TWIN instructions.\n\n", coset_out );
    if( NULL != duped ) {
        print_2_symm_ops( coset_out, "Untransformed Supergroup Matricies",
                          "Transformed to Subgroup's Lattice",
                           duped, t->super );
    }

    free( duped );

/* read a SHELX .ins file if it has been specified. */
	    if( NULL != t->shelx_ins_file ) {
            errno = 0;
            orig_ins_file = read_shelx_ins_file( t->shelx_ins_file );
        if( NULL == orig_ins_file ) {
            fprintf( stderr,"%s:%d: %s\n", __FILE__,__LINE__,
            0 != errno ? strerror(errno) : "read_shelx_ins_file() returned NULL" );
            fclose( coset_out );
            return;
        }
    }
    else {  /* we are done */
       fputs( "### End of COSET Output ###\n", coset_out );
       fclose( coset_out );
       return;
    }

/* prepare the SHELX .ins files which contain the new TWIN instructions
 * for a subdequent least-squares job(s).
 */
    if( (NULL != t->new_base_name)  ) {
            twin_shelx_instr = twin_ins_list( t->super );
            if( NULL == twin_shelx_instr ) {
                fprintf( stderr,"%s:%d: %s\n", __FILE__,__LINE__,
                0 != errno ? strerror(errno) : "twin_ins_list() returned NULL" );
                fclose( coset_out );
                return;
            }
            new_ins_file_list = write_new_ins_files( t->new_base_name, 
                                twin_shelx_instr, orig_ins_file );
            if( NULL == new_ins_file_list ) {
                fprintf( stderr, "%s:%d: %s\n", __FILE__,__LINE__, 
                0 != errno ? strerror(errno) : "list of new .ins filename couldn't be written." );
                fclose( coset_out );
                dealloc_list( twin_shelx_instr );
                dealloc_list( twin_shelx_instr );
                dealloc_list( orig_ins_file );
                return;
            }
    }

    if( NULL != t->shelx_executable ) {
        errno = 0;
        job_list = setup_shelx_jobs(new_ins_file_list, t->shelx_ins_file );
        if( NULL != job_list ) {
            jobs_run = spawn_shelx_jobs( job_list, t->shelx_executable );
        }
        else {
           fprintf(stderr, "%s:%d: Couldn't setup SHELX jobs: %s\n", 
                   __FILE__, __LINE__, errno != 0 ? strerror(errno) : 
                   "setup_shelx_jobs() returned NULL, but errno not set." );
           return;
        }
        fprintf( coset_out, "%d SHELX jobs were run.  Please examine .res and .lst files.\n",
                            jobs_run );
    }

    dealloc_list( twin_shelx_instr );
    dealloc_list( new_ins_file_list );
    dealloc_list( job_list );
    dealloc_list( orig_ins_file );

    fputs( "### End of COSET Output ###\n", coset_out );
    fclose( coset_out );
    return;
}

