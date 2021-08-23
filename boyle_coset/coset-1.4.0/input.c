/* read input file using a Finite State Machine (FSM) for the 
 * Flack 'coset' program.
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
#ifdef NEED_C89_COMPATIBILITY
#undef _ISOC99_SOURCE
#else
#define _ISOC99_SOURCE
#endif


/* if we are compiling this to use the system() function, we guess that
 * we also need to use the nonstandard "t" mode in the fopen() statement.
 * (see the fopen() statement in the open_file() fsm function.
 */
#ifdef USE_SYSTEM_FUNCTION
#define USE_NONSTANDARD_FOPEN 
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "queue.h"
#include "symm_mat.h"
#include "task.h"
#include "coset.h"
#include "matrix.h"
#include "dupstr.h"
#include "input.h"

/* for compilers which don't have C99 functions available */
#ifdef NEED_C89_COMPATIBILITY
#include "c89_util.h"
#endif

#define COMMENT_CHAR '#'

#define DELIMITER '\n'  /* for get_line */

/* bit operations for 'flags' */
#define HAS_TITLE       (1 << 0)
#define HAS_ALGORITHM   (1 << 1)
#define HAS_SUPERGROUP  (1 << 2)
#define HAS_SUBGROUP    (1 << 3)
#define HAS_RMAT        (1 << 4)
#define HAS_TRANS       (1 << 5)
#define HAS_INSFILE     (1 << 6)
#define HAS_OUTFILE     (1 << 7)
#define HAS_EXEC        (1 << 8)
#define HAS_NEWINS      (1 << 9)
#define HAS_END         (1 << 10)
#define N_FLAGS         11
#define ALL_FLAGS       0x7ff /* 00000111 11111111 */

#define NEWINS_REQUIRES (HAS_INSFILE|HAS_TRANS)
#define EXEC_REQUIRES   (HAS_TRANS|HAS_NEWINS)


/* prototypes for the finite state machines state functions. */
static fsm *read_line( struct fsm *f );
static fsm *open_file( struct fsm *f );
static fsm *close_file( struct fsm *f );
static fsm *title( struct fsm *f );
static fsm *algorithm( struct fsm *f );
static fsm *supergroup( struct fsm *f );
static fsm *subgroup( struct fsm *f );
static fsm *rmat( struct fsm *f );
static fsm *trans( struct fsm *f );
static fsm *insfile( struct fsm *f );
static fsm *outfile( struct fsm *f );
static fsm *exec( struct fsm *f );
static fsm *newins( struct fsm *f );
static fsm *end( struct fsm *f );

/* Borrow get_line from my user_dialog utility, but just incorporate get_line()
 * as a static function.
*/

/* get_line() is a utility function which fills the 'line' variable with
 * a well formed C string.  If fgets() returns NULL, 'line' is assigned
 * a value of "".  The function returns the length of the string.
 */
static size_t get_line( FILE *f , char *line, size_t len )
{
    char  *ptr;
    int   c;

    if( !line || 0 == len ) {
        fputs( "Invalid buffer to hold input.\n", stderr );
        return 0;
    }

    ptr = fgets( line, (int)len, f );
    if( NULL == ptr ) {
        line[0] = '\0';
        return 0;
    }

    if( NULL != (ptr = strchr(line, DELIMITER)) )
        *ptr = '\0';
     else {
        line[len-1] = '\0';
        while( (c = fgetc(f)) != '\n' && c != EOF );
        fputs( "Warning: Input truncated and flushed.\n", stderr );
     }

     return strlen(line);
}


/* a little utility function to convert the entire string from
 * case to case.  the 'convert' pointer to function should be 
 * <ctype.h> 'toupper' or 'tolower' as appropriate.  's' *must*
 * be NUL terminated or mucho problems ensue.
 */

static void convert_string( char *s, int (*convert)(int) )
{
    char *p = s;
    while( '\0' != *p ) {
        *p = convert( *p );
         p++;
    }
}

/* skip_keyword(): finds the first blank in the line and advances the
 * pointer until the first non-blank character.
 */
static char *skip_keyword( char *s )
{
    char *p0 = strchr( s, ' ' );
    while( isblank( *++p0 ) );
    return p0;
}

/* for the directive which have only a single file name as an argument
 * we write a general function which gets wrapped by the directive specific
 * state functions.
 */
static char *get_filename( char *input_line )
{
   char *p;
   size_t len;

   p = skip_keyword( input_line );
   len = strlen( p );
   return duppsubstr( input_line, p, len );

}

/* gen_error_message(): generates an error message which can be printed
 * out at the appropriate time.
 */
static void gen_error_message(char *buf, size_t buflen, const char *fmt, ...  )
{
    va_list ap;

    va_start( ap, fmt );
    vsnprintf( buf, buflen, fmt, ap );
    va_end( ap );

    return;
}

#ifdef USE_SKIP_UNTIL
/* skip_until(): read the input file until keyword 'target' is reached, and then
 * 'targ_fun' is assigned to f->next.
 */

static fsm *skip_until( const char *target, struct fsm *f, struct fsm *(*targ_fun)(struct fsm *) )
{
    size_t len;
    char nibble[NIBBLE_LEN];

    f->next = NULL;
    while( NULL == f->next ) {
        len = get_line(f->inp, f->line, sizeof(f->line) );
        f->line_num++;
        strncpy( nibble, f->line, NIBBLE_LEN - 1 );
        nibble[NIBBLE_LEN-1] = '\0';
        convert_string( nibble, toupper );
        if( 0 == strncmp(target, nibble, NIBBLE_LEN-1) ) {
            f->next = targ_fun;
        }

    }
    return f;

}
#endif


/* these are the finite state machines's state functions */

static fsm *read_line( struct fsm *f )
{
    size_t len;
    char nibble[NIBBLE_LEN] = {0};  /* for future reference: nibble should never 
                                     * be made a pointer for dynamic memory allocation.
                                     * It is unnecessary, and it will throw off the
                                     * sizeof(nibble) expressions in this function.
                                     */


/* this is the table of keyword and associated state functions.  Should be
 * NULL terminated.
 */
    struct keyword_table tabl[] = {
           { "TITLE", title },
           { "ALGORITHM", algorithm },
           { "SUPERGROUP", supergroup },
           { "SUBGROUP", subgroup },
           { "RMAT", rmat },
           { "TRANS", trans },
           { "INSFILE", insfile },
           { "OUTFILE", outfile },
           { "EXEC", exec },
           { "NEWINS", newins },
           { "END", end },
           {  NULL, NULL }
           };

    struct keyword_table *tp;

/* read the line from the file and increment line counter, deal with pathological conditions */
    len = get_line(f->inp, f->line, sizeof(f->line) );
    if( 0 == len ) { /* there is some issue */
        if( feof(f->inp) ) {  /* someone probably forgot to put an 'END' statement at the end of their input file */
            fputs(  "??? Missing END statement at end of input file ???\n", stderr );
            f->next = close_file;
            return f;
        }
        else if( ferror(f->inp) ) {
            fprintf( stderr, "%s:%d: %s: %s\n", __FILE__,__LINE__, f->input_filename, 
                     0 != errno ? strerror(errno) : "ferror() was set" );
            exit( EXIT_FAILURE );
        }
    }
    f->line_num++;

/* check to see if it is a comment (line starts with '#' character). Skip
 * to next line.
 */
    if( COMMENT_CHAR == f->line[0] ) {
        f->next = read_line;
        return f;
    }

/* taste (nibble) the first 3 characters and convert to upper case to
 * see which directive is in the line.
 */
    strncpy( nibble, f->line, sizeof(nibble) - 1 );
    nibble[NIBBLE_LEN-1] = '\0';
    convert_string( nibble, toupper );

/* lookup the keyword and make appropriate assignments */
    for( tp = tabl; NULL != tp->func; tp++ ) {
         if( 0 == strncmp(tp->keywd, nibble, sizeof(nibble) - 1) ) {
             f->next = tp->func;
             break;
         } 
    }

/* if nothing matches, just read the next line.  We may want to 
 * change this in the future.
 */
    if( NULL == f->next ) {
        f->next = read_line;
    } 

    return f;
}


/* open_file(): opens the input file and initiializes the task queue */        

static fsm *open_file( struct fsm *f )
{
#ifdef USE_NONSTANDARD_FOPEN
   const char *mode = "rt";
#else  /* we use only strings defined by the ANSI C standard. */
   const char *mode = "r";
#endif
   errno = 0;
   f->inp = fopen( f->input_filename, mode );
   if( NULL == f->inp ) {
       char *etmp = errno != 0 ? strerror(errno) : "couldn't open input file";
       gen_error_message( f->err_msg, sizeof(f->err_msg), "%s:%d: %s",__FILE__, __LINE__,
                          etmp );
                 
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   errno = 0;
   f->task_queue = malloc( sizeof(*f->task_queue) );
   if( NULL == f->task_queue ) {
       char *etmp = errno != 0 ? strerror(errno): "could't allocate task_queue";
       gen_error_message( f->err_msg, sizeof(f->err_msg), "%s:%d: %s",__FILE__,__LINE__,
                          etmp );
       f->last_err = errno;
       f->next = NULL;
       fclose( f->inp );
   }
   else { /* we have memory allocated, so init the queue */
      queue_init( f->task_queue, dealloc_task );
      f->next = read_line;
   }  

   return f;
}

static fsm *close_file( struct fsm *f )
{
   fclose( f->inp );
   if( NULL != f->tsk ) { /* load the last task into the queue */
       queue_enqueue( f->task_queue, f->tsk );
       f->tsk = NULL;
   }

   f->last_err = 0;
   memset(f->err_msg, 0, sizeof(f->err_msg) );
   f->next = NULL;
   return f;
}


static fsm *title( struct fsm *f )
{
   size_t len;
   char *p;

   if( NULL != f->tsk ) { /* this is not the first task in the file. Take
                           * old task and place it on the queue.
                           */
       queue_enqueue( f->task_queue, f->tsk );
       f->tsk = NULL;
   }

/* for a new task, reset the FSM's 'flags' and 'rmats_read' members. */
   f->flags = 0;  
   f->rmats_read = 0;
   errno = 0;
   f->tsk = malloc( sizeof(*f->tsk) );
   if( NULL == f->tsk ) {
       char *etmp = errno != 0 ? strerror(errno) : "couldn't allocate struct task";
       gen_error_message( f->err_msg, sizeof(f->err_msg), "%s:%d: %s", __FILE__, __LINE__, etmp );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   init_task( f->tsk );

/* set the book keeping flag */
   f->flags |= HAS_TITLE;

/* skip over the TITLE keyword and advance pointer until first nonblank
 * line.  
 */

   p = skip_keyword( f->line );

/* allocate memory for title, dealing with error.  If everything goes OK, then
 * copy string pointed to by 'p0' to 'title'.
 */
   len = strlen( p );
   errno = 0;
   f->tsk->title = dupnstr( p, len );
   if( NULL == f->tsk->title ) {
       char *etmp = errno != 0 ? strerror(errno) : "couldn't allocate title buffer";
       gen_error_message( f->err_msg, sizeof(f->err_msg), "%s:%d: %s", __FILE__, __LINE__, etmp );
       f->last_err = errno;
       f->tsk->title = NULL;
   }
   f->next = read_line;
   return f;
}

static fsm *algorithm( struct fsm *f )
{
   char *p;

   p = skip_keyword(f->line );

   p[0] = toupper( p[0] );
              if( 'A' == p[0] ) {
                  f->tsk->coset_decomp = coset_decomposition_A;
              }
              else if( 'B' == p[0] ) {
                  f->tsk->coset_decomp = coset_decomposition_B;
              }
              else {
                 f->tsk->coset_decomp = NULL;
                 gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s",
                                   f->input_filename, f->line_num, "bad input line, algorithm not set" );
                 f->last_err = -3;
                 f->next = NULL;
                 return f;
              }
   f->tsk->algorithm_name = p[0];
   f->flags |= HAS_ALGORITHM;
   f->next = read_line;
   return f;
}
    
static fsm *supergroup( struct fsm *f )
{
   static const char *fmt = "%s%s";  /* format: SUPERGROUP <name> */
   char keyword[12];
   int n_scanned = 0,
       point_group_num,
       er;

   n_scanned = sscanf(f->line, fmt, keyword, f->tsk->super_name);
   if( 2 != n_scanned ) {  
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s",
                         f->input_filename, f->line_num, "bad input line, supergroup not set" );
       f->last_err = -3;
       f->next = NULL;
       return f;
   }

  
   point_group_num = lookup_supergroup( f->tsk->super_name );
   errno = 0;
   f->tsk->super = select_symm_ops( point_group_num, &er ); 
   if( NULL == f->tsk->super ) { /* error */
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "select_symm_ops() returned NULL" );
       f->last_err = er;
       f->next = NULL;
       return f;
   }

   f->flags |= HAS_SUPERGROUP;
   f->next = read_line;
   return f;
}    


static fsm *subgroup( struct fsm *f )
{
   static const char *fmt = "%s%s%d"; /* format: SUBGROUP <name> <number of matricies> */
   char keyword[12];
   int i,
       n_scanned,
       n_mat;


   n_scanned = sscanf( f->line, fmt, keyword, f->tsk->sub_name, &n_mat );
   if( 3 != n_scanned ) {
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s",
                         f->input_filename, f->line_num, "bad input line, subgroup not set" );
       f->last_err = -4;
       f->next = NULL;
       return f;
   }

   f->tsk->sub = malloc( (n_mat+1) * sizeof(*f->tsk->sub) );
   if( NULL == f->tsk->sub ) {
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "malloc() returned NULL" );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

/* zero bcm and truefalse members the struct symm_op array including the sentinel value at [n_mat] */
   for( i = 0; i <= n_mat; i++ ) {
        f->tsk->sub[i].bcm = 0;
        f->tsk->sub[i].truefalse = False;
   }
   f->tsk->n_subgroup_mats = n_mat;
   f->flags |= HAS_SUBGROUP;
   f->next = read_line;
   return f;
}

static fsm *rmat( struct fsm *f )
{
   static const char *fmt = "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf";
   int i, j, k, n_scanned;
   char keyword[12];
   double mat[3][3] = { {0.0} };
   double itm[3][3] = { {0.0} }; /* inverse transpose matrix */

   if( !(f->flags & HAS_SUBGROUP) ) {
       gen_error_message( f->err_msg, sizeof(f->err_msg), "%s:%d: input error: SUBGROUP must precede RMAT directive.",
                          f->input_filename, f->line_num );
       f->last_err = -5;
       f->next = NULL;
       return f; 
   }

/* start with some bookeeping.  we skip any RMAT lines after
 * the proper number of matricies have been read in.
 */
   if( f->rmats_read == f->tsk->n_subgroup_mats ) {
       f->next = read_line;
       return f;
   }

/* read the matrix elements from the line in the file to a temporary
 * 3x3 matrix and calculate it's inverse transpose.
 */  
    n_scanned = sscanf( f->line, fmt, keyword,
                       &mat[0][0], &mat[0][1], &mat[0][2],
                       &mat[1][0], &mat[1][1], &mat[1][2],
                       &mat[2][0], &mat[2][1], &mat[2][2] );
    if( 10 != n_scanned ) {
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s",
                         f->input_filename, f->line_num, "bad input line for RMAT" );
       f->last_err = -3;
       f->next = NULL;
       return f;
    }

    calculate_inverse_transpose( itm, mat );

   

/* search for first 'False' element, and read data into that one */
   for( k = 0; k < f->tsk->n_subgroup_mats; k++ ) {
        if( False == f->tsk->sub[k].truefalse ) {
            break; 
        }
   } 

/* copy the matrix */
   for( i = 0; i < 3; i++ ) {
        for( j = 0; j < 3; j++ ) {
             f->tsk->sub[k].mat[i][j] = itm[i][j];
        }
   }
/* mark the element 'True' so it doesn't get overwritten
 * and calculate Bit enCoded Matrix. Initialize other struct 
 * members to zero values.
 */
   f->tsk->sub[k].truefalse = True;
   f->tsk->sub[k].bcm = encode_matrix(f->tsk->sub[k].mat );
   f->tsk->sub[k].n_fold = 0;
   f->tsk->sub[k].rotation_angle = 0.0;

   f->rmats_read++;

   f->flags |= HAS_RMAT;
   f->next = read_line;
   return f;
}

static fsm *trans( struct fsm *f )
{
   static const char *fmt = "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf";
   char keyword[12];
   int n_scanned;

   n_scanned = sscanf( f->line, fmt, keyword,
                       &f->tsk->trans_mat[0][0], &f->tsk->trans_mat[0][1], &f->tsk->trans_mat[0][2],
                       &f->tsk->trans_mat[1][0], &f->tsk->trans_mat[1][1], &f->tsk->trans_mat[1][2],
                       &f->tsk->trans_mat[2][0], &f->tsk->trans_mat[2][1], &f->tsk->trans_mat[2][2] );
   if( 10 != n_scanned ) {
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s",
                         f->input_filename, f->line_num, "bad input line for TRANS" );
       f->last_err = -4;
       f->next = NULL;
       return f;
   }
   f->tsk->trans_mat_bcm = encode_matrix( f->tsk->trans_mat );
   f->flags |= HAS_TRANS;
   f->next = read_line;
   return f;
}

static fsm *insfile( struct fsm *f )
{
   f->tsk->shelx_ins_file = get_filename( f->line );
   if( NULL == f->tsk->shelx_ins_file ) {
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "get_filename() returned NULL" );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   f->flags |= HAS_INSFILE;
   f->next = read_line;
   return f;
}

static fsm *outfile( struct fsm *f )
{
   f->tsk->outfile = get_filename( f->line );
   if( NULL == f->tsk->outfile ) {
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "get_filename() returned NULL" );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   f->flags |= HAS_OUTFILE;
   f->next = read_line;
   return f;

   return f;
}

static fsm *exec( struct fsm *f )
{
   if(EXEC_REQUIRES != (EXEC_REQUIRES & f->flags)) {
      gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", f->input_filename, f->line_num,
                        "EXEC requires TRANS and NEWINS to precede it" );
      f->last_err = -6;
      f->next = NULL;
   }   
   f->tsk->shelx_executable = get_filename( f->line );
   if( NULL == f->tsk->shelx_executable ) {
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "get_filename() returned NULL" );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   f->flags |= HAS_EXEC;
   f->next = read_line;
   return f;
}

static fsm *newins( struct fsm *f )
{

   if(NEWINS_REQUIRES != (NEWINS_REQUIRES & f->flags)) {
      gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", f->input_filename, f->line_num,
                        "NEWINS requires INSFILE and TRANS to precede it" );
      f->last_err = -6;
      f->next = NULL;
   }   
   f->tsk->new_base_name = get_filename( f->line );
   if( NULL == f->tsk->new_base_name ) {
       int src_line = __LINE__;
       char *file = __FILE__;
       gen_error_message(f->err_msg, sizeof(f->err_msg), "%s:%d: %s", file, src_line,
                0 != errno ? strerror(errno) : "get_filename() returned NULL" );
       f->last_err = errno;
       f->next = NULL;
       return f;
   }

   f->flags |= HAS_NEWINS;
   f->next = read_line;
   return f;
}

static fsm *end( struct fsm *f )
{
   f->next = close_file;
   f->flags |= HAS_END;
   return f;
}


/* finally, here are the public functions for read_input_file module. */
Queue *read_input_file( char *fname )
{
    Queue *tq = NULL;
    fsm state_machine;
    fsm *state = &state_machine;


    fsm_init( &state_machine );
    state_machine.next = open_file;
    state_machine.input_filename = fname;

    do {
       state = (*(state->next))(&state_machine);
    } while(NULL != state->next );
    
    if( 0 != strcmp( "", state_machine.err_msg ) ) {
        fprintf( stderr, "State machine error code: %d: %s\n", state_machine.last_err,
                  state_machine.err_msg );
    }

    tq = fsm_pass_task_queue( &state_machine );

 
    return tq;
}

void fsm_init( struct fsm *f )
{
   f->next = NULL;
   f->flags = 0;
   f->line_num = 0;
   f->input_filename = NULL;
   f->inp = NULL;
   memset( f->line, 0, sizeof(f->line));
   f->rmats_read = 0;
   f->last_err = 0;
   memset( f->err_msg, 0, sizeof(f->err_msg) );
   f->tsk = NULL;
   f->task_queue = NULL;

   return;
}

Queue *fsm_pass_task_queue( struct fsm *f )
{
    return f->task_queue;
}

