/* implementation of functions used for reading, editing, and writing
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
#ifdef NEED_C89_COMPATIBILITY
#undef _ISOC99_SOURCE
#else
#define _ISOC99_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "symm_mat.h"
#include "coset.h"
#include "dupstr.h"
#include "shelx.h"
#include "sll.h"
#include "dynamic_sll.h"

/* for compilers which don't have C99 functions available */
#ifdef NEED_C89_COMPATIBILITY
#include "c89_util.h"
#endif

#ifdef USE_SYSTEM_FUNCTION
#define USE_NONSTANDARD_FOPEN  /* for nonconforming C implementations of fopen() */
#endif



/* first are some static functions for private use in this file scope */
static char *create_basf_instruction( int n )
{
    static char *basf = "BASF ";
    char  *p, *ret;
    float starting_value = 0.0;
    char buf[SHELX_LINE_LEN] = {0};
    int n_comp = 0, n_abs = 0;
    int i, k_val = 0;
    int i_ovr;  /* tests for snprintf() overrun */

    n_abs = abs( n );

    if( (1 == n_abs) || (2 == n_abs) ) { /* for identity, inversion, and 2 fold twins */
         k_val = 1;
         n_comp = 2;
    }
    else { /* for trilled and higher order multiple domains */
         k_val = n_abs - 1;
         n_comp = n_abs;
    }
    starting_value = 1.0 / (float)n_comp;

    p = buf;
    sprintf( p, "%s", basf );

/* loop over k_val filling buf with snprintf() */
    for( i = 0; i < k_val; i++ ) {
         p = buf + strlen( buf );
         if( strlen(buf) + BASF_BUF_LEN <  SHELX_LINE_LEN ) {
             i_ovr = snprintf( p, BASF_BUF_LEN, BASF_SCALE_FORMAT, starting_value );
         }
         else {
             fprintf( stderr, "%s:%d -- Potential buffer overrun on BASF instruction.", 
                      __FILE__, __LINE__ );
             fprintf( stderr, "snprintf() reports writing %d bytes\n", i_ovr ); 
         }
    }
    buf[strlen(buf)] = '\n';

    errno = 0;
    ret = dupstr( buf );
    return ret;  /*caller must free() the returned pointer. */

}
static char *format_shelx_twin_instruction( double twin_law[3][3], int nc )
{
    int n_comp, ret;
    size_t buf_size = TWIN_INS_BUF_LEN;
    char tmp[SHELX_LINE_LEN] = {0}, *basf, *twin_ins;

    switch( nc ) {
         case -1:
            /* fall through */
         case  1:
            /* fall through */
         case  2:
            n_comp = 2;
            break;
         default:
            n_comp = abs( nc );
            break;
    }

    basf = create_basf_instruction( nc );
    if( NULL == basf ) {
        fprintf( stderr, "%s:%d -- %s\n", __FILE__,__LINE__, errno != 0 ? strerror(errno) : "NULL pointer returned from create_basf_instruction()" );
    }

   
    errno = 0;
    twin_ins = malloc( buf_size );
    if( NULL == twin_ins )
        return NULL;
    ret = snprintf( tmp, sizeof(tmp), TWIN_INS_FORMAT, twin_law[0][0], twin_law[0][1], twin_law[0][2],
                                                       twin_law[1][0], twin_law[1][1], twin_law[1][2],
                                                       twin_law[2][0], twin_law[2][1], twin_law[2][2], n_comp );
    if( (size_t)ret >=  sizeof(tmp) )
        fprintf( stderr, "%s:%d -- TWIN instruction truncated.\n", __FILE__, __LINE__ );

    ret = snprintf( twin_ins, buf_size, "%s%s", basf, tmp );
    if( (size_t)ret >=  buf_size )
        fprintf( stderr, "%s:%d -- BASF/TWIN instruction truncated.\n", __FILE__, __LINE__ );

    free( basf );
    return twin_ins;  /* caller should free() this */

}

/* write_shelx_ins_file(): write a SHELX .ins file with a filename of 'name' and with the
 * contents stored in the 'ins' linked list.
 */
static void write_shelx_ins_file( char *name, SLinkedList *ins )
{
    FILE *ins_fp;
    SLinkedListElem *el;
#ifdef USE_NONSTANDARD_FOPEN
    const char *mode = "wt";
#else  /* we use only modes defined by the ANSI C Standard */
    const char *mode = "w";
#endif
 

    ins_fp = fopen( name, mode );
    if( NULL == ins_fp ) {
	fprintf( stderr, "Error: %s: %s\n", name, errno != 0 ? strerror(errno): "couldn't open file" );
    }

    fprintf( stdout, "Writing new SHELX .ins file %s ...\n", name );
    for( el = sll_list_head(ins); NULL != el; el = sll_list_next(el) ) {
	 fputs( sll_list_data(el), ins_fp );
    }

    fclose( ins_fp );
    return;
}

static SLinkedList *shelx_insert_twin_ins( char *twin, SLinkedList *s )
{

   const char *match_this = "FVAR";  /* Doesn't have to be FVAR, can be any SHELX instruction
				      * which has to be in .ins file.  Could be "UNIT", for example.
				      */
   int match_flag = 0;
   SLinkedList *edited_list = NULL;
   SLinkedListElem *el;
   size_t match_len;
   
   match_len = strlen( match_this );
   edited_list = alloc_list_init( edited_list, free );
   if( NULL == edited_list )
       return NULL;


   for( el = sll_list_head(s); NULL != el; el = sll_list_next(el) ) {
	size_t len;
	char *tmp; 
	len = strlen(sll_list_data(el) ) + 1;
        errno = 0;
	tmp = malloc( len );   /* this probably isn't the fastest way to do this, but
                                * doing a 'deep copy' probably will make cleaning up
                                * memory allocations at the end of the program less messy.
                                */
	if( NULL == tmp ) {  /* deal with error */
            fprintf( stderr, "%s:%d: %s\n", __FILE__, __LINE__, errno != 0 ? strerror(errno) :
                             "malloc() returned NULL" );
	}
	strncpy(tmp, sll_list_data(el), len );
	sll_insert_next( edited_list, sll_list_tail(edited_list), tmp );
	
	if( 0 == strncmp(match_this, sll_list_data(el), match_len) ) {
	    sll_insert_next( edited_list, sll_list_tail(edited_list), twin );
	    match_flag = 1;
	}
   } 
   if( 0 == match_flag )
       fprintf( stderr, "SHELX instruction \"%s\" not found in this list.\n", match_this);

   return edited_list;

}



SLinkedList *twin_ins_list( struct symm_op *s )
{
    SLinkedList *list = NULL;
    char *tp;    

    int i;

    list = alloc_list_init( list, free );
    if( NULL == list )
        return NULL;

    for( i = 0; 0 != s[i].bcm; i++ ) {
         if( True == s[i].truefalse ) {
             tp = format_shelx_twin_instruction( s[i].mat, s[i].n_fold );
             sll_insert_next( list, sll_list_tail(list), tp );
         }
    }

    return list;
}

SLinkedList *read_shelx_ins_file( char *ins_file_name )
{
    FILE *ins;
    SLinkedList *ins_list = NULL;  /* lines of the file in a single linked list */

    char line[SHELX_LINE_LEN];
    char *scpy;
#ifdef USE_NONSTANDARD_FOPEN
    const char *mode = "rt";
#else  /* use only ANSI C Standard mode */
    const char *mode = "r";
#endif

    errno = 0;
    ins = fopen( ins_file_name, mode );
    if( NULL == ins ) {
        fprintf( stderr, "Error: %s: %s\n", ins_file_name, errno != 0 ? strerror(errno): "couldn't open file" );
        return NULL;
    }
    errno = 0; 
    ins_list = alloc_list_init( ins_list, free );
    if( NULL == ins_list ) {
        fclose( ins );
        return NULL;
    }

/* copy and duplicate the strings (dupstr()) and insert the line into this list */
    while( NULL != fgets( line, (int)sizeof(line), ins )) {  /* don't bother taking off '\n' characters */
           scpy = dupstr( line );
           if( NULL == scpy ) { /* error */
               fclose( ins );
               sll_destroy( ins_list );
               free( ins_list );
               return NULL;
           }
           sll_insert_next( ins_list, sll_list_tail(ins_list), scpy );
    }

    fclose( ins );
    return ins_list;
}

char *get_basename( char *filename, int delim_char )
{
    char *p0, *p1, *p2;
    size_t len;

    p0 = filename;
    p1 = strrchr( p0, delim_char );
    len = p1 - p0;

    p2 = dupnstr( p0, len );
    return p2;
}

SLinkedList *write_new_ins_files( char *base_name, SLinkedList *twin_laws, SLinkedList *ins )
{
#define OUTPUT_FILENAME_FORMAT "%s_%02d.ins"
    int i = 0,
        ret;
    size_t len;
    SLinkedList *twin_list_ins;
    char *ins_file_name = NULL,  /* the SHELX .ins file name */
         *basf_twin;             /* pointer char buffer which contains BASF and 
                                  * TWIN instructions.
                                  */
    SLinkedList *new_ins_file_names = NULL;

    new_ins_file_names = alloc_list_init( new_ins_file_names, free );
    errno = 0;
    if( NULL == new_ins_file_names ) {
        return NULL;
    }

    while( -1 != sll_remove_next( twin_laws, NULL, (void **)&basf_twin ) ) {
         i++;
         len = strlen(base_name) + 8;
         errno = 0;
         ins_file_name = malloc(len);
         if( NULL == ins_file_name ) {
             fprintf( stderr, "%s:%d: malloc(): %s\n", __FILE__, __LINE__, 
                      errno != 0 ? strerror(errno) : "could not allocate memory!\n" );
             dealloc_list( new_ins_file_names );
             return NULL;
         }
         twin_list_ins = shelx_insert_twin_ins( basf_twin, ins );
         ret = snprintf( ins_file_name, len, OUTPUT_FILENAME_FORMAT, base_name, i );
         if( (size_t)ret >= len ) {
             fprintf( stderr, "Truncated filename: %s\n", ins_file_name );
         }
         write_shelx_ins_file( ins_file_name, twin_list_ins );
         ret = sll_insert_next( new_ins_file_names, sll_list_tail(new_ins_file_names),ins_file_name );
         if( -1 == ret ) {
             fprintf( stderr, "%s:%d -- Couldn't insert \"%s\" into \"new_ins_file_names\"\n", 
                      __FILE__, __LINE__, ins_file_name ); 
         }     
         dealloc_list( twin_list_ins );
    }

    return new_ins_file_names;
#undef OUTPUT_FILENAME_FORMAT
}
