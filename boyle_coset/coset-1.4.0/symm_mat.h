/* contains public interface for symmetry matrix routines 
 * for the Flack left coset decomposition program which uses
 * uses alogorithms outlined in Acta Cryst. (1987), A43,
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

#ifdef  NEED_C89_COMPATIBILITY
#undef _ISOC99_SOURCE 
#else
#define _ISOC99_SOURCE
#endif

#ifndef SYMM_MAT_H
#define SYMM_MAT_H 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "sll.h"

#ifdef  NEED_C89_COMPATIBILITY
#include "c89_util.h"
#else
#include <stdbool.h>
#endif

#define True  1
#define False 0

#define N_SYMM_OPS_432      24
#define N_SYMM_OPS_622      12
#define N_SYMM_OPS_422       8
#define N_SYMM_OPS_32        6  /* point group 32 in hexagonal setting */
#define N_SYMM_OPS_222       4
#define N_SYMM_OPS_2         2
#define N_SYMM_OPS_1         1

#define BCM_ERROR            0x2aaaa  /* gives a matrix with all elements as -1 */
#define INVERSION_BCM        0x20202
#define Centric  1
#define Acentric 0

/* for read_symm_file() and write_symm_file():
 * for both read and write statements, the order of the matrix elements
 * should be:
 * 11 12 13 21 22 23 31 32 33
 */

#define READ_FORMAT "%4s%lf%lf%lf%lf%lf%lf%lf%lf%lf"
#define WRITE_FORMAT "RMAT %7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f"
#define READ_EXPECT 10

#define READ_LINE_LEN 132


/* set bits for the 'flags' parameter for print_symm_ops() */
#define TRUTH_VALUE  (1 << 0)
#define BIT_PATTERN  (1 << 1)
#define HEX_PATTERN  (1 << 2)



/* basic data structure for storing symmetry operations. */
struct symm_op {
       bool truefalse;
       double mat[3][3];
       unsigned int bcm;  /* bcm => Binary (en)coded matrix 
                           * see the implementation file for
                           * details.
                           */
       int   n_fold;       /* identifies the type of rotation.  Used for information as well as
                            * the number of BASF parameters for the SHELX refinement.
                            */
       float rotation_angle; /* (in degrees) angle of rotation for symm element */
       double eig_val;
       double eig_vec[3];
       };

/* public prototypes */

/* caller must free() pointer returned by 
 * select_symm_ops()
 */
struct symm_op *select_symm_ops( const int pt_group, int *ierr );
int lookup_supergroup( const char *s );
void transform_group( struct symm_op *g, double tm[3][3] );
unsigned int encode_matrix( double fm[3][3] );
void decode_matrix( double fm[3][3], unsigned int cmx );
void unitize_eigen_vector( struct symm_op *s );
void analyze_symm_op( struct symm_op *s );
void analyze_symm_group( struct symm_op *g );
int is_centric( struct symm_op *s );
int is_symmetric( struct symm_op *s );
int count_ops( struct symm_op *s );
struct symm_op *duplicate_ops( struct symm_op *s );
int set_truth_value( struct symm_op *s, int value, int total, ... );
struct symm_op *read_symm_file( const char *fname );
int write_symm_file(const char *fname, struct symm_op *s);
void print_2_symm_ops( FILE *out, const char *header1, const char *header2,
                       struct symm_op *s1, struct symm_op *s2 );

/* for debugging only: */
void print_symm_ops( FILE *, struct symm_op *, unsigned int  );
int symm_op_diagnostic( FILE * );

#endif
