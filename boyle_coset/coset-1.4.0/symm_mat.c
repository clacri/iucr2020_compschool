/* contains implementation for symmetry matrix interface (symm_mat.h) 
 * for the Flack left coset decomposition program which uses
 * uses alogorithms outlined in Acta Cryst. (1987), A43,
 * 564-568, by H. D. Flack.
 *
 * Copyright (C) 2008, 2012, 2013, 2017 Paul D. Boyle
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

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>

#ifdef NEED_C89_COMPATIBILITY
#include "c89_util.h"
#else
#include <stdbool.h>
#endif


#include "symm_mat.h"
#include "matrix.h"
#include "float_util.h"
#include "eigen.h"

struct symm_op_trans_tabl {
       char *name;
       int  id;
       };

/* First some private functions for this file (mostly for diagnostic uses)*/

static void matrix_print( FILE *out, double m[3][3] )
{
#define FMT "%6.2f%6.2f%6.2f\n"
   int k;
   for( k = 0; k < 3; k++ )
        fprintf( out, FMT, m[k][0], m[k][1], m[k][2] );
#undef FMT
}

/* In case we want to output the bit representation of an unsigned int */
static void print_bits( FILE *out, unsigned int s )
{
    int j;

    for( j = CHAR_BIT * sizeof(unsigned int) - 1; j >= 0; j-- ) {
         fprintf( out, "%d", s & (1 << j) ? 1 : 0 );
         if( 0 == j % CHAR_BIT ) fprintf( out, " " );
         fflush( out );
    }
    fputc( '\n', out );

    return;
}

static void print_symm_op_info( FILE *out, struct symm_op *s, int idx, unsigned int flags )
{
   if( flags & TRUTH_VALUE ) {
       fprintf( out, "Symm operator (%d): %s\n", idx+1, s->truefalse == True ? "True" : "False" );
   }
   else {
       fprintf( out, "Symm operator (%d):\n", idx+1 );
   }

   if( flags & HEX_PATTERN ) {
       fprintf( out, "Hex Value of Encoded Matrix: %#x\n", s->bcm );
   }
   if( flags & BIT_PATTERN ) {
       print_bits( out, s->bcm );
   }
   fprintf( out, "Stored Matrix:\n" );
   matrix_print( out, s->mat );
   return;
}

/*
 * The symmetry matricies are encoded as a single 'unsigned int' to facilitate
 * matrix comparisons.  Rather than checking each matrix element individually
 * (a total of 9 double point comparisons), only a single integer value needs
 * to be compared to determine whether a pair of matricies are identical.
 *
 * In crystallographic point group symmetry matricies, the only values which
 * occur are -1, 0, and 1.  These are encoded as follows:
 *		Decimal		Hexadecimal	Binary (Bits)
 *		-1		0x2		10
 *		 0		0x0		00
 *		 1		0x1		01
 * 
 * The elements of the matrix are encoded using 2 bits per element
 * at the following bit positions in an unsigned int.
 *
 * Matrix Element	Starting Bit Offset
 * [0][0]		 0
 * [0][1]		 2
 * [0][2]		 4
 * [1][0]		 6
 * [1][1]		 8
 * [1][2]		10
 * [2][0]		12
 * [2][1]		14
 * [2][2]		16
 *
 * The offsets are tabulated in the offset table matrix (see offset_table).
 */

static unsigned int encode_value( double val, int row, int col, unsigned int codex )
{
     int offset_table[3][3] = {{  0,  2,  4 },
                               {  6,  8, 10 },
                               { 12, 14, 16 }};

     unsigned int cv = 0;
     if( is_equal( val, 1.0 ) ) {
         cv = 0x1;
     }
     else if(is_zero( val ) ) {
         cv = 0x0;
     }
     else if( is_equal( val, -1.0 ) ) {
         cv = 0x2;
     }
     else { /* we have a value which is not 1, 0, or -1, return an error value. */
        return BCM_ERROR;
     }

     codex |= ( cv << offset_table[row][col]);
     return codex;
}

static double decode_value(unsigned int cs, int row, int col )
{
     int offset_table[3][3] = {{  0,  2,  4 },
                               {  6,  8, 10 },
                               { 12, 14, 16 }};

      double ret_val;
      unsigned int extracted = 0;

      extracted = 0x3 & ( cs >> offset_table[row][col] );

      switch( extracted ) {
          case 0x0:
             ret_val = 0.0;
             break;
          case 0x1:
             ret_val = 1.0;
             break;
          case 0x2:
             ret_val = -1.0;
             break;
          default:
              fprintf( stderr, "Unallowed value extracted: %#x\n", extracted );
      }

      return ret_val;

}

/*** Symmetry matricies for the supergroups ***/

/* numbers in comments refer to the symmetry operation list order given
 * Flack's paper p. 567 for spacegroup P 432 (#207)
 */

static struct symm_op ops_432[N_SYMM_OPS_432] = {
              {True, 
                {{ 1.0,  0.0,  0.0},   /* 1 */
                 { 0.0,  1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0,  0.0,  0.0},   /* 3 */
                 { 0.0,  1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0,  0.0,  0.0},   /* 2 */
                 { 0.0, -1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  0.0,  0.0},   /* 4 */
                 { 0.0, -1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 13 */
                 { 1.0,  0.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 14 */  
                 {-1.0,  0.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0,  0.0,  0.0},   /* 18 */
                 { 0.0,  0.0,  1.0},
                 { 0.0,  1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0,  0.0,  0.0},   /* 19 */
                 { 0.0,  0.0, -1.0},
                 { 0.0, -1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0,  1.0},   /* 22 */
                 { 0.0, -1.0,  0.0},
                 { 1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0, -1.0},   /* 24 */
                 { 0.0, -1.0,  0.0},
                 {-1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 16 */
                 { 1.0,  0.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 15 */
                 {-1.0,  0.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0,  1.0},   /* 5 */
                 { 1.0,  0.0,  0.0},
                 { 0.0,  1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0,  1.0},   /* 6 */
                 {-1.0,  0.0,  0.0},
                 { 0.0, -1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0, -1.0},   /* 7 */
                 {-1.0,  0.0,  0.0},
                 { 0.0,  1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0, -1.0},   /* 8 */
                 { 1.0,  0.0,  0.0},
                 { 0.0, -1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 9 */
                 { 0.0,  0.0,  1.0},
                 { 1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 10 */
                 { 0.0,  0.0,  1.0},
                 {-1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 11 */
                 { 0.0,  0.0, -1.0},
                 {-1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 12 */
                 { 0.0,  0.0, -1.0},
                 { 1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  0.0,  0.0},   /* 17 */
                 { 0.0,  0.0,  1.0},
                 { 0.0, -1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  0.0,  0.0},   /* 20 */
                 { 0.0,  0.0, -1.0},
                 { 0.0,  1.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0,  1.0},   /* 21 */
                 { 0.0,  1.0,  0.0},
                 {-1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  0.0, -1.0},   /* 23 */
                 { 0.0,  1.0,  0.0},
                 { 1.0,  0.0,  0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}} };
              

/* numbers in comments refer to the symmetry operation list order given
 * Flack's paper p. 567 for spacegroup P 622 (#177)
 */

static struct symm_op ops_622[N_SYMM_OPS_622] = {
              {True, 
                {{ 1.0,  0.0,  0.0},   /* 1 */
                 { 0.0,  1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 7 */
                 { 1.0,  0.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  0.0,  0.0},   /* 8 */
                 {-1.0, -1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0, -1.0,  0.0},   /* 9 */
                 { 0.0,  1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True, 
                {{-1.0,  0.0,  0.0},   /* 4 */
                 { 0.0, -1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 10 */
                 {-1.0,  0.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0,  0.0,  0.0},   /* 11 */
                 { 1.0,  1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  1.0,  0.0},   /* 12 */
                 { 0.0, -1.0,  0.0},
                 { 0.0,  0.0, -1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{-1.0, -1.0,  0.0},   /* 2 */
                 { 1.0,  0.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0,  1.0,  0.0},   /* 3 */
                 {-1.0, -1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 1.0,  1.0,  0.0},   /* 5 */
                 {-1.0,  0.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}},
              {True,
                {{ 0.0, -1.0,  0.0},   /* 6 */
                 { 1.0,  1.0,  0.0},
                 { 0.0,  0.0,  1.0}},0,0,0.0,0.0,{0.0,0.0,0.0}} };

/**************************************************************/
               
/* implementations for the public functions */

/* encode_matrix(): creates a Bit enCoded Matrix (BCM) from an input 3x3 matrix. */
unsigned int encode_matrix( double fm[3][3] )
{
    int i, j; 
    unsigned int cmx = 0; /* cmx -> Coded Matrix */
 
    for( i = 0; i < 3; i++ ) {
         for( j = 0; j < 3; j++ ) {
              cmx = encode_value( fm[i][j], i, j, cmx );
              if( BCM_ERROR == cmx ) {   /* we have an error */
                  return cmx;
              }
         }
    }
    return cmx;
}

/* decode_matrix(): decodes a Bit enCoded Matrix (BCM) and writes the 3x3 matrix to
 * the 'fm' parameter.
 */
void decode_matrix( double fm[3][3], unsigned int cmx )
{
    int i, j;

    for( i = 0; i < 3; i++ )
         for( j = 0; j < 3; j++ )
              fm[i][j] = decode_value( cmx, i, j );

   
    return;
}

/* unitize_eigen_vector():  function finds the smallest absolute value of an eigen vector component
 * and scales it to a value of '1'.  Other components are scaled similarly.
 */

void unitize_eigen_vector( struct symm_op *s )
{
    int i;
    double smallest = 0.0;
    double arr[3];

    for( i = 0; i < 3; i++ ) {
         arr[i] = fabs( s->eig_vec[i] );
         if( is_zero( arr[i] ) ) {  /* get rid of the 0.0 components.  Hopefully 99.99 is impossibly large  */
             arr[i] = 99.99;
         }
    }

    smallest = arr[0];
    for( i = 1; i < 3; i++ ) {
         if( arr[i] < smallest )
             smallest = arr[i];
    }
    

    for( i = 0; i < 3; i++ ) { /*scale everything to the smallest magnitude */
         s->eig_vec[i] /= smallest;
    }
    return;
}


/* analyze_symm_op():  identifies the type of symmetry element the symm_op
 * matrix represents, and calculates the rotation angle for n_fold rotations.
 */

void analyze_symm_op( struct symm_op *s )
{
    double tr, det, phi;
    int i_tr, i_det;  /* int values for tr and det */
    int i, j;

/* the eigen related arrays declared below are where we store the results of 
 * eigen value, eigen vector calculations.
 */
    struct real_eigen_data real_eig[3];
    struct complex_eigen_data cmplx_eig[3];


/* the formula for phi coded below comes from an old version of the 
 * International Tables.  Unfortunately, I can't find it the volume and page reference.
 * Fortunately, however, this formula is given in the "Rotation Matrix" article on
 * Wikipedia in Section 2.3.2 (as of Nov. 2017). Here is the URL:
 * https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_angle
 */

    det = determinant( s->mat );
    tr  = trace( s->mat );

    phi = acos( (tr/det - 1.0) / 2.0 );
    s->rotation_angle = (float)RAD2DEG( phi );

/* use the trace and determinant to characterize the symmetry element. See 
 * Giacovazzo, "Fundamentals of Crystallography" 1994, p. 43
 */
    i_tr = round_to_nearest_int( tr );
    i_det = round_to_nearest_int( det );

    if( 1 == i_det ) {
        switch( i_tr ) {
           case -1:  /* 2-fold proper rotation */
              s->n_fold = 2;
              break;
           case  0:  /* 3-fold proper rotation */
              s->n_fold = 3;
              break;
           case  1:  /* 4-fold proper rotation */
              s->n_fold = 4;
              break;
           case  2:  /* 6-fold proper rotation */
              s->n_fold = 6;
              break;
           case  3:  /* 1-fold proper rotation */
              s->n_fold = 1;
              break;
           default:  /* we should never get here */
              fprintf( stderr, "%s:%d: Error: i_det = %d, i_tr = %d\n", __FILE__, __LINE__,
                                i_det, i_tr );    
              fprintf( stderr, "These values do not correspond to any crystallographic point symmetry operation.\n" );
              exit( EXIT_FAILURE );
              break;
        }
    }
    else if( -1 == i_det ) {
        switch( i_tr ) {
           case -3:  /* 1-fold improper rotation (inversion center) */
              s->n_fold = -1;
              break;
           case -2:  /* 6-fold improper rotation */
              s->n_fold = -6;
              break;
           case -1:  /* 4-fold improper rotation */
              s->n_fold = -4;
              break;
           case  0:  /* 3-fold improper rotation */
              s->n_fold = -3;
              break;
           case  1:  /* 2-fold improper rotation (mirror plane) */
              s->n_fold = -2;
              break;
           default:  /* we should never get here */
              fprintf( stderr, "%s:%d: Error: i_det = %d, i_tr = %d\n", __FILE__, __LINE__,
                                i_det, i_tr );    
              fprintf( stderr, "These values do not correspond to any crystallographic point symmetry operation.\n" );
              exit( EXIT_FAILURE );
              break;
        }
    }
    else {  /* we have a problem */
       fprintf( stderr, "%s:%d: Error: symmetry matrix determinant is neither +1 or -1\n",
                __FILE__, __LINE__ );
       exit( EXIT_FAILURE );
    }

/* now we work out the eigenvalues and eigenvectors for the symmetry matrix.
 * First thing to do is to determine whether matrix is symmetric or nonsymmetric and
 * solve using the appropriate eigen solver.
 */

/* Then we pick the correct eigenvalue and eigenvector.  The eigen value must be real and the
 * sign of the eigen value should be the same as the sign of the determinant of the symmetry
 * matrix.
 */

    if( is_symmetric( s ) ) {
        symm_eigen_solve( s->mat, real_eig );
        for( i = 0; i < 3; i++ ) {
             if( match_signs(det, real_eig[i].eig_value) ) {
                 s->eig_val = real_eig[i].eig_value;
                 for( j = 0; j < 3; j++ ) { /* copy the eigen data to symm_op */
                      s->eig_vec[j] = real_eig[i].eig_vector[j];
                 }
             }
        }
    }
    else {
        non_symm_eigen_solve( s->mat, cmplx_eig );
        for( i = 0; i < 3; i++ ) {
             if( is_real_eigen_value(&cmplx_eig[i]) ) {
                 int conv_errs = 0;
                 struct real_eigen_data tmp;
                 conv_errs = convert_complex_eig_to_real_eig( &cmplx_eig[i], &tmp );
                 if( 0 != conv_errs ) {
                     print_eigen_results( stderr, "Bad Eigen results!", cmplx_eig, 3 );
                     fprintf( stderr, "%s: %d:  Eigenvector index (%d): Invalid conversion of complex to real. Exiting.\n",
                                       __FILE__, __LINE__, i );
                     exit( EXIT_FAILURE );
                 }
                 if( match_signs( det, tmp.eig_value ) ) {
                     s->eig_val = tmp.eig_value;
                     for( j = 0; j < 3; j++ ) {
                          s->eig_vec[j] = tmp.eig_vector[j];
                     }
                 }
             }
        }
    }

/* make smallest non-zero eigenvector component(s) equal to 1 */
    unitize_eigen_vector( s );

    return;
}

/* transform_group(): transforms an array of struct symm_op from one
 * basis to another.
 */

void transform_group( struct symm_op *g, double tm[3][3] )
{
    int i;

    for( i = 0; 0 != g[i].bcm; i++ ) {
         double tmp[3][3] = { {0.0} };
         similarity_transform( tmp, g[i].mat, tm );
         copy_matrix(g[i].mat, tmp );
         g[i].bcm = encode_matrix( g[i].mat );
    }
    return;
}

void analyze_symm_group( struct symm_op *g ) {
     int i;

     for( i = 0; 0 != g[i].bcm; i++ ) {
          analyze_symm_op( &g[i] );
     }
   
     return;
}


/* lookup_supergroup():  takes a character string designating the supergroup 
 * and returns an integer representation of that group.
 */

int lookup_supergroup( const char *s )
{
    int supergroup_point_group = -1;
    int i,
        n;
    struct symm_op_trans_tabl stab[] = {
           {"-1", 1 },
           {"2/m", 2 },
           {"mmm", 222 },
           {"4/mmm", 422 },
           {"-3m", 32 },
           {"6/mmm", 622 },
           {"m-3m", 432 }
           };

   n = sizeof( stab ) / sizeof( struct symm_op_trans_tabl );
   for( i = 0; i < n; i++ ) {
        if( 0 == strcmp( stab[i].name, s ) ) {
            supergroup_point_group = stab[i].id;
            break;
        }
   }

   return supergroup_point_group;
}


struct symm_op *select_symm_ops( const int pt_group, int *ierr )
{
    struct symm_op *ret = NULL,
                   *sp;
    int n_elem = 0,  /* number of symmetry elements */
        k_elem,      /* k_elem length of array including sentinel value */
        i;           /* counting index */
    int *ip;          

/* The following declarations are arrays of indicies of the symmetry
 * operators for the holoaxial point groups to be used for the G (super)
 * group symmetry.
 */

/*first do the ones for cubic, tetragonal, orthorhombic, monoclinic,
 * and triclinic.
 */
    int m432[N_SYMM_OPS_432] = {0,1,2,3,4,5,6,7,8,9,10,
                                11,12,13,14,15,16,17,18,19,20,
                                21,22,23};
    int m422[N_SYMM_OPS_422] = {0,1,2,3,4,5,10,11};
    int m222[N_SYMM_OPS_222] = {0,1,2,3};
    int m2[N_SYMM_OPS_2] = {0,1};
    int m1[N_SYMM_OPS_1] = {0};

/* here are the ones for hexangonal, trigonal, and rhombohedral groups */
    int m622[N_SYMM_OPS_622] = {0,1,2,3,4,5,6,7,8,9,10,11};
    int m32[N_SYMM_OPS_32] = {0,1,2,3,8,9};

/* here is the inversion center matrix used for expanding the 
 * supergroup symmetry operators.
 */

    double invcen[3][3] = {
          {-1.0,  0.0,  0.0 },
          { 0.0, -1.0,  0.0 },
          { 0.0,  0.0, -1.0 }
                         };


    switch( pt_group ) {
      case 432:
           n_elem = sizeof( m432 ) / sizeof( int );
           ip = m432;
           sp = ops_432;
           *ierr = 0;
           break;
      case 622:
           n_elem = sizeof( m622 ) / sizeof( int );
           ip = m622;
           sp = ops_622;
           *ierr = 0;
           break;
      case 32:
           n_elem = sizeof( m32 ) / sizeof( int );
           ip = m32;
           sp = ops_622;
           *ierr = 0;
           break;
      case 422:
           n_elem = sizeof( m422 ) / sizeof( int );
           ip = m422;
           sp = ops_432;
           *ierr = 0;
           break;
      case 222:
           n_elem = sizeof( m222 ) / sizeof( int );
           ip = m222;
           sp = ops_432;
           *ierr = 0;
           break;
      case 2:
           n_elem = sizeof( m2 ) / sizeof( int );
           ip = m2;
           sp = ops_432;
           *ierr = 0;
         break;
      case 1:
           n_elem = sizeof( m1 ) / sizeof( int );
           ip = m1;
           sp = ops_432;
           *ierr = 0;
         break;
      default:
         *ierr = -1;  /* not an allowed point group, caller deals with error */
          return ret;
         break;
    }

/* now allocate memory for symmetry ops and their inversion related
 * pairs.  Caller must free() this allocation!!!
 */
    k_elem = 1 + 2 * n_elem;  /* +1 to store sentinel value */
    ret = malloc( k_elem * sizeof( *ret ) );
    if( NULL == ret ) {
        *ierr = -2;  /* memory allocation error */
        return NULL;
    }
    ret[k_elem - 1].bcm = 0;  /* sentinel value */
    ret[k_elem - 1].truefalse = False;

/* now copy symmetry matricies to the allocated array */
    for( i = 0; i < n_elem; i++ ) {
         ret[i] = sp[ ip[i] ];
         ret[i].truefalse = True;
         ret[i].bcm = encode_matrix( ret[i].mat );
         matrix_multiply3x3( ret[n_elem + i].mat, ret[i].mat, invcen  );
         ret[n_elem+i].bcm = encode_matrix( ret[n_elem+i].mat );
         ret[n_elem+i].truefalse = True;
    }
   
    return ret;

}

/* set_truth_value(): sets truefalse value to 'value' for 'total',
 * selected indicies. Returns 0 if everything OK, -1 otherwise.
 *
 * if 'total' is <= 0 then all the truefalse values are set to 'value'
 */
int set_truth_value( struct symm_op *s, int value, int total, ... )
{
   int i, idx, ret_code;
   int c; 
   va_list ap;

/* first do a sanity check on 'value' */
   switch( value ) {
       case True:
          ret_code = 0;
          break;
       case False:
          ret_code = 0;
          break;
       default:
            return -1;
   }

/* count number of symmetry ops for sanity check for array bounds */
   c = count_ops( s );
   c = c - 1;

/* for total > 0, do the 'pick and choose' mode */
   if( total >  0 ) {
       while( i < total ) {
           idx = va_arg( ap, int );
           if( (idx < 0) || (idx > c) )  {
               return -1;
           }
           s[idx].truefalse = value;
           i++;
       }

       va_end( ap );
   }
   else {  /* else -- set them all the same value */
      for( i = 0; 0 != s[i].bcm; i++ ) {
           s[i].truefalse = value;
      }
   }
   return ret_code;
   
}


/* count_ops() -- returns number of symmetry operations
 * in an array of struct symm_op pointed to by 's'.
 */
int count_ops( struct symm_op *s )
{
    int count;

    if( NULL == s )
        return 0;

    for( count = 0; 0 != s[count].bcm; count++ );

    return count;
}

/* duplicate_ops(): like dupstr, except for an array of struct symm_op.
 * caller must free() the allocated memory.
 */

struct symm_op *duplicate_ops( struct symm_op *s )
{
    struct symm_op *ret = NULL;
    int n_ops;
    size_t sz;

    n_ops = count_ops( s ) + 1;
    sz = n_ops * sizeof(*ret);
    errno = 0;
    ret = malloc( sz );
    if( NULL != ret )
        memcpy( ret, s, sz );

    return ret;
}


/* is_centric() return 1 if one of the members of the 
 * symm_op array is an inversion operator.  Otherwise,
 * the function returns 0.
 */

int is_centric( struct symm_op *s )
{
    int i,
        truth_val = 0;

    for( i = 0; 0 != s[i].bcm; i++ ) {
         if( INVERSION_BCM == s[i].bcm ) {
             truth_val = 1; 
             break;
         }
    }
    
    return truth_val;
}

/* is_symmetric() return 1 if a symmetry matrix is symmetric (i.e. the 
 * transpose of the input matrix is equal to the input matrix). 
 * The function returns 0 if the matrix is not symmetric and -1 if there
 * has been some error.
 */


int is_symmetric ( struct symm_op *s )
{

    struct symm_op tr = {0,
                         { { 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0 }
                         }, 0, 0, 0.0,0.0,{0.0,0.0,0.0}
                        }; 

    memcpy( &tr, s, sizeof(tr) );
    transpose_matrix( tr.mat );
    tr.bcm = encode_matrix( tr.mat );
    return tr.bcm == s->bcm ? 1 : 0 ;
}


/* print_symm_ops() -- publically callable version of 
 * print_symm_opt_info()
 */
void print_symm_ops( FILE *out, struct symm_op *s, unsigned int flags )
{
    int i;
    for( i = 0; 0 != s[i].bcm; i++ ) {
         print_symm_op_info( out, &s[i], i, flags );
    }
    return; 
}

/* print_2_symm_ops(): an unfortunately complicated function parameter list
 * for printing out two arrays of struct symm_op's in tandem. The array 
 * dimensions must be equal.
 */
void print_2_symm_ops( FILE *out, const char *header1, const char *header2, 
                       struct symm_op *s1, struct symm_op *s2 )
{
#define PRINT2_FORMAT "%6.2f%6.2f%6.2f\t\t\t%6.2f%6.2f%6.2f\n"
    int i, j,
        sc1,
        sc2;

    sc1 = count_ops( s1 );
    sc2 = count_ops( s2 );
    if( sc1 != sc2 ) {
        fprintf( out, "Arrays not equal length: s1 has %d elements, s2 has %d elements\n",
                       sc1, sc2 );
        return;
    }

    fprintf( out, "%s\t%s\n", header1, header2 );
    for( i = 0; 0 != s2[i].bcm; i++ ) {
         if( True == s2[i].truefalse ) {
             if( INVERSION_BCM == s2[i].bcm ) { /* deal with special case for inversion twinning */
                  fprintf( out, "Twin domains related by inversion\n" );
             }
             else if( 0 != s2[i].n_fold ) {  /* deal with rotation twin laws */
                  fprintf( out, "** %d-fold (%s) rotation axis about the [%4.2f %5.2f %5.2f] direction. **\n",
                           s2[i].n_fold, s2[i].n_fold > 0 ? "proper" : "improper",
                           s2[i].eig_vec[0], s2[i].eig_vec[1], s2[i].eig_vec[2] );
             }
             for( j = 0; j < 3; j++ ) {
                  fprintf( out, PRINT2_FORMAT, s1[i].mat[j][0],s1[i].mat[j][1],s1[i].mat[j][2],
                                               s2[i].mat[j][0],s2[i].mat[j][1],s2[i].mat[j][2] );
             }
             fputs( "\n\n", out );
         }
    }
    fputc( '\n', out );
    fflush( out );
    return;
#undef PRINT2_FORMAT
}


/* this function checks the encoding/decoding routines for
 * the symmetry matricies.
 */
int symm_op_diagnostic( FILE *output )
{
    int i, k,
        n_elem = 0;
    int err = 0;
    double fmat[3][3] = { {0.0} };

    int groups[7] = {1, 2, 222, 422,  32, 622, 432 };
                      
    struct symm_op *s = NULL;

/* if no file handle is given, then write to stdout */
    if( NULL == output )
        output = stdout;

    for( i = 0; i < 7; i++ ) {
         switch( groups[i] ) {
             case 432:
             n_elem = 2 * N_SYMM_OPS_432;
             break;
             case 622:
             n_elem = 2 * N_SYMM_OPS_622;
             break;
             case 32:
             n_elem = 2 * N_SYMM_OPS_32;
             break;
             case 422:
             n_elem = 2 * N_SYMM_OPS_422;
             break;
             case 222:
             n_elem = 2 * N_SYMM_OPS_222;
             break;
             case 2:
             n_elem = 2 * N_SYMM_OPS_2;
             break;
             case 1:
             n_elem = 2 * N_SYMM_OPS_1;
             break;
         }
         s = select_symm_ops( groups[i], &err );
         fprintf( output, "\n\nEncodings and Matricies for Group: %d\n",
                          groups[i] );

         for( k = 0; k < n_elem; k++ ) {
             double tmp[3][3];
             print_symm_op_info( output, &s[k], k, TRUTH_VALUE|BIT_PATTERN|HEX_PATTERN );
             decode_matrix( fmat, s[k].bcm);
             fprintf( output, "Decoded Matrix:\n" );
             matrix_print( output, fmat );
             fprintf(output, "Difference Matrix:\n" );
             matrix_subtract3x3( tmp, s[k].mat, fmat );
             matrix_print( output, tmp );
             fputs("", output);
         }
         free( s );

    }     
    return 0;
}
