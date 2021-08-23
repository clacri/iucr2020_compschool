/* implementation of matrix operations for 3x3 matricies 
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


#define _ISOC99_SOURCE
#include <stdio.h>
#include <string.h>

#include "matrix.h"
#include "float_util.h"
#ifdef USER_DIALOG
#include "user_dialog.h"
#endif

/***********************************************************************/
/* determinant() */
/* This is a function which calculates the determinant
 * for a N_DIMxN_DIM matrix.  Nothing fancy here.  Maybe this
 * should be rewritten with pointer arithmetic, but the
 * hell with that for now.
 */

double determinant( double m[N_DIM][N_DIM] )
{
    double determ;

    determ = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) +
             m[1][0] * (m[2][1] * m[0][2] - m[0][1] * m[2][2]) +
             m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);

    return determ;
}

/***********************************************************************/
double trace( double m[N_DIM][N_DIM] )
{
    double tr = 0.0;
    int i;
    
    for( i = 0; i < N_DIM; i++ )
         tr += m[i][i];

    return tr;
}
/***********************************************************************/
/* invert_matrix() */
void invert_matrix( double det, double m[N_DIM][N_DIM], double inv[N_DIM][N_DIM])
{
    inv[0][0] =   (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
    inv[0][1] = - (m[0][1] * m[2][2] - m[0][2] * m[2][1]) / det;
    inv[0][2] =   (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
    inv[1][0] = - (m[1][0] * m[2][2] - m[2][0] * m[1][2]) / det;
    inv[1][1] =   (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    inv[1][2] = - (m[0][0] * m[1][2] - m[0][2] * m[1][0]) / det;
    inv[2][0] =   (m[1][0] * m[2][1] - m[2][0] * m[1][1]) / det;
    inv[2][1] = - (m[0][0] * m[2][1] - m[0][1] * m[2][0]) / det;
    inv[2][2] =   (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;

    return;
}
/***********************************************************************/
/* negate_matrix() switches the signs on any terms in matrix which != 0
 * terms are copied to 'negated' function parameter.
 */
void negate_matrix( double negated[N_DIM][N_DIM], double m[N_DIM][N_DIM] )
{
   int i, j;

   for( i = 0; i < N_DIM; i++ ) {
        for( j = 0; j < N_DIM; j++ ) {
             if( !is_zero( m[i][j] ) ) {
                 negated[i][j] = -m[i][j];
             }
             else {
                 negated[i][j] = m[i][j];
             }
        }
   }
   return;
}

/***********************************************************************/
void transpose_matrix( double m[N_DIM][N_DIM] )
{
    double tmp;
    int i, j;

    for( i = 0; i < N_DIM; i++ )
       for( j = 0; j < i; j++ ) {
            tmp = m[i][j];
            m[i][j] = m[j][i];
            m[j][i] = tmp;
       }

}

/***********************************************************************/
void matrix_multiply3x3( double prod[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM] )
{
    int i, j, k;
    for( i = 0; i < N_DIM; i++ )
         for( j = 0; j < N_DIM; j++ ) {
              prod[i][j] = 0.0;
              for( k = 0; k < N_DIM; k++ )
                   prod[i][j] += a[i][k] * b[k][j];
         }
}

/***********************************************************************/
void matrix_add3x3( double res[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM] )
{
     int i, j;

     for( i = 0; i < N_DIM; i++ )
          for( j = 0; j < N_DIM; j++ )
             res[i][j] = a[i][j] + b[i][j];

     return;
}
/***********************************************************************/
void matrix_subtract3x3( double res[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM] )
{
     int i, j;

     for( i = 0; i < N_DIM; i++ )
          for( j = 0; j < N_DIM; j++ )
             res[i][j] = a[i][j] - b[i][j];

     return;
}
/***********************************************************************/
void calculate_inverse_transpose( double inv_tr[N_DIM][N_DIM], double matrx[N_DIM][N_DIM] )
{
    double det;

    det = determinant( matrx );
    invert_matrix( det, matrx, inv_tr );
    transpose_matrix( inv_tr );

    return;
}

/***********************************************************************/
void similarity_transform( double transformed[N_DIM][N_DIM], double mat[N_DIM][N_DIM], double trans_mat[N_DIM][N_DIM] )
{
    double det;
    double invrt[N_DIM][N_DIM] = { { 0.0 } };
    double tmp[N_DIM][N_DIM] = { { 0.0 } };

    det = determinant( trans_mat );
    invert_matrix( det, trans_mat, invrt );
    matrix_multiply3x3( tmp, mat, invrt );
    matrix_multiply3x3( transformed, trans_mat, tmp );

    return;
}
/***********************************************************************/
void copy_matrix( double dest[N_DIM][N_DIM], double src[N_DIM][N_DIM] )
{
    int i, j;
    for( i = 0; i < N_DIM; i++ )
         for( j = 0; j < N_DIM; j++ )
              dest[i][j] = src[i][j];

    return;
}
/***********************************************************************/

#ifdef USE_DIALOG
void get_matrix( double m[N_DIM][N_DIM] )
{
#define N_EXPECTED 3
   int i,
       n_scanned;
   int read_status = 0;
   size_t len;


   char line[LINE_LEN],
        *ptr;

   for( i = 0; i < N_DIM; i++ ) {
      int tmp_i = i + 1;
      do {
      read_status = 0;
      fprintf( stdout, "Enter elements %d1 %d2 %d3: ",
                       tmp_i, tmp_i, tmp_i );
      fflush( stdout );
      len = get_line( stdin, line, LINE_LEN );
      if( 0 == len ) {
          fputs( "Error in getting line.", stderr );
          return;
      }
      if( NULL != (ptr = strchr(line, '\n')) ) *ptr = '\0';
         n_scanned = sscanf( line, "%lf%*[\t ,]%lf%*[\t ,]%lf",
                             &m[i][0], &m[i][1], &m[i][2] );
         if( N_EXPECTED != n_scanned ) {
             read_status = 1; /* try again */
             fprintf( stderr, "Expected %d items, you typed %d, try again.\n",
                               N_EXPECTED, n_scanned );
         }
      } while( read_status );
  }

#undef N_EXPECTED
}
#endif

/***********************************************************************/
void print_matrix( FILE *out, double m[N_DIM][N_DIM], const char *format )
{
    int i;
    for( i = 0; i < N_DIM; i++ ) {
         fprintf( out, format, m[i][0],m[i][1],m[i][2] );
    }

    return;
}

