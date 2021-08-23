/* contains implementation for COSET's eigen interface (eigen.h)
 * Copyright (C) 2008, 2012, 2013, 2017, 2021 Paul D. Boyle
 *
 * Written by:
 *      Paul D. Boyle
 *      Department of Chemistry
 *      University of Western Ontario
 *      London, Ontario CANADA N6A 5B7
 *
 *      May 2021
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
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <complex.h>

#include <lapacke.h>

#include "matrix.h"
#include "eigen.h"
#include "float_util.h"

#define STORE_FULL_MATRIX             0
#define STORE_UPPER_TRIANGULAR_MATRIX 1
#define STORE_LOWER_TRIANGULAR_MATRIX 2

#define N_DIM 3 
#define LDA N_DIM
#define LDVR N_DIM

/* start with some routines which will bring the matrix forms used by COSET into a
 * form usable by LAPACKE.
 */

/* fill_full_matrix(): duplicates a matrix in its entirety and returns a dynamically
 * allocated array containing the matrix's contents.
 *
 * Called by convert_3x3_matrix()
 */
static double *fill_full_matrix( double mt[][N_DIM] )
{
    int i, j;
    int n_row = 3;
    int n_elements = 0;
    double *mp = NULL;

    errno = 0;
    mp = malloc( (n_row * N_DIM) * sizeof(*mp) );
    if( NULL == mp )
        return NULL;

    n_elements = n_row * N_DIM;

    /* initialize array to all 0.0 */
    for( i = 0; i < n_elements; i++ )
         mp[i] = 0.0;

    for( i = 0; i < n_row; i++ ) {
         for( j = 0; j < n_row; j++ ) {
              mp[i * n_row + j] = mt[i][j];
         }
    }

    return mp;
}

/* fill_upper_matrix(): duplicates the upper triangular part of a symmetric matrix
 * and returns a dynamically allocated array containing the matrix's contents.
 *
 * Should only be used for symmetric matrices
 *
 * Called by convert_3x3_matrix()
 */
static double *fill_upper_matrix( double mt[][N_DIM] )
{
    int i, j;
    int n_row = 3;
    int n_elements = 0;
    double *mp = NULL;

    errno = 0;
    mp = malloc( (n_row * N_DIM) * sizeof(*mp) );
    if( NULL == mp )
        return NULL;

    n_elements = n_row * N_DIM;

    /* initialize array to all 0.0 */
    for( i = 0; i < n_elements; i++ )
         mp[i] = 0.0;

    for( i = 0; i < n_row; i++ ) {
         for( j = i; j < n_row; j++ ) {
              mp[i * n_row + j] = mt[i][j];
         }
    }

    return mp;
}

/* fill_lower_matrix(): duplicates the lower triangular part of a symmetric matrix
 * and returns a dynamically allocated array containing the matrix's contents.
 *
 * Should only be used for symmetric matrices
 *
 * Called by convert_3x3_matrix()
 */
static double *fill_lower_matrix( double mt[][N_DIM] )
{
    int i, j;
    int n_row = 3;
    int n_elements = 0;
    double *mp = NULL;

    errno = 0;
    mp = malloc( (n_row * N_DIM) * sizeof(*mp) );
    if( NULL == mp )
        return NULL;

    n_elements = n_row * N_DIM;

    /* initialize array to all 0.0 */
    for( i = 0; i < n_elements; i++ )
         mp[i] = 0.0;

    for( i = 0; i < n_row; i++ ) {
         for( j = 0; j <= i ; j++ ) {
              mp[i * n_row + j] = mt[i][j];
         }
    }

    return mp;
}
/* convert_3x3_matrix():  takes a statically allocated 3x3 matrix and creates a
 * dynamically allocated array suitable for input into LAPACKE functions.  Allowed
 * values of 'flag' are specified above.
 *
 *  ** Caller should free() returned pointer **
 */
static double *convert_3x3_matrix( double m[N_DIM][N_DIM], int flag )
{
    double *(*f)( double [N_DIM][N_DIM] );

    switch( flag ) {
       case STORE_FULL_MATRIX:
          f = fill_full_matrix;
          break;
       case STORE_UPPER_TRIANGULAR_MATRIX:
          f = fill_upper_matrix;
          break;
       case STORE_LOWER_TRIANGULAR_MATRIX:
          f = fill_lower_matrix;
          break;
       default:
          fprintf( stderr, "Unrecognized switch() option: %s: %d\n", __FILE__, __LINE__ );
          exit( EXIT_FAILURE );
          break;
    }

    return f(m);
}
/* reform_3x3_matrix(): takes dynamically allocated array and fills a statically
 * allocated 3x3 matrix with its contents.
 */

static void reform_3x3_matrix( double *arr, double m[N_DIM][N_DIM] )
{
    int i, j;
    for( i = 0; i < N_DIM; i++ ) {
         for( j = 0; j < N_DIM; j++ ) {
              m[i][j] = arr[ i * N_DIM + j];
         }
    }
    return;
}

/* marshal_eigen_results(): used to extract results from LAPACKE calculation and bring into into an
 * array of 'struct complex_eigen_data' based on the 'print_eigenvectors function found at this link:
 * https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/lapacke_dgeev_row.c.htm
 * 
 * and make appropriate changes
 */

static void marshal_eigen_results( struct complex_eigen_data *eig, int n_dim, double *val_re, double *val_im, double *vecs, int lenv )
{   
    int i, j;
    double complex vc[N_DIM][N_DIM] = {{0.0}};

    
    /* first collect the eigenvalues */
    for( i = 0; i < n_dim; i++ ) {
         eig[i].eig_value = val_re[i] + val_im[i] * I;
    }
    /* second, get the eigenvector components into complex numbers */
    for( i = 0; i < n_dim; i++ ) {
         j = 0; 
         while( j < n_dim ) {
             if(cimag(eig[j].eig_value)  == (0.0 * I) ) {
                 vc[i][j] = vecs[i*lenv + j] + 0.00 * I;
                 j++;
             }
             else {
                 vc[i][j] = vecs[i*lenv + j] + vecs[i*lenv+(j+1)] * I;
                 vc[i][j+1] = conj( vc[i][j] );
                 j += 2;
             }
         }
    }
    
    /* finally, copy the eigenvector components to the complex_eigen_data struct 'eig' */
    for( i = 0; i < n_dim; i++ ) {
         for( j = 0; j < n_dim; j++ ) {
              eig[i].eig_vector[j] = vc[j][i];
         }
    }

    return;
}


/* below are the implementations for the public functions. */

/* is_real_eigen_value(): returns 1 if eigenvalue has no imaginary component or 0 if it does */
int is_real_eigen_value( struct complex_eigen_data *d )
{
    return is_zero( (double)cimag(d->eig_value) );
}

/* convert_complex_eig_to_real_eig():  this function is a convenience function to bring complex eigen values
 * which happen to have all zeros for the imaginary components into a struct real_eigen_data structure.
 */
int convert_complex_eig_to_real_eig( struct complex_eigen_data *in, struct real_eigen_data *out )
{
    int i;
    int err_code = 0;

    if( is_real_eigen_value( in ) ) {
        out->eig_value = (double)creal(in->eig_value);
    }
    else {
         err_code  = -1;
    }
    for( i = 0; i < 3; i++ ) {
         if( is_zero( cimag( in->eig_vector[i] ) ) ) {
             out->eig_vector[i] = (double)creal(in->eig_vector[i]);
         }
         else {
             err_code--;  /* decrement the err_code when non-zero imaginary components are encountered */
         }
    }

    return err_code; /* valid conversions return 0. Anything else is an error */
}

/* a macro for printing complex numbers see:
 * stackoverflow.com/questions/4099433/c-complex-number-and-printf
 *  used in the print_eigen_results() function below.
 */
#define printfc(f,c) (fprintf(f, "%6.2f%+6.2fi ",creal(c),cimag(c)))

/* print_eigen_results() included for debugging */

void print_eigen_results( FILE *out, char *desc, struct complex_eigen_data *eig, int n_dim )
{
    int i, j;
    printf( "%s\n", desc );

    fprintf( out, "\tEigenvalue\t\tEigenvector Components\n" );
    for( i = 0; i < n_dim; i++ ) {
         fprintf( out, "[%d]  %6.2f%+6.2fi: ", i, creal(eig[i].eig_value),cimag(eig[i].eig_value) );
         for( j = 0; j < 3; j++ ) {
              printfc(out, eig[i].eig_vector[j]);
         }
         fprintf( out, "\n" );
    }
    return;
}

/* symm_eigen_solve():  function used for determining the eigenvalues and eigenvectors
 * of symmetric matrices.
 */
int symm_eigen_solve( double mat[][3], struct real_eigen_data *ret )
{
    int i, j;
    int n = N_DIM, lda = LDA, info;

    int storage_flag = STORE_UPPER_TRIANGULAR_MATRIX;
    char lapacke_storage_flag;

    double eig_values[N_DIM] = {0.0};
    double eig_vectors[N_DIM][N_DIM] = {{0.0}};

    double *a= NULL;

/*  Convert the 3x3 matrix form used by COSET into a form usable
 *  by the LAPACKE functions.  Either of the macros (see above):
 *  STORE_UPPER_TRIANGULAR_MATRIX  -or-
 *  STORE_LOWER_TRIANGULAR_MATRIX
 *
 *  may be used. We choose to use the upper triangular matrix form.
 */
    errno = 0;
    a = convert_3x3_matrix( mat, storage_flag );
    if( NULL == a ) {
        return errno;
    }
/*  LAPACKE specific code used to solve for eigenvalues and eigenvectors
 */
    if( STORE_UPPER_TRIANGULAR_MATRIX == storage_flag ) {
        lapacke_storage_flag = 'U';
    }
    else {
        lapacke_storage_flag = 'L';
    }

    /* OK, we solve for the eigenvalues and eigenvectors! */
    info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', lapacke_storage_flag, n, a, lda, eig_values );

    /* if 'info' is greater than zero, then the eigen solver didn't work correctly.
     * We crash hard if this doesn't work.
     */
    if( info > 0 ) {
        fprintf( stderr, "LAPACKE_dsyev() failed to compute eigenvalues\n" );
        exit( EXIT_FAILURE );
    }    


/*  After solution of the eigen problem, we now bring the eigenvalues and
 *  eigenvectors into a form usable by the COSET program.
 *
 *  In 'a', the eigenvectors are stored as columns.  We bring into the 
 *  static 'eig_vector' matrix and then transpose it. Only then do we
 *  copy the data into the 'ret' struct.
 */
    reform_3x3_matrix( a, eig_vectors );
    transpose_matrix( eig_vectors );

    for( i = 0; i < 3; i++ ) {
         ret[i].eig_value = eig_values[i];
         for( j = 0; j < 3; j++ )
              ret[i].eig_vector[j] = eig_vectors[i][j];
    }

/* Deallocate the dynamically allocated memory used.
 */
    free( a );

    return 0;
}


/* non_symm_eigen_solve(): function used for determining the eigenvalues and eigenvectors
 * of nonsymmetric matrices.
 */
int non_symm_eigen_solve( double mat[][3], struct complex_eigen_data *ret )
{
    int n = N_DIM, lda = LDA, ldvr = LDVR, info = 0;
    double wr[N_DIM] = {0.0}, /* wr, wi used for storing the real and imaginary components */
           wi[N_DIM] = {0.0}; /* of the eigenvalues produced by dgeev() */
    double vr[N_DIM * LDVR] = {0.0}; /* vr used for storing right eigenvectors produced by dgeev() */
    double *a = NULL; /* the dynamically allocated array used by LAPACKE */

   
    errno = 0;
/*  Convert the 3x3 matrix form used by COSET into a form usable
 *  by the LAPACKE function.
 */

    a = convert_3x3_matrix( mat, STORE_FULL_MATRIX );
    if( NULL == a ) {
        fprintf( stderr, "%s: %d: %s\n", __FILE__, __LINE__, errno != 0 ? strerror(errno) : "memory allocation error." );
        exit( EXIT_FAILURE );
    }

/*  LAPACKE specific code used to solve for eigenvalues and eigenvectors.
 *  Use LAPACKE_dgeev() function.
 * calculate only right eigenvectors
 */
    info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'V', n, a, lda, wr, wi, NULL, 3, vr, ldvr  );

    /* we crash hard if this function gives erroneous results (i.e info > 0 ) */
    if( info > 0 ) {
        fprintf( stderr, "LAPACKE_dgeev() failed to compute eigenvalues\n" );
        exit( EXIT_FAILURE );
    }

    /* marshal the eigenvalues and eigenvectors returned by dgeev() and bring them into a form usable by the 
     * COSET program. The 'ret' variable are where the eigenvalues and eigenvectors are returned
     * back to the calling function.
     */
    marshal_eigen_results( ret, n, wr, wi,  vr, ldvr );

/* Deallocate the dynamically allocated memory used.
 */
    free( a );

    return 0;
}
