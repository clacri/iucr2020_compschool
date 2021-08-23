/* contains implementation for COSET's eigen interface (eigen.h)
 * Copyright (C) 2008, 2012, 2013, 2017 Paul D. Boyle
 *
 * Written by:
 *      Paul D. Boyle
 *      Department of Chemistry
 *      University of Western Ontario
 *      London, Ontario CANADA N6A 5B7
 *
 *      November 2017
 * with modifications made June 2021     
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>

#include "eigen.h"
#include "float_util.h"

/* we write and use custom error handler to avoid GSL default behaviour of calling abort()
 * and dumping a core file.
 */

static void coset_gsl_error_handler( const char *reason, const char *src_file, int src_line, int gsl_err_num )
{
    fprintf( stderr, "Error: %s (near line %d): %s: %s\n", src_file, src_line, reason, gsl_strerror(gsl_err_num) );
    exit( EXIT_FAILURE );
}


/* reform_matrix(): takes a 3x3 symmetry matrix and puts into a form (gsl_matrix *) suitable for being
 * passed to GSL functions.  Caller must gsl_matrix_free() to free the allocated memory.
 */
static gsl_matrix *reform_matrix( double m[][3] )
{

    int i, j;
    gsl_matrix *mp = NULL;

    errno = 0;
    mp = gsl_matrix_alloc( 3, 3 );
    if( NULL == mp )
        return NULL;

    for( i = 0; i < 3; i++ ) {
         for( j = 0; j < 3; j++ ) {
              gsl_matrix_set( mp, i, j, m[i][j] );
         }
    }

    return mp;
}


/* *** Public Functions *** */

/* a macro for printing complex numbers see:
 * stackoverflow.com/questions/4099433/c-complex-number-and-printf
 * used in the print_eigen_results() function below.
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



/* symm_eigen_solve():  function used for determining the eigen values and eigen vectors
 * of symmetric matrices.
 */


int symm_eigen_solve( double mat[][3], struct real_eigen_data *ret )
{
    int i, j;
    int gsl_status = GSL_SUCCESS;

    gsl_matrix *m = NULL;

/* GSL data types and allocations */
    gsl_vector *eval = gsl_vector_alloc( 3 );
    gsl_vector *evec_col = gsl_vector_alloc( 3 );
    gsl_matrix *evec = gsl_matrix_alloc( 3, 3 );
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc( 3 );

/* install our own GSL error handler callback function */

    gsl_set_error_handler( &coset_gsl_error_handler );

/*  Convert the 3x3 matrix form used by COSET into a form usable
 *  by the GSL functions.
 */
    m = reform_matrix( mat );

/*  GSL specific code used to solve for eigenvalues and eigenvectors
 */

    gsl_status = gsl_eigen_symmv( m, eval, evec, w );
    if(GSL_SUCCESS != gsl_status ) {
          GSL_ERROR( "gsl_eigen_symmv() failed", gsl_status );
    }

    gsl_eigen_symmv_free( w );
    gsl_eigen_symmv_sort( eval, evec, GSL_EIGEN_SORT_ABS_ASC );


/*  After solution of the eigen problem, we now bring the eigenvalues and
 *  eigenvectors into a form usable by the COSET program.
 */
    for( i = 0; i < 3; i++ ) {
         ret[i].eig_value = gsl_vector_get( eval, i );
         gsl_matrix_get_col(evec_col, evec, i );
         for( j = 0; j < 3; j++ )
              ret[i].eig_vector[j] = gsl_vector_get( evec_col, j);
    }

/* Deallocate the dynamically allocated memory used.
 */
    gsl_vector_free( eval );
    gsl_vector_free( evec_col );
    gsl_matrix_free( evec );
    gsl_matrix_free( m );

    return 0;
}


/* non_symm_eigen_solve(): function used for determining the eigenvalues and eigenvectors
 * of nonsymmetric matrices.
 */
int non_symm_eigen_solve( double mat[][3], struct complex_eigen_data *ret )
{
    int i, j;
    int gsl_status = GSL_SUCCESS;
    gsl_matrix *m = NULL;

/* GSL data types used and allocated below */
    gsl_complex eig_val;
    gsl_vector_complex_view eig_vec;
    gsl_complex z;
   
    gsl_vector_complex *eval = gsl_vector_complex_alloc( 3 );
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc( 3, 3 );
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc( 3 );

/* install our own custom GSL error handler function */

    gsl_set_error_handler( &coset_gsl_error_handler );

/*  Convert the 3x3 matrix form used by COSET into a form usable
 *  by the GSL functions.
 */
    m = reform_matrix( mat );

/*  GSL specific code used to solve for eigenvalues and eigenvectors
 */

    gsl_status = gsl_eigen_nonsymmv( m, eval, evec, w );
    if(GSL_SUCCESS != gsl_status ) {
          GSL_ERROR( "gsl_eigen_nonsymmv() failed", gsl_status );
    }

    gsl_eigen_nonsymmv_free( w );
    gsl_eigen_nonsymmv_sort( eval, evec, GSL_EIGEN_SORT_ABS_DESC );

/*  After solution of the eigen problem, we now bring the eigenvalues and
 *  eigenvectors into a form usable by the COSET program.
 */
    for( i = 0; i < 3; i++ ) {
         eig_val = gsl_vector_complex_get( eval, i );
         ret[i].eig_value = GSL_REAL(eig_val) + GSL_IMAG(eig_val) * I;
         eig_vec = gsl_matrix_complex_column( evec, i );
         for( j = 0; j < 3; j++ ) {
              z = gsl_vector_complex_get( &eig_vec.vector, j );
              ret[i].eig_vector[j] = GSL_REAL(z) + GSL_IMAG(z) * I;
         }
    }

/* Deallocate the dynamically allocated memory used.
 */
    gsl_vector_complex_free( eval );
    gsl_matrix_complex_free( evec );
    gsl_matrix_free( m );

    return 0;
}
