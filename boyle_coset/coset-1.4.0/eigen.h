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
#ifndef EIGEN_H
#define EIGEN_H

#include <complex.h>

/* the two structs below are for storing the eigenvalues and eigenvectors
 * calculated by the third party linear algebra libraries.
 */

struct real_eigen_data {
       double eig_value;
       double eig_vector[3];
       };

struct complex_eigen_data {
       double complex eig_value;
       double complex eig_vector[3];
       };


/* is_real_eigen_value(): returns 1 if eigenvalue or has no imaginary component or 0 if they do */
int is_real_eigen_value( struct complex_eigen_data *d );


/* convert_complex_eig_to_real_eig():  this function is a convenience function to bring complex eigenvalues
 * which happen to have all zeros for the imaginary components into a struct real_eigen_data structure.
 */
int convert_complex_eig_to_real_eig( struct complex_eigen_data *in, struct real_eigen_data *out );

/* symm_eigen_solve():  function used for determining the eigenvalues and eigenvectors
 * of symmetric matrices.
 */
int symm_eigen_solve( double mat[][3], struct real_eigen_data *ret );

/* non_symm_eigen_solve(): function used for determining the eigenvalues and eigenvectors
 * of nonsymmetric matrices.
 */
int non_symm_eigen_solve( double mat[][3], struct complex_eigen_data *ret );


/* for printing out the raw eigenvalues and eigenvectors.  Mostly useful for debugging */
void print_eigen_results( FILE *out, char *desc, struct complex_eigen_data *eig, int n_dim );


#endif
