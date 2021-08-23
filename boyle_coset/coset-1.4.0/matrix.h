/* public interface for matrix operations for 3x3 matricies
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

#ifndef MATRIX_H
#define MATRIX_H 
#include <stdio.h>

#define N_DIM 3

double determinant(double m[N_DIM][N_DIM]);
double trace( double m[N_DIM][N_DIM] );
void invert_matrix(double det, double m[N_DIM][N_DIM], double inv[N_DIM][N_DIM]);
void negate_matrix( double negated[N_DIM][N_DIM], double m[N_DIM][N_DIM] );
void transpose_matrix(double m[N_DIM][N_DIM]);
void matrix_multiply3x3(double prod[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM]);
void matrix_add3x3(double prod[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM]);
void matrix_subtract3x3(double prod[N_DIM][N_DIM], double a[N_DIM][N_DIM], double b[N_DIM][N_DIM]);
void calculate_inverse_transpose(double inv_tr[N_DIM][N_DIM], double matrx[N_DIM][N_DIM]);
void similarity_transform( double transformed[3][3], double mat[3][3], double trans_mat[3][3] );
void copy_matrix( double dest[N_DIM][N_DIM], double src[N_DIM][N_DIM] );
#ifdef USE_DIALOG
void get_matrix(double m[N_DIM][N_DIM]);
#endif
void print_matrix(FILE *out, double m[N_DIM][N_DIM], const char *format);

#endif
