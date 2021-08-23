/* contains implementation for used for the Flack left coset decomposition
 * algorithms as outlined in Acta Cryst. (1987), A43, 564-568, by H. D.
 * Flack.
 *
 * If you want to use the extended B algorithm compile with 
 * USE_EXTENDED_B_ALGORITHM defined 
 * Copyright (C) 2008, 2012, 2013 Paul D. Boyle
 *
 * Written by:
 *      Paul D. Boyle
 *      Department of Chemistry
 *      University of Western Ontario
 *      London, Ontario CANADA N6A 5B7
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

#include "symm_mat.h"
#include "matrix.h"
#include "coset.h"

#ifdef USE_ELEMENT_COMPARISON  /* normally not used, we compare Bit eCoded Matrcies (bcm) */
/* if you want to use an element by element comparison, you will need integrate the 
 * is_equal_matrices() function below into the algorithms yourself.
 */
#include "float_util.h"   /*  you must add float_util.c to Makefile if not already there. */

#define PERFECT_SCORE 9  /* All 9 elements must match up */

/* this function is used when we want to do an element by element comparison of each matrix.
 * the function pointer 'compare' should be a function like the 'is_equal()' function
 * in float_util.c, or you can code your own according to the prototype.  
 */
int is_equal_matrices( double m1[3][3], double m2[3][3], int (*compare)( double, double ) )
{
    int score = 0;  /* will equal 9 if matricies are identical */
    int i, j;
    for( i = 0; i < 3; i++ )
         for( j = 0; j < 3; j++ )
              score += compare( m1[i][j], m2[i][j] );
 
    return PERFECT_SCORE == score ? 1 : 0; 
}
#endif /* end of USE_ELEMENT_COMPARISON conditional compilation directive */

/* Here is Flack's algorithm A
 */
void coset_decomposition_A( struct symm_op *G, struct symm_op *H )
{
   int i, j;
   struct symm_op *G_k;

   for( i = 0;  0 != G[i].bcm; i++ ) {
        if( True == G[i].truefalse ) {
            for( j = 1; 0 != H[j].bcm; j++ ) {
                 struct symm_op prod =  {False,{{0.0}},0,0,0.0,0.0,{0.0,0.0,0.0}};
                 matrix_multiply3x3( prod.mat, G[i].mat, H[j].mat ); 
                 prod.bcm = encode_matrix( prod.mat );
                 for( G_k = &G[i+1]; 0 != G_k->bcm; G_k++ ) {
                      if( prod.bcm == G_k->bcm ) {  /* bit coded matricies compared */
                          G_k->truefalse = False;
                      }
                 }
            }
        }
   }
   return;

}

/* Here is Flack's algorithm B
 */
void coset_decomposition_B( struct symm_op *G, struct symm_op *H )
{
   int i, j,
       centric_flag = 0,
       G_count,
       H_count;
   struct symm_op *G_k;

   G_count = count_ops( G );
   G_count /= 2;  /* we don't need to loop over the centrically related operations */

   H_count = count_ops( H );
   centric_flag = is_centric( H );
   if( Centric  == centric_flag ) {
       H_count /= 2;
   }

/* turn off the centrically related elements of G */
   for( i = G_count; 0 != G[i].bcm; i++ ) {
        G[i].truefalse = False;
   }

/* now do Algorithm B (Flack, p.567) */
   for( i = 0; i < G_count; i++ ) {
        if( True == G[i].truefalse ) {
            for( j = 1; j < H_count; j++ ) {
                 struct symm_op prod_pos = {False,{ {0.0} }, 0,0,0.0,0.0,{0.0,0.0,0.0}};
                 struct symm_op prod_neg = {False,{ {0.0} }, 0,0,0.0,0.0,{0.0,0.0,0.0}};
                 matrix_multiply3x3( prod_pos.mat, G[i].mat, H[j].mat );
                 prod_pos.bcm = encode_matrix( prod_pos.mat );
                 negate_matrix( prod_neg.mat, prod_pos.mat );
                 prod_neg.bcm = encode_matrix( prod_neg.mat );
                 for( G_k = &G[i+1]; G_k < &G[G_count]; G_k++ ) {
                      if( (prod_pos.bcm == G_k->bcm) ||
                          (prod_neg.bcm == G_k->bcm) ) {
                          G_k->truefalse = False;
                      }
                 }
                     
            }
        }
   } /* end of Flack algorithm B */
#ifdef USE_EXTENDED_B_ALGORITHM

/* this next part isn't part of the original algorithm, but I think it is
 * more convenient for the end user to list out explicitly those matricies which are 
 * implied centrically related twin laws for non-centrosymmetric
 * crystals.
 *
 * Loop over G, if the element is marked 'True', then find
 * it's centrically related pair, and mark that one true too.
 */
   if( Acentric  == centric_flag ) {
       for( i = 0; i < G_count; i++ ) {
            if( True == G[i].truefalse ) {
                double neg[3][3] = { {0.0} };
                unsigned int bcm;
                negate_matrix( neg, G[i].mat );
                bcm = encode_matrix( neg );
                for( j = G_count; 0 != G[j].bcm; j++ ) {
                     if( bcm == G[j].bcm ) {
                         G[j].truefalse = True;
                     }
                }
            
            }
       }
   }
#endif  /* end of USE_EXTENDED_B_ALGORITHM conditional compilation directive */

   return;

}

/* function for printing out the system of representatives.
 * The 'start' parameter is the index we want to start 
 * the output from.
 */
void output_cosets( FILE *out, struct symm_op *s, int start )
{
#define FMT "%6.2f%6.2f%6.2f\n"

    int i;
    fputs("\n*** Possible Twin Laws for this Subgroup-Supergroup Relationship ***\n", out );
    for( i = start; 0 != s[i].bcm; i++ ) {
         if( True == s[i].truefalse ) {
             fprintf( out, "Possible Twin Law (%d):\n", i );
             print_matrix( out, s[i].mat, FMT );
         }
    }
    fflush( out );
    return;
#undef FMT
}

