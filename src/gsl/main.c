#include <stdio.h>
#include <gsl/gsl_linalg.h>
     
int main (void)
     {
     double a_data[] = { 1.0, 0.6, 0.0,
                         0.0, 1.5, 1.0,
                         0.0, 1.0, 1.0 };
     /*
      * Inverse is
      *    1  -1.2   1.2
      *    0   2.0  -2.0
      *    0  -2.0   3.0
      */
     double inva[9];
     int s, i, j;
     
     gsl_matrix_view m
          = gsl_matrix_view_array(a_data, 3, 3);
     gsl_matrix_view inv
          = gsl_matrix_view_array(inva,3,3);
     gsl_permutation * p = gsl_permutation_alloc (3);
 
     printf("The matrix is\n");
     for (i = 0; i < 3; ++i)
         for (j = 0; j < 3; ++j)
             printf(j==2?"%6.3f\n":"%6.3f ", gsl_matrix_get(&m.matrix,i,j));
 
     gsl_linalg_LU_decomp (&m.matrix, p, &s);    
     gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
 
     printf("The inverse is\n");
     for (i = 0; i < 3; ++i)
         for (j = 0; j < 3; ++j)
             printf(j==2?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
 
     gsl_permutation_free (p);
     return 0;
     }
