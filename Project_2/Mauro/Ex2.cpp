#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

using namespace std;

int main (void)
{
    const int n = 3; // Dimensionality
    double alpha[n] = {0.112618, 0.399711, 2.00251}; // parameters
    double mat[n][n] = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}; // Just some test matrices (mat MUST BE POSITIVE DEFINITE)
    double mat2[n][n] = {{1, 3, 0}, {3, -1, 2}, {0, 2, -2}};
    gsl_matrix *Sqp = gsl_matrix_alloc(n, n); // Overlap
    gsl_matrix *Hqp = gsl_matrix_alloc(n, n); // Hamiltonian
    gsl_matrix *evecs = gsl_matrix_alloc(n, n); // Sqp eigenvectors
    gsl_vector *evals = gsl_vector_alloc(n); // Sqp eigenvalues
    double val;
    for (int i = 0; i < n; i++) // Just set matrix elements
    {
        for (int j = 0; j < n; j++)
        {
            //val = M_PI/(alpha[i]+alpha[j]);
            val = mat[i][j];
            gsl_matrix_set(Sqp, i, j, pow((M_PI/(alpha[i]+alpha[j])), 1.5));
            val = mat2[i][j];
            gsl_matrix_set(Hqp, i, j, 3.0*pow(M_PI, 1.5)*alpha[i]*alpha[j] / pow(alpha[i]+alpha[j], 2.5) - 2.0 * M_PI / (alpha[i] + alpha[j]));
            cout << val << " ";
        }
        cout << endl;
    }

    // Diagonalize S
    gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(Sqp, evals, evecs, ws);
    gsl_eigen_symmv_free(ws);
    gsl_eigen_symmv_sort(evals, evecs, GSL_EIGEN_SORT_VAL_ASC);

    // Build Sm12
    gsl_matrix *Sm12 = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(Sm12);
    for (int i = 0; i < n; i++)
    {
        val = 1. / sqrt(abs(gsl_vector_get(evals, i)));
        gsl_matrix_set(Sm12, i, i, val);
    }

    //Product between U (evecs) and Sm12
    gsl_matrix *V = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(V);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., evecs, Sm12, 0., V); // Not needed, can just do manually but it was for practise

    //Product between H and V
    gsl_matrix *tempH = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(tempH);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Hqp, V, 0., tempH);

    //Product between Vt and tempH
    gsl_matrix *newH = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(newH);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., V, tempH, 0., newH);

    /*for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << gsl_matrix_get(newH, i, j) << " ";
        }
        cout << endl;
    }*/

    // Diagonalize newH
    gsl_matrix *Cp = gsl_matrix_alloc(n, n);
    gsl_vector *epsilon = gsl_vector_alloc(n);
    gsl_eigen_symmv_workspace *ws2 = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv(newH, epsilon, Cp, ws2);
    gsl_eigen_symmv_free(ws2);
    //gsl_eigen_symmv_sort(epsilon, Cp, GSL_EIGEN_SORT_VAL_ASC); // Sort them NOT! THIS IS WRONG

    // Calculate true eigenvectors
    gsl_matrix *C = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(C); // For good measure
    //Product between V and Cp
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Cp, 0., C);

    double norm; // Normalize columns
    for (int i = 0; i < n; i++)
    {
        norm = 0.;
        for (int j = 0; j < n; j++)
        {
            val = gsl_matrix_get(C, j, i);
            norm += val*val;
        }
        norm = sqrt(norm);
        for (int j = 0; j < n; j++)
        {
            val = gsl_matrix_get(C, j, i);
            gsl_matrix_set(C, j, i, val/norm);
        }
    }

    cout << "With long procedure:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << gsl_vector_get(epsilon, i) << endl;
        for (int j = 0; j < n; j++)
        {
            cout << gsl_matrix_get(C, i, j) << " ";
        }
        cout << endl;
    }

    // Using dedicated subroutine
    gsl_matrix *genEvecs = gsl_matrix_alloc(n, n);
    gsl_vector *genEvals = gsl_vector_alloc(n);
    gsl_eigen_gensymmv_workspace *genws = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(Hqp, Sqp, genEvals, genEvecs, genws);
    gsl_eigen_gensymmv_free(genws);

    cout << "With dedicated subroutine" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << gsl_vector_get(genEvals, i) << endl;
        for (int j = 0; j < n; j++)
        {
            cout << gsl_matrix_get(genEvecs, i, j) << " ";
        }
        cout << endl;
    }

    /*// Some check
    gsl_matrix *test = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(test);
    gsl_blas_dgemm(CblasConjTrans, CblasNoTrans, 1., V, V, 0., test);
    cout << "What he calls orthogonal:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << gsl_matrix_get(test, i, j) << " ";
        }
        cout << endl;
    }*/
}
