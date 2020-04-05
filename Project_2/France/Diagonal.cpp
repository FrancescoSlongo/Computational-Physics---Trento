// Computational Physics 2020 - Project 1: Hydrogen Atom with a Gaussian Basis
// Master course in Physics - Università di Trento

// We use some matricial diagonalization in order to estimate the groundstate energy of the hydrogen atom

#include <iostream>	        // input and output
#include <fstream>	        // Library for writing in a file (ofile etc.)
#include <cmath>		    // math functions
#include <string>		    // for using strings
#include <iomanip>	        // for setprecision, setw

#include <gsl/gsl_math.h>   // gsl library for math operations
#include <gsl/gsl_eigen.h>  // gsl library for eigenvalues
#include <gsl/gsl_blas.h>   // gsl library for linear algebra

using namespace std;

double pi = 4*atan(1);

// Function to cout the vectors in terminal
void vec_out(int d, gsl_vector * A)
{
    for(int i = 0; i < d; i++)
    {
        cout << gsl_vector_get(A, i) << "  ";
    }
    cout << endl;
}

// Function to cout the matrices in terminal
void mat_out(int d, gsl_matrix * A)
{
    for(int i = 0; i < d; i++)
    {
        for(int j = 0; j < d; j++)
        {
            cout << gsl_matrix_get(A, i, j) << "  ";
        }
        cout << endl;
    }
}

int main()
{
    // Dimensionality of the problem
    int d = 3;
    double alpha[d] = {0.112618, 0.399711, 2.00251};

    // I allocate the memory for the matrices
    gsl_matrix *S = gsl_matrix_alloc(d, d); // Overlap
    gsl_matrix *H = gsl_matrix_alloc(d, d); // Hamiltonian

    // I define the matrices that I want to use
    double mat1[d][d] = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    double mat2[d][d] = {{1, 3, 0}, {3, -1, 2}, {0, 2, -2}};

    // I set the values of the matrices
    for(int j = 0; j < d; j++)
    {
        for(int i = 0; i < d; i++)
        {
            // If you want to test with matrices
            //gsl_matrix_set(S, i, j, mat1[i][j]);
            //gsl_matrix_set(H, i, j, mat2[i][j]);
            // If you want to solve real problem
            gsl_matrix_set(S, i, j, pow((pi/(alpha[i]+alpha[j])), 1.5));
            gsl_matrix_set(H, i, j, 3.0*pow(pi, 1.5)*alpha[i]*alpha[j] / pow(alpha[i]+alpha[j], 2.5) - 2.0 * pi / (alpha[i] + alpha[j]));
        }
    }

    // I plot the matrices
    cout << "Matrix S" << endl;
    mat_out(d, S);
    cout << "Matrix H" << endl;
    mat_out(d, H);

    // I get the eigenvalues of the matrix S
    // I use the gsl subroutines for symmetric matrices

    // I allocate the space for the eigenvectors and the eigenvalues
    gsl_matrix *eigVec = gsl_matrix_alloc(d, d); // Eigenvectors matrix
    gsl_vector *eigVal = gsl_vector_alloc(d);    // Eigenvalues vector
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(d); // Workspace required for the eigenvalues

    // I get the eigenvalues of the eigenvalues
    gsl_eigen_symmv(S, eigVal, eigVec, w);

    // Now I build the matrix V, the one needed for the generalised eigenvalue problem
    gsl_matrix *V = gsl_matrix_alloc(d, d);

    // I set the values of the matrices
    for(int j = 0; j < d; j++)
    {
        for(int i = 0; i < d; i++)
        {
            gsl_matrix_set(V, i, j, gsl_matrix_get(eigVec, i, j)/sqrt(gsl_vector_get(eigVal,j)));
        }
    }
    cout << "V" << endl;
    mat_out(d, V);

    gsl_matrix_free(eigVec);     // I do not need the matrix for eigenvectors anymore

    // Now I calculate the matrix H'
    // I allocate the memory for the first
    gsl_matrix *H1 = gsl_matrix_alloc(d, d);        // New H matrix
    gsl_matrix *Htmp = gsl_matrix_alloc(d, d);      // Temporary H matrix

    // Matrix multiplication
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, H, 0, Htmp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Htmp, V, 0, H1);

    cout << "H1" << endl;
    mat_out(d, H1);

    gsl_matrix_free(Htmp);   // Free Htmp

    // Now I get the eigenvalues and the eigenvectors of H'
    gsl_matrix *C = gsl_matrix_alloc(d, d); // Eigenvectors matrix
    gsl_vector *E = gsl_vector_alloc(d);    // Eigenvalues vector
    gsl_matrix *Ctmp = gsl_matrix_alloc(d, d);      // Temporary C matrix

    // I get the eigenvalues of the eigenvalues
    gsl_eigen_symmv(H1, E, Ctmp, w);

    // Eigenvalues are not sorted from lower to bigger, I sort them because it is easier in this way
    // (and I need to sort the eigenvector matrix as well)
    gsl_eigen_symmv_sort(E, Ctmp, GSL_EIGEN_SORT_VAL_ASC);

    // Now I do not need the workspace anymore
    gsl_eigen_symmv_free(w);

    // Now I get the true eigenvectors
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, Ctmp, 0, C);

    cout << "Our algorithm" << endl;
    // Now I out the eigenvalues
    cout << "Energies" << endl;
    vec_out(d, E);

    // Now I out the eigenvectors
    // First I have to normalise them
    double A; // Normalization factor
    for(int j = 0; j < d; j++)
    {
        A = 0;
        for(int i = 0; i < d; i++)
        {
            A += gsl_matrix_get(C, i, j)*gsl_matrix_get(C, i, j);
        }
        for(int i = 0; i < d; i++)
        {
            gsl_matrix_set(C, i, j, gsl_matrix_get(C, i, j)/sqrt(A));
        }
    }
    cout << "Eigenvectors" << endl;
    mat_out(d, C);


    // I set the values of the matrices
    for(int j = 0; j < d; j++)
    {
        for(int i = 0; i < d; i++)
        {
            // If you want to test with matrices
            //gsl_matrix_set(S, i, j, mat1[i][j]);
            //gsl_matrix_set(H, i, j, mat2[i][j]);
            // If you want to solve real problem
            gsl_matrix_set(S, i, j, pow((pi/(alpha[i]+alpha[j])), 1.5));
            gsl_matrix_set(H, i, j, 3.0*pow(pi, 1.5)*alpha[i]*alpha[j] / pow(alpha[i]+alpha[j], 2.5) - 2.0 * pi / (alpha[i] + alpha[j]));
        }
    }

    // Now I use the dedicated subroutines of gsl
    // I allocate the space for the eigenvectors and the eigenvalues
    gsl_matrix *CGsl = gsl_matrix_alloc(d, d); // Eigenvectors matrix
    gsl_vector *EGsl = gsl_vector_alloc(d);    // Eigenvalues vector
    gsl_eigen_gensymmv_workspace * wGsl = gsl_eigen_gensymmv_alloc(d); // Workspace required for the eigenvalues

    // I get the eigenvalues of the eigenvalues
    gsl_eigen_gensymmv(H, S, EGsl, CGsl, wGsl);

    // Now I do not need the workspace anymore
    gsl_eigen_gensymmv_free(wGsl);

    // Now I out the result
    cout << "GSL" << endl;
    // Now I out the eigenvalues
    cout << "Energies" << endl;
    vec_out(d, EGsl);
    // Now I out the eigenvectors
    cout << "Eigenvectors" << endl;
    mat_out(d, CGsl);

    return 0;
}
