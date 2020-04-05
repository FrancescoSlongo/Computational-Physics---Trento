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

double diag(int d, double * alpha)
{
    // I allocate the memory for the matrices
    gsl_matrix *S = gsl_matrix_alloc(d, d); // Overlap
    gsl_matrix *H = gsl_matrix_alloc(d, d); // Hamiltonian
    // Now I get the eigenvalues and the eigenvectors of H'
    gsl_matrix *C = gsl_matrix_alloc(d, d); // Eigenvectors matrix
    gsl_vector *E = gsl_vector_alloc(d);    // Eigenvalues vector

    // I set the values of the matrices
    for(int j = 0; j < d; j++)
    {
        for(int i = 0; i < d; i++)
        {
            // If you want to solve real problem
            gsl_matrix_set(S, i, j, pow((pi/(alpha[i]+alpha[j])), 1.5));
            gsl_matrix_set(H, i, j, 3.0*pow(pi, 1.5)*alpha[i]*alpha[j] / pow(alpha[i]+alpha[j], 2.5) - 2.0 * pi / (alpha[i] + alpha[j]));
        }
    }

    // Now I use the dedicated subroutines of gsl
    // I allocate the space for the eigenvectors and the eigenvalue
    gsl_eigen_gensymmv_workspace * wGsl = gsl_eigen_gensymmv_alloc(d); // Workspace required for the eigenvalues

    // I get the eigenvalues of the eigenvalues
    gsl_eigen_gensymmv(H, S, E, C, wGsl);

    // Now I free S and H
    gsl_matrix_free(S);
    gsl_matrix_free(H);
    gsl_matrix_free(C);
    // Now I do not need the workspace anymore
    gsl_eigen_gensymmv_free(wGsl);
    return gsl_vector_get(E, 0);
}

double steep_desc(int d, double h, double gamma, double resolution, double* alpha)
{
    double alpha0[d] = {0, 0, 0};
    double der[d];
    // Variables for derivative
    double alpha1[d];
    double alpha2[d];
    double alpha3[d];
    double alpha4[d];
    int MaxCyc = 100000;

    double grad;        // Gradient for if condition

    int m = 0;          // Counter for cycle
    while(m < MaxCyc)
    {
        gamma = gamma/(m/10000 + 1.0);
        for(int i = 0; i < d; i++)
        {
            alpha0[i] = alpha[i];
        }
        // I set the derivative at 5 points
        grad = 0;
        for(int i = 0; i < d; i++)
        {
            for(int j = 0; j < d; j++)
            {
                if(i != j)
                {
                    alpha1[j] = alpha0[j];
                    alpha2[j] = alpha0[j];
                    alpha3[j] = alpha0[j];
                    alpha4[j] = alpha0[j];
                } else {
                    alpha1[j] = alpha0[j] - 2*h;
                    alpha2[j] = alpha0[j] - h;
                    alpha3[j] = alpha0[j] + h;
                    alpha4[j] = alpha0[j] + 2*h;
                }
            }

            der[i] = (diag(d, alpha1) - 8.0 * diag(d, alpha2) + 8.0 * diag(d, alpha3) - diag(d, alpha4)) / (12.0 * h);
            alpha[i] -= gamma*der[i];
            grad += der[i]*der[i];
        }
        m++;
    }
    cout << "Energies" << endl;
    cout << diag(d, alpha) << endl;
    cout << "Alpha" << endl;
    for(int i = 0; i < d; i++)
    {
        cout << alpha[i] << "\t";
    }
    return diag(d, alpha);
}

int main()
{
    // Dimensionality of the problem
    int d = 3;
    double alpha[d] = {0.120377, 0.4455, 2.30887};
    double h = 0.000001;
    double gamma = 0.001;
    double resolution = 0.0000001;
    // Now I get the eigenvalues and the eigenvectors of H'
    double E;    // Eigenvalues vector


    // Now I out the result
    cout << "GSL" << endl;
    // Now I out the eigenvalues
    cout << "Energies" << endl;
    E = diag(d, alpha);
    cout << E << endl;
    cout << "Steepest descent" << endl;
    E = steep_desc(d, h, gamma, resolution, alpha);

    return 0;
}
