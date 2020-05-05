#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#include <gsl/gsl_complex.h> // Not necessary
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_multimin.h> // Not necessary actually

using namespace std;

const int n = 3;
gsl_matrix *genEvecs = gsl_matrix_alloc(n, n);
gsl_vector *genEvals = gsl_vector_alloc(n);
gsl_matrix *Sqp = gsl_matrix_alloc(n, n); // Overlap
gsl_matrix *Hqp = gsl_matrix_alloc(n, n); // Hamiltonian
double alphas[n], val, val2;

double my_f(const gsl_vector *v, void *params) // Solves generalized eigenvalue problem
{
    double *p = (double *)params; // Reminescent of using minimization functions of gsl library
    for (int i = 0; i < n; i++)
    {
        alphas[i] = gsl_vector_get(v, i) + p[i + 1]; // This is for derivatives (see my_df below)
    }

    for (int i = 0; i < n; i++) // Just set matrix elements
    {
        for (int j = 0; j < n; j++)
        {
            /*if (i < 3 && j < 3) // This is for mixing bases
            {
                val = M_PI/(alphas[i]+alphas[j]) * sqrt(M_PI/(alphas[i]+alphas[j]));
                val2 = 3.*pow(M_PI, 1.5)*alphas[i]*alphas[j]/pow(alphas[i]+alphas[j], 2.5) - 2*M_PI/(alphas[i]+alphas[j]);
            }
            else if (j == i && i >= 3)
            {
                val = M_PI * sqrt(M_PI) / 2. / sqrt(alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
                val2 = 5. * M_PI * sqrt(M_PI) / 2. * alphas[i]*alphas[j] / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / sqrt(alphas[i]+alphas[j]);
                val2 -= 2.*M_PI/3. / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]); // Just because it's more readable
            }
            else
            {
                val = 0;
                val2 = 0;
            }*/

            // First is for s-waves, second for p-waves
            val = M_PI/(alphas[i]+alphas[j]) * sqrt(M_PI/(alphas[i]+alphas[j]));
            //val = M_PI * sqrt(M_PI) / 2. / sqrt(alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
            gsl_matrix_set(Sqp, i, j, val);
            // First is for s-waves, second two are for p-waves
            val2 = 3.*pow(M_PI, 1.5)*alphas[i]*alphas[j]/pow(alphas[i]+alphas[j], 2.5) - 2*M_PI/(alphas[i]+alphas[j]);
            //val2 = 5. * M_PI * sqrt(M_PI) / 2. * alphas[i]*alphas[j] / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / sqrt(alphas[i]+alphas[j]);
            //val2 -= 2.*M_PI/3. / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
            gsl_matrix_set(Hqp, i, j, val2);
        }
    }

    gsl_eigen_gensymmv_workspace *genws = gsl_eigen_gensymmv_alloc(n); // Solve generalized problem, also with vectors
    gsl_eigen_gensymmv(Hqp, Sqp, genEvals, genEvecs, genws);
    gsl_eigen_gensymmv_free(genws);

    return gsl_vector_min(genEvals); // Return the smallest eigenvalue
}

void my_df(const gsl_vector *v, void *params, double *df) // 5 point derivative
{
    double *p = (double*)params; // In here there's derivation step

    int i = 0;
    while (i < n)
    {
        alphas[i] = gsl_vector_get(v, i);
        i++;
    }

    double h = p[0];

    double der;
    i = 0;
    while (i < n) // For every dimension calculate values of function at the 4 different points by modifying the parameter array
    {
        der = 0;
        p[1 + i] = -2*h;
        der += my_f(v, (void*)p);
        p[1 + i] = -h;
        der -= 8.*my_f(v, (void*)p);
        p[1 + i] = h;
        der += 8.*my_f(v, (void*)p);
        p[1 + i] = 2*h;
        der -= my_f(v, (void*)p);
        p[1 + i] = 0.; // Restore original value for parameter array

        der /= (12.*h);

        df[i] = der;

        i++;
    }
}

int main()
{
    double p[n + 1];
    int i = 0;
    while (i < n)
    {
        p[i + 1] = 0.; // To be used for derivatives (see above)
        i++;
    }
    p[0] = 1e-8; // Derivation step

    // Starting point
    gsl_vector *x = gsl_vector_alloc(n);
    // Just a bunch of initializations for different dimensions
    //double alphav[n] = {.2};
    //double alphav[n] = {0.2, 1.3};
    double alphav[n] = {.4, .9, 1.4};
    //double alphav[n] = {0.3, 0.6, 0.8, 1.6};
    //double alphav[n] = {0.3, 0.6, 0.8, 1.6, 1.2, 1.3};
    i = 0;
    while (i < n)
    {
        gsl_vector_set(x, i, alphav[i]);
        i++;
    }

    double step = 1e-1; // To multiply by the gradient
    const int iter = 1e6; // Maximum iterations of minimization
    int j = n;
    double gradThre = 1e-10; // Threshold for modulus of gradient
    double gradAbs = 1; // Modulus of gradient
    double grad[n];

    i = 0;
    while (i < iter && gradAbs > gradThre) // Minimization loop
    {
        while (j < n) // Update position (skip i = 0 because j = n, j = 0 is set at the end of this cycle)
        {
            alphav[j] -= step*grad[j];
            gsl_vector_set(x, j, alphav[j]);
            j++;
        }
        gradAbs = 0;
        my_df(x, (void*)p, grad); // Compute gradient

        j = 0;
        for (j = 0; j < n; j++) // Calculate modulus of gradient
        {
            gradAbs += grad[j]*grad[j];
        }
        gradAbs = sqrt(gradAbs);
        //cout << gradAbs << endl;
        //cout << "Energy: " << my_f(x, (void*)p) << endl;
        i++;

        j = 0;
    }
    int itDone = i--; // Iterations done

    cout << "Values: ";
    i = 0;
    while (i < n)
    {
        cout << setprecision(20) << alphav[i] - 8./9./M_PI << " ";
        i++;
    }
    cout << endl;
    cout << "With an energy of: " << setprecision(20) << my_f(x, (void*)p) + 4./3./M_PI << endl;
    cout << "And gradient: " << gradAbs << endl;
    cout << "After " << itDone << " passages." << endl;
    int idx = gsl_vector_min_index(genEvals);
    cout << "Coefficients:" << endl;
    i = 0;
    while (i < n)
    {
        cout << setprecision(20) << gsl_matrix_get(genEvecs, i, idx) << " ";
        i++;
    }
    cout << endl;
}
