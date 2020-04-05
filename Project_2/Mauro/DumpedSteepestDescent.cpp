#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_multimin.h>

using namespace std;

const int n = 3;
gsl_matrix *genEvecs = gsl_matrix_alloc(n, n);
gsl_vector *genEvals = gsl_vector_alloc(n);
gsl_matrix *Sqp = gsl_matrix_alloc(n, n); // Overlap
gsl_matrix *Hqp = gsl_matrix_alloc(n, n); // Hamiltonian
double alphas[n], val, val2;

double my_f(const gsl_vector *v, void *params)
{
    double *p = (double *)params;
    int i = 0;
    while (i < n)
    {
        alphas[i] = gsl_vector_get(v, i) + p[i + 1];
        i++;
    }

    i = 0;
    int j;
    while (i < n) // Just set matrix elements
    {
        j = 0;
        while (j < n)
        {
            /*if (i < 3 && j < 3)
            {
                val = M_PI/(alphas[i]+alphas[j]) * sqrt(M_PI/(alphas[i]+alphas[j]));
                val2 = 3.*pow(M_PI, 1.5)*alphas[i]*alphas[j]/pow(alphas[i]+alphas[j], 2.5) - 2*M_PI/(alphas[i]+alphas[j]);
            }
            else if (i >= 3 && j >= 3)
            {
                val2 = M_PI * sqrt(M_PI) / 2. / sqrt(alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
                val2 = 5. * M_PI * sqrt(M_PI) / 2. * alphas[i]*alphas[j] / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / sqrt(alphas[i]+alphas[j]);
                val2 -= 2.*M_PI/3. / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
            }
            else
            {
                val = 0;
                val2 = 0;
            }*/
            //val = M_PI/(alphas[i]+alphas[j]) * sqrt(M_PI/(alphas[i]+alphas[j]));
            val = M_PI * sqrt(M_PI) / 2. / sqrt(alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
            gsl_matrix_set(Sqp, i, j, val);
            //val2 = 3.*pow(M_PI, 1.5)*alphas[i]*alphas[j]/pow(alphas[i]+alphas[j], 2.5) - 2*M_PI/(alphas[i]+alphas[j]);
            val2 = 5. * M_PI * sqrt(M_PI) / 2. * alphas[i]*alphas[j] / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]) / sqrt(alphas[i]+alphas[j]);
            val2 -= 2.*M_PI/3. / (alphas[i]+alphas[j]) / (alphas[i]+alphas[j]);
            gsl_matrix_set(Hqp, i, j, val2);
            j++;
        }
        i++;
    }

    gsl_eigen_gensymmv_workspace *genws = gsl_eigen_gensymmv_alloc(n); // Don't actually care of eigenvectors but you never know
    gsl_eigen_gensymmv(Hqp, Sqp, genEvals, genEvecs, genws);
    gsl_eigen_gensymmv_free(genws);

    return gsl_vector_min(genEvals);
}

void my_df(const gsl_vector *v, void *params, double *df)
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
    while (i < n)
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
        p[1 + i] = 0.;

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
        p[i + 1] = 0.;
        i++;
    }
    p[0] = 1e-8;

    // Starting point
    gsl_vector *x = gsl_vector_alloc(n);
    double alphav[n] = {0.3, 1.6, 1.};
    //double alphav[n] = {0.3, 1.6};
    alphav[0] = 1.5;
    alphav[1] = 1.;
    alphav[2] = .5;
    //double alphav[n] = {1.};
    //double alphav[n] = {0.3, 0.6, 1., 1.2, 1.4, 1.6};
    i = 0;
    while (i < n)
    {
        gsl_vector_set(x, i, alphav[i]);
        i++;
    }

    double step = 1e-3;
    const int iter = 1e6;
    int j;
    static double stepVec[iter];
    double gradThre = 1e-10;
    double gradAbs = 1;
    double grad[n];

    i = 0;
    while (i < 1e4)
    {
        stepVec[i] = 1e-3;
        i++;
    }
    while (i < 1e5)
    {
        stepVec[i] = 1e-4;
        i++;
    }
    while (i < iter)
    {
        stepVec[i] = 1e-5;
        i++;
    }

    i = 0;
    while (i < iter && gradAbs > gradThre)
    {
        step = stepVec[i];
        gradAbs = 0;
        my_df(x, (void*)p, grad);

        // Calculate immediately abs of grad
        j = 0;
        while (j < n)
        {
            gradAbs += grad[j]*grad[j];
            j++;
        }
        gradAbs = sqrt(gradAbs);

        j = 0;
        while (j < n)
        {
            alphav[j] -= step*grad[j]/gradAbs;
            gsl_vector_set(x, j, alphav[j]);
            //gradAbs += grad[j]*grad[j];
            j++;
        }
        //gradAbs = sqrt(gradAbs);
        //cout << "Energy: " << my_f(x, (void*)p) << endl;
        i++;
    }
    int itDone = i--;

    cout << "Values: ";
    i = 0;
    while (i < n)
    {
        cout << alphav[i] << " ";
        i++;
    }
    cout << endl;
    cout << "With an energy of: " << setprecision(15) << my_f(x, (void*)p) << endl;
    cout << "And gradient: " << gradAbs << endl;
    cout << "After " << itDone << " passages." << endl;
}
