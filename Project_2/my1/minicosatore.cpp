#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<cstring>
#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_siman.h>

using namespace std;

const int N = 3;

/* set up parameters for this simulated annealing run */

/* how many points do we try before stepping */
#define N_TRIES 200

/* how many iterations for each T? */
#define ITERS_FIXED_T 1000

/* max step size in random walk */
#define STEP_SIZE 1.0

/* Boltzmann constant */
#define K 1.0

/* initial temperature */
#define T_INITIAL 0.008

/* damping factor for temperature */
#define MU_T 1.003
#define T_MIN 2.0e-6

gsl_siman_params_t params  = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

void print_gsl_matrix(gsl_matrix * A, string name)
{
    cout<<"matrix "<<name<< ":"<<endl;
        for(int i=0;i<(*A).size2;i++)
    {
        for(int j=0;j<(*A).size1;j++)
        {
            cout<<gsl_matrix_get(A, i, j)<<" ";
        }
        cout<<endl;
    }

}


void print_gsl_vector(gsl_vector * v, string name)
{
    cout<<"vector "<<name<< ":"<<endl;
        for(int i=0;i<(*v).size;i++)
    {
        cout<<gsl_vector_get(v, i)<<" ";
    }
    cout<<endl;
}

//energy function
double gs_E(void * b)
{
    double * a = ((double *) b);
    gsl_matrix * S = gsl_matrix_alloc(N,N);
    gsl_matrix * H = gsl_matrix_alloc(N,N);
    gsl_eigen_gensymm_workspace *w = gsl_eigen_gensymm_alloc(N);
    gsl_vector * val = gsl_vector_alloc(N);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            gsl_matrix_set(S, i, j, pow(M_PI/(a[i]+a[j]),3./2.));
            gsl_matrix_set(H, i, j, 3.*pow(M_PI,3./2.)*a[i]*a[j]/pow(a[i]+a[j],5./2.)-2*M_PI/(a[i]+a[j]));
        }
    }


    gsl_eigen_gensymm(H, S, val, w);

    return gsl_vector_min(val);

    gsl_eigen_gensymm_free(w);
    gsl_matrix_free(S);
    gsl_vector_free(val);
}

//take step function
void t_step(const gsl_rng * r, void *xp, double step_size)
{
    double * old_x = ((double *) xp);
    double new_x[N];

    double u;
    for(int i=0; i<N;i++)
    {
        u = gsl_rng_uniform(r);
        new_x[i] = u * 2 * step_size - step_size + old_x[i];
    }

    memcpy(xp, &new_x, sizeof(new_x));
}

//distance function
double dist(void *xp, void *yp)
{
    double * x = ((double *) xp);
    double * y = ((double *) yp);
    double dd=0;

    for(int i=0;i<N;i++)
    {
        dd+=(x[i]-y[i])*(x[i]-y[i]);
    }

    return sqrt(dd);
}

//print position func
void p_pos(void *xp)
{
    for(int i=0;i<N;i++)
    {
        printf(" %.6lf  ", ((double *) xp)[i]);
    }
}

int main()
{
    /*double a[N];

    for(int i=0;i<N;i++)
    {
        cout<<"alpha i "<<i<<" -> ";
        cin>>a[i];
    }

    cout<< gs_E(a)<<endl;*/

    const gsl_rng_type * T;
    gsl_rng * r;

    double x_initial[N];
    for(int i=0;i<N;i++)
    {
        cout<<"alpha i "<<i<<" -> ";
        cin>>x_initial[i];
    }

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_siman_solve(r, x_initial, gs_E, t_step, dist, p_pos, NULL, NULL, NULL, sizeof(double[N]), params);

    gsl_rng_free (r);

    printf("best params: ");
    for(int i=0;i<N;i++)
    {
        printf("%.12lf ",x_initial[i]);
    }
    printf("\nwith energy: %.12lf", gs_E(x_initial));

    return 0;

}