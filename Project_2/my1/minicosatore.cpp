#include<iostream>
#include<cmath>
#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

using namespace std;

const int N = 3;

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

int main()
{
    gsl_matrix * S = gsl_matrix_alloc(N,N);
    double a[N];
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(N);
    gsl_vector * val = gsl_vector_alloc(N);

    /*for(int i=0;i<N;i++)
    {
        cout<<"alpha i "<<i<<" -> ";
        cin>>a[i];
    }*/

    for(a[0]=0;a[0]<3;a[0]+=0.01)
    {
        for(a[1]=0;a[1]<3;a[1]+=0.01)
        {
            for(a[2]=0;a[2]<3;a[2]+=0.01)
            {
                for(int i=0;i<N;i++)
                {
                    for(int j=0;j<N;j++)
                    {
                        gsl_matrix_set(S, i, j, 0.5*pow(3.14,3./2.)*pow(1/(a[i]+a[j]),5./2.));

                    }

                }
                //print_gsl_matrix(S, "S");

                gsl_eigen_symm(S, val, w);

                //print_gsl_vector(val, "eigenvalues");
                if((gsl_vector_min(val)<0)&&(a[0]!=a[1])&&(a[1]!=a[2])&&(a[0]!=a[2]))
                {
                    cout<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
                    print_gsl_vector(val, "eigenvalues");
                }
            }
        }
    }



}
