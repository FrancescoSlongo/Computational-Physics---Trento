#include<iostream>
#include<cmath>
#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

using namespace std;

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
    //define matrix S
    gsl_matrix * S = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(S);
    gsl_matrix_set(S, 1, 2, .5);
    gsl_matrix_set(S, 2, 1, .5);

    //diagonalize S
    gsl_vector * val = gsl_vector_alloc(3);
    gsl_matrix * U = gsl_matrix_alloc(3,3);

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(3);

    gsl_matrix * V = gsl_matrix_alloc(3,3);
    gsl_matrix * tmp = gsl_matrix_calloc(3,3);

    print_gsl_matrix(S,"S");

    gsl_eigen_symmv(S, val, U, w);
    print_gsl_vector( val, "eigenvalues");
    print_gsl_matrix(U, "eigenvectors");
    //print_gsl_matrix(S,"S");
    gsl_eigen_symmv_sort(val, U, GSL_EIGEN_SORT_VAL_ASC);

    //Calculate s^-1/2 and then V
    for(int i=0;i<(*U).size1;i++)
    {
        gsl_matrix_set(tmp, i, i, 1./sqrt(gsl_vector_get(val,i)));
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0, U, tmp, 0.0, V);
    print_gsl_matrix(V, "V");

    //define matrix H
    gsl_matrix * H = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(H);
    gsl_matrix_set(H, 0, 1, .5);
    gsl_matrix_set(H, 1, 0, .5);

    gsl_matrix * Hp = gsl_matrix_alloc(3,3);
    gsl_matrix_set_zero(tmp);

    //find the tranformed H'
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., H, V, 0., tmp);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., V, tmp, 0., Hp);

    print_gsl_matrix(Hp, "Hp");

    //diagonlize H'
    gsl_matrix * Cp = gsl_matrix_alloc(3,3);
    gsl_vector * E = gsl_vector_alloc(3);

    gsl_eigen_symmv(Hp, E, Cp, w);
    print_gsl_vector( E, "eigenvalues");
    print_gsl_matrix(Cp, "eigenvectors");

    //find the eigenvectors of the original H
    gsl_matrix * C = gsl_matrix_alloc(3,3);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Cp, 0., C);
    print_gsl_matrix(C, "C");


    //compare with the builtin solver
    gsl_matrix * X = gsl_matrix_alloc(3,3);
    gsl_vector * Y = gsl_vector_alloc(3);
    gsl_eigen_gensymmv_workspace * ww = gsl_eigen_gensymmv_alloc(3);

    gsl_matrix_set_identity(S);
    gsl_matrix_set(S, 1, 2, .5);
    gsl_matrix_set(S, 2, 1, .5);

    gsl_matrix_set_identity(H);
    gsl_matrix_set(H, 0, 1, .5);
    gsl_matrix_set(H, 1, 0, .5);

    gsl_eigen_gensymmv(H, S, Y, X, ww);
    print_gsl_vector(Y, "gsl eigenvalues");
    print_gsl_matrix(X, "gsl eigenvectors");
}
