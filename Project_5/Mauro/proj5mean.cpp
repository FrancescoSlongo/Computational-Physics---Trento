#include <iostream>
#include "mpi.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

using namespace std;

double orbs(double x, double y, double omega, int sel){ // Calculates the orbitals in (x, y)
    if (sel == 1) // 1s
        return exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2) // 1px
        return x*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3) // 1py
        return y*exp(-0.5*omega*(x*x+y*y));
}

double DXorbs(double x, double y, double omega, int sel){ // Calculates the x component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -x*omega*exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2)
        return (1-omega*x*x)*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3)
        return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
}

double DYorbs(double x, double y, double omega, int sel){ // Calculates the y component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -y*omega*exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2)
        return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3)
        return (1-omega*y*y)*exp(-0.5*omega*(x*x+y*y));
}

double DDorbs(double x, double y, double omega, int sel){ // Calculates the laplacian of the orbitals in (x, y)
    if (sel == 1)
        return omega*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-2);
    else if (sel == 2)
        return omega*x*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
    else// if (sel == 3)
        return omega*y*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
}

/// This fcns calculate the single orbitals, components of the gradients and laplacians
double wfS(double x, double y, double omega){
    return exp(-0.5*omega*(x*x+y*y));
}
double DXwfS(double x, double y, double omega){
    return -x*omega*exp(-0.5*omega*(x*x+y*y));
}
double DYwfS(double x, double y, double omega){
    return -y*omega*exp(-0.5*omega*(x*x+y*y));
}
double DDwfS(double x, double y, double omega){
    return omega*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-2);
}
double wfPx(double x, double y, double omega){
    return x*exp(-0.5*omega*(x*x+y*y));
}
double DXwfPx(double x, double y, double omega){
    return (1-omega*x*x)*exp(-0.5*omega*(x*x+y*y));
}
double DYwfPx(double x, double y, double omega){
    return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
}
double DDwfPx(double x, double y, double omega){
    return omega*x*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
}
double wfPy(double x, double y, double omega){
    return y*exp(-0.5*omega*(x*x+y*y));
}
double DXwfPy(double x, double y, double omega){
    return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
}
double DYwfPy(double x, double y, double omega){
    return (1-omega*y*y)*exp(-0.5*omega*(x*x+y*y));
}
double DDwfPy(double x, double y, double omega){
    return omega*y*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
}

double **pposU;
double **pposD;
int *orbSelectorU; // 1->s, 2->px, 3->py
int *orbSelectorD;
gsl_matrix *AijU;
gsl_matrix *DXAijU;
gsl_matrix *DYAijU;
gsl_matrix *DDAijU;
gsl_matrix *AtempU, *AinvU;

gsl_matrix *AijD;
gsl_matrix *DXAijD;
gsl_matrix *DYAijD;
gsl_matrix *DDAijD;
gsl_matrix *AtempD, *AinvD;
double omega;

int szU, szD; // Sizes of Slater det for up and down spins

void initPos(void)
{
    double offset = M_PI/12;
    double step = 2*M_PI/(szU+szD);
    for (int i = 0; i < szU; i++)
    {
        pposU[i][0] = cos(offset + i*step);
        pposU[i][1] = sin(offset + i*step);
        cout << pposU[i][0] << ", " << pposU[i][1] << endl;
    }
    for (int i = szU; i < szU+szD; i++)
    {
        pposD[i-szU][0] = cos(offset + i*step);
        pposD[i-szU][1] = sin(offset + i*step);
        cout << pposD[i-szU][0] << ", " << pposD[i-szU][1] << endl;
    }
}

void initMatsUD(int size1, int size2)
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size1; j++)
        {
            gsl_matrix_set(AijU, j, i, orbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i])); // On row j put orbital i (column) for particle j
            gsl_matrix_set(DXAijU, j, i, DXorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
            gsl_matrix_set(DYAijU, j, i, DYorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
            gsl_matrix_set(DDAijU, j, i, DDorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(AijD, j, i, orbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
            gsl_matrix_set(DXAijD, j, i, DXorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
            gsl_matrix_set(DYAijD, j, i, DYorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
            gsl_matrix_set(DDAijD, j, i, DDorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
        }
    }
    /*cout << gsl_matrix_get(AijD, 0, 0) << endl;
    cout << gsl_matrix_get(DXAijD, 0, 0) << endl;
    cout << gsl_matrix_get(DYAijD, 0, 0) << endl;
    cout << gsl_matrix_get(DDAijD, 0, 0) << endl;*/
}

gsl_permutation *permU;
gsl_permutation *permD;
int sgnU, sgnD;
void LUdecomp(void) // Saves LU decomposition in Atemp, permutation in perm and signum in sgn
{
    // Copy matrices cause decomp changes them
    gsl_matrix_memcpy(AtempU, AijU);
    //cout << gsl_matrix_get(AtempU, 0, 0) << " " << gsl_matrix_get(AtempU, 1, 0) << " " << gsl_matrix_get(AtempU, 0, 1) << " " << gsl_matrix_get(AtempU, 1, 1) << endl;
    gsl_linalg_LU_decomp(AtempU, permU, &sgnU);

    gsl_matrix_memcpy(AtempD, AijD);
    gsl_linalg_LU_decomp(AtempD, permD, &sgnD);
}

double mfGradU[2], mfGradD[2];
double mfLapU, mfLapD;
double stKin, jfKin;
void evalKin(void)
{
    gsl_linalg_LU_invert(AtempU, permU, AinvU);
    gsl_linalg_LU_invert(AtempD, permD, AinvD);

    // Use the LU decomposed matrix
    mfGradU[0] = 0;
    mfGradD[0] = 0;
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DXAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DXAijD, AinvD, 0., AtempD);
    // Traces
    for (int i = 0; i < szU; i++)
        mfGradU[0] += gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGradD[0] += gsl_matrix_get(AtempD, i, i);

    mfGradU[1] = 0;
    mfGradD[1] = 0;
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DYAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DYAijD, AinvD, 0., AtempD);
    for (int i = 0; i < szU; i++)
        mfGradU[1] += gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGradD[1] += gsl_matrix_get(AtempD, i, i);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DDAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., DDAijD, AinvD, 0., AtempD);
    mfLapU = 0;
    mfLapD = 0;
    for (int i = 0; i < szU; i++)
        mfLapU += gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfLapD += gsl_matrix_get(AtempD, i, i);

    //cout << mfGradU[0] << ", " << mfGradU[1] << ", " << mfGradD[0] << ", " << mfGradD[1] << ", " << mfLapD << ", " << mfLapU << endl;
    stKin = -.5*(mfLapU + mfLapD) - mfGradU[0]*mfGradD[0] - mfGradU[1]*mfGradD[1];
    jfKin = -.25*(mfLapU + mfLapD - mfGradU[0]*mfGradU[0] - mfGradU[1]*mfGradU[1] - mfGradD[0]*mfGradD[0] - mfGradD[1]*mfGradD[1]);
}

double detU, detD;
void evalDetUD(void)
{
    detU = gsl_linalg_LU_det(AtempU, sgnU);
    detD = gsl_linalg_LU_det(AtempD, sgnD);
}

double pot;
void evalPot(void)
{
    pot = 0;
    for (int i = 0; i < szU; i++)
    {
        pot += .5*omega*omega*(pposU[i][0]*pposU[i][0] + pposU[i][1]*pposU[i][1]);
    }
    for (int i = 0; i < szD; i++)
    {
        pot += .5*omega*omega*(pposD[i][0]*pposD[i][0] + pposD[i][1]*pposD[i][1]);
    }
}

int main(int argc, char **argv)
{
    szU = 3;
    szD = 2;
    pposU = new double*[szU];
    for (int i = 0; i < szU; i++)
        pposU[i] = new double[2];
    orbSelectorU = new int[szU];
    pposD = new double*[szD];
    for (int i = 0; i < szD; i++)
        pposD[i] = new double[2];
    orbSelectorD = new int[szD];

    pposU[0][0] = 1;
    pposU[0][1] = 0;
    orbSelectorU[0] = 1;
    pposU[1][0] = 0;
    pposU[1][1] = 1;
    orbSelectorU[1] = 2;
    pposU[2][0] = 2;
    pposU[2][1] = 1;
    orbSelectorU[2] = 3;
    pposD[0][0] = 1;
    pposD[0][1] = 1;
    orbSelectorD[0] = 1;
    pposD[1][0] = 0;
    pposD[1][1] = 1;
    orbSelectorD[1] = 2;
    initPos();

    permU = gsl_permutation_alloc(szU);
    AijU = gsl_matrix_alloc(szU, szU);
    DXAijU = gsl_matrix_alloc(szU, szU);
    DYAijU = gsl_matrix_alloc(szU, szU);
    DDAijU = gsl_matrix_alloc(szU, szU);
    AtempU = gsl_matrix_alloc(szU, szU);
    AinvU = gsl_matrix_alloc(szU, szU);

    permD = gsl_permutation_alloc(szD);
    AijD = gsl_matrix_alloc(szD, szD);
    DXAijD= gsl_matrix_alloc(szD, szD);
    DYAijD = gsl_matrix_alloc(szD, szD);
    DDAijD = gsl_matrix_alloc(szD, szD);
    AtempD = gsl_matrix_alloc(szD, szD);
    AinvD = gsl_matrix_alloc(szD, szD);

    omega = 1;

    initMatsUD(szU, szD);

    LUdecomp();
    evalDetUD();
    evalKin();
    evalPot();
    cout << stKin << ", " << jfKin << endl;
    cout << pot << endl;
}
