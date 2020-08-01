#include <iostream>
#include "mpi.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <random>
#include <vector>
#include <fstream>

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
/*double wfS(double x, double y, double omega){
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
}*/

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

string orbName(int sel)
{
    if (sel == 1)
        return "1s";
    else if (sel == 2)
        return "1px";
    else// if (sel == 3)
        return "1py";
}

void initPos(double offset, double r, bool verbose) // Initialize particles equidistant on circle of radius 1
{
    double step = 2*M_PI/(szU+szD);
    for (int i = 0; i < szU; i++)
    {
        pposU[i][0] = r*cos(offset + i*step);
        pposU[i][1] = r*sin(offset + i*step);
        if (verbose)
            cout << "UP   Particle " << i << ": (" << pposU[i][0] << ", " << pposU[i][1] << ") in orbital " << orbName(orbSelectorU[i]) << endl;
    }
    for (int i = szU; i < szU+szD; i++)
    {
        pposD[i-szU][0] = r*cos(offset + i*step);
        pposD[i-szU][1] = r*sin(offset + i*step);
        if (verbose)
            cout << "DOWN Particle " << i-szU << ": (" << pposD[i-szU][0] << ", " << pposD[i-szU][1] << ") in orbital " << orbName(orbSelectorD[i-szU]) << endl;
    }
}

void initMatsUD(int size1, int size2) // Fill matrices
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

void dersMatsUD(int size1, int size2) // Fill matrices
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size1; j++)
        {
            gsl_matrix_set(DXAijU, j, i, DXorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
            gsl_matrix_set(DYAijU, j, i, DYorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
            gsl_matrix_set(DDAijU, j, i, DDorbs(pposU[j][0], pposU[j][1], omega, orbSelectorU[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(DXAijD, j, i, DXorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
            gsl_matrix_set(DYAijD, j, i, DYorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
            gsl_matrix_set(DDAijD, j, i, DDorbs(pposD[j][0], pposD[j][1], omega, orbSelectorD[i]));
        }
    }
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
void evalKin(void) // Kinetic energies
{
    gsl_linalg_LU_invert(AtempU, permU, AinvU);
    gsl_linalg_LU_invert(AtempD, permD, AinvD);

    // Use the LU decomposed matrix
    mfGradU[0] = 0;
    mfGradD[0] = 0;
    double gradXt, gradYt;
    /*for (int t = 0; t < szU; t++)
    {
        gradXt = 0;
        gradXt = 0;
        for (int i = 0; i < szU; i++)
        {
            gradXt += gsl_matrix_get(DXAijU, t, i)*gsl_matrix_get(AinvU, t, i);
            gradYt += gsl_matrix_get(DYAijU, t, i)*gsl_matrix_get(AinvU, t, i);
        }

        mfGradU[0] += gradXt;
        mfGradU[1] += gradYt;
    }
    for (int t = 0; t < szD; t++)
    {
        gradXt = 0;
        gradXt = 0;
        for (int i = 0; i < szD; i++)
        {
            gradXt += gsl_matrix_get(DXAijD, t, i)*gsl_matrix_get(AinvD, t, i);
            gradYt += gsl_matrix_get(DYAijD, t, i)*gsl_matrix_get(AinvD, t, i);
        }

        mfGradD[0] += gradXt;
        mfGradD[1] += gradYt;
    }

    double lapt;
    for (int t = 0; t < szU; t++)
    {
        lapt = 0;
        for (int i = 0; i < szU; i++)
            lapt += gsl_matrix_get(DDAijU, t, i)*gsl_matrix_get(AinvU, t, i);

        mfLapU += lapt;
    }
    for (int t = 0; t < szD; t++)
    {
        lapt = 0;
        for (int i = 0; i < szU; i++)
            lapt += gsl_matrix_get(DDAijD, t, i)*gsl_matrix_get(AinvD, t, i);

        mfLapD += lapt;
    }*/
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijD, AinvD, 0., AtempD);
    // Traces
    for (int i = 0; i < szU; i++)
        mfGradU[0] += gsl_matrix_get(AtempU, i, i)*gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGradD[0] += gsl_matrix_get(AtempD, i, i)*gsl_matrix_get(AtempD, i, i);

    mfGradU[1] = 0;
    mfGradD[1] = 0;
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijD, AinvD, 0., AtempD);
    for (int i = 0; i < szU; i++)
        mfGradU[1] += gsl_matrix_get(AtempU, i, i)*gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGradD[1] += gsl_matrix_get(AtempD, i, i)*gsl_matrix_get(AtempD, i, i);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijD, AinvD, 0., AtempD);
    mfLapU = 0;
    mfLapD = 0;
    for (int i = 0; i < szU; i++)
        mfLapU += gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfLapD += gsl_matrix_get(AtempD, i, i);

    //cout << mfGradU[0] << ", " << mfGradU[1] << ", " << mfGradD[0] << ", " << mfGradD[1] << ", " << mfLapD << ", " << mfLapU << endl;
    stKin = -.5*(mfLapU + mfLapD);
    jfKin = -.25*(mfLapU + mfLapD - mfGradU[0] - mfGradU[1] - mfGradD[0] - mfGradD[1]);
}

double pot;
void evalPot(void) // Calculate potential
{
    pot = 0;
    for (int i = 0; i < szU; i++)
        pot += .5*(pposU[i][0]*pposU[i][0] + pposU[i][1]*pposU[i][1]);

    for (int i = 0; i < szD; i++)
        pot += .5*(pposD[i][0]*pposD[i][0] + pposD[i][1]*pposD[i][1]);
}

double evalDetUD(gsl_matrix *MAT, int sgnTemp) // Calculate Slater determinants
{
    return gsl_linalg_LU_det(MAT, sgnTemp);
}

long double MTrand(int seed) // Uniformly distributed numbers using MT19937
{
    static mt19937 twister(seed);
    static uniform_real_distribution<long double> dist(0.0, 1.0);

    return dist(twister);
}

gsl_matrix *AnewU;
gsl_matrix *AnewD;
gsl_permutation *newPermU;
gsl_permutation *newPermD;
int newSgnU, newSgnD;
void MCproc(int Ntherm, int Niter, double delta, double omega)
{
    double newpos[10][2];
    double newDetU, newDetD;
    double Pacc, eta;
    vector <double> stKinArr, jfKinArr, potArr;
    stKinArr.reserve(Niter+1);
    jfKinArr.reserve(Niter+1);
    potArr.reserve(Niter+1);
    double stKinCum = 0, jfKinCum = 0, potCum = 0;

    LUdecomp();
    double detU = evalDetUD(AtempU, sgnU);
    double detD = evalDetUD(AtempD, sgnD);
    evalKin();
    evalPot();
    stKinArr.push_back(stKin);
    jfKinArr.push_back(jfKin);
    potArr.push_back(pot);
    stKinCum += stKin;
    jfKinCum += jfKin;
    potCum += pot;
    int acc = 0;

    for (int k = 0; k < Niter; k++)
    {
        gsl_matrix_memcpy(AnewU, AijU);
        // UP spin
        for (int i = 0; i < szU; i++)
        {
            newpos[i][0] = pposU[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposU[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szU; j++)
                gsl_matrix_set(AnewU, i, j, orbs(newpos[i][0], newpos[i][1], omega, orbSelectorU[j]));
        }

        gsl_matrix_memcpy(AnewD, AijD);
        // DOWN spin
        for (int i = 0; i < szD; i++)
        {
            newpos[i+szU][0] = pposD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i+szU][1] = pposD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szD; j++)
                gsl_matrix_set(AnewD, i, j, orbs(newpos[i+szU][0], newpos[i+szU][1], omega, orbSelectorD[j]));
        }

        gsl_matrix_memcpy(AtempU, AnewU);
        gsl_matrix_memcpy(AtempD, AnewD);
        gsl_linalg_LU_decomp(AtempU, newPermU, &newSgnU);
        gsl_linalg_LU_decomp(AtempD, newPermD, &newSgnD);
        newDetU = evalDetUD(AtempU, newSgnU);
        newDetD = evalDetUD(AtempD, newSgnD);

        Pacc = newDetU*newDetU*newDetD*newDetD/detD/detD/detU/detU;
        eta = MTrand(0);
        if (eta < Pacc)
        {
            for (int i = 0; i < szU; i++)
            {
                pposU[i][0] = newpos[i][0];
                pposU[i][1] = newpos[i][1];
            }
            for (int i = 0; i < szD; i++)
            {
                pposD[i][0] = newpos[i+szU][0];
                pposD[i][1] = newpos[i+szU][1];
            }
            gsl_matrix_memcpy(AijU, AnewU);
            gsl_matrix_memcpy(AijD, AnewD);
            gsl_permutation_memcpy(permU, newPermU);
            gsl_permutation_memcpy(permD, newPermD);
            dersMatsUD(szU, szD);
            detU = newDetU;
            detD = newDetD;
            evalKin();
            evalPot();

            acc++;
        }
        stKinArr.push_back(stKin);
        jfKinArr.push_back(jfKin);
        potArr.push_back(pot);
        stKinCum += stKin;
        jfKinCum += jfKin;
        potCum += pot;
    }
    cout << "Accepted steps: " << 100*(double)acc/Niter << "%." << endl;
    cout << "Mean energy: std " << stKinCum/Niter << ", KFC: " << jfKinCum/Niter << endl;
    cout << "Mean potential: " << potCum/Niter << endl;
    cout << "Total E: " << stKinCum/Niter + potCum/Niter << endl;

    /*ofstream outEs;
    outEs.open("./Es.txt");
    for (int i = 0; i < Niter; i++)
        outEs << i << " " << stKinArr[i] << " " << jfKinArr[i] << " " << potArr[i] << endl;*/
}

int main(int argc, char **argv)
{
    szU = 1;
    szD = 1;
    pposU = new double*[szU];
    for (int i = 0; i < szU; i++)
        pposU[i] = new double[2];
    orbSelectorU = new int[szU];
    pposD = new double*[szD];
    for (int i = 0; i < szD; i++)
        pposD[i] = new double[2];
    orbSelectorD = new int[szD];

    //pposU[0][0] = 1;
    //pposU[0][1] = 0;
    orbSelectorU[0] = 1;
    //orbSelectorU[1] = 2;
    //orbSelectorU[2] = 3;
    //pposU[1][0] = 2;
    //pposU[1][1] = -1;
    //pposU[2][0] = 2;
    //pposU[2][1] = 1;
    //pposD[0][0] = 0;
    //pposD[0][1] = 2;
    orbSelectorD[0] = 1;
    //orbSelectorD[1] = 2;
    //orbSelectorD[2] = 3;
    //pposD[1][0] = -2;
    //pposD[1][1] = 3;
    initPos(0, 1., 1);

    permU = gsl_permutation_alloc(szU);
    AijU = gsl_matrix_alloc(szU, szU);
    DXAijU = gsl_matrix_alloc(szU, szU);
    DYAijU = gsl_matrix_alloc(szU, szU);
    DDAijU = gsl_matrix_alloc(szU, szU);
    AtempU = gsl_matrix_alloc(szU, szU);
    AinvU = gsl_matrix_alloc(szU, szU);

    AnewU = gsl_matrix_alloc(szU, szU);
    newPermU = gsl_permutation_alloc(szU);

    permD = gsl_permutation_alloc(szD);
    AijD = gsl_matrix_alloc(szD, szD);
    DXAijD= gsl_matrix_alloc(szD, szD);
    DYAijD = gsl_matrix_alloc(szD, szD);
    DDAijD = gsl_matrix_alloc(szD, szD);
    AtempD = gsl_matrix_alloc(szD, szD);
    AinvD = gsl_matrix_alloc(szD, szD);

    AnewD = gsl_matrix_alloc(szD, szD);
    newPermD = gsl_permutation_alloc(szD);

    omega = 1;

    initMatsUD(szU, szD);

    srand(53);
    MCproc(0, 1000000, 1, 1);
    //LUdecomp();
    //cout << evalDetUD(AtempU, sgnU)*evalDetUD(AtempD, sgnD) << endl;
    //evalKin();
    //cout << stKin << ", " << jfKin << endl;
}
