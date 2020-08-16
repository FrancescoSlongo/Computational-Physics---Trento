// ASSUMING det1*det2 + det1*det2, CASES OF N = 4, N = 5
#include <iostream>
#include "mpi.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <random>
#include <vector>
#include <fstream>
#include <iomanip>
#include <time.h>

using namespace std;

double orbs(double x, double y, int sel){ // Calculates the orbitals in (x, y)
    if (sel == 1) // 1s
        return exp(-0.5*(x*x+y*y));
    else if (sel == 2) // 1px
        return x*exp(-0.5*(x*x+y*y));
    else// if (sel == 3) // 1py
        return y*exp(-0.5*(x*x+y*y));
}

double DXorbs(double x, double y, int sel){ // Calculates the x component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -x*exp(-0.5*(x*x+y*y));
    else if (sel == 2)
        return (1-x*x)*exp(-0.5*(x*x+y*y));
    else// if (sel == 3)
        return -x*y*exp(-0.5*(x*x+y*y));
}

double DYorbs(double x, double y, int sel){ // Calculates the y component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -y*exp(-0.5*(x*x+y*y));
    else if (sel == 2)
        return -x*y*exp(-0.5*(x*x+y*y));
    else// if (sel == 3)
        return (1-y*y)*exp(-0.5*(x*x+y*y));
}

double DDorbs(double x, double y, int sel){ // Calculates the laplacian of the orbitals in (x, y)
    if (sel == 1)
        return exp(-0.5*(x*x+y*y))*(x*x+y*y-2);
    else if (sel == 2)
        return x*exp(-0.5*(x*x+y*y))*(x*x+y*y-4);
    else// if (sel == 3)
        return y*exp(-0.5*(x*x+y*y))*(x*x+y*y-4);
}

double **pposUD; // Positions of particles
int *orbSelectorUD1; // Associates orbital to particle 1->s, 2->px, 3->py
int *orbSelectorUD2; // Same as above but for second determinant
gsl_matrix *AijU1; // Slater matrix UP SPIN
gsl_matrix *DXAijU1; // X component of gradient of all orbitals in Slater matrix UP SPIN
gsl_matrix *DYAijU1; // Y component UP SPIN
gsl_matrix *DDAijU1; // Laplacian UP SPIN
gsl_matrix *AtempU1, *AinvU1; // Multi-purpose temporary matrix, inverse of Slater matrix UP SPIN
gsl_matrix *AijU2; // Same as above, second determinant
gsl_matrix *DXAijU2;
gsl_matrix *DYAijU2;
gsl_matrix *DDAijU2;
gsl_matrix *AtempU2, *AinvU2;

gsl_matrix *AijD1; // All same as above, but DOWN SPIN
gsl_matrix *DXAijD1;
gsl_matrix *DYAijD1;
gsl_matrix *DDAijD1;
gsl_matrix *AtempD1, *AinvD1;
gsl_matrix *AijD2;
gsl_matrix *DXAijD2;
gsl_matrix *DYAijD2;
gsl_matrix *DDAijD2;
gsl_matrix *AtempD2, *AinvD2;

int szU, szD; // Sizes of Slater matrices for up and down spins

string orbName(int sel) // Just for the next function, kinda useless
{
    if (sel == 1)
        return "1s";
    else if (sel == 2)
        return "1px";
    else// if (sel == 3)
        return "1py";
}

void initPos(double offset, double r, bool verbose) // Initialize particles equidistant on circle of radius r with offset angle wrt 0
{
    double step = 2*M_PI/(szU+szD);
    for (int i = 0; i < szU; i++)
    {
        pposUD[i][0] = r*cos(offset + i*step);
        pposUD[i][1] = r*sin(offset + i*step);
        if (verbose)
        {
            if (orbSelectorUD1[i] != orbSelectorUD2[i])
                cout << "UP   Particle " << i << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbitals " << orbName(orbSelectorUD1[i]) << " and " << orbName(orbSelectorUD2[i]) << endl;
            else
                cout << "UP   Particle " << i << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbital " << orbName(orbSelectorUD1[i]) << endl;
        }
    }
    for (int i = szU; i < szU+szD; i++)
    {
        pposUD[i][0] = r*cos(offset + i*step);
        pposUD[i][1] = r*sin(offset + i*step);
        if (verbose)
        {
            if (orbSelectorUD1[i] != orbSelectorUD2[i])
                cout << "DOWN Particle " << i-szU << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbitals " << orbName(orbSelectorUD1[i]) << " and " << orbName(orbSelectorUD2[i]) << endl;
            else
                cout << "DOWN Particle " << i-szU << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbital " << orbName(orbSelectorUD1[i]) << endl;
        }
    }
}

double** ajtw; // (szU+szD)*(szU+szD) Block matrix with aupup on block diagonals and aupdown off diagonal
double** bjtw; // Same but coefficients b
void initJtw(double auu, double aud, double buu, double bud) // Fill those matrices
{
    for (int i = 0; i < szU; i++)
    {
        for (int j = 0; j < szU; j++)
        {
            ajtw[i][j] = auu;
            bjtw[i][j] = buu;
        }
        for (int j = 0; j < szD; j++)
        {
            ajtw[i][j+szU] = aud;
            bjtw[i][j+szU] = bud;
        }
    }
    for (int i = 0; i < szD; i++)
    {
        for (int j = 0; j < szU; j++)
        {
            ajtw[i+szU][j] = aud;
            bjtw[i+szU][j] = bud;
        }
        for (int j = 0; j < szD; j++)
        {
            ajtw[i+szU][j+szU] = auu;
            bjtw[i+szU][j+szU] = buu;
        }
    }
}

void initMatsUD(int size1, int size2) // Fill Slater and derivatives matrices
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size1; j++)
        {
            gsl_matrix_set(AijU1, j, i, orbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i])); // On row j put orbital i (column) for particle j
            gsl_matrix_set(DXAijU1, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(DYAijU1, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(DDAijU1, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(AijU2, j, i, orbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
            gsl_matrix_set(DXAijU2, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
            gsl_matrix_set(DYAijU2, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
            gsl_matrix_set(DDAijU2, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(AijD1, j, i, orbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DXAijD1, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DYAijD1, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DDAijD1, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(AijD2, j, i, orbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
            gsl_matrix_set(DXAijD2, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
            gsl_matrix_set(DYAijD2, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
            gsl_matrix_set(DDAijD2, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
        }
    }
}

void dersMatsUD(int size1, int size2) // Fill only derivatives matrices (needed during MC)
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size1; j++)
        {
            gsl_matrix_set(DXAijU1, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(DYAijU1, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(DDAijU1, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD1[i]));
            gsl_matrix_set(DXAijU2, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i])); //
            gsl_matrix_set(DYAijU2, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
            gsl_matrix_set(DDAijU2, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD2[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(DXAijD1, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DYAijD1, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DDAijD1, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD1[i+size1]));
            gsl_matrix_set(DXAijD2, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1])); //
            gsl_matrix_set(DYAijD2, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
            gsl_matrix_set(DDAijD2, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD2[i+size1]));
        }
    }
}

gsl_permutation *permU1;
gsl_permutation *permU2;
gsl_permutation *permD1;
gsl_permutation *permD2;
int sgnU1, sgnU2, sgnD1, sgnD2;
void LUdecomp(void) // Saves LU decompositions in AtempU/D, permutations in permU/D and signi in sgnU/D
{
    // Copy matrices cause decomp changes them
    gsl_matrix_memcpy(AtempU1, AijU1);
    gsl_linalg_LU_decomp(AtempU1, permU1, &sgnU1);
    gsl_matrix_memcpy(AtempU2, AijU2);
    gsl_linalg_LU_decomp(AtempU2, permU2, &sgnU2);

    gsl_matrix_memcpy(AtempD1, AijD1);
    gsl_linalg_LU_decomp(AtempD1, permD1, &sgnD1);
    gsl_matrix_memcpy(AtempD2, AijD2);
    gsl_linalg_LU_decomp(AtempD2, permD2, &sgnD2);
}

double mfDetU1, mfDetU2, mfDetD1, mfDetD2; // Save determinants as global cause needed in multiple functions
double mfLapU, mfLapD; // Laplacians of mean field term
double stKin, jfKin; // Final kinetic energies
double mfGrad[10][2]; // Gradient of mean field
double jtwGrad[10][2]; // Gradient of Jastrow
void evalKin(void) // Evaluate kinetic energies (SEE PDF)
{
    // Jastrow term
    double jtwLap1 = 0, jtwLap2 = 0;
    double dx, dy, dr;
    double uprime; // U'(r)/r
    double ratio; // Ratio of determinants (see PDF)
    for (int l = 0; l < szU+szD; l++)
    {
        jtwGrad[l][0] = 0;
        jtwGrad[l][1] = 0;
    }
    for (int l = 0; l < szU+szD; l++)
    {
        for (int i = 0; i < l; i++)
        {
            dx = pposUD[l][0] - pposUD[i][0];
            dy = pposUD[l][1] - pposUD[i][1];
            dr = sqrt(dx*dx+dy*dy);
            uprime = ajtw[l][i]/(1+bjtw[l][i]*dr)/(1+bjtw[l][i]*dr)/dr;
            jtwLap1 += 2*uprime;
            jtwLap1 -= 4*ajtw[l][i]*bjtw[l][i]/(1+bjtw[l][i]*dr)/(1+bjtw[l][i]*dr)/(1+bjtw[l][i]*dr);
            jtwGrad[l][0] += uprime*dx;
            jtwGrad[l][1] += uprime*dy;
            jtwGrad[i][0] -= uprime*dx;
            jtwGrad[i][1] -= uprime*dy;
        }
    }
    for (int l = 0; l < szU+szD; l++)
        jtwLap2 += jtwGrad[l][0]*jtwGrad[l][0] + jtwGrad[l][1]*jtwGrad[l][1];

    ratio = mfDetU2*mfDetD2/mfDetU1/mfDetD1;
    // Mean field part
    gsl_linalg_LU_invert(AtempU1, permU1, AinvU1);
    gsl_linalg_LU_invert(AtempD1, permD1, AinvD1);
    gsl_linalg_LU_invert(AtempU2, permU2, AinvU2);
    gsl_linalg_LU_invert(AtempD2, permD2, AinvD2);

    // Use the LU decomposed matrix to store the product
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijU1, AinvU1, 0., AtempU1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijD1, AinvD1, 0., AtempD1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijU2, AinvU2, 0., AtempU2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijD2, AinvD2, 0., AtempD2);
    // Traces
    for (int i = 0; i < szU; i++)
        mfGrad[i][0] = gsl_matrix_get(AtempU1, i, i)/(1 + ratio) + gsl_matrix_get(AtempU2, i, i)/(1 + 1./ratio);
    for (int i = 0; i < szD; i++)
        mfGrad[i+szU][0] = gsl_matrix_get(AtempD1, i, i)/(1 + ratio) + gsl_matrix_get(AtempD2, i, i)/(1 + 1./ratio);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijU1, AinvU1, 0., AtempU1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijD1, AinvD1, 0., AtempD1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijU2, AinvU2, 0., AtempU2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijD2, AinvD2, 0., AtempD2);
    for (int i = 0; i < szU; i++)
        mfGrad[i][1] = gsl_matrix_get(AtempU1, i, i)/(1 + ratio) + gsl_matrix_get(AtempU2, i, i)/(1 + 1./ratio);
    for (int i = 0; i < szD; i++)
        mfGrad[i+szU][1] = gsl_matrix_get(AtempD1, i, i)/(1 + ratio) + gsl_matrix_get(AtempD2, i, i)/(1 + 1./ratio);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijU1, AinvU1, 0., AtempU1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijD1, AinvD1, 0., AtempD1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijU2, AinvU2, 0., AtempU2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijD2, AinvD2, 0., AtempD2);
    mfLapU = 0;
    mfLapD = 0;
    for (int i = 0; i < szU; i++)
        mfLapU += gsl_matrix_get(AtempU1, i, i)/(1 + ratio) + gsl_matrix_get(AtempU2, i, i)/(1 + 1./ratio);
    for (int i = 0; i < szD; i++)
        mfLapD += gsl_matrix_get(AtempD1, i, i)/(1 + ratio) + gsl_matrix_get(AtempD2, i, i)/(1 + 1./ratio);

    //cout << mfGradU[0] << ", " << mfGradU[1] << ", " << mfGradD[0] << ", " << mfGradD[1] << ", " << mfLapD << ", " << mfLapU << endl;

    double jtwMfUD = 0, mfUD2 = 0;
    for (int l = 0; l < szU+szD; l++)
    {
        jtwMfUD += jtwGrad[l][0]*mfGrad[l][0] + jtwGrad[l][1]*mfGrad[l][1];
        mfUD2 += mfGrad[l][0]*mfGrad[l][0] + mfGrad[l][1]*mfGrad[l][1];
    }
    stKin = -.5*(mfLapU + mfLapD + jtwLap1 + jtwLap2 + 2*jtwMfUD);
    jfKin = -.25*(mfLapU + mfLapD + jtwLap1 - mfUD2);
}

double pot;
void evalPot(void) // Calculate potential
{
    pot = 0;
    for (int i = 0; i < szU+szD; i++)
    {
        pot += .5*(pposUD[i][0]*pposUD[i][0] + pposUD[i][1]*pposUD[i][1]);
        for (int j = 0; j < i; j++)
            pot += 1/sqrt((pposUD[i][0]-pposUD[j][0])*(pposUD[i][0]-pposUD[j][0]) + (pposUD[i][1]-pposUD[j][1])*(pposUD[i][1]-pposUD[j][1]));
    }
}

long double MTrand(int seed) // Uniformly distributed numbers using MT19937
{
    static mt19937 twister(seed);
    static uniform_real_distribution<long double> dist(0.0, 1.0);

    return dist(twister);
}

double evalJtw(double **pos)
{
    double ftemp = 1.;
    double dr;
    for (int i = 0; i < szU+szD; i++)
    {
        for (int j = 0; j < i; j++)
        {
            dr = sqrt((pos[j][0] - pos[i][0])*(pos[j][0] - pos[i][0]) + (pos[j][1] - pos[i][1])*(pos[j][1] - pos[i][1]));
            ftemp *= exp(ajtw[j][i]*dr/(1+bjtw[j][i]*dr));
        }
    }
    return ftemp;
}

double evalJtwBs(double **pos, double buu, double bud)
{
    double ftemp = 1.;
    double dr;
    for (int i = 0; i < szU+szD; i++)
    {
        for (int j = 0; j < i; j++)
        {
            dr = sqrt((pos[j][0] - pos[i][0])*(pos[j][0] - pos[i][0]) + (pos[j][1] - pos[i][1])*(pos[j][1] - pos[i][1]));
            if ((i < szU && j < szU) || (i >= szU && j >= szU))
                ftemp *= exp(ajtw[j][i]*dr/(1+buu*dr));
            else
                ftemp *= exp(ajtw[j][i]*dr/(1+bud*dr));
        }
    }
    return ftemp;
}

gsl_matrix *AnewU1; // Slater matrices for proposed positions in Metropolis
gsl_matrix *AnewD1;
gsl_matrix *AnewU2;
gsl_matrix *AnewD2;
gsl_permutation *newPermU1; // Permutations for proposed positions
gsl_permutation *newPermD1;
gsl_permutation *newPermU2;
gsl_permutation *newPermD2;
int newSgnU1, newSgnD1, newSgnU2, newSgnD2; // Signs for proposed positions
ofstream varEfid;
double** newpos; // Proposed positions
void MCproc(int Ntherm, int Niter, double delta, bool verbose, bool save1, bool save2)
{
    // MC without regeighting. Ntherm thermalization steps, Niter iterations, step delta.
    // Verbose outputs all energies with uncertainties, save1 saves all energies for each iteration, save2 saves total energy and uncertainty for 1 MC run on external file
    double newDetU1, newDetD1, newDetU2, newDetD2; // Determinants for proposed positions
    double Pacc, eta; // Acceptance probability, randomly generated number
    vector <double> stKinArr, jfKinArr, potArr;
    stKinArr.reserve(Niter);
    jfKinArr.reserve(Niter);
    potArr.reserve(Niter);
    double stKinCum, jfKinCum, potCum;
    double stKinCum2, jfKinCum2, potCum2;
    double sigmaSt, sigmaJf, sigmaPot;
    double Eavg, Ecum2, sigma;
    double jtwTerm, newJtwTerm;

    LUdecomp();
    mfDetU1 = gsl_linalg_LU_det(AtempU1, sgnU1);
    mfDetD1 = gsl_linalg_LU_det(AtempD1, sgnD1);
    mfDetU2 = gsl_linalg_LU_det(AtempU2, sgnU2);
    mfDetD2 = gsl_linalg_LU_det(AtempD2, sgnD2);
    evalKin();
    evalPot();
    jtwTerm = evalJtw(pposUD);

    // Thermalization steps
    for (int k = 0; k < Ntherm; k++)
    {
        // UP spin
        for (int i = 0; i < szU; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szU; j++)
            {
                gsl_matrix_set(AnewU1, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD1[j]));
                gsl_matrix_set(AnewU2, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD2[j]));
            }
        }

        // DOWN spin
        for (int i = szU; i < szU+szD; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = szU; j < szU+szD; j++)
            {
                gsl_matrix_set(AnewD1, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD1[j]));
                gsl_matrix_set(AnewD2, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD2[j]));
            }
        }

        gsl_matrix_memcpy(AtempU1, AnewU1);
        gsl_matrix_memcpy(AtempD1, AnewD1);
        gsl_matrix_memcpy(AtempU2, AnewU2);
        gsl_matrix_memcpy(AtempD2, AnewD2);
        gsl_linalg_LU_decomp(AtempU1, newPermU1, &newSgnU1);
        gsl_linalg_LU_decomp(AtempD1, newPermD1, &newSgnD1);
        gsl_linalg_LU_decomp(AtempU2, newPermU2, &newSgnU2);
        gsl_linalg_LU_decomp(AtempD2, newPermD2, &newSgnD2);
        newDetU1 = gsl_linalg_LU_det(AtempU1, newSgnU1);
        newDetD1 = gsl_linalg_LU_det(AtempD1, newSgnD1);
        newDetU2 = gsl_linalg_LU_det(AtempU2, newSgnU2);
        newDetD2 = gsl_linalg_LU_det(AtempD2, newSgnD2);
        newJtwTerm = evalJtw(newpos);

        Pacc = (newDetU1*newDetD1 + newDetU2*newDetD2)*(newDetU1*newDetD1 + newDetU2*newDetD2)/(mfDetU1*mfDetD1 + mfDetU2*mfDetD2)/(mfDetU1*mfDetD1 + mfDetU2*mfDetD2);
        Pacc *= (newJtwTerm*newJtwTerm/jtwTerm/jtwTerm);
        eta = MTrand(0);
        if (eta < Pacc)
        {
            for (int i = 0; i < szU+szD; i++)
            {
                pposUD[i][0] = newpos[i][0];
                pposUD[i][1] = newpos[i][1];
            }
            gsl_matrix_memcpy(AijU1, AnewU1);
            gsl_matrix_memcpy(AijD1, AnewD1);
            gsl_matrix_memcpy(AijU2, AnewU2);
            gsl_matrix_memcpy(AijD2, AnewD2);
            gsl_permutation_memcpy(permU1, newPermU1);
            gsl_permutation_memcpy(permD1, newPermD1);
            gsl_permutation_memcpy(permU2, newPermU2);
            gsl_permutation_memcpy(permD2, newPermD2);
            dersMatsUD(szU, szD);
            mfDetU1 = newDetU1;
            mfDetD1 = newDetD1;
            mfDetU2 = newDetU2;
            mfDetD2 = newDetD2;
            jtwTerm = newJtwTerm;
            evalKin();
            evalPot();
        }
    }

    stKinArr.push_back(stKin);
    jfKinArr.push_back(jfKin);
    potArr.push_back(pot);
    stKinCum = stKin;
    jfKinCum = jfKin;
    potCum = pot;
    stKinCum2 = stKin*stKin;
    jfKinCum2 = jfKin*jfKin;
    potCum2 = pot*pot;
    Ecum2 = (stKin + pot)*(stKin + pot);
    int acc = 0;

    for (int k = 0; k < Niter; k++)
    {
        // UP spin
        for (int i = 0; i < szU; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szU; j++)
            {
                gsl_matrix_set(AnewU1, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD1[j]));
                gsl_matrix_set(AnewU2, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD2[j]));
            }
        }

        // DOWN spin
        for (int i = szU; i < szU+szD; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = szU; j < szU+szD; j++)
            {
                gsl_matrix_set(AnewD1, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD1[j]));
                gsl_matrix_set(AnewD2, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD2[j]));
            }
        }

        gsl_matrix_memcpy(AtempU1, AnewU1);
        gsl_matrix_memcpy(AtempD1, AnewD1);
        gsl_matrix_memcpy(AtempU2, AnewU2);
        gsl_matrix_memcpy(AtempD2, AnewD2);
        gsl_linalg_LU_decomp(AtempU1, newPermU1, &newSgnU1);
        gsl_linalg_LU_decomp(AtempD1, newPermD1, &newSgnD1);
        gsl_linalg_LU_decomp(AtempU2, newPermU2, &newSgnU2);
        gsl_linalg_LU_decomp(AtempD2, newPermD2, &newSgnD2);
        newDetU1 = gsl_linalg_LU_det(AtempU1, newSgnU1);
        newDetD1 = gsl_linalg_LU_det(AtempD1, newSgnD1);
        newDetU2 = gsl_linalg_LU_det(AtempU2, newSgnU2);
        newDetD2 = gsl_linalg_LU_det(AtempD2, newSgnD2);
        newJtwTerm = evalJtw(newpos);

        Pacc = (newDetU1*newDetD1 + newDetU2*newDetD2)*(newDetU1*newDetD1 + newDetU2*newDetD2)/(mfDetU1*mfDetD1 + mfDetU2*mfDetD2)/(mfDetU1*mfDetD1 + mfDetU2*mfDetD2);
        Pacc *= (newJtwTerm*newJtwTerm/jtwTerm/jtwTerm);
        eta = MTrand(0);
        if (eta < Pacc)
        {
            for (int i = 0; i < szU+szD; i++)
            {
                pposUD[i][0] = newpos[i][0];
                pposUD[i][1] = newpos[i][1];
            }
            gsl_matrix_memcpy(AijU1, AnewU1);
            gsl_matrix_memcpy(AijD1, AnewD1);
            gsl_matrix_memcpy(AijU2, AnewU2);
            gsl_matrix_memcpy(AijD2, AnewD2);
            gsl_permutation_memcpy(permU1, newPermU1);
            gsl_permutation_memcpy(permD1, newPermD1);
            gsl_permutation_memcpy(permU2, newPermU2);
            gsl_permutation_memcpy(permD2, newPermD2);
            dersMatsUD(szU, szD);
            mfDetU1 = newDetU1;
            mfDetD1 = newDetD1;
            mfDetU2 = newDetU2;
            mfDetD2 = newDetD2;
            jtwTerm = newJtwTerm;
            evalKin();
            evalPot();
            /*if (jfKin > 1000)
            {
                cout << "Ops" << endl;
                for (int i = 0; i < szU; i++)
                {
                    cout << pposU[i][0] << " " << pposU[i][1] << endl;
                }
                for (int i = 0; i < szD; i++)
                {
                    cout << pposD[i][0] << " " << pposD[i][1] << endl;
                }
            }*/

            acc++;
        }
        stKinArr.push_back(stKin);
        jfKinArr.push_back(jfKin);
        potArr.push_back(pot);
        stKinCum += stKin;
        jfKinCum += jfKin;
        potCum += pot;
        stKinCum2 += stKin*stKin;
        jfKinCum2 += jfKin*jfKin;
        potCum2 += pot*pot;
        Ecum2 += (stKin + pot)*(stKin + pot);
    }
    stKinCum /= (Niter+1);
    jfKinCum /= (Niter+1);
    potCum /= (Niter+1);
    stKinCum2 /= (Niter+1);
    jfKinCum2 /= (Niter+1);
    potCum2 /= (Niter+1);
    sigmaSt = sqrt((stKinCum2-stKinCum*stKinCum)/Niter);
    sigmaJf = sqrt((jfKinCum2-jfKinCum*jfKinCum)/Niter);
    sigmaPot = sqrt((potCum2-potCum*potCum)/Niter);
    Eavg = stKinCum + potCum;
    Ecum2 /= (Niter+1);
    sigma = sqrt((Ecum2-Eavg*Eavg)/Niter);
    if (verbose)
    {
        cout << "Accepted steps: " << 100*(double)acc/Niter << "%." << endl;
        cout << "Mean energy: std: " << setprecision(10) << stKinCum << " pm " << sigmaSt << ", KFC: " << jfKinCum << " pm " << sigmaJf << endl;
        cout << "Mean potential: " << setprecision(10) << potCum << " pm " << sigmaPot << endl;
    }
    cout << "Total E: " << setprecision(20) << Eavg << " pm " << sigma << endl;

    if (save1)
    {
        ofstream outEs;
        outEs.open("./Es.txt");
        for (int i = 0; i < Niter; i++)
            outEs << i << " " << stKinArr[i] << " " << jfKinArr[i] << " " << potArr[i] << endl;
        outEs.close();
    }
    stKinArr.clear();
    jfKinArr.clear();
    potArr.clear();

    if (save2)
        varEfid << setprecision(20) << Eavg << " " << sigma << endl;
}

/*double *rwEs, *rwSigmas;
void MCprocRW(int Ntherm, int Niter, double delta, double params[5][2], bool verbose)
{
    // params is a 5x2 matrix with variations on rows and buu, bud on column (0 is reference)
    double newDetU, newDetD;
    double Pacc, eta;
    double* stKinCum = new double[5]; // Assuming 0 values, + h1, -h1, +h2, -h2
    double* jfKinCum = new double[5];
    double* potCum = new double[5];
    double* Eavg = new double[5];
    double* Ecum2 = new double[5];
    double* sigma = new double[5];
    double cWeight;
    double* weights = new double[4];
    double jtwTerm, newJtwTerm;
    double rwJtwTerm;

    LUdecomp();
    double detU = gsl_linalg_LU_det(AtempU, sgnU);
    double detD = gsl_linalg_LU_det(AtempD, sgnD);
    evalKin();
    evalPot();
    jtwTerm = evalJtwBs(pposUD, params[0][0], params[0][1]);

    // Thermalization steps
    for (int k = 0; k < Ntherm; k++)
    {
        // UP spin
        for (int i = 0; i < szU; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szU; j++)
                gsl_matrix_set(AnewU, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD[j]));
        }

        // DOWN spin
        for (int i = szU; i < szU+szD; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = szU; j < szU+szD; j++)
                gsl_matrix_set(AnewD, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD[j]));
        }

        gsl_matrix_memcpy(AtempU, AnewU);
        gsl_matrix_memcpy(AtempD, AnewD);
        gsl_linalg_LU_decomp(AtempU, newPermU, &newSgnU);
        gsl_linalg_LU_decomp(AtempD, newPermD, &newSgnD);
        newDetU = gsl_linalg_LU_det(AtempU, newSgnU);
        newDetD = gsl_linalg_LU_det(AtempD, newSgnD);
        newJtwTerm = evalJtwBs(newpos, params[0][0], params[0][1]);

        Pacc = newDetU*newDetU*newDetD*newDetD/detD/detD/detU/detU;
        Pacc *= (newJtwTerm*newJtwTerm/jtwTerm/jtwTerm);
        eta = MTrand(0);
        if (eta < Pacc)
        {
            for (int i = 0; i < szU+szD; i++)
            {
                pposUD[i][0] = newpos[i][0];
                pposUD[i][1] = newpos[i][1];
            }
            gsl_matrix_memcpy(AijU, AnewU);
            gsl_matrix_memcpy(AijD, AnewD);
            gsl_permutation_memcpy(permU, newPermU);
            gsl_permutation_memcpy(permD, newPermD);
            dersMatsUD(szU, szD);
            detU = newDetU;
            detD = newDetD;
            jtwTerm = newJtwTerm;
            evalKin();
            evalPot();
        }
    }

    stKinCum[0] = stKin;
    jfKinCum[0] = jfKin;
    potCum[0] = pot;
    Ecum2[0] = (stKin + pot)*(stKin + pot);
    for (int s = 1; s < 5; s++) // Loop over reweighting parameters
    {
        rwJtwTerm = evalJtwBs(pposUD, params[s][0], params[s][1]);
        cWeight = rwJtwTerm*rwJtwTerm/jtwTerm/jtwTerm;
        stKinCum[s] = stKin*cWeight;
        jfKinCum[s] = jfKin*cWeight;
        potCum[s] = pot*cWeight;
        Ecum2[s] = (stKin + pot)*(stKin + pot)*cWeight;
        weights[s-1] = cWeight;
    }

    int acc = 0;

    // True loop
    for (int k = 0; k < Niter; k++)
    {
        // UP spin
        for (int i = 0; i < szU; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = 0; j < szU; j++)
                gsl_matrix_set(AnewU, i, j, orbs(newpos[i][0], newpos[i][1], orbSelectorUD[j]));
        }

        // DOWN spin
        for (int i = szU; i < szU+szD; i++)
        {
            newpos[i][0] = pposUD[i][0] + delta*(rand() / (double)RAND_MAX - 0.5);
            newpos[i][1] = pposUD[i][1] + delta*(rand() / (double)RAND_MAX - 0.5);
            for (int j = szU; j < szU+szD; j++)
                gsl_matrix_set(AnewD, i-szU, j-szU, orbs(newpos[i][0], newpos[i][1], orbSelectorUD[j]));
        }

        gsl_matrix_memcpy(AtempU, AnewU);
        gsl_matrix_memcpy(AtempD, AnewD);
        gsl_linalg_LU_decomp(AtempU, newPermU, &newSgnU);
        gsl_linalg_LU_decomp(AtempD, newPermD, &newSgnD);
        newDetU = gsl_linalg_LU_det(AtempU, newSgnU);
        newDetD = gsl_linalg_LU_det(AtempD, newSgnD);
        newJtwTerm = evalJtwBs(newpos, params[0][0], params[0][1]);

        Pacc = newDetU*newDetU*newDetD*newDetD/detD/detD/detU/detU;
        Pacc *= (newJtwTerm*newJtwTerm/jtwTerm/jtwTerm);
        eta = MTrand(0);
        if (eta < Pacc)
        {
            for (int i = 0; i < szU+szD; i++)
            {
                pposUD[i][0] = newpos[i][0];
                pposUD[i][1] = newpos[i][1];
            }
            gsl_matrix_memcpy(AijU, AnewU);
            gsl_matrix_memcpy(AijD, AnewD);
            gsl_permutation_memcpy(permU, newPermU);
            gsl_permutation_memcpy(permD, newPermD);
            dersMatsUD(szU, szD);
            detU = newDetU;
            detD = newDetD;
            jtwTerm = newJtwTerm;
            evalKin();
            evalPot();

            acc++;
        }

        stKinCum[0] += stKin;
        jfKinCum[0] += jfKin;
        potCum[0] += pot;
        Ecum2[0] += (stKin + pot)*(stKin + pot);
        for (int s = 1; s < 5; s++) // Loop over reweighting parameters
        {
            rwJtwTerm = evalJtwBs(pposUD, params[s][0], params[s][1]);
            cWeight = rwJtwTerm*rwJtwTerm/jtwTerm/jtwTerm;
            stKinCum[s] += stKin*cWeight;
            jfKinCum[s] += jfKin*cWeight;
            potCum[s] += pot*cWeight;
            Ecum2[s] += (stKin + pot)*(stKin + pot)*cWeight;
            weights[s-1] += cWeight;
        }
    }

    stKinCum[0] /= (Niter+1);
    jfKinCum[0] /= (Niter+1);
    potCum[0] /= (Niter+1);
    Eavg[0] = stKinCum[0] + potCum[0];
    Ecum2[0] /= (Niter+1);
    sigma[0] = sqrt((Ecum2[0]-Eavg[0]*Eavg[0])/Niter);
    rwEs[0] = Eavg[0];
    rwSigmas[0] = sigma[0];
    if (verbose)
    {
        cout << "Accepted steps: " << 100*(double)acc/Niter << "%." << endl;
        cout << "REFERENCE PARAMS = " << setprecision(5)<< params[0][0] << ", " << params[0][1] << "; total E: " << setprecision(20) << Eavg[0] << " pm " << sigma[0] << endl;
    }

    for (int s = 1; s < 5; s++)
    {
        stKinCum[s] /= weights[s-1];
        jfKinCum[s] /= weights[s-1];
        potCum[s] /= weights[s-1];
        Eavg[s] = stKinCum[s] + potCum[s];
        Ecum2[s] /= weights[s-1];
        sigma[s] = sqrt((Ecum2[s]-Eavg[s]*Eavg[s])/Niter);
        rwEs[s] = Eavg[s];
        rwSigmas[s] = sigma[s];
        if (verbose)
        {
            cout << "PARAMS = " << setprecision(5)<< params[s][0] << ", " << params[s][1] << "; total E: " << setprecision(20) << Eavg[s] << " pm " << sigma[s] << endl;
            cout << "Weights sum = " << weights[s-1] << endl;
        }
    }

    delete[] stKinCum;
    delete[] jfKinCum;
    delete[] potCum;
    delete[] Eavg;
    delete[] Ecum2;
    delete[] sigma;
    delete[] weights;
}*/

int main(int argc, char **argv)
{
    szU = 2;
    szD = 2;
    pposUD = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        pposUD[i] = new double[2];
    orbSelectorUD1 = new int[szU+szD];
    orbSelectorUD2 = new int[szU+szD];
    newpos = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        newpos[i] = new double[2];

    ajtw = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        ajtw[i] = new double[szU+szD];

    bjtw = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        bjtw[i] = new double[szU+szD];

    orbSelectorUD1[0] = 1;
    orbSelectorUD1[1] = 2;
    orbSelectorUD1[2] = 1;
    orbSelectorUD1[3] = 3;
    orbSelectorUD2[0] = 1;
    orbSelectorUD2[1] = 3;
    orbSelectorUD2[2] = 1;
    orbSelectorUD2[3] = 2;

    initPos(M_PI/12, 1., 1);

    permU1 = gsl_permutation_alloc(szU);
    AijU1 = gsl_matrix_alloc(szU, szU);
    DXAijU1 = gsl_matrix_alloc(szU, szU);
    DYAijU1 = gsl_matrix_alloc(szU, szU);
    DDAijU1 = gsl_matrix_alloc(szU, szU);
    AtempU1 = gsl_matrix_alloc(szU, szU);
    AinvU1 = gsl_matrix_alloc(szU, szU);
    permU2 = gsl_permutation_alloc(szU); //
    AijU2 = gsl_matrix_alloc(szU, szU);
    DXAijU2 = gsl_matrix_alloc(szU, szU);
    DYAijU2 = gsl_matrix_alloc(szU, szU);
    DDAijU2 = gsl_matrix_alloc(szU, szU);
    AtempU2 = gsl_matrix_alloc(szU, szU);
    AinvU2 = gsl_matrix_alloc(szU, szU);

    AnewU1 = gsl_matrix_alloc(szU, szU);
    newPermU1 = gsl_permutation_alloc(szU);
    AnewU2 = gsl_matrix_alloc(szU, szU);
    newPermU2 = gsl_permutation_alloc(szU);

    permD1 = gsl_permutation_alloc(szD);
    AijD1 = gsl_matrix_alloc(szD, szD);
    DXAijD1 = gsl_matrix_alloc(szD, szD);
    DYAijD1 = gsl_matrix_alloc(szD, szD);
    DDAijD1 = gsl_matrix_alloc(szD, szD);
    AtempD1 = gsl_matrix_alloc(szD, szD);
    AinvD1 = gsl_matrix_alloc(szD, szD);
    permD2 = gsl_permutation_alloc(szD); //
    AijD2 = gsl_matrix_alloc(szD, szD);
    DXAijD2 = gsl_matrix_alloc(szD, szD);
    DYAijD2 = gsl_matrix_alloc(szD, szD);
    DDAijD2 = gsl_matrix_alloc(szD, szD);
    AtempD2 = gsl_matrix_alloc(szD, szD);
    AinvD2 = gsl_matrix_alloc(szD, szD);

    AnewD1 = gsl_matrix_alloc(szD, szD);
    newPermD1 = gsl_permutation_alloc(szD);
    AnewD2 = gsl_matrix_alloc(szD, szD);
    newPermD2 = gsl_permutation_alloc(szD);

    srand(time(NULL));
    varEfid.open("./varE.txt");
    initMatsUD(szU, szD);
    initJtw(1./3, 1, .3, .4);
    MCproc(100, 1000000, 1, 1, 0, 0); // save1 is for energies at every step, save2 for total energy & uncertainty for multiple runs

    /*double buu = 1, bud = 0;
    double buuStep = .01, budStep = .001;
    rwEs = new double[5];
    rwSigmas = new double[5];
    double bvs[5][2] = {{buu, bud}, {buu-buuStep, bud}, {buu+buuStep, bud}, {buu, bud-budStep}, {buu, bud+budStep}};
    initJtw(1./3., 1, buu, bud);*/

    /*MCprocRW(1000, 100000, 1.5, bvs, 1);
    cout << endl;
    initJtw(2, 2./3., buu, bud);
    MCproc(100, 100000, 1.5, 0, 0, 0);
    initJtw(2, 2./3., buu, bud-budStep);
    MCproc(100, 100000, 1.5, 0, 0, 0);
    initJtw(2, 2./3., buu, bud+budStep);
    MCproc(100, 100000, 1.5, 0, 0, 0);*/

    /*double grad[2] = {10, 10}, gradCoeff = .05;
    double gradThres = 1e-11;
    double gradAbs = 10;
    while (gradAbs > gradThres)
    {
        MCprocRW(0, 100000, 1.5, bvs, 0);
        grad[0] = (rwEs[2]-rwEs[1])/2/buuStep;
        grad[1] = (rwEs[4]-rwEs[3])/2/budStep;
        gradAbs = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
        cout << setprecision(10) << bvs[0][0] << ", " << bvs[0][1] << ": (" << grad[0] << ", " << grad[1] << ") abs: " << gradAbs << "   E = " << rwEs[0] << " pm " << rwSigmas[0] << endl;
        buu -= grad[0]*gradCoeff;
        bud -= grad[1]*gradCoeff;
        initJtw(1./3., 1, buu, bud);
        for (int i = 0; i < 5; i++)
        {
            bvs[i][0] = buu;
            bvs[i][1] = bud;
        }
        bvs[1][0] -= buuStep;
        bvs[2][0] += buuStep;
        bvs[3][1] -= budStep;
        bvs[4][1] += budStep;
        //cout << rwEs[0] << " " << rwEs[1] << " " << rwEs[2] << endl;
    }
    buu += grad[0]*gradCoeff;
    bud += grad[1]*gradCoeff;
    cout << buu << ", " << bud << endl;*/

    /*for (int z = 0; z < 20; z++)
    {
        bud = .1+z*.01;
        initJtw(2, 2./3., 1, bud);
        cout << z << endl;
        varEfid << bud << " ";
        MCproc(0, 100000, 1.5, 1, 0, 1);
    }*/

    /*LUdecomp();
    double detU = gsl_linalg_LU_det(AtempU, sgnU);
    double detD = gsl_linalg_LU_det(AtempD, sgnD);
    evalKin();
    evalPot();
    cout << setprecision(20) << stKin << ", " << jfKin << "; " << pot << endl;*/
    varEfid.close();
}
