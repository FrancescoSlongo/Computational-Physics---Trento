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

double omega;
double orbs(double x, double y, int sel){ // Calculates the orbitals in (x, y)
    if (sel == 1) // 1s
        return exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2) // 1px
        return x*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3) // 1py
        return y*exp(-0.5*omega*(x*x+y*y));
}

double DXorbs(double x, double y, int sel){ // Calculates the x component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -x*omega*exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2)
        return (1-omega*x*x)*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3)
        return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
}

double DYorbs(double x, double y, int sel){ // Calculates the y component of the gradient of the orbitals in (x, y)
    if (sel == 1)
        return -y*omega*exp(-0.5*omega*(x*x+y*y));
    else if (sel == 2)
        return -omega*x*y*exp(-0.5*omega*(x*x+y*y));
    else// if (sel == 3)
        return (1-omega*y*y)*exp(-0.5*omega*(x*x+y*y));
}

double DDorbs(double x, double y, int sel){ // Calculates the laplacian of the orbitals in (x, y)
    if (sel == 1)
        return omega*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-2);
    else if (sel == 2)
        return omega*x*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
    else// if (sel == 3)
        return omega*y*exp(-0.5*omega*(x*x+y*y))*(omega*x*x+omega*y*y-4);
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
gsl_matrix *AluU1; // Stores LU decomposition UP SPIN
gsl_matrix *AluU2;

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
gsl_matrix *AluD1;
gsl_matrix *AluD2;

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
    gsl_matrix_memcpy(AluU1, AijU1);
    gsl_linalg_LU_decomp(AluU1, permU1, &sgnU1);
    gsl_matrix_memcpy(AluU2, AijU2);
    gsl_linalg_LU_decomp(AluU2, permU2, &sgnU2);

    gsl_matrix_memcpy(AluD1, AijD1);
    gsl_linalg_LU_decomp(AluD1, permD1, &sgnD1);
    gsl_matrix_memcpy(AluD2, AijD2);
    gsl_linalg_LU_decomp(AluD2, permD2, &sgnD2);
}

double mfDetU1, mfDetU2, mfDetD1, mfDetD2; // Save determinants as global cause needed in multiple functions
double mfLapU, mfLapD; // Laplacians of mean field term
double mfGrad[10][2]; // Gradient of mean field
double jtwGrad[10][2]; // Gradient of Jastrow
double *stKinRW, *jfKinRW; // Arrays for reweighted kinetic energies
int pointsRW; // Points to reweight (considering starting point too)
void evalKinRW(double **params) // Evaluate kinetic energies (SEE PDF)
{
    double ratio = mfDetU2*mfDetD2/mfDetU1/mfDetD1;
    // Mean field part
    gsl_linalg_LU_invert(AluU1, permU1, AinvU1);
    gsl_linalg_LU_invert(AluD1, permD1, AinvD1);
    gsl_linalg_LU_invert(AluU2, permU2, AinvU2);
    gsl_linalg_LU_invert(AluD2, permD2, AinvD2);

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

    double mfUD2 = 0;
    for (int l = 0; l < szU+szD; l++)
        mfUD2 += mfGrad[l][0]*mfGrad[l][0] + mfGrad[l][1]*mfGrad[l][1];

    // Jastrow term
    double jtwLap1, jtwLap2;
    double dx, dy, dr;
    double uprime; // U'(r)/r
    double jtwMfUD;
    double buu, bud;
    for (int t = 0; t < pointsRW; t++)
    {
        buu = params[t][0];
        bud = params[t][1];
        jtwMfUD = 0;
        jtwLap1 = 0;
        jtwLap2 = 0;
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
                if ((i < szU && l < szU) || (i >= szU && l >= szU))
                {
                    uprime = ajtw[l][i]/(1+buu*dr)/(1+buu*dr)/dr;
                    jtwLap1 += 2*uprime;
                    jtwLap1 -= 4*ajtw[l][i]*buu/(1+buu*dr)/(1+buu*dr)/(1+buu*dr);
                }
                else
                {
                    uprime = ajtw[l][i]/(1+bud*dr)/(1+bud*dr)/dr;
                    jtwLap1 += 2*uprime;
                    jtwLap1 -= 4*ajtw[l][i]*bud/(1+bud*dr)/(1+bud*dr)/(1+bud*dr);
                }
                jtwGrad[l][0] += uprime*dx;
                jtwGrad[l][1] += uprime*dy;
                jtwGrad[i][0] -= uprime*dx;
                jtwGrad[i][1] -= uprime*dy;
            }
        }
        for (int l = 0; l < szU+szD; l++)
            jtwLap2 += jtwGrad[l][0]*jtwGrad[l][0] + jtwGrad[l][1]*jtwGrad[l][1];

        for (int l = 0; l < szU+szD; l++)
            jtwMfUD += jtwGrad[l][0]*mfGrad[l][0] + jtwGrad[l][1]*mfGrad[l][1];

        stKinRW[t] = -.5*(mfLapU + mfLapD + jtwLap1 + jtwLap2 + 2*jtwMfUD);
        jfKinRW[t] = -.25*(mfLapU + mfLapD + jtwLap1 - mfUD2);
    }
}

double pot;
void evalPot(void) // Calculate potential
{
    pot = 0;
    for (int i = 0; i < szU+szD; i++)
    {
        pot += .5*omega*omega*(pposUD[i][0]*pposUD[i][0] + pposUD[i][1]*pposUD[i][1]);
        for (int j = 0; j < i; j++)
            pot += 1/sqrt((pposUD[i][0]-pposUD[j][0])*(pposUD[i][0]-pposUD[j][0]) + (pposUD[i][1]-pposUD[j][1])*(pposUD[i][1]-pposUD[j][1]));
    }
}

long double MTrand(int seed) // Uniformly distributed numbers using MT19937
{
    static mt19937 twister(rand());
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

double evalJtwRW(double **pos, double buu, double bud)
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
double** newpos; // Proposed positions
double *rwEs, *rwSigmas;
void MCprocRW(int Ntherm, int Niter, int Nvoid, double delta, double **params, bool verbose)
{
    // params is a 5x2 matrix with variations on rows and buu, bud on column (0 is reference)
    double newDetU1, newDetD1, newDetU2, newDetD2; // Determinants for proposed positions
    double Pacc, eta;
    double* stKinCum = new double[pointsRW]; // Assuming 0 values, + h1, -h1, +h2, -h2
    double* jfKinCum = new double[pointsRW];
    double* potCum = new double[pointsRW];
    double* Eavg = new double[pointsRW];
    double* Ecum2 = new double[pointsRW];
    double* sigma = new double[pointsRW];
    double cWeight;
    double* weights = new double[pointsRW - 1];
    double jtwTerm, newJtwTerm;
    double rwJtwTerm;

    LUdecomp();
    mfDetU1 = gsl_linalg_LU_det(AluU1, sgnU1);
    mfDetD1 = gsl_linalg_LU_det(AluD1, sgnD1);
    mfDetU2 = gsl_linalg_LU_det(AluU2, sgnU2);
    mfDetD2 = gsl_linalg_LU_det(AluD2, sgnD2);
    evalKinRW(params);
    evalPot();
    jtwTerm = evalJtwRW(pposUD, params[0][0], params[0][1]);

    /// Thermalization steps
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

        gsl_matrix_memcpy(AluU1, AnewU1);
        gsl_matrix_memcpy(AluD1, AnewD1);
        gsl_matrix_memcpy(AluU2, AnewU2);
        gsl_matrix_memcpy(AluD2, AnewD2);
        gsl_linalg_LU_decomp(AluU1, newPermU1, &newSgnU1);
        gsl_linalg_LU_decomp(AluD1, newPermD1, &newSgnD1);
        gsl_linalg_LU_decomp(AluU2, newPermU2, &newSgnU2);
        gsl_linalg_LU_decomp(AluD2, newPermD2, &newSgnD2);
        newDetU1 = gsl_linalg_LU_det(AluU1, newSgnU1);
        newDetD1 = gsl_linalg_LU_det(AluD1, newSgnD1);
        newDetU2 = gsl_linalg_LU_det(AluU2, newSgnU2);
        newDetD2 = gsl_linalg_LU_det(AluD2, newSgnD2);
        newJtwTerm = evalJtwRW(newpos, params[0][0], params[0][1]);

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
            evalKinRW(params);
            evalPot();
        }
    }

    stKinCum[0] = stKinRW[0];
    jfKinCum[0] = jfKinRW[0];
    potCum[0] = pot;
    Ecum2[0] = (stKinRW[0] + pot)*(stKinRW[0] + pot);
    for (int s = 1; s < pointsRW; s++) // Loop over reweighting parameters
    {
        rwJtwTerm = evalJtwRW(pposUD, params[s][0], params[s][1]);
        cWeight = rwJtwTerm*rwJtwTerm/jtwTerm/jtwTerm;
        stKinCum[s] = stKinRW[s]*cWeight;
        jfKinCum[s] = jfKinRW[s]*cWeight;
        potCum[s] = pot*cWeight;
        Ecum2[s] = (stKinRW[s] + pot)*(stKinRW[s] + pot)*cWeight;
        weights[s-1] = cWeight;
    }

    int acc = 0;

    /// True loop
    for (int k = 0; k < Niter; k++)
    {
        for (int z = 0; z < Nvoid; z++)
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
            newJtwTerm = evalJtwRW(newpos, params[0][0], params[0][1]);

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
                gsl_matrix_memcpy(AluU1, AtempU1);
                gsl_matrix_memcpy(AluU2, AtempU2);
                gsl_matrix_memcpy(AluD1, AtempD1);
                gsl_matrix_memcpy(AluD2, AtempD2);
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
                evalKinRW(params);
                evalPot();

                acc++;
            }
        }

        stKinCum[0] += stKinRW[0];
        jfKinCum[0] += jfKinRW[0];
        potCum[0] += pot;
        Ecum2[0] += (stKinRW[0] + pot)*(stKinRW[0] + pot);
        for (int s = 1; s < pointsRW; s++) // Loop over reweighting parameters
        {
            rwJtwTerm = evalJtwRW(pposUD, params[s][0], params[s][1]);
            cWeight = rwJtwTerm*rwJtwTerm/jtwTerm/jtwTerm;
            stKinCum[s] += stKinRW[s]*cWeight;
            jfKinCum[s] += jfKinRW[s]*cWeight;
            potCum[s] += pot*cWeight;
            Ecum2[s] += (stKinRW[s] + pot)*(stKinRW[s] + pot)*cWeight;
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

    for (int s = 1; s < pointsRW; s++)
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
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    omega = 1;

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
    //orbSelectorUD1[4] = 2;

    orbSelectorUD2[0] = 1;
    orbSelectorUD2[1] = 3;
    orbSelectorUD2[2] = 1;
    orbSelectorUD2[3] = 2;
    //orbSelectorUD1[4] = 3;

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
    AluU1 = gsl_matrix_alloc(szU, szU);
    AluU2 = gsl_matrix_alloc(szU, szU);

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
    AluD1 = gsl_matrix_alloc(szD, szD);
    AluD2 = gsl_matrix_alloc(szD, szD);

    initMatsUD(szU, szD);

    /// Steepest descent
    pointsRW = 5;
    stKinRW = new double[5];
    jfKinRW = new double[5];
    rwEs = new double[pointsRW];
    rwSigmas = new double[pointsRW];
    double buu = .5, bud = .51;
    double buuStep = .001, budStep = .001;
    double **bvs;
    bvs = new double*[5];
    for (int i = 0; i < 5; i++)
        bvs[i] = new double[2];
    for (int i = 0; i < 5; i++)
    {
        bvs[i][0] = buu;
        bvs[i][1] = bud;
    }
    bvs[1][0] -= buuStep;
    bvs[2][0] += buuStep;
    bvs[3][1] -= budStep;
    bvs[4][1] += budStep;

    double grad[2] = {10, 10}, gradCoeff = .05;
    double gradThres = 1e-11;
    double gradAbs = 10;

    initJtw(1./3., 1, buu, bud);
    while (gradAbs > gradThres)
    {
        MCprocRW(10, 100000, 1, 1.5, bvs, 0);
        grad[0] = (rwEs[2]-rwEs[1])/2/buuStep;
        grad[1] = (rwEs[4]-rwEs[3])/2/budStep;
        gradAbs = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
        cout << setprecision(6) << "buu = " << bvs[0][0] << ", bud = " << bvs[0][1] << ": grad = (" << grad[0] << ", " << grad[1] << ") (abs: " << gradAbs << ")  E = " << rwEs[0] << " pm " << rwSigmas[0] << endl;
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
    cout << "Final results: buu = " << buu << ", bud = " << bud << endl;
}
