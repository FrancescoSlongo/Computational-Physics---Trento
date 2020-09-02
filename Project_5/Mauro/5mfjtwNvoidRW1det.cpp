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
int *orbSelectorUD; // Associates orbital to particle 1->s, 2->px, 3->py
gsl_matrix *AijU; // Slater matrix UP SPIN
gsl_matrix *DXAijU; // X component of gradient of all orbitals in Slater matrix UP SPIN
gsl_matrix *DYAijU; // Y component UP SPIN
gsl_matrix *DDAijU; // Laplacian UP SPIN
gsl_matrix *AtempU, *AinvU; // Multi-purpose temporary matrix, inverse of Slater matrix UP SPIN
gsl_matrix *AluU; // Stores LU decomposition UP SPIN

gsl_matrix *AijD; // All same as above, but DOWN SPIN
gsl_matrix *DXAijD;
gsl_matrix *DYAijD;
gsl_matrix *DDAijD;
gsl_matrix *AtempD, *AinvD;
gsl_matrix *AluD;

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
            cout << "UP   Particle " << i << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbital " << orbName(orbSelectorUD[i]) << endl;
    }
    for (int i = szU; i < szU+szD; i++)
    {
        pposUD[i][0] = r*cos(offset + i*step);
        pposUD[i][1] = r*sin(offset + i*step);
        if (verbose)
            cout << "DOWN Particle " << i-szU << ": (" << pposUD[i][0] << ", " << pposUD[i][1] << ") in orbital " << orbName(orbSelectorUD[i]) << endl;
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
            gsl_matrix_set(AijU, j, i, orbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i])); // On row j put orbital i (column) for particle j
            gsl_matrix_set(DXAijU, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
            gsl_matrix_set(DYAijU, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
            gsl_matrix_set(DDAijU, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(AijD, j, i, orbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
            gsl_matrix_set(DXAijD, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
            gsl_matrix_set(DYAijD, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
            gsl_matrix_set(DDAijD, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
        }
    }
    /*cout << gsl_matrix_get(AijD, 0, 0) << endl;
    cout << gsl_matrix_get(DXAijD, 0, 0) << endl;
    cout << gsl_matrix_get(DYAijD, 0, 0) << endl;
    cout << gsl_matrix_get(DDAijD, 0, 0) << endl;*/
}

void dersMatsUD(int size1, int size2) // Fill only derivatives matrices (needed during MC)
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size1; j++)
        {
            gsl_matrix_set(DXAijU, j, i, DXorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
            gsl_matrix_set(DYAijU, j, i, DYorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
            gsl_matrix_set(DDAijU, j, i, DDorbs(pposUD[j][0], pposUD[j][1], orbSelectorUD[i]));
        }
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size2; j++)
        {
            gsl_matrix_set(DXAijD, j, i, DXorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
            gsl_matrix_set(DYAijD, j, i, DYorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
            gsl_matrix_set(DDAijD, j, i, DDorbs(pposUD[j+size1][0], pposUD[j+size1][1], orbSelectorUD[i+size1]));
        }
    }
}

gsl_permutation *permU;
gsl_permutation *permD;
int sgnU, sgnD;
void LUdecomp(void) // Saves LU decompositions in AtempU/D, permutations in permU/D and signi in sgnU/D
{
    // Copy matrices cause decomp changes them
    gsl_matrix_memcpy(AluU, AijU);
    gsl_linalg_LU_decomp(AluU, permU, &sgnU);

    gsl_matrix_memcpy(AluD, AijD);
    gsl_linalg_LU_decomp(AluD, permD, &sgnD);
}

double mfLapU, mfLapD; // Laplacians of mean field term
double mfGrad[10][2]; // Gradient of mean field
double jtwGrad[10][2]; // Gradient of Jastrow
double *stKinRW, *jfKinRW; // Arrays for reweighted kinetic energies
int pointsRW; // Points to reweight (considering starting point too)
void evalKinRW(double **params) // Evaluate kinetic energies (SEE PDF)
{
    // Mean field part
    gsl_linalg_LU_invert(AluU, permU, AinvU);
    gsl_linalg_LU_invert(AluD, permD, AinvD);

    // Use the LU decomposed matrix to store the product
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DXAijD, AinvD, 0., AtempD);
    // Traces
    for (int i = 0; i < szU; i++)
        mfGrad[i][0] = gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGrad[i+szU][0] = gsl_matrix_get(AtempD, i, i);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DYAijD, AinvD, 0., AtempD);
    for (int i = 0; i < szU; i++)
        mfGrad[i][1] = gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfGrad[i+szU][1] = gsl_matrix_get(AtempD, i, i);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijU, AinvU, 0., AtempU);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., DDAijD, AinvD, 0., AtempD);
    mfLapU = 0;
    mfLapD = 0;
    for (int i = 0; i < szU; i++)
        mfLapU += gsl_matrix_get(AtempU, i, i);
    for (int i = 0; i < szD; i++)
        mfLapD += gsl_matrix_get(AtempD, i, i);

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
    //static mt19937 twister(seed);
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

gsl_matrix *AnewU; // Slater matrices for proposed positions in Metropolis
gsl_matrix *AnewD;
gsl_permutation *newPermU; // Permutations for proposed positions
gsl_permutation *newPermD;
int newSgnU, newSgnD; // Signs for proposed positions
double** newpos; // Proposed positions
double *rwEs, *rwSigmas;
void MCprocRW(int Ntherm, int Niter, int Nvoid, double delta, double **params, bool verbose)
{
    // params is a 5x2 matrix with variations on rows and buu, bud on column (0 is reference)
    double newDetU, newDetD;
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
    double detU = gsl_linalg_LU_det(AluU, sgnU);
    double detD = gsl_linalg_LU_det(AluD, sgnD);
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

        gsl_matrix_memcpy(AluU, AnewU);
        gsl_matrix_memcpy(AluD, AnewD);
        gsl_linalg_LU_decomp(AluU, newPermU, &newSgnU);
        gsl_linalg_LU_decomp(AluD, newPermD, &newSgnD);
        newDetU = gsl_linalg_LU_det(AluU, newSgnU);
        newDetD = gsl_linalg_LU_det(AluD, newSgnD);
        newJtwTerm = evalJtwRW(newpos, params[0][0], params[0][1]);

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
        for (int z = 0; z < Nvoid + 1; z++)
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
            newJtwTerm = evalJtwRW(newpos, params[0][0], params[0][1]);

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
                gsl_matrix_memcpy(AluU, AtempU);
                gsl_matrix_memcpy(AluD, AtempD);
                gsl_permutation_memcpy(permU, newPermU);
                gsl_permutation_memcpy(permD, newPermD);
                dersMatsUD(szU, szD);
                detU = newDetU;
                detD = newDetD;
                jtwTerm = newJtwTerm;

                acc++;
            }
        }
        evalKinRW(params);
        evalPot();

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

    szU = 1;
    szD = 1;
    pposUD = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        pposUD[i] = new double[2];
    orbSelectorUD = new int[szU+szD];
    newpos = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        newpos[i] = new double[2];

    ajtw = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        ajtw[i] = new double[szU+szD];

    bjtw = new double*[szU+szD];
    for (int i = 0; i < szU+szD; i++)
        bjtw[i] = new double[szU+szD];

    orbSelectorUD[0] = 1;
    orbSelectorUD[1] = 1;
    //orbSelectorUD[1] = 2;
    //orbSelectorUD[2] = 3;
    //orbSelectorUD[3] = 1;
    //orbSelectorUD[4] = 2;
    //orbSelectorUD[5] = 3;

    initPos(M_PI/12, 1., 1);

    permU = gsl_permutation_alloc(szU);
    AijU = gsl_matrix_alloc(szU, szU);
    DXAijU = gsl_matrix_alloc(szU, szU);
    DYAijU = gsl_matrix_alloc(szU, szU);
    DDAijU = gsl_matrix_alloc(szU, szU);
    AtempU = gsl_matrix_alloc(szU, szU);
    AinvU = gsl_matrix_alloc(szU, szU);
    AnewU = gsl_matrix_alloc(szU, szU);
    newPermU = gsl_permutation_alloc(szU);
    AluU = gsl_matrix_alloc(szU, szU);

    permD = gsl_permutation_alloc(szD);
    AijD = gsl_matrix_alloc(szD, szD);
    DXAijD = gsl_matrix_alloc(szD, szD);
    DYAijD = gsl_matrix_alloc(szD, szD);
    DDAijD = gsl_matrix_alloc(szD, szD);
    AtempD = gsl_matrix_alloc(szD, szD);
    AinvD = gsl_matrix_alloc(szD, szD);
    AnewD = gsl_matrix_alloc(szD, szD);
    newPermD = gsl_permutation_alloc(szD);
    AluD = gsl_matrix_alloc(szD, szD);

    initMatsUD(szU, szD);

    /// Steepest descent
    pointsRW = 5;
    stKinRW = new double[5];
    jfKinRW = new double[5];
    rwEs = new double[pointsRW];
    rwSigmas = new double[pointsRW];
    double buu = 1, bud = .3853;
    double buuStep = .001, budStep = .0001;
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

    double grad[2] = {10, 10}, gradCoeff = .03;
    double gradThres = 1e-11;
    double gradAbs = 10;

    initJtw(1./3., 1, buu, bud);
    while (gradAbs > gradThres)
    {
        MCprocRW(10, 1e7, 0, 2, bvs, 0);
        grad[0] = (rwEs[2]-rwEs[1])/2/buuStep;
        grad[1] = (rwEs[4]-rwEs[3])/2/budStep;
        gradAbs = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
        cout << setprecision(10) << "buu = " << bvs[0][0] << ", bud = " << bvs[0][1] << ": grad = (" << grad[0] << ", " << grad[1] << ") (abs: " << gradAbs << ")  E = " << rwEs[0] << " pm " << rwSigmas[0] << endl;
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
