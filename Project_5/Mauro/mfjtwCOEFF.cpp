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
int *orbSelectorUD; // Associates orbital to particle 1->s, 2->px, 3->py
gsl_matrix *AijU; // Slater matrix UP SPIN
gsl_matrix *DXAijU; // X component of gradient of all orbitals in Slater matrix UP SPIN
gsl_matrix *DYAijU; // Y component UP SPIN
gsl_matrix *DDAijU; // Laplacian UP SPIN
gsl_matrix *AtempU, *AinvU; // Multi-purpose temporary matrix, inverse of Slater matrix UP SPIN

gsl_matrix *AijD; // All same as above, but DOWN SPIN
gsl_matrix *DXAijD;
gsl_matrix *DYAijD;
gsl_matrix *DDAijD;
gsl_matrix *AtempD, *AinvD;

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
    gsl_matrix_memcpy(AtempU, AijU);
    //cout << gsl_matrix_get(AtempU, 0, 0) << " " << gsl_matrix_get(AtempU, 1, 0) << " " << gsl_matrix_get(AtempU, 0, 1) << " " << gsl_matrix_get(AtempU, 1, 1) << endl;
    gsl_linalg_LU_decomp(AtempU, permU, &sgnU);

    gsl_matrix_memcpy(AtempD, AijD);
    gsl_linalg_LU_decomp(AtempD, permD, &sgnD);
}

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

    // Mean field part
    gsl_linalg_LU_invert(AtempU, permU, AinvU);
    gsl_linalg_LU_invert(AtempD, permD, AinvD);

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

double evalJtw(double **pos) // Calculate jastrow factor given positions
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

double evalJtwBs(double **pos, double buu, double bud) // Calculate jastrow given positions and bupup and bupdown
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
ofstream varEfid;
double** newpos; // Proposed positions

int main()
{
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

    permD = gsl_permutation_alloc(szD);
    AijD = gsl_matrix_alloc(szD, szD);
    DXAijD= gsl_matrix_alloc(szD, szD);
    DYAijD = gsl_matrix_alloc(szD, szD);
    DDAijD = gsl_matrix_alloc(szD, szD);
    AtempD = gsl_matrix_alloc(szD, szD);
    AinvD = gsl_matrix_alloc(szD, szD);

    AnewD = gsl_matrix_alloc(szD, szD);
    newPermD = gsl_permutation_alloc(szD);

    srand(time(NULL));
    initMatsUD(szU, szD);
    initJtw(1./3., 1, .4, .5);
    //initJtw(0, 0, .4, .5);

    int Niter = 1e5; // MC steps for check
    int Ntau = 100; // Maximum possible Nvoid to check
    int Ntherm = 100; // Thermalization steps
    double delta = 1.8; // Delta to check
	ofstream fidC, fidV;
	fidC.open("./coeffs.txt"); // Save coefficients c(k)
	fidV.open("./pots.txt"); // Save potentials for all steps (not needed)
    double newDetU, newDetD; // Determinants for proposed positions
    double Pacc, eta; // Acceptance probability, randomly generated number
    vector <double> stKinArr, jfKinArr, potArr, coeffs; // Vectors of energies and c(k)
    stKinArr.reserve(Niter);
    jfKinArr.reserve(Niter);
    potArr.reserve(Niter);
    coeffs.reserve(Ntau);
    double stKinCum, jfKinCum, potCum;
    double stKinCum2, jfKinCum2, potCum2;
    double sigmaSt, sigmaJf, sigmaPot;
    double Eavg, Ecum2, sigma;
    double jtwTerm, newJtwTerm;

    LUdecomp();
    double detU = gsl_linalg_LU_det(AtempU, sgnU); // Just the beginning of the standard MC procedure
    double detD = gsl_linalg_LU_det(AtempD, sgnD);
    evalKin();
    evalPot();
    jtwTerm = evalJtw(pposUD);
    int lps, acc;

    for (int z = 0; z < 10; z++) // Consider multiple delta, save coefficients as "k c(k)", Ntau rows repeated for these 10 values of delta
    {
        delta = 1.5 + z*.1;
        //initJtw(1./3., 1, .1+z*.05, .5); // Can also check with fixed delta and b varying
        lps = 0;
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
            newJtwTerm = evalJtw(newpos);

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

            fidV << lps << " " << stKin + pot << endl; // Can comment if you don't want to save potentials
            lps++;
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

        acc = 0;

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
            newJtwTerm = evalJtw(newpos);

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
            fidV << lps << " " << stKin + pot << endl; // Can comment if you don't want to save potentials
            lps++;

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
            //if ((k % (int)1e5) == 0)
            //    cout << "Coefficient progress: " << (double)k / Niter * 100 << "%.       " << '\r';
        }

        cout << endl;
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

        cout << "Accepted steps: " << 100*(double)acc/Niter << "%." << endl;
        cout << "Mean energy: std: " << setprecision(10) << stKinCum << " pm " << sigmaSt << ", KFC: " << jfKinCum << " pm " << sigmaJf << endl;
        cout << "Mean potential: " << setprecision(10) << potCum << " pm " << sigmaPot << endl;
        cout << "Total E: " << setprecision(20) << Eavg << " pm " << sigma << endl;

        // COEFFICIENT CALCULATION (trust, copypaste from last year's code)
        int k = 0;
        double num, den = 0;
        for (int k = 0; k < Niter + 1; k++)
        {
            den += (stKinArr[k] + potArr[k] - Eavg)*(stKinArr[k] + potArr[k] - Eavg);
        }

        k = 0;
        for (int k = 0; k < Ntau; k++)
        {
            num = 0;
            for (int j = 0; j < Niter + 1 - k; j++)
            {
                num += (stKinArr[j] + potArr[j] - Eavg)*(stKinArr[j+k] + potArr[j+k] - Eavg);
            }
            coeffs[k] = num/den;
            fidC << k << " " << coeffs[k] << endl; // Save coefficients
        }
        stKinArr.clear();
        jfKinArr.clear();
        potArr.clear();
    }
    fidV.close();
    fidC.close();
}
