// TO USE THIS, YOU NEED TO COMPILE IT TOGETHER WITH THE FILES
// THAT YOU CAN DOWNLOAD FROM https://www.alglib.net/download.php
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ctime>

#include "linalg.h"
#include "ap.h"

using namespace std;
using namespace alglib;

ofstream ofid;
ofstream ofidfcn;
ifstream initfid;

double Ne, h, r0, rmax; // Number of electrons, mesh step, initial point, final point
int Nmesh; // Number of points in mesh
double *pPsi, *hRho, *hRhoPart, *Ucoul; // WFs, density, density cumulator, coulomb potential
// DENSITIES ARE DEFINED AS RHOTILDE IN REPORT, THAT IS RHO*4*M_PI*r*r
double *ddPsi, *psiddPsi, *angPsi2; // Second derivative of WF, WF*ddWF cumulator, WF^2/r^2 cumulator
double rs, Rc, Rc2, Rc3; // WS radius and Jellium radius to the pow of 1, 2 and 3
int levels[10][2], chlvs[10][2]; // Fixed levels to occupy, levels to check energy of (for interchanges)
double energies[10]; // Energies of levels to occupy
int nlevels, nchlvs; // # of occupied levels, # of levels to check
double lnum; // Maximum angular momentum
double minen, maxen; // Parameters for prelimiary check

double potJel(double r, double l) // Jellium potential with angular momentum
{
    if (r < Rc)
        return 0.5*Ne/Rc3*(r*r-3*Rc2) + 0.5*l*(l+1)/r/r;
    else
        return -Ne/r + 0.5*l*(l+1)/r/r;
}

double potExt(double r) // Jellium potential
{
    if (r < Rc)
        return 0.5*Ne/Rc3*(r*r-3*Rc2);
    else
        return -Ne/r;
}

double ers, rp; // Just some useful global definitions
double gamKS = -0.103756, beta1 = 0.56371, beta2 = 0.27358; // Ortiz-Ballone parameters
double potXC(int idx) // Exchange + correlation potential
{
    rp = r0 + idx*h;
    ers = cbrt(3.*rp*rp/hRho[idx]);
    return gamKS*(1.+7./6.*beta1*sqrt(ers)+4./3.*beta2*ers)/(1.+beta1*sqrt(ers)+beta2*ers)/(1.+beta1*sqrt(ers)+beta2*ers) - cbrt(.75*hRho[idx]/rp/rp/M_PI/M_PI);
}

void evalUcoul(double *dst) // COULOMB CALCULATOR
{
    double ans = 0;
    double r;
    for (int index = 0; index < Nmesh; index+=2) // Odd number of points for both integrals (including r = 0)
    {
        // pPsi[index] at r = r0 + index*h, need to get up to this coordinate for first integral
        ans = 0;
        if (index == 0)
        {
            ans = dst[0]/r0;
        }
        else
        {
            for (int i = 2; i < index+1; i++) // Simpson without initial 3 points (including r = 0)
            {
                r = r0 + i*h;
                (i % 2) ? ans += 4*dst[i] : ans += 2*dst[i];
                //             = 4*realRho*r*r*4*M_PI
            }
            ans -= dst[2];
            ans += 1.125*(3*dst[0] + 3*dst[1] + dst[2]);
            ans /= (r0+index*h);
            ans -= dst[index]/(r0+index*h);
        }

        if (index == Nmesh - 2)
        {
            ans += 1.5*(dst[index]/(r0 + index*h) + dst[index + 1]/(r0 + (index+1)*h));
        }
        else
        {
            for (int i = index; i < Nmesh - 3; i++) // Simpson without last 3 points
            {
                r = r0 + i*h;
                (i % 2) ? ans += 4*dst[i]/r : ans += 2*dst[i]/r;
            }
            ans -= dst[index]/(r0+index*h);
            ans -= dst[Nmesh-4]/(r0+(Nmesh-4)*h);

            // Simpson 3/8 rule, though with 9/8 cause at the end do *h/3
            ans += 1.125*(dst[Nmesh-1]/(r0+(Nmesh-1)*h) + 3*dst[Nmesh-2]/(r0+(Nmesh-2)*h) + 3*dst[Nmesh-3]/(r0+(Nmesh-3)*h) + dst[Nmesh-4]/(r0+(Nmesh-4)*h));

        }

        ans *= h/3;
        Ucoul[index] = ans;
    }

    for (int index = 1; index < Nmesh; index+=2) // Odd number of points for both integrals, including r = 0
    {
        ans = 0;
        for (int i = 0; i < index+1; i++)
        {
            r = r0 + i*h;
            (i % 2) ? ans += 2*dst[i] : ans += 4*dst[i];
        }
        ans /= (r0+index*h);
        ans -= dst[index]/(r0+index*h); // Same as moving this one line above and add r^2 instead of r only

        if (index == Nmesh - 1)
        {
            ans += 3*dst[Nmesh-1]/(r0+(Nmesh-1)*h);
        }
        else
        {
            for (int i = index; i < Nmesh; i++)
            {
                r = r0 + i*h;
                (i % 2) ? ans += 2*dst[i]/r : ans += 4*dst[i]/r;
            }
            ans -= dst[index]/(r0+index*h);
            ans -= dst[Nmesh-1]/(r0+(Nmesh-1)*h);
        }

        ans *= h/3;
        Ucoul[index] = ans;
    }
}

double tpz(int index) // Trapezoidal rule for Ucoulomb (just as test)
{
    double ans = 0;
    double r;
    for (int i = 0; i < index; i++) // Consider also r = 0
    {
        r = r0 + i*h;
        ans += 0.5*hRho[i];
    }
    ans /= (r0+index*h);
    ans += hRho[index]/4/(r0+index*h);

    ans += hRho[index]/4/(r0+index*h);
    for (int i = index+1; i < Nmesh-1; i++)
    {
        r = r0 + i*h;
        ans += 0.5*hRho[i]/r;
    }
    ans += hRho[Nmesh-1]/4/(r0+(Nmesh-1)*h);

    ans *= h/2;
    return 4*ans;
}

void evalUcoulT(void) // For trapezoidal rule
{
    for (int i = 0; i < Nmesh; i++)
    {
        Ucoul[i] = tpz(i);
    }
}

double reFcn(bool verbose) // Calculate potentials of E[rho]
{
    double ansdd = 0, ansll = 0, ansExt = 0, ansCoul = 0, ansX = 0, ansC = 0;
    for (int i = 0; i < Nmesh; i++)
    {
        rp = r0 + h*i;
        ers = cbrt(3.*rp*rp/hRhoPart[i]);
        if (i % 2)
        {
            // Kinetic
            ansll += 2*angPsi2[i];
            ansdd -= 2*psiddPsi[i]; // note minus
            // Potentials
            ansExt += 2*potExt(rp)*hRhoPart[i];
            ansCoul += Ucoul[i]*hRhoPart[i];
            ansX -= 1.5*cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i]; // note minus and coeff 2*3/4
            ansC += 2*gamKS / (1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];
        }
        else
        {
            // Kinetic
            ansll += 4*angPsi2[i];
            ansdd -= 4*psiddPsi[i];
            // Potentials
            ansExt += 4*potExt(rp)*hRhoPart[i];
            ansCoul += 2*Ucoul[i]*hRhoPart[i];
            ansX -= 3*cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i];
            ansC += 4*gamKS / (1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];
        }
    }
    int i = Nmesh-1;
    ansll -= angPsi2[i];
    ansdd += psiddPsi[i];

    ansExt -= potExt(rp)*hRhoPart[i];
    ansCoul -= .5*Ucoul[i]*hRhoPart[i];
    ansX += .75*cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i]; // minus 2*3/4
    ansC -= gamKS / (1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];

    if (verbose)
    {
        cout << "Angular kinetic: " << 0.5*ansll*h/3 << endl;
        cout << "Radial kinetic: " << 0.5*ansdd*h/3 << endl;
        cout << "External potential: " << ansExt*h/3 << endl;
        cout << "Coulomb potential: " << ansCoul*h/3 << endl;
        cout << "Exchange potential: " << ansX*h/3 << endl;
        cout << "Correlation potential: " << ansC*h/3 << endl;
    }

    return (0.5*ansll+0.5*ansdd+ansExt+ansCoul+ansX+ansC)*h/3;
}

double epsFcn(bool verbose) // Evaluate functional E_epsilon[rho]
{
    double ansEs = 0, ansCoul = 0, ansX = 0, ansC = 0;
    for (int i = 0; i < nlevels; i++)
    {
        ansEs += 2*(2*levels[i][1]+1)*energies[i];
    }

    for (int i = 0; i < Nmesh; i++)
    {
        rp = r0 + i*h;
        ers = cbrt(3.*rp*rp/hRhoPart[i]);
        if (i % 2)
        {
            ansCoul += Ucoul[i]*hRhoPart[i];
            ansX -= .5*cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i];
            ansC += 2*gamKS*(beta1*sqrt(ers)/6+beta2*ers/3)/(1.+beta1*sqrt(ers)+beta2*ers)/(1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];
        }
        else
        {
            ansCoul += 2*Ucoul[i]*hRhoPart[i];
            ansX -= cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i];
            ansC += 4*gamKS*(beta1*sqrt(ers)/6+beta2*ers/3)/(1.+beta1*sqrt(ers)+beta2*ers)/(1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];
        }
    }
    int i = Nmesh-1;
    ansCoul -= .5*Ucoul[i]*hRhoPart[i];
    ansX += .25*cbrt(.75*hRhoPart[i]/rp/rp/M_PI/M_PI)*hRhoPart[i];
    ansC -= gamKS*(beta1*sqrt(ers)/6+beta2*ers/3)/(1.+beta1*sqrt(ers)+beta2*ers)/(1.+beta1*sqrt(ers)+beta2*ers)*hRhoPart[i];

    if (verbose)
    {
        cout << "SP energies sum: " << ansEs << endl;
        cout << "Coulomb correction: " << ansCoul*h/3 << endl;
        cout << "Exchange correction: " << ansX*h/3 << endl;
        cout << "Correlation correction: " << ansC*h/3 << endl;
    }

    return ansEs - (ansCoul + ansX + ansC)*h/3;
}

double smps2(double* y1, double* y2, int N, double h) // General Simpson integration considering y(r = 0) = 0
{
    double ans = 0;
    for (int i = 0; i < N; i++)
        (i % 2) ? ans += 2*y1[i]*y2[i] : ans += 4*y1[i]*y2[i];

    ans -= y1[N-1]*y2[N-1];
    ans *= h/3;
    return ans;
}

double smps(double* y1, int N, double h) // Same as above, with just a fcn
{
    double ans = 0;
    for (int i = 0; i < N; i++)
        (i % 2) ? ans += 2*y1[i]: ans += 4*y1[i];

    ans -= y1[N-1];
    ans *= h/3;
    return ans;
}

void ddf(double* f, double* otp, int N, double h) // O(h^4) second derivative
{
    int i = 0;
    otp[i] = 45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5];
    otp[i] /= (12*h*h);
    i++;
    otp[i] = 45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5];
    otp[i] /= (12*h*h);
    for (i = 2; i < N-2; i++)
    {
        otp[i] = -f[i-2]+16*f[i-1]-30*f[i]+16*f[i+1]-f[i+2];
        otp[i] /= (12*h*h);
    }
    i = N-2;
    otp[i] = 45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5];
    otp[i] /= (12*h*h);
    i++;
    otp[i] = 45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5];
    otp[i] /= (12*h*h);
}

static real_1d_array diag; // For the diagonalizer
static real_1d_array offd;
static real_2d_array evec;

void findif(double Emin, double l, int nneed, int disc) // FD eqn solver
{
    double V, r;
    diag.setlength(Nmesh);
    for (int j = 0; j < Nmesh-2; j++)
    {
        V = 1./h/h + potJel(r0+j*h, l) + disc*(potXC(j) + Ucoul[j]); // Diagonal elements
        diag[j] = V;
        offd[j] = -0.5/h/h; // Off diagonal elements
    }

    V = 1./h/h + potJel(r0+(Nmesh-1)*h, l) + disc*(potXC(Nmesh-1) + Ucoul[Nmesh-1]);
    diag[Nmesh-1] = V;

    smatrixtdevdi(diag, offd, Nmesh, 2, 0, 4, evec); // Diagonalizer finding eigenvalues with indexes 0, 1, 2, 3, 4
    //smatrixtdevd(diag, offd, Nmesh, 2, evec);
    int sel = 0;
    while (diag[sel] < Emin) // This is because the eigensolver gives the 1st eigenvalue as a large negative number
    {
        sel++;
    }

    double norm;
    for (int i = 0; i < nneed; i++) // Loop over the needed energy levels (principal quantum number)
    {
        norm = 0;
        for (int j = 0; j < Nmesh; j++) // Save eigenvector and normalization
        {
            pPsi[j] = evec[j][sel+i];
            (j % 2) ? norm += 2*pPsi[j]*pPsi[j] : norm += 4*pPsi[j]*pPsi[j];
        }
        norm -= pPsi[Nmesh-1]*pPsi[Nmesh-1];

        norm *= h/3;

        for (int j = 0; j < Nmesh - 1; j++) // Accumulate angular kinetic term and density
        {
            pPsi[j] /= sqrt(norm);
            r = r0 + j*h;
            hRhoPart[j] += 2*(2*l+1)*pPsi[j]*pPsi[j];
            angPsi2[j] += 2*(2*l+1)*l*(l + 1.)*pPsi[j]*pPsi[j]/r/r;
        }
        ddf(pPsi, ddPsi, Nmesh, h);
        for (int k = 0; k < Nmesh; k++) // Accumulate derivative term
        {
            psiddPsi[k] += 2*(2*l+1)*pPsi[k]*ddPsi[k];
        }
        energies[(int)(lnum*i + l)] = diag[sel+i]; // Save and output energy
        cout << setprecision(15) << "Level n: " << i << ", l: " << l << " has E = " << energies[(int)(lnum*i + l)] << endl;
    }
}

double findifCheck(double Emin, double l, int n, int disc) // Check of levels to see interchanges (only eigenvalues are needed)
{
    double V, r;
    diag.setlength(Nmesh);
    for (int j = 0; j < Nmesh-2; j++)
    {
        V = 1./h/h + potJel(r0+j*h, l) + disc*(potXC(j) + Ucoul[j]);
        diag[j] = V;
        offd[j] = -0.5/h/h;
    }

    V = 1./h/h + potJel(r0+(Nmesh-1)*h, l) + disc*(potXC(Nmesh-1) + Ucoul[Nmesh-1]);
    diag[Nmesh-1] = V;

    smatrixtdevdi(diag, offd, Nmesh, 0, 0, n+1, evec);
    int sel = 0;
    while (diag[sel] < Emin)
    {
        sel++;
    }

    return diag[sel];
}

void KScalc(double Emin, bool checkSphere, int disc) // KS eqns solver
{
    for (int k = 0; k < Nmesh; k++) // Begin by setting cumulators to 0
    {
        hRhoPart[k] = 0;
        psiddPsi[k] = 0;
        angPsi2[k] = 0;
    }
    double r;
    int nneed;

    for (int nl = 0; nl < lnum; nl++) // Loop over requested leves
    {
        nneed = 1;
        if (nl < nlevels - lnum)
        {
            nneed = 2;
        }
        findif(Emin, nl, nneed, disc);
    }
    hRhoPart[Nmesh-1] = hRhoPart[Nmesh-2];

    /*if (checkSphere) This checks levels (nmax+1, lmax) and (nmax, lmax+1) (NOT WORKING HERE)
    {
        int lch;
        double Eerr;
        for (int nch = 0; nch < (nlevels / 4)+2; nch++)
        {
            lch = 0;
            for (int k = 0; k < nlevels; k++)
            {
                if (levels[k][0] == nch && levels[k][1] == lch)
                {
                    Emin = energies[k] + 1e-8;
                    lch = levels[k][1] + 1;
                }
            }
            Eerr = fpmCheck(Emin, Estep, Et, nch, lch);
            for (int k = 0; k < nlevels; k++)
            {
                if (Eerr < energies[k])
                {
                    cout << "!-- WARNING: E" << nch << lch << " = " << Eerr << " < E" << levels[k][0] << levels[k][1] << " = " << energies[k] << " ---!" << endl;
                }
            }
        }
    }*/

    if (checkSphere) // Checks energies of requested levels
    {
        int nch, lch;
        double Eerr;
        for (int i = 0; i < nchlvs; i++)
        {
            nch = chlvs[i][0];
            lch = chlvs[i][1];
            for (int k = 0; k < nlevels; k++)
            {
                if (levels[k][0] == nch && levels[k][1]+1 == lch)
                {
                    Emin = energies[k] + 1e-8;
                }
            }
            Eerr = findifCheck(Emin, lch, nch, disc);
            for (int k = 0; k < nlevels; k++)
            {
                if (Eerr < energies[k]) // Output warning if levels interchange
                {
                    cout << "!-- WARNING: E" << nch << lch << " = " << Eerr << " < E" << levels[k][0] << levels[k][1] << " = " << energies[k] << " ---!" << endl;
                }
            }
        }
    }

    cout << "Density normalization: " << smps(hRhoPart, Nmesh, h) << endl; // Just some checks
    cout << "Density in last point: " << setprecision(15) << hRhoPart[Nmesh-1]/4/M_PI/(r0+(Nmesh-1)*h)/(r0+(Nmesh-1)*h) << endl;
}

double evalSplt(void) // Spillout evaluation (acts differently based on even or odd number of points)
{
    int fst = (int)floor(Rc/h);
    double ans = 0;
    if (fst % 2)
    {
        for (int i = fst; i < Nmesh; i++)
        {
            (i % 2) ? ans += 2*hRhoPart[i] : ans += 4*hRhoPart[i];
        }
        ans -= hRhoPart[fst];
        ans -= hRhoPart[Nmesh-1];
    }
    else
    {
        for (int i = fst; i < Nmesh-3; i++)
        {
            (i % 2) ? ans += 4*hRhoPart[i] : ans += 2*hRhoPart[i];
        }
        ans -= hRhoPart[fst];
        ans -= hRhoPart[Nmesh-4];
        ans += 1.125*(hRhoPart[Nmesh-4]+3*hRhoPart[Nmesh-3]+3*hRhoPart[Nmesh-2]+3*hRhoPart[Nmesh-1]);
    }
    return ans*h/3;
}

int main()
{
	rmax = 40;
	double Estep = .1; // Not needed
	double Ethre = 1e-15; // Not needed
	double alpha = .3; // .3 For Ne = 40, rmax = 40 2e-3 // .4 for Ne=20, rmax= 40 2e-3 // .4 for Ne=8 rmax = 40 2e-3
	int scIter = 100; // Iterations of SC procedure

    ofid.open("./KSrhoFDo.txt"); // Output files
	ofidfcn.open("./KSfcnFDo.txt");
	Ne = 40;
	nlevels = 6;
	lnum = 4;
	for (int i = 0; i < nlevels; i++) // Fixed energy levels
    {
        levels[i][0] = i / (int)lnum;
        levels[i][1] = i % (int)lnum;
    }
    chlvs[0][0] = 0; // Energy levels to check
    chlvs[0][1] = 4;
    nchlvs = 1;

	rs = 3.93; // WS radius of Jellium
	//rs = 4.86;
	Rc = rs*cbrt(Ne);
	Rc2 = Rc*Rc;
	Rc3 = rs*rs*rs*Ne;
	h = 2e-3;

    r0 = h;
    Nmesh = (int)ceil(rmax / h);
    (Nmesh % 2) ? Nmesh++ : true; // Set even number of points
    pPsi = new double[Nmesh];
    hRho = new double[Nmesh];
    hRhoPart = new double[Nmesh];
    ddPsi = new double[Nmesh];
    psiddPsi = new double[Nmesh];
    angPsi2 = new double[Nmesh];
    Ucoul = new double[Nmesh]; // Direct term, it's pretty convenient to have it already cause it's bad to calculate
    double *tempCoul = new double[Nmesh];
    diag.setlength(Nmesh);
    offd.setlength(Nmesh-1);
    evec.setlength(Nmesh, 3);
    double Efcn1, Efcn2;
    cout << "Rmax = " << rmax << ", h = " << h << ", Nmesh = " << Nmesh << ", Rc = " << Rc << ", Ne = " << Ne << endl;

    double tmp, r;
    for (int i = 0; i < Nmesh; i++) // Initialize some stuff
    {
        r = r0 + i*h;
        hRho[i] = 0.;
        hRhoPart[i] = 0;
        tempCoul[i] = 0;
    }
    evalUcoul(hRhoPart);
    int dst;
    bool parm;

    for (int i = 0; i < scIter; i++) // SC procedure
    {
        for (int k = 0; k < Nmesh; k++) // Mixing
        {
            hRho[k] = alpha*hRhoPart[k] + (1.-alpha)*hRho[k];
            Ucoul[k] = alpha*Ucoul[k] + (1.-alpha)*tempCoul[k];
            if (i==0) // 0 iteration
            {
                hRho[k] = 1.;
                Ucoul[k] = 0;
            }
        }
        minen = -4; // Preliminary check parameters (not really needed with FD)
        maxen = 1;
        dst = 0;
        if (i > 0)
        {
            minen = energies[0];
            maxen = energies[0];
            for (int k = 1; k < nlevels; k++)
            {
                if (energies[k] < minen)
                {
                    minen = energies[k];
                }
                if (energies[k] > maxen)
                {
                    maxen = energies[k];
                }
            }
            dst = 1; // This parameter distinguishes between using all potentials or Jellium potential only
        }
        else
        {
            for (int i = 0; i < nlevels; i++) // Output HO parameters for a check
            {
                cout << setprecision(15) << "H.O. level n: " << levels[i][0] << ", l: " << levels[i][0] << " has E = " << -1.5*Ne/Rc+sqrt(Ne/Rc3)*(1.5+2*levels[i][0]+levels[i][1]) << endl;
            }
        }
        cout << minen*1.1<< endl;
        cout << abs(maxen)*.1 << endl;

        cout << "Step " << i << ":" << endl;
        KScalc(minen*1.1, 0, dst);

        parm = 0;
        for (int k = 0; k < Nmesh; k++)
        {
            r = r0 + k*h;
            tempCoul[k] = Ucoul[k]; // Save previous Coulomb for mixing
            if (cbrt(3*r*r/hRhoPart[k]) < 1) // Check if rs > 1 always
                parm = 1;
        }
        if (parm)
            cout << "!--- Warning: rs < 1 ---!" << endl;

        evalUcoul(hRhoPart);

        // Evaluate and save functionals
        Efcn1 = reFcn(0);
        Efcn2 = epsFcn(0);
        cout << "Functionals diff: " << abs(Efcn1 - Efcn2) << ", fcn: " << Efcn1 << endl;
        cout << endl;
        ofidfcn << setprecision(15) << i << " " << abs(Efcn1 - Efcn2) << " " << Efcn1 << endl;
    }

    for (int i = 0; i < Nmesh; i++) // Save final density and potentials
    {
        r = r0 + i*h;
        ofid << setprecision(20) << hRhoPart[i]/4/M_PI/r/r << " " << potJel(r, 0) << " " << potXC(i) << " " << Ucoul[i] << endl;
    }

    cout << "Total energy of system: " << 3.*Ne*Ne/5/Rc + Efcn1 << endl; // Total energy considering Jellium self interaction

    double spillout = evalSplt();
    cout << "Spillout: " << spillout << ", polarizability: " << Rc3*(1.+spillout/Ne) << endl; // Spillout and polarizability

    Efcn1 = reFcn(1); // Output functionals at the end as checks
    Efcn2 = epsFcn(1);

	ofid.close();
	ofidfcn.close();
	// Bye bye
}

