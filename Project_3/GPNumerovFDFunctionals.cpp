#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ctime>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

//#include <Eigenvalues> // For Eigen, just need to add libraries in a folder and add folder to include path
//#include <Dense>

using namespace std;
//using namespace Eigen; // For Eigen library

ofstream ofid;
ofstream ofidAlpha;

double Na, h, r0;
int Nmesh;
double *pPsi, *hRho;
// pPsi is wf output of Numerov/FD at previous step, hRho is the density that goes into the Hartree potential,
// which is equal to the "old" in alpha*pPsi*pPsi + (1-alpha)*hRho at next step
double *ddPsi; // Double derivative

double pot(int i) // Potential
{
    double r2 = (r0+i*h)*(r0+i*h);
	return 0.5*r2 + Na*hRho[i]/r2;
}

double smps2(double* y1, double* y2, int N, double h) // Simpson integration with product of 2 functions
{
    double ans = y1[0]*y2[0];
    for (int i = 1; i < N - 1; i++)
        (i % 2) ? ans += 4*y1[i]*y2[i] : ans += 2*y1[i]*y2[i];

    ans += y1[N-1]*y2[N-1];
    ans *= h/3;
    return ans;
}

double smpsEfcn(double* Vifcn, int N, double h) // Functionals evaluation with simpson, saves also Vint functional to use it in mu-Vint
{
    double T, Ve, Vi, r = r0;
    double factor;
    T = pPsi[0]*ddPsi[0];
    Ve = pPsi[0]*pPsi[0]*r*r;
    Vi = pPsi[0]*pPsi[0]*pPsi[0]*pPsi[0]/r/r;
    for (int i = 1; i < N - 1; i++)
    {
        (i % 2) ? factor = 4 : factor = 2;
        r = r0 + i*h;
        T += factor*pPsi[i]*ddPsi[i];
        Ve += factor*pPsi[i]*pPsi[i]*r*r;
        Vi += factor*pPsi[i]*pPsi[i]*pPsi[i]*pPsi[i]/r/r;
    }
    r = r0 + (N-1)*h;
    T += pPsi[N-1]*ddPsi[N-1];
    Ve += pPsi[N-1]*pPsi[N-1]*r*r;
    Vi += pPsi[N-1]*pPsi[N-1]*pPsi[N-1]*pPsi[N-1]/r/r;
    T *= h/3 * (-.5);
    Ve *= h/3 * .5;
    Vi *= h/3 * .5*Na;

    *Vifcn = Vi;

    return T+Ve+Vi;
}

double smps4(double* y1, double* y2, double* y3, double* y4, int N, double h) // Simpson integration with 4 fcns multiplied
{
    double ans = y1[0]*y2[0]*y3[0]*y4[0];
    for (int i = 1; i < N - 1; i++)
        (i % 2) ? ans += 4*y1[i]*y2[i]*y3[i]*y4[i] : ans += 2*y1[i]*y2[i]*y3[i]*y4[i];

    ans += y1[N-1]*y2[N-1]*y3[N-1]*y4[N-1];
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

double kPrev2, kPrev1, kCurr, uPrev2, uPrev1, uCurr;

double numerov(double E) // Numerov calculator
{
	double norm = 0;
	kPrev2 = 2 * (E - pot(0));
	kPrev1 = 2 * (E - pot(1));

	pPsi[0] = r0;
	pPsi[1] = r0+h;

	norm += pPsi[0]*pPsi[0];
	norm += 4*pPsi[1]*pPsi[1];

	int i = 2;
	while (i < Nmesh - 1)
	{
		kCurr = 2 * (E - pot(i));
		pPsi[i] = (pPsi[i-1] * (2. - 5*h*h/6.*kPrev1) - pPsi[i-2]*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
		(i % 2) ? norm += 4*pPsi[i]*pPsi[i] : norm += 2*pPsi[i]*pPsi[i];

		kPrev2 = kPrev1;
		kPrev1 = kCurr;
		i++;
	}

	kCurr = 2 * (E - pot(i));
	pPsi[i] = (pPsi[i-1] * (2. - 5*h*h/6.*kPrev1) - pPsi[i-2]*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
	norm += pPsi[i]*pPsi[i];

	norm *= h/3;

	for (i = 0; i < Nmesh; i++)
    {
        pPsi[i] /= sqrt(norm);
    }

    return pPsi[Nmesh - 1];
}

double tridet(double t) // Determinant of tridiagonal matrix of this problem with our method
{
    double fm2 = 0.;
    double fm1 = 1.;
    double fm0;
    for (int i = 0; i < Nmesh; i++)
    {
        fm0 = (2 + 2*h*h*pot(i) - 2*h*h*t)*fm1 - fm2;
        fm2 = fm1;
        fm1 = fm0;
    }
    return fm0;
}

gsl_vector *diag, *offd, *zs, *wf;

double findif(double Et, double p0, bool my1) // FD with our methods (1) (my1 = 1) and (2) (my1 = 0)
{
    double Ec = 0.;
    double Es = .1;
    double dlo = tridet(Ec);
    Ec += Es;
    double dhi = tridet(Ec);

    while (dlo*dhi > 0) // Just false position method with determinant of characteristic polynomial
    {
        dlo = dhi;
        Ec += Es;
        dhi = tridet(Ec);
    }

    double Ep = Ec - Es;

    int i = 0;
    double En;
    double dtemp;
    int maxCycles = 1e6;
    while (i < maxCycles && abs(Ec - Ep) > Et)
    {
        En = Ec - dhi*(Ec - Ep)/(dhi - dlo);
        dtemp = tridet(En);
        if (dtemp - dlo == 0 || dhi - dlo == 0)
        {
            i = maxCycles;
        }
        else if (dtemp * dlo < 0)
        {
            dhi = dtemp;
            Ec = En;
        }
        else if (dtemp * dhi < 0)
        {
            dlo = dtemp;
            Ep = En;
        }

        i++;
        if (i >= maxCycles)
            cout << "Overflow" << endl;
    }

    double norm = 0;
    if (my1)
    {
        // Our 1st implementation of FD
        pPsi[Nmesh-1] = p0; // Last point
        double V = 1./h/h + pot(Nmesh-2);
        pPsi[Nmesh-2] = 2*h*h*(V-En)*pPsi[Nmesh-1];

        norm += 4*pPsi[Nmesh-1]*pPsi[Nmesh-1];
        norm += pPsi[Nmesh-2]*pPsi[Nmesh-2];
        for (int j = Nmesh-3; j > -1; j--)
        {
            V = 1./h/h + pot(j+1);
            pPsi[j] = 2*h*h*(V-En)*pPsi[j+1] - pPsi[j+2];
            (j % 2) ? norm += 4*pPsi[j]*pPsi[j] : norm += 2*pPsi[j]*pPsi[j];
        }
        norm -= pPsi[0]*pPsi[0]; // Cause it added it twice in prev loop, should be added once cause it's last point
    }
    else
    {
    	// Our 2nd implementation of FD
        gsl_vector_set_zero(zs);
        double V;
        for (int j = 0; j < Nmesh-1; j++)
        {
            V = 1./h/h + pot(j) - En;
            gsl_vector_set(diag, j+1, V); // Diagonals
            gsl_vector_set(offd, j+1, -0.5/h/h); // Off diagonals
        }
        V = 1./h/h + pot(Nmesh-1);
        gsl_vector_set(diag, Nmesh, V);

        gsl_vector_set(diag, 0, 1.);
        gsl_vector_set(offd, 0, -1.);
        gsl_vector_set(zs, 0, -h);

        gsl_linalg_solve_symm_tridiag(diag, offd, zs, wf); // Solve linear system

        pPsi[0] = gsl_vector_get(wf, 1);
        norm += pPsi[0]*pPsi[0];
        for (int j = 1; j < Nmesh - 1; j++)
        {
            pPsi[j] = gsl_vector_get(wf, j+1);
            (j % 2) ? norm += 4*pPsi[j]*pPsi[j] : norm += 2*pPsi[j]*pPsi[j];
        }
        pPsi[Nmesh-1] = gsl_vector_get(wf, Nmesh);
        norm += pPsi[Nmesh-1]*pPsi[Nmesh-1];
    }

    norm *= h/3;

    for (int j = 0; j < Nmesh - 1; j++) // Normalize wf
    {
        pPsi[j] /= sqrt(norm);
    }

    return En;
}

double fpm(double Emin, double Estep, double Ethre) // False position method for Numerov
{
    double Ecurr = Emin;
    double ulo = numerov(Ecurr);
    Ecurr = Emin + Estep;
    double uhi = numerov(Ecurr);

    while (ulo*uhi > 0)
    {
        ulo = uhi;
        Ecurr += Estep;
        uhi = numerov(Ecurr);
    }

    double Eprev = Ecurr - Estep;

    int i = 0;
    double Enew;
    double utemp;
    int maxCycles = 1e5;
    while (i < maxCycles && abs(Ecurr - Eprev) > Ethre)
    {
        Enew = Ecurr - uhi*(Ecurr - Eprev)/(uhi - ulo);
        utemp = numerov(Enew);
        if (utemp - ulo == 0 || uhi - ulo == 0)
        {
            i = maxCycles;
        }
        else if (utemp * ulo < 0)
        {
            uhi = utemp;
            Ecurr = Enew;
        }
        else if (utemp * uhi < 0)
        {
            ulo = utemp;
            Eprev = Enew;
        }

        i++;
        if (i >= maxCycles)
            cout << "Overflow" << endl;
    }

    return Enew;
}

/*gsl_matrix *Evecs; // For GSL diagonalization
gsl_vector *Evals;
gsl_matrix *M;

double GSLfindif(void)
{
    gsl_matrix_set_zero(M); // Initialize matrices
    double V;
    double od = -.5/h/h;
    for (int i = 0; i < Nmesh; i++)
    {
        V = 1./h/h + pot(i);
        gsl_matrix_set(M, i, i, V);
        if (i > 0)
        {
            gsl_matrix_set(M, i, i-1, od);
        }
        if (i < Nmesh-2)
        {
            gsl_matrix_set(M, i, i+1, od);
        }
    }

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(Nmesh); // Diagonalize
    gsl_eigen_symmv(M, Evals, Evecs, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(Evals, Evecs, GSL_EIGEN_SORT_VAL_ASC);

    double norm = 0; // Calculate normalization
    pPsi[0] = gsl_matrix_get(Evecs, 0, 0);
    norm += pPsi[0]*pPsi[0];
    for (int j = 1; j < Nmesh - 1; j++)
    {
        pPsi[j] = gsl_matrix_get(Evecs, j, 0);
        (j % 2) ? norm += 4*pPsi[j]*pPsi[j] : norm += 2*pPsi[j]*pPsi[j];
    }
    pPsi[Nmesh-1] = gsl_matrix_get(Evecs, Nmesh-1, 0);
    norm += pPsi[Nmesh-1]*pPsi[Nmesh-1];

    norm *= h/3;

    for (int j = 0; j < Nmesh - 1; j++) // Normalize wf
    {
        pPsi[j] /= sqrt(norm);
    }

    return gsl_vector_min(Evals);
}*/

/*double EIGENfindif(void) // EIGEN library diagnoalizer (same structure as above)
{
    MatrixXd mat(Nmesh, Nmesh);
    double V;
    double od = -.5/h/h;
    for (int i = 0; i < Nmesh; i++)
    {
        V = 1./h/h + pot(i);
        mat(i, i) = V;
        if (i > 0)
        {
            mat(i, i-1) = od;
        }
        if (i < Nmesh-2)
        {
            mat(i, i+1) = od;
        }
    }

    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(mat);

    double norm = 0;
    pPsi[0] = es.eigenvectors().col(0)[0];
    norm += pPsi[0]*pPsi[0];
    for (int j = 1; j < Nmesh - 1; j++)
    {
        pPsi[j] = es.eigenvectors().col(0)[j];
        (j % 2) ? norm += 4*pPsi[j]*pPsi[j] : norm += 2*pPsi[j]*pPsi[j];
    }
    pPsi[Nmesh-1] = es.eigenvectors().col(0)[Nmesh-1];
    norm += pPsi[Nmesh-1]*pPsi[Nmesh-1];

    norm *= h/3;

    for (int j = 0; j < Nmesh - 1; j++)
    {
        pPsi[j] /= sqrt(norm);
    }

    return es.eigenvalues()[0];
}*/

int main()
{
	double rmax = 6;
	double Emin = 0.01;
	double Estep = .1;
	double Ethre = 1e-13;
	int maxCycles = 1e5;
	double alpha = 0;
	int scIter = 1;

	Na = 1;
	h = 1e-3;

	ofid.open("./GPwf.txt");

    r0 = h;
    Nmesh = (int)ceil(rmax / h);
    !(Nmesh % 2) ? Nmesh++ : true; // Odd number of points in mesh (due to using simpson)
    cout << Nmesh << " points with h = " << h << endl;
    pPsi = new double[Nmesh];
    hRho = new double[Nmesh];
    diag = gsl_vector_alloc(Nmesh+1);
    offd = gsl_vector_alloc(Nmesh);
    zs = gsl_vector_alloc(Nmesh+1);
    wf = gsl_vector_alloc(Nmesh+1);

    ddPsi = new double[Nmesh]; // Second deriavtive of wf

    double Vifcn, mu, Efcn1, Efcn2; // Functionals initialization
    double etime; // Time taken

	// For GSL diagonalization
	// LEAVE COMMENTED IF GSL DIAGONALIZATION IS UNUSED, for small h may have problems allocating memory
    /*M = gsl_matrix_alloc(Nmesh, Nmesh);
    Evecs = gsl_matrix_alloc(Nmesh, Nmesh);
    Evals = gsl_vector_alloc(Nmesh);*/

    double r;
    for (int i = 0; i < Nmesh; i++) // Wf ansatz
    {
        r = i*h+r0;
        pPsi[i] = 2./pow(M_PI, .25)*r*exp(-r*r/2);
        hRho[i] = 0.;
    }

    for (int i = 0; i < scIter; i++)
    {
        for (int k = 0; k < Nmesh; k++) // Set density to be used in Hartree potential (at 1st step, hRho = 0)
        {
            r = k*h+r0;
            hRho[k] = alpha*pPsi[k]*pPsi[k] + (1.-alpha)*hRho[k];
        }

        clock_t begin = clock();
        mu = fpm(Emin, Estep, Ethre); // Various minimization algorithms
        //mu = findif(1e-10, 1e-20, 1);
        //mu = GSLfindif();
        //mu = EIGENfindif();
        etime = (double)(clock() - begin) / CLOCKS_PER_SEC;
        cout << "Step " << i+1 << "; mu = " << setprecision(15) << mu << endl;

        ddf(pPsi, ddPsi, Nmesh, h); // Compute derivative
        Efcn1 = smpsEfcn(&Vifcn, Nmesh, h); // Compute functionals
        Efcn2 = mu - Vifcn;
        cout << "Functionals diff: " << abs(Efcn1 - Efcn2) << endl;
        ofidAlpha << fixed << setprecision(15) << i << "\t" << Efcn1 << "\t" << abs(Efcn1 - Efcn2) << endl; // Save iteration, functional and difference
    }

    for (int i = 0; i < Nmesh; i++) // Save wf and potential
    {
        r = i*h+r0;
        ofid << fixed << setprecision(15) << r << " " << pPsi[i] << " " << 0.5*r*r + Na*hRho[i]/r/r << "\n";
    }

	ofid.close();
	ofidAlpha.close();
}
