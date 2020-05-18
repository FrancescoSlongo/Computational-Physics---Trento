// Computational Physics 2020 - Project 4: DFT Theory
// Master course in Physics - Università di Trento

// Kohm Sham Equation - Jellium

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <string.h>  // For memset
#include <iomanip>	// for setprecision, setw

using namespace std;

ofstream ofile;

// Parameters of the mesh
const double h = pow(10.0, -3);         // Mesh step
const double h2 = h*h;                  // Useful parameter
const double rmax = 25.0;                // Max value of the mesh
const int N = int(ceil(rmax/h));    // Number of mesh points

// Parameters of Jellium
const int Ne = 40;                   // Number of electrons
const double rs = 3.93;             // Wigner-Seitz radius 3.93 - 4.86
const double rhoI = 3.0/(4.0 * M_PI * pow(rs,3));
const double Rc = rs * pow(Ne, 1.0/3.0);
const double Rc2 = pow(Rc, 2);
const double Rc3 = rs*rs*rs*Ne;

// Parameters of correlation energy
const double gam = -0.103756;
const double beta1 = 0.56371;
const double beta2 = 0.27358;

// Other useful parameters
const double alpha = 0.1;

double Vext(int i)
{
    if(i*h <= Rc)
    {
        return 0.5* Ne *(i*i*h2 - 3.0 * Rc2)/Rc3;
    }
    else
    {
        return - Ne / (i*h);
    }
}

double Vcoul(int i, double * rho)
{
    double U1, U2;
    U1 = 0.0;
    U2 = 0.0;

    // Notice: at the end I have to divide by h
    if(i == 1)
    {
        // First integral
        U1 += rho[1]/2;
        // Second integral
        U2 += (rho[1] + rho[N - 1]*(N-1));
        for(int j = 1; j < N/2 - 2; j++)
        {
            U2 += 4.0*rho[2*j]*2*j;
            U2 += 2.0*rho[2*j+1]*(2*j+1);
        }
        U2 += 4.0*rho[N-2]*(N-2);
        U2 /= 3.0;
        U2 += (rho[N - 1]*(N - 1) + rho[N]*N)/2.0;
    }
    else if(i == N - 1)
    {
        // First integral
        U1 += rho[1] + rho[i]*i*i;
        for(int j = 1; j < N/2 - 2; j++)
        {
            U1 += 4.0*rho[2*j]*4*j*j;
            U1 += 2.0*rho[2*j+1]*(2*j+1)*(2*j+1);
        }
        U1 += 4.0*rho[N-2]*(N-2)*(N-2);
        U1 /= 3.0;
        U1 += rho[1];     // I add the first point in a different way
        U1 /= i;

        U2 += (rho[N - 1]*(N - 1) + rho[N]*N)/2.0;
    }
    else if(i == N)
    {
        // First integral
        U1 += rho[i]*i*i;
        for(int j = 1; j < N/2; j++)
        {
            U1 += 4.0*rho[2*j-1]*(2*j-1)*(2*j-1);
            U1 += 2.0*rho[2*j]*4*j*j;
        }
        U1 += 4.0*rho[N-1]*(N-1)*(N-1);
        U1 /= 3.0;
        U1 /= i;
    }
    else if((i%2 == 0) && (i != 1 && i != N - 1 && i != N))
    {
        // First integral
        U1 += rho[i]*i*i;
        for(int j = 1; j < i/2; j++)
        {
            U1 += 4.0*rho[2*j-1]*(2*j-1)*(2*j-1);
            U1 += 2.0*rho[2*j]*4*j*j;
        }
        U1 += 4.0*rho[i-1]*(i-1)*(i-1);
        U1 /= 3.0;
        U1 /= i;
        // Second integral
        U2 += (rho[i]*i + rho[N]*N);
        for(int j = i/2 + 1; j < N/2; j++)
        {
            U2 += 4.0*rho[2*j-1]*(2*j-1);
            U2 += 2.0*rho[2*j]*2.0*j;
        }
        U2 += 4.0*rho[N-1]*(N-1);
        U2 /= 3.0;
    } else if((i%2 == 1) && (i != 1 && i != N - 1 && i != N))
    {
        // First integral
        U1 += rho[1] + rho[i]*i*i;
        for(int j = 1; j < (i - 1)/2; j++)
        {
            U1 += 2.0*rho[2*j]*4*j*j;
            U1 += 4.0*rho[2*j+1]*(2*j+1)*(2*j+1);
        }
        U1 += 4.0*rho[i-1]*(i-1)*(i-1);
        U1 /= 3.0;
        U1 += rho[1];
        U1 /= i;
        // Second integral
        U2 += rho[i]*i + rho[N-1]*(N-1);
        for(int j = (i + 1)/2; j < N/2 -1; j++)
        {
            U2 += 4.0*rho[2*j]*2*j;
            U2 += 2.0*rho[2*j+1]*(2*j+1);
        }
        U2 += 4.0*rho[N - 2]*(N - 2);
        U2 /= 3.0;
        U2 += (rho[N - 1]*(N - 1) + rho[N]*N)/2.0;
    }
    return 4.0*M_PI*h2*(U1 + U2);
}

const double M_PI23 = 3.0/(M_PI * M_PI);
double Vx(double rho)
{
    return -pow((M_PI23 * rho), 1.0/3.0);
}

double Vc(double rho)
{
    double rstmp;
    rstmp = pow( (3.0 / (4.0 * M_PI * rho)) , 1.0/3.0);
    return gam / (1.0 + beta1*sqrt(rstmp) + beta2*rstmp);
}

double VcDe(double rho)
{
    double rstmp;
    rstmp = pow( (3.0 / (4.0 * M_PI * rho)) , 1.0/3.0);
    return gam * (beta1/6.0 * sqrt(rstmp)+ beta2/3.0 * rstmp) / (1.0 + beta1*sqrt(rstmp) + beta2*rstmp)/ (1.0 + beta1*sqrt(rstmp) + beta2*rstmp);
}

double Vc12(double rho)
{
    double rstmp;
    rstmp = pow( (3.0 / (4.0 * M_PI * rho)) , 1.0/3.0);
    return gam * (1.0 + 7.0/6.0*beta1*sqrt(rstmp) + 4.0/3.0*beta2*rstmp) / (1.0 + beta1*sqrt(rstmp) + beta2*rstmp)/ (1.0 + beta1*sqrt(rstmp) + beta2*rstmp);
}

double Amp(double * y) // Given a wave function, it calculate the normalization amplitude
{
	double A;
	// I take the square of the wave function for each point
	A = y[0]*y[0]+y[N]*y[N];
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        A += 4*y[2*m - 1]*y[2*m - 1];
        A += 2*y[2*m]*y[2*m];
	}
	A += 4*y[N-1]*y[N-1];
	A *= h/3.0;
	// I give back the amplitude
	return 1.0/sqrt(A);
}

double FunT(double * y)
{
    double T;
	// I calculate the V due to the harmonic potential
	T = 1.0 / (12.0 * h2) * ((45.0 * y[1] - 154.0 * y[2] + 214.0 * y[3] - 156.0 * y[4] + 61.0 * y[5] - 10.0 * y[6]) * y[1] + (45.0 * y[N - 1] - 154.0 * y[N - 2] + 214.0 * y[N - 3] - 156.0 * y[N - 4] + 61.0 * y[N - 5] - 10.0 * y[N - 6]) * y[N - 1]);
	for(int m = 1; m < N/2 - 1; m++)
	{
	    T += 1.0 / (3.0 * h2) * (-30.0 * y[2*m] + 16.0 * y[2*m + 1] + 16.0 * y[2*m - 1] - y[2*m + 2] - y[2*m - 2]) * y[2*m];
	    T += 1.0 / (6.0 * h2) * (-30.0 * y[2*m+1] + 16.0 * y[2*m + 2] + 16.0 * y[2*m] - y[2*m+3] - y[2*m - 1]) * y[2*m+1];
	}
	T += 1.0 / (3.0 * h2) * (-30.0 * y[N - 2] + 16.0 * y[N - 1] + 16.0 * y[N - 3] - y[N] - y[N - 4]) * y[N - 2];
	T *= -h/3.0;
	T += -h*(1.0 / (12.0 * h2) * ((45.0 * y[1] - 154.0 * y[2] + 214.0 * y[3] - 156.0 * y[4] + 61.0 * y[5] - 10.0 * y[6]) * y[1] + (45.0 * y[N - 1] - 154.0 * y[N - 2] + 214.0 * y[N - 3] - 156.0 * y[N - 4] + 61.0 * y[N - 5] - 10.0 * y[N - 6]) * y[N - 1]));; // I add another time the smallest derivative at the extremes
	// I give back the potential
	return T;
}

double FunVcen(double * y)
{
    double V;
	// I calculate the V due to the external potential
    V = y[1]*y[1]+y[N]*y[N]/(N*N);
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*y[2*m-1]*y[2*m-1]/((2*m-1)*(2*m-1));
        V += 2*y[2*m]*y[2*m]/(4*m*m);
	}
	V += 4*y[N - 1]*y[N - 1]/((N-1)*(N-1));
	V /= 3.0*h;
	// I give back the potential
	return V;
}

double FunVext(double * rho, double * VextVec)
{
    double V;
	// I calculate the V due to the external potential
    V = VextVec[N]* rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*VextVec[2*m-1]*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*VextVec[2*m]*rho[2*m]*(4*m*m);
	}
	V += 4*VextVec[N-1]*rho[N - 1]*(N-1)*(N-1);
	V *= 4.0 * M_PI * h * h2/3.0;
	// I give back the potential
	return V;
}

double FunVcoul(double * rho, double * VcoulVec)
{
    double V;
	// I calculate the V due to the external potential
    V = VcoulVec[N]*rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*VcoulVec[2*m-1]*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*VcoulVec[2*m]*rho[2*m]*(4*m*m);
	}
	V += 4*VcoulVec[N-1]*rho[N - 1]*(N-1)*(N-1);
	V *= 2.0 * M_PI * h * h2/3.0;
	// I give back the potential
	return V;
}

double FunVx(double * rho, double * VxVec)
{
    double V;
	// I calculate the V due to the external potential
    V = VxVec[N] * rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*VxVec[2*m-1]*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*VxVec[2*m]*rho[2*m]*(4*m*m);
	}
	V += 4*VxVec[N-1]*rho[N - 1]*(N-1)*(N-1);
	V *= M_PI * h * h2;
	// I give back the potential
	return V;
}

double FunVDex(double * rho, double * VxVec)
{
    double V;
	// I calculate the V due to the external potential
    V = VxVec[N] * rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*VxVec[2*m-1]*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*VxVec[2*m]*rho[2*m]*(4*m*m);
	}
	V += 4*VxVec[N-1]*rho[N - 1]*(N-1)*(N-1);
	V *= M_PI * h * h2/3.0;
	// I give back the potential
	return V;
}

double FunVc(double * rho)
{
    double V;
	// I calculate the V due to the external potential
    V = Vc(rho[N]) * rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*Vc(rho[2*m-1])*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*Vc(rho[2*m])*rho[2*m]*(4*m*m);
	}
	V += 4*Vc(rho[N-1])*rho[N - 1]*(N-1)*(N-1);
	V *= 4.0 * M_PI * h * h2 / 3.0;
	// I give back the potential
	return V;
}

double FunVDec(double * rho)
{
    double V;
	// I calculate the V due to the external potential
    V = VcDe(rho[N]) * rho[N]*N * N;
	for(int m = 1; m < N/2.0 - 1; m++)
	{
        V += 4*VcDe(rho[2*m-1])*rho[2*m-1]*(2*m-1)*(2*m-1);
        V += 2*VcDe(rho[2*m])*rho[2*m]*(4*m*m);
	}
	V += 4*VcDe(rho[N-1])*rho[N - 1]*(N-1)*(N-1);
	V *= 4.0 * M_PI * h * h2 / 3.0;
	// I give back the potential
	return V;
}

void PropInd(double E,  int l, double * y, double * VcentrVec, double * VextVec) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	// I calculate k2
	double k0; // Value of k2
	double k1; // Value of k2
	double k_1; // Value of k2
	// Initial conditions
	y[0] = 0.0;
	y[1] = pow(h, l + 1);
	y[2] = pow(2.0*h, l+1);
    k0 = 2.0*(E - VextVec[2]) - VcentrVec[2];
    k_1 = 2.0*(E - VextVec[1]) - VcentrVec[1];
    k1 = 2.0*(E - VextVec[3]) - VcentrVec[3];
	// Calculate next step of wf
	for(int i = 2; i < N - 1; i++)
    {
        y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k0) - y[i-1]*(1.0 + h2 /12.0* k_1)) / (1.0 + h2 / 12.0 * k1);
        k_1 = k0;
        k0 = k1;
        k1 = 2.0*(E - VextVec[i+2]) - VcentrVec[i+2];
    }
    y[N] =(y[N-1]*(2.0 - 5.0/6.0 * h2 * k0) - y[N-2]*(1.0 + h2 /12.0* k_1)) / (1.0 + h2 / 12.0 * k1);
}

void NumerovInd(int l, int nmax, double * y, double * rho,  double * VextVec, double * Esaved) // Method that implements the secant method with Numerov
{
	// Parameters of the secant method in the Numerov algorithm
	double DE = 0.01;      // Step of E

	// Auxiliary variables useful for Numerov
	double E, E1, E2, Eaux;
	double y1, y2, ytmp, tmp;

	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-13;
	int maxCycNum = 1e6;
	int sol = 0;

	// I set some of the variables
	E = -5.0;
	ytmp = 0.0;
	tmp = 0.0;
	int j; 	//Counter for the while cycle
	int m;
	double A; // Normalization amplitude

	// I precalculate the potentials that do not change
	double * VcentrVec = new double[N + 1];
	VcentrVec[0] = 0.0;
	for(int i = 1; i < N + 1; i++)
    {
        VcentrVec[i] = double(l * (l +1.0))/ (i*i*h2);
    }

	// Loop for energies
	while(sol < nmax){
		// I use Numerov
		PropInd(E, l, y, VcentrVec, VextVec);
		// Condition if the sign of the wavefunction changes at the extremity
		if(y[N]*tmp < 0)
		{
			// I prepare the variables for the method and I save the variables to continue the while loop
			E1 = E - DE;
			E2 = E;
			y1 = tmp;
			y2 = y[N];
			tmp = y[N];
			j = 0;
			// Loop for secant method
			while(j < maxCycNum && abs(E2 - E1) > resolution)
			{
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				PropInd(Eaux, l, y, VcentrVec, VextVec);

				// Condition to choose the next point
				if(y[N]*y1 > 0)
				{
					y1 = y[N];
					E1 = Eaux;
				}
				else if(y[N]*y1 < 0)
				{
					y2 = y[N];
					E2 = Eaux;
				}
				else if(y[N]*y1 == 0)
				{
					j = maxCycNum;
				}
				j++;
			}
			Esaved[0] += 2.0*(2.0*l+1.0)*Eaux;
			sol += 1;
			// Normalise the wavefunction and add the new density
			A = Amp(y);
			for(int m = 1; m < N + 1; m++)
			{
				y[m] = y[m]*A;
				rho[m] += (2.0*l + 1.0)/(2.0 * M_PI * m * m * h2)*y[m]*y[m]; // Update the density
			}
            Esaved[1] += (2.0*l+1)*FunT(y);
			Esaved[2] += (2.0*l + 1.0)*l*(l + 1.0)*FunVcen(y);
		} else {
		tmp = y[N];
		}
		E = E + DE;
	}
}

void PropKS(double E,  int l, double * y, double * VcentrVec, double * VextVec, double * VcoulVec, double * VxVec, double * Vc12Vec) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	// I calculate k2
	double k0; // Value of k2
	double k1; // Value of k2
	double k_1; // Value of k2
	// Initial conditions
	y[0] = 0.0;
	y[1] = pow(h, l + 1);
	y[2] = pow(2.0*h, l+1);
    k0 = 2.0*(E - VextVec[2] - VcoulVec[2] - VxVec[2] - Vc12Vec[2]) - VcentrVec[2];
    k_1 = 2.0*(E - VextVec[1] - VcoulVec[1] - VxVec[1] - Vc12Vec[1]) - VcentrVec[1];
    k1 = 2.0*(E - VextVec[3] - VcoulVec[3] - VxVec[3] - Vc12Vec[3]) - VcentrVec[3];
	// Calculate next step of wf
	for(int i = 2; i < N - 1; i++)
    {
        y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k0) - y[i-1]*(1.0 + h2 /12.0* k_1)) / (1.0 + h2 / 12.0 * k1);
        k_1 = k0;
        k0 = k1;
        k1 = 2.0*(E - VextVec[i+2] - VcoulVec[i+2] - VxVec[i+2] - Vc12Vec[i+2]) - VcentrVec[i+2];
    }
    y[N] =(y[N-1]*(2.0 - 5.0/6.0 * h2 * k0) - y[N-2]*(1.0 + h2 /12.0* k_1)) / (1.0 + h2 / 12.0 * k1);
}

void NumerovKS(int l,  int nmax, double * y, double * rho, double * VextVec, double * VcoulVec, double * VxVec, double * Vc12Vec, double * Esaved) // Method that implements the secant method with Numerov
{
	// Parameters of the secant method in the Numerov algorithm
	double DE = 0.01;      // Step of E

	// Auxiliary variables useful for Numerov
	double E, E1, E2, Eaux;
	double y1, y2, ytmp, tmp;

	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-13;
	int maxCycNum = 1e6;
	int sol = 0;

	// I set some of the variables
	E = -5.5;
	ytmp = 0.0;
	tmp = 0.0;
	int j; 	//Counter for the while cycle
	int m;
	double A; // Normalization amplitude

	// I precalculate the potentials that do not change
	double * VcentrVec = new double[N + 1];
	VcentrVec[0] = 0.0;
	for(int i = 1; i < N + 1; i++)
    {
        VcentrVec[i] = double(l * (l +1.0))/ (i*i*h2);
    }

	// Loop for energies
	while(sol < nmax){
		// I use Numerov
		PropKS(E, l, y, VcentrVec, VextVec, VcoulVec, VxVec, Vc12Vec);
		// Condition if the sign of the wavefunction changes at the extremity
		if(y[N]*tmp < 0)
		{
			// I prepare the variables for the method and I save the variables to continue the while loop
			E1 = E - DE;
			E2 = E;
			y1 = tmp;
			y2 = y[N];
			tmp = y[N];
			j = 0;
			// Loop for secant method
			while(j < maxCycNum && abs(E2 - E1) > resolution)
			{
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				PropKS(Eaux, l, y, VcentrVec, VextVec, VcoulVec, VxVec, Vc12Vec);

				// Condition to choose the next point
				if(y[N]*y1 > 0)
				{
					y1 = y[N];
					E1 = Eaux;
				}
				else if(y[N]*y1 < 0)
				{
					y2 = y[N];
					E2 = Eaux;
				}
				else if(y[N]*y1 == 0)
				{
					j = maxCycNum;
				}
				j++;
			}
			Esaved[0] += 2.0*(2.0*l+1.0)*Eaux;
			sol += 1;
			// Normalise the wavefunction and add the new density
			A = Amp(y);
			for(int m = 1; m < N + 1; m++)
			{
				y[m] = y[m]*A;
				rho[m] += (2.0*l + 1.0)/(2.0 * M_PI * m * m * h2)*y[m]*y[m]; // Update the density
			}
            Esaved[1] += (2.0*l+1)*FunT(y);
			Esaved[2] += (2.0*l + 1.0)*l*(l + 1.0)*FunVcen(y);
		} else {
		tmp = y[N];
		}
		E = E + DE;
	}
}

void Mix(double * rho, double * rhoPot, double * VcoulVec, double * VcoulVecOld)
{
    // I mix the new wave function with the previous one
    for(int i = 0; i < N + 1; i ++)
    {
        rhoPot[i] = alpha*rho[i] + (1.0 - alpha) *rhoPot[i];
        VcoulVec[i] = alpha*VcoulVec[i] + (1.0 - alpha) *VcoulVecOld[i];
    }
}

void PrintWf(double * y)
{
    for(int i = 0; i < N + 1; i ++)
    {
        ofile << fixed << setprecision(20) << y[i] << "\t";
    }
}

void Set_zero(double * rho)
{
    for(int i = 0; i < N + 1; i ++)
    {
        rho[i] = 0.0;
    }
}

int main()
{
    cout << "Mesh parameters" << endl;
    cout << "h:= " << h << " rmax:= " << rmax << " N:= " << N << " alpha:= " << alpha <<endl;

    cout << "Jellium parameters" << endl;
    cout << "Ne:= " << Ne << " rs:= " << rs << " rhoI:= " << rhoI << " Rc:= " << Rc << endl;

    // Parameters of Numerov
    int lmax;
    int * nmax = new int[4];
    if(Ne == 2)
    {
        lmax = 1;
        nmax[0] = 1;
        nmax[1] = 0;
        nmax[2] = 0;
        nmax[3] = 0;
    }
    if(Ne == 40)
    {
        lmax = 4;
        nmax[0] = 2;
        nmax[1] = 2;
        nmax[2] = 1;
        nmax[3] = 1;
    }

    int MaxCyc = 1e3;
    double Et = 1.0e-13;

    // Array for the important quantities
    double * y = new double[N + 1]; // Where I save the wave function
    double * rho = new double[N + 1]; // Where I save the density of the system
    double * rhoPot = new double[N + 1]; // Where I save the density used in the potentials
    double * Esaved = new double[3]; // Place where I save numerov eigenvalues, the kinetic energy and the centrifugal energy
    for(int i = 0; i < 3; i++)
    {
        Esaved[i] = 0;
    }
    double * VextVec = new double[N + 1];
    double * VcoulVec = new double[N + 1];
    double * VcoulVecOld = new double[N + 1];
    double * VxVec = new double[N + 1];
    double * Vc12Vec = new double[N + 1];
    for(int i = 0; i < N + 1; i++)
    {
        rho[i] = 0.0;
        rhoPot[i] = 0.0;
        VextVec[i] = Vext(i);
    }
    VcoulVec[0] = 0.0;
    rhoPot[0] = 0.0; // In the cycles I always do not consider the first point

    // Parameters for the functionals
    double Uext = 0.0;
    double Ucoul = 0.0;
    double Ux = 0.0;
    double UDex = 0.0;
    double Uc = 0.0;
    double UDec = 0.0;
    double E1 = 0.0;
    double E2 = 0.0;

    // First I calculate the independent-electron density
    cout << endl;
    cout << "Independent electrons" << endl;
    for(int l = 0; l < lmax; l++)
    {
        NumerovInd(l, nmax[l], y, rho,  VextVec, Esaved);
    }
    rho[0] = rho[1];
    VxVec[0] = Vx(rho[0]);
    Vc12Vec[0] = Vc12(rho[0]);
    for(int i = 0; i < N + 1; i++)
    {
        VcoulVec[i] = Vcoul(i, rho);
        VxVec[i] = Vx(rho[i]);
        Vc12Vec[i] = Vc12(rho[i]);
    }
    Uext = FunVext(rho, VextVec);
    Ucoul = FunVcoul(rho, VcoulVec);
    Ux = FunVx(rho, VxVec);
    UDex = FunVDex(rho, VxVec);
    Uc = FunVc(rho);
    UDec = FunVDec(rho);
    E1 = Esaved[1] + Esaved[2] + Ucoul + Uext + Ux + Uc;
    E2 = Esaved[0] - Ucoul - UDex - UDec;
    // Now I start the self consistent procedure
    // I initialize the density for the potential
    cout << endl;
    cout << "Kohm-Sham equation" << endl;


    minen = -4;
        maxen = 1;
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
        }

    int k = 0;
    while((k < MaxCyc) && (abs(E2 - E1) > Et))
    {
        Mix(rho, rhoPot, VcoulVec, VcoulVecOld);
        for(int i = 0; i < N + 1; i++)
        {
            VxVec[i] = Vx(rhoPot[i]);
            Vc12Vec[i] = Vc12(rhoPot[i]);
        }
        Set_zero(rho);
        for(int i = 0; i < 3; i++)
        {
            Esaved[i] = 0;
        }
        for(int l = 0; l < lmax; l++)
        {
            NumerovKS(l, nmax[l], y, rho,  VextVec, VcoulVec, VxVec, Vc12Vec, Esaved);
        }
        for(int i = 0; i < N + 1; i++)
        {
            VcoulVecOld[i] = Vcoul(i, rho);
            VxVec[i] = Vx(rho[i]);
            Vc12Vec[i] = Vc12(rho[i]);
        }

        Uext = FunVext(rho, VextVec);
        Ucoul = FunVcoul(rho, VcoulVecOld);
        Ux = FunVx(rho, VxVec);
        UDex = FunVDex(rho, VxVec);
        Uc = FunVc(rho);
        UDec = FunVDec(rho);

        E1 = Esaved[1] + Esaved[2] + Ucoul + Uext + Ux + Uc;
        E2 = Esaved[0] - Ucoul - UDex - UDec;
        cout << "k:= " << k << " Func diff:= " << abs(E1 - E2) << endl;
        k++;
    }
    cout << endl;
    cout << "Final functionals" << endl;
    cout << fixed << setprecision(20) << "Numerov:= " << Esaved[0] << endl;
    cout << fixed << setprecision(20) << "T:= " << Esaved[1] << endl;
    cout << fixed << setprecision(20) << "Ucen:= " << Esaved[2] << endl;
    cout << fixed << setprecision(20) << "Uext:= " << Uext << endl;
    cout << fixed << setprecision(20) << "Ucoul:= " << Ucoul << endl;
    cout << fixed << setprecision(20) << "Ux:= " << Ux << endl;
    cout << fixed << setprecision(20) << "UDex:= " << UDex << endl;
    cout << fixed << setprecision(20) << "Uc:= " << Uc << endl;
    cout << fixed << setprecision(20) << "UDec:= " << UDec << endl;
    cout << fixed << setprecision(20) << "E:= " << Esaved[1] + Esaved[2] + Uext + Ucoul + Ux + Uc << endl;
    cout << fixed << setprecision(20) << "Eepsilon:= " << Esaved[0] - Ucoul - UDex - UDec << endl;

    return 0;
}
