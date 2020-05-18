#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <string.h>  // For memset
#include <iomanip>	// for setprecision, setw

using namespace std;

// Parameters of the mesh
const double h = pow(10.0, -3);         // Mesh step
const double h2 = h*h;                  // Useful parameter
const double rmax = 15.0;                // Max value of the mesh
const int N = int(ceil(rmax/h));    // Number of mesh points


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
	V *= 4.0 * M_PI * h * h2/3.0;
	// I give back the potential
	return V;
}

int main()
{
    double * rho = new double[N + 1]; // Where I save the density of the system
    double * UculoA = new double[N + 1];
    double * Uculo = new double[N + 1];
    double err = 0.0;
    rho[0] = 1.0/(4.0*M_PI);
    UculoA[0] = 0.5;
    Uculo[0] = 0.5;
    for(int i = 1; i < N + 1; i++)
    {
        rho[i] = exp(-i*i*h2)/(4.0*M_PI);
        UculoA[i] = sqrt(M_PI) * erf(i*h)/(4.0*i*h);
    }
    cout << "Relative error" << endl;
    for(int i = 1; i < N + 1; i++)
    {
        Uculo[i] = Vcoul(i, rho);
        cout << abs(Uculo[i] - UculoA[i])/UculoA[i] << endl;
        err += abs(Uculo[i] - UculoA[i])/UculoA[i];
    }
    cout << fixed << setprecision(20) << err << endl;

    cout << "integral difference" << endl;
    cout << fixed << setprecision(20) << abs(FunVcoul(rho, Uculo) - sqrt(M_PI/2.0)/8.0) << endl;

    return 0;
}
