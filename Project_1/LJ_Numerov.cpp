// Computational Physics 2020 - Project 1: Numerov Algorithm
// Master course in Physics - Università di Trento

// 3D Schrodinger equation solver 
// We use the Numerov Algorithm to solve the Schrodinger for the LJ potential

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <iomanip>	// for setprecision, setw

using namespace std;

ofstream ofile;

double mu = 1.6536699e-27; // reduced mass in Kg
double eps = 5.9e-3*1.60218e-19; // J
double sigma = 3.18e-10; // m
double hbar = 0.26524240934; // hbar in natural units
double b5 = sqrt(8.0/25.0)/hbar; // Constant that enters in the inital conditions
double pi = 4*atan(1);

double VLJ(double x) // Function that gives back the Leanard-Johnes potential in one point
{
	double x6 = pow(x,6);
	return 4.0*(1.0/(x6*x6) - 1.0/x6);
}

void Prop(int N, int l, double E, double h, double * y, double (*V) (double)) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	double h2; // Auxiliary variable to save computation time
	h2 = h*h;
	
	// I calculate k2
	double * k2 = new double[N + 1]; // Value of k2
	for(int i = 1; i < N+1; i++)
	{
		k2[i] = 2.0/pow(hbar,2)*(E - V(i*h)) - double(l * (l +1.0))/ (i*i*h2);
	}
	
	// Initial conditions
	y[0] = 0.0;
	for(int i = 1; i < N/20+1; i++)
	{
		y[i] = exp(-b5/pow(i*h,5));
	}
	// Apply Numerov Algorithm
	for(int i= 501; i < N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
	delete[] k2;
}

double Bessel(int N, double h, int l, double x)
{
	// I define the vectors that will contain my functions
	double j1, j2, j3; // Auxiliary variables
	double y; //Output of the function 
	
	// I initialize the know functions
	j1 = cos(x)/x;
	j2 = sin(x)/x;
	// Conditions in the case it is not needed to apply the recursive
	if(l == -1)
	{
		y = j1;
	} else if(l == 0)
	{
		y = j2;
	} else if(l > 0) // Case when we have to apply the recursive formula
	{
		int m = 0; // Auxiliary variable for while
		// Loop for the recursive formula
		while(m < l)
		{	
			j3 = (2.0*m + 1.0) / x * j2 - j1;
			j1 = j2;
			j2 = j3;
			m++;
		}
		y = j3;
	}
	return y;
}

double Neumann(int N, double h, int l, double x)
{
	// I define the vectors that will contain my functions
	double n1, n2, n3; // Auxiliary variables
	double y; //Output of the function 
	
	// I initialize the know functions
	n1 = sin(x)/x;
	n2 = -cos(x)/x;
	// Conditions in the case it is not needed to apply the recursive
	if(l == -1)
	{
		y = n1;
	} else if(l == 0)
	{
		y = n2;
	} else if(l > 0) // Case when we have to apply the recursive formula
	{
		int m = 0; // Auxiliary variable for while
		// Loop for the recursive formula
		while(m < l)
		{	
			n3 = (2.0*m + 1.0) / x * n2 - n1;
			n1 = n2;
			n2 = n3;
			m++;
		}
		y = n3;
	}
	return y;
}

int main(){
	// Variables about the discretization of the space
	double rmax = 7; // Maximum point for the check of the energy
	int N; // Number of mesh point
	double h; // Length of Delta x
	N = pow(10, 4); // Defining number of mesh points
	h = rmax/double(N);
	
	cout << "rmax:= " << rmax << endl;
	cout << "Number of mesh points:= " << N << endl;
	
	// Energy settings for our problem
	double Emax = 0.6; // Maximum value of energy scanned
	double deltaE = 0.01; // Resolution of energy
	double E = 0;
	int lmax = 6; // Value of angular momentum
	cout << "E:= " << E << endl;
	cout << "Maximum value of the angular momentum to check:= " << lmax << endl;
	
	// I calculate the value of k from the energy
	double k = sqrt(2*E/pow(hbar,2));
	
	// Initialization of some variables
	double * y = new double[N + 1]; // Array where I save the wf
	// Points where I want to calculate the Bessel and the Neumann functions
	int n1 = N/10 *(rmax - 1);
	int n2 = N/10*(rmax - 0.5);
	
	double kappa; // Variable defined to calculate phase shift
	double deltal; // Phase shift
	double sigmatot = 0; // Total cross section
	
	ofile.open("Cross.txt");
	while(E < Emax)
	{	
		for(int l = 0; l < lmax+1; l++)
		{
			// Application of the algorithm for propagation of the wave function
			Prop(N, l, E, h, y, VLJ);
	
			// Print the wave function
			//ofile.open("LJ.txt");
			//for(int m = 0; m < N+1; m++)
			//{
			//	ofile << fixed<< setprecision(20) << y[m] << "\t";
			//}
			//ofile.close();
	
			// Now I calculate the phase shift
			kappa = y[n1]*n2/(y[n2]*n1);
	
			deltal = atan2(kappa*Bessel(N, h, l, k*h*n2) - Bessel(N, h, l, k*h*n1), kappa*Neumann(N, h, l, k*h*n2) - Neumann(N, h, l, k*h*n1));
			sigmatot += 4.0*pi/pow(k,2)*(2.0*l + 1.0)*pow(sin(deltal),2);
		}
		ofile << fixed << setprecision(5) << E << "\t" << sigmatot << endl;
		E += deltaE;
	}
	ofile.close();
	return 0;
}