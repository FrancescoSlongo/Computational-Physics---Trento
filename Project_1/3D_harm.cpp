// Computational Physics 2020 - Project 1: Numerov Algorithm
// Master course in Physics - Università di Trento

// 3D Schrodinger equation solver 
// We use the Numerov Algorithm

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <iomanip>	// for setprecision, setw

using namespace std;

ofstream ofile;

double Vh(double x) // Function that gives back the harmonic potential in one point
{
	return 0.5*x*x;
}

void Prop(int N, int l, double E, double h, double * y, double (*V) (double)) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	double h2; // Auxiliary variable to save computation time
	h2 = h*h;
	
	// I calculate k2
	double * k2 = new double[N + 1]; // Value of k2
	for(int i = 1; i < N+1; i++)
	{
		k2[i] = 2.0*(E - V(i *h)) - double(l * (l +1.0))/ (i*i*h2);
	}
	
	// Initial conditions
	y[0] = 0.0;
	y[1] = pow(h, l + 1);
	y[2] = pow(2.0 * h, l+1);
	
	// Apply Numerov Algorithm
	for(int i= 2; i < N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
	delete[] k2;
}

double Amp(int N, double * y) // Given a wave function, it calculate the normalization amplitude
{	
	double A;
	// I take the square of the wave function for each point
	for(int m = 0; m < N + 1; m++)
	{
		A += y[m]*y[m];
	}
	// I give back the amplitude
	return 1.0/sqrt(A);
}

void Numerov(int N, double h, int l, double E0, double Emax, double DeltaE, double (*V) (double)) // Method that implements Numerov
{
	// Secant method
	// Variables used for the secant method
	double * y = new double[N + 1]; // Value of the wave function at the point xi
	double E, E1, E2, Eaux; 
	double y1, y2, ytmp, tmp;
	
	E = E0;
	ytmp = 0;
	tmp = 0;
	
	int j; 	//Counter for the while cycle
	double A; // Normalization amplitude
	
	// Do you want to plot the eigenfunctios?
	char answ;
	cout << "Do you want to print the eigenfunctions (y/n)? " ;
	cin >> answ;
	if(answ == 'y')
	{
			ofile.open("Eig_3Dh.txt");
	}
	
	// Loop for energies
	while(E < Emax + DeltaE){
	// I use Numerov
	Prop(N, l, E, h, y, V);
		
	// Condition if the sign of the wavefunction changes at the extremity
	if(y[N]*tmp < 0)
	{		
		// I prepare the variables for the method and I save the variables to continue the while loop
		E1 = E - DeltaE;
		E2 = E;
		y1 = tmp;
		y2 = y[N];
		tmp = y[N];
		j = 0;
		
		// Loop for secant method
		while(j < 1000 && abs(E2 - E1) > 0.000001)
		{	
			// Formula to calculate next Energy
			Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
			// I use Numerov
			Prop(N, l, Eaux, h, y, V);
			
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
				E2 = Eaux;
				E1 = E2;
			}
		}
		// If I have to print something I renormalize the eigenfunction
		if(answ == 'y')
		{
			A = Amp(N, y);
			for(int m = 0; m < N+1; m++)
			{
				ofile << fixed<< setprecision(5) << y[m]*A << "\t";
			}
			ofile << endl;
		}
		cout <<  fixed << setprecision(8) << "Energy value is: " << Eaux << endl;
		} else {
			tmp = y[N];
		}
		E = E + DeltaE;
	}
	ofile.close();
	delete[] y;
}

int main(){
	// Variables about the discretization of the space
	double rmax = 7; // Maximum point for the check of the energy
	int N; // Number of mesh point
	int Nexp = 4; // Exponent for defining number of mesh points
	double h; // Length of Delta x
	N = pow(10, Nexp); // Defining number of mesh points
	h = rmax/double(N);
	
	// Energy settings for our problem
	double DeltaE = 0.5; // Difference in energy in cheching solutions
	double E0 = 0.3; // Initial value for checking energy
	double Emax = 5.4; // Maximum value of the energy scanned
	int l = 0; // Value of angular momentum
	
	// Application of the algorithm
	Numerov(N, h, l, E0, Emax, DeltaE, Vh);
	
	return 0;
}