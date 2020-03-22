// Computational Physics 2020 - Project 1: Numerov Algorithm
// Master course in Physics - Università di Trento

// 1D Schrodinger equation solver 
// We use the Numerov Algorithm

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <iomanip>	// for setprecision, setw

using namespace std;

ofstream ofile;

double pi = 4.0*atan(1.0);

inline double Vh(double x) // Function that gives back the harmonic potential in one point
{
	return 0.5*x*x;
}

void Prop_ext(int N, double E, double h, double * y) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	double h2; // Auxiliary variable to save computation time
	h2 = h*h;
	
	// I calculate k2
	double * k2 = new double[2*N + 1]; // Value of k2
	for(int i = 0; i < 2*N+1; i++)
	{
		k2[i] = 2.0*(E - Vh((i - N)*h));
	}
	
	// Initial conditions
	y[0] = 0.0;
	y[1] = h;
	
	// Apply Numerov Algorithm
	for(int i= 1; i < 2*N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
	delete[] k2;
}

void Prop_sym(int N, double E, double h, double * y) // Algorithm that implements Numerov and propagates the wave function, starting from the center and supposing a symmetric wave function
{
	double h2;  // Auxiliary variable to save computation time
	h2 = h*h;
	
	// I calculate k2
	double * k2 = new double[2*N + 1]; // Value of k2
	for(int i = 0; i < 2*N+1; i++)
	{
		k2[i] = 2.0*(E - Vh((i - N)*h));
	}
	
	// Initial conditions
	y[N] = 1.0;
	y[N+1] = y[N]*(2.0 - 5.0/6.0 * h2 * k2[N]) / (2.0 + h2 / 6.0 * k2[N+1]);
	
	// Apply Numerov Algorithm
	for(int i= N+1; i < 2*N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
		
	// Copy the mirrored part
	for(int i= 0; i < N; i++)
		{
			y[i] = y[2*N-i];
		}
	delete[] k2;
}

void Prop_asym(int N, double E, double h, double * y) // Algorithm that implements Numerov and propagates the wave function, starting from the center and supposing a asymmetric wave function
{
	double h2;  // Auxiliary variable to save computation time
	h2 = h*h;
	
	// I calculate k2
	double * k2 = new double[2*N + 1]; // Value of k2
	for(int i = N; i < 2*N+1; i++)
	{
		k2[i] = 2.0*(E - Vh((i - N)*h));
	}
	
	// Initial conditions
	y[N] = 0.0;
	y[N+1] = h;
	
	// Apply Numerov Algorithm
	for(int i= N+1; i < 2*N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
		
	// Copy the mirrored part
	for(int i= 0; i < N; i++)
		{
			y[i] = -y[2*N-i];
		}
	delete[] k2;
}

double Amp(int N, double h, double * y) // Given a wave function, it calculate the normalization amplitude
{	
	double A;
	// I take the square of the wave function for each point
	A = 0;
	for(int m = 1; m < 2*N; m++)
	{
		A += y[m]*y[m];
	}
	A *= h;
	// I give back the amplitude
	return 1.0/sqrt(A);
}

void Numerov(int N, double h, double E0, double Emax, double DeltaE, void (*Prop) (int, double, double, double *)) // Method that implements the secant method with Numerov
{
	// Secant method
	// Variables used for the secant method
	double * y = new double[2*N + 1]; // Value of the wave function at the point xi
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
			ofile.open("Eig_1Dh.txt");
	}
	
	// Required resolution for the secant method
	double resolution;
	cout << "What resolution do you want for the secant method?";
	cin >> resolution;
	
	// Loop for energies
	while(E < Emax + DeltaE){
		// I use Numerov
		Prop(N, E, h, y);
		
		// Condition if the sign of the wavefunction changes at the extremity
		if(y[2*N]*tmp < 0)
		{		
			// I prepare the variables for the method and I save the variables to continue the while loop
			E1 = E - DeltaE;
			E2 = E;
			y1 = tmp;
			y2 = y[2*N];
			tmp = y[2*N];
			j = 0;
			
			// Loop for secant method
			while(j < 1000 && abs(E2 - E1) > resolution)
			{	
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				Prop(N, Eaux, h, y);
				
				// Condition to choose the next point
				if(y[2*N]*y1 > 0)
				{
					y1 = y[2*N];
					E1 = Eaux;
				}		
				else if(y[2*N]*y1 < 0) 
				{
					y2 = y[2*N];
					E2 = Eaux;
				}
				else if(y[2*N]*y1 == 0)
				{
					E2 = Eaux;
					E1 = E2;
				}
			}
			// If I have to print something I renormalize the eigenfunction
			if(answ == 'y')
			{
				A = Amp(N, h,y);
				for(int m = 0; m < 2*N+1; m++)
				{
					ofile << fixed<< setprecision(5) << y[m]*A << "\t";
				}
				ofile << endl;
			}
			cout <<  fixed << setprecision(8) << "Energy value is: " << Eaux << endl;
		} else {
		tmp = y[2*N];
		}
		E = E + DeltaE;
	}
	ofile.close();
	delete[] y;
}

void Exact1D(int N, int n, double h, double * y)
{
	double x = 0; // Auxiliary variable
	if(n == 0)
	{
		for(int i = 0; i < 2*N+1; i++)
		{	
			x = (i-N)*h;
			y[i] = 1.0/sqrt(sqrt(pi))*exp(-x*x/2);	
		}
	}
	if(n == 1)
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1/(sqrt(sqrt(pi))*sqrt(2.0))*(x)*exp(-x*x/2);	
		}
	}
	if(n == 2)
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1/(sqrt(sqrt(pi))*sqrt(2.0))*(2*x*x - 1.0)*exp(-x*x/2);	
		}
	}
	if(n == 3)
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1/(sqrt(sqrt(pi))*sqrt(3.0))*(2.0*x*x*x - 3.0*x)*exp(-x*x/2);	
		}
	}
	if(n == 4)
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1/(sqrt(sqrt(pi))*sqrt(16.0*6.0))*(4.0*x*x*x*x - 12.0*x*x + 3.0)*exp(-x*x/2);	
		}
	}
}

int main()
{		
	// Variables about the discretization of the space
	double rmax = 7.0; // Maximum point for the check of the energy
	int N ; // Number of mesh point
	int Nexp = 4; // Exponent for defining number of mesh points
	double h; // Length of Delta x
	
	// Energy settings for our problem
	double DeltaE = 0.1; // Difference in energy in cheching solutions
	double E0 = 0.1; // Initial value for checking energy
	double Emax = 0.2;  // Maximum value of the energy scanned

	// I apply the method that I want
	//int count;
	//cout << "What kind of method do you want to apply?" << endl;
	//cout << "Press 0 for the general one (starting from one extreme), 1 for only the symmetric eigenfunctions and 2 for the antisymmetric eigenfunctions.  ";
	//cin >> count;
	//if (count == 0)
	//{
	//	Numerov(N, h, E0, Emax, DeltaE, Prop_ext);
	//} else if (count == 1) 
	//{
	//	Numerov(N, h, E0, Emax, DeltaE, Prop_sym);
	//} else if (count == 2) 
	//{
	//	Numerov(N, h, E0, Emax, DeltaE, Prop_asym);
	//} else 
	//{
	//	cout << "You put the wrong input!!!";
	//}
	
	// Check for the errors due to the discretization and Numerov
	double * Eexact = new double[5]; // Exact values of the energy
	Eexact[0] = 0.5;
	Eexact[1] = 1.5;
	Eexact[2] = 2.5;
	Eexact[3] = 3.5;
	Eexact[4] = 4.5;
	
	int expmax = 5;   	// Maximum exponent for the 
	double A, B; 	// Amplitude for the normalization
	double err, errtmp; // Error between Numerov and exact wavefunction
	
	ofile.open("Errors_1D.txt");
	// Loop to check the error on the wavefunction
	for(int l = 0; l < 5; l++)
	{
		for(int s = 1; s <= expmax; s++)
		{
			N = pow(10,s);
			h = 7.0/double(N);
			double * y  = new double[2*N + 1]; // Array where I save Numerov
			double * z  = new double[2*N + 1]; // Array where I save the exact solution
			cout << "Doing Energy and Number of mesh points(exponent) equal to " << Eexact[l] << "  " << s << endl;
			// Numerov
			Prop_ext(N, Eexact[l], h, y);
			A = Amp(N,h, y);
			for(int m = 0; m < 2*N+1; m++)
			// Exact solution
			Exact1D(N, l, h, z);
			B = Amp(N, h, z);
			// I calculate the error
			err = 1.0;
			errtmp = 0.0;
			
			for(int m = 0; m < 2*N+1; m++)
			{
				errtmp += y[m]*z[m];
			}
			errtmp *= h*A*B;
			err -= pow(errtmp,2);
			ofile << fixed << setprecision(15) << err << "\t";
			cout << "err:= " << err << endl;
			delete[] y, z;
		}
		ofile << endl;
	}
	ofile.close();
	return 0;
}

