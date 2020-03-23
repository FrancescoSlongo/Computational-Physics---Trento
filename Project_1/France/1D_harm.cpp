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

// Algorithm that implements Numerov and propagates the wave function, starting from the center and supposing a symmetric or antisymmetric wave function 
void Prop_as(int N, double E, char c, double h, double * y) 
{
	double h2 = h*h;// Auxiliary variable to save computation time
	
	// I calculate k2
	double * k2 = new double[2*N + 1]; // Value of k2
	for(int i = 0; i < 2*N+1; i++)
	{
		k2[i] = 2.0*(E - Vh((i - N)*h));
	}
	
	// Initial conditions
	if(c == 's')
	{
		y[N] = 1.0;
		y[N+1] = 1.0-h2/2.0;
	} else if(c == 'a') 
	{
		y[N] = 0.0;
		y[N+1] = h;
	}
	
	// Apply Numerov Algorithm
	for(int i= N+1; i < 2*N; i++)
		{
			y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
		}
		
	// Copy the mirrored part
	if(c == 's')
	{
		for(int i= 0; i < N; i++)
		{
			y[i] = -y[2*N-i];
		}
	} else if(c == 'a')
	{
		for(int i= 0; i < N; i++)
		{
			y[i] = -y[2*N-i];
		}
	}
	delete[] k2;
}

double Amp(int N, double h, double * y) // Given a wave function, it calculate the normalization amplitude
{	
	double A;
	// I take the square of the wave function for each point
	A = 0;
	A = y[0]*y[0]+y[2*N]*y[2*N];
	for(int m = 1; m < 2*N; m++)
	{
		if(m%2 ==1)
		{
			A += 4*y[m]*y[m];
		} else {
			A += 2*y[m]*y[m];
		}
		
	}
	A *= h/3.0;
	// I give back the amplitude
	return 1.0/sqrt(A);
}


// Given a wave function, it calculate the normalization amplitude only on half of the wave function and then it multiplies for the analytical result
double AmpHalf(int N, double h, double * y)
{	
	double A;
	// I take the square of the wave function for each point
	A = 0;
	A = y[0]*y[0]+y[N]*y[N];
	for(int m = 1; m < N; m++)
	{
		if(m%2 ==1)
		{
			A += 4.0*y[m]*y[m];
		} else {
			A += 2.0*y[m]*y[m];
		}
		
	}
	A *= h/3.0;
	// I give back the amplitude
	return 1.0/sqrt(2.0*A);
}

void Numerov(int N, double h, double E0, double Emax, double DeltaE, double * Esaved, int * s) // Method that implements the secant method with Numerov
{
	// Secant method
	// Variables used for the secant method
	double * y = new double[2*N + 1]; // Value of the wave function at the point xi
	double E, E1, E2, Eaux; 
	double y1, y2, ytmp, tmp;
	double * Eout = new double[5];
	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-12; 
	int maxCyc = 1e5;
	
	E = E0;
	ytmp = 0;
	tmp = 0;
	
	int j; 	//Counter for the while cycle
	double A; // Normalization amplitude
	
	
	// Loop for energies
	while(E < Emax + DeltaE){
		// I use Numerov
		Prop_ext(N, E, h, y);
		
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
			while(j < maxCyc && abs(E2 - E1) > resolution)
			{	
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				Prop_ext(N, Eaux, h, y);
				
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
				j++;
			}
			// If I have to print something I renormalize the eigenfunction
			A = AmpHalf(N, h, y);
			for(int m = 0; m < 2*N+1; m++)
			{
				ofile << fixed<< setprecision(13) << y[m]*A << "\t";
			}
			ofile << endl;
			cout <<  fixed << setprecision(10) << "Energy value is: " << Eaux << endl;
			Esaved[s[0]] = Eaux;
			s[0] = s[0] + 1;
		} else {
		tmp = y[2*N];
		}
		E = E + DeltaE;
	}
	delete[] y;
}

void Exact1D(int N, int n, double h, double * y, char s)
{
	double x; // Auxiliary variables
	if(n == 0) // E = 0.5
	{
		for(int i = 0; i < 2*N+1; i++)
		{	
			x = (i-N)*h;
			y[i] = 1.0/sqrt(sqrt(pi))*exp(-x*x/2);	
			if (s == 'y')
			{
				ofile << fixed<< setprecision(13) << y[i] << "\t";
			}
		}
	}
	if(n == 1) // E = 1.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1.0/sqrt(sqrt(pi))*sqrt(2.0)*(x)*exp(-x*x/2);	
			if (s == 'y')
			{
				ofile << fixed<< setprecision(13) << y[i] << "\t";
			}
		}
	}
	if(n == 2) // E = 2.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1.0/(sqrt(sqrt(pi))*sqrt(2.0))*(2*x*x - 1.0)*exp(-x*x/2);	
			if (s == 'y')
			{
				ofile << fixed<< setprecision(13) << y[i] << "\t";
			}
		}
	}
	if(n == 3) // E = 3.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1.0/(sqrt(sqrt(pi))*sqrt(3.0))*(2.0*x*x*x - 3.0*x)*exp(-x*x/2);	
			if (s == 'y')
			{
				ofile << fixed<< setprecision(13) << y[i] << "\t";
			}
		}
	}
	if(n == 4) // E = 4.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1.0/(sqrt(sqrt(pi))*sqrt(24.0))*(4.0*x*x*x*x - 12.0*x*x + 3.0)*exp(-x*x/2);	
			if (s == 'y')
			{
				ofile << fixed<< setprecision(13) << y[i] << "\t";
			}
		}
	}
	if (s == 'y')
	{
		ofile << endl;
	}
}

int main()
{		
	// Variables about the discretization of the space
	double rmax = 7.0; // Maximum point for the check of the energy
	int N[10] = {101, 215, 465, 1001, 2155, 4641, 10001, 21545, 46417, 100001}; // Number of mesh point
	double h[10]; // Length of Delta x
	int MAX = 10; // Variable to change when I want to run in my PC or Mauro's one
	for(int i = 0; i < 10; i++)
	{
		h[i] = rmax/double(N[i]);	
	}
	
	// Plotting analytical wf
	cout << "Hey! I am plotting the analitical wavefunctions!" << endl;
	string nameAnal1, nameAnal2;
	nameAnal1 = "Wf_Anal_";
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		nameAnal2 = nameAnal1;
		ofile.open(nameAnal2.append(to_string(i)).append(".txt"));
		double * y  = new double[2*N[i] + 1];
		for(int j = 0; j < 5; j++)
		{
			Exact1D(N[i], j, h[i], y, 'y');
		}
		ofile.close();
		delete[] y;
	}
	
	// Energy settings for Numerov
	double DeltaE = 0.01; // Difference in energy in cheching solutions
	double E0 = 0.105; // Initial value for checking energy
	double Emax = 4.6;  // Maximum value of the energy scanned
	
	// I keep saved the energy values that I got with Numerov
	double * Esaved = new double [5*MAX];
	for(int i = 0; i < 15; i++)
	{
		Esaved[i] = 0.0;
	}
	// First I check the values that I get out from the Numerov algorithm, in order to make a comparison later
	string nameNumerov1, nameNumerov2;
	nameNumerov1 =  "Wf_Numerov_";
	cout << "Hey! I am doing Numerov! " << endl;
	int * s = new int[1]; // Counter to keep track where I save energies with Numerov
	s[0] = 0;
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		nameNumerov2 = nameNumerov1;
		ofile.open(nameNumerov2.append(to_string(i)).append(".txt"));
		Numerov(N[i], h[i], E0, Emax, DeltaE, Esaved, s);
		ofile.close();
	}

	// Setting for Numerov but exact energies
	double Eexact[5] = {0.5, 1.5, 2.5, 3.5, 4.5}; // Exact values of the energy
	double A; // Norm of Wf 
	
	cout << "Hey! I am plotting the wavefunctions with Numerov but exact energies!" << endl;
	string nameExact1, nameExact2;
	nameExact1 =  "Wf_Exact_";
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		double * y  = new double[2*N[i] + 1]; // Array where I save Numerov
		nameExact2 = nameExact1;
		ofile.open(nameExact2.append(to_string(i)).append(".txt"));
		for(int j = 0; j < 5; j++)
		{
			Prop_ext(N[i], Eexact[j], h[i], y);
			A = AmpHalf(N[i], h[i], y);
			for(int m = 0; m < 2*N[i]+1; m++)
			{
				ofile << fixed<< setprecision(13) << y[m]*A << "\t";
			}
			ofile << endl;
		}
		ofile.close();
	}
	
	
	// Now I calculate the error between the analytical wave functions and the results we get
	double * erry = new double [5*MAX];
	double * errz = new double [5*MAX];
	double tmpy, tmpz;
	double B;
	cout << "Hey! I am plotting the errors!" << endl;
	string nameYo1, nameYo2;
	nameYo1 =  "Wf_Yo_";
	for(int i = 0; i < 3; i++)
	{
		cout << "Index i :=" <<  i << endl;
		double * x = new double[2*N[i] + 1]; // Array where I save the exact solution
		double * y = new double[2*N[i] + 1]; // Array where I save Numerov
		double * z = new double[2*N[i] + 1]; // Array where I save Numerov with exact energies
		//nameYo2 = nameYo1;
		//ofile.open(nameYo2.append(to_string(i)).append(".txt"));
		for(int j = 0; j < 5; j++)
		{
			// I assign the exact solution to a vector
			Exact1D(N[i], j, h[i], x, 'n');
			// I assign the exact solution to Numerov
			Prop_ext(N[i], Esaved[5*i+j], h[i], y);
			// I assign the exact solution to Numerov with exact energies
			Prop_ext(N[i], Eexact[j], h[i], z);
			// I calculate the amplitudes
			A = AmpHalf(N[i], h[i], y);
			B = AmpHalf(N[i], h[i], z);
			
			//for(int m = 0; m < 2*N[i]+1; m++)
			//{
			//	ofile << fixed<< setprecision(13) << y[m]*A << "\t";
			//}
			//ofile << endl;
			
			// Now I calculate the error on the w.f. using Cavalieri-Simson
			tmpy = x[0]*y[0]+x[2*N[i]]*y[2*N[i]];
			tmpz = x[0]*z[0]+x[N[i]]*z[N[i]];
			erry[5*i+j] = 1.0;
			errz[5*i+j] = 1.0;
			for(int m = 1; m < N[i]; m++)
			{
				if(m%2 ==1)
				{
					tmpy += 4.0*x[m]*y[m];
					tmpz += 4.0*x[m]*z[m];
				} else {
					tmpy += 2.0*x[m]*y[m];
					tmpz += 2.0*x[m]*z[m];
				}
			}
			// For the Numerov I continue after the center
			for(int m = N[i]+1; m < 2*N[i]; m++)
			{
				if(m%2 ==1)
				{
					tmpy += 4.0*x[m]*y[m];
				} else {
					tmpy += 2.0*x[m]*y[m];
				}
			}
			tmpy *= A*h[i]/3.0;
			tmpz *= B*h[i]/3.0;
			erry[5*i+j] -= pow(tmpy,2);
			errz[5*i+j] -= 4.0*pow(tmpz,2);
		}
		ofile.close();
		delete[] x,y,z;
	}
	
	// I print the errors of Numerov on a txt
	cout << "Errors about Numerov" << endl;
	ofile.open("Err_Numerov.txt");
	for(int i = 0; i < MAX; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			cout << fixed<< setprecision(13) << erry[5*i+j] << "\t";
			ofile << fixed<< setprecision(13) << erry[5*i+j] << "\t";
		}
		cout << endl;
		ofile << endl;
	}
	ofile.close();
	
	// I print the data of Numerov with exact energies on a txt
	cout << "Errors about exact Numerov" << endl;
	ofile.open("Err_Exact.txt");
	for(int i = 0; i < MAX; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			cout << fixed<< setprecision(13) << errz[5*i+j] << "\t";
			ofile << fixed<< setprecision(13) << errz[5*i+j] << "\t";
		}
		cout << endl;
		ofile << endl;
	}
	ofile.close();
	
	return 0;
}

