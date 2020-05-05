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
		y[N+1] = (2.0 - 5.0/6.0*h2*k2[N]) / (1.0 + h2/12.0 *k2[N+1]) / 2.0;
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
			y[i] = y[2*N-i];
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

double Amp(int N, double h, double * y) // Given a wave function, it calculate the normalization amplitude using the Cavalieri-Simpson method
{	
	double A;
	// I take the square of the wave function for each point
	A = 0.0;
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

void Numerov(int N, double h, double E0, double Emax, double DeltaE, double * Esaved, int * s) // Method that implements the secant method with Numerov
{
	// False point method method
	// Variables used for the secant method
	double * y = new double[2*N + 1]; // Where I save the wave function
	double E, E1, E2, Eaux; 
	double y1, y2, ytmp, tmp;
	double * Eout = new double[5];
	int counter1 = 0; //Counter used to use symmetric or asymmetric boundary conditions
	int counter2 = 0;
	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-12; 
	int maxCyc = 1e5;
	
	E = E0;
	ytmp = 0;
	tmp = 0;
	
	int j; 	//Counter for the while cycle
	double A; // Normalization amplitude
	
	
	// Loop forN-1 energies
	while(E < Emax + DeltaE){
		// I use Numerov
		if(counter1 == 0)
		{
			Prop_as(N, E, 's', h, y);
		}
		else 
		{
			Prop_as(N, E, 'a', h, y);
		}
		
		// Condition if the sign of the wavefunction changes at the extremity
		if(y[2*N]*tmp < 0 && counter2 != 0)
		{		
			// I prepare the variables for the method and I save the variables to continue the while loop
			E1 = E - DeltaE;
			E2 = E;
			y1 = tmp;
			y2 = y[2*N];
			tmp = y[2*N];
			j = 0;
			// Loop for secant method
			while(j < maxCyc && abs(E2 - E1)/(E1 + E2) > resolution)
			{	
				
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				if(counter1 == 0)
				{
					Prop_as(N, Eaux, 's', h, y);
				}
				else 
				{
					Prop_as(N, Eaux, 'a', h, y);
				}
				
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
					j = maxCyc;
				}
				j++;
			}
			// If I have to print something I renormalize the eigenfunction
			//A = Amp(N, h, y);
			//for(int m = 0; m < 2*N+1; m++)
			//{
			//	y[m] = y[m]*A;
			//	ofile << fixed<< setprecision(20) << y[m] << "\t";
			//}
			//A = Amp(N, h, y);
			//cout << "A:=" << 1.0 - A << endl; 
			//ofile << endl;
			cout <<  fixed << setprecision(10) << "Energy value is: " << Eaux << endl;
			// I save the energies
			Esaved[s[0]] = Eaux;
			s[0] = s[0] + 1;
			counter1 = (counter1 + 1)%2;
			counter2 = 0;
		} else {
		tmp = y[2*N];
		counter2 = 1;
		}
		E = E + DeltaE;
	}
	delete[] y;
}

// Analytical wave functions for the 1D harmonic oscillator
void Exact1D(int N, int n, double h, double * y)
{
	double x; // Auxiliary variables
	if(n == 0) // E = 0.5
	{
		for(int i = 0; i < 2*N+1; i++)
		{	
			x = (i-N)*h;
			y[i] = 1.0/sqrt(sqrt(pi))*exp(-x*x/2);	
		}
	}
	if(n == 1) // E = 1.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1.0/sqrt(sqrt(pi))*sqrt(2.0)*(x)*exp(-x*x/2);	
		}
	}
	if(n == 2) // E = 2.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1.0/(sqrt(sqrt(pi))*sqrt(2.0))*(2*x*x - 1.0)*exp(-x*x/2);	
		}
	}
	if(n == 3) // E = 3.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = -1.0/(sqrt(sqrt(pi))*sqrt(3.0))*(2.0*x*x*x - 3.0*x)*exp(-x*x/2);	
		}
	}
	if(n == 4) // E = 4.5
	{
		for(int i = 0; i < 2*N + 1; i++)
		{
			x = (i-N)*h;
			y[i] = 1.0/(sqrt(sqrt(pi))*sqrt(24.0))*(4.0*x*x*x*x - 12.0*x*x + 3.0)*exp(-x*x/2);	
		}
	}
}

int main()
{		
	// Variables about the discretization of the space
	double rmax = 7.0; // Maximum point for the check of the energy
	int N[60] = {100,    113,    127,    143,    160,    180,    202,
          227,    256,    287,    323,    363,    408,    459,
          516,    580,    651,    732,    823,    925,   1040,
         1169,   1315,   1478,   1661,   1868,   2100,   2360,
         2653,   2983,   3353,   3770,   4238,   4764,   5356,
         6021,   6769,   7610,   8555,   9618,  10812,  12155,
        13665,  15362,  17270,  19415,  21827,  24538,  27586,
        31012,  34864,  39195,  44063,  49536,  55689,  62606,
        70382,  79124,  88952, 100000}; // Number of mesh point
	double h[60]; // Length of Delta x
	int MAX = 60; // Variable to change when I want to run in my PC or Mauro's one
	ofile.open("h.txt");
	
	for(int i = 0; i < MAX; i++)
	{
		h[i] = rmax/double(N[i]);
		ofile << fixed<< setprecision(20) << h[i] << "\t";
	}
	ofile.close();
	
	// Plotting analytical wf
	double A, B; // Normalization constant
	cout << "Hey! I am plotting the analitical wavefunctions!" << endl;
	string nameAnal1, nameAnal2;
	nameAnal1 = "Wf_Anal_";
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		//nameAnal2 = nameAnal1;
		//ofile.open(nameAnal2.append(to_string(i)).append(".txt"));
		double * y  = new double[2*N[i] + 1];
		for(int j = 0; j < 5; j++)
		{
			Exact1D(N[i], j, h[i], y);
			//A = Amp(N[i], h[i], y);
			//for(int m = 0; m < 2*N[i]+1; m++)
			//{
			//	ofile << fixed<< setprecision(20) << y[m] << "\t";
			//}
			//ofile << endl;
		}
		//ofile.close();
		delete[] y;
	}
	
	
	// Setting for Numerov but exact energies
	double Eexact[5] = {0.5, 1.5, 2.5, 3.5, 4.5}; // Exact values of the energy
	
	cout << "Hey! I am plotting the wavefunctions with Numerov but exact energies!" << endl;
	string nameExact1, nameExact2;
	nameExact1 =  "Wf_Exact_";
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		double * y  = new double[2*N[i] + 1]; // Array where I save Numerov
		//nameExact2 = nameExact1;
		//ofile.open(nameExact2.append(to_string(i)).append(".txt"));
		for(int j = 0; j < 5; j++)
		{
			if(j%2 == 0)
			{
				Prop_as(N[i], Eexact[j], 's', h[i], y);
			} 
			else 
			{
				Prop_as(N[i], Eexact[j], 'a', h[i], y);
			}
			//A = Amp(N[i], h[i], y);
			//for(int m = 0; m < 2*N[i]+1; m++)
			//{
				//y[m] = y[m]*A;
				//ofile << fixed<< setprecision(20) << y[m] << "\t";
			//}
			//A = Amp(N[i], h[i], y);
			//cout << "A:=" << 1.0 - A << endl; 
			ofile << endl;
		}
		ofile.close();
	}
	
	// Energy settings for Numerov
	double DeltaE = 0.01; // Difference in energy in cheching solutions
	double E0 = 0.105; // Initial value for checking energy
	double Emax = 4.6;  // Maximum value of the energy scanned
	
	// I keep saved the energy values that I got with Numerov
	double * Esaved = new double [5*MAX];
	for(int i = 0; i < 5*MAX; i++)
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
		//ofile.open(nameNumerov2.append(to_string(i)).append(".txt"));
		Numerov(N[i], h[i], E0, Emax, DeltaE, Esaved, s);
		//ofile.close();
	}
	
	// Now I calculate the error between the analytical function and the Numerov one
	cout << "Hey! I am plotting the errors between Numerov and the anlytical ones!" << endl;
	double err, errtmp;
	ofile.open("Wf_Err.txt");
	for(int i = 0; i < MAX; i++)
	{
		cout << "Index i :=" <<  i << endl;
		double * y  = new double[2*N[i] + 1]; // Array where I save Numerov
		double * z  = new double[2*N[i] + 1]; // Array where I save the anlytical wave functions
		for(int j = 0; j < 5; j++)
		{
			Exact1D(N[i], j, h[i], z);
			if(j%2 == 0)
			{
				Prop_as(N[i], Esaved[5*i+j], 's', h[i], y);
			} 
			else 
			{
				Prop_as(N[i], Esaved[5*i+j], 'a', h[i], y);
			}
			// I calculate both normalisation amplitudes
			A = Amp(N[i], h[i], y);
			B = Amp(N[i], h[i], z);
			// I calculate the error
			errtmp = y[0]*z[0]+y[2*N[i]]*z[2*N[i]];
			for(int m = 1; m < 2*N[i]; m++)
			{
				if(m%2 ==1)
				{
					errtmp += 4*y[m]*z[m];
				} else {
					errtmp += 2*y[m]*z[m];
				}
			}
			errtmp *= A*B*h[i]/3.0;
			err = 1.0 - pow(errtmp,2);
			cout << "N := " << N[i] << " Energy := " << Esaved[5*i+j] <<" err:=" << err << endl; 
			ofile << fixed << setprecision(15) << err << "\t";
		}
		ofile << endl;
		delete[] y,z;
	}
	ofile.close();
	
	return 0;
}

