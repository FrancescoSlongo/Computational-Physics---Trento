// Computational Physics 2020 - Project 1: Numerov Algorithm
// Master course in Physics - Università di Trento

// Spherical Bessel and Neumann functions calculator

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <iomanip>		// for setprecision, setw

using namespace std;

ofstream ofile;

void Bessel(int N, double h, int l, double * y)
{
	// I define the vectors that will contain my functions
	double * j1 = new double[N + 1];
	double * j2 = new double[N + 1];
	double * j3 = new double[N + 1];
	double * x = new double[N + 1]; //Auxiliary variable
	
	// The function in zero si often divergent in 0, I initilize the value and keep it to 0
	j1[0] = 0.0;
	j2[0] = 0.0;
	j3[0] = 0.0;
	x[0] = 0.0;
	
	// I initialize the know functions
	for(int i = 1; i < N + 1; i++)
	{	
		x[i] = i*h;
		j1[i] = cos(x[i])/x[i];
		j2[i] = sin(x[i])/x[i];
	}		
	// Conditions in the case it is not needed to apply the recursive
	if(l == -1)
	{
		for(int i =0; i < N + 1; i++)
		{
			y[i] = j1[i];
		}
	} else if(l == 0)
	{
		for(int i =0; i < N + 1; i++)
		{
			y[i] = j2[i];
		}
	} else if(l > 0) // Case when we have to apply the recursive formula
	{
		int m = 0; // Auxiliary variable for while
		// Loop for the recursive formula
		while(m < l)
		{
			for(int i = 1; i < N + 1; i++)
			{	
				j3[i] = (2.0*m + 1.0) / x[i] * j2[i] - j1[i];
				j1[i] = j2[i];
				j2[i] = j3[i];
			}
			m++;
		}
		for(int i =0; i < N + 1; i++)
		{
			y[i] = j3[i];
		}
	}
	for(int m = 0; m < N+1; m++)
	{
		ofile << fixed<< setprecision(5) << y[m] << "\t";
	}
	ofile << endl;
	delete[] j1, j2, j3;
}

void Neuman(int N, double h, int l, double * y)
{
	// I define the vectors that will contain my functions
	double * n1 = new double[N + 1];
	double * n2 = new double[N + 1];
	double * n3 = new double[N + 1];
	double * x = new double[N + 1]; //Auxiliary variable
	
	// The function in zero si often divergent in 0, I initilize the value and keep it to 0
	n1[0] = 0.0;
	n2[0] = 0.0;
	n3[0] = 0.0;
	x[0] = 0.0;
	
	// I initialize the know functions
	for(int i = 1; i < N + 1; i++)
	{	
		x[i] = i*h;
		n1[i] = sin(x[i])/x[i];
		n2[i] = -cos(x[i])/x[i];
	}		
	// Conditions in the case it is not needed to apply the recursive
	if(l == -1)
	{
		for(int i =0; i < N + 1; i++)
		{
			y[i] = n1[i];
		}
	} else if(l == 0)
	{
		for(int i =0; i < N + 1; i++)
		{
			y[i] = n2[i];
		}
	} else if(l > 0) // Case when we have to apply the recursive formula
	{
		int m = 0; // Auxiliary variable for while
		// Loop for the recursive formula
		while(m < l)
		{
			for(int i = 1; i < N + 1; i++)
			{	
				n3[i] = (2.0*m + 1.0) / x[i] * n2[i] - n1[i];
				n1[i] = n2[i];
				n2[i] = n3[i];
			}
			m++;
		}
		for(int i =0; i < N + 1; i++)
		{
			y[i] = n3[i];
		}
	}
	for(int m = 0; m < N+1; m++)
	{
		ofile << fixed<< setprecision(5) << y[m] << "\t";
	}
	ofile << endl;
	delete[] n1, n2, n3;
}

int main()
{
	int l = 2; // Value of the angular momentum for the spherical bessel function
	int N = pow(10,4); // Number of points for the Bessel
	double h = 10.0/double(N);
	
	// I define the vectors that will contain my functions
	double * y = new double[N + 1];
	
	ofile.open("Spherical.txt");
	Bessel(N, h, l, y);
	Neuman(N, h, l, y);
	ofile.close();
}