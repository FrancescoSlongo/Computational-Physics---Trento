// Computational Physics 2020 - Project 4: DFT Theory
// Master course in Physics - Università di Trento

// Kohm Sham Equation - Jellium

#include <iostream>	// input and output
#include <fstream>	// Library for writing in a file (ofile etc.)
#include <cmath>		// math functions
#include <string>		// for using strings
#include <iomanip>	// for setprecision, setw

using namespace std;

ofstream ofile;

// Parameters of the mesh
const double h = pow(10.0, -3);         // Mesh step
const double h2 = h*h;                  // Useful parameter
const double rmax = 15.0;                // Max value of the mesh
const int N = int(ceil(rmax/h));    // Number of mesh points

// Parameters of Jellium
const int Ne = 20;                   // Number of electrons
const double rs = 3.93;             // Wigner-Seitz radius 3.93 - 4.86
const double rhoI = 3.0/(4.0 * M_PI * pow(rs,3));
const double Rc = rs * pow(Ne, 1.0/3.0);
const double Rc2 = pow(Rc, 2);
const double Rc3 = rs*rs*rs*Ne;

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

double Vh(double x) // Function that gives back the harmonic potential in one point
{
	return 0.5*x*x;
}

double Amp(double * y) // Given a wave function, it calculate the normalization amplitude
{
	double A;
	// I take the square of the wave function for each point
	A = y[0]*y[0]+y[N]*y[N];
	for(int m = 1; m < N; m++)
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

void Prop(double E,  int l, double * y) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	// I calculate k2
	double * k2 = new double[N + 1]; // Value of k2
	for(int i = 1; i < N + 1; i++)
	{
		k2[i] = 2.0*(E - Vext(i)) - double(l * (l +1.0))/ (i*i*h2);;
	}
	// Initial conditions
	y[0] = 0.0;
	y[1] = pow(h, l + 1);
	y[2] = pow(2.0*h, l+1);

	// Calculate next step of wf
	for(int i = 2; i < N; i++)
    {
        y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
    }
    delete k2;
}

void Numerov(int l, double E0,  double * y, int nmax, double * Esaved) // Method that implements the secant method with Numerov
{
	// False point method method
	// Variables used for the false point method

	// Parameters of the secant method in the Numerov algorithm
	double DE = 0.01;      // Step of E

	// Auxiliary variables useful for Numerov
	double E, E1, E2, Eaux;
	double y1, y2, ytmp, tmp;

	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-13;
	int maxCyc = 1e6;
	int sol = 0;

	// I set some of the variables
	E = E0;
	ytmp = 0.0;
	tmp = 0.0;
	int j; 	//Counter for the while cycle
	int m;
	double A; // Normalization amplitude

	// Loop for energies
	while(sol < nmax){
		// I use Numerov
		Prop(E, l, y);
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
			while(j < maxCyc && abs(E2 - E1) > resolution)
			{
				// Formula to calculate next Energy
				Eaux = E1 - y1 * (E2 - E1) / (y2 - y1);
				// I use Numerov
				Prop(Eaux, l, y);

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
					j = maxCyc;
				}
				j++;
			}
			sol += 1;
			// Normalise the wavefunction
			A = Amp(y);
			for(int m = 0; m < N + 1; m++)
			{
				y[m] = y[m]*A;
			}
			m = 0;
			while(Esaved[m] != 0)
            {
                m++;
            }
            Esaved[m] = Eaux;
		} else {
		tmp = y[N];
		}
		E = E + DE;
	}

}

void PrintWf(double * y)
{
    for(int i = 0; i < N + 1; i ++)
    {
        ofile << fixed << setprecision(20) << y[i] << "\t";
    }
}


int main()
{
    cout << "Mesh parameters" << endl;
    cout << "h:= " << h << " rmax:= " << rmax << " N:= " << N << endl;

    cout << "Jellium parameters" << endl;
    cout << "Ne:= " << Ne << " rs:= " << rs << " rhoI:= " << rhoI << " Rc:= " << Rc << endl;

    // Parameters of Numerov
    double E;
    double E0 = -5.0;
    int lmax = 2;
    int nmax = 2;
    string name1, name2;
    name1 = "Wf";

    // Place where to save the wave function and the energies
    double * Esaved = new double[lmax * nmax];
    int * lsaved = new int[lmax * nmax];
    int * nsaved = new int[lmax * nmax];

    for(int i = 0; i < lmax * nmax; i++)
    {
        Esaved[i] = 0;
        nsaved[i] = 0;
        lsaved[i] = 0;
    }


    cout << "Results" << endl;
    int j = 0;
    for(int l = 0; l < lmax; l++)
    {
        name2 = name1;
        name2.append(to_string(Ne));
        name2.append("_");
        name2.append(to_string(l));
        name2.append(".txt");
        ofile.open(name2);
        // I initialize the arrays for the wave functions
        double * y = new double[N + 1]; // Where I save the wavefunction
        j = 0;
        while(Esaved[j] != 0)
        {
            j++;
        }
        for(int m = 0; m < nmax; m++)
        {
            lsaved[j+m] = l;
            nsaved[j+m] = m + 1;
        }

        Numerov(l, E0, y, nmax, Esaved);
        PrintWf(y);
        ofile.close();
        delete y;
    }

    // Reorder the energies
    for(int i = 0; i < lmax * nmax; i++)
    {
        cout << "l:= " << lsaved[i] << " n:= " << nsaved[i] << " E:= " << fixed << setprecision(10) << Esaved[i] << endl;
    }
    return 0;
}
