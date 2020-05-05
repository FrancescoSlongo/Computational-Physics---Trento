// Computational Physics 2020 - Project 3: Mean Field methods
// Master course in Physics - Università di Trento

// 3D Gross-Pitaevsky solver

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
const double rmax = 7.0;                // Max value of the mesh
const int N = int(ceil(rmax/h));    // Number of mesh points

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

void Prop(double mu, double Na, double * yHar, double * y) // Algorithm that implements Numerov and propagates the wave function, starting from one extremity
{
	// I calculate k2
	double * k2 = new double[N + 1]; // Value of k2
	for(int i = 1; i < N + 1; i++)
	{
		k2[i] = 2.0*(mu  - Na * yHar[i]*yHar[i]/(i*i*h2)) - i*i*h2;
	}
	// Initial conditions
	y[0] = 0.0;
	y[1] = h;
	y[2] = 2*h;

	// Calculate next step of wf
	for(int i = 2; i < N; i++)
    {
        y[i+1] =(y[i]*(2.0 - 5.0/6.0 * h2 * k2[i]) - y[i-1]*(1.0 + h2 /12.0* k2[i-1])) / (1.0 + h2 / 12.0 * k2[i+1]);
    }
    delete k2;
}

double Numerov(double Na, double * yHar, double * y) // Method that implements the secant method with Numerov
{
	// False point method method
	// Variables used for the false point method

	// Parameters of the secant method in the Numerov algorithm
	double mu0 = 1.401;       // Starting value of mu
	double Dmu = 0.01;      // Step of mu

	// Auxiliary variables useful for Numerov
	double mu, mu1, mu2, muaux;
	double y1, y2, ytmp, tmp;

	// Required resolution for the secant method and maximum number of iterations
	double resolution = 1.0e-12;
	int maxCyc = 1e5;
	int sol = 0;

	// I set some of the variables
	mu = mu0;
	ytmp = 0.0;
	tmp = 0.0;
	int j; 	//Counter for the while cycle
	double A; // Normalization amplitude

	// Loop forN-1 energies
	while(sol == 0){
		// I use Numerov
		Prop(mu, Na, yHar, y);
		// Condition if the sign of the wavefunction changes at the extremity
		if(y[N]*tmp < 0)
		{
			// I prepare the variables for the method and I save the variables to continue the while loop
			mu1 = mu - Dmu;
			mu2 = mu;
			y1 = tmp;
			y2 = y[N];
			tmp = y[N];
			j = 0;
			// Loop for secant method
			while(j < maxCyc && abs(mu2 - mu1) > resolution)
			{
				// Formula to calculate next Energy
				muaux = mu1 - y1 * (mu2 - mu1) / (y2 - y1);
				// I use Numerov
				Prop(muaux, Na, yHar, y);

				// Condition to choose the next point
				if(y[N]*y1 > 0)
				{
					y1 = y[N];
					mu1 = muaux;
				}
				else if(y[N]*y1 < 0)
				{
					y2 = y[N];
					mu2 = muaux;
				}
				else if(y[N]*y1 == 0)
				{
					j = maxCyc;
				}
				j++;
			}
			sol = 1;
			// Normalise the wavefunction
			A = Amp(y);
			for(int m = 0; m < N + 1; m++)
			{
				y[m] = y[m]*A;
			}
		} else {
		tmp = y[N];
		}
		mu = mu + Dmu;
	}

	return muaux;
}

void Init_Guess(double alpha, double * y) // Function that initialise the guess wf
{
    // I initialize the guess wave function as the groundstate of an harmonic oscillator
    for(int i = 0; i < N + 1; i ++)
    {
        y[i] = 2.0*i*h*sqrt(alpha)/pow(M_PI,0.25)*exp(-i*i*h2 / 2.0);
    }
}

double Tkin(double * y)
{
    double T;
	// I calculate the V due to the harmonic potential
	T = 1.0 / (12.0 * h2) * ((45.0 * y[1] - 154.0 * y[2] + 214.0 * y[3] - 156.0 * y[4] + 61.0 * y[5] - 10.0 * y[6]) * y[1] + (45.0 * y[N - 1] - 154.0 * y[N - 2] + 214.0 * y[N - 3] - 156.0 * y[N - 4] + 61.0 * y[N - 5] - 10.0 * y[N - 6]) * y[N - 1]);
	for(int m = 2; m < N - 1; m++)
	{
		if(m%2 ==1)
		{
			T += 1.0 / (3.0 * h2) * (-30.0 * y[m] + 16.0 * y[m + 1] + 16.0 * y[m - 1] - y[m + 2] - y[m - 2]) * y[m];
		} else {
			T += 1.0 / (6.0 * h2) * (-30.0 * y[m] + 16.0 * y[m + 1] + 16.0 * y[m - 1] - y[m + 2] - y[m - 2]) * y[m];
		}
	}
	T *= -0.5*h/3.0;
	T += -0.5*h*(1.0 / (12.0 * h2) * ((45.0 * y[1] - 154.0 * y[2] + 214.0 * y[3] - 156.0 * y[4] + 61.0 * y[5] - 10.0 * y[6]) * y[1] + (45.0 * y[N - 1] - 154.0 * y[N - 2] + 214.0 * y[N - 3] - 156.0 * y[N - 4] + 61.0 * y[N - 5] - 10.0 * y[N - 6]) * y[N - 1]));; // I add another time the smallest derivative at the extremes
	// I give back the potential
	return T;
}

double Vext(double * y)
{
    double V;
	// I calculate the V due to the harmonic potential
	V = N * N * h2 *y[N]*y[N]; // The term y[0]*y[0] is multiplied by 0
	for(int m = 1; m < N; m++)
	{
		if(m%2 ==1)
		{
			V += 4.0*m*m*h2*y[m]*y[m];
		} else {
			V += 2.0*m*m*h2*y[m]*y[m];
		}
	}
	V *= 0.5*h/3.0;
	// I give back the potential
	return V;
}

double Vint(double Na, double * yHar, double * y)
{
    double V;
	// I calculate the V due to the interaction potential
	V = y[N] * y[N] * yHar[N] * yHar[N] / (N * N * h2);  // I do not count y[0] because the wf goes to zero in zero
	for(int m = 1; m < N; m++)
	{
		if(m%2 ==1)
		{
			V += 4.0 * pow(y[m], 4) / (m * m * h2);
		} else {
			V += 2.0 * pow(y[m], 4) / (m * m * h2);
		}
	}
	V *= 0.5*Na*h/3.0;
	// I give back the potential
	return V;
}

void Mix(double alpha, double * yHar, double * y)
{
    // I mix the new wave function with the previous one
    for(int i = 0; i < N + 1; i ++)
    {
        yHar[i] = sqrt(alpha*y[i]*y[i] + (1.0 - alpha) *yHar[i]*yHar[i]);
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

    // Gross-Pitaevsky
    double Na[3] = {0.01, 1, 100};     // Parameter in front of Harmonic potential

    // Parameter of the mixing
    double alpha[4] = {0.01, 0.1, 0.5, 1};

    // I initialize the arrays for the wave functions
    double * yHar = new double[N + 1]; // Guess wave function
    double * y = new double[N + 1]; // Where I save the wavefunction

    // Functionals of the wave function
    double mu, T, Ve, Vi, E1, E2;
    string name1, name2;
    name1 = "Wf";

    cout << "Results" << endl;
    for(int i = 0; i < 3; i++)
    {
        name2 = name1;
        name2.append(to_string(i));
        name2.append(".txt");
        ofile.open(name2);

        for(int j = 0; j < 4; j++)
        {
            Init_Guess(alpha[j], yHar);

            // Apply Numerov
            mu = Numerov(Na[i], yHar, y);
            cout << "Na:= " << Na[i] << " Alpha:= " << alpha[j] << " Mu:= " << mu << endl;

            PrintWf(y);

            // Calculate some functionals
            T = Tkin(y);
            Ve = Vext(y);
            Vi = Vint(Na[i], yHar, y);
            E1 = mu - Vi;
            E2 = T + Ve + Vi;

            cout << "T:= " << T << " Vext:= " << Ve << " Vin:= " << Vi<< endl;
            cout << "E1:= " << E1 << " E2:= " << E2  << " Diff:= " << abs(E1 - E2) << endl;
            cout << endl;

            // Mix the wave function
            Mix(alpha[j], yHar, y);

            // Apply Numerov
            mu = Numerov(Na[i], yHar, y);

            ofile << endl;
        }
        ofile.close();
    }
    return 0;
}
