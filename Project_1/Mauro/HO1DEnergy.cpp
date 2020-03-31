#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

double pot(double x) // Potential for 1D (and 3D HO)
{
	return 0.5*x*x;
}

double kPrev2, kPrev1, kCurr, uPrev2, uPrev1, uCurr, r;

double numerovTot(double E, double r0, double h, int meshpts, bool disc) // A whole evolution with Numerov
{
	double norm = 0; // Normalization variable
	kPrev2 = 2 * (E - pot(r0));
	kPrev1 = 2 * (E - pot(r0 + h));
	
	if (disc) // Depending on even/odd we have different initial condition for the wf
	{
		uPrev2 = 2.;
		uPrev1 = (2. - 5./6.*h*h*kPrev2)/(1. + 1./12.*h*h*kPrev1);
	}
	else
	{
		uPrev2 = 0.;
		uPrev1 = h;
	}
	
	norm += uPrev2*uPrev2; // Cavalieri-Simpson integration
	norm += 4*uPrev1*uPrev1;
	
	r = r0 + 2*h;
	
	int i = 0;
	while (i < meshpts - 3) // Just evolution with Numerov
	{
		kCurr = 2 * (E - pot(r));
		uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
		(i % 2) ? norm += 4*uCurr*uCurr : norm += 2*uCurr*uCurr;
		
		uPrev2 = uPrev1;
		kPrev2 = kPrev1;
		uPrev1 = uCurr;
		kPrev1 = kCurr;
		r += h;
		i++;
	}
	
	kCurr = 2 * (E - pot(r));
	uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
	norm += uCurr*uCurr;
	
	norm *= h/3;
	
	uCurr /= sqrt(norm); // Only normalize last point because I don't really need anything else
	
	return uCurr;
}

int main()
{
	double rmax = 10;
	double hVec[20] = {0.1, 0.07, 4.655e-2, 3.2e-2, 2.167e-2, 1.478e-2, 1e-2, 6.9e-3, 4.5e-3, 3.2e-3, 2.2e-3, 1.5e-3, 1e-3, 6.94e-4, 4.7e-4, 3.23e-4, 2.2e-4, 1.5e-4, 1e-4, 7e-05};
	double h;
	double Emin0 = 0.01, Emin;
	double Ecurr, Enew, Eprev;
	double Estep = .1, Ethreshold = 1e-13;
	int maxCycles = 1e5; // To stop false position method after too many iterations (otherwise it can take a lot ot time to converge especially for small rmax)
	double ulo, uhi, utemp;
	double r0 = 0.; // Initial condition
	double* uvals[2];
	const int nhs = 20;
	const int nEs = 3;
	double Es[nhs][nEs]; // Values of h on rows, values of E on columns
	
	Ecurr = Emin;
	
	bool disc;
	cout << "Even/odd [1/0]? ";
	cin >> disc;
	
	FILE* fid = fopen("HO1DE_1e-13_h_rmax10_ciao.txt", "w");
	
	int i, j = 0, k, n;
	double meshPts;
	while (j < 10)
	{
		h = hVec[j];
		meshPts = ceil(rmax / h); // Get number of mesh points
		!((int)meshPts % 2) ? meshPts++ : true; // Se to odd number (for Cavalieri-Simpson)
		cout << meshPts << " points with h = " << h << endl;
		
		Emin = Emin0;
		Ecurr = Emin;
		k = 0;
		n = 0;
		while (n < nEs)
		{
			ulo = numerovTot(Ecurr, r0, h, meshPts, disc);
			Ecurr = Emin + Estep;
			uhi = numerovTot(Ecurr, r0, h, meshPts, disc);
			
			while (ulo*uhi > 0) // Start looking for interval in which y(rmax, E) changes sign
			{
				ulo = uhi;
				Ecurr += Estep;
				uhi = numerovTot(Ecurr, r0, h, meshPts, disc);
			}
			
			Eprev = Ecurr - Estep;
			Emin = Ecurr;
			
			i = 0;
			while (i < maxCycles && abs(Ecurr - Eprev) > Ethreshold) // False position method
			{	
				Enew = Ecurr - uhi*(Ecurr - Eprev)/(uhi - ulo);
				utemp = numerovTot(Enew, r0, h, meshPts, disc);
				if (utemp * ulo < 0)
				{
					uhi = utemp;
					Ecurr = Enew;
				}
				else if (utemp * uhi < 0)
				{
					ulo = utemp;
					Eprev = Enew;
				}
				else if (utemp - ulo == 0 || uhi - ulo == 0)
				{
					i = maxCycles;
				}
	
				i++;
				if (i == maxCycles) // Stop if surpasses max number of cycles
					cout << "Over" << endl;
			}
			
			/*while (i < maxCycles && abs(Ecurr - Eprev) > Ethreshold) // Normal secant method (converges less fast)
			{	
				Enew = Ecurr - uhi*(Ecurr - Eprev)/(uhi - ulo);
				ulo = uhi;
				Eprev = Ecurr;
				Ecurr = Enew;
				uhi = numerovTot(Ecurr, r0, h, meshPts, disc);
				i++;
				if (uhi - ulo == 0)
				{
					i = maxCycles;
				}
			}*/
		
			cout << setprecision(15) << Ecurr << endl; // Return and save energy, start from above new energy for a new research
			Es[j][k] = Ecurr;
			Ecurr = Emin;
			k++;
			n++;
		}
		j++;
	}
	
	for (j = 0; j < nhs; j++) // Save into file
	{
		fprintf(fid, "%d ", hVec[j]);
		for (k = 0; k < nEs; k++)
		{
			fprintf(fid, "%0.20lf ", Es[j][k]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
}

