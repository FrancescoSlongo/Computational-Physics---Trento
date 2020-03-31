#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

double pot(double x)
{
	return 0.5*x*x;
}

double kPrev2, kPrev1, kCurr, uPrev2, uPrev1, uCurr, r;

double numerovTot(double E, double r0, double h, double l, int meshpts) // Same as 1D, just pass angular momentum and have effective potential
{
	double norm = 0;
	kPrev2 = 2 * (E - pot(r0)) - l*(l+1)/r0/r0;
	kPrev1 = 2 * (E - pot(r0 + h)) - l*(l+1)/(r0+h)/(r0+h);
	
	uPrev2 = pow(r0, l+1); // Different initial conditions
	uPrev1 = pow(r0+h, l+1);
	
	norm += uPrev2*uPrev2;
	norm += 4*uPrev1*uPrev1;
	
	r = r0 + 2*h;
	
	int i = 0;
	while (i < meshpts - 3)
	{
		kCurr = 2 * (E - pot(r)) - l*(l+1)/r/r;
		uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
		(i % 2) ? norm += 4*uCurr*uCurr : norm += 2*uCurr*uCurr;
		
		uPrev2 = uPrev1;
		kPrev2 = kPrev1;
		uPrev1 = uCurr;
		kPrev1 = kCurr;
		r += h;
		i++;
	}
	
	kCurr = 2 * (E - pot(r)) - l*(l+1)/r/r;
	uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
	norm += uCurr*uCurr;
	
	norm *= h/3;
	
	uCurr /= sqrt(norm);
	
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
	int maxCycles = 1e5;
	double ulo, uhi, utemp;
	double r0 = 0.;
	double* uvals[2];
	const int nhs = 20;
	const int nEs = 3;
	double Es[nhs][nEs];
	
	Ecurr = Emin;
	
	double l;
	cout << "Ang. mom.: "; // Pass angular momentum
	cin >> l;
	
	FILE* fid = fopen("HO3DE_1e-13_h_rmax10_1.txt", "w");
	
	int i, j = 0, k, n;
	double meshPts;
	while (j < nhs)
	{
		h = hVec[j];
		meshPts = ceil(rmax / h);
		!((int)meshPts % 2) ? meshPts++ : true;
		cout << meshPts << " points with h = " << h << endl;
		r0 = h;
		
		Emin = Emin0;
		Ecurr = Emin;
		k = 0;
		n = 0;
		while (n < nEs)
		{
			ulo = numerovTot(Ecurr, r0, h, l, meshPts);
			Ecurr = Emin + Estep;
			uhi = numerovTot(Ecurr, r0, h, l, meshPts);
			
			while (ulo*uhi > 0)
			{
				ulo = uhi;
				Ecurr += Estep;
				uhi = numerovTot(Ecurr, r0, h, l, meshPts);
			}
			
			Eprev = Ecurr - Estep;
			Emin = Ecurr;
			
			i = 0;
			while (i < maxCycles && abs(Ecurr - Eprev) > Ethreshold)
			{	
				Enew = Ecurr - uhi*(Ecurr - Eprev)/(uhi - ulo);
				utemp = numerovTot(Enew, r0, h, l, meshPts);
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
				if (i == maxCycles)
					cout << "Over" << endl;
			}
		
			cout << setprecision(15) << Ecurr << endl;
			Es[j][k] = Ecurr;
			Ecurr = Emin;
			k++;
			n++;
		}
		j++;
	}
	
	for (j = 0; j < nhs; j++)
	{
		fprintf(fid, "%lf ", hVec[j]);
		for (k = 0; k < nEs; k++)
		{
			fprintf(fid, "%0.20lf ", Es[j][k]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
}

