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

double numerovTot(double E, double r0, double h, int meshpts, bool disc)
{
	double norm = 0;
	kPrev2 = 2 * (E - pot(r0));
	kPrev1 = 2 * (E - pot(r0 + h));
	
	if (disc)
	{
		uPrev2 = 2.;
		uPrev1 = (2. - 5./6.*h*h*kPrev2)/(1. + 1./12.*h*h*kPrev1);
	}
	else
	{
		uPrev2 = 0.;
		uPrev1 = h;
	}
	
	norm += uPrev2*uPrev2;
	norm += 4*uPrev1*uPrev1;
	
	r = r0 + 2*h;
	
	int i = 0;
	while (i < meshpts - 3)
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
	
	uCurr /= sqrt(norm);
	
	return uCurr;
}

int main()
{
	double rmax = 10;
	int meshVec[10] = {101, 215, 465, 1001, 2155, 4641, 10001, 21545, 46417, 100001};
	double h;
	double Emin0 = 0.01, Emin;
	double Emax = 10.;
	double Ecurr, Enew, Eprev;
	double Estep = .1, Ethreshold = 1e-14;
	int maxCycles = 1e5;
	double ulo, uhi, utemp;
	double r0 = 0.;
	double* uvals[2];
	double Es[10][6];
	
	Ecurr = Emin;
	
	bool disc;
	cout << "Even/odd [1/0]? ";
	cin >> disc;
	
	//FILE* fid = fopen("HO1DE_E_1e-12_10_nonorm_wiki.txt", "w");
	
	int i, j = 0, k;
	double meshPts;
	while (j < 10)
	{
		meshPts = meshVec[j];
		h = (rmax - r0)/meshPts;
		cout << meshPts << " points with h = " << h << endl;
		
		Emin = Emin0;
		Ecurr = Emin;
		k = 0;
		while (Ecurr < Emax)
		{
			ulo = numerovTot(Ecurr, r0, h, meshPts, disc);
			Ecurr = Emin + Estep;
			uhi = numerovTot(Ecurr, r0, h, meshPts, disc);
			
			while (ulo*uhi > 0)
			{
				ulo = uhi;
				Ecurr += Estep;
				uhi = numerovTot(Ecurr, r0, h, meshPts, disc);
			}
			
			Eprev = Ecurr - Estep;
			Emin = Ecurr;
			
			i = 0;
			while (i < maxCycles && abs(Ecurr - Eprev) > Ethreshold)
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
				if (i == maxCycles)
					cout << "Ghei" << endl;
			}
			cout << ulo + uhi << endl;
			
			/*while (i < maxCycles && abs(Ecurr - Eprev) > Ethreshold)
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
		
			cout << setprecision(15) << Ecurr << endl;
			Es[j][k] = Ecurr;
			Ecurr = Emin;
			k++;
		}
		j++;
	}
	
	/*for (j = 0; j < 10; j++)
	{
		fprintf(fid, "%d ", meshVec[j]);
		for (k = 0; k < 6; k++)
		{
			fprintf(fid, "%0.20lf ", Es[j][k]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);*/
}

