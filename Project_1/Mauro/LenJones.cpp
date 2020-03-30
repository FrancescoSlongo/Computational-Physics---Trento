#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

double epsilon = 0.0059;
double unitE = 1.60217662e-19;
double sigma = 3.18e-10;
double mu = 0.9958627415;
double unitU =  1.660539e-27;
double hbar = 1.0545718e-34;
double bconst;
double kPrev2, kPrev1, kCurr, uPrev2, uPrev1, uCurr, r;
FILE* fid = fopen("./LJwavefunc2.txt", "w");
double M_PI = 3.141592653589793;

double pot(double x)
{
	return 4./pow(x, 6.)*(1./pow(x, 6.) - 1.);
}

double psi0(double x)
{
	//cout << 1. / x << ", " << pow(1. / x, 5.) << ", " << bconst*pow(1. / x, 5.) << ", " << exp(-bconst*pow(1. / x, 5.)) << endl;
	return exp(-bconst*pow(1. / x, 5.));
}

double psi0Log(double x)
{
	return -bconst*pow(1. / x, 5.);
}

double bessel(double x, double l)
{
	double j1 = sin(x)/x, j0 = cos(x)/x, jm;
	for (int i = 0; i < l; i++)
	{
		jm = (2*i+1)/x*j1 - j0;
		j0 = j1;
		j1 = jm;
	}
	return j1;
}

double neumann(double x, double l)
{
	double j1 = -cos(x)/x, j0 = sin(x)/x, jm;
	for (int i = 0; i < l; i++)
	{
		jm = (2*i+1)/x*j1 - j0;
		j0 = j1;
		j1 = jm;
	}
	return j1;
}

void numerovLog(double E, double r0, double h, double cond, double l, double* vals, double* uvals)
{
	r = r0;
	kPrev2 = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/(r*r);
	uPrev2 = psi0Log(r);
	fprintf(fid, "%lf %0.12lf %0.12lf\n", r, exp(uPrev2), pot(r));
	
	r += h;
	kPrev1 = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/(r*r);
	uPrev1 = psi0Log(r);
	fprintf(fid, "%lf %0.12lf %0.12lf\n", r, exp(uPrev1), pot(r));
	
	uCurr = uPrev1;
	
	r += h;
	while (exp(uCurr) < cond)
	{
		kCurr = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/r/r;
		
		uCurr = uPrev1 + log(2. - 5*h*h/6.*kPrev1 - exp(uPrev2 - uPrev1)*(1. + h*h/12.*kPrev2)) - log(1. + h*h/12.*kCurr);
		fprintf(fid, "%lf %0.12lf %0.12lf\n", r, exp(uCurr), pot(r));
		vals[0] = r - h; // Zero is i-2 last
		uvals[0] = exp(uPrev1);
		vals[1] = r;
		uvals[1] = exp(uCurr);
		
		uPrev2 = uPrev1;
		kPrev2 = kPrev1;
		uPrev1 = uCurr;
		kPrev1 = kCurr;
		
		r += h;
	}
}

void numerovTot(double E, double r0, double h, double l, double rmax1, double rmax2, double* vals, double* uvals)
{
	r = r0;
	kPrev2 = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/(r*r);
	uPrev2 = psi0(r);
	//fprintf(fid, "%lf %0.12lf %0.12lf\n", r, uPrev2, pot(r));
	
	r += h;
	kPrev1 = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/(r*r);
	uPrev1 = psi0(r);
	//fprintf(fid, "%lf %0.12lf %0.12lf\n", r, uPrev1, pot(r));
	int i = 0;
	
	r += h;
	while (r < rmax1)
	{
		kCurr = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/r/r;
		
		uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
		//fprintf(fid, "%lf %0.12lf %0.12lf\n", r, uCurr, pot(r));
		//cout << uCurr << endl;
		
		uPrev2 = uPrev1;
		kPrev2 = kPrev1;
		uPrev1 = uCurr;
		kPrev1 = kCurr;
		r += h;
		i++;
	}
	
	vals[0] = r-h;
	uvals[0] = uCurr;
	
	while (r < rmax2)
	{
		kCurr = 2/(hbar*hbar)*(E - pot(r)) - l*(l+1)/r/r;
		
		uCurr = (uPrev1 * (2. - 5*h*h/6.*kPrev1) - uPrev2*(1. + h*h/12.*kPrev2))/(1. + h*h/12.*kCurr);
		//fprintf(fid, "%lf %0.12lf %0.12lf\n", r, uCurr, pot(r));
		//cout << uCurr << endl;
		
		uPrev2 = uPrev1;
		kPrev2 = kPrev1;
		uPrev1 = uCurr;
		kPrev1 = kCurr;
		r += h;
		i++;
	}
	
	vals[1] = r-h;
	uvals[1] = uCurr;
	
	//cout << "Steps: " << i << endl;;
}

int main()
{
	double l = 0.;
	cout << "Max ang. mom.: ";
	cin >> l;
	double h = 1e-6;
	cout << "Step(mom): ";
	cin >> h;
	double r0;
	double E;
	cout << "Energy: ";
	cin >> E;
	double rmax1 = 7., rmax2 = 8.;
	double vals[2], uvals[2];
	double k;
	double kappa;
	double phL;
	double totSigma = 0;
	
	hbar /= (sigma*sqrt(mu*epsilon*unitE*unitU));
	bconst = 2.*sqrt(2.)/hbar/5.;
	k = sqrt(2*E)/hbar;
	
	int j = 0;
	while (j <= l)
	{
		numerovTot(E, .35, h, j, rmax1, rmax2, vals, uvals);
		//cout << "(" << vals[0] << ", " << uvals[0] << "); (" << vals[1] << ", " << uvals[1] << ")" << endl;
		kappa = uvals[0]*vals[1]/(uvals[1]*vals[0]);
		phL = atan2(kappa*bessel(k*vals[1], j) - bessel(k*vals[0], j), kappa*neumann(k*vals[1], j) - neumann(k*vals[0], j));
		cout << "Phase shift for l = " << j << ": " << setprecision(10) << phL << endl;
		totSigma += (2.*j+1.)*sin(phL)*sin(phL);
		j++;
	}
	totSigma *= 4*M_PI/k/k;
	cout << "Total cross section: " << totSigma << endl;
	
	fclose(fid);
}

