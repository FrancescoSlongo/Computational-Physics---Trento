#include<iostream>
#include<cmath>

using namespace std;

double *x;
double *y;
int MESH_SIZE;
double X_MIN = 0.4;
double X_MAX = 15;
double h = 1e-3;

//LJ potential
inline double pot_lj(double x)
{
    double x6=pow(x,6);
    return 4*(1./x6/x6-1./x6); //hbar = 0.265243
}

//Numerov evolution for an LJ potential
double numerov_lj(double E, int l)
{
    int i;
    double k0, k1, k_1;

    k_1 = 2./0.265243/0.265243*( E - pot_lj(x[1])) - l*(l+1)/x[1]/x[1];
    k0 = 2./0.265243/0.265243*( E - pot_lj(x[2])) - l*(l+1)/x[2]/x[2];

    for(i=2;i<MESH_SIZE;i++)
    {
	//evaluate the k^2
        k1 = 2./0.265243/0.265243*( E - pot_lj(x[i+1])) - l*(l+1)/x[i+1]/x[i+1];
	//Numerov recursive formula
        y[i+1] = (y[i]*(2-5*h*h/6*k0)-y[i-1]*(1+h*h/12*k_1))/(1+h*h/12*k1);
        k_1 = k0;
        k0 = k1;
    }

    return y[MESH_SIZE-1];
}

//generator of Bessel functions
double bessel_j(double x, int l)
{
    double j_1,j0,j1;

    //initial conditions
    j_1=cos(x)/x;
    j0=sin(x)/x;

    if(l>0)
    {
        for(int i=0;i<l;i++)
        {
	    //recursive formula
            j1=(2*i+1)/x*j0-j_1;
            j_1=j0;
            j0=j1;
        }
    }
    return j0;
}

//Generator of Neumann functions
double bessel_n(double x, int l)
{
    double j_1,j0,j1;

    //initial conditions
    j_1=sin(x)/x;
    j0=-cos(x)/x;

    if(l>0)
    {
        for(int i=0;i<l;i++)
        {
	    //recursive formula
            j1=(2*i+1)/x*j0-j_1;
            j_1=j0;
            j0=j1;
        }
    }
    return j0;
}

double find_phaseshift(double E, int beg, double span, FILE *out)
{
    double rat, tand_n, tand_d;
    double stot=0;
	//beg is the position of X_1 and span is the distance between X_1 and X_2
    int x1 = beg ;
    int x2 = beg + int(span/h);

    fprintf(out, "%10.12lf %10.12lf ", (X_MIN+x1*h), (X_MIN + x2*h));
    for(int l=0;l<7;l++) //l might be set higher l=15 to see effect of truncating the partial wave sum on the total cross section
    {
        numerov_lj(E, l);
        rat = x[x2]* y[x1]/ x[x1]/ y[x2];
        tand_n = (rat* bessel_j(x[x2]*sqrt(E*2/0.265243/0.265243), l) - bessel_j(x[x1]*sqrt(E*2/0.265243/0.265243),l));
        tand_d = (rat* bessel_n(x[x2]*sqrt(E*2/0.265243/0.265243), l) - bessel_n(x[x1]*sqrt(E*2/0.265243/0.265243),l));
	//save the phase shfts
        fprintf(out, " %10.12lf", atan(tand_n/tand_d));
	//increment the cross section by the partial wave
        stot += (2*l+1)* sin(atan(tand_n/tand_d)) * sin(atan(tand_n/tand_d));
    }
    //divide by the prefactors to obtain total cross section 
    stot = stot * 4 *3.14 /(E*2/0.265243/0.265243);

    //also save cross section
    fprintf(out, " %10.12lf\n", stot);

    return stot;
}



int main()
{
    cout<<"X_MIN= "<<X_MIN<<" h= "<<h<<" X_MAX = "<<X_MAX<<endl;

    FILE *out;
    out = fopen("test_e.txt", "w+");
    //now the only free parameter is E so fix the others and generate mesh
    double span = 0.5;
    MESH_SIZE = ceil((X_MAX-X_MIN)/h);
    x=new double[MESH_SIZE];
    y=new double[MESH_SIZE];
    int beg = ceil((5.5 - X_MIN)/h);


    for(int i=0;i<MESH_SIZE;i++)
    {
        x[i]= X_MIN +i*h ;
    }
    //initial conditions
    y[0] = 0.;
    y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
    y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

    for(double E=0;E<0.5932;E+=1e-3) //test for each energy int the requested interval
    {
        find_phaseshift(E, beg, span, out);
    }
}
