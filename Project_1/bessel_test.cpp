#include<iostream>
#include<cmath>

using namespace std;

double h = 1e-3;
double X_MIN = 0.;
double X_MAX = 15.;
const int MESH_SIZE = ceil( (X_MAX - X_MIN)/h);

//bessel function generator
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
	    //recursion formula for bessel functions
            j1=(2*i+1)/x*j0-j_1;
            j_1=j0;
            j0=j1;
        }
    }
    return j0;
}

//Neumann functions generator
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
	    //recursion formula for bessel functions
            j1=(2*i+1)/x*j0-j_1;
            j_1=j0;
            j0=j1;
        }
    }
    return j0;
}



int main()
{
    double analyt, calcul;
    double *x = new double[MESH_SIZE];
    double *y = new double[MESH_SIZE];
    FILE *out;
    out = fopen("test_b.txt","w+");

    cout<<"step size h= "<<h<<endl;
    
    //analitycal functions from wikipedia>spherical Bessel functions
    int i;
    for(i=0;i<MESH_SIZE;i++)
    {	
	//for each point save the recursively generated function, the analytical one and their difference for l=1,2,3
        x[i]= X_MIN +i*h ;
        fprintf(out, "%10.12lf ",x[i]);
        calcul = bessel_j(x[i], 1);
        analyt = sin(x[i])/x[i]/x[i] - cos(x[i])/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf ", calcul, analyt, calcul-analyt);
        calcul = bessel_n(x[i], 1);
        analyt = - cos(x[i])/x[i]/x[i]-sin(x[i])/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf ", calcul, analyt, calcul-analyt);
        calcul = bessel_j(x[i], 2);
        analyt = (3./x[i]/x[i]-1)*sin(x[i])/x[i] - 3*cos(x[i])/x[i]/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf ", calcul, analyt, calcul-analyt);
        calcul = bessel_n(x[i], 2);
        analyt = (-3./x[i]/x[i]+1)* cos(x[i])/x[i]-3*sin(x[i])/x[i]/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf ", calcul, analyt, calcul-analyt);
        calcul = bessel_j(x[i], 3);
        analyt = (15./pow(x[i],3)-6/x[i])*sin(x[i])/x[i] - (15/x[i]/x[i]-1)*cos(x[i])/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf ", calcul, analyt, calcul-analyt);
        calcul = bessel_n(x[i], 3);
        analyt = (-15/pow(x[i],3)+6/x[i])*cos(x[i])/x[i] -(15/x[i]/x[i]-1)* sin(x[i])/x[i];
        fprintf(out, "%10.12lf %10.12lf %10.12lf\n", calcul, analyt, calcul-analyt);
    }

    fclose(out);
}
