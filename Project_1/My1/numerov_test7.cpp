#include<iostream>
#include<cmath>

using namespace std;

double *x;
double *y;
int MESH_SIZE;
double X_MIN = 0.4;
double X_MAX = 7;
double h = 1e-3;

inline double pot_lj(double x)
{
    double x6=pow(x,6);
    return 4*(1./x6/x6-1./x6); //hbar = 0.265243
}

double numerov_lj(double E, int l)
{
    int i;
    double k0, k1, k_1;

    k_1 = 2./0.265243/0.265243*( E - pot_lj(x[1])) - l*(l+1)/x[1]/x[1];
    k0 = 2./0.265243/0.265243*( E - pot_lj(x[2])) - l*(l+1)/x[2]/x[2];

    for(i=2;i<MESH_SIZE;i++)
    {
        k1 = 2./0.265243/0.265243*( E - pot_lj(x[i+1])) - l*(l+1)/x[i+1]/x[i+1];
        y[i+1] = (y[i]*(2-5*h*h/6*k0)-y[i-1]*(1+h*h/12*k_1))/(1+h*h/12*k1);
        //cout<<y[i+1]<<endl;
        k_1 = k0;
        k0 = k1;
    }

    return y[MESH_SIZE-1];
}

double bessel_j(double x, int l)
{
    double j_1,j0,j1;

    j_1=cos(x)/x;
    j0=sin(x)/x;

    if(l>0)
    {
        for(int i=0;i<l;i++)
        {
            j1=(2*i+1)/x*j0-j_1;
            j_1=j0;
            j0=j1;
        }
    }
    return j0;
}

double bessel_n(double x, int l)
{
    double j_1,j0,j1;

    j_1=sin(x)/x;
    j0=-cos(x)/x;

    if(l>0)
    {
        for(int i=0;i<l;i++)
        {
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
    int x1 = beg ;
    int x2 = beg + int(span/h);
    //cout<<"x1 = "<<x1<<" x2 = "<<x2<<endl;
    fprintf(out, "%10.12lf ", E);
    for(int l=0;l<15;l++)
    {
        numerov_lj(E, l);
        rat = x[x2]* y[x1]/ x[x1]/ y[x2];
        tand_n = (rat* bessel_j(x[x2]*sqrt(E*2/0.265243/0.265243), l) - bessel_j(x[x1]*sqrt(E*2/0.265243/0.265243),l));
        tand_d = (rat* bessel_n(x[x2]*sqrt(E*2/0.265243/0.265243), l) - bessel_n(x[x1]*sqrt(E*2/0.265243/0.265243),l));
        //cout<<"d[l= "<<l<<"] = "<<atan(tand)<<endl;
        fprintf(out, " %10.12lf", atan(tand_n/tand_d));
        stot += (2*l+1)* sin(atan(tand_n/tand_d)) * sin(atan(tand_n/tand_d));
    }
    stot = stot * 4 *3.14 /(E*2/0.265243/0.265243);
    //cout<<"total cross section = "<<stot<<endl;

    fprintf(out, " %10.12lf\n", stot);

    return stot;
}





int main()
{
    cout<<"X_MIN= "<<X_MIN<<" h= "<<h<<" X_MAX = "<<X_MAX<<endl;

    FILE *out;
    out = fopen("test_e.txt", "w+");
    double span = 0.5;

    MESH_SIZE = ceil((X_MAX-X_MIN)/h);
    x=new double[MESH_SIZE];
    y=new double[MESH_SIZE];
    int beg = ceil((5.5 - X_MIN)/h);


    for(int i=0;i<MESH_SIZE;i++)
    {
        x[i]= X_MIN +i*h ;
    }

    y[0] = 0.;
    y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
    y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

    for(double E=0;E<0.5932;E+=1e-3)
    {
        find_phaseshift(E, beg, span, out);
        //fprintf(out, "%10.12lf ", zz);
    }
}
