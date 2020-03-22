#include<iostream>
#include<cmath>

using namespace std;

double *x;
double *y;
int MESH_SIZE;
double X_MIN = 0.;
double X_MAX = 7;
double h = 1e-4;

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
    fprintf(out, "%10.12lf %10.12lf ", (X_MIN+x1*h), (X_MIN + x2*h));
    for(int l=0;l<7;l++)
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
    double E = 0.3;
    char tmp;

    cout<<"E= "<<E<<" h= "<<h<<" X_MAX = "<<X_MAX<<endl;
    cout<<"test: starting point / Rmax / distance between points [s/r/d]-> ";
    cin>>tmp;


    if(tmp=='s') //Tests the difference in the phaseshifts in function of the starting 0<x<1 of the numerov algorithm
    {
        FILE *out;
        out = fopen("test_s.txt", "w+");
        double span = 0.5;

        for(X_MIN=0;X_MIN<1;X_MIN+=h*10)
        {
            MESH_SIZE = ceil((X_MAX-X_MIN)/h);
            x=new double[MESH_SIZE];
            y=new double[MESH_SIZE];
            int beg = ceil((5. - X_MIN)/h);


            for(int i=0;i<MESH_SIZE;i++)
            {
                x[i]= X_MIN +i*h ;
            }

            y[0] = 0.;
            y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
            y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

            find_phaseshift(E, beg, span, out);
            //fprintf(out, "%10.12lf ", zz);
        }
        fclose(out);
    }
    else if(tmp=='r') //Test the difference in the paseshifts in function of 4<Rmax<6 sigma
    {
        FILE *out;
        out = fopen("test_r.txt", "w+");
        X_MIN = 0.4;
        double span = 0.5;

        MESH_SIZE = ceil((X_MAX-X_MIN)/h);
        x=new double[MESH_SIZE];
        y=new double[MESH_SIZE];

        for(int i=0;i<MESH_SIZE;i++)
        {
            x[i]= X_MIN +i*h ;
        }

        y[0] = 0.;
        y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
        y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

        for(int beg = ceil((4. - X_MIN)/h);beg<ceil((6.4 - X_MIN)/h);beg+=10)
        {
            find_phaseshift(E, beg, span, out);
            //fprintf(out, "%10.12lf ", zz);
        }
        fclose(out);
    }
    else if(tmp=='d') //Test the difference in the paseshifts in function of spacing between the two control points
    {
        FILE *out;
        out = fopen("test_d.txt", "w+");
        X_MIN = 0.4;
        int beg = ceil((5.5 - X_MIN)/h);

        MESH_SIZE = ceil((X_MAX-X_MIN)/h);
        x=new double[MESH_SIZE];
        y=new double[MESH_SIZE];

        for(int i=0;i<MESH_SIZE;i++)
        {
            x[i]= X_MIN +i*h ;
        }

        y[0] = 0.;
        y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
        y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

        for(double span = 0.9;span<1.5;span+=h)
        {
            find_phaseshift(E, beg, span, out);
            //fprintf(out, "%10.12lf ", zz);
        }
        fclose(out);
    }

}
