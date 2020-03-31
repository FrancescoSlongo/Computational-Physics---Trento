#include<iostream>
#include<cmath>

using namespace std;

double *x;
double *y;
int MESH_SIZE;
double X_MIN = 0.;
double X_MAX = 15.;
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
    for(int l=0;l<7;l++)
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
    double E = 0.3; //fixed energy
    char tmp;

    cout<<"E= "<<E<<" h= "<<h<<" X_MAX = "<<X_MAX<<endl;
    cout<<"test: starting point / Rmax / distance between points [s/r/d]-> ";
    cin>>tmp;


    if(tmp=='s') //Tests the difference in the phaseshifts in function of the starting 0<x<1 of the numerov algorithm
    {
        FILE *out;
        out = fopen("test_s.txt", "w+");
        double span = 0.5;  //distance between x_1 and x_2 is fixed

        for(X_MIN=0;X_MIN<1;X_MIN+=h*10)
        {	
	    //set mash and position of x_1 and x_2 for each iteration
            MESH_SIZE = ceil((X_MAX-X_MIN)/h);
            x=new double[MESH_SIZE];
            y=new double[MESH_SIZE];
            int beg = ceil((11. - X_MIN)/h);


            for(int i=0;i<MESH_SIZE;i++)
            {
                x[i]= X_MIN +i*h ;
            }
	    //set initial conditions
            y[0] = 0.;
            y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
            y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

            find_phaseshift(E, beg, span, out);
        }
        fclose(out);
    }
    else if(tmp=='r') //Test the difference in the paseshifts in function of 4<Rmax<15 sigma
    {
        FILE *out;
        out = fopen("test_r.txt", "w+");
        X_MIN = 0.4; //initial point is fixed
        double span = 0.5; //distance between points is fixed

        MESH_SIZE = ceil((X_MAX-X_MIN)/h);
        x=new double[MESH_SIZE];
        y=new double[MESH_SIZE];

        for(int i=0;i<MESH_SIZE;i++)
        {
            x[i]= X_MIN +i*h ;
        }
	//initial conditions
        y[0] = 0.;
        y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
        y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

        for(int beg = ceil((4. - X_MIN)/h);beg<ceil((14.4 - X_MIN)/h);beg+=10)
        {
            find_phaseshift(E, beg, span, out);
        }
        fclose(out);
    }
    else if(tmp=='d') //Test the difference in the paseshifts in function of spacing between the two control points
    {
        FILE *out;
        out = fopen("test_d.txt", "w+");
        X_MIN = 0.4; //initial point is fixed
        int beg = ceil((11. - X_MIN)/h); //x_1 is fixed

        MESH_SIZE = ceil((X_MAX-X_MIN)/h);
        x=new double[MESH_SIZE];
        y=new double[MESH_SIZE];

        for(int i=0;i<MESH_SIZE;i++)
        {
            x[i]= X_MIN +i*h ;
        }
	//initial conditions
        y[0] = 0.;
        y[1] = exp(-2.*sqrt(2)/5./0.265243/pow(x[1],5));
        y[2] = exp(-2*sqrt(2)/5./0.265243/pow(x[2],5));

        for(double span = 0.1;span<3;span+=h)
        {
            find_phaseshift(E, beg, span, out);
        }
        fclose(out);
    }

}
