#include<iostream>
#include<cmath>
#include<iomanip>
using namespace std;

double X_MAX;
double h;
int MESH_SIZE;
double *y, *rho;

int Ne;
double rs;
double Rc;
int n=0;

inline double pot(int i)
{
    double r2 = h*(i+1)*h*(i+1);
    if(h*(i+1)<Rc)
    {
        return 0.5*Ne/(Rc*Rc*Rc)*(r2 - 3.*Rc*Rc);
    }
    else
    {
        return -Ne/h/(i+1);
    }
}

/*inline double pot(int i)
{
    double r2 = h*(i+1)*h*(i+1);
    return 0.5*r2;
}*/


double numerov(double E, double l)
{
    double k_1 = 2*(E - pot(0))- l*(l+1)/h/h;
    double k0 = 2*(E - pot(1))- l*(l+1)/h/h/4;
    double k1;

    y[0]=h;
    y[1]=pow(2.*h,l);


    for(int i=1;i<MESH_SIZE-1;i++)
    {
        k1 = 2*(E - pot(i+1)) - l*(l+1)/h/(i+1)/h/(i+1);

        y[i+1] = (y[i]*(2-5*h*h/6*k0)-y[i-1]*(1+h*h/12*k_1))/(1+h*h/12*k1);
        k_1 = k0;
        k0 = k1;
    }

    return y[MESH_SIZE-1];
}

double bisezione(double E0, double E1, double precision, double l)
{
    double x1, x2, xm, y1, y2, ym;
    int jj=0;

    x1 = E0;
    x2 = E1;
    y1 = numerov(x1, l);
    y2 = numerov(x2, l);
    while(abs(x2-x1)>precision&jj<10000)
    {
        xm = x2 - y2*(x1-x2)/(y1-y2);
        //xm = (x1+x2)/2.;
        ym = numerov(xm, l);
        //cout<<"x1 "<<x1<<" x2 "<<x2<<" xm "<<xm;
        //cout<<"   y1 "<<y1<<" y2 "<<y2<<" ym "<<ym<<endl;
        if(ym*y1<0.)
        {
            x2=xm;
            y2=ym;
        }
        else if(ym*y2<0.)
        {
            x1 = xm;
            y1 = ym;
        }
        else if(ym==0.)
        {
            x2 = xm;
            x1 = x2;
        }
        else if(ym==y1 || ym==y2)
        {
            break;
        }
        jj++;
    }

    return xm;
}

double findE(double E0, double Emax, double dE, double precision, int lmax, double *Energ, int *L)
{
    int Sscan = 0;
    int n=0;
    int *Vscan = new int[lmax];
    double *x1 = new double[lmax];
    double *x2 = new double[lmax];
    double E=E0;
    double tmp;

    for(int l=0;l<=lmax;l++)
    {
        x1[l]=numerov(E, l);
        Vscan[l] = 0;
    }

    while(n<Ne&&E<Emax)
    {
        for(int l=0;l<=lmax;l++)
        {
            x2[l]=numerov(E+dE, l);
            if(x1[l]*x2[l]<=0)
            {
                Vscan[l]=1;
                Sscan+=1;
            }
        }

        if(Sscan==1)
        {
            for(int i=0;i<=lmax;i++)
            {
                if(Vscan[i]==1)
                {
                    tmp = bisezione(E, E+dE, precision, i);
                    cout<<"l= "<<i<<" E= "<<tmp<<endl;
                    n+=2*(2*i+1);

                    for(int j=0;j<10;j++)
                    {
                        if(!Energ[j])
                        {
                            Energ[j]=tmp;
                            L[j]=i;
                            break;
                        }
                    }
                }
            }
        }
        else if(Sscan>1)
        {
            /*for(int i=0;i<=lmax;i++)
            {
                cout<<Vscan[i];
            }
            cout<<endl;*/
            n+=findE(E, E+dE, dE/10., precision, lmax, Energ, L);
        }

        for(int l=0;l<=lmax;l++)
        {
            x1[l]=x2[l];
            Vscan[l]=0;
        }
        E=E+dE;
        Sscan=0;
    }
    //cout<<"holes-> "<<n-Ne<<endl;
    return n;
    delete[] Vscan;
    delete[] x1;
    delete[] x2;
}

void printWF(double *E, int *L)
{
    FILE *out;
    out = fopen("output.txt", "w+");

    for(int i=0;i<MESH_SIZE;i++)
    {
        rho[i]=0;
    }

    for(int i=0;i<10;i++)
    {
        if(E[i]!=0)
        {
            numerov(E[i], L[i]);

            double norm = y[0]*y[0]+4*y[1]*y[1];
            //cout<<"norm "<<norm<<endl;
            for(int j=2;j<MESH_SIZE-1;j++)
            {
                (j % 2) ? norm += 4*y[j]*y[j] : norm += 2*y[j]*y[j];
            }
            norm+=y[MESH_SIZE-1]*y[MESH_SIZE-1];
            norm*=h/3;

            for(int j=0;j<MESH_SIZE;j++)
            {
                y[j]/=sqrt(abs(norm));
                rho[j]+=2*(2*L[i]+1)*y[j]*y[j];

                fprintf(out, "%.20lf ", y[j]);

            }

            fprintf(out, "\n");

        }
    }

    for(int i=0;i<MESH_SIZE;i++)
    {
        fprintf(out, "%.20lf ", rho[i]);
    }
    fclose(out);
}


int main()
{
    X_MAX = 15.;
    h = 0.0001;
    MESH_SIZE = 2*ceil(X_MAX/h/2)+1;
    y = new double[MESH_SIZE];
    rho = new double[MESH_SIZE];


    Ne = 20;
    rs = 3.93; //Na
    //rs = 4.86; //K
    Rc = rs*pow(Ne,1./3.);

    cout<<"h= "<<h<<" X_MAX= "<<X_MAX<<" MESH= "<<MESH_SIZE<<" rs= "<<rs<<" Ne= "<<Ne<<" Rc= "<<Rc<<endl;

    double E[10];
    int L[10];

    for(int i=0;i<10;i++)
    {
        E[i]=0;
        L[i]=0;
    }

    int n=findE(-10, 0, 0.1, 1e-10, 4, E, L);
    cout<<"total # electrons= "<<n<<endl;

    printWF(E, L);
    /*int i=0;
    while(E[i]!=0)
    {
        cout<<E[i]<<" "<<L[i]<<endl;
        i++;
    }*/

}
