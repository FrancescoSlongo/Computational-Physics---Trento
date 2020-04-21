#include<iostream>
#include<cmath>

using namespace std;

const double X_MAX = 7.;
const double h = 0.001;
const int MESH_SIZE = ceil(X_MAX/h);
double *yold = new double[MESH_SIZE];
double *y = new double[MESH_SIZE];

inline double pot(int i, double Na)
{
    return (h*(i+1)*h*(i+1)/2. + Na *yold[i]*yold[i]);
}

double numerov(double E)
{
    double Na = 1.;
    double k_1 = E - pot(0, Na);
    double k0 = E - pot(1, Na);
    double k1;
    double norm= h*h*h+2*h*2*h*h;

    y[0]=h;
    y[1]=2.*h;


    for(int i=1;i<MESH_SIZE-1;i++)
    {
        k1 = E - pot(i+1, Na);

        y[i+1] = (y[i]*(2-5*h*h/6*2*k0)-y[i-1]*(1+h*h/12*2*k_1))/(1+h*h/12*2*k1);
        k_1 = k0;
        k0 = k1;
        norm += y[i+1]*y[i+1]*h;
    }
    //cout<<"norm "<<norm<<endl;
    for(int i=0;i<MESH_SIZE;i++)
    {
        y[i]/=sqrt(abs(norm));
    }

    return y[MESH_SIZE-1];
}

double findGS(double precision)
{
    double tmp1, tmp2;
    double x1, x2, xm, y1, y2, ym;
    int jj=0;
    double E;
    double E0 = 0.;
    double dE = .01;
    double Emax = 6.;
    double n;

    tmp1 = numerov(E0);
    for(n=E0+dE;n<Emax;n+=dE)
    {
        tmp2 =numerov(n);
        if(tmp1*tmp2<0)
        {
            x1 = n-dE;
            x2 = n;
            y1 = numerov(x1);
            y2 = numerov(x2);
            while(abs(x2-x1)>precision&jj<10000)
            {
                //xm = x2 - y2*(x1-x2)/(y1-y2);
                xm = (x1+x2)/2.;
                ym = numerov(xm);
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
                jj++;
            }
            E=xm;
            cout<<"E = ";
            printf("%5.20f  %d\n", E, jj);
            break;
        }
        else if(tmp1*tmp2==0.)
        {
            if(tmp1==0.) E = n-dE;
            if(tmp2==0.) E = n;
            cout<<"E = ";
            printf("%5.10f\n", E);
            break;
        }
        tmp1 = tmp2;
    }
    if(n>=Emax){cout<<"No Energy Found"<<endl;}
    return E;
}

void GPsolve(double alpha)
{
    findGS(1e-12);
    for(int i=0;i<1000;i++)
    {
        for(int j=0;j<MESH_SIZE;j++)
        {
            yold[j] = alpha*y[j]+(1-alpha)*yold[j];
        }

        findGS(1e-12);
    }
}

inline double FDdet(double E)
{
    double Na = 1.;
    double f_1 = 0.;
    double f0 = 1.;
    double f1;


    for(int i=0;i<MESH_SIZE;i++)
    {
        f1 = (2. + 2.*h*h*pot(i, Na)-2.*h*h*E)*f0 - f_1;
        //cout<<"f_1 "<<f_1<<" f0 "<<f0<<" f1 "<<f1<<" pot "<<pot(i, Na)<<endl;
        f_1 = f0;
        f0 = f1;
    }
    //cout<<"f1 "<<f1<<endl;
    return f1;
}

double findGS_FD(double precision)
{
    double tmp1, tmp2;
    double x1, x2, xm, y1, y2, ym;
    int jj=0;
    double E;
    double E0 = 0.;
    double dE = .01;
    double Emax = 6.;
    double n;
    double Na = 1.;

    tmp1 = FDdet(E0);
    for(n=E0+dE;n<Emax;n+=dE)
    {
        tmp2 =FDdet(n);
        if(tmp1*tmp2<0)
        {
            x1 = n-dE;
            x2 = n;
            y1 = FDdet(x1);
            y2 = FDdet(x2);
            while(abs(x2-x1)>precision&jj<10000)
            {
                xm = x2 - y2*(x1-x2)/(y1-y2);
                //xm = (x2+x1)/2.;
                ym = FDdet(xm);
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
                jj++;
            }
            E=xm;
            cout<<"E = ";
            printf("%5.20f  %d\n", E, jj);
            break;
        }
        else if(tmp1*tmp2==0.)
        {
            if(tmp1==0.) E = n-dE;
            if(tmp2==0.) E = n;
            cout<<"E = ";
            printf("%5.10f\n", E);
            break;
        }
        tmp1 = tmp2;
    }
    if(n>=Emax){cout<<"No Energy Found"<<endl;}

    //evaluate the WF
    y[0] = h;
    y[1] = y[0]*(2. + 2.*h*h*pot(1, Na)-2.*h*h*E);
    for(int i=2;i<MESH_SIZE;i++)
    {
       y[i]= y[i-1]*(2. + 2.*h*h*pot(i, Na)-2.*h*h*E)-y[i-2];
    }

    return E;
}

void GPsolve_FD(double alpha)
{
    findGS_FD(1e-12);
    for(int i=0;i<10;i++)
    {
        for(int j=0;j<MESH_SIZE;j++)
        {
            yold[j] = alpha*y[j]+(1-alpha)*yold[j];
        }

        findGS_FD(1e-12);
    }
}




int main()
{
    FILE *out;
    out = fopen("output.txt", "w+");
    double *tmp;
    double alpha = 0.1;

    for(int i=0;i<MESH_SIZE;i++)
    {
        //yold[i]=alpha*2./pow(M_PI,1./4.)*exp(-0.5*h*(i+1)*h*(i+1));
        yold[i]=0;
    }
    cout<<"h= "<<h<<" X_MAX = "<<X_MAX<<" alpha = "<<alpha<<endl;
    //cout<<numerov(2.5);

    //GPsolve_FD(alpha);
    findGS_FD(1e-17);

    for(int i=0;i<MESH_SIZE;i++)
    {
        fprintf(out, "%.12lf ", y[i]);
    }

    fclose(out);
    return 0;
}