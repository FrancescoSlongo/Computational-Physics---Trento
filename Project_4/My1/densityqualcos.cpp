#include<iostream>
#include<cmath>
#include<iomanip>
#include<vector>
using namespace std;

double X_MAX;
double h;
int MESH_SIZE;
double *y, *ddy, *rho, *Ucoul, *UcoulOld, *rhoOld;

int Ne;
double rs;
double Rc;
int n=0;

/*inline double pot(int i)
{
    double r2 = h*(i+1)*h*(i+1);
    double p1 = UcoulOld[i] - pow(3./M_PI/M_PI*rhoOld[i], 1./3.);
    if(rhoOld[i]>1e-20)
    {
        double rs = pow(3/4/M_PI/rhoOld[i], 1./3.);
        double gamma = -0.103756;
        double b1 = 0.56371;
        double b2 = 0.27358;
        p1 += gamma*(1.+7.*b1/6.*sqrt(rs)+4*b2/3.*rs)/(1+b1*sqrt(rs)+b2*rs)/(1+b1*sqrt(rs)+b2*rs);
    }
    if(h*(i+1)<Rc)
    {
        return 0.5*Ne/(Rc*Rc*Rc)*(r2 - 3.*Rc*Rc)+p1;
    }
    else
    {
        return -Ne/h/(i+1)+p1;
    }
}*/

inline double pot(int i)
{
    double r2 = h*(i+1)*h*(i+1);
    return 0.5*r2;
}

void ddf() // O(h^4) second derivative
{
    int i = 0;
    ddy[i] = 45*y[i]-154*y[i+1]+214*y[i+2]-156*y[i+3]+61*y[i+4]-10*y[i+5];
    ddy[i] /= (12*h*h);
    i++;
    ddy[i] = 45*y[i]-154*y[i+1]+214*y[i+2]-156*y[i+3]+61*y[i+4]-10*y[i+5];
    ddy[i] /= (12*h*h);
    for (i = 2; i < MESH_SIZE-2; i++)
    {
        ddy[i] = -y[i-2]+16*y[i-1]-30*y[i]+16*y[i+1]-y[i+2];
        ddy[i] /= (12*h*h);
    }
    i = MESH_SIZE-2;
    ddy[i] = 45*y[i]-154*y[i-1]+214*y[i-2]-156*y[i-3]+61*y[i-4]-10*y[i-5];
    ddy[i] /= (12*h*h);
    i++;
    ddy[i] = 45*y[i]-154*y[i-1]+214*y[i-2]-156*y[i-3]+61*y[i-4]-10*y[i-5];
    ddy[i] /= (12*h*h);
}

double numerov(double E, double l)
{
    double k_1 = 2*(E - pot(0))- l*(l+1)/h/h;
    double k0 = 2*(E - pot(1))- l*(l+1)/h/h/4.;
    double k1;

    y[0]=pow(h,l+1);
    y[1]=pow(2.*h,l+1);

    for(int i=1;i<MESH_SIZE-1;i++)
    {
        k1 = 2.*(E - pot(i+1)) - l*(l+1)/h/(i+2)/h/(i+2);

        y[i+1] = (y[i]*(2.-5.*h*h/6.*k0)-y[i-1]*(1.+h*h/12.*k_1))/(1.+h*h/12.*k1);
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
    while(abs(x2-x1)>precision && jj<10000)
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

double findE(double E0, double Emax, double dE, double precision, int lmax, vector<double> &Energ, vector<int> &L)
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

    while(n<Ne && E<Emax)
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
                    //cout<<E<<" "<<E+dE<<endl;
                    //cout<<"l= "<<i<<" E= "<<setprecision(12)<<tmp<<endl;
                    n+=2*(2*i+1);

                    /*for(int j=0;j<10;j++)
                    {
                        if(!Energ[j])
                        {
                            Energ[j]=tmp;
                            L[j]=i;
                            break;
                        }
                    }*/

                    Energ.push_back(tmp);
                    L.push_back(i);
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

void calculateWF(const vector<double> &E, const vector<int> &L)
{
    for(int i=0;i<MESH_SIZE;i++)
    {
        rho[i]=0;
    }

    for(int i=0;i<E.size();i++)
    {
        numerov(E[i], L[i]);

        double norm = y[0]*y[0];
        //cout<<"norm "<<norm<<endl;
        for(int j=1;j<MESH_SIZE-1;j++)
        {
            (j % 2) ? norm += 4*y[j]*y[j] : norm += 2*y[j]*y[j];
        }
        norm+=y[MESH_SIZE-1]*y[MESH_SIZE-1];
        norm*=h/3;

        for(int j=0;j<MESH_SIZE;j++)
        {
            y[j]/=sqrt(abs(norm));
            rho[j]+=2*(2*L[i]+1)*y[j]*y[j]/4/M_PI/h/h/(j+1)/(j+1);
        }

    }
}

void printWF(const vector<double> &E, const vector<int> &L)
{
    FILE *out;
    out = fopen("output.txt", "w+");

    for(int i=0;i<MESH_SIZE;i++)
    {
        rho[i]=0;
    }

    for(int i=0;i<E.size();i++)
    {
        numerov(E[i], L[i]);

        double norm = y[0]*y[0];
        //cout<<"norm "<<norm<<endl;
        for(int j=1;j<MESH_SIZE-1;j++)
        {
            (j % 2) ? norm += 4*y[j]*y[j] : norm += 2*y[j]*y[j];
        }
        norm+=y[MESH_SIZE-1]*y[MESH_SIZE-1];
        norm*=h/3;

        for(int j=0;j<MESH_SIZE;j++)
        {
            y[j]/=sqrt(abs(norm));
            rho[j]+=2*(2*L[i]+1)*y[j]*y[j]/4/M_PI/h/h/(i+1)/(i+1);

            fprintf(out, "%.20lf ", y[j]);

        }

        fprintf(out, "\n");
    }

    for(int i=0;i<MESH_SIZE;i++)
    {
        fprintf(out, "%.20lf ", rho[i]);
    }
    fclose(out);
}

void calcCoul(double *Ucoul)
{
    double int1, int2;
    /*for(int i=0;i<MESH_SIZE;i++)  CASI SPECIALI i=0,1, MS-1, MS  e pari
    {
        int1=rho[0]*h*h+4.*rho[1]*4.*h*h;

        for(int j=2;j<=i-1;j++)
        {
            (j % 2) ? int1 += 4*rho[j]*(j+1)*h*(j+1)*h : int1 += 2*rho[j]*(j+1)*h*(j+1)*h;
        }

        int1 += rho[i]*h*h*(i+1)*(i+1);
        int1 /= 3*(i+1);


        int2=rho[i+1]*h*(i+2)+4.*rho[i+2]*(i+3)*h;

        for(int j=i+3;j<=MESH_SIZE-1;j++)
        {
            (j % 2) ? int2 += 4*rho[j]*(j+1)*h : int2 += 2*rho[j]*(j+1)*h;
        }

        int2 += rho[MESH_SIZE]*h*(MESH_SIZE+1);
        int2 *= h/3;

        Ucoul[i]=4*M_PI*(int1+int2);
    }*/

    for(int i=0;i<MESH_SIZE;i++)
    {
        int1 = 0;
        for(int j=0;j<=i;j++)
        {
            int1 += h*rho[j]*h*(j+1)*h*(j+1);
        }
        int1/=h*(i+1);

        int2=0;
        for(int j=i+1;j<MESH_SIZE;j++)
        {
            int2 += h*rho[j]*h*(j+1);
        }

        Ucoul[i]=4*M_PI*(int1+int2);
    }
}

double functionals(const vector<double> &E,const vector<int> &L)
{
    //first calculate E in order to use numerov with old Ucoul for kinetic term
    //calculate kinetic and correction functionals
    double intKin=0;
    double intCorr=0;
    double intkk=0;
    double intcc=0;
    double norm=0;
    for(int i=0;i<E.size();i++)
    {
        numerov(E[i],L[i]);

        norm = y[0]*y[0];
        //cout<<"norm "<<norm<<endl;
        for(int j=1;j<MESH_SIZE-1;j++)
        {
            (j % 2) ? norm += 4*y[j]*y[j] : norm += 2*y[j]*y[j];
        }
        norm+=y[MESH_SIZE-1]*y[MESH_SIZE-1];
        norm*=h/3;

        for(int j=0;j<MESH_SIZE;j++)
        {
            y[j]/=sqrt(abs(norm));
        }

        ddf();

        intkk = y[0]*ddy[0];
        intcc = y[0]*y[0]/h/h;
        for(int j=1;j<MESH_SIZE-1;j++)
        {
            (j % 2) ? intkk += 4*y[j]*ddy[j] : intkk += 2*y[j]*ddy[j];
            (j % 2) ? intcc += 4*y[j]*y[j]/h/h/(j+1)/(j+1) : intcc += 2*y[j]*y[j]/h/h/(j+1)/(j+1);
        }
        intkk += y[MESH_SIZE-1]*ddy[MESH_SIZE-1];
        intcc += y[MESH_SIZE-1]*y[MESH_SIZE-1]/h/h/MESH_SIZE/MESH_SIZE;
        intkk *= h/3;
        intcc *= h/3;

        intKin -= (2*L[i]+1)*intkk;
        //cout<<intkk<<endl;
        intCorr += L[i]*(2*L[i]+1)*(L[i]+1)*intcc;
    }
    cout<<"Kinetic func= "<<setprecision(15)<<intKin<<endl;
    cout<<"Correct func= "<<setprecision(15)<<intCorr<<endl;

    //calculate direct coulomb functional
    double intCoul=Ucoul[0]*h*h*rho[0];
    for(int i=1;i<MESH_SIZE-1;i++)
    {
        (i % 2) ? intCoul += 4*Ucoul[i]*rho[i]*h*h*(i+1)*(i+1) : intCoul += 2*Ucoul[i]*rho[i]*h*h*(i+1)*(i+1);
    }
    intCoul += Ucoul[MESH_SIZE-1]*rho[MESH_SIZE-1]*h*h*MESH_SIZE*MESH_SIZE;
    intCoul *= h/3.*2.*M_PI;
    cout<<"Coulomb func= "<<setprecision(15)<<intCoul<<endl;

    //calculate external potential functional
    double intExt=Ne/2/Rc/Rc/Rc*(h*h-3.*Rc*Rc)*rho[0]*h*h; //vext =-Ne/h/(i+1) for r>R =-Ne/2/Rc/Rc/Rc*(h*(i+1)*h*(i+1)*-3.*Rc*Rc)
    for(int i=1;h*(i+1)<Rc;i++)
    {
        if(i % 2)
        {
            intExt+=2*Ne/Rc/Rc/Rc*(h*(i+1)*h*(i+1)-3.*Rc*Rc)*rho[i]*h*(i+1)*(i+1)*h;
        }
        else
        {
            intExt+=Ne/Rc/Rc/Rc*(h*(i+1)*h*(i+1)-3.*Rc*Rc)*rho[i]*h*(i+1)*(i+1)*h;
        }
    }
    for(int i=floor(Rc/h-1);i<MESH_SIZE-1;i++)
    {
        (i % 2) ? intExt-=4*Ne*rho[i]*h*(i+1) : intExt-=2*Ne*rho[i]*h*(i+1);
    }
    intExt -= Ne*rho[MESH_SIZE-1]*h*MESH_SIZE;
    intExt *= h/3*4*M_PI;
    cout<<"externa func= "<<setprecision(15)<<intExt<<endl;

    //calculate exchange functional
    double intEpsX = -3./4.*pow(3./M_PI/M_PI*rho[0], 1./3.)*rho[0]*h*h;
    for(int i=1;i<MESH_SIZE-1;i++)
    {
        if(i % 2)
        {
            intEpsX -= 3.*pow(3./M_PI/M_PI*rho[i], 1./3)*rho[i]*h*h*(i+1)*(i+1);
        }
        else
        {
            intEpsX -= 3./2*pow(3./M_PI/M_PI*rho[i], 1./3)*rho[i]*h*h*(i+1)*(i+1);
        }
    }
    intEpsX -= 3./4*pow(3./M_PI/M_PI*rho[MESH_SIZE-1], 1./3)*rho[MESH_SIZE-1]*h*h*MESH_SIZE*MESH_SIZE;
    intEpsX *= h/3*4*M_PI;
    cout<<"Exchang func= "<<setprecision(15)<<intEpsX<<endl;

    //calculate the correlation functional
    double gamma = -0.103756;
    double b1 = 0.56371;
    double b2 = 0.27358;
    double rs = pow(3./4/M_PI/rho[0],1./3.);
    double intEpsC = gamma/(1+b1*sqrt(rs)+b2*rs)*rho[0]*h*h;
    for(int i=1;i<MESH_SIZE-1;i++)
    {
        rs = pow(3./4./M_PI/rho[i], 1./3);
        (i % 2) ? intEpsC += 4.*gamma/(1+b1*sqrt(rs)+b2*rs)*rho[i]*h*h*(i+1)*(i+1) : intEpsC += 2*gamma/(1+b1*sqrt(rs)+b2*rs)*rho[i]*h*h*(i+1)*(i+1);
    }
    rs = pow(3./4/M_PI/rho[MESH_SIZE-1], 1./3.);
    intEpsC += gamma/(1+b1*sqrt(rs)+b2*rs)*rho[MESH_SIZE-1]*h*h*MESH_SIZE*MESH_SIZE;
    intEpsC *= h/3.*4*M_PI;
    cout<<"Correla func= "<<setprecision(15)<<intEpsC<<endl;

    double Efunc = intKin + intCorr + intExt + intCoul + intEpsX +intEpsC;
    cout<<"totEner func= "<<setprecision(15)<<Efunc<<endl;


    //now calculate the functional E_epsilon
    //sum the numerov eigenvalues with correct muliplicity
    double E_eps = 0;
    for(int i=0;i<E.size();i++)
    {
        E_eps += 2*(2*L[i]+1)*E[i];
    }
    cout<<"numerov func= "<<setprecision(15)<<E_eps<<endl;

    //calculate derivative of exchange functional
    double intDEpsX = - pow(rho[0], 4./3.)/4.*pow(3/M_PI/M_PI, 1./3.)*h*h;
    for(int i=1;i<MESH_SIZE-1;i++)
    {
        if(i % 2)
        {
            intDEpsX -= 4*pow(rho[i],4./3.)/4.*pow(3/M_PI/M_PI, 1./3.)*h*h*(i+1)*(i+1);
        }
        else
        {
            intDEpsX -= 2*pow(rho[i],4./3.)/4.*pow(3/M_PI/M_PI, 1./3.)*h*h*(i+1)*(i+1);
        }
    }
    intDEpsX -=pow(rho[MESH_SIZE-1],4./3.)/4.*pow(3/M_PI/M_PI, 1./3.)*h*h*MESH_SIZE*MESH_SIZE;
    intDEpsX *= h/3.*4*M_PI;
    cout<<"DeExcha func= "<<setprecision(15)<<intDEpsX<<endl;

    //calculate correlation functional
    rs = pow(3./4/M_PI/rho[0],1./3.);
    double intDEpsC = gamma/(1+b1*sqrt(rs)+b2*rs)*(b1/6.*sqrt(rs)+b2/3.*rs)/(1+b1*sqrt(rs)+b2*rs)*rho[0]*h*h;
    for(int i=1;i<MESH_SIZE-1;i++)
    {
        rs = pow(3./4./M_PI/rho[i], 1./3.);
        if(i%2)
        {
            intDEpsC += 4*gamma/(1+b1*sqrt(rs)+b2*rs)*(b1/6.*sqrt(rs)+b2/3.*rs)/(1+b1*sqrt(rs)+b2*rs)*rho[i]*h*h*(i+1)*(i+1);
        }
        else
        {
            intDEpsC += 2*gamma/(1+b1*sqrt(rs)+b2*rs)*(b1/6.*sqrt(rs)+b2/3.*rs)/(1+b1*sqrt(rs)+b2*rs)*rho[i]*h*h*(i+1)*(i+1);
        }
    }
    rs = pow(3./4./M_PI/rho[MESH_SIZE-1], 1./3.);
    intDEpsC += gamma/(1+b1*sqrt(rs)+b2*rs)*(b1/6.*sqrt(rs)+b2/3.*rs)/(1+b1*sqrt(rs)+b2*rs)*rho[MESH_SIZE-1]*h*h*MESH_SIZE*MESH_SIZE;
    intDEpsC *=h/3.*4*M_PI;
    cout<<"DeCorre func= "<<setprecision(15)<<intDEpsC<<endl;

    //sum up to find energy functional from eigenvalues of numerov
    E_eps -= intCoul + intDEpsX + intDEpsC;
    cout<<"TotNuEn func= "<<setprecision(15)<<E_eps<<endl;

    return Efunc-E_eps;
}

int main()
{
    X_MAX = 15.;
    h = 0.001;
    MESH_SIZE = 2*ceil(X_MAX/h/2)+1;
    y = new double[MESH_SIZE];
    ddy = new double[MESH_SIZE];
    rho = new double[MESH_SIZE];
    rhoOld = new double[MESH_SIZE];
    Ucoul = new double[MESH_SIZE];
    UcoulOld = new double[MESH_SIZE];
    double *tmp;
    double diff;

    for(int i=0;i<MESH_SIZE;i++)
    {
        Ucoul[i]=0;
        UcoulOld[i]=0;
        rho[i]=0;
        rhoOld[i]=0;
    }


    Ne = 40;
    rs = 3.93; //Na
    //rs = 4.86; //K
    Rc = rs*pow(Ne,1./3.);
    double alpha = 0.2;

    cout<<"h= "<<h<<" X_MAX= "<<X_MAX<<" MESH= "<<MESH_SIZE<<" rs= "<<rs<<" Ne= "<<Ne<<" Rc= "<<Rc<<" alpha= "<<alpha<<endl;

    vector<double> E;
    E.reserve(10);
    vector<int> L;
    L.reserve(10);

    //Self consistent procedure
    for(int i=0;i<1;i++)
    {
        //reset the energy levels
        E.clear();
        L.clear();

        //calculate new rho
        findE(-10, 10, 0.1, 1e-17, 4, E, L);
        calculateWF(E, L);

        //mixing
        //for(int j=0;j<MESH_SIZE;j++)
        //{
        //    rho[j]=alpha*rho[j]+(alpha-1)*rhoOld[j];
        //}

        //calculate functionals
        calcCoul(Ucoul);
        diff = functionals(E, L);
        cout<<"func diff= "<<diff<<endl;

        //subtitute the old stuff with the new stuff for a new step
        /*tmp = UcoulOld;
        UcoulOld = Ucoul;
        Ucoul = tmp;

        tmp = rhoOld;
        rhoOld = rho;
        rho = tmp;*/
    }


    //final energies and save WFs on file
    for(int i=0;i<E.size();i++)
    {
        cout<<"l= "<<L[i]<<" E= "<<setprecision(15)<<E[i]<<endl;
    }
    cout<<"total # electrons= "<<n<<endl;

    printWF(E, L);
}
