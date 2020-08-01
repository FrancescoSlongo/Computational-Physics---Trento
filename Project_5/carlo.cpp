#include<iostream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<random>
#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
using namespace std;

/*
1| inizializzo le sigle particle WFs con loro parametri
2| inizializzo le posizioni x,y                                                                             <___
3| calcolo (H Y)/Y e suo quadrato                                                                               |
    3.1| calcolo matrice con derivate prime e seconde delle single particle WFs                                 |
    3.2| calcolo matrice inversa                                                                                |
    3.3| calcolo derivata determinante con formula matrice delle derivate * matrice inversa                     |   ripeto per stimare <E>
    3.?| se cambio una sola posizione devo ricalcolare una riga e basta o anche la matrice inversa?             |
4| calcolo nuova posizione con x(i+1)=x(i) + D*(rand()-.5)                                                      |
5| calcolo la nuova Y(x(i+1)) e valuto P=|Y(x(i+1)|^2/|Y(x(i))|^2 e scelgo se accettare o no con metropolis ____|

6| cambio parametri delle single particle WFs
7| cerco minimo di <E> dai parametri
*/

int dim1, dim2;
int *orba, *orbb;
gsl_matrix *A   ;
gsl_matrix *dxA ;
gsl_matrix *dyA ;
gsl_matrix *ddA ;
gsl_matrix *B   ;
gsl_matrix *dxB ;
gsl_matrix *dyB ;
gsl_matrix *ddB ;
double *xa, *xb, *ya, *yb;



double state(double x, double y, double alpha, int which){
    if(which==0){
        return exp(-0.5*alpha*(x*x+y*y));
    } else if(which==1){
        return x*exp(-0.5*alpha*(x*x+y*y));
    } else if(which==2){
        return y*exp(-0.5*alpha*(x*x+y*y));
    }
}

double dxstate(double x, double y, double alpha, int which){
    if(which==0){
        return -x*alpha*exp(-0.5*alpha*(x*x+y*y));
    } else if(which==1){
        return (1-alpha*x*x)*exp(-0.5*alpha*(x*x+y*y));
    } else if(which==2){
        return (-alpha*x*y)*exp(-0.5*alpha*(x*x+y*y));
    }
}

double dystate(double x, double y, double alpha, int which){
    if(which==0){
        return -y*alpha*exp(-0.5*alpha*(x*x+y*y));
    } else if(which==1){
        return (-alpha*x*y)*exp(-0.5*alpha*(x*x+y*y));
    } else if(which==2){
        return (1-alpha*y*y)*exp(-0.5*alpha*(x*x+y*y));
    }
}

double ddstate(double x, double y, double alpha, int which){
    if(which==0){
        return alpha*exp(-0.5*alpha*(x*x+y*y))*(alpha*x*x+alpha*y*y-2);
    } else if(which==1){
        return alpha*x*exp(-0.5*alpha*(x*x+y*y))*(alpha*x*x+alpha*y*y-4);
    } else if(which==2){
        return alpha*y*exp(-0.5*alpha*(x*x+y*y))*(alpha*x*x+alpha*y*y-4);
    }
}


void init_matrix(){

    for(int i=0;i<dim1;i++){
        for(int j=0;j<dim1;j++){
            gsl_matrix_set(A,   j, i,   state(xa[j], ya[j], 1, orba[i]));
            gsl_matrix_set(dxA, j, i, dxstate(xa[j], ya[j], 1, orba[i]));
            gsl_matrix_set(dyA, j, i, dystate(xa[j], ya[j], 1, orba[i]));
            gsl_matrix_set(ddA, j, i, ddstate(xa[j], ya[j], 1, orba[i]));
        }
    }
    //cout<<gsl_matrix_get(A,0,0)<<" "<<gsl_matrix_get(dxA,0,0)<<" "<<gsl_matrix_get(dyA,0,0)<<" "<<gsl_matrix_get(ddA,0,0)<<endl;

    for(int i=0;i<dim2;i++){
        for(int j=0;j<dim2;j++){
            gsl_matrix_set(B,   j, i,   state(xb[j], yb[j], 1, orbb[i]));
            gsl_matrix_set(dxB, j, i, dxstate(xb[j], yb[j], 1, orbb[i]));
            gsl_matrix_set(dyB, j, i, dystate(xb[j], yb[j], 1, orbb[i]));
            gsl_matrix_set(ddB, j, i, ddstate(xb[j], yb[j], 1, orbb[i]));
        }
    }
    //cout<<gsl_matrix_get(B,0,0)<<" "<<gsl_matrix_get(dxB,0,0)<<" "<<gsl_matrix_get(dyB,0,0)<<" "<<gsl_matrix_get(ddB,0,0)<<endl;
}

void calculateT(double *output, gsl_matrix *Maa, gsl_matrix *Mbb, gsl_permutation *pa, gsl_permutation *pb){
    gsl_matrix *Ma = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *Mb = gsl_matrix_alloc(dim2, dim2);

    gsl_matrix_memcpy(Ma, Maa);
    gsl_matrix_memcpy(Mb, Mbb);


    gsl_matrix *Ia = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *Ib = gsl_matrix_alloc(dim2, dim2);

    gsl_linalg_LU_invert(Ma, pa, Ia);
    gsl_linalg_LU_invert(Mb, pb, Ib);

    for(int i=0;i<dim1;i++){
        for(int j=0;j<dim1;j++){
            gsl_matrix_set(dxA, j, i, dxstate(xa[j], ya[j], 1, orba[i]));
            gsl_matrix_set(dyA, j, i, dystate(xa[j], ya[j], 1, orba[i]));
            gsl_matrix_set(ddA, j, i, ddstate(xa[j], ya[j], 1, orba[i]));
        }
    }

    for(int i=0;i<dim2;i++){
        for(int j=0;j<dim2;j++){
            gsl_matrix_set(dxB, j, i, dxstate(xb[j], yb[j], 1, orbb[i]));
            gsl_matrix_set(dyB, j, i, dystate(xb[j], yb[j], 1, orbb[i]));
            gsl_matrix_set(ddB, j, i, ddstate(xb[j], yb[j], 1, orbb[i]));
        }
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, dxA, Ia, 0, Ma);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, dxB, Ib, 0, Mb);

    double trdxa = 0, trdxb = 0, trdxa2 = 0, trdxb2 = 0;

    for(int i=0;i<dim1;i++){
        trdxa += gsl_matrix_get(Ma, i, i);
        trdxa2 += gsl_matrix_get(Ma, i, i) * gsl_matrix_get(Ma, i, i);
    }

    for(int i=0;i<dim2;i++){
        trdxb += gsl_matrix_get(Mb, i, i);
        trdxb2 += gsl_matrix_get(Mb, i, i) * gsl_matrix_get(Mb, i, i);
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, dyA, Ia, 0, Ma);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, dyB, Ib, 0, Mb);

    double trdya = 0, trdyb = 0, trdya2 = 0, trdyb2 = 0;

    for(int i=0;i<dim1;i++){
        trdya += gsl_matrix_get(Ma, i, i);
        trdya2 += gsl_matrix_get(Ma, i, i) * gsl_matrix_get(Ma, i, i);
    }

    for(int i=0;i<dim2;i++){
        trdyb += gsl_matrix_get(Mb, i, i);
        trdyb2 += gsl_matrix_get(Mb, i, i) * gsl_matrix_get(Mb, i, i);
    }

    gsl_blas_dgemm(CblasTrans, CblasTrans, 1, ddA, Ia, 0, Ma);
    gsl_blas_dgemm(CblasTrans, CblasTrans, 1, ddB, Ib, 0, Mb);

    double trdda = 0, trddb = 0;

    for(int i=0;i<dim1;i++){
        trdda += gsl_matrix_get(Ma, i, i);
    }

    for(int i=0;i<dim2;i++){
        trddb += gsl_matrix_get(Mb, i, i);
    }

    //cout<<trdxa<<" "<<trdya<<" "<<trdxb<<" "<<trdyb<<" "<<trdda<<" "<<trddb<<endl;

    output[0] = -0.5*(trdda + trddb /*+ 2*(trdxa*trdxb + trdya*trdyb)*/);
    //output[1] = 0.25*((trdxa*trdxa + trdya*trdya) + (trdxb*trdxb + trdyb*trdyb) - trdda - trddb);
    output[1] = 0.25*((trdxa2 + trdya2) + (trdxb2 + trdyb2) - trdda - trddb);
}

double pot(){
    double p=0;

    for(int i=0;i<dim1;i++){
        p += 0.5*(xa[i]*xa[i] + ya[i]*ya[i]);
    }
    for(int i=0;i<dim2;i++){
        p += 0.5*(xb[i]*xb[i] + yb[i]*yb[i]);
    }

    return p;
}

double monte_carlo_out(int steps, double delta){
    double out[3];
    int acc = 0;
    vector <double> T;
    T.reserve(steps);
    vector <double> Tjf;
    Tjf.reserve(steps);
    vector <double> V;
    V.reserve(steps);
    vector <double> X;
    X.reserve(steps+1);
    vector <double> Y;
    Y.reserve(steps+1);
    double Tcum = 0, Tjfcum = 0, Vcum = 0, Tcum2 = 0, Tjfcum2 = 0, Vcum2 = 0, Ecum2 = 0;

    X.push_back(xa[1]);
    Y.push_back(ya[1]);
    init_matrix();

    gsl_permutation *pa = gsl_permutation_alloc(dim1);
    gsl_permutation *pb = gsl_permutation_alloc(dim2);
    gsl_matrix *Ma = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *Mb = gsl_matrix_alloc(dim2, dim2);
    int sa, sb;

    gsl_matrix_memcpy(Ma, A);
    gsl_matrix_memcpy(Mb, B);
    //cout<<gsl_matrix_get(A,0,0)<<" "<<gsl_matrix_get(A,1,0)<<endl;
    //cout<<gsl_matrix_get(A,0,1)<<" "<<gsl_matrix_get(A,1,1)<<endl;
    gsl_linalg_LU_decomp(Ma, pa, &sa);
    //cout<<gsl_matrix_get(Ma,0,0)<<" "<<gsl_matrix_get(Ma,1,0)<<endl;
    //cout<<gsl_matrix_get(Ma,0,1)<<" "<<gsl_matrix_get(Ma,1,1)<<endl;
    gsl_linalg_LU_decomp(Mb, pb, &sb);

    double detai, detbi, detf;

    detai = gsl_linalg_LU_det(Ma, sa);
    detbi = gsl_linalg_LU_det(Mb, sb);

    gsl_matrix *tempa = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *tempb = gsl_matrix_alloc(dim2, dim2);
    gsl_permutation *temppa = gsl_permutation_alloc(dim1);
    gsl_permutation *temppb = gsl_permutation_alloc(dim2);

    int temps;
    double xtmp, ytmp, prob;

    for(int k=0;k<steps;k++){
        for(int i=0;i<dim1;i++){
            xtmp = xa[i] + delta*(rand()/(double)RAND_MAX -0.5);
            ytmp = ya[i] + delta*(rand()/(double)RAND_MAX -0.5);

            gsl_matrix_memcpy(tempa, A);
            for(int j=0;j<dim1;j++){
                gsl_matrix_set(tempa, i, j, state(xtmp, ytmp, 1, orba[j]));
            }

            gsl_matrix_memcpy(Ma, tempa);
            gsl_linalg_LU_decomp(Ma, temppa, &temps);

            detf = gsl_linalg_LU_det(Ma, temps);
            prob = detf*detf/detai/detai;
            //cout<<"detf= "<<detf<<" deta= "<<detai<<" detb= "<<detbi<<" prob= "<<prob<<endl;
            if(prob > rand()/(double)RAND_MAX){
                xa[i] = xtmp;
                ya[i] = ytmp;
                gsl_matrix_memcpy(A, tempa);
                gsl_permutation_memcpy(pa, temppa);
                detai = detf;

                acc+=1;
            }
        }
        gsl_matrix_memcpy(Ma, A);
        gsl_linalg_LU_decomp(Ma, pa, &sa);

        for(int i=0;i<dim2;i++){
            xtmp = xb[i] + delta*(rand()/(double)RAND_MAX -0.5);
            ytmp = yb[i] + delta*(rand()/(double)RAND_MAX -0.5);

            gsl_matrix_memcpy(tempb, B);
            for(int j=0;j<dim2;j++){
                gsl_matrix_set(tempb, i, j, state(xtmp, ytmp, 1, orbb[j]));
            }
            gsl_matrix_memcpy(Mb, tempb);
            gsl_linalg_LU_decomp(Mb, temppb, &temps);

            detf = gsl_linalg_LU_det(Mb, temps);
            prob = detf*detf/detbi/detbi;
            //cout<<"detf= "<<detf<<" deta= "<<detai<<" detb= "<<detbi<<" prob= "<<prob<<endl;
            if(prob > rand()/(double)RAND_MAX){
                xb[i] = xtmp;
                yb[i] = ytmp;
                gsl_matrix_memcpy(B, tempb);
                gsl_permutation_memcpy(pb, temppb);
                detbi = detf;
                //cout<<"A "<<gsl_matrix_get(A,0,0)<<" "<<gsl_matrix_get(A,1,0)<<" "<<gsl_matrix_get(A,0,1)<<" "<<gsl_matrix_get(A,1,1)<<endl;
                //cout<<"B "<<gsl_matrix_get(tempb,0,0)<<" "<<gsl_matrix_get(tempb,1,0)<<" "<<gsl_matrix_get(tempb,0,1)<<" "<<gsl_matrix_get(tempb,1,1)<<endl;

                acc+=1;
            }
        }
        gsl_matrix_memcpy(Mb, B);
        gsl_linalg_LU_decomp(Mb, pb, &sb);

        calculateT(out, Ma, Mb, pa, pb);

        T.push_back(out[0]);
        Tjf.push_back(out[1]);
        V.push_back(pot());
        Tcum += out[0];
        Tjfcum += out[1];
        Vcum += V.back();
        Tcum2 += out[0]*out[0];
        Tjfcum2 += out[1]*out[1];
        Vcum2 += V.back()*V.back();
        Ecum2 += (out[0] + V.back())*(out[0]+V.back());
        X.push_back(xa[1]);
        Y.push_back(ya[1]);
    }
    double sigma1, sigma2;
    cout<<"acceptance rate= "<<(double)acc/steps/(dim1+dim2)<<endl;
    sigma1 = sqrt((Tcum2/steps-Tcum*Tcum/steps/steps)/(steps-1));
    sigma2 = sqrt((Tjfcum2/steps-Tjfcum*Tjfcum/steps/steps)/(steps-1));
    cout<<"T= "<<Tcum/steps<<" +\\- "<< sigma1<< " Tjf= "<<Tjfcum/steps<<" +\\- "<< sigma2<<endl;
    sigma1 = sqrt((Vcum2/steps-Vcum*Vcum/steps/steps)/(steps-1));
    sigma2 = sqrt((Ecum2/steps-(Vcum+Tcum)/steps*(Vcum+Tcum)/steps)/(steps+1));
    cout<<"V= "<<Vcum/steps<<" +\\- "<<sigma1<<" E= "<<setprecision(20)<<(Vcum+Tcum)/steps<<" +\\- "<<sigma2<<endl;

    FILE *oof;
    oof = fopen("output.txt", "w+");
    for(int i=0;i<T.size(); i++){
        fprintf(oof, "%.20lf %.20lf %.20lf %20lf %20lf\n", T[i], Tjf[i], V[i], X[i], Y[i]);
    }

    return (Vcum+Tcum)/steps;

}

double monte_carlo(int steps, double delta){
    double out[3];
    int acc = 0;
    double Tcum = 0, Tjfcum = 0, Vcum = 0;

    init_matrix();

    gsl_permutation *pa = gsl_permutation_alloc(dim1);
    gsl_permutation *pb = gsl_permutation_alloc(dim2);
    gsl_matrix *Ma = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *Mb = gsl_matrix_alloc(dim2, dim2);
    int sa, sb;

    gsl_matrix_memcpy(Ma, A);
    gsl_matrix_memcpy(Mb, B);

    gsl_linalg_LU_decomp(Ma, pa, &sa);
    gsl_linalg_LU_decomp(Mb, pb, &sb);

    double detai, detbi, detf;

    detai = gsl_linalg_LU_det(Ma, sa);
    detbi = gsl_linalg_LU_det(Mb, sb);


    gsl_matrix *tempa = gsl_matrix_alloc(dim1, dim1);
    gsl_matrix *tempb = gsl_matrix_alloc(dim2, dim2);
    gsl_permutation *temppa = gsl_permutation_alloc(dim1);
    gsl_permutation *temppb = gsl_permutation_alloc(dim2);

    int temps;
    double xtmp, ytmp, prob;

    for(int k=0;k<steps;k++){
        for(int i=0;i<dim1;i++){
            xtmp = xa[i] + delta*(rand()/(double)RAND_MAX -0.5);
            ytmp = ya[i] + delta*(rand()/(double)RAND_MAX -0.5);

            gsl_matrix_memcpy(tempa, A);
            for(int j=0;j<dim1;j++){
                gsl_matrix_set(tempa, i, j, state(xtmp, ytmp, 1, orba[j]));
            }

            gsl_matrix_memcpy(Ma, tempa);
            gsl_linalg_LU_decomp(Ma, temppa, &temps);

            detf = gsl_linalg_LU_det(Ma, temps);
            prob = detf*detf/detai/detai;
            //cout<<"detf= "<<detf<<" deta= "<<detai<<" detb= "<<detbi<<" prob= "<<prob<<endl;
            if(prob > rand()/(double)RAND_MAX){
                xa[i] = xtmp;
                ya[i] = ytmp;
                gsl_matrix_memcpy(A, tempa);
                gsl_permutation_memcpy(pa, temppa);
                detai = detf;

                acc+=1;
            }
        }
        gsl_matrix_memcpy(Ma, A);
        gsl_linalg_LU_decomp(Ma, pa, &sa);

        for(int i=0;i<dim2;i++){
            xtmp = xb[i] + delta*(rand()/(double)RAND_MAX -0.5);
            ytmp = yb[i] + delta*(rand()/(double)RAND_MAX -0.5);

            gsl_matrix_memcpy(tempb, B);
            for(int j=0;j<dim2;j++){
                gsl_matrix_set(tempb, i, j, state(xtmp, ytmp, 1, orbb[j]));
            }
            gsl_matrix_memcpy(Mb, tempb);
            gsl_linalg_LU_decomp(Mb, temppb, &temps);

            detf = gsl_linalg_LU_det(Mb, temps);
            prob = detf*detf/detbi/detbi;
            //cout<<"detf= "<<detf<<" deta= "<<detai<<" detb= "<<detbi<<" prob= "<<prob<<endl;
            if(prob > rand()/(double)RAND_MAX){
                xb[i] = xtmp;
                yb[i] = ytmp;
                gsl_matrix_memcpy(B, tempb);
                gsl_permutation_memcpy(pb, temppb);
                detbi = detf;
                //cout<<"A "<<gsl_matrix_get(A,0,0)<<" "<<gsl_matrix_get(A,1,0)<<" "<<gsl_matrix_get(A,0,1)<<" "<<gsl_matrix_get(A,1,1)<<endl;
                //cout<<"B "<<gsl_matrix_get(tempb,0,0)<<" "<<gsl_matrix_get(tempb,1,0)<<" "<<gsl_matrix_get(tempb,0,1)<<" "<<gsl_matrix_get(tempb,1,1)<<endl;

                acc+=1;
            }
        }
        gsl_matrix_memcpy(Mb, B);
        gsl_linalg_LU_decomp(Mb, pb, &sb);

        calculateT(out, Ma, Mb, pa, pb);

        Tcum += out[0];
        Tjfcum += out[1];
        Vcum += pot();
    }

    return (Vcum+Tcum)/steps;

}

int main(){
    dim1 = 3;
    dim2 = 3;

    orba = new int[dim1] {0, 1, 2};
    orbb = new int[dim2] {0, 1, 2};

    A   = gsl_matrix_alloc(dim1, dim1);
    dxA = gsl_matrix_alloc(dim1, dim1);
    dyA = gsl_matrix_alloc(dim1, dim1);
    ddA = gsl_matrix_alloc(dim1, dim1);
    B   = gsl_matrix_alloc(dim2, dim2);
    dxB = gsl_matrix_alloc(dim2, dim2);
    dyB = gsl_matrix_alloc(dim2, dim2);
    ddB = gsl_matrix_alloc(dim2, dim2);


    xa = new double[dim1] {1, .7,.3};
    ya = new double[dim1] {1, 0, 2};
    xb = new double[dim2] {0, 2, 1};
    yb = new double[dim2] {1, .3, .5};

    srand(29876527);

    monte_carlo_out(10000, 1.);

}
