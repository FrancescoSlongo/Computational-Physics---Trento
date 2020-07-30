#include<iostream>
#include<cmath>
#include<iomanip>
#include<vector>
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

int main(){

}
