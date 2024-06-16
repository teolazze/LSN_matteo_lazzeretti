#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../rangen/random.h"

using namespace std;

double err(double media, double mediaquad, int n){
    if(n == 0) return 0;
    else{ 
        return sqrt((mediaquad - media * media) / n);
    }
}

int main (int argc, char *argv[]){

    ofstream flusso_out;
    flusso_out.open("risultati.dat");

    
    int M = 1E9;
    int N = 1000;
    int L = M/N;

    double d = 1;
    double l = 0.8;


    Random rnd = generaRand();

    double sumquadAi = 0;
    double sumAi = 0;
    double mediaAi = 100;

    for(int i = 0; i < N; i++){
        int c = 0;

        for (int j = 0; j < L; j++){
            float theta_rand = rnd.Rannyu(0,2*mediaAi);
            float z_rand = rnd.Rannyu();
            double z1 = z_rand;
            double z2 = z_rand+l*cos(theta_rand);
            if((z1 < 0 || z1 > 1) || (z2 < 0 || z2 > 1)) c++;
        }

        double Ai = 2*l*L/(c*d);
        sumAi += Ai;
        mediaAi = sumAi/(i+1);        
        sumquadAi += Ai*Ai;
        double mediaquadAi = sumquadAi/(i+1);
        double sigma1 = err(mediaAi, mediaquadAi, i);

        flusso_out << i+1 << " " << mediaAi << " " << sigma1 << endl;

    }

    flusso_out.close();
    return 0;
}
