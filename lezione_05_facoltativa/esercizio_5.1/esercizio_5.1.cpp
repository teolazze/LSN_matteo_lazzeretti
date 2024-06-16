#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"
#include "system.h"

using namespace std;

double lin_beta(double N){
    return N;
}


int main (int argc, char *argv[]){

    ofstream f1out, f2out, f3out, f4out;
    ofstream flusso_out1, flusso_out2, flusso_out3, flusso_out4;

    f1out.open("OUTPUT/GS/uniform/expr.dat");
    f2out.open("OUTPUT/GS/gauss/expr.dat");
    f3out.open("OUTPUT/ES/uniform/expr.dat");
    f4out.open("OUTPUT/ES/gauss/expr.dat");

    flusso_out1.open("OUTPUT/GS/uniform/sampling.dat");
    flusso_out2.open("OUTPUT/GS/gauss/sampling.dat");
    flusso_out3.open("OUTPUT/ES/uniform/sampling.dat");
    flusso_out4.open("OUTPUT/ES/gauss/sampling.dat");

    vector<double> p_start = {15,0,0};
    double delta_gs = 1.23, delta_es = 3;
    double sigma_gs = 0.8, sigma_es = 2;
    Random rnd = generaRand();


    Metropolis M(p_start, delta_gs, 1000, 100, 1, 0, 0, &rnd);
    M.expected_r(f1out, flusso_out1);

    M.set_sampling(1, sigma_gs);
    M.expected_r(f2out, flusso_out2);

    M.set_state(1);
    M.set_sampling(1, sigma_es);
    M.expected_r(f4out, flusso_out4);  

    M.set_sampling(0,delta_es);
    M.expected_r(f3out, flusso_out3);  


    return 0;
}