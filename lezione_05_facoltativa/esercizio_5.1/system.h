#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include <vector>
#include "../../rangen/random.h"

using namespace std;
using namespace arma;

// Dichiarazione della classe "Metropolis" che rappresenta una città
class Metropolis {
private:
    vector<double> _p_start;
    vector<double> _p;
    double _a0;
    int _Ntot, _Nacc;
    int _type, _type_sampling;
    double _acceptance;
    double _delta;
    int _Niter, _Nblocks, _Nskip;
    Random* _rnd;          // Generatore di numeri casuali
    

public:
    Metropolis(vector<double> p_start, double delta, int Niter, int Nblocks, int Nskip, int type, int type_sampling, Random* rnd) :_Nskip(Nskip), _type_sampling(type_sampling), _type(type), _Nblocks(Nblocks), _p_start(p_start), _p(p_start), _delta(delta), _Niter(Niter), _rnd(rnd) {_Ntot = 0; _Nacc = 0; _acceptance = 0; _a0 = 1;}
    double getacceptance(){return _acceptance; }
    double ground_state(vector<double> x);
    double excited_state(vector<double> x);
    void set_sampling(int i, double delta){_type_sampling = i; _delta = delta;};
    void set_state(int i){_type = i;}
    void Metr();
    void expected_r(ofstream &f1out, ofstream &f2out);
    double err(double media, double mediaquad, int n);

};

void start(int* type, int* n_block, int* n_iter);

#endif // __System__