#ifndef __System__
#define __System__

// Standard C++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Standard C library for some system functions
#include <stdlib.h>

// Custom random number generator header
#include "../rangen/random.h"

// Using directives to avoid repetitive namespace qualification
using namespace std;

// Declaration of the "Metropolis" class which represents a Metropolis algorithm simulation
class Metropolis {
private:
    // Private members representing the state of the simulation
    double _xstart;         // Initial position
    double _x;              // Current position
    int _Ntot, _Nacc;       // Total and accepted moves
    double _acceptance;     // Acceptance ratio
    double _sigma, _m;      // Parameters for the wavefuction
    double _delta_h;        // Step size for the Metropolis algorithm
    int _Niter, _Nblocks;   // Number of iterations and blocks
    double _E, _errE;       // Energy and error in energy
    Random* _rnd;           // Pointer to random number generator instance

public:
    // Constructor with initialization list for member variables
    Metropolis(double xstart, double m, double sigma, double delta_h, int Niter, int Nblocks, Random* rnd)
        : _xstart(xstart), _x(xstart), _sigma(sigma), _m(m), _delta_h(delta_h),
          _Niter(Niter), _Nblocks(Nblocks), _rnd(rnd), _E(0), _errE(0) {}

    // Setters for mean and sigma parameters
    void setmean(double mean) { _m = mean; }
    void setsigma(double sigma) { _sigma = sigma; }

    // Getters for various properties of the simulation state
    double getacceptance() const { return _acceptance; }
    double getmean() const { return _m; }
    double getsigma() const { return _sigma; }
    double getE() const { return _E; }
    double geterrE() const { return _errE; }

    // Methods related to the Metropolis algorithm
    double Ham();                      // Calculate Hamiltonian
    double prob(double x);             // Calculate probability
    void Metr();                       // Perform Metropolis algorithm
    void computeEnergy(ofstream &fout, bool printx); // Compute energy with output to file
    void computeEnergy(ofstream &fout, ofstream &flusso_out, bool printx); // Compute energy with output to multiple files
    void computeEnergy();              // Compute energy without file output
    double err(double media, double mediaquad, int n); // Calculate error
};

// Declaration of the "SA" (Simulated Annealing) class
class SA {
private:
    // Private members representing the state of the simulated annealing process
    double _sigmastart, _mstart; // Initial values for sigma and mean
    double _sigma, _m;           // Current values for sigma and mean
    double _delta_s, _delta_m;   // Step sizes for sigma and mean
    double _betastart, _beta, _beta_supp; // Initial and current values for beta
    double _E, _errE, _deltaE;   // Energy, error in energy, and change in energy
    int _N, _Nacc, _Ncheck;      // Iteration count, accepted moves count, and check
    double _acceptance, _e;      // Acceptance ratio and adjust interval
    Metropolis _M;               // Metropolis instance used within SA
    Random* _rnd;                // Pointer to random number generator instance
    int _type;                   // Type of simulated annealing
    bool _change;                // Flag to indicate if parameters changed
    double _mean_m, _mean_s;     // Mean values for sigma and mean

public:
    // Constructor for Simulated Annealing class
    SA(Metropolis M, double betastart, Random* rnd, int type)
        : _M(M), _betastart(betastart), _rnd(rnd), _type(type),
          _beta(betastart), _m(M.getmean()), _mstart(M.getmean()),
          _sigma(M.getsigma()), _sigmastart(M.getsigma()) {
        _delta_s = 1.0 / _betastart;
        _delta_m = 1.0 / _betastart;
        _M.computeEnergy();
        _E = _M.getE();
        _errE = _M.geterrE();
        _N = 0;
        _Nacc = 0;
        _beta_supp = _beta;
        _change = false;
        _e = 1;
    }

    // Getters for various properties of the simulated annealing state
    double getE() const { return _E; }
    double geterrE() const { return _errE; }
    double getsigma() { return _sigma; }
    double getmean() { return _m; }
    double getmean_sigma() { return _mean_s; }
    double getmean_mean() { return _mean_m; }
    double getbeta() { return _beta; }

    // Methods related to Simulated Annealing process
    void evolve();       // Evolve the system
    void sample_at_fixed_params(); // Sampling with fixed beta
    void Metr_SA();      // Perform Metropolis algorithm for Simulated Annealing
    void change_beta();  // Change the beta parameter
};

// Function declaration to start the simulation
void start(int* type, int* n_block, int* n_iter);

#endif // __System__
