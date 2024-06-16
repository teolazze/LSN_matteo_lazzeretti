#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"
#include <vector>
#include <unordered_set>

using namespace std;

//////////////////////////////////////
// METROPOLIS CLASS
//////////////////////////////////////

// Calculate the Hamiltonian (energy) of the system
double Metropolis::Ham() {
    // Kinetic and potential terms of the Hamiltonian
    double kinetic = -0.5 * (pow(_x, 2) + pow(_m, 2) - pow(_sigma, 2) - 2 * _x * _m * tanh((_x * _m) / pow(_sigma, 2))) / pow(_sigma, 4);
    double potential = pow(_x, 4) - 5.0 / 2 * pow(_x, 2);
    return kinetic + potential;
}

// Calculate the probability density function
double Metropolis::prob(double x) {
    double psi = exp(-pow((x - _m), 2) / (2 * pow(_sigma, 2))) + exp(-pow((x + _m), 2) / (2 * pow(_sigma, 2)));
    return psi * psi;
}

// Perform a Metropolis move
void Metropolis::Metr() {
    // Propose a new position
    double x = _x + _rnd->Rannyu(-_delta_h, _delta_h);
    // Acceptance probability
    double p = this->prob(x) / this->prob(_x);
    // Accept or reject the move
    if (p > 1 || _rnd->Rannyu() < p) {
        _Nacc++;
        _x = x;
    }
    _Ntot++;
    _acceptance = static_cast<double>(_Nacc) / _Ntot;
}

// Calculate the energy of the system and write it to file
void Metropolis::computeEnergy(ofstream &fout, bool printx) {
    _acceptance = 0;
    _Nacc = 0;
    _Ntot = 0;
    double sum_quad_block = 0;
    double sum_block = 0;

    // Loop over blocks
    for (int k = 0; k < _Nblocks; k++) {
        double sum_parz = 0;
        // Loop over iterations within a block
        for (int i = 0; i < _Niter; i++) {
            if (printx) fout << _x << endl; // Print current position if required
            this->Metr(); // Perform a Metropolis move
            sum_parz += this->Ham(); // Accumulate energy contribution
        }
        sum_block += sum_parz / _Niter; // Average energy per iteration in the block
        sum_quad_block += pow(sum_parz / _Niter, 2); // Square of the above for error calculation
    }

    _E = sum_block / _Nblocks; // Mean energy over all blocks
    _errE = this->err(sum_block / _Nblocks, sum_quad_block / _Nblocks, _Nblocks); // Statistical error in mean energy
}

// Calculate the energy of the system and write it to two different files
void Metropolis::computeEnergy(ofstream &f1out, ofstream &f2out, bool printx) {
    _acceptance = 0;
    _Nacc = 0;
    _Ntot = 0;
    double sum_quad_block = 0;
    double sum_block = 0;

    // Loop over blocks
    for (int k = 0; k < _Nblocks; k++) {
        double sum_parz = 0;
        // Loop over iterations within a block
        for (int i = 0; i < _Niter; i++) {
            if (printx) f1out << _x << endl; // Print current position if required
            this->Metr(); // Perform a Metropolis move
            sum_parz += this->Ham(); // Accumulate energy contribution
        }
        sum_block += sum_parz / _Niter; // Average energy per iteration in the block
        sum_quad_block += pow(sum_parz / _Niter, 2); // Square of the above for error calculation
        f2out << sum_block / (k + 1) << "   " << this->err(sum_block / (k + 1), sum_quad_block / (k + 1), k + 1) << "   " << _acceptance << endl; // Write intermediate results to f2out
    }

    _E = sum_block / _Nblocks; // Mean energy over all blocks
    _errE = this->err(sum_block / _Nblocks, sum_quad_block / _Nblocks, _Nblocks); // Statistical error in mean energy
}

// Calculate the energy of the system without writing it to file
void Metropolis::computeEnergy() {
    _acceptance = 0;
    _Nacc = 0;
    _Ntot = 0;
    double sum_quad_block = 0;
    double sum_block = 0;

    // Loop over blocks
    for (int k = 0; k < _Nblocks; k++) {
        double sum_parz = 0;
        // Loop over iterations within a block
        for (int i = 0; i < _Niter; i++) {
            this->Metr(); // Perform a Metropolis move
            sum_parz += this->Ham(); // Accumulate energy contribution
        }
        sum_block += sum_parz / _Niter; // Average energy per iteration in the block
        sum_quad_block += pow(sum_parz / _Niter, 2); // Square of the above for error calculation
    }

    _E = sum_block / _Nblocks; // Mean energy over all blocks
    _errE = this->err(sum_block / _Nblocks, sum_quad_block / _Nblocks, _Nblocks); // Statistical error in mean energy
}

// Calculate the statistical error
double Metropolis::err(double mean, double mean_squared, int n) {
    if (n == 1) return 0; // No error if there's only one sample
    return sqrt((mean_squared - mean * mean) / (n - 1)); // Calculate standard error using mean and mean of squares
}

//////////////////////////////////////
// SIMULATED ANNEALING
//////////////////////////////////////

// Perform a Metropolis move in the context of Simulated Annealing
void SA::Metr_SA() {
    // Propose new values for mean and sigma
    double mean = _m + _rnd->Rannyu(-_delta_m, _delta_m);
    double sigma = _sigma + _rnd->Rannyu(-_delta_s, _delta_s);

    _M.setmean(mean); // Set new mean
    _M.setsigma(sigma); // Set new sigma
    _M.computeEnergy(); // Calculate energy with new parameters

    double newEnergy = _M.getE(); // Get new energy
    double p = exp(-_beta * (newEnergy - _E)); // Metropolis acceptance probability

    // Accept or reject the move
    if (p > 1 || _rnd->Rannyu() < p) {
        _change = true; // Accepted move flag
        _deltaE = newEnergy - _E; // Energy change
        _E = newEnergy; // Update energy
        _errE = _M.geterrE(); // Update energy error
        _m = mean; // Update mean
        _sigma = sigma; // Update sigma
        _Nacc++; // Increase acceptance count
    } else {
        _deltaE = 0; // No energy change if rejected
        _change = false; // Rejected move flag
    }

    _M.setmean(_m); // Reset mean in Metropolis object
    _M.setsigma(_sigma); // Reset sigma in Metropolis object
    _N++; // Increase step count
    _acceptance = static_cast<double>(_Nacc) / _N; // Acceptance ratio
}

// Evolve the system using Simulated Annealing
void SA::evolve() {
    double lim = 10000; // Arbitrary limit for beta or cooling parameter
    string filename; // Filename for output data

    // Determine output file name based on cooling type
    ofstream fout;
    if (_type == 0) filename = "OUTPUT_LIN/data.dat";
    if (_type == 1) filename = "OUTPUT_QDR/data.dat";
    if (_type == 2) filename = "OUTPUT_GEO/data.dat";

    fout.open(filename); // Open output file
    fout << "           beta           mean          sigma           E           errE" << endl; // Output file header

    // Perform Simulated Annealing
    vector<double> energy_history; // Vector to store recent mean energy values
    const int group_size = 10; // Number of groups for stability check
    int stable_count = 0; // Counter for stable energy values
    double last_mean_energy = 0.0; // Last computed mean energy
    double last_beta = 0.0; // Last used beta value

    do {
        this->change_beta(); // Adjust beta according to chosen cooling scheme
        this->Metr_SA(); // Perform Metropolis move in Simulated Annealing context

        if (_change) {
            // Write current state to file
            fout << setw(15) << _beta << setw(15) << _m << setw(15) << _sigma << setw(15) << _E << setw(15) << _errE << endl; // Write current state to file

            energy_history.push_back(_E); // Store current mean energy

            // Check for stability if enough data points are available
            if (energy_history.size() >= group_size) {
                double sum_energy = 0.0;
                // Calculate mean energy over last `group_size` entries
                for (int i = energy_history.size() - group_size; i < energy_history.size(); i++) {
                    sum_energy += energy_history[i];
                }
                double mean_energy = sum_energy / group_size;

                // Check stability condition
                if (abs(mean_energy - last_mean_energy) < 3 * _errE) {
                    stable_count++; // Increment stable count
                } else {
                    stable_count = 0; // Reset stable count if condition is not met
                }

                // Update last computed mean energy
                last_mean_energy = mean_energy;

                energy_history.clear(); // Clear energy history for next round
                last_beta = _beta; // Update last beta value used
            }
        }

        // Continue loop until energy values stabilize for 10 iterations or beta exceeds limit
    } while (stable_count < 10 && _beta < 10000);

    fout.close(); // Close output file
}

// Change the beta parameter based on the chosen cooling scheme
void SA::change_beta() {
    if (_type == 0) {
        _beta += 1; // Linear cooling scheme
    } else if (_type == 1) {
        _beta_supp += 0.01; // Quadratic cooling scheme
        _beta = _beta_supp * _beta_supp;
    } else if (_type == 2) {
        _beta *= 1.005; // Geometric cooling scheme
    }

    if (!_change) _e += 1; // Adjust counter based on acceptance
    if (_change) _e = 1;

    _delta_s = 0.01 / _e; // Update step sizes
    _delta_m = 0.01 / _e;
}

// Perform sampling at fixed parameters
void SA::sample_at_fixed_params() {
    ofstream f1out, f2out; // Output file streams
    string filename1, filename2; // Filenames for output

    int N = 5000; // Number of iterations
    int block_size = 100; // Block size for averaging
    int num_blocks = 50; // Number of blocks

    _Nacc = 0; // Reset acceptance count
    _N = 0; // Reset step count

    vector<double> block_means_m(num_blocks, 0.0); // Vector for mean values of mean
    vector<double> block_means_s(num_blocks, 0.0); // Vector for mean values of sigma

    _delta_s = 0.01; // Initial step size for sigma
    _delta_m = 0.01; // Initial step size for mean

    // Determine output filenames based on cooling type
    if (_type == 0) filename1 = "OUTPUT_LIN/params.dat";
    if (_type == 1) filename1 = "OUTPUT_QDR/params.dat";
    if (_type == 2) filename1 = "OUTPUT_GEO/params.dat";

    if (_type == 0) filename2 = "OUTPUT_LIN/mu&sigma.dat";
    if (_type == 1) filename2 = "OUTPUT_QDR/mu&sigma.dat";
    if (_type == 2) filename2 = "OUTPUT_GEO/mu&sigma.dat";

    f1out.open(filename1); // Open first output file
    f1out << " beta mean sigma E errE acc" << endl; // Output file header
    f2out.open(filename2); // Open second output file
    f2out << " beta mu errmu sigma errsigma" << endl; // Output file header

    // Loop over blocks
    for (int i = 0; i < num_blocks; i++) {
        double sum_m = 0.0, sum_s = 0.0;

        // Loop over iterations within a block
        for (int j = 0; j < block_size; j++) {
            this->Metr_SA(); // Perform Metropolis move in Simulated Annealing context

            sum_m += _m; // Accumulate mean values
            sum_s += _sigma; // Accumulate sigma values

            // Write current state to file
            f1out << setw(15) << _beta << setw(15) << _m << setw(15) << _sigma << setw(15) << _E << setw(15) << _errE << setw(15) << _acceptance << endl;
        }

        block_means_m[i] = sum_m / block_size; // Calculate mean of means for block
        block_means_s[i] = sum_s / block_size; // Calculate mean of sigmas for block
    }

    double mean_m = 0.0, mean_s = 0.0;
    // Calculate overall mean of means and sigmas
    for (int i = 0; i < num_blocks; i++) {
        mean_m += block_means_m[i];
        mean_s += block_means_s[i];
    }

    mean_m /= num_blocks; // Mean of mean
    mean_s /= num_blocks; // Mean of sigma

    double var_m = 0.0, var_s = 0.0;
    // Calculate variance of means and sigmas
    for (int i = 0; i < num_blocks; i++) {
        var_m += (block_means_m[i] - mean_m) * (block_means_m[i] - mean_m);
        var_s += (block_means_s[i] - mean_s) * (block_means_s[i] - mean_s);
    }

    var_m /= (num_blocks - 1); // Variance of mean
    var_s /= (num_blocks - 1); // Variance of sigma

    double err_m = sqrt(var_m / num_blocks); // Standard error of mean
    double err_s = sqrt(var_s / num_blocks); // Standard error of sigma

    // Write mean and error values to second output file
    f2out << setw(15) << _beta << setw(15) << mean_m << setw(15) << var_m << setw(15) << mean_s << setw(15) << var_s << setw(15) << endl;

    _mean_m = mean_m; // Store mean of mean
    _mean_s = mean_s; // Store mean of sigma

    f1out.close(); // Close first output file
    f2out.close(); // Close second output file
}

// Function to read input parameters and initialize simulation
void start(int* type, int* n_block, int* n_iter) {
    ifstream input("input.dat"); // Open input file
    string property; // Property read from input file
    *type = -1; // Initialize simulation type
    *n_block = -1; // Initialize number of blocks
    *n_iter = -1; // Initialize number of iterations

    // Read input properties from file
    while (!input.eof()) {
        input >> property; // Read property
        if (property == "COOLING") {
            string sim_type; // Simulation type
            input >> sim_type; // Read simulation type
            bool check = false; // Check flag for valid simulation type
            if (sim_type == "lin") {
                *type = 0; // Linear cooling
                check = true; // Set check flag
            } else if (sim_type == "qdr") {
                *type = 1; // Quadratic cooling
                check = true; // Set check flag
            } else if (sim_type == "geo") {
                *type = 2; // Geometric cooling
                check = true; // Set check flag
            }
            if (!check) {
                throw invalid_argument("Invalid argument SIMULATION TYPE."); // Throw exception for invalid simulation type
            }
        } else if (property == "N_BLOCK") {
            input >> *n_block; // Read number of blocks
        } else if (property == "N_ITERATION") {
            input >> *n_iter; // Read number of iterations
        } else {
            cerr << "PROBLEM: unknown input" << endl; // Print error message for unknown input
        }
    }

    // Validate input parameters
    if (*type == -1) {
        throw invalid_argument("Missing SIMULATION TYPE input."); // Throw exception for missing simulation type
    }
    if (*n_block < 0) {
        throw invalid_argument("Wrong N_BLOCK input."); // Throw exception for invalid number of blocks
    }
    if (*n_iter < 0) {
        throw invalid_argument("Wrong N_ITERATION input."); // Throw exception for invalid number of iterations
    }

    input.close(); // Close input file
    return; // Return from function
}


///////////////////////////////////////////////////////////////////
//Prima bozza
///////////////////////////////////////////////////////////////////
// // Evolve the system using Simulated Annealing
// void SA::evolve() {
// double lim = 10000;
// string filename;

// // Determine output filename based on cooling type
// ofstream fout;
// if (_type == 0) filename = "OUTPUT_LIN/data.dat";
// if (_type == 1) filename = "OUTPUT_QDR/data.dat";
// if (_type == 2) filename = "OUTPUT_GEO/data.dat";

// fout.open(filename);
// fout << " beta mean sigma E errE" << endl;

// // Perform simulated annealing
// do{
// this->change_beta(); // Change the temperature (beta) according to the chosen schedule
// this->Metr_SA(); // Perform a Metropolis step in the Simulated Annealing context

// // Output the current state to the file
// fout << setw(15) << _beta << setw(15) << _m << setw(15) << _sigma << setw(15) << _E << setw(15) << _errE << endl;

// // Check if the error in energy is greater than the absolute change in energy
// if (_errE > abs(_deltaE)) {
// if (_change) _Ncheck++; // Increment the check counter if a change was accepted
// } else {
// _Ncheck = 0; // Reset the check counter if the error condition is not met
// }

// // Continue the loop under the following conditions:
// // 1. The check counter is less than 10, or no moves have been accepted, or the change in energy is greater than 1.
// // 2. The beta value is less than the limit (10000).
// // 3. The error count (_e) is less than 15.
// }while((_Ncheck < 10 || _Nacc == 0 || abs(_deltaE) > 1) && _beta < lim && _e < 15);
// }