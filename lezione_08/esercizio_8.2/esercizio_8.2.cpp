#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"  // Include for random number generator
#include "../system.h"  // Include for Metropolis and SA classes

using namespace std;


int main (int argc, char *argv[]){
    // Declare variables to hold simulation parameters
    int type, n_block, n_iter;

    // Read and initialize the simulation parameters from input file
    start(&type, &n_block, &n_iter);

    ofstream flusso_out;  // Output file stream
    string filename1;  // Filename for output

    // Determine output filename based on simulation type
    if(type == 0) filename1 = "OUTPUT_LIN/sampling.dat";
    if(type == 1) filename1 = "OUTPUT_QDR/sampling.dat";
    if(type == 2) filename1 = "OUTPUT_GEO/sampling.dat";

    flusso_out.open(filename1);  // Open output file for writing

    ofstream fout;  // Output file stream
    string filename2;  // Filename for output

    // Determine output filename based on simulation type
    if(type == 0) filename2 = "OUTPUT_LIN/block_mean.dat";
    if(type == 1) filename2 = "OUTPUT_QDR/block_mean.dat";
    if(type == 2) filename2 = "OUTPUT_GEO/block_mean.dat";

    fout.open(filename2);  // Open output file for writing

    // Initialize parameters for Metropolis algorithm
    double mean = 0.6;
    double sigma = 0.6;
    double x_start = mean;
    double delta_h = 3;

    // Generate random number generator
    Random rnd = generaRand();

    // Initialize Metropolis object
    Metropolis M(x_start, mean, sigma, delta_h, n_iter, n_block, &rnd);

    // Initialize Simulated Annealing (SA) system with the Metropolis object and starting beta of 10
    SA system(M, 10, &rnd, type);

    // Perform the simulated annealing process
    system.evolve();

    system.sample_at_fixed_params();

    // Create a new Metropolis object with optimized mean and sigma from SA and new parameters
    Metropolis M1(x_start, system.getmean_mean(), system.getmean_sigma(), delta_h, 10000, 100, &rnd);

    // Compute energy using the new Metropolis object and write the results to the output file
    M1.computeEnergy(flusso_out, fout, true);

    return 0;
}
