#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"  // Assuming this is where Random and generaRand() are defined
#include "system.h"

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd = generaRand();  // Initialize random number generator

    vector<city> cities;  // Vector to store cities

    int N_city;  // Number of cities
    int N_iter;  // Number of iterations
    int N_path;  // Number of paths
    int type;    // Simulation type

    start(&type, &N_city, &N_iter, &N_path);  // Initialize simulation parameters from input

    ofstream f1out, f2out, fcout, ffout, fmout;  // Output file streams

    // Depending on the simulation type, generate cities and open output files
    if (type == 0) {  // Circular arrangement of cities
        for (int k = 0; k < N_city; k++) {
            double theta = rnd.Rannyu(0, 2 * M_PI);  // Random angle
            city c(cos(theta), sin(theta));  // Create city at (cos(theta), sin(theta))
            cities.push_back(c);  // Add city to vector
        }
        f1out.open("../OUTPUT_CIRCLE/starting_pop.dat");  // Open output file for starting population
        f2out.open("../OUTPUT_CIRCLE/final_pop.dat");     // Open output file for final population
        fcout.open("../OUTPUT_CIRCLE/city.dat");          // Open output file for city coordinates
        fmout.open("../OUTPUT_CIRCLE/mean.dat");          // Open output file for mean distances
        ffout.open("../OUTPUT_CIRCLE/firstpath.dat");     // Open output file for first path
    }
    if (type == 1) {  // Square arrangement of cities
        for (int k = 0; k < N_city; k++) {
            double x = rnd.Rannyu(-1, 1);  // Random x coordinate in [-1, 1)
            double y = rnd.Rannyu(-1, 1);  // Random y coordinate in [-1, 1)
            city c(x, y);  // Create city at (x, y)
            cities.push_back(c);  // Add city to vector
        }
        f1out.open("../OUTPUT_SQUARE/starting_pop.dat");  // Open output file for starting population
        f2out.open("../OUTPUT_SQUARE/final_pop.dat");     // Open output file for final population
        fcout.open("../OUTPUT_SQUARE/city.dat");          // Open output file for city coordinates
        fmout.open("../OUTPUT_SQUARE/mean.dat");          // Open output file for mean distances
        ffout.open("../OUTPUT_SQUARE/firstpath.dat");     // Open output file for first path
    }
    if (type == 2) {  // Read cities from a configuration file
        ifstream infile("config.dat");  // Open configuration file
        double x, y;
        while (infile >> x >> y) {
            cities.push_back(city(x, y));  // Read city coordinates and add to vector
        }
        infile.close();  // Close configuration file

        f1out.open("../OUTPUT_SQUARE/starting_pop.dat");  // Open output file for starting population
        f2out.open("../OUTPUT_SQUARE/final_pop.dat");     // Open output file for final population
        fcout.open("../OUTPUT_SQUARE/city.dat");          // Open output file for city coordinates
        fmout.open("../OUTPUT_SQUARE/mean.dat");          // Open output file for mean distances
        ffout.open("../OUTPUT_SQUARE/firstpath.dat");     // Open output file for first path
    }

    vector<int> v;
    for (int i = 0; i < N_city; i++) v.push_back(i);  // Initialize a vector with city indices

    path my_path(v, &rnd);  // Create a path object with the vector of city indices
    my_path.check();  // Validate the path

    System my_syst(my_path, &rnd, N_path);  // Initialize the system with the path and random generator

    my_syst.stamp_city(cities, fcout);  // Write city coordinates to file

    my_syst.initialize();  // Initialize the system with random mutations

    my_syst.check();  // Validate the initial paths

    my_syst.stamp_path(f1out);  // Write starting paths to file

    my_syst.sortL1(cities);  // Sort paths based on L1 distance

    // Perform evolutionary iterations
    for (int i = 0; i < N_iter; i++) {
        my_syst.crossover2(0.7, 1.6);  // Perform crossover with specified probability and potential
        my_syst.mutation(0.1, 0.1, 0.1, 0.1);  // Perform mutation with specified probabilities
        my_syst.sortL2(cities);  // Sort paths based on L2 distance
        my_syst.check();  // Validate paths after operations
        my_syst.restart();  // Restart paths
        my_syst.printfirstpath(ffout, cities);  // Print the first path to file
        my_syst.printmean(fmout, cities, N_path/2, 1, i);  // Print mean distances to file

        printProgressBar(i + 1, N_iter);  // Print progress bar
    }

    cout << endl;  // Print newline after progress bar completion

    my_syst.stamp_path(f2out);  // Write final paths to file

    f1out.close();  // Close output files
    f2out.close();
    fcout.close();
    fmout.close();
    ffout.close();
}
