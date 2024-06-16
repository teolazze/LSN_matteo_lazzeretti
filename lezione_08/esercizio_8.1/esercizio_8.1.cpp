#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"
#include "../system.h"

using namespace std;

int main (int argc, char *argv[]){

    int n_block = 100; // Numero di blocchi

    ofstream f1out, f2out, f3out; // Dichiarazione degli stream di output per i file di campionamento
    ofstream flusso_out1, flusso_out2, flusso_out3; // Dichiarazione degli stream di output per i file di media

    // Apertura dei file di output del campionamento
    f1out.open("OUTPUT/sampling1.dat");
    f2out.open("OUTPUT/sampling2.dat");
    f3out.open("OUTPUT/sampling3.dat");

    // Apertura dei file di output della media in funzione del blocco
    flusso_out1.open("OUTPUT/mean1.dat");
    flusso_out2.open("OUTPUT/mean2.dat");
    flusso_out3.open("OUTPUT/mean3.dat");
 
    double mean = 0.8; // Valore medio
    double sigma = 0.4; // Deviazione standard
    double x_start = mean; // Valore iniziale
    double delta_h = 3; // Passo di Metropolis
    Random rnd = generaRand(); // Generazione dell'oggetto Random

    // Creazione del primo oggetto Metropolis con 1000 passi
    Metropolis M1(x_start, mean, sigma, delta_h, 1000, n_block, &rnd);
    M1.computeEnergy(f1out, flusso_out1, true); // Calcolo dell'energia e scrittura nei file di output

    // Creazione del secondo oggetto Metropolis con 1000 passi
    Metropolis M2(x_start, mean, sigma, delta_h, 1000, n_block, &rnd);
    M2.computeEnergy(f2out, flusso_out2, true); // Calcolo dell'energia e scrittura nei file di output

    // Creazione del terzo oggetto Metropolis con 10000 passi
    Metropolis M3(x_start, mean, sigma, delta_h, 10000, n_block, &rnd);
    M3.computeEnergy(f3out, flusso_out3, true); // Calcolo dell'energia e scrittura nei file di output

    return 0; // Fine del programma
}
