#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../rangen/random.h"

using namespace std;


double calculateError(double mean, double meanOfSquares, int n);  // Funzione per calcolare l'errore

int main(int argc, char *argv[]) {

    // Apri il file di output
    ofstream output_file;
    output_file.open("results.dat");

    // Parametri della simulazione
    int totalIter = 1E9;  // Numero di iterazioni totali
    int numBlocks = 10000;  // Numero blocchi totali
    int IterPerBlock = totalIter / numBlocks;
    int numBins = 100;

    // Inizializza il generatore di numeri casuali
    Random rnd = generaRand();

    // Variabili per accumulare statistiche
    double sumOfSquaresMeans = 0, sumOfSquaresVariances = 0;
    double sumOfMeans = 0, sumOfVariances = 0;

    int* binCounts;  // Vettore per contenere il numero di dati in ogni bin
    binCounts = new int[numBins];
    for (int j = 0; j < numBins; j++) {
        binCounts[j] = 0;
    }

    // Ciclo sul numero di blocchi
    for (int i = 0; i < numBlocks; i++) {
        double blockSum = 0, blockVarianceSum = 0;

        // Ciclo sui campioni in ogni blocco
        for (int j = 0; j < IterPerBlock; j++) {
            float randomNum = rnd.Rannyu();
            blockSum += randomNum;
            blockVarianceSum += (randomNum - 0.5) * (randomNum - 0.5);

            binCounts[int(randomNum * numBins)]++;
        }

        // Calcola le medie in funzione del blocco
        double blockMean = blockSum / IterPerBlock;
        sumOfMeans += blockMean;
        double cumulativeMean = sumOfMeans / (i + 1); // Valore medio al blocco i-esimo
        
        sumOfSquaresMeans += blockMean * blockMean;
        double cumulativeMeanOfSquares = sumOfSquaresMeans / (i + 1);
        double meanError = calculateError(cumulativeMean, cumulativeMeanOfSquares, i + 1);  // Errore sul valore medio al blocco i-esimo

        // Calcola le varianze in funzione del blocco
        double blockVariance = blockVarianceSum / IterPerBlock;
        sumOfVariances += blockVariance;
        double cumulativeVariance = sumOfVariances / (i + 1);  // Varianza al blocco i-esimo
        
        sumOfSquaresVariances += blockVariance * blockVariance;
        double cumulativeVarianceOfSquares = sumOfSquaresVariances / (i + 1);
        double varianceError = calculateError(cumulativeVariance, cumulativeVarianceOfSquares, i + 1);  // Errore sulla varianza al blocco i-esimo

        // Calcola il chi-quadro per ogni blocco
        double chiSquared = 0;
        for (int j = 0; j < numBins; j++) {
            chiSquared += pow((binCounts[j] - (double)IterPerBlock / numBins), 2) / ((double)IterPerBlock / numBins);
            binCounts[j] = 0;
        }

        // Scrivi i risultati nel file
        output_file << (i + 1) << " " << cumulativeMean << " " << meanError << " " 
                    << cumulativeVariance << " " << varianceError << " " << chiSquared << endl;
    }

    delete[] binCounts;

    // Chiudi il file di output
    output_file.close();
    return 0;
}




// Funzione per calcolare l'errore
double calculateError(double mean, double meanOfSquares, int n) {
    if (n == 1) return 0;
    else {
        return sqrt((meanOfSquares - mean * mean) / (n-1));
    }
}