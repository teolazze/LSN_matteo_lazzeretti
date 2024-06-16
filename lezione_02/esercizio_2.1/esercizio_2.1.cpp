#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"  // Include la libreria per la generazione di numeri casuali

using namespace std;

// Funzione per il campionamento tramite il metodo Accept-Reject
float AccRej_3(Random* rnd);

// Funzione per calcolare l'errore statistico
double err(double media, double mediaquad, int n);

int main (int argc, char *argv[]) {

    /////////////////////////////////////////////////////////////////////////////////////
    // Questo codice calcola lo stesso integrale campionando con tre diverse funzioni:
    // 1. distribuzione uniforme
    // 2. p(x) = 2(1-x)
    // 3. p(x) = 3/2(1-x^2)
    /////////////////////////////////////////////////////////////////////////////////////

    ofstream flusso_out;  // Stream di output per scrivere i risultati su file
    flusso_out.open("risultati.dat");  // Apre il file "risultati.dat" per scrivere i dati

    int M = 5E4;  // Numero totale di campioni
    int N = 500;  // Numero di blocchi
    int L = M/N;  // Numero di campioni per blocco

    Random rnd = generaRand();  // Inizializza il generatore di numeri casuali

    // Variabili per sommare i valori degli integrali e i loro quadrati
    float somma_I1 = 0, somma_I2 = 0, somma_I3 = 0;
    float somma_I1_quad = 0, somma_I2_quad = 0, somma_I3_quad = 0;

    for(int i = 0; i < N; i++) {
        float sum_parz1 = 0, sum_parz2 = 0, sum_parz3 = 0;
        
        // Esegue il calcolo dell'integrale su ciascun blocco
        for(int j = 0; j < L; j++) {
            float rand1 = rnd.Rannyu();  // Numero casuale uniforme tra 0 e 1
            float rand2 = 1-sqrt(1-rnd.Rannyu());  // Numero casuale con distribuzione p(x) = 2(1-x)
            float rand3 = AccRej_3(&rnd);  // Numero casuale con distribuzione p(x) = 3/2(1-x^2)

            // Calcola l'integrale per ciascun campione usando le diverse distribuzioni
            sum_parz1 += M_PI/2 * cos(M_PI/2 * rand1);
            sum_parz2 += M_PI/2 * cos(M_PI/2 * rand2) / (2 * (1 - rand2));
            sum_parz3 += M_PI/2 * cos(M_PI/2 * rand3) / (3./2 * (1 - rand3 * rand3));
        }

        // Media degli integrali per il blocco corrente

        //Campionamento 1
        float I1 = sum_parz1 / L;
        somma_I1 += I1;
        somma_I1_quad += I1 * I1;
        float sigma_I1 = err(somma_I1 / (i + 1), somma_I1_quad / (i + 1), i + 1);

        //Campionamento 2
        float I2 = sum_parz2 / L;
        somma_I2 += I2;
        somma_I2_quad += I2 * I2;
        float sigma_I2 = err(somma_I2 / (i + 1), somma_I2_quad / (i + 1), i + 1);

        //Campionamento 3
        float I3 = sum_parz3 / L;
        somma_I3 += I3;
        somma_I3_quad += I3 * I3;
        float sigma_I3 = err(somma_I3 / (i + 1), somma_I3_quad / (i + 1), i + 1);

        // Scrive i risultati su file
        flusso_out << i + 1 << " " << somma_I1 / (i + 1) << " " << sigma_I1 << " "
                   << somma_I2 / (i + 1) << " " << sigma_I2 << " "
                   << somma_I3 / (i + 1) << " " << sigma_I3 << endl;
    }

    flusso_out.close();  // Chiude il file di output
    return 0;
}

// Funzione per il campionamento tramite il metodo di Accettazione-Rifiuto
float AccRej_3(Random* rnd) {
    double x = rnd->Rannyu();  // Numero casuale uniforme tra 0 e 1
    double y = rnd->Rannyu(0, 3./2);  // Numero casuale uniforme tra 0 e 3/2
    if(y < 3./2 * (1 - x * x)) {
        return x;  // Accetta il campione
    } else {
        return AccRej_3(rnd);  // Rifiuta il campione e riprova
    }
    return 0;
}

// Funzione per calcolare l'errore statistico
double err(double media, double mediaquad, int n) {
    if(n == 1) return 0;
    else {
        return sqrt((mediaquad - media * media) / (n - 1));
    }
}
