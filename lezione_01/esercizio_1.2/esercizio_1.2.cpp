#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../rangen/random.h"  // Include la libreria per la generazione di numeri casuali

using namespace std;

int main (int argc, char *argv[]) {

    // Array di ofstream per gestire i tre flussi di output
    ofstream flusso_out[3];
    string name[3] = {"std", "exp", "lor"};  // Nomi dei file di output

    // Crea e apre i file di output
    for (int i = 0; i < 3; i++) {
        string namefile = "risultati_" + name[i] + ".dat";
        flusso_out[i].open(namefile);
    }

    // Numero di campioni e array di valori N
    int M = 1E5;  // Numero di iterazioni
    int N[4] = {1, 2, 10, 100};  // Array dei diversi valori di N

    // Inizializza il generatore di numeri casuali
    Random rnd = generaRand();
    
    // Ciclo principale per generare i campioni
    for (int k = 0; k < M; k++) {
        // Per ogni valore di N
        for (int i = 0; i < 4; i++) {
            double sum_std = 0, sum_exp = 0, sum_lor = 0;  // Somme per le diverse distribuzioni

            // Genera N[i] campioni e somma i valori per ogni distribuzione
            for (int j = 0; j < N[i]; j++) {
                sum_std += rnd.Rannyu();  // Somma campioni uniformi
                sum_exp += rnd.Exp(1);  // Somma campioni esponenziali
                sum_lor += rnd.Lorentz(1, 0);  // Somma campioni lorentziani
            }

            // Scrive le medie nei rispettivi file di output
            flusso_out[0] << sum_std / N[i] << " ";
            flusso_out[1] << sum_exp / N[i] << " ";
            flusso_out[2] << sum_lor / N[i] << " ";
        }

        // Inserisce una nuova linea alla fine di ogni riga nei file di output
        for (int j = 0; j < 3; j++) {
            flusso_out[j] << endl;
        }
    }

    // Chiude tutti i file di output
    for (int j = 0; j < 3; j++) {
        flusso_out[j].close();
    }

    return 0;
}
