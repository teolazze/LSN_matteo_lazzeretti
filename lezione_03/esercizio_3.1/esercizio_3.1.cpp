#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include "../../rangen/random.h"  // Includi la libreria per il generatore di numeri casuali

using namespace std;

// Dichiarazione delle funzioni
float max(float x1, float x2);  // Funzione per calcolare il massimo tra due numeri float
double err(double media, double mediaquad, int n);  // Funzione per calcolare l'errore standard della media

int main(int argc, char *argv[]) {

    ofstream flusso_out;
    flusso_out.open("risultati.dat");  // Apre il file "risultati.dat" per l'output

    // Definizione delle costanti
    int M = 1E6;  // Numero totale di simulazioni
    int N = 500;  // Numero di blocchi
    int L = M / N;  // Numero di simulazioni per blocco
    int H = 100;  // Numero di intervalli di tempo discreti

    float r = 0.1;  // Tasso di interesse
    float sigma = 0.25;  // Volatilità
    float T = 1;  // Tempo di scadenza
    float K = 100;  // Prezzo di esercizio

    Random rnd = generaRand();  // Inizializzazione del generatore di numeri casuali

    // Variabili per accumulare i risultati
    float sum_dir_C = 0, sum_disc_C = 0;  // Accumulatore per opzioni di tipo Call, per stima diretta e discreta
    float sum_dir_quad_C = 0, sum_disc_quad_C = 0;  // Quadrato dell'accumulatore per opzioni di tipo Call, per stima diretta e discreta
    float sum_dir_P = 0, sum_disc_P = 0;  // Accumulatore per opzioni di tipo Put, per stima diretta e discreta
    float sum_dir_quad_P = 0, sum_disc_quad_P = 0;  // Quadrato dell'accumulatore per opzioni di tipo Put, per stima diretta e discreta

    // Loop sui blocchi
    for (int i = 0; i < N; i++) {
        float sum_parz_dir_C = 0, sum_parz_disc_C = 0;  // Accumulatori parziali per opzioni di tipo Call, per stima diretta e discreta
        float sum_parz_dir_P = 0, sum_parz_disc_P = 0;  // Accumulatori parziali per opzioni di tipo Put, per stima diretta e discreta

        // Loop sulle simulazioni all'interno di ogni blocco
        for (int j = 0; j < L; j++) {
            // Simulazione per opzioni di tipo Call e Put con metodo diretto
            float S_i_dir = 100 * exp((r - sigma * sigma / 2) * T + sigma * rnd.Gauss(0, 1) * sqrt(T));
            float C_i_dir = exp(-r * T) * max(0., S_i_dir - K);
            float P_i_dir = exp(-r * T) * max(0., K - S_i_dir);
            sum_parz_dir_C += C_i_dir;
            sum_parz_dir_P += P_i_dir;

            // Simulazione per opzioni di tipo Call e Put con metodo discreto
            float S_i_disc = 100;
            for (int k = 0; k < H; k++) {
                S_i_disc = S_i_disc * exp((r - sigma * sigma / 2) * T / float(H) + sigma * rnd.Gauss(0, 1) * sqrt(T / float(H)));
            }
            float C_i_disc = exp(-r * T) * max(0., S_i_disc - K);
            float P_i_disc = exp(-r * T) * max(0., K - S_i_disc);
            sum_parz_disc_C += C_i_disc;
            sum_parz_disc_P += P_i_disc;
        }

        // Aggiornamento degli accumulatori totali per i risultati con metodo diretto
        sum_dir_C += sum_parz_dir_C / L;
        sum_dir_quad_C += pow(sum_parz_dir_C / L, 2);
        sum_dir_P += sum_parz_dir_P / L;
        sum_dir_quad_P += pow(sum_parz_dir_P / L, 2);

        // Aggiornamento degli accumulatori totali per i risultati con metodo discreto
        sum_disc_C += sum_parz_disc_C / L;
        sum_disc_quad_C += pow(sum_parz_disc_C / L, 2);
        sum_disc_P += sum_parz_disc_P / L;
        sum_disc_quad_P += pow(sum_parz_disc_P / L, 2);

        // Scrittura dei risultati su file
        flusso_out << i + 1 << " " << sum_dir_C / (i + 1) << " " << err(sum_dir_C / (i + 1), sum_dir_quad_C / (i + 1), i + 1) << " "
                   << sum_dir_P / (i + 1) << " " << err(sum_dir_P / (i + 1), sum_dir_quad_P / (i + 1), i + 1) << " "
                   << sum_disc_C / (i + 1) << " " << err(sum_disc_C / (i + 1), sum_disc_quad_C / (i + 1), i + 1) << " "
                   << sum_disc_P / (i + 1) << " " << err(sum_disc_P / (i + 1), sum_disc_quad_P / (i + 1), i + 1) << endl;
    }

    flusso_out.close();  // Chiude il file di output
    return 0;  // Fine del programma
}

// Definizione della funzione per calcolare il massimo tra due numeri float
float max(float x1, float x2) {
    if (x1 > x2)
        return x1;
    return x2;
}

// Definizione della funzione per calcolare l'errore standard della media
double err(double media, double mediaquad, int n) {
    if (n == 1)
        return 0;
    else {
        return sqrt((mediaquad - media * media) / (n - 1));
    }
}
