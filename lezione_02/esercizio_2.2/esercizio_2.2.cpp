#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h" // Includi la libreria per la generazione di numeri casuali
#include "libes_2.2.h" // Includi il file di intestazione per le funzioni definite in "libes_2.2.h"

using namespace std;

int main(int argc, char* argv[]) {
    
    int usr = -1;

    // Chiedi all'utente di scegliere il tipo di Random Walk
    cout << "Premere 0 per avere Random Walk a diversa lunghezza correlati" << endl;
    cout << "Premere 1 per avere Random Walk a diversa lunghezza scorrelati" << endl;
    cin >> usr;

    // Verifica se l'input è valido
    if(usr == 0 || usr == 1){}
    else{
        cout << "Numero inserito non corretto" << endl;
        return -1;
    }

    string namefile;
    // Scegli il nome del file in base alla scelta dell'utente
    if(usr == 0){
        namefile = "risultati_corr.dat";
    }
    else{
        namefile = "risultati_scorr.dat";
    }
    // Apre il file di output
    ofstream flusso_out;
    flusso_out.open(namefile);

    int M = 1000000;  // numero di iterazioni totali
    int N = 100; // numero blocchi 
    int L = M / N;  // lunghezza blocco
    int lunghezza = 1E2;  // lunghezza del Random Walk

    Random rnd = generaRand(); // Inizializza il generatore di numeri casuali

    vector<double> somma_distanza_continua(lunghezza, 0); // Vettore per la somma delle distanze nel caso continuo
    vector<double> sommaquad_distanza_continua(lunghezza, 0); // Vettore per la somma dei quadrati delle distanze nel caso continuo

    vector<double> somma_distanza_discreta(lunghezza, 0); // Vettore per la somma delle distanze nel caso discreto
    vector<double> sommaquad_distanza_discreta(lunghezza, 0); // Vettore per la somma dei quadrati delle distanze nel caso discreto
    
    // Ciclo sui blocchi
    for(int i = 0; i < N; i++) {
        vector<double> distanza_continua(lunghezza, 0); // Vettore per il calcolo delle distanze nel caso continuo
        vector<double> distanza_discreta(lunghezza, 0); // Vettore per il calcolo delle distanze nel caso discreto

        vector<vector<double>> rw_continuo; // Vettore per il Random Walk nel caso continuo
        vector<vector<double>> rw_discreto; // Vettore per il Random Walk nel caso discreto
        // Ciclo sugli step all'interno di un blocco
        for(int j = 0; j < L; j++) {
            // Genera il Random Walk a seconda del tipo scelto
            if(usr == 0){
                rw_continuo = generaRWc_corr(&rnd, lunghezza);
                rw_discreto = generaRWd_corr(&rnd, lunghezza);
            }
            else{
                rw_continuo = generaRWc_scorr(&rnd, lunghezza);
                rw_discreto = generaRWd_scorr(&rnd, lunghezza);
            }
            // Calcola la somma delle distanze dei vari RW del blocco
            for(int k = 0; k < lunghezza; k++) {
                distanza_continua[k] += moduloquad(&rw_continuo)[k];
                distanza_discreta[k] += moduloquad(&rw_discreto)[k]; 
            }
        }

        // Calcola la media delle distanze e il loro quadrato per entrambi i casi
        for(int k = 0; k < lunghezza; k++) { 
            distanza_continua[k] = sqrt(distanza_continua[k]/L);
            distanza_discreta[k] = sqrt(distanza_discreta[k]/L);
        }

        // Aggiorna le somme e le somme dei quadrati delle distanze
        for(int k = 0; k < lunghezza; k++) {
            somma_distanza_continua[k] += distanza_continua[k];
            sommaquad_distanza_continua[k] += distanza_continua[k]*distanza_continua[k];

            somma_distanza_discreta[k] += distanza_discreta[k];
            sommaquad_distanza_discreta[k] += distanza_discreta[k]*distanza_discreta[k];
        } 
    }

    // Calcola la media e l'errore delle medie per entrambi i casi e stampa i risultati su file
    for(int k = 0; k < lunghezza; k++) {
        flusso_out << k << " " << somma_distanza_discreta[k]/N << " " << errore(somma_distanza_discreta[k], sommaquad_distanza_discreta[k], N, k) 
                        << " " << somma_distanza_continua[k]/N << " " << errore(somma_distanza_continua[k], sommaquad_distanza_continua[k], N, k) << endl;
    }

    // Chiude il file di output
    flusso_out.close();

    cout << "Risultati stampati in: " << namefile << endl;

    return 0;
}
