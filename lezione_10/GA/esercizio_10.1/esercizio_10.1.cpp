#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "mpi.h"  // Libreria MPI per la comunicazione tra processi
#include "../../../rangen/random.h"  // Libreria per la generazione di numeri casuali
#include "system.h"  // Header del sistema che gestisce il problema

using namespace std;

int main(int argc, char *argv[]) {

    int size, rank;

    MPI_Init(&argc, &argv);  // Inizializzazione MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Numero totale di processi MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Rank del processo corrente

    // Controllo per non usare più di 8 cores (rank va da 0 a size-1)
    if (rank > 8) throw invalid_argument("Non usare più di 7 cores.");

    Random rnd = generaRand(rank + 1);  // Inizializzazione del generatore di numeri casuali

    vector<city> cities;  // Vettore di città

    int N_city;  // Numero di città
    int N_iter;  // Numero di iterazioni
    int N_path;  // Numero di percorsi
    int type;    // Tipo (non specificato nel codice fornito)
    bool paral;  // Flag per il tipo di esecuzione (parallela o singola)

    // Funzione non definita nel codice fornito, presa come inizializzazione del sistema
    start(&type, &N_city, &N_iter, &N_path, &paral);

    // Lettura delle coordinate delle città da file
    ifstream infile("cap_prov_ita.dat");
    double x, y;
    while (infile >> x >> y) {
        cities.push_back(city(x, y));  // Creazione delle città nel vettore
    }
    infile.close();

    // Scelta della cartella di output in base al tipo di esecuzione (parallela o singola)
    string name_folder;
    if (paral)
        name_folder = "PARALL";
    else
        name_folder = "SINGLE";

    // Creazione dei nomi dei file di output specifici per il processo corrente
    string starting_pop = "OUTPUT_" + name_folder + "/P" + to_string(rank + 1) + "/starting_pop.dat";
    string final_pop = "OUTPUT_" + name_folder + "/P" + to_string(rank + 1) + "/final_pop.dat";
    string mean = "OUTPUT_" + name_folder + "/P" + to_string(rank + 1) + "/mean.dat";
    string firstpath = "OUTPUT_" + name_folder + "/P" + to_string(rank + 1) + "/firstpath.dat";
    string best = "OUTPUT_" + name_folder + "/best";

    // Apertura dei file di output per il processo corrente
    ofstream f1out, f2out, ffout, fmout, fbout;
    f1out.open(starting_pop);
    f2out.open(final_pop);
    fmout.open(mean);
    ffout.open(firstpath);
    fbout.open(best);

    // Creazione di un vettore di città ordinate per indice
    vector<int> v;
    for (int i = 0; i < N_city; i++) v.push_back(i);

    // Creazione di un percorso casuale iniziale usando il generatore di numeri casuali
    path my_path(v, &rnd);
    my_path.check();  // Controllo dell'integrità del percorso

    // Creazione del sistema con il percorso iniziale, il generatore di numeri casuali e il numero di percorsi desiderato
    System my_syst(my_path, &rnd, N_path);

    my_syst.initialize();  // Inizializzazione del sistema
    my_syst.check();       // Controllo dell'integrità del sistema

    my_syst.stamp_path(f1out);  // Scrittura del percorso iniziale nel file di output

    // Ordinamento delle città nel sistema rispetto a una certa metrica (presumo)
    my_syst.sortL1(cities);

    // Probabilità di migrazione tra processi (attivata solo se paral è true)
    double p_migr = 0.05;
    if (!paral) p_migr = 0;

    // Loop principale delle iterazioni
    for (int i = 0; i < N_iter; i++) {
        my_syst.crossover2(0.8, 2.2);  // Operazione di crossover sui percorsi
        my_syst.mutation(0.25, 0.1, 0.1, 0.2);  // Operazione di mutazione sui percorsi
        my_syst.sortL2(cities);  // Riordinamento dei percorsi secondo un'altra metrica
        my_syst.check();  // Controllo dell'integrità del sistema dopo le operazioni
        my_syst.migration(rank, size, p_migr);  // Operazione di migrazione tra processi
        my_syst.restart();  // Riavvio del sistema per la prossima iterazione
        my_syst.printfirstpath(ffout, cities);  // Scrittura del primo percorso nel file di output
        my_syst.printmean(fmout, cities, N_path / 2, 1, i);  // Calcolo e scrittura della media nel file di output

        if (rank == 0) printProgressBar(i + 1, N_iter);  // Stampa della barra di avanzamento (presumo)
    }

    if (rank == 0) cout << endl;  // Fine della stampa della barra di avanzamento

    my_syst.sortL2(cities);  // Riordinamento finale dei percorsi
    my_syst.stamp_path(f2out);  // Scrittura del percorso finale nel file di output
    my_syst.stamp_best(fbout, rank, size, cities);  // Scrittura del miglior percorso nel file di output

    // Chiusura di tutti i file di output
    f1out.close();
    f2out.close();
    fmout.close();
    ffout.close();
    fbout.close();

    MPI_Finalize();  // Finalizzazione MPI

    return 0;
}
