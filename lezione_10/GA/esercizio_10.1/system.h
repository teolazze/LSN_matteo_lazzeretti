#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "../../../rangen/random.h"
#include "mpi.h"

using namespace std;
using namespace arma;

// Dichiarazione della classe "city" che rappresenta una città
class city {
private:
    double _x, _y;

public:
    city(double x, double y) : _x(x), _y(y) {}
    double getx() const { return _x; }
    double gety() const { return _y; }
    double distance_L1(const city& c) const {
        return abs(c.getx() - _x) + abs(c.gety() - _y);
    }
    double distance_L2(const city& c) const {
        return sqrt(pow(c.getx() - _x, 2) + pow(c.gety() - _y, 2));
    }
};

// Dichiarazione della classe "path" che rappresenta un percorso
class path {
private:
    vector<int> _path;
    int _ncities;
    Random* _rnd;          // Generatore di numeri casuali

    void swap(int c1, int c2);

public:
    // Costruttore
    path(vector<int> path, Random* rnd) : _path(path), _rnd(rnd), _ncities(path.size()) {}

    // Metodo per verificare la correttezza del percorso
    void check();

    // Metodi per accedere ai membri della classe
    int get_ncities() { return _ncities; }
    int get_pos(int k) const { return _path[k]; }
    vector<int> get_vec() const { return _path; }
    void write_pos(int k, int city) { _path[k] = city; }

    // Metodo per aggiornare il percorso
    void update_path(const vector<int>& new_path);

    // Metodi per calcolare la lunghezza del percorso
    double length_L1(const vector<city>& cities) const;
    double length_L2(const vector<city>& cities) const;

    // Metodi per le mutazioni del percorso
    void mut_pairperm();
    void mut_shift();
    void mut_permute();
    void mut_inversion();

    // Metodo per il crossover tra due percorsi
    void crossover(path* my_path);
};

// Dichiarazione della classe "System" che gestisce il sistema di percorsi
class System {
private:
    vector<path> _syst;
    vector<int> _son;
    vector<int> _parents;
    int _npop;
    int _ncities_sist;
    Random* _rnd;

public:
    // Costruttore
    System(path start_path, Random* rnd, int npop) : _npop(npop), _rnd(rnd), _ncities_sist(start_path.get_ncities()) {
        for (int i = 0; i < npop; i++) _syst.push_back(start_path);
    }

    // Metodo per inizializzare il sistema
    void initialize();

    // Metodi per stampare le informazioni del sistema
    void stamp_path(ofstream& fout);
    void stamp_best(ofstream &fout, int rank, int size, const vector<city>& cities);
    void sortL1(const vector<city>& cities);
    void sortL2(const vector<city>& cities);
    void stamp_city(const vector<city>& cities, ofstream& fout) const;
    void printfirstpath(ofstream& fout, const vector<city> &cities) const;
    void printmean(ofstream& fout, const vector<city>& cities, int N_mean, int d, int generation);
    void check();

    // Metodi per le operazioni di mutazione e crossover
    void mutation(double p1, double p2, double p3, double p4);
    void crossover(double p, double pot);
    void crossover2(double p, double pot);
    void restart();

    void migration(int rank, int size, double p);

    // Metodo per selezionare un elemento in base a una legge di potenza
    int select(double p){return int((_npop)*pow(_rnd->Rannyu(),p));};
};

// Funzione per leggere le impostazioni iniziali dal file "input.dat"
void start(int* type, int* n_city, int* n_iter, int* n_path, bool* paral);

// Funzione per stampare una barra di avanzamento
void printProgressBar(int current, int total);

#endif // __System__