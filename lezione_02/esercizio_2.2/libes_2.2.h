#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include "../../rangen/random.h"

using namespace std;


///////////////////////////////////////////////////////////////////////////////
// Calcola il vettore somma di due vettori
///////////////////////////////////////////////////////////////////////////////
vector<double> somma(vector<double>* v1, vector<double>* v2) {
    // Assicurati che i due vettori abbiano la stessa dimensione
    if (v1->size() != v2->size()) {
        cerr << "Errore: i due vettori hanno dimensioni diverse." << endl;
        return {};
    }

    vector<double> v3;
    v3.reserve(v1->size()); // Riserva la memoria per evitare riallocazioni

    // Somma elemento per elemento
    for(size_t i = 0; i < v1->size(); ++i) {
        v3.push_back((*v1)[i] + (*v2)[i]);
    }

    return v3;
}

///////////////////////////////////////////////////////////////////////////////
// Effettua la divisione tra vettore e scalare
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>>* divide(vector<vector<double>>* v, int L) {
    double Ld = static_cast<double>(L);
    for(size_t j = 0; j < v->size(); j++) {
        for(int i = 0; i < 3; i++) {
            (*v)[j][i] = (*v)[j][i] / Ld;
        }
    }
    return v;
}

///////////////////////////////////////////////////////////////////////////////
// Calcola il vettore dei moduli dato un vettore di posizioni (x,y,z)
///////////////////////////////////////////////////////////////////////////////
vector<double> moduloquad(vector<vector<double>>* v) {
    vector<double> v2;
    for(size_t k = 0; k < v->size(); k++) {
        double somma = 0;
        for(int i = 0; i < 3; i++){
            somma += ((*v)[k][i]) * ((*v)[k][i]);
        }
        v2.push_back(somma);
    }
    return v2;
}

///////////////////////////////////////////////////////////////////////////////
// Genera -1 o 1 casualmente
///////////////////////////////////////////////////////////////////////////////
int passo(Random* rnd) {
    float rand = rnd->Rannyu();
    return 2*int(2 * rand) - 1;
}

///////////////////////////////////////////////////////////////////////////////
// Genera un angolo polare theta in tre dimensioni
///////////////////////////////////////////////////////////////////////////////
double theta(Random* rnd) {
    float rand = rnd->Rannyu();
    return acos(1-2*rand);
}

///////////////////////////////////////////////////////////////////////////////
// Genera un Random Walk (caso discreto) in cui la posizione i-esima è ottenuta evolvendo la (i-1)-esima
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> generaRWd_corr(Random* rnd, int lungh) {
    vector<double> coord(3, 0);
    vector<vector<double>> v;  // Vettore contenente le posizioni (x,y,z) al variare dello step
    v.push_back(coord);  // Parto da (0,0,0)
    for(int j = 0; j < lungh; j++) {
        vector<double> appo(3,0);
        int i = int(3 * rnd->Rannyu()); // Decido la componente x, y o z che si muove di +1 o -1
        appo[i]=passo(rnd);
        coord = somma(&coord, &appo);  // Sommo al RW precedente il passo successivo
        v.push_back(coord);  // Inserisco nel RW la nuova posizione
    }
    return v;
}

///////////////////////////////////////////////////////////////////////////////
// Genera un Random Walk (caso discreto); i RW a differente lunghezza sono indipendenti
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> generaRWd_scorr(Random* rnd, int lungh) {
    vector<vector<double>> v;  // Vettore contenente le posizioni (x,y,z) al variare dello step

    for(int k = 0; k < lungh; k++){  // Per ogni step del RW riparto da capo ovvero da (0,0,0)
        vector<double> coord(3, 0);  // Parto da (0,0,0)
        for(int j = 0; j < k; j++) {
            vector<double> appo(3,0);
            int i = int(3 * rnd->Rannyu());  // Decido la componente x, y o z che si muove di +1 o -1
            appo[i]=passo(rnd);
            coord = somma(&coord, &appo);
        }
        v.push_back(coord);  // Inserisco la posizione finale del nuovo RW
    }
    return v; 
}

///////////////////////////////////////////////////////////////////////////////
// Genera un Random Walk (caso continuo) in cui la posizione i-esima è ottenuta evolvendo la (i-1)-esima
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> generaRWc_corr(Random* rnd, int lungh) {
    vector<double> coord(3, 0);
    vector<vector<double>> v;  // Vettore contenente le posizioni (x,y,z) al variare dello step
    v.push_back(coord);  // Parto da (0,0,0)
    for(int j = 0; j < lungh; j++) {
        vector<double> appo;
        double thet = theta(rnd);  // Estraggo l'angolo polare
        double phi = rnd->Rannyu(0,2*M_PI);  // Estraggo l'angolo phi
        appo.push_back(sin(thet)*cos(phi));  // Aggiorno la componente x
        appo.push_back(sin(thet)*sin(phi));  // Aggiorno la componente y
        appo.push_back(cos(thet));  // Aggiorno la componente z
        coord = somma(&coord, &appo);  // Sommo al RW precedente il passo successivo
        v.push_back(coord);  // Inserisco nel RW la nuova posizione
    }
    return v;
}

///////////////////////////////////////////////////////////////////////////////
// Genera un Random Walk (caso continuo); i RW a differente lunghezza sono indipendenti
///////////////////////////////////////////////////////////////////////////////
vector<vector<double>> generaRWc_scorr(Random* rnd, int lungh) {
    vector<vector<double>> v;  // Vettore contenente le posizioni (x,y,z) al variare dello step

    for(int k = 0; k < lungh; k++){
        vector<double> coord(3, 0);   // Parto da (0,0,0)
        for(int j = 0; j < k; j++) {  // Per ogni step del RW riparto da capo ovvero da (0,0,0)
            vector<double> appo;
            double thet = theta(rnd);  // Estraggo l'angolo polare
            double phi = rnd->Rannyu(0,2*M_PI);  // Estraggo l'angolo phi
            appo.push_back(sin(thet)*cos(phi));  // Aggiorno la componente x
            appo.push_back(sin(thet)*sin(phi));  // Aggiorno la componente y
            appo.push_back(cos(thet));  // Aggiorno la componente z
            coord = somma(&coord, &appo);  // Sommo al RW il passo successivo
        }
        v.push_back(coord);
    }
    return v;
}

///////////////////////////////////////////////////////////////////////////////
// Calcola la deviazione standard
///////////////////////////////////////////////////////////////////////////////
double errore(double somma, double sommaquad, int N, int passi){
    if(passi == 0 || passi == 1) return 0;
    else{
        return sqrt(1./(N-1)*(sommaquad/N-pow(somma/N,2)));
    }


}