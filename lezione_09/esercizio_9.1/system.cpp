#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h" 
#include <vector>
#include <unordered_set>

using namespace std;
using namespace arma;

/////////////////////////////////////////
// Classe path
/////////////////////////////////////////

void path::swap(int c1, int c2){
  int appo = _path[c1];  // Salva temporaneamente l'elemento in posizione c1
  _path[c1] = _path[c2];  // Sposta l'elemento in posizione c2 a c1
  _path[c2] = appo;       // Sposta l'elemento salvato temporaneamente in c2
}

void path::check(){
  if (_path.size() != _ncities) {  // Controlla se la lunghezza del percorso è diversa da quella attesa
    throw invalid_argument("Vettore di lunghezza diversa da quella spettata.");
  }

  if(_path[0] != 0){  // Controlla se la prima città nel percorso non è quella corretta (la città 0)
    throw invalid_argument("Città iniziale non giusta.");
  }

  unordered_set<int> labels;  // Utilizza un set non ordinato per controllare la ripetizione delle città

  for (int num = 0; num < _ncities; num++) {
    if (!labels.insert(_path[num]).second){  // Inserisce nel set e controlla se c'è già la città nel percorso
      throw invalid_argument("Città ripetuta nel percorso.");
    }
    if(_path[num] < 0 || _path[num] > _ncities-1){  // Controlla se l'indice della città è valido
      throw invalid_argument("Indice città non valido.");
    }
  }
}

double path::length_L1(const vector<city> &cities) const {
  double sum = 0;
  for(int i = 0; i < _ncities-1; i++){
    sum += cities[_path[i]].distance_L1(cities[_path[i+1]]);  // Calcola la lunghezza totale utilizzando la distanza L1 tra città consecutive
  }
  sum += cities[_path[_ncities-1]].distance_L1(cities[_path[0]]);  // Aggiunge la distanza dall'ultima città alla prima

  return sum;
}

double path::length_L2(const vector<city> &cities) const {
  double sum = 0;
  for(int i = 0; i < _ncities-1; i++){
    sum += cities[_path[i]].distance_L2(cities[_path[i+1]]);  // Calcola la lunghezza totale utilizzando la distanza L2 tra città consecutive
  }
  sum += cities[_path[_ncities-1]].distance_L2(cities[_path[0]]);  // Aggiunge la distanza dall'ultima città alla prima

  return sum;
}

void path::mut_pairperm(){
  int pos1 = int((_ncities-1)*_rnd->Rannyu()) + 1;  // Seleziona casualmente una posizione nel percorso, escludendo la prima città
  int pos2 = int((_ncities-1)*_rnd->Rannyu()) + 1;  // Seleziona casualmente un'altra posizione nel percorso, escludendo la prima città

  if(pos1 == 0 || pos2 == 0){
    cout << "stai facendo lo swap con la prima posizione" << endl;  // Stampa un messaggio se si sta facendo lo swap con la prima posizione
  }
  swap(pos1, pos2);  // Esegue lo swap delle due posizioni nel percorso
}

void path::mut_shift(){
  int m = int((_ncities-1)*_rnd->Rannyu()) + 1;  // Seleziona casualmente una posizione nel percorso, escludendo la prima città
  int n = int(_rnd->Rannyu(m, _ncities));  // Seleziona casualmente un'altra posizione nel percorso, tra m e l'ultima città

  if (m < 1 || n >= _ncities) {
    cerr << "Indice non valido o conteggio troppo grande!" << endl;  // Stampa un messaggio di errore se l'indice non è valido
    return;
  }

  // Shift degli elementi dal m-esimo al n-esimo nel percorso
  rotate(_path.begin() + 1, _path.begin() + m, _path.begin() + n + 1);
}

void path::mut_permute() {
    int m = int((_ncities / 2 - 1) * _rnd->Rannyu()) + 1;  // Genera m tale che m < N/2

    // Genera un indice di inizio tale che i due blocchi contigui di lunghezza m non includano l'elemento zero
    int start1 = int((_ncities - 2 * m - 1) * _rnd->Rannyu()) + 1;

    // Indice iniziale del secondo blocco contiguo
    int start2 = start1 + m;

    // Scambia gli elementi dei due blocchi contigui di lunghezza m
    std::swap_ranges(_path.begin() + start1, _path.begin() + start1 + m, _path.begin() + start2);
}

void path::mut_inversion(){
    // Genera due indici casuali per delimitare il sottovettore da invertire
    int start = int((_ncities - 1) * _rnd->Rannyu()) + 1;
    int end = int((_ncities - 1) * _rnd->Rannyu()) + 1;

    // Assicurati che start sia minore di end
    if (start > end) {
      int appo = end;
      end = start;
      start = appo;
    }

    // Inverti il sottovettore
    std::reverse(_path.begin() + start, _path.begin() + end + 1);
}

void path::update_path(const vector<int>& new_path) {
    _path = new_path;
}

void path::crossover(path* my_path) {
    int cut = int((_ncities - 1) * _rnd->Rannyu()) + 1;

    vector<int> parent1 = _path;
    vector<int> parent2 = my_path->get_vec();

    vector<int> child1 = parent1;
    vector<int> child2 = parent2;

    // Inizializza i figli con la prima parte dei rispettivi genitori
    for (int i = cut; i < _ncities; ++i) {
        child1[i] = 0;
        child2[i] = 0;
    }

    // Riempie la seconda parte di child1 con le città di parent2
    int index1 = cut;
    for (int i = 0; i < _ncities; ++i) {
        if (find(child1.begin(), child1.end(), parent2[i]) == child1.end()) {
            child1[index1++] = parent2[i];
        }
    }

    // Riempie la seconda parte di child2 con le città di parent1
    int index2 = cut;
    for (int i = 0; i < _ncities; ++i) {
        if (find(child2.begin(), child2.end(), parent1[i]) == child2.end()) {
            child2[index2++] = parent1[i];
        }
    }

    // Aggiorna i percorsi con i nuovi figli
    _path = child1;
    my_path->update_path(child2);
}


/////////////////////////////////////////
// Classe system
/////////////////////////////////////////

void System::initialize(){
  for(int i = 0; i < _npop; i++){
    for(int j = 0; j < 100*_ncities_sist; j++){
      _syst[i].mut_pairperm();
    }
    _syst[i].check();
  }
}

void System::stamp_path(ofstream &fout){
  fout << "Path:" << endl;
  for(int i = 0; i < _npop; i++){
    for(int j = 0; j < _ncities_sist; j++){
      fout << _syst[i].get_pos(j) << "  " ;
    }
    fout << endl;
  }
}

void System::sortL1(const vector<city> &cities){
  sort(_syst.begin(), _syst.end(), [&cities](const path &a, const path &b) {
    return a.length_L1(cities) < b.length_L1(cities);
  });
}

void System::sortL2(const vector<city> &cities){
  sort(_syst.begin(), _syst.end(), [&cities](const path &a, const path &b) {
    return a.length_L2(cities) < b.length_L2(cities);
  });
}

void System::mutation(double p1, double p2, double p3, double p4){
  if (p1 < 0.0 || p1 > 1.0 || p2 < 0.0 || p2 > 1.0 || p3 < 0.0 || p3 > 1.0 || p4 < 0.0 || p4 > 1.0) {
    throw invalid_argument("Le probabilità di mutazione devono essere comprese tra zero e uno.");    
  }

  for(int i = 0; i < _son.size(); i++){
    if(_rnd->Rannyu() < p1) _syst[_son[i]].mut_pairperm();  // Applica la mutazione di scambio con probabilità p1
    if(_rnd->Rannyu() < p2) _syst[_son[i]].mut_shift();     // Applica la mutazione di shift con probabilità p2
    if(_rnd->Rannyu() < p3) _syst[_son[i]].mut_permute();   // Applica la mutazione di permutazione con probabilità p3
    if(_rnd->Rannyu() < p4) _syst[_son[i]].mut_inversion(); // Applica la mutazione di inversione con probabilità p4
  }
}

void System::stamp_city(const vector<city> &cities, ofstream &fout) const {

  if (cities.size() != _ncities_sist) {
    throw invalid_argument("Le città non sono in numero giusto.");    
  }

  fout << "    x      y" << endl;
  for(int i = 0; i < _ncities_sist; i++){
    fout << " " << cities[i].getx() << "  " << cities[i].gety() << endl;  // Stampa le coordinate x e y di tutte le città
  }
}

void System::check(){
  for(int i = 0; i < _npop; i++){
    _syst[i].check();  // Esegue il controllo di validità per ogni percorso nel sistema
  }
}

void System::crossover(double p, double pot) {
    int c = 0;
    int target = round(p * _npop);

    while (c < target) {
        int pa1;
        int pa2;

        // Trova un path1 unico
        bool unique1 = false;
        while (!unique1) {
            pa1 = this->select(pot);  // Seleziona un genitore potenziato casualmente
            unique1 = true;
            for (int k = 0; k < _son.size(); k++) {
                if (pa1 == _son[k]) {
                    unique1 = false;
                    break;
                }
            }
        }
        _son.push_back(pa1);

        // Trova un path2 unico
        bool unique2 = false;
        while (!unique2) {
            pa2 = this->select(pot);  // Seleziona un altro genitore potenziato casualmente
            unique2 = true;
            for (int k = 0; k < _son.size(); k++) {
                if (pa2 == _son[k]) {
                    unique2 = false;
                    break;
                }
            }
        }
        _son.push_back(pa2);

        // Esegue l'operazione di crossover
        _syst[pa1].crossover(&_syst[pa2]);
        c += 2;
    }
}

void System::crossover2(double p, double pot) {
    int c = 0;
    int target = round(p * _npop);

    while (c < target) {
        int pa1;
        int pa2;

        // Trova un path1 unico
        bool unique1 = false;
        while (!unique1) {
            pa1 = this->select(pot);  // Seleziona un genitore potenziato casualmente
            unique1 = true;
            for (int k = 0; k < _son.size(); k++) {
                if (pa1 == _son[k]) {
                    unique1 = false;
                    break;
                }
            }
        }
        _parents.push_back(pa1);

        // Trova un path2 unico
        bool unique2 = false;
        while (!unique2) {
            pa2 = this->select(pot);  // Seleziona un altro genitore potenziato casualmente
            unique2 = true;
            for (int k = 0; k < _son.size(); k++) {
                if (pa2 == _son[k]) {
                    unique2 = false;
                    break;
                }
            }
        }
        _parents.push_back(pa2);

        // Esegue l'operazione di crossover
        path child1 = _syst[pa1];
        path child2 = _syst[pa2];
        child1.crossover(&child2);

        int a1 = this->select(1/pot);  // Seleziona un antenato casualmente potenziato
        int a2 = this->select(1/pot);  // Seleziona un altro antenato casualmente potenziato

        _son.push_back(a1);
        _son.push_back(a2);

        _syst[a1] = child1;
        _syst[a2] = child1;

        c += 2;
    }

}

void System::restart(){
  vector<int> v;
  _son = v;
}

void System::printfirstpath(ofstream &fout, const vector<city> &cities) const {

  for(int i = 0; i < _ncities_sist; i++){
    fout << " " << _syst[0].get_pos(i);  // Stampa il primo percorso nel sistema
  }

  fout << "       " << _syst[0].length_L2(cities) << endl;  // Stampa la lunghezza del primo percorso utilizzando la distanza L2
}

void System::printmean(ofstream &fout, const vector<city> &cities, int N_mean, int d, int generation){

  if (d != 0 & d != 1) {
    throw invalid_argument("L'indice per la stampa della distanza deve essere: 0-L1, 1-L2.");    
  }

  double sum = 0;
  if(d == 0){
    for(int i = 0; i < N_mean; i++){
      sum += _syst[i].length_L1(cities);  // Calcola la somma delle lunghezze dei percorsi utilizzando la distanza L1
    }
  }
  if(d == 1){
    for(int i = 0; i < N_mean; i++){
      sum += _syst[i].length_L2(cities);  // Calcola la somma delle lunghezze dei percorsi utilizzando la distanza L2
    }
  }

  fout << generation << "   " << sum/N_mean << endl;  // Stampa la media delle lunghezze dei percorsi per la generazione corrente
}


/////////////////////////////////////////
// Funzioni varie
/////////////////////////////////////////


void printProgressBar(int current, int total) {
    int barWidth = 70; // Larghezza della barra di avanzamento
    float progress = (float)current / total;
    int pos = barWidth * progress;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";  // Segna il progresso fino alla posizione corrente
        else if (i == pos) std::cout << ">";  // Indica la posizione corrente sulla barra di avanzamento
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";  // Stampa percentuale di avanzamento
    std::cout.flush();
}

void start(int* type, int* n_city, int* n_iter, int* n_path){

  ifstream input("input.dat");  // Apre il file di input
  string property;
  *type = -1;
  *n_city = -1;
  *n_iter = -1;
  *n_path = -1;

  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      string sim_type;
      input >> sim_type;
      bool check = false;
      if(sim_type == "circle"){
        *type = 0;  // Imposta il tipo di simulazione a cerchio
        check = true;
      }
      if(sim_type == "square"){
        *type = 1;  // Imposta il tipo di simulazione a quadrato
        check = true;
      }
      if(sim_type == "file"){
        *type = 2;  // Legge le città da un file di input
        check = true;
      }
      if(check == false){
        throw invalid_argument("Invalid argument SIMULATION TYPE.");  // Tipo di simulazione non valido
      }
    } else if( property == "N_CITY" ){
      input >> *n_city;  // Legge il numero di città dal file di input
    } else if( property == "N_ITERATION" ){
      input >> *n_iter;  // Legge il numero di iterazioni dal file di input
    } else if( property == "N_PATH" ){
      input >> *n_path;  // Legge il numero di percorsi dal file di input
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  if(*type == -1){
    throw invalid_argument("Missing SIMULATION TYPE input.");  // Tipo di simulazione non specificato nel file di input
  }
  if(*n_city < 0){
    throw invalid_argument("Wrong N_CITY input.");  // Numero di città non valido nel file di input
  }
  if(*n_iter < 0){
    throw invalid_argument("Wrong N_ITERATION input.");  // Numero di iterazioni non valido nel file di input
  }
  if(*n_path < 0){
    throw invalid_argument("Wrong N_PATH input.");  // Numero di percorsi non valido nel file di input
  }
  input.close();  // Chiude il file di input
  return;
}

