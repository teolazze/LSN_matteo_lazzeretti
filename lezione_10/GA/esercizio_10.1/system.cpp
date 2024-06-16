#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"
#include <vector>
#include <unordered_set>
#include "mpi.h"

using namespace std;
using namespace arma;

// Swappa due città
void path::swap(int c1, int c2){
  int appo = _path[c1];
  _path[c1] = _path[c2];
  _path[c2] = appo;
}

// Check su correttezza del path
void path::check(){
  if (_path.size() != _ncities) {
    throw invalid_argument("Vettore di lunghezza diversa da quella spettata.");
  }

  if(_path[0] != 0){
    throw invalid_argument("Città iniziale non giusta.");
  }

  unordered_set<int> labels;

  for (int num = 0; num < _ncities; num++) {
    if (!labels.insert(_path[num]).second){
      throw invalid_argument("Città ripetuta nel percorso.");
    }
    if(_path[num] < 0 || _path[num] > _ncities-1){
      throw invalid_argument("Indice città non valido.");
    }
  }
}

// Calcola lunghezza L1 dato il path
double path::length_L1(const vector<city> &cities) const {
  double sum = 0;
  for(int i = 0; i < _ncities-1; i++){
    sum += cities[_path[i]].distance_L1(cities[_path[i+1]]);
  }
  sum += cities[_path[_ncities-1]].distance_L1(cities[_path[0]]);

  return sum;
}

// Calcola lunghezza L2 dato il path
double path::length_L2(const vector<city> &cities) const {
  double sum = 0;
  for(int i = 0; i < _ncities-1; i++){
    sum += cities[_path[i]].distance_L2(cities[_path[i+1]]);
  }
  sum += cities[_path[_ncities-1]].distance_L2(cities[_path[0]]);

  return sum;
}

// Mutazione pair permutation
void path::mut_pairperm(){
  int pos1 = int((_ncities-1)*_rnd->Rannyu()) + 1;
  int pos2 = int((_ncities-1)*_rnd->Rannyu()) + 1;

  if(pos1 == 0 || pos2 == 0){
    cout << "stai facendo lo swap con la prima posizione" << endl;
  }
  swap(pos1, pos2);
}

// Mutazione shift
void path::mut_shift(){
  int m = int((_ncities-1)*_rnd->Rannyu()) + 1;
  int n = int(_rnd->Rannyu(m, _ncities));

  if (m < 1 || n >= _ncities) {
    cerr << "Indice non valido o conteggio troppo grande!" << endl;
    return;
  }
  rotate(_path.begin() + 1, _path.begin() + m, _path.begin() + n + 1);
}

// Mutazione permute
void path::mut_permute() {
  int m = int((_ncities / 2 - 1) * _rnd->Rannyu()) + 1; // Genera m tale che m < N/2

  int start1 = int((_ncities - 2 * m - 1) * _rnd->Rannyu()) + 1;

  int start2 = start1 + m;

  std::swap_ranges(_path.begin() + start1, _path.begin() + start1 + m, _path.begin() + start2);
}

// Mutazione inversione
void path::mut_inversion(){
  int start = int((_ncities - 1) * _rnd->Rannyu()) + 1;
  int end = int((_ncities - 1) * _rnd->Rannyu()) + 1;

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

// Funzione di crossover
void path::crossover(path* my_path) {
  int cut = int((_ncities - 1) * _rnd->Rannyu()) + 1;


  vector<int> parent1 = _path;
  vector<int> parent2 = my_path->get_vec();

  vector<int> child1 = parent1;
  vector<int> child2 = parent2;

  for (int i = cut; i < _ncities; ++i) {
      child1[i] = 0;
      child2[i] = 0;
  }

  int index1 = cut;
  for (int i = 0; i < _ncities; ++i) {
    if (find(child1.begin(), child1.end(), parent2[i]) == child1.end()) {
      child1[index1++] = parent2[i];
    }
  }

  int index2 = cut;
  for (int i = 0; i < _ncities; ++i) {
    if (find(child2.begin(), child2.end(), parent1[i]) == child2.end()) {
      child2[index2++] = parent1[i];
    }
  }

  _path = child1;
  my_path->update_path(child2);
}

////////////////////////////////////
////////////////////////////////////

// Inizializza
void System::initialize(){
  for(int i = 0; i < _npop; i++){
    for(int j = 0; j < 100*_ncities_sist; j++){
      _syst[i].mut_pairperm();
      _syst[i].mut_inversion();
      _syst[i].mut_shift();
      _syst[i].mut_permute();
    }
    _syst[i].check();
  }
}

// stampa tutti i path
void System::stamp_path(ofstream &fout){
  fout << "Path:" << endl;
  for(int i = 0; i < _npop; i++){
    for(int j = 0; j < _ncities_sist; j++){
      fout << _syst[i].get_pos(j) << "  " ;
    }
    fout << endl;
  }
}

// riordina con L1
void System::sortL1(const vector<city> &cities){
  sort(_syst.begin(), _syst.end(), [&cities](const path &a, const path &b) {
    return a.length_L1(cities) < b.length_L1(cities);
  });
}

// riordina con L2
void System::sortL2(const vector<city> &cities){
  sort(_syst.begin(), _syst.end(), [&cities](const path &a, const path &b) {
    return a.length_L2(cities) < b.length_L2(cities);
  });
}

// applica mutazione sui figli
void System::mutation(double p1, double p2, double p3, double p4){
  if (p1 < 0.0 || p1 > 1.0 || p2 < 0.0 || p2 > 1.0 || p3 < 0.0 || p3 > 1.0 || p4 < 0.0 || p4 > 1.0) {
    throw invalid_argument("Le probabilità di mutazione devono essere compresi tra zero e uno.");    
  }

  for(int i = 0; i < _son.size(); i++){
    if(_rnd->Rannyu() < p1) _syst[_son[i]].mut_pairperm();
    if(_rnd->Rannyu() < p2) _syst[_son[i]].mut_shift();
    if(_rnd->Rannyu() < p3) _syst[_son[i]].mut_permute();
    if(_rnd->Rannyu() < p4) _syst[_son[i]].mut_inversion();
  }
}

// stampa posizioni città
void System::stamp_city(const vector<city> &cities, ofstream &fout) const {

  if (cities.size() != _ncities_sist) {
    throw invalid_argument("Le città non sono in numero giusto.");    
  }

  fout << "    x      y" << endl;
  for(int i = 0; i < _ncities_sist; i++){
    fout << " " << cities[i].getx() << "  " << cities[i].gety() << endl;
  }
}

// funzione check su tutti i path
void System::check(){
  for(int i = 0; i < _npop; i++){
    _syst[i].check();
  }
}

// Funzione di crossover 1 non utilizzata
void System::crossover(double p, double pot) {
  int c = 0;
  int target = round(p * _npop);

    
  while (c < target) {
    int pa1;
    int pa2;

    bool unique1 = false;
    while (!unique1) {
      pa1 = this->select(pot);
      unique1 = true;
      for (int k = 0; k < _son.size(); k++) {
        if (pa1 == _son[k]) {
          unique1 = false;
          break;
        }
      }
    }
    _son.push_back(pa1);

    bool unique2 = false;
    while (!unique2) {
      pa2 = this->select(pot);
      unique2 = true;
      for (int k = 0; k < _son.size(); k++) {
        if (pa2 == _son[k]) {
          unique2 = false;
          break;
        }
      }
    }
    _son.push_back(pa2);

    _syst[pa1].crossover(&_syst[pa2]);
    c += 2;
  }
}

// Funzione di crossover utilizzata (e commentata)
void System::crossover2(double p, double pot) {
  int c = 0;
  int target = round(p * _npop);  // Numero di crossover da eseguire

  while (c < target) {
    int pa1;
    int pa2;

    // Trova un percorso pa1
    bool unique1 = false;
    while (!unique1) {
      pa1 = this->select(pot);  // Seleziona un percorso basato sulla probabilità potenziata
      unique1 = true;
      for (int k = 0; k < _son.size(); k++) {
        if (pa1 == _son[k]) {
          unique1 = false;  // Verifica se è già stato selezionato come figlio (_son)
          break;
        }
      }
    }
    _parents.push_back(pa1);  // Memorizza pa1 come genitore

    // Trova un percorso pa2
    bool unique2 = false;
    while (!unique2) {
      pa2 = this->select(pot);  // Seleziona un altro percorso basato sulla probabilità potenziata
      unique2 = true;
      for (int k = 0; k < _son.size(); k++) {
        if (pa2 == _son[k]) {
          unique2 = false;  // Verifica se è già stato selezionato come figlio (_son)
          break;
        }
      }
    }
    _parents.push_back(pa2);  // Memorizza pa2 come genitore

    // Esegue l'operazione di crossover
    path child1 = _syst[pa1];  // Crea un figlio inizializzato con il percorso pa1
    path child2 = _syst[pa2];  // Crea un figlio inizializzato con il percorso pa2
    child1.crossover(&child2);  // Esegue il crossover tra child1 e child2

    // Seleziona due indici casuali a1 e a2 basati sulla probabilità inversa di pot
    int a1 = this->select(1. / pot);
    int a2 = this->select(1. / pot);

    _son.push_back(a1);  // Memorizza a1 come figlio (_son)
    _son.push_back(a2);  // Memorizza a2 come figlio (_son)

    _syst[a1] = child1;  // Assegna child1 al percorso a1
    _syst[a2] = child1;  // Assegna child2 al percorso a2

    c += 2;  // Incrementa il contatore dei crossover eseguiti
  }
}

// restart dopo una generazione
void System::restart(){
  vector<int> v;
  _son = v;
}

// stampa il miglior elemento
void System::printfirstpath(ofstream &fout, const vector<city> &cities) const {

  for(int i = 0; i < _ncities_sist; i++){
    fout << " " << _syst[0].get_pos(i);
  }

  fout << "       " << _syst[0].length_L2(cities) << endl;

}

// stampa il valore medio delle lunghezze dei primi N path
void System::printmean(ofstream &fout, const vector<city> &cities, int N_mean, int d, int generation){

  if (d != 0 & d != 1) {
    throw invalid_argument("L'indice per la stampa della distanza deve essere: 0-L1, 1-L2.");    
  }

  double sum = 0;
  if(d == 0){
    for(int i = 0; i < N_mean; i++){
      sum += _syst[i].length_L1(cities);
    }
  }
    if(d == 1){
    for(int i = 0; i < N_mean; i++){
      sum += _syst[i].length_L2(cities);
    }
  }

  fout << generation << "   " << sum/N_mean << endl;
}

// Metodo che gestisce la migrazione dei percorsi tra i processi MPI
void System::migration(int rank, int size, double p) {

  bool change = false;  // Flag per indicare se la migrazione avviene o meno
  int ex1 = -1, ex2 = -1;  // Indici dei processi destinatari per la migrazione

  // Solo il processo con rank 0 decide se avviare la migrazione
  if (rank == 0) {
    // Decide casualmente se avviare la migrazione basandosi sulla probabilità p
    if (_rnd->Rannyu() < p) {
      change = true;  // Abilita la migrazione
      ex1 = int(size * _rnd->Rannyu());  // Seleziona casualmente il primo processo destinatario
      do {
        ex2 = int(size * _rnd->Rannyu());  // Seleziona casualmente il secondo processo destinatario
      } while (ex2 == ex1);  // Assicura che i due processi destinatari siano diversi
    }
  }

  // Trasmette ex1 e ex2 a tutti i processi MPI
  MPI_Bcast(&ex1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ex2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  // Se ex1 e ex2 sono stati selezionati correttamente
  if (ex1 != -1 && ex2 != -1) {
    MPI_Status stat1, stat2;
    MPI_Request req;

    // Allocazione di array temporanei per i messaggi da inviare e ricevere
    int* imesg1 = new int[_ncities_sist];
    int* imesg2 = new int[_ncities_sist];
    int itag1 = 1;
    int itag2 = 2;

    int pos_path1 = 0;  // Indice del primo percorso da scambiare
    if (rank == ex1) {
      // Se il processo corrente è ex1, copia il percorso da migrare in imesg1
      for (int j = 0; j < _ncities_sist; j++) {
        imesg1[j] = _syst[pos_path1].get_pos(j);  // Salva la posizione delle città nel percorso
      }
    }

    int pos_path2 = _npop - 1;  // Indice del secondo percorso da scambiare
    if (rank == ex2) {
      // Se il processo corrente è ex2, copia il percorso da migrare in imesg2
      for (int j = 0; j < _ncities_sist; j++) {
        imesg2[j] = _syst[pos_path2].get_pos(j);  // Salva la posizione delle città nel percorso
      }
    }

    // Se il processo corrente è ex1, invia il percorso imesg1 a ex2 e ricevi imesg2 da ex2
    if (rank == ex1) {
      MPI_Send(&imesg1[0], _ncities_sist, MPI_INTEGER, ex2, itag1, MPI_COMM_WORLD);
      MPI_Recv(&imesg2[0], _ncities_sist, MPI_INTEGER, ex2, itag2, MPI_COMM_WORLD, &stat2);
    }
    // Se il processo corrente è ex2, invia il percorso imesg2 a ex1 e ricevi imesg1 da ex1
    else if (rank == ex2) {
      MPI_Send(&imesg2[0], _ncities_sist, MPI_INTEGER, ex1, itag2, MPI_COMM_WORLD);
      MPI_Recv(&imesg1[0], _ncities_sist, MPI_INTEGER, ex1, itag1, MPI_COMM_WORLD, &stat1);
    }

    // Se il processo corrente è ex1, aggiorna il percorso locale con imesg2 ricevuto da ex2
    if (rank == ex1) {
      for (int j = 0; j < _ncities_sist; j++) {
        _syst[pos_path1].write_pos(j, imesg2[j]);  // Scrive nel percorso locale
      }
    }

    // Se il processo corrente è ex2, aggiorna il percorso locale con imesg1 ricevuto da ex1
    if (rank == ex2) {
      for (int j = 0; j < _ncities_sist; j++) {
        _syst[pos_path2].write_pos(j, imesg1[j]);  // Scrive nel percorso locale
      }
    }

    // Scambia i percorsi di partenza e arrivo per il secondo scambio
    pos_path1 = _npop - 1;  // Indice del primo percorso da scambiare
    if (rank == ex1) {
      // Se il processo corrente è ex1, copia il percorso da migrare in imesg1
      for (int j = 0; j < _ncities_sist; j++) {
        imesg1[j] = _syst[pos_path1].get_pos(j);  // Salva la posizione delle città nel percorso
      }
    }

    pos_path2 = 0;  // Indice del secondo percorso da scambiare
    if (rank == ex2) {
      // Se il processo corrente è ex2, copia il percorso da migrare in imesg2
      for (int j = 0; j < _ncities_sist; j++) {
        imesg2[j] = _syst[pos_path2].get_pos(j);  // Salva la posizione delle città nel percorso
      }
    }

    // Se il processo corrente è ex1, invia il percorso imesg1 a ex2 e ricevi imesg2 da ex2
    if (rank == ex1) {
      MPI_Send(&imesg1[0], _ncities_sist, MPI_INTEGER, ex2, itag1, MPI_COMM_WORLD);
      MPI_Recv(&imesg2[0], _ncities_sist, MPI_INTEGER, ex2, itag2, MPI_COMM_WORLD, &stat2);
    }
    // Se il processo corrente è ex2, invia il percorso imesg2 a ex1 e ricevi imesg1 da ex1
    else if (rank == ex2) {
      MPI_Send(&imesg2[0], _ncities_sist, MPI_INTEGER, ex1, itag2, MPI_COMM_WORLD);
      MPI_Recv(&imesg1[0], _ncities_sist, MPI_INTEGER, ex1, itag1, MPI_COMM_WORLD, &stat1);
    }

    // Se il processo corrente è ex1, aggiorna il percorso locale con imesg2 ricevuto da ex2
    if (rank == ex1) {
      for (int j = 0; j < _ncities_sist; j++) {
        _syst[pos_path1].write_pos(j, imesg2[j]);  // Scrive nel percorso locale
      }
    }

    // Se il processo corrente è ex2, aggiorna il percorso locale con imesg1 ricevuto da ex1
    if (rank == ex2) {
      for (int j = 0; j < _ncities_sist; j++) {
        _syst[pos_path2].write_pos(j, imesg1[j]);  // Scrive nel percorso locale
      }
    }

    // Deallocazione della memoria allocata per gli array temporanei
    delete[] imesg1;
    imesg1 = nullptr;

    delete[] imesg2;
    imesg2 = nullptr;
  }
}

// stampa il migliore tra tutti i percorsi
void System::stamp_best(ofstream &fout, int rank, int size, const vector<city>& cities){
  double* irecv = nullptr;
  double* isend = new double[2];
  isend[0] = (double)rank;
  isend[1] = _syst[0].length_L2(cities);

  if (rank == 0) {
    irecv = new double[2 * size];
  }

  MPI_Gather(isend, 2, MPI_DOUBLE, irecv, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int pos_min = 0;

  if (rank == 0) {
    double min = irecv[1];
    for (int i = 0; i < 2 * size; i += 2) {
      if (irecv[i + 1] < min) {
        min = irecv[i + 1];
        pos_min = (int)irecv[i];
      }
    }
    cout << "Best path length: " << min << endl;
  }

  MPI_Bcast(&pos_min, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == pos_min) {
    cout << "Best path for rank: " << rank << endl;
    fout << _ncities_sist << endl;
    for (int j = 0; j < _ncities_sist; j++) {
      fout << _syst[0].get_pos(j) << "  ";
    }
    fout << endl;
    fout << _syst[0].length_L2(cities) << endl;
  }

  delete[] isend;
  if (irecv != nullptr) {
    delete[] irecv;
  }

}

///////////////////////////////////////////////////
///////////////////////////////////////////////////

// Funzione di start
void start(int* type, int* n_city, int* n_iter, int* n_path, bool* paral){

  ifstream input("input.dat");
  string property;
  *type = -1;
  *n_city = -1;
  *n_iter = -1;
  *n_path = -1;
  *paral = 0;

  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
      string sim_type;
      input >> sim_type;
      bool check = false;
      if(sim_type == "circle"){
        *type = 0;
        check = true;
      }
      if(sim_type == "square"){
        *type = 1;
        check = true;
      }
      if(sim_type == "file"){
        *type = 2;
        check = true;
      }
      if(check == false){
        throw invalid_argument("Invalid argument SIMULATION TYPE."); 
      }
    } else if( property == "N_CITY" ){
      input >> *n_city;
    } else if( property == "N_ITERATION" ){
      input >> *n_iter;
    } else if( property == "N_PATH" ){
      input >> *n_path;
    } else if( property == "PARALLEL" ){
      input >> *paral;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  if(*type == -1){
    throw invalid_argument("Missing SIMULATION TYPE input."); 
  }
  if(*n_city < 0){
    throw invalid_argument("Wrong N_CITY input."); 
  }
  if(*n_iter < 0){
    throw invalid_argument("Wrong N_ITERATION input."); 
  }
  if(*n_path < 0){
    throw invalid_argument("Wrong N_PATH input."); 
  }
  input.close();
  return;
}

// Barra di progressione
void printProgressBar(int current, int total) {
    int barWidth = 70;
    float progress = (float)current / total;
    int pos = barWidth * progress;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}