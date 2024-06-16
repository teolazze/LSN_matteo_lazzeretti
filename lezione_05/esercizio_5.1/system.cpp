#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"
#include <vector>
#include <unordered_set>

using namespace std;
using namespace arma;



double Metropolis::ground_state(vector<double> p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  double psi = pow(_a0,-3./2)/sqrt(M_PI)*exp(-r/_a0);
  return psi*psi;
}

double Metropolis::excited_state(vector<double> p){
  double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  double psi = pow(_a0,-5./2)/8*sqrt(2/M_PI)*r*exp(-r/(2*_a0))*p[2]/r;
  return psi*psi;
}

void Metropolis::Metr(){
  vector<double> p = {0,0,0};
  if(_type_sampling == 0){
    for(int i = 0; i < 3; i++){
      p[i] = _p[i] + _rnd->Rannyu(-_delta,_delta);
    }
  }
  if(_type_sampling == 1){
    for(int i = 0; i < 3; i++){
      p[i] = _p[i] + _rnd->Gauss(0, _delta);
    }
  }
  //cout << p[0] << endl;

  double prob = -1;
  if(_type == 0) prob = this->ground_state(p)/this->ground_state(_p);
  if(_type == 1) prob = this->excited_state(p)/this->excited_state(_p);
  if(prob == -1) throw invalid_argument("type non giusto");

  //cout << prob << endl;

  if(prob > 1){
    _Nacc++;
    _p = p;
  }else{
    if(_rnd->Rannyu() < prob ){
      _Nacc++;
      _p = p;
    }
  }
  _Ntot++;
  _acceptance = (double) _Nacc/_Ntot;
}

void Metropolis::expected_r(ofstream &f1out, ofstream &f2out){
  double sum = 0;
  double sum_quad = 0;

  for(int k = 0; k < _Nskip; k++){
    for(int i = 0; i < _Niter; i++){
      this->Metr();
      f2out << _p[0] << "  " << _p[1] << "  " << _p[2] << "  " << _acceptance << endl;
    }
  }

  for(int k = 0; k < _Nblocks; k++){
    double sum_parz = 0;
    for(int i = 0; i < _Niter; i++){
      this->Metr();
      sum_parz += sqrt(_p[0]*_p[0]+_p[1]*_p[1]+_p[2]*_p[2]);
      f2out << _p[0] << "  " << _p[1] << "  " << _p[2] << "  " << _acceptance << endl;
    }
    //f2out << sum_parz/_Niter << endl;
    sum += sum_parz/_Niter;
    sum_quad += pow(sum_parz/_Niter,2);
    double incertezza = this->err(sum/(k+1), sum_quad/(k+1), k+1);
    f1out << sum/(k+1) << "  " << incertezza << endl;
  }
  _p = _p_start;
  _acceptance = 0;
  _Nacc = 0;
  _Ntot = 0;
}

double Metropolis::err(double media, double mediaquad, int n){
    if(n == 1) return 0;
    else{ 
        return sqrt((mediaquad - media * media) / (n-1));
    }
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void start(int* type, int* n_block, int* n_iter){

  ifstream input("input.dat");
  string property;
  *type = -1;
  *n_block = -1;
  *n_iter = -1;

  while ( !input.eof() ){
    input >> property;
    if( property == "COOLING" ){
      string sim_type;
      input >> sim_type;
      bool check = false;
      if(sim_type == "lin"){
        *type = 0;
        check = true;
      }
      if(sim_type == "qdr"){
        *type = 1;
        check = true;
      }
      if(sim_type == "geo"){
        *type = 2;
        check = true;
      }
      if(check == false){
        throw invalid_argument("Invalid argument SIMULATION TYPE."); 
      }
    } else if( property == "N_BLOCK" ){
      input >> *n_block;
    } else if( property == "N_ITERATION" ){
      input >> *n_iter;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  if(*type == -1){
    throw invalid_argument("Missing SIMULATION TYPE input."); 
  }
  if(*n_block < 0){
    throw invalid_argument("Wrong N_block input."); 
  }
  if(*n_iter < 0){
    throw invalid_argument("Wrong N_ITERATION input.");
  }
  input.close();
  return;
}