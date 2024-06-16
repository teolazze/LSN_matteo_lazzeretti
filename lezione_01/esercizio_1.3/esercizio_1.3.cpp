#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../rangen/random.h" // Includi la libreria per generare numeri casuali

using namespace std;

// Funzione per calcolare l'errore della media
double err(double media, double mediaquad, int n);

int main (int argc, char *argv[]){

    // Apertura del file di output
    ofstream flusso_out;
    flusso_out.open("risultati.dat");
    
    // Definizione dei parametri
    int M = 1E8; // Numero totale di punti generati
    int N = 1000; // Numero di blocchi
    int L = M/N; // Numero di punti in ciascun blocco

    double d = 1; // Distanza tra le linee
    double l = 0.8; // Lunghezza dell'ago

    // Inizializzazione del generatore di numeri casuali
    Random rnd = generaRand(); 

    // Variabili per il calcolo di PI
    double sumquad_measurePI = 0;
    double sum_measurePI = 0;

    // Loop sui blocchi
    for(int i = 0; i < N; i++){
        int c = 0;

        // Loop sui punti all'interno di ciascun blocco
        for (int j = 0; j < L; j++){
            float z_rand = rnd.Rannyu(0,d);
            double x, z;
            do{
                x = rnd.Rannyu(-1,1);
                z = rnd.Rannyu(-1,1);
            }while(x*x+z*z > 1); // Accetto solo se il punto è nella circonferenza di raggio unitario
            double z1 = z_rand;
            double z2 = z_rand+l/(sqrt(x*x+z*z))*z;
            if((z1 < 0 || z1 > 1) || (z2 < 0 || z2 > 1)) c++;
        }

        // Calcolo di PI e aggiornamento delle somme
        double measurePI = 2*l*L/(c*d);  // Misura del blocco
        sum_measurePI += measurePI;  // Somma delle misure dei blocchi
        double media_measurePI = sum_measurePI/(i+1);  // Media delle misure dei blocchi
        sumquad_measurePI += measurePI*measurePI;  // Somma delle misure dei blocchi al quadrato
        double mediaquad_measurePI = sumquad_measurePI/(i+1);  // Media delle misure dei blocchi al quadrato
        double sigma = err(media_measurePI, mediaquad_measurePI, i+1);  // Deviazione standard delle misure dei blocchi

        // Scrittura dei risultati sul file di output
        flusso_out << i+1 << " " << media_measurePI << " " << sigma << endl;

    }

    // Chiusura del file di output e terminazione del programma
    flusso_out.close();
    return 0;
}


double err(double media, double mediaquad, int n){
    if(n == 1) return 0;
    else{ 
        return sqrt((mediaquad - media * media) / (n-1));
    }
}