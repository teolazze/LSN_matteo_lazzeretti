#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int isend[2], irecv[2];
    for(int i=0;i<3;i++) isend[i]=rank+i+1;

    

    MPI_Finalize();
    return 0;

}