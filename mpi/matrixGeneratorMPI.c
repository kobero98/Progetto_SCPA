#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <mpi.h>



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int M;
    int N;
    int K;
    int seed;
    if(argc!=6){
        printf("eseguire con ./<nomeEseguibile> <seed> <nomeCartela> <M> <N> <K>\n");
        return -1;
    }
    seed=atoi(argv[1]);  
    M=atoi(argv[3]);
    N=atoi(argv[4]);
    K=atoi(argv[5]);
    if(M == 0 || N==0|| K==0){
        printf("malloc per la creazione dei file fallita\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
        MPI_Finalize();
        return -1;
    }
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_File file;
    MPI_File fileShadow;
    char *filename=(char*)malloc(sizeof(char)*(strlen(argv[2])+2));
    char *filenameShadow=(char*)malloc(sizeof(char)*(strlen(argv[2])+6));
    if(filename==NULL || filenameShadow==NULL ){
        printf("malloc per la creazione dei file fallita\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
        MPI_Finalize();
        return -1;
    }
    // Specifica il nome del file
    if(rank==0){
        sprintf(filename,"%s/A\0",argv[2]);
        sprintf(filenameShadow,"%s/A.txt\0",argv[2]);
        N=K;
    }
    if(rank==1){
        sprintf(filename,"%s/B\0",argv[2]);
        sprintf(filenameShadow,"%s/B.txt\0",argv[2]);
        M=K;
    }
    if(rank==2){
        sprintf(filename,"%s/C\0",argv[2]);
        sprintf(filenameShadow,"%s/C.txt\0",argv[2]);
    }
    // Inizializzazione di dati locali con valori casuali
    srand(seed);  // Inizializza il seed del generatore di numeri casuali basato sul rank
    // Apertura del file per la scrittura
    // Processo 0 scrive le dimensioni della matrice nel file
    MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_File_open(MPI_COMM_SELF, filenameShadow, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fileShadow);

    if (rank < 3) {
        MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
        MPI_File_open(MPI_COMM_SELF, filenameShadow, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fileShadow);
        MPI_File_write(file, &M, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write(file, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
        char buffer[100];
        sprintf(buffer,"%d,%d\n\0",M,N);
        MPI_File_write(fileShadow, &buffer,strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
        for(int i=0;i<M;i++){
            int s;
            char buffer2[20];
            for(int j=0;j<N-1;j++){
                int s=rand()%100;
                sprintf(buffer2,"%d,\0",s);
                MPI_File_write(file, &s, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_write(fileShadow, &buffer2,strlen(buffer2), MPI_CHAR, MPI_STATUS_IGNORE);
            }
            s=rand()%100;
            sprintf(buffer2,"%d\n\0",s);
            MPI_File_write(file, &s, 1, MPI_INT, MPI_STATUS_IGNORE);
            MPI_File_write(fileShadow, &buffer2, strlen(buffer2), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        MPI_File_close(&file);
        MPI_File_close(&fileShadow);
    }
    // Chiusura del file
    MPI_Finalize();

    return 0;
}