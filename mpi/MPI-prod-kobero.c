#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


#define DIMS 2  //number of dimensions of a matrix
#define NO (0)
#define YES NO+1
#define DEBUG

//va valutato il giusto ordine di operazioni
//e si possono fare meglio di n^3 operazioni?
void calcolo_Computazionale(float* localA,float* localB,float* localC,int m,int n,int k){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0;z<k;z++){
                localC[i*n+j]+=localA[i*k+z]*localB[j*k+z];
            }
        }
    }
}
void init_localC(float * C,int n,int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            C[i*m+j]=1.;
        }
    }
}
void init_Vector(float *vet,int dim){
    for(int i=0;i<dim;i++){
        vet[i]=1.;
    }
}
/*
    questa é da modificare
*/
void init_SingleCase(float* A,float* B,float* C,int M,int N,int K){
    //init A matrice M*K
    for(int i=0;i<M;i++){
        for (int j = 0; j < K; j++) A[i*K+j]=1.;
    }
    //init B matrice K*N
    for(int i=0;i<K;i++){
        for (int j = 0; j < N; j++) B[i*N+j]=1.;
    }
    //init C matrice M*N
    for(int i=0;i<M;i++){
        for (int j = 0; j < N; j++) C[i*N+j]=1.;
    }
}
int main(int argc, char **argv) {

    //informazioni che il processo 0 deve inviare agli altri processi
    int K; //Numero Colonne A //Numero righe B
    int M;  //Numero Righe di A //Numero righe di C 
    int N;  //Numero Colonne di B //Numero colonne di C
    int nb =1;      //Iper-parametro della mattonella
    int mb =1;      //iper-parametro della mttonella
    //porzioni locali delle matrici
    float *localA;
    float *localB;
    float *localC;
    float *C;   //variabile utilizzata solo dal processo 0
    //variabili MPI
    int my_rank;        //il rank personale all'interno dell comunicatore globale
    int my_coord[DIMS]; //le cordinate sulla griglia all'interno del comunicatore cartesiano 
    int p;              //Numero di processi utilizzati
    MPI_Comm comm_world_copy;   //copia del comunicatore MPI_COMM_WORLD (è sempre buona norma averla)
    MPI_Comm comm_cart;         //nuovo comunicatore relativo alla griglia di processi da associare alla matrice C
    //variabili utili per suddividere il lavoro dei processi
    int n;                      //numero di colonne per la sottomatrice (di C) da assegnare a ciascun processo
    int m;                      //numero di righe per la sottomatrice (di C) da assegnare a ciascun processo
    int proc_dims[DIMS];        //[0]: indica le righe della griglia dei processi [1]: indica le colonne
    //variabili di appoggio
    int periods[DIMS];          //array che indica se ciascuna dimensione della matrice deve essere periodica (i.e. circolare) o meno
    //indici ciclo for
    int i;
    int j;
    int k;
    //variabili per misurare le prestazioni
    double start;
    double startAftearCreate;
    double end;
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_copy); //duplicazione del comunicatore MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(argc != 4){
        printf("%d input non valido \n",my_rank);
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }
    K=atoi(argv[1]);  
    M=atoi(argv[2]);   
    N=atoi(argv[3]);  
    if(K==0 || M==0 || N==0){
        printf("%d input non valido \n",my_rank);
        MPI_Abort(comm_world_copy, EXIT_FAILURE);    
    }
    float* C_SingoleCase;
    //Elaborazione matrice singola elaborazione fatta solo dal processo 0
    //serve per misurare la correttezza del programma 
    double startSingleExecution;
    double endSingleExecution;
    double startSingleComputation;
    if(my_rank == 0){
        startSingleExecution=MPI_Wtime();
        float* A_SingoleCase=(float*) malloc(sizeof(float)*M*K);
        float* B_SingoleCase=(float*) malloc(sizeof(float)*K*N);
        C_SingoleCase=(float*) malloc(sizeof(float)*M*N);
        init_SingleCase(A_SingoleCase,B_SingoleCase,C_SingoleCase,M,N,K);
        startSingleComputation=MPI_Wtime();
        calcolo_Computazionale(A_SingoleCase,B_SingoleCase,C_SingoleCase,M,N,K);
        endSingleExecution=MPI_Wtime();
        free(A_SingoleCase);
        free(B_SingoleCase);
    }

    MPI_Barrier(comm_world_copy);
    start=MPI_Wtime();
    //inizializzazione di periods a soli false + inizializzazione di proc_dims a soli 0
    for(i=0; i<DIMS; i++) {
        periods[i] = NO;
        proc_dims[i] = 0;   //proc_dims[0] = #righe della mesh di processi; proc_dims[1] = #colonne della mesh di processi
    }
    //definizione (e successiva creazione) della topologia più quadrata possibile per i p processi
    MPI_Dims_create(p, DIMS, proc_dims);
    MPI_Cart_create(comm_world_copy, DIMS, proc_dims, periods, NO, &comm_cart);
    MPI_Cart_coords(comm_cart, my_rank, DIMS,my_coord);
    
    MPI_Group totalGroup;
    //TODO: controllare gli errori
    MPI_Comm_group(comm_world_copy, &totalGroup);
    
    //calcolo delle dimensioni di base delle sottomatrici (di C) da assegnare a ciascun processo
    n = N/(proc_dims[1]*nb);
    int resto_colonne = N-n*(proc_dims[1]*nb);
    m = M/(proc_dims[0]*mb);
    int resto_righe = N-n*(proc_dims[0]*mb);
    if(my_coord[1]<resto_colonne) n++;
    if(my_coord[0]<resto_righe) m++;

    int* my_row_rank=(int*)malloc(sizeof(int)*proc_dims[1]);
    int* my_column_rank=(int*)malloc(sizeof(int)*proc_dims[0]);
    if(my_row_rank == NULL || my_column_rank==NULL)
    {
        printf("%d Unable to allocate array of int for column or row \n",my_rank);
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }

    int coords[DIMS];
    coords[0]=my_coord[0];
    //mi metto i processi inerenti alla mia riga
    for(i=0;i<proc_dims[1];i++){
        coords[1]=i;
        MPI_Cart_rank(comm_cart ,coords ,&(my_row_rank[i]));
    }
    MPI_Group my_group_row;
    MPI_Comm my_comm_row;
    //TODO: check errori
    MPI_Group_incl(totalGroup, proc_dims[1], my_row_rank,&my_group_row);
    MPI_Comm_create_group(comm_world_copy,my_group_row,0,&my_comm_row);//Lo 0 é il tag
    if(my_comm_row == MPI_COMM_NULL)
    {
        printf("Unable to allocate comunicatore per righe.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }
    free(my_row_rank);

    coords[1]=my_coord[1];
    //mi segno i processi inerenti alla mia colonna
    for(i=0;i<proc_dims[0];i++){
        coords[0]=i;
        MPI_Cart_rank(comm_cart ,coords,&(my_column_rank[i]));
    }
    //TODO: check errori
    MPI_Group my_group_column;
    MPI_Comm my_comm_column;
    MPI_Group_incl(totalGroup, proc_dims[0], my_column_rank,&my_group_column);
    MPI_Comm_create_group(comm_world_copy, my_group_column,0,&my_comm_column);//Lo 0 é il tag
    if(my_comm_column==MPI_COMM_NULL)
    {
        printf("Unable to allocate comunicatore per colonne.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }
    free(my_column_rank);

    localA =(float*) malloc((sizeof(float)*K*m));
    localB =(float*) malloc((sizeof(float)*K*n));
    localC =(float*) malloc((sizeof(float)*n*m));
    if(my_rank==0) {
        C=(float*) malloc((sizeof(int)*N*M));
        if(C==NULL){
            printf("processo 0 Unable to allocate total Matrix\n");
            MPI_Abort(comm_world_copy, EXIT_FAILURE);
        }
    }
    if(!(localA && localB && localC)) {
        printf("Unable to allocate local_A, local_B, local_C and intermediate_C.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }
    //Condivido della matrice A
    for(i=0;i<m;i++){
        if(my_coord[1]==0){
           init_Vector(&(localA[i*K]),K);
        }    
        MPI_Bcast(&(localA[i*K]), K, MPI_FLOAT,0, my_comm_row);
    }
    //La matrice B viene creata per colonne
    //Condivisione della matrice B
    for(i=0;i<n;i++){
        if(my_coord[0]==0){
           init_Vector(&(localB[i*K]),K);
        }    
        MPI_Bcast(&(localB[i*K]), K, MPI_FLOAT,0, my_comm_column);
    }
    //La matrice C ognuno si alloca la sua
    init_localC(localC,n,m);

    MPI_Barrier(comm_world_copy);
    startAftearCreate=MPI_Wtime();
    //Calcolo Cooutazionale
    calcolo_Computazionale(localA,localB,localC,m,n,K);

    //stampa del risultato finale
    if(my_rank==0){
        //tratto i dati del processo 0 in modo differente
        for(j=0;j<m;j++){
                for(k=0;k<n;k++){
                    int row=j*proc_dims[0];
                    int column=k*proc_dims[1];
                    C[row*m+column]=localC[j*m+k];
                }
            }
        //copio tutte le matrici locali all'interno del processo 0
        for(i=1;i<p;i++){
            int size[2];
            MPI_Status status;
            MPI_Recv(&size,2, MPI_FLOAT,i, 0, comm_world_copy, &status);
            //chech status if(status.MPI_ERROR<0)
            float * buffer=(float*) malloc(sizeof(int)*size[0]*size[1]);
            MPI_Recv(buffer, size[0]*size[1], MPI_FLOAT,i,0, comm_world_copy,&status);
            int coord_i[DIMS];
            MPI_Cart_coords(comm_cart, i, DIMS,coord_i);
            for(j=0;j<size[0];j++){
                for(k=0;k<size[1];k++){
                    int row=j*proc_dims[0]+coord_i[0];
                    int column=k*proc_dims[1]+coord_i[1];
                    C[row*M+column]=buffer[j*size[1]+k];
                }
            }
            free(buffer);
        }
    }
    else{
        int dim[2]={m,n};
        MPI_Ssend(dim,2,MPI_FLOAT,0,0,comm_world_copy);
        MPI_Ssend(localC,m*n,MPI_FLOAT,0,0,comm_world_copy);
    }
    MPI_Barrier(comm_world_copy);
    end=MPI_Wtime();
    //calcolo MAX DIFF
    if(my_rank==0){
        double maxErr=0.0;
        int countt=0;
        for(i=0; i<M; i++) {
            for(j=0; j<N; j++) {
                if(maxErr<(C[i*N+j] - C_SingoleCase[i*N+j])){
                    maxErr= C[i*N+j] - C_SingoleCase[i*N+j];
                    countt++;
                }
            }
        }
        printf("Max Error=%f Number error=%d\n",maxErr,countt);
        free(C);
        free(C_SingoleCase);
        printf("tempo esecuzione caso singolo tempo totale=%fs tempo senza creazione=%fs",endSingleExecution-startSingleExecution,endSingleExecution-startSingleComputation);
        printf("speed up: %f/%f=%f ",endSingleExecution-startSingleComputation,end-startAftearCreate,(endSingleExecution-startSingleComputation)/(end-startAftearCreate))
    }
    free(localA);
    free(localB);
    free(localC);
    printf("%d)tempo senza creazione=%f  tempo Totale=%fs\n",my_rank,end-startAftearCreate,end-start);
    MPI_Finalize();
    return 0;
}
