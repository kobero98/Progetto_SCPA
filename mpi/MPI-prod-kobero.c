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
void calcolo_Computazionale(int* localA,int* localB,int* localC,int m,int n,int k){
    printf("inizio a fare i cicli\n");
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int z=0;z<k;z++){
                localC[i*n+j]+=localA[i*k+z]*localB[j*k+z];
                // printf("A[%d][%d] = %d\n", i, z, localA[i*k+z]);
                // printf("B[%d][%d] = %d\n", i, z, localB[j*k+z]);
                // printf("C[%d][%d] = %d\n", i, j, localC[i*n+j]);
            }
        }
    }
}
void init_localC(int * C,int n,int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            C[i*m+j]=1;
        }
    }
}
void init_Vector(int *vet,int dim){
    for(int i=0;i<dim;i++){
        vet[i]=1;
    }
}
/*
    questa é da modificare
*/
void func_gen_data(int my_rank,int  my_coordinate[2],int sizeA,int sizeB,int sizeC,float * A,float * B,float * C){
        int i=0;
        for(;i<sizeA;i++){
            A[i]=1.0;
        }
        i=0;
        for(;i<sizeB;i++){
            B[i]=1.0;
        }
        i=0;
        for(;i<sizeC;i++){
            C[i]=1.0;
        }
}

int main(int argc, char **argv) {
    //informazioni che il processo 0 deve inviare agli altri processi
    int K=80;  //Numero Colonne A //Numero righe B
    int M=80;  //Numero Righe di A //Numero righe di C 
    int N=80;  //Numero Colonne di B //Numero colonne di C
    int nb =1;      //Iper-parametro della mattonella
    int mb =1;      //iper-parametro della mttonella
    //porzioni locali delle matrici
    int *localA;
    int *localB;
    int *localC;
    //variabili MPI
    int my_rank;
    int my_coord[DIMS];
    int p;
    MPI_Comm comm_world_copy;   //copia del comunicatore MPI_COMM_WORLD (è sempre buona norma averla)
    MPI_Comm comm_cart;         //nuovo comunicatore relativo alla griglia di processi da associare alla matrice C
    //variabili utili per suddividere il lavoro dei processi
    int n;                      //numero di righe per la sottomatrice (di C) da assegnare a ciascun processo
    int m;                      //numero di colonne per la sottomatrice (di C) da assegnare a ciascun processo
    int proc_dims[DIMS];        //array che indica il numero di processi che va a finire in ciascuna dimensione della matrice C (--> la mesh coinvolgerà la matrice C)
    int proc_dims_hat[DIMS];    //array che indica il numero di righe e colonne che, nella mesh di processi, sono composte da una riga/colonmna di C in più (le entry sono risp. N%n, M%m)
    //variabili di appoggio
    int periods[DIMS];          //array che indica se ciascuna dimensione della matrice deve essere periodica (i.e. circolare) o meno
    int mesh_row_index;         //variabile che tiene traccia dell'indice riga all'interno della mesh di processi
    int mesh_col_index;         //variabile che tiene traccia dell'indice colonna all'interno della mesh di processi
    int my_mesh_row;            //il PROPRIO indice riga all'interno della mesh di processi
    int my_mesh_col;            //il PROPRIO indice colonna all'interno della mesh di processi
    //indici ciclo for
    int i;
    int j;
    int k;

    MPI_Init(NULL, NULL);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_copy); //duplicazione del comunicatore MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //printf("Hello from process %d of %d!\n", my_rank, p);
    
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

    //TODO: Check malloc
    int* my_row_rank=(int*)malloc(sizeof(int)*proc_dims[1]);
    int* my_column_rank=(int*)malloc(sizeof(int)*proc_dims[0]);
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
    localA =(int*) malloc((sizeof(int)*K*m));
    localB =(int*) malloc((sizeof(int)*K*n));
    localC =(int*) malloc((sizeof(int)*n*m));
    if(!(localA && localB && localC)) {
        printf("Unable to allocate local_A, local_B, local_C and intermediate_C.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }
    //Condivido della matrice A
    for(i=0;i<m;i++){
        if(my_coord[1]==0){
           init_Vector(&(localA[i*K]),K);
        }    
        MPI_Bcast(&(localA[i*K]), K, MPI_INT,0, my_comm_row);
    }
    //La matrice B viene creata per colonne
    //Condivisione della matrice B
    for(i=0;i<n;i++){
        if(my_coord[0]==0){
           init_Vector(&(localB[i*K]),K);
        }    
        MPI_Bcast(&(localB[i*K]), K, MPI_INT,0, my_comm_column);
    }
    //La matrice C ognuno si alloca la sua
    init_localC(localC,n,m);
    //Calcolo Cooutazionale
    calcolo_Computazionale(localA,localB,localC,m,n,K);
    //printf("%d finito la computazion\n",my_rank);
    // if(col_comms[0] != MPI_COMM_NULL)    //check su se il processo appartiene effettivamente al comunicatore col_comms[0]
    //     MPI_Gatherv(intermediate_C, block_sizes_C[my_mesh_row], MPI_INT, total_C, block_sizes_C, block_displs_C, MPI_INT, 0, col_comms[0]);
    //stampa del risultato finale
    printf("RESULT:\n");
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            printf("%d) C[%d][%d] = %d\n",my_rank,i*proc_dims[0]+my_coord[0], j*proc_dims[1]+my_coord[1], localC[i*m+j]);
        }
    }
    MPI_Finalize();
    return 0;
}
