//Seconda implementazione prodotto matrice vettore
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#define isNull(X) \
        if(X==NULL) return NULL;
float * init_vet(int n,int rank){
    float * b =(float*) malloc((n)*sizeof(float));
    isNull(b);
    if(b==NULL){
        return NULL;

    }
    for(int i=0;i<n;i++){
        b[i]=i+rank;
    }
    return b;
}
float * func_prod(float **A,float * b,int n,int m){
    float * vet_result =(float*) malloc((m)*sizeof(float));
    isNull(vet_result);
    for(int i=0;i<n;i++){
        float result=0.0;
        for(int j=0;j<m;j++){
            result=result+(A[i][j]*b[i]);
        }
        vet_result[i]=result;
    }
    return vet_result;
}
//questa secondo me va più veloce
float * func_prod2(float **A,float * b,int n,int m){
    float * vet_result =(float*) calloc(sizeof(float),n);
    isNull(vet_result);
    for(int j=0;j<m;j++){
        for(int i=0;i<n;i++){
            vet_result[j]=vet_result[j]+(A[i][j]*b[i]);
        }
    }
    return vet_result;
}
void free_matrix(float ** A,int n){
    int i=0;
    for(;i<n;i++){
        free(A[i]);
    }
    free(A);
}
float ** init_matrix(int n,int m,int rank){
    float ** A =(float**) malloc((n)*sizeof(float*));
    isNull(A);
    for(int i=0;i<n;i++){
        A[i]=init_vet(m,rank);
        isNull(A[i]);
    }
    return A;
}

int main(int argc,char* argv[]){
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    //Aquisisco i dati
    int n,m;
    if(rank==0){
        printf("Inserire quante righe(N):\n");
        scanf("%d",&n);
        printf("Inserire quante colonne(M):\n");
        scanf("%d",&m);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);
    //tutti scoprono tutto
    int ndims[2]={0,0};
    MPI_Dims_create(size,2,ndims);
    if(rank==0) printf("%d %d\n",ndims[0],ndims[1]);
    int periods[2] = {0,0};
    MPI_Comm newCommunicator;
    MPI_Cart_create(MPI_COMM_WORLD,2,ndims,periods,0,&newCommunicator);
    int coords[2];
    MPI_Cart_coords(newCommunicator,rank,2,coords);

    //calcolo le due lunghezze
    int len_n,len_m;
    if(coords[1] < n%ndims[1]) len_n=n/ndims[1]+1;
    else len_n=n/ndims[1];

    if(coords[0] < m%ndims[0]) len_m=n/ndims[0]+1;
    else len_m=m/ndims[0];
    printf("%d] (%d,%d)\n",rank,coords[0],coords[1]);
    float ** A=init_matrix(len_n,len_m,rank);
    if(A==NULL){
        printf("Matrice A is NULL\n");
        return 0;
    }
    float * b=init_vet(len_n,rank);
    if(b==NULL){
        printf("vettore b is NULL\n");
        return 0;
    }
    float * c=func_prod2(A,b,len_n,len_m);
    if(c==NULL){
        printf("vettore c is NULL\n");
        return 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free_matrix(A,len_n);
    free(b);
    int tag=50;
    if(coords[1]==0){
        for(int i=1;i<ndims[1];i++){
                float *vet_arrived=(float*) malloc(sizeof(float)*len_n);
                MPI_Recv(vet_arrived,len_n,MPI_FLOAT,MPI_ANY_SOURCE,tag,newCommunicator,NULL);
                for(int j=0;j<len_n;j++) c[j]=c[j]+vet_arrived[j];
                free(vet_arrived);
        }
        for(int j=0;j<len_n;j++) printf("%d]%f\n",j,c[j]);
    }else{
        int dest;
        int destCoords[2]={coords[0],0};
        MPI_Cart_rank(newCommunicator, destCoords, &dest);
        MPI_Send(c,len_n,MPI_FLOAT,dest,tag,newCommunicator);
    }
    free(c);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
