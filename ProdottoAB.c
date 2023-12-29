//prima implementazione prodotto matrice vettore
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>


float * func_prod(float **A,float * b,int n,int m){
    float * vet_result =(float*) malloc((m)*sizeof(float));
    for(int i=0;i<n;i++){
        float result=0.0;
        for(int j=0;j<m;j++){
            result=result+(A[i][j]*b[i]);
        }
        vet_result[i]=result;
    }
    return vet_result;
}
float * init_vet(int n,int rank){
    float * b =(float*) malloc((n)*sizeof(float));
    for(int i=0;i<n;i++){
        b[i]=i+rank;
    }
    return b;
}
float ** init_matrix(int n,int m,int rank){
    float ** A =(float**) malloc((n)*sizeof(float*));
    for(int i=0;i<n;i++){
        A[i]=init_vet(m,rank);
    }
    return A;
}
int ScomposizioneInFattori(int p){
    int lato_M,lato_m;
    int start = (int)sqrt(p);
     for(int i=start;i>=1;i--){
         if(p%i==0){
            lato_m=i;
            lato_M=p/i;
            printf("%d %d %d\n",p,lato_m,lato_M);
            break;
        }
    }
    return lato_m;
}
int main(int argc,char*argv[]){
    int rank,size;
    int tag=50;
    int i,j;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    //calcolo della latenza
    int local_l1,local_l2,local_n;
    int n,m;
    int rest;
    if(rank==0){
        printf("Inserire quante righe(N):\n");
        scanf("%d",&n);
        printf("Inserire quante colonne(M):\n");
        scanf("%d",&m);
        for(int i=1;i<size;i++){
            MPI_Send(&n,1,MPI_INT,i,tag,MPI_COMM_WORLD);
            MPI_Send(&m,1,MPI_INT,i,tag,MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(&n,1,MPI_INT,0,tag,MPI_COMM_WORLD,NULL);
        MPI_Recv(&m,1,MPI_INT,0,tag,MPI_COMM_WORLD,NULL);
        //printf("%d]n=%d,m=%d\n",rank,n,m);
    }
    //sezione che devono fare tutti i processi
    int latoCorto = ScomposizioneInFattori(size);
    int latoLungo = size/latoCorto;
    int div_m,div_n,num_X,num_Y;
    int resto_m,resto_n;
    //qui devo fare un analisi più attenta della condizione di vantaggio
    if(n>sqrt(m)){
        num_X=latoCorto;
        num_Y=latoLungo;
        div_m=m/latoCorto;
        div_n=n/latoLungo;
        resto_m=m%latoCorto;
        resto_n=n%latoLungo;
    }else{
        num_X=latoLungo;
        num_Y=latoCorto;
        div_m=m/latoLungo;
        div_n=n/latoCorto;
        resto_m=m%latoLungo;
        resto_n=n%latoCorto;
    }

    int rank_j=rank%num_X;
    int rank_i=rank/num_X;

    //deve essere vero questa cosa:-> rank = rank_i*div_m + rank_j
    //qui é da pensarci un attimo ma non dovrebbe essere difficile
    int local_start_x= div_m*(rank_j)+((rank_j<resto_m)?rank_j:resto_m);
    int local_start_y= div_n*(rank_i)+((rank_i<resto_n)?rank_i:resto_n); //qui é da dove partire per le y

    int len_n = div_n + ((rank_i < (n%div_n))?1:0); //quanto devo eseguire per le y
    int len_m = div_m + ((rank_j < (m%div_m))?1:0); //quanto devo eseguire per le x
    float ** A=init_matrix(len_n,len_m,rank);
    float * b=init_vet(len_n,rank);
    float * c=func_prod(A,b,len_n,len_m);
    int k=0;
    if(rank_j==0){
        for(int i=0;i<num_X-1;i++){
            float *vet_arrived=(float*) malloc(sizeof(float)*len_n);
            MPI_Recv(vet_arrived,len_n,MPI_FLOAT,MPI_ANY_SOURCE,tag+10,MPI_COMM_WORLD,NULL);
            for(int j=0;j<len_n;j++) c[j]=c[j]+vet_arrived[j];
        }
        for(int j=0;j<len_n;j++) printf("%d]%f\n",local_start_y+j,c[j]);
    }   
    else{
        //printf("S-%d] rank_i=%d rank_j=%d\n",rank,rank_i,rank_j);
        MPI_Ssend(c,len_n,MPI_FLOAT,rank-(rank_j),tag+10,MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
