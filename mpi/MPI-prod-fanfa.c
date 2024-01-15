#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define DIMS 2  //number of dimensions of a matrix
#define RAND_DIV 1000000000.0   //divisor of rand() when matrixes components are generated

#define NO 0
#define YES NO+1

#define DEBUG



/* function which has the responsibility to allocate the array representing the association rank - mesh coordinates for all the processes */
int **create_all_cart_coords(int p, MPI_Comm comm_world_copy, MPI_Comm comm_cart) {
    //local auxiliary variables
    int i;  //loop index
    //function output
    int **all_cart_coords;

    //allocation of all_cart_coords
    all_cart_coords = (int **)malloc(p*sizeof(int *));
    if(!all_cart_coords)
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    for(i=0; i<p; i++) {
        //allocation of all_cart_coords[i]
        all_cart_coords[i] = (int *)malloc(DIMS*sizeof(int));
        if(!all_cart_coords[i])
            MPI_Abort(comm_world_copy, EXIT_FAILURE);

        //get the coordinates (mesh indexes) for process i
        MPI_Cart_coords(comm_cart, i, DIMS, all_cart_coords[i]);

    }

    return all_cart_coords;

}





/* function which has the responsibility to allocate the array pointing to all the arrays of process ranks of each mesh row */
int **create_row_ranks_list(int num_mesh_rows, int num_mesh_cols, int p, int **all_cart_coords, MPI_Comm comm_world_copy) {
    //local auxiliary variables (loop indexes)
    int i;
    int j;
    int l;
    //function output
    int **row_ranks_list;

    //allocation of row_ranks_list[i]
    row_ranks_list = (int **)malloc(num_mesh_rows*sizeof(int *));
    if(!row_ranks_list)
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    for(i=0; i<num_mesh_rows; i++) {
        //allocation of row_ranks_list
        row_ranks_list[i] = (int *)malloc(num_mesh_cols*sizeof(int));
        if(!row_ranks_list[i])
            MPI_Abort(comm_world_copy, EXIT_FAILURE);

        /* Here we are populating the ranks list associated to i-th row of the mesh.
         * Index j represents the j-th entry of this rank list.
         * Index l is only useful to iterate on all_cart_coords array to check which processes are exactly in the i-th row of the mesh.
         */
        for(j=0; j<num_mesh_cols; j++) {
            for(l=0; l<p; l++) {
                if(all_cart_coords[l][0] == i) {
                    row_ranks_list[i][j] = l;   //l = rank of the found process having row index i in the mesh
                    break;  //once we have populated the j-th entry, we can stop to scan all_cart_coords.
                }

            }

        }

    }

    return row_ranks_list;

}





/* function which has the responsibility to allocate the array pointing to all the arrays of process ranks of each mesh column */
int **create_col_ranks_list(int num_mesh_rows, int num_mesh_cols, int p, int **all_cart_coords, MPI_Comm comm_world_copy) {
    //local auxiliary variables (loop indexes)
    int i;
    int j;
    int l;
    //function output
    int **col_ranks_list;

    //allocation of col_ranks_list
    col_ranks_list = (int **)malloc(num_mesh_cols*sizeof(int *));
    if(!col_ranks_list)
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    for(i=0; i<num_mesh_cols; i++) {
        //allocation of col_ranks_list[i]
        col_ranks_list[i] = (int *)malloc(num_mesh_cols*sizeof(int));
        if(!col_ranks_list[i])
            MPI_Abort(comm_world_copy, EXIT_FAILURE);

        /* Here we are populating the ranks list associated to i-th column of the mesh.
         * Index j represents the j-th entry of this rank list.
         * Index l is only useful to iterate on all_cart_coords array to check which processes are exactly in the i-th column of the mesh.
         */
        for(j=0; j<num_mesh_rows; j++) {
            for(l=0; l<p; l++) {
                if(all_cart_coords[l][1] == i) {
                    col_ranks_list[i][j] = l;   //l = rank of the found process having column index i in the mesh
                    break;  //once we have populated the j-th entry, we can stop to scan all_cart_coords.
                }

            }

        }

    }

    return col_ranks_list;

}





/* function which has the responsibility to calculate the number of rows (or, alternatively, the number of columns) of matrix A 
 * ruled by the current process.
 *
 * IDEA:
 * rows_superblock = mb*proc_dims[0]
 * cols_superblock = nb*proc_dims[1]
 * 
 * superblocks_in_mesh_col = M/rows_superblock
 * superblocks in mesh_row = K/cols_superblock
 *
 * exceeding_rows = M%rows_superblock
 * exceeding_cols = K%cols_superblock
 *
 * k = kb*superblocks_in_mesh_row
 * m = mb*superblocks_in_mesh_col
 */
int get_local_dim_A(int matrix_dim, int block_dim, int *proc_dims, int **all_cart_coords, int my_rank, int dim) {
    /* LOCAL AUXILIARY VARIABLES */
    //number of rows and columns in each (complete) superblock
    int rows_or_cols_superblock;
    //number of superblocks in each mesh row and in each mesh column
    int superblocks_in_mesh_col_or_row;
    //number of exceeding rows and exceeding columns with respect to the suddivision of the matrix given by the mesh
    int exceeding_rows_or_cols;
    int exceeding_cols;
    //minimum number of rows (m) and minimum number of columns (k) of local_A for all the processes
    int m_or_k;

    //get the number of rows or columns in each (complete) superblock
    rows_or_cols_superblock = block_dim*proc_dims[dim];
    //get the number of superblocks in each mesh row or in each mesh column
    superblocks_in_mesh_col_or_row = matrix_dim/rows_or_cols_superblock;
    //get the number of exceeding rows or exceeding columns with respect to the suddivision of the matrix given by the mesh
    exceeding_rows_or_cols = matrix_dim % rows_or_cols_superblock;
    //get the minimum number of rows (m) or the minimum number of columns (k) of local_A foreach process
    m_or_k = block_dim * superblocks_in_mesh_col_or_row;

    //return the number of rows or columns of local_A for the current process
    if(all_cart_coords[my_rank][dim] < exceeding_rows_or_cols / block_dim) //exceeding_rows/mb = #processes which will have mb extra rows
        return m_or_k + block_dim;                                         //exceeding_cols/kb = #processes which will have kb extra columns
    else if(all_cart_coords[my_rank][dim] == exceeding_rows_or_cols / block_dim)    //exceeding_rows%mb = #extra rows that process with rank == exceeding_rows/mb has
        return m_or_k + (exceeding_rows_or_cols % block_dim);                       //exceeding_cols%kb = #extra columns that process with rank == exceeding_cols/kb has
    else
        return m_or_k;

}





int main(int argc, char **argv) {
    /* INFORMATION THAT IS PROVIDED BY THE USER AND PROCESS 0 HAS TO SEND TO OTHER PROCESSES */
    //total dimensions of the matrixes
    int K;  //#rows of matrix B and #columns of matrix A (intially try with 12)
    int N;  //#columns of matrixes B, C (initially try with 8)
    int M;  //#rows of matrixes A, C (initially try with 10)
    //dimensions of matrix blocks
    int kb; //#columns of one block in matrix A (initially try with 2)
    int mb; //#rows of one block in matrix A (initially try with 2)

    /* VARIABLES WHERE THE MATRIXES ARE DEFINED */
    //total matrixes
    float *A;
    float *B;
    float *C;
    //local portions of the matrixes
    float *local_A;
    float *local_B;
    float *local_C;

    /* MPI VARIABLES */
    int my_rank;
    int p;                          //#processes
    MPI_Comm comm_world_copy;       //copy of MPI_COMM_WORLD
    MPI_Comm comm_cart;    //new communicator related to processes topology, which is related to a super-block in matrix A
    //communicators related to the single mesh rows and to the single mesh columns
    MPI_Comm *row_comms;    //communicators related to the mesh rows
    MPI_Comm *col_comms;    //communicators related to the mesh columns
    //processes groups
    MPI_Group total_group;  //group associated to comm_world_copy communicator
    MPI_Group *row_groups;  //groups associated to the single communicators related to the mesh rows
    MPI_Group *col_groups;  //groups associated to the single communicators related to the mesh columns 

    /* VARIABLES USEFUL TO DIVIDE THE WORK BETWEEN THE PROCESSES */
    //array representing the number of processes on each dimension of super-block in matrix A (super-block = set of p blocks where each block is assigned to a different process)
    int proc_dims[DIMS];
    //number of rows and columns of local_A, local_B and local_C for the current process
    int my_rows_A;
    int my_cols_A;
    int my_rows_B;
    int my_cols_B;
    int my_rows_C;
    int my_cols_C;

    /* AUXILIARY VARIABLES */
    int periods[DIMS];      //array indicating if each dimension of processes topology has to be periodic (i.e. circular) or not
    //data structures
    int **all_cart_coords;  //array associating the rank with the corresponding mesh indexes foreach process (where rank = index of all_cart_coords)
    int **row_ranks_list;   //array representing all the process ranks foreach row of the mesh (where row of the mesh = index of row_ranks_list)
    int **col_ranks_list;   //array representing all the process ranks foreach column of the mesh (where column of the mesh = index of col_ranks_list)

    /* LOOP INDEXES */
    int i;
    int j;
    int k;



    /* MPI INITIALIZATION */
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_copy); //duplication of MPI_COMM_WORLD communicator
    MPI_Comm_size(comm_world_copy, &p);
    MPI_Comm_rank(comm_world_copy, &my_rank);

    /* CHECK OF #ARGUMENTS PASSED BY THE USER */
    if(argc < 6) {
        if(my_rank == 0) {
            printf("Usage: MPI-prod K N M kb mb\n");
            fflush(stdout);
        }
        return 0;
    }



    /* CREATION OF ALL COMMUNICATORS */
    //initialization of periods with only false + initialization of proc_dims with only 0
    for(i=0; i<DIMS; i++) {
        periods[i] = NO;
        proc_dims[i] = 0;   //proc_dims[0] = #rows of processes mesh; proc_dims[1] = #columns of processes mesh
    }

    //definition (and following creation) of the most squared possible topology (mesh) for the p processes
    MPI_Dims_create(p, DIMS, proc_dims);    //after here, proc_dims = [#processes as rows of the mesh; #processes as columns of the mesh]
    MPI_Cart_create(comm_world_copy, DIMS, proc_dims, periods, NO, &comm_cart);
    //definition of the processes group which partecipate to comm_world_copy communicator (i.e. group of all the processes)
    MPI_Comm_group(comm_world_copy, &total_group);

    //initialization of the array associating the rank with the corresponding mesh indexes foreach process
    all_cart_coords = create_all_cart_coords(p, comm_world_copy, comm_cart);

    //initialization of the array of the processes rows (where each row is represented by an array of ranks)
    row_ranks_list = create_row_ranks_list(proc_dims[0], proc_dims[1], p, all_cart_coords, comm_world_copy);
    //initialization of the array of the processes columns (where each column is represented by an array of ranks)
    col_ranks_list = create_col_ranks_list(proc_dims[0], proc_dims[1], p, all_cart_coords, comm_world_copy);

    //memory allocation for row groups, col_groups, row_comms and col_comms
    row_groups = (MPI_Group *)malloc(proc_dims[0]*sizeof(MPI_Group));   //one group foreach mesh row
    col_groups = (MPI_Group *)malloc(proc_dims[1]*sizeof(MPI_Group));   //one group foreach mesh column
    row_comms = (MPI_Comm *)malloc(proc_dims[0]*sizeof(MPI_Comm));      //one communicator foreach mesh row
    col_comms = (MPI_Comm *)malloc(proc_dims[1]*sizeof(MPI_Comm));      //one communicator foreach mesh column
    if(!(row_comms && col_comms))
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    //get all the sub-groups related to the single mesh rows 
    for(i=0; i<proc_dims[0]; i++) { //proc_dims[0] = #mesh rows
        //remember that we have a number of sub-groups equal to the number of mesh rows and a number of processes in each sub-group equal to the number of mesh columns.
        MPI_Group_incl(total_group, proc_dims[1], row_ranks_list[i], &row_groups[i]);
        MPI_Comm_create(comm_world_copy, row_groups[i], &row_comms[i]);   //communicators row_comms[i] will be used in MPI_Scatterv().
    }
    //get all the sub-groups related to the singole mesh columns
    for(i=0; i<proc_dims[1]; i++) { //proc_dims[1] = #mesh columns
        //remember that we have a number of sub-groups equal to the number of mesh columns and a number of processes in each sub-group equal to the number of mesh rows.
        MPI_Group_incl(total_group, proc_dims[0], col_ranks_list[i], &col_groups[i]);
        MPI_Comm_create(comm_world_copy, col_groups[i], &col_comms[i]);   //communicators col_comms[i] will be used in MPI_Scatterv().
    }



    /* INFORMATION EXCHANGE */
    //we assume that only process 0 initially knows the values of K, N, M, kb, mb.
    if(my_rank == 0) {
        //user input acquisition
        K = atoi(argv[1]);
        N = atoi(argv[2]);
        M = atoi(argv[3]);
        kb = atoi(argv[4]);
        mb = atoi(argv[5]);
    }

    //process 0 sends to all the other processes the acquired information.
    MPI_Bcast(&K, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&N, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&M, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&kb, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&mb, 1, MPI_INT, 0, comm_world_copy);

    /* DIVISION, ALLOCATION AND INITIALIZATION OF THE MATRIXES */
    //get the number of rows and columns of local_A for the current process
    my_rows_A = get_local_dim_A(M, mb, proc_dims, all_cart_coords, my_rank, 0);
    my_cols_A = get_local_dim_A(K, kb, proc_dims, all_cart_coords, my_rank, 1);
    //get also the number of rows and columns of local_B and local_C for the current process
    my_rows_B = my_cols_A;
    my_cols_B = N;
    my_rows_C = my_rows_A;
    my_cols_C = N;

    //allocate the local matrices
    local_A = (float *)malloc(my_rows_A*my_cols_A*sizeof(float));
    local_B = (float *)malloc(my_rows_B*my_cols_B*sizeof(float));
    local_C = (float *)malloc(my_rows_C*my_cols_C*sizeof(float));
    if(!(local_A && local_B && local_C))
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    //local_A has to be initialized by all the processes.
    for(i=0; i<my_rows_A*my_cols_A; i++) {
        local_A[i] = rand() / RAND_DIV;
    }

    //local_B has to be initialized by the processes at first row of the mesh. Its content has to be sent in broadcast on col_comms[i] (foreach i).
    if(all_cart_coords[my_rank][0] == 0) {  //i.e. if I am a process at first row of the mesh
        for(i=0; i<my_rows_B*my_cols_B; i++) {
            local_B[i] = rand() / RAND_DIV;
        }
        
    }
    for(j=0; j<proc_dims[1]; j++) {
        if(col_comms[j] != MPI_COMM_NULL)
            MPI_Bcast(local_B, my_rows_B*my_cols_B, MPI_FLOAT, 0, col_comms[j]);
    }

    //local_C has to be initialized by the processes at first column of the mesh. Its content has to be sent in broadcast on row_comms[i] (foreach i).
    if(all_cart_coords[my_rank][1] == 0) {  //i.e. if I am a process at first row of the mesh
        for(i=0; i<my_rows_C*my_cols_C; i++) {
            local_C[i] = rand() / RAND_DIV;
        }
        
    }
    for(j=0; j<proc_dims[0]; j++) {
        if(row_comms[j] != MPI_COMM_NULL)
            MPI_Bcast(local_C, my_rows_C*my_cols_C, MPI_FLOAT, 0, row_comms[j]);
    }


    /* END OF THE EXECUTION */
    MPI_Finalize();
    return 0;

}






/*
int main(int argc, char **argv) {
    //informazioni che il processo 0 deve inviare agli altri processi
    int big_k;  //K
    int big_n;  //N
    int big_m;  //M
    //matrici complete
    float total_A[M][K];
    float total_B[K][N];
    float total_C[M][N];
    //porzioni locali delle matrici
    float *local_A;
    float *local_B;
    float *local_C;
    float *intermediate_C;      //porzione della matrice C comprendente un'intera riga della mesh di processi; fa da intermediaria nella fase di gather dei risultati.
    //variabili MPI
    int my_rank;
    int p;
    MPI_Comm comm_world_copy;   //copia del comunicatore MPI_COMM_WORLD (è sempre buona norma averla)
    MPI_Comm comm_cart;         //nuovo comunicatore relativo alla griglia di processi da associare alla matrice C
    //variabili utili per suddividere il lavoro dei processi
    int m;                      //numero di righe per la sottomatrice (di C) da assegnare a ciascun processo
    int n;                      //numero di colonne per la sottomatrice (di C) da assegnare a ciascun processo
    int proc_dims[DIMS];        //array che indica il numero di processi che va a finire in ciascuna dimensione della matrice C (--> la mesh coinvolgerà la matrice C)
    int proc_dims_hat[DIMS];    //array che indica il numero di righe e colonne che, nella mesh di processi, sono composte da una riga/colonmna di C in più (le entry sono risp. M%m, N%n)
    //variabili di appoggio
    int periods[DIMS];          //array che indica se ciascuna dimensione della matrice deve essere periodica (i.e. circolare) o meno
    int mesh_row_index;         //variabile che tiene traccia dell'indice riga all'interno della mesh di processi
    int mesh_col_index;         //variabile che tiene traccia dell'indice colonna all'interno della mesh di processi
    int my_mesh_row;            //il PROPRIO indice riga all'interno della mesh di processi
    int my_mesh_col;            //il PROPRIO indice colonna all'interno della mesh di processi
    //indici ciclo for
    int i;
    int j;
    int l;

    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_copy); //duplicazione del comunicatore MPI_COMM_WORLD
    MPI_Comm_size(comm_world_copy, &p);
    MPI_Comm_rank(comm_world_copy, &my_rank);
    printf("Hello from process %d of %d!\n", my_rank, p);

    //inizializzazione di periods a soli false + inizializzazione di proc_dims a soli 0
    for(i=0; i<DIMS; i++) {
        periods[i] = NO;
        proc_dims[i] = 0;   //proc_dims[0] = #righe della mesh di processi; proc_dims[1] = #colonne della mesh di processi

    }

    //definizione (e successiva creazione) della topologia più quadrata possibile per i p processi
    MPI_Dims_create(p, DIMS, proc_dims);
    MPI_Cart_create(comm_world_copy, DIMS, proc_dims, periods, NO, &comm_cart);

    //definizione dei comunicatori relativi alle singole righe della mesh e alle singole colonne della mesh
    MPI_Comm *row_comms;
    MPI_Comm *col_comms;
    //definizione dei gruppi di processi (quello relativo al comunicatore comm_world_copy, quelli relativi alle singole righe della mesh e quelli relativi alle singole colonne della mesh)
    MPI_Group total_group;  //gruppo associato al comunicatore comm_world_copy
    MPI_Group *row_groups;  //gruppi associati ai singoli comunicatori relativi alle righe della mesh
    MPI_Group *col_groups;  //gruppi associati ai singoli comunicatori relativi alle colonne della mesh

    //allocazione di memoria per tutte le variabili sopra definite
    row_comms = (MPI_Comm *)malloc(proc_dims[0]*sizeof(MPI_Comm));
    col_comms = (MPI_Comm *)malloc(proc_dims[1]*sizeof(MPI_Comm));
    if(!(row_comms && col_comms)) {
        printf("Unable to allocate row_comms and col_comms.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }
    row_groups = (MPI_Group *)malloc(proc_dims[0]*sizeof(MPI_Group));
    col_groups = (MPI_Group *)malloc(proc_dims[1]*sizeof(MPI_Group));
    if(!(row_groups && col_groups)) {
        printf("Unable to allocate row_groups and col_groups.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }

    //definizione di tutti gli array di processi per ciascuna riga e per ciascuna colonna della mesh (i.e. per ciascun sotto-comunicatore)
    int **row_ranks_list;
    int **col_ranks_list;
    //array di appoggio che tengono traccia del prossimo indice da popolare per ciascun array contenuto nelle due ranks_list
    int *row_indexes_list;
    int *col_indexes_list;
    //variabile d'appoggio: ospita le coordinate di tutti i processi nell'ambito della griglia cartesiana
    int **all_cart_coords;

    //allocazione di memoria per tutte le variabili sopra definite
    //1) allocazione degli array
    row_indexes_list = (int *)malloc(proc_dims[0]*sizeof(int));
    col_indexes_list = (int *)malloc(proc_dims[1]*sizeof(int));
    if(!(row_indexes_list && col_indexes_list)) {
        printf("Unable to allocate row_indexes_list and col_indexes_list.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }
    //2) allocazione dei doppi puntatori per le matrici
    row_ranks_list = (int **)malloc(proc_dims[0]*sizeof(int *));
    col_ranks_list = (int **)malloc(proc_dims[1]*sizeof(int *));
    all_cart_coords = (int **)malloc(p*sizeof(int *));
    if(!(row_ranks_list && col_ranks_list && all_cart_coords)) {
        printf("Unable to allocate row_ranks_list, col_ranks_list and all_cart_coords.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }
    //3) allocazione delle entry delle matrici vere e proprie (con inizializzazione del puntatore relativo alla prima riga)
    row_ranks_list[0] = (int *)malloc(proc_dims[0]*proc_dims[1]*sizeof(int));
    col_ranks_list[0] = (int *)malloc(proc_dims[1]*proc_dims[0]*sizeof(int));
    all_cart_coords[0] = (int *)malloc(p*DIMS*sizeof(int));
    if(!(row_ranks_list[0] && col_ranks_list[0] && all_cart_coords[0])) {
        printf("Unable to allocate row_ranks_list[0], col_ranks_list[0] and all_cart_coords[0].\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }
    //4) inizializzazione dei puntatori relativi alle righe successive alla prima
    for(i=1; i<proc_dims[0]; i++)
        row_ranks_list[i] = row_ranks_list[0] + i*proc_dims[1];
    for(i=1; i<proc_dims[1]; i++)
        col_ranks_list[i] = col_ranks_list[0] + i*proc_dims[0];
    for(i=1; i<p; i++)
        all_cart_coords[i] = all_cart_coords[0] + i*DIMS;

    //inizializzazione a 0 di tutti i campi delle due indexes_list
    for(i=0; i<proc_dims[0]; i++) {
        row_indexes_list[i] = 0;
    }
    for(i=0; i<proc_dims[1]; i++) {
        col_indexes_list[i] = 0;
    }

    //ottenimento del gruppo dei processi che partecipano al comunicatore comm_world_copy (che sarebbero tutti i processi)
    MPI_Comm_group(comm_world_copy, &total_group);

    //ottenimento di tutti i gruppi associati alle singole righe della mesh
    for(i=0; i<p; i++) {
        MPI_Cart_coords(comm_cart, i, DIMS, all_cart_coords[i]);  //calcolo di delle cordinate per il processo i
        
        //popolamento dell'opportuno array contenuto in row_ranks_list; all_cart_coords[i][0] è l'indice riga del processo all'interno della mesh
        row_ranks_list[all_cart_coords[i][0]][row_indexes_list[all_cart_coords[i][0]]] = i;     //i = rank corrente
        row_indexes_list[all_cart_coords[i][0]]++;

        //popolamento dell'opportuno array contenuto in col_ranks_list; all_cart_coords[i][1] è l'indice colonna del processo all'interno della mesh
        col_ranks_list[all_cart_coords[i][1]][col_indexes_list[all_cart_coords[i][1]]] = i;     //i = rank corrente
        col_indexes_list[all_cart_coords[i][1]]++;

    }

    //ottenimento di tutti i sottogruppi legati alle singole righe della mesh e ai relativi sotto-comunicatori
    for(i=0; i<proc_dims[0]; i++) {
        //ricordiamo che abbiamo un numero di sottogruppi pari al numero di righe della mesh e un numero di processi in ogni sottogruppo pari al numero di colonne della mesh
        MPI_Group_incl(total_group, proc_dims[1], row_ranks_list[i], &row_groups[i]);
        MPI_Comm_create(comm_world_copy, row_groups[i], &row_comms[i]);   //i comunicatori qui definiti (in row_comms) verranno utilizzati nella funzione MPI_Scatterv().
    
    }
    //ottenimento di tutti i sottogruppi legati alle singole colonne della mesh e ai relativi sotto-comunicatori
    for(i=0; i<proc_dims[1]; i++) {
        //ricordiamo che abbiamo un numero di sottogruppi pari al numero di colonne della mesh e un numero di processi in ogni sottogruppo pari al numero di righe della mesh
        MPI_Group_incl(total_group, proc_dims[0], col_ranks_list[i], &col_groups[i]);
        MPI_Comm_create(comm_world_copy, col_groups[i], &col_comms[i]);   //i comunicatori qui definiti (in col_comms) verranno utilizzati nella funzione MPI_Scatterv().
    
    }
    //è possibile supporre che la prima entry di row_comms contenga il comunicatore relativo ai rappresentanti delle colonne della mesh,
    //mentre la prima entry di col_comms contenga il comunicatore relativo ai rappresentanti delle righe della mesh.


    //si assume che solo il processo 0 inizialmente conosca i valori di K, N, M e le tre matrici.
    if(my_rank == 0) {
        big_k = K;
        big_n = N;
        big_m = M;

        //inizializzazione randomica delle matrici
        for(i=0; i<big_m; i++) {    //matrice A
            for(j=0; j<big_k; j++) {
                total_A[i][j] = rand() % 10;

                #ifdef DEBUG
                printf("A[%d][%d] = %f\n", i, j, total_A[i][j]);
                #endif

            }

        }
        for(i=0; i<big_k; i++) {    //matrice B
            for(j=0; j<big_n; j++) {
                total_B[i][j] = rand() % 10;

                #ifdef DEBUG
                printf("B[%d][%d] = %f\n", i, j, total_B[i][j]);
                #endif

            }

        }
        for(i=0; i<big_m; i++) {    //matrice C
            for(j=0; j<big_n; j++) {
                total_C[i][j] = rand() % 10;

                #ifdef DEBUG
                printf("C[%d][%d] = %f\n", i, j, total_C[i][j]);
                #endif

            }

        }

    }

    //invio dei messaggi a tutti gli altri processi da parte del processo 0
    MPI_Bcast(&big_k, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&big_n, 1, MPI_INT, 0, comm_world_copy);
    MPI_Bcast(&big_m, 1, MPI_INT, 0, comm_world_copy);

    //calcolo delle dimensioni di base delle sottomatrici (di C) da assegnare a ciascun processo
    m = big_m/proc_dims[0];
    n = big_n/proc_dims[1];

    //calcolo del numero di righe e di colonne che, nella mesh di processi, sono composte da una riga/colonna di matrice in più
    proc_dims_hat[0] = big_m - m*proc_dims[0];
    proc_dims_hat[1] = big_n - n*proc_dims[1];


    //array con numero righe, numero colonne, spiazzamento righe e spiazzamento colonne della sottomatrice di ciascun processo
    int *row_sizes;
    int *row_displs;
    int *col_sizes;
    int *col_displs;
    //row_sizes e row_displs rispettivamente moltiplicate per il numero di colonne totale di A, per il numero di righe totale di B e per il numero di colonne totale di C
    int *block_sizes_A;
    int *block_displs_A;
    int *block_sizes_C;
    int *block_displs_C;

    //allocazione di memoria per tutte le variabili sopra definite
    row_sizes = (int *)malloc(proc_dims[0]*sizeof(int));
    row_displs = (int *)malloc(proc_dims[0]*sizeof(int));
    col_sizes = (int *)malloc(proc_dims[1]*sizeof(int));
    col_displs = (int *)malloc(proc_dims[1]*sizeof(int));
    if(!(row_sizes && row_displs && col_sizes && col_displs)) {
        printf("Unable to allocate row_sizes, row_displs, col_sizes and col_displs.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }
    block_sizes_A = (int *)malloc(proc_dims[0]*sizeof(int));
    block_displs_A = (int *)malloc(proc_dims[0]*sizeof(int));
    block_sizes_C = (int *)malloc(proc_dims[0]*sizeof(int));
    block_displs_C = (int *)malloc(proc_dims[0]*sizeof(int));
    if(!(block_sizes_A && block_displs_A && block_sizes_C && block_displs_C)) {
        printf("Unable to allocate block_sizes_A, block_displs_A, block_sizes_C and block_displs_C.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }

    //calcolo delle entry dei 4 array suddetti che indicano come verrà effettivamente suddiviso il lavoro (potrebbero esservi dei processi che prendono una riga e/o una colonna in più).
    mesh_row_index = 0;
    mesh_col_index = 0;

    for(i=0; i<p; i++) {
        //riga base + num righe
        if(all_cart_coords[i][0] < proc_dims_hat[0] && i%proc_dims[0] == 0) {    //caso in cui il processo i prende una riga di C in più
            row_displs[mesh_row_index] = all_cart_coords[i][0]*(m+1);
            row_sizes[mesh_row_index] = m+1;
            mesh_row_index++;

        } else if(all_cart_coords[i][0] >= proc_dims_hat[0] && i%proc_dims[0] == 0) {
            row_displs[mesh_row_index] = proc_dims_hat[0]*(m+1) + (all_cart_coords[i][0]-proc_dims_hat[0])*m;
            row_sizes[mesh_row_index] = m;
            mesh_row_index++;
        
        }
            
        //colonna base + num colonne
        if(all_cart_coords[i][1] < proc_dims_hat[1] && i < proc_dims[1]) {   //caso in cui il processo i prende una colonna di C in più
            col_displs[mesh_col_index] = all_cart_coords[i][1]*(n+1);
            col_sizes[mesh_col_index] = n+1;
            mesh_col_index++;
        
        } else if(all_cart_coords[i][1] >= proc_dims_hat[1] && i < proc_dims[1]) {
            col_displs[mesh_col_index] = proc_dims_hat[1]*(n+1) + (all_cart_coords[i][1]-proc_dims_hat[1])*n;
            col_sizes[mesh_col_index] = n;
            mesh_col_index++;

        }    

    }

    //inizializzazione del proprio indice riga e del proprio indice colonna all'interno della mesh di processi
    my_mesh_row = all_cart_coords[my_rank][0];
    my_mesh_col = all_cart_coords[my_rank][1];


    //allocazione della porzione locale delle tre matrici e di intermediate_C
    local_A = (float *)malloc(row_sizes[my_mesh_row]*big_k*sizeof(float));
    local_B = (float *)malloc(big_k*col_sizes[my_mesh_col]*sizeof(float));
    local_C = (float *)malloc(row_sizes[my_mesh_row]*col_sizes[my_mesh_col]*sizeof(float));
    intermediate_C = (float *)malloc(row_sizes[my_mesh_row]*big_n*sizeof(float));
    if(!(local_A && local_B && local_C && intermediate_C)) {
        printf("Unable to allocate local_A, local_B, local_C and intermediate_C.\n");
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    }

    //azzeramento esplicito di local_A, local_B, local_C e intermediate_C
    memset(local_A, 0, row_sizes[my_mesh_row]*big_k*sizeof(float));
    memset(local_B, 0, big_k*col_sizes[my_mesh_col]*sizeof(float));
    memset(local_C, 0, row_sizes[my_mesh_row]*col_sizes[my_mesh_col]*sizeof(float));
    memset(intermediate_C, 0, row_sizes[my_mesh_row]*big_n*sizeof(float));

    
    //il processo 0 invia al rappresentante di ciascuna riga di mesh l'intera matrice B e l'intera matrice C
    if(col_comms[0] != MPI_COMM_NULL) { //check su se il processo appartiene effettivamente al comunicatore col_comms[0]
        MPI_Bcast(total_B, big_k*big_n, MPI_INT, 0, col_comms[0]);
        MPI_Bcast(total_C, big_m*big_n, MPI_INT, 0, col_comms[0]);
    
    }
    //il processo 0 invia al rappresentante di ciascuna colonna di mesh l'intera matrice A
    if(row_comms[0] != MPI_COMM_NULL) { //check su se il processo appartiene effettivamente al comunicatore row_comms[0]
        MPI_Bcast(total_A, big_m*big_k, MPI_INT, 0, row_comms[0]);

    }

    //split della matrice C tra tutti i processi
    for(i=0; i<big_m; i++) { //nel caso della mesh deve essere distribuita una riga per volta.
        //calcolo del comunicatore (i.e. dell'insieme di processi) a cui deve essere distribuita la riga i-esima
        for(j=0; j<proc_dims[0]; j++) {
            if(i >= row_displs[j] && (j == proc_dims[0]-1 || i < row_displs[j+1])) { //la riga esatta della matrice deve appartenere al range di indici della corrispettiva riga di mesh.
                mesh_row_index = j;
                break;
            
            }

        }    

        if(row_comms[mesh_row_index] != MPI_COMM_NULL)  //check su se il processo appartiene effettivamente al comunicatore row_comms[mesh_row_index]
            MPI_Scatterv(total_C[i], col_sizes, col_displs, MPI_INT, &local_C[(i-row_displs[my_mesh_row])*col_sizes[my_mesh_col]], col_sizes[my_mesh_col], MPI_INT, 0, row_comms[mesh_row_index]);
                                                                
    }

    //split della matrice A tra tutti i processi
    //per far ciò, prima si inizializzano gli array block_sizes_A e block_displs_A.
    for(i=0; i<proc_dims[0]; i++) {
        block_sizes_A[i] = row_sizes[i]*big_k;
        block_displs_A[i] = row_displs[i]*big_k;

    }
    for(i=0; i<proc_dims[1]; i++) { //una Scatterv() per ogni sotto-comunicatore (i.e. per ogni colonna di mesh)
        if(col_comms[i] != MPI_COMM_NULL) //check su se il processo appartiene effettivamente al comunicatore col_comms[i]
            MPI_Scatterv(total_A, block_sizes_A, block_displs_A, MPI_INT, local_A, block_sizes_A[my_mesh_row], MPI_INT, 0, col_comms[i]);

    }

    //split della matrice B tra tutti i processi
    for(i=0; i<big_k; i++) {    //anche una matrice da suddividere vericalmente deve essere distribuita una riga per volta
        for(j=0; j<proc_dims[0]; j++) { //una Scatterv() per ogni sotto-comunicatore (i.e. per ogni riga di mesh)
            if(row_comms[j] != MPI_COMM_NULL)   //check su se il processo appartiene effettivamente al comunicatore row_comms[j]
                MPI_Scatterv(total_B[i], col_sizes, col_displs, MPI_INT, &local_B[i*col_sizes[my_mesh_col]], col_sizes[my_mesh_col], MPI_INT, 0, row_comms[j]);
        
        }

    }

    //calcolo dei risultati parziali di C <-- A*B + C
    for(i=0; i<row_sizes[my_mesh_row]; i++) {
        for(j=0; j<col_sizes[my_mesh_col]; j++) {
            for(l=0; l<big_k; l++) {
                local_C[i*col_sizes[my_mesh_col]+j] += local_A[i*big_k+l]*local_B[l*col_sizes[my_mesh_col]+j];

            }

        }

    }

    //ricollezionamento dei risultati parziali contestualmente alle singole righe della mesh di processi
    for(i=0; i<proc_dims[0]; i++) {
        for(j=0; j<row_sizes[i]; j++) {
            if(row_comms[i] != MPI_COMM_NULL)   //check su se il processo appartiene effettivamente al comunicatore row_comms[i]
                MPI_Gatherv(&local_C[j*col_sizes[my_mesh_col]], col_sizes[my_mesh_col], MPI_INT, &intermediate_C[j*big_n], col_sizes, col_displs, MPI_INT, 0, row_comms[i]);

        }
    
    }

    //a questo punto si mettono insieme i risultati parziali per ottenere il risultato finale dato da C <-- A*B + C
    //per far ciò, prima si inizializzano gli array block_sizes_C e block_displs_C.
    for(i=0; i<proc_dims[0]; i++) {
        block_sizes_C[i] = row_sizes[i]*big_n;
        block_displs_C[i] = row_displs[i]*big_n;

    }
    if(col_comms[0] != MPI_COMM_NULL)    //check su se il processo appartiene effettivamente al comunicatore col_comms[0]
        MPI_Gatherv(intermediate_C, block_sizes_C[my_mesh_row], MPI_INT, total_C, block_sizes_C, block_displs_C, MPI_INT, 0, col_comms[0]);

    //stampa del risultato finale
    if(my_rank == 0) {
        printf("RESULT:\n");
        for(i=0; i<big_m; i++) {
            for(j=0; j<big_n; j++) {
                printf("C[%d][%d] = %f\n", i, j, total_C[i][j]);

            }

        }

    }

    MPI_Finalize();
    return 0;

}
*/
