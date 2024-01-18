#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define DIMS 2  //number of dimensions of a matrix
#define RAND_DIV 1000000000.0   //divisor of rand() when matrixes components are generated

#define NO 0
#define YES NO+1

#define DEBUG



//TODO: execution of code by varying num of processes on course server

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
    //local auxiliary variables (loop indexes + max_proc_taken)
    int i;
    int j;
    int l;
    int max_proc_taken; //variable useful to keep track of the last process registered in i-th row
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
        max_proc_taken = -1;
        for(j=0; j<num_mesh_cols; j++) {
            for(l=max_proc_taken+1; l<p; l++) { //we iterate from the process succeeding max_proc_taken to not find the same process foreach iteration.
                if(all_cart_coords[l][0] == i) {
                    row_ranks_list[i][j] = l;   //l = rank of the found process having row index i in the mesh
                    max_proc_taken = l;
                    break;  //once we have populated the j-th entry, we can stop to scan all_cart_coords.
                }
            }
        }
    }

    return row_ranks_list;

}





/* function which has the responsibility to allocate the array pointing to all the arrays of process ranks of each mesh column */
int **create_col_ranks_list(int num_mesh_rows, int num_mesh_cols, int p, int **all_cart_coords, MPI_Comm comm_world_copy) {
    //local auxiliary variables (loop indexes + max_proc_taken)
    int i;
    int j;
    int l;
    int max_proc_taken; //variable useful to keep track of the last process registered in i-th column
    //function output
    int **col_ranks_list;

    //allocation of col_ranks_list
    col_ranks_list = (int **)malloc(num_mesh_cols*sizeof(int *));
    if(!col_ranks_list)
        MPI_Abort(comm_world_copy, EXIT_FAILURE);

    for(i=0; i<num_mesh_cols; i++) {
        //allocation of col_ranks_list[i]
        col_ranks_list[i] = (int *)malloc(num_mesh_rows*sizeof(int));
        if(!col_ranks_list[i])
            MPI_Abort(comm_world_copy, EXIT_FAILURE);

        /* Here we are populating the ranks list associated to i-th column of the mesh.
         * Index j represents the j-th entry of this rank list.
         * Index l is only useful to iterate on all_cart_coords array to check which processes are exactly in the i-th column of the mesh.
         */
        max_proc_taken = -1;
        for(j=0; j<num_mesh_rows; j++) {
            for(l=max_proc_taken+1; l<p; l++) { //we iterate from the process succeeding max_proc_taken to not find the same process foreach iteration.
                if(all_cart_coords[l][1] == i) {
                    col_ranks_list[i][j] = l;   //l = rank of the found process having column index i in the mesh
                    max_proc_taken = l;
                    break;  //once we have populated the j-th entry, we can stop to scan all_cart_coords.
                }
            }
        }
    }
    
    return col_ranks_list;

}





/* function which has the responsibility to calculate the number of rows (or, alternatively, the number of columns) of matrix A 
 * ruled by the current process. If dim == 0, we are calculating the number of rows; else, we are calculating the number of columns.
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

    //check important in order to be sure the code works fine
    if(dim != 0)
        dim = 1;

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





/* function which returns the minimum between two integers */
int get_min(int a, int b) {
    if(a > b)
        return b;
    else
        return a;
    
}





/* function which frees a buffer (passed as input) which points to num_entries other buffers */
void free_all(int **buffer, int num_entries) {
    //local auxiliary variable (loop index)
    int i;

    for(i=0; i<num_entries; i++) {
        free(buffer[i]);
    }

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
    //total matrix C (it has to be allocated only by process with rank 0)
    float *C;
    //local portions of the matrixes
    float *local_A;
    float *local_B;
    float *local_C;
    //auxiliary buffer for matrix C exchange (it has to be allocated only by the processes at column 0 of the mesh)
    float *other_C;

    /* MPI VARIABLES */
    int my_rank;
    int p;                      //#processes
    int tag;                    //tag for MPI_Send() and MPI_Recv() operations
    MPI_Status status;          //status of msg reception during MPI_Recv()
    MPI_Comm comm_world_copy;   //copy of MPI_COMM_WORLD
    MPI_Comm comm_cart;         //new communicator related to processes topology, which is related to a super-block in matrix A
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
    //counters
    int num_recv;   //variable counting the number of times process 0 is invoking MPI_Recv() during the final gather operation. It helps preventing deadlocks.
    int max_recv;   //variable indicating the max number of times process 0 has to invoke MPI_Recv() during the final gather operation.

    /* LOOP INDEXES */
    int i;
    int j;
    int l;



    /* MPI INITIALIZATION */
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_copy); //duplication of MPI_COMM_WORLD communicator
    MPI_Comm_size(comm_world_copy, &p);
    MPI_Comm_rank(comm_world_copy, &my_rank);
    tag = 0;
    //plant seed for successive rand() invokations
    srand(my_rank+1);   //srand(0) == srand(1)



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
    my_rows_A = get_local_dim_A(M, mb, proc_dims, all_cart_coords, my_rank, 0); //last parameter == 0: we are counting rows number.
    my_cols_A = get_local_dim_A(K, kb, proc_dims, all_cart_coords, my_rank, 1); //last parameter == 1: we are counting columns number.
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

    //allocate auxiliary buffer for matrix C exchange; only processes at column 0 of the mesh have to perform this allocation.
    if(all_cart_coords[my_rank][1] == 0) {
        other_C = (float *)malloc(my_rows_C*my_cols_C*sizeof(float));
        if(!other_C)
            MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }

    //allocate total matrix C; only process 0 has to perform this allocation.
    if(my_rank == 0) {
        C = (float *)malloc(N*M*sizeof(float));
        if(!C)
            MPI_Abort(comm_world_copy, EXIT_FAILURE);
    }

    //local_A has to be initialized by all the processes.
    for(i=0; i<my_rows_A*my_cols_A; i++) {
        local_A[i] = rand() / RAND_DIV;
    }

    //local_B has to be initialized by the processes at first row of the mesh. Its content has to be sent in broadcast on col_comms[i] (foreach i).
    if(all_cart_coords[my_rank][0] == 0) {  //i.e. if I am a process at first row of the mesh --> i.e. if my_mesh_row == 0
        for(i=0; i<my_rows_B*my_cols_B; i++) {
            local_B[i] = rand() / RAND_DIV;
        }
        
    }
    for(j=0; j<proc_dims[1]; j++) {
        if(col_comms[j] != MPI_COMM_NULL)
            MPI_Bcast(local_B, my_rows_B*my_cols_B, MPI_FLOAT, 0, col_comms[j]);
    }

    //local_C has to be initialized by the processes at first column of the mesh. Its content has NOT to be sent in broadcast on row_comms[i] (foreach i)
    //because in the end the different outputs of local_C for the different processes have to be summed and intial value of matrix C shall be counted only once.
    if(all_cart_coords[my_rank][1] == 0) {  //i.e. if I am a process at first row of the mesh --> i.e. if my_mesh_col == 0
        for(i=0; i<my_rows_C*my_cols_C; i++) {
            local_C[i] = rand() / RAND_DIV;
        }
    } else {
        memset(local_C, 0, my_rows_C*my_cols_C*sizeof(float));
    }



    /* DEBUG PRINTS */
    #ifdef DEBUG
    if(my_rank == 0) {
        printf("PROCESS %d\n", my_rank);
        printf("\nMATRIX A\n");
        for(i=0; i<my_rows_A; i++) {
            for(j=0; j<my_cols_A; j++) {
                printf("%f  ", local_A[i*my_cols_A+j]);
            }
            printf("\n");
        }
        printf("\nMATRIX B\n");
        for(i=0; i<my_rows_B; i++) {
            for(j=0; j<my_cols_B; j++) {
                printf("%f  ", local_B[i*my_cols_B+j]);
            }
            printf("\n");
        }
        printf("\nMATRIX C (INPUT)\n");
        for(i=0; i<my_rows_C; i++) {
            for(j=0; j<my_cols_C; j++) {
                printf("%f  ", local_C[i*my_cols_C+j]);
            }
            printf("\n");
        }
    }
    #endif



    /* CALCULATE PARTIAL RESULTS OF C <-- A*B + C */
    for(i=0; i<my_rows_A; i++) {   //TODO: which loop ordering is the most efficient?
        for(j=0; j<my_cols_B; j++) {
            for(l=0; l<my_cols_A; l++) {
                local_C[i*my_cols_C+j] += local_A[i*my_cols_A+l]*local_B[l*my_cols_B+j];
            }
        }
    }



    /* DEBUG PRINTS */
    #ifdef DEBUG
    if(my_rank == 0) {
        printf("\nMATRIX C (LOCAL OUTPUT)\n");
        for(i=0; i<my_rows_C; i++) {
            for(j=0; j<my_cols_C; j++) {
                printf("%f  ", local_C[i*my_cols_C+j]);
            }
            printf("\n");
        }
    }
    #endif


    
    /* GATHER OF MATRIX C (OUTPUT) */
    //PHASE 1: processes on the same mesh row have the outputs of matrix C related to the same matrix entries: they simply have to be summed.
    for(i=0; i<proc_dims[0]; i++) { //loop on all the mesh rows
        if(all_cart_coords[my_rank][0] == i) {  //condition: belonging to i-th mesh row
            for(j=1; j<proc_dims[1]; j++) { //loop on all the processes belonging to the same mesh row (except itself)

                if(all_cart_coords[my_rank][1] == j && row_comms[i] != MPI_COMM_NULL) {   //one process at time has to send its output to the row mesh representative.
                    MPI_Send(local_C, my_rows_C*my_cols_C, MPI_FLOAT, 0, tag, row_comms[i]);
                }
                else if(all_cart_coords[my_rank][1] == 0 && row_comms[i] != MPI_COMM_NULL) {
                    //row mesh representative receives the output and sums it with local_C (iterating on all the processes of mesh row).
                    MPI_Recv(other_C, my_rows_C*my_cols_C, MPI_FLOAT, j, tag, row_comms[i], &status);
                    //sum: local_C = local_C + other_C
                    for(l=0; l<my_rows_C*my_cols_C; l++) {
                        local_C[l] += other_C[l];
                    }
                }
            }
        }
    }

    //PHASE 2: processes on different mesh rows have the outputs of matrix C related to different matrix entries: they have to be gathered (again with send & recv).
    max_recv = ceil(1.0*M/mb);  //max_recv corresponds to the number of strips in which matrix C is divided.
    num_recv = 0;
    for(i=0; i<proc_dims[0]; i++) {  //loop on all the processes belonging to the first mesh column (i.e. loop on all the mesh rows)
        if(all_cart_coords[my_rank][1] == 0) {  //condition: belonging to first mesh column
            
            for(j=0; j<ceil(1.0*my_rows_C/mb); j++) {   //loop on all the strips of matrix C belonging to the process on which we are iterating (ceil==roof)
                if(all_cart_coords[my_rank][0] == i && col_comms[0] != MPI_COMM_NULL) {   //one process at time sends its local_C to process 0 (included process 0 itself).
                    MPI_Send(&local_C[j*N*mb], N*get_min(mb, my_rows_C-j*mb), MPI_FLOAT, 0, tag, col_comms[0]); //my_rows_C-j*mb = "remains"
                } 
                if(all_cart_coords[my_rank][0] == 0 && col_comms[0] != MPI_COMM_NULL) { //process 0 receives local_C from the sender and gathers it into matrix C.
                    if(num_recv == max_recv)    //if process 0 has received all the matrix C strips, it shall not invoke MPI_Recv() anymore.
                        break;
                    MPI_Recv(&C[(j*proc_dims[0]+i)*mb*N], N*get_min(mb, my_rows_C-j*mb), MPI_FLOAT, i, tag, col_comms[0], &status); //my_rows-j*mb = "remains"
                    //(j*proc_dims[0]+i)*mb*N = initial offset of matrix C from which we have to copy strip j of local copy of C associated to i-th process
                    num_recv++;
                }
            }
        }
    }



    /* DEBUG PRINTS */
    #ifdef DEBUG
    if(my_rank == 0) {
        printf("\nMATRIX C (GLOBAL OUTPUT)\n");
        for(i=0; i<M; i++) {
            for(j=0; j<N; j++) {
                printf("%f  ", C[i*N+j]);
            }
            printf("\n");
        }
    }
    #endif



    /* END OF THE EXECUTION */
    free_all(all_cart_coords, p);
    free_all(row_ranks_list, proc_dims[0]);
    free_all(col_ranks_list, proc_dims[1]);

    free(row_groups);
    free(col_groups);
    free(row_comms);
    free(col_comms);

    free(local_A);
    free(local_B);
    free(local_C);

    if(all_cart_coords[my_rank][1] == 0)
        free(other_C);
    if(my_rank == 0)
        free(C);

    MPI_Finalize();
    return 0;

}
