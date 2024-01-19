#! /bin/bash
mkdir target
mkdir result

module load mpi;
mpicc matrixGeneratorMPI.c -o target/mgmpi
mpicc MPI-prod-kobero.c -o target/MPI-prod-kobero
cd target
seed=123456789
for K in {10,50,100,1000,2000,3000,5000}
do 
    dir="matrix_Seed_$seed-$K"
    echo ciao $dir $seed
    mkdir $dir
    mpirun -n 3 mgmpi $seed $dir $K $K $K
    for p in {1,5,10,15,20}
    do
        for n in {1}
        do
            echo "dimension problema: $K numero processi: $p ";  
            mpirun -n $p MPI-prod-kobero $dir > ../result/result_Kobero_$K-$p
        done
    done
done
