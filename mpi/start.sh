#! /bin/bash
mkdir target
mkdir result

module load mpi;
mpicc matrixGeneratorMPI.c -o target/mgmpi
mpicc MPI-prod-kobero.c -o target/MPI-prod-kobero
cd target
seed=123456789
for K in {10,50,100,1000,2000}
do 
    dir="matrix_Seed_$seed-$K"
    echo ciao $dir $seed
    mkdir $dir
    mpirun -n 3 mgmpi $seed $dir $K $K $K
    mpirun -n 1 MPI-prod-kobero $dir > ../result/result_Kobero_$K-1
    for p in {5,10,15,20}
    do
        echo "[" > ../result/result_Kobero_$K-$p
        for n in {1..19}
        do
            echo "dimension problema: $K numero processi: $p ";
            echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result_Kobero_$K-$p
            mpirun -n $p MPI-prod-kobero $dir >> ../result/result_Kobero_$K-$p
            echo "{}]}," >> ../result/result_Kobero_$K-$p
        done
        echo "{\"iterazione\":20,\"vettore\":[" >> ../result/result_Kobero_$K-$p
        mpirun -n $p MPI-prod-kobero $dir >> ../result/result_Kobero_$K-$p
        echo "{}]}]" >> ../result/result_Kobero_$K-$p
    done

done
