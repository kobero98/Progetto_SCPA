#! /bin/bash
mkdir target
mkdir result

module load mpi/openmpi-x86_64;
mpicc MPI-matrix-gen.c -o target/MPI-matrix-gen -O3
mpicc MPI-prod-$1.c -o target/MPI-prod-$1 -O3
cd target

seed=123456789
for K in {10,50,100,1000,2000,3000,4000}
do 
    dir="matrix-seed-$seed-$K"
    echo ciao $dir $seed
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir $K $K $K
    mpirun -n 1 MPI-prod-$1 $dir > ../result/result-$1-$K-1

    for p in {5,10,15,20}
    do
        echo "[" > ../result/result-$1-$K-$p

        for n in {1..19}
        do
            echo "dimension problema: $K numero processi: $p ";
            echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$1-$K-$p
            mpirun -n $p MPI-prod-$1 $dir >> ../result/result-$1-$K-$p
            echo "{}]}," >> ../result/result-$1-$K-$p
        done
        echo "{\"iterazione\":20,\"vettore\":[" >> ../result/result-$1-$K-$p
        mpirun -n $p MPI-prod-$1 $dir >> ../result/result-$1-$K-$p
        echo "{}]}]" >> ../result/result-$1-$K-$p

    done
done

for K in {5000,8000,10000,11000,12000,15000}
do
    dir="matrix-seed-$seed-$K"
    echo ciao $dir $seed
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir $K $K $K

    for p in {4,8,16,20}
    do
        echo "[" > ../result/result-$1-$K-$p

        for n in {1..19}
        do
            echo "dimension problema: $K numero processi: $p ";
            echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$1-$K-$p
            mpirun -n $p MPI-prod-$1 $dir >> ../result/result-$1-$K-$p
            echo "{}]}," >> ../result/result-$1-$K-$p
        done
        echo "{\"iterazione\":20,\"vettore\":[" >> ../result/result-$1-$K-$p
        mpirun -n $p MPI-prod-$1 $dir >> ../result/result-$1-$K-$p
        echo "{}]}]" >> ../result/result-$1-$K-$p

    done
done
