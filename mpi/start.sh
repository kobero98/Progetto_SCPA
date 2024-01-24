#! /bin/bash
mkdir target
mkdir result

module load mpi/openmpi-x86_64;
mpicc MPI-matrix-gen.c -o target/MPI-matrix-gen -O3
mpicc MPI-prod-kobero.c -o target/MPI-prod-kobero -O3
mpicc MPI-prod-fanfa.c -o target/MPI-prod-fanfa -lm -O3
cd target

seed=123456789

for K in {500,1000,5000,10000,15000}
do 
    dir="kobero-squared-$K"
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir $K $K $K
    mpirun -n 1 MPI-prod-kobero $dir > ../result/result-$K-1

    for p in {2,4,6,8,10,12,14,16,18,20}
    do
        echo "[" > ../result/result-$K-$p

        for n in {1..9}
        do
            echo "dimensione problema: $K numero processi: $p ";
            echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$K-$p
            mpirun -n $p MPI-prod-kobero $dir >> ../result/result-$K-$p
            echo "{}]}," >> ../result/result-$K-$p
        done
        echo "{\"iterazione\":10,\"vettore\":[" >> ../result/result-$K-$p
        mpirun -n $p MPI-prod-kobero $dir >> ../result/result-$K-$p
        echo "{}]}]" >> ../result/result-$K-$p

    done
done



for K in {32,64,128,256}
do
    dir="kobero-rectangular-$K"
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir 15000 15000 $K
    mpirun -n 1 MPI-prod-kobero $dir > ../result/result-$K-1

    for p in {2,4,6,8,10,12,14,16,18,20}
    do
        echo "[" > ../result/result-$K-$p

        for n in {1..9}
        do
            echo "dimensione problema: $K numero processi: $p ";
            echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$K-$p
            mpirun -n $p MPI-prod-kobero $dir >> ../result/result-$K-$p
            echo "{}]}," >> ../result/result-$K-$p
        done
        echo "{\"iterazione\":10,\"vettore\":[" >> ../result/result-$K-$p
        mpirun -n $p MPI-prod-kobero $dir >> ../result/result-$K-$p
        echo "{}]}]" >> ../result/result-$K-$p

    done
done



for K in {500,1000,5000,10000,15000}
do 
    dir="fanfa-squared-$K"
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir $K $K $K
    mpirun -n 1 MPI-prod-fanfa $dir 1 1 > ../result/result-$K-1-1

    for kb in {1,2,4,8,16}
    do

        for p in {2,4,6,8,10,12,14,16,18,20}
        do
            echo "[" > ../result/result-$K-$kb-$p

            for n in {1..9}
            do
                echo "dimensione problema: $K numero processi: $p ";
                echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$K-$kb-$p
                mpirun -n $p MPI-prod-fanfa $dir $kb $kb >> ../result/result-$K-$kb-$p
                echo "{}]}," >> ../result/result-$K-$kb-$p
            done
            echo "{\"iterazione\":10,\"vettore\":[" >> ../result/result-$K-$kb-$p
            mpirun -n $p MPI-prod-fanfa $dir $kb $kb >> ../result/result-$K-$kb-$p
            echo "{}]}]" >> ../result/result-$K-$kb-$p

        done
    done
done



for K in {32,64,128,256}
do 
    dir="fanfa-rectangular-$K"
    mkdir $dir
    mpirun -n 3 MPI-matrix-gen $seed $dir 15000 15000 $K
    mpirun -n 1 MPI-prod-fanfa $dir 1 1 > ../result/result-$K-1-1

    for kb in {1,2,4,8,16}
    do

        for p in {2,4,6,8,10,12,14,16,18,20}
        do
            echo "[" > ../result/result-$K-$kb-$p

            for n in {1..9}
            do
                echo "dimensione problema: $K numero processi: $p ";
                echo "{\"iterazione\":$n,\"vettore\":[" >> ../result/result-$K-$kb-$p
                mpirun -n $p MPI-prod-fanfa $dir $kb $kb >> ../result/result-$K-$kb-$p
                echo "{}]}," >> ../result/result-$K-$kb-$p
            done
            echo "{\"iterazione\":10,\"vettore\":[" >> ../result/result-$K-$kb-$p
            mpirun -n $p MPI-prod-fanfa $dir $kb $kb >> ../result/result-$K-$kb-$p
            echo "{}]}]" >> ../result/result-$K-$kb-$p

        done
    done
done
