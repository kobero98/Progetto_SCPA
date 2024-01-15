#! /bin/bash
cd target;
for K in {100,1000,2000,3000,5000}
do 
    echo $K+" tentativo";  
    mpirun MPI-prod-kobero $K $K $K > ../result/result_Kobero_$K
donels
