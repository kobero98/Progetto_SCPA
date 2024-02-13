#! /bin/bash

make

echo "["  >  file_cuda_prod.json
echo "["  >  file_cuda_prod_ms.json
for XBD in {1,2,4,16,32,64,128,256,512,1024}
do    
    YBD=$((1024 / $XBD ))
    echo "$XBD, $YBD"
    for K in {32,64,128,256,512,1024}
    do
        for M in {512,1024,2048,4096}
        do
            echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, "  >>  file_cuda_prod.json
            echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, "  >>  file_cuda_prod_ms.json
            ./Cuda-prod $M $K $M n $XBD $YBD >> file_cuda_prod.json
            ./Cuda-prod-ms-read $M $K $M "n" $XBD $YBD >>file_cuda_prod_ms.json
            echo "},"  >>  file_cuda_prod.json
            echo "},"  >>  file_cuda_prod_ms.json
        done
    done
done
echo "]"  >>  file_cuda_prod.json
echo "]"  >>  file_cuda_prod_ms.json