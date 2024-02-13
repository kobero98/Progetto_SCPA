#! /bin/bash
make

mkdir squared
mkdir rectangular
mkdir special

# echo "[" > squared/file_cuda_prod.json
# echo "[" > squared/file_cuda_prod_v2.json
# echo "[" > squared/file_cuda_prod_ms.json
# for XBD in {1,2,4,16,32,64,128,256,512,1024}
# do    
#     YBD=$((1024 / $XBD))
#     echo "$XBD, $YBD"
#     for K in {500,750,1000,2000,5000,10000}
#     do
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K, " >> squared/file_cuda_prod.json
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K, " >> squared/file_cuda_prod_v2.json
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K, " >> squared/file_cuda_prod_ms.json
#         ./Cuda-prod $K $K $K n $XBD $YBD >> squared/file_cuda_prod.json
#         ./Cuda-prod-v2 $K $K $K n $XBD $YBD >> squared/file_cuda_prod_v2.json
#         ./Cuda-prod-ms-read $K $K $K n $XBD $YBD >> squared/file_cuda_prod_ms.json
#         echo "}," >> squared/file_cuda_prod.json
#         echo "}," >> squared/file_cuda_prod_v2.json
#         echo "}," >> squared/file_cuda_prod_ms.json
#     done
# done
# echo "]" >> squared/file_cuda_prod.json
# echo "]" >> squared/file_cuda_prod_v2.json
# echo "]" >> squared/file_cuda_prod_ms.json



# echo "[" > rectangular/file_cuda_prod.json
# echo "[" > rectangular/file_cuda_prod_v2.json
# echo "[" > rectangular/file_cuda_prod_ms.json
# for XBD in {1,2,4,16,32,64,128,256,512,1024}
# do    
#     YBD=$((1024 / $XBD))
#     echo "$XBD, $YBD"
#     for K in {32,64,128,256}
#     do
#         M=15000 #scegliamolo multiplo di due per comodità
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> rectangular/file_cuda_prod.json
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> rectangular/file_cuda_prod_v2.json
#         echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> rectangular/file_cuda_prod_ms.json
#         ./Cuda-prod $M $K $M n $XBD $YBD >> rectangular/file_cuda_prod.json
#         ./Cuda-prod-v2 $M $K $M n $XBD $YBD >> rectangular/file_cuda_prod_v2.json
#         ./Cuda-prod-ms-read $M $K $M n $XBD $YBD >> rectangular/file_cuda_prod_ms.json
#         echo "}," >> rectangular/file_cuda_prod.json
#         echo "}," >> rectangular/file_cuda_prod_v2.json
#         echo "}," >> rectangular/file_cuda_prod_ms.json
#     done
# done
# echo "]" >> rectangular/file_cuda_prod.json
# echo "]" >> rectangular/file_cuda_prod_v2.json
# echo "]" >> rectangular/file_cuda_prod_ms.json



echo "[" > special/file_cuda_prod.json
echo "[" > special/file_cuda_prod_v2.json
echo "[" > special/file_cuda_prod_ms.json
for XBD in {1,2,4,16,32,64,128,256,512,1024}
do    
    YBD=$((1024 / $XBD))
    echo "$XBD, $YBD"
    for K in {32,64,128,256,512,1024}
    do
        for M in {512,1024,2048,4096}
        do
            echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> special/file_cuda_prod.json
            echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> special/file_cuda_prod_v2.json
            echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$M,\"K\":$K, " >> special/file_cuda_prod_ms.json
            ./Cuda-prod $M $K $M n $XBD $YBD >> special/file_cuda_prod.json
            ./Cuda-prod-v2 $M $K $M n $XBD $YBD >> special/file_cuda_prod_v2.json
            ./Cuda-prod-ms-read $M $K $M n $XBD $YBD >> special/file_cuda_prod_ms.json
            echo "}," >> special/file_cuda_prod.json
            echo "}," >> special/file_cuda_prod_v2.json
            echo "}," >> special/file_cuda_prod_ms.json
        done
    done
done
echo "]" >> special/file_cuda_prod.json
echo "]" >> special/file_cuda_prod_v2.json
echo "]" >> special/file_cuda_prod_ms.json
