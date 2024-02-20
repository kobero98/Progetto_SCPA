#! /bin/bash
make

mkdir Error_Speed_up

echo "[" > squared/file_cuda_prod.json
echo "[" > squared/file_cuda_prod_v2.json
echo "[" > squared/file_cuda_prod_ms.json
echo "[" > squared/file_cuda_prod_ms_v2.json
for XBD in {512,1024} #{1,2,4,16,32,64,128,256,512,1024}
do    
    YBD=$((1024 / $XBD))
    echo "$XBD, $YBD"
    for K in {500,750,1000,2000,3000,4000}
    do
        echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K,\"result\":[ " >> squared/file_cuda_prod.json
        echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K,\"result\":[ " >> squared/file_cuda_prod_v2.json
        echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K,\"result\":[ " >> squared/file_cuda_prod_ms.json
        echo "{ \"XDB\":$XBD,\"YDB\":$YBD,\"M\":$K,\"K\":$K,\"result\":[ " >> squared/file_cuda_prod_ms_v2.json
        for i in {1..9}
        do
            echo "{" >> squared/file_cuda_prod.json
            echo "{" >> squared/file_cuda_prod_v2.json
            echo "{" >> squared/file_cuda_prod_ms.json
            echo "{" >> squared/file_cuda_prod_ms_v2.json
            
            ./Cuda-prod $K $K $K y $XBD $YBD >> squared/file_cuda_prod.json
            ./Cuda-prod-v2 $K $K $K y $XBD $YBD >> squared/file_cuda_prod_v2.json
            ./Cuda-prod-ms-read $K $K $K y $XBD $YBD >> squared/file_cuda_prod_ms.json
            ./Cuda-prod-ms-v2 $K $K $K y >> squared/file_cuda_prod_ms_v2.json
            
            echo "}," >> squared/file_cuda_prod.json
            echo "}," >> squared/file_cuda_prod_v2.json
            echo "}," >> squared/file_cuda_prod_ms.json
            echo "}," >> squared/file_cuda_prod_ms_v2.json

        done
        echo "{" >> squared/file_cuda_prod.json
        echo "{" >> squared/file_cuda_prod_v2.json
        echo "{" >> squared/file_cuda_prod_ms.json
        echo "{" >> squared/file_cuda_prod_ms_v2.json
        
        ./Cuda-prod $K $K $K n $XBD $YBD >> squared/file_cuda_prod.json
        ./Cuda-prod-v2 $K $K $K n $XBD $YBD >> squared/file_cuda_prod_v2.json
        ./Cuda-prod-ms-read $K $K $K n $XBD $YBD >> squared/file_cuda_prod_ms.jsonv
        ./Cuda-prod-ms-v2 $K $K $K n >> squared/file_cuda_prod_ms_v2.json

        echo "}]" >> squared/file_cuda_prod.json
        echo "}]" >> squared/file_cuda_prod_v2.json
        echo "}]" >> squared/file_cuda_prod_ms.json
        echo "}]" >> squared/file_cuda_prod_ms_v2.json
        
        
        echo "}," >> squared/file_cuda_prod.json
        echo "}," >> squared/file_cuda_prod_v2.json
        echo "}," >> squared/file_cuda_prod_ms.json
        echo "}," >> squared/file_cuda_prod_ms_v2.json
    done
done
truncate -s -2 squared/file_cuda_prod.json
truncate -s -2 squared/file_cuda_prod_v2.json
truncate -s -2 squared/file_cuda_prod_ms.json
truncate -s -2 squared/file_cuda_prod_ms_v2.json

echo "]" >> squared/file_cuda_prod.json
echo "]" >> squared/file_cuda_prod_v2.json
echo "]" >> squared/file_cuda_prod_ms.json
echo "]" >> squared/file_cuda_prod_ms_v2.json

