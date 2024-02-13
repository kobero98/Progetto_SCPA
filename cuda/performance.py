import pandas as pd
import sys

from pandasql import sqldf

def main():
    # Check su se il programma è stato invocato con il parametro adeguato
    if len(sys.argv) < 2 or not (sys.argv[1] == "squared" and sys.argv[1] == "rectangular" and sys.argv[1] == "special"):
        print("Usage: python performance.py <folder_name>")
        exit(-1)

    df1 = pd.read_json(sys.argv[1] + "/file_cuda_prod.json")    #df1 è relativo all'esecuzione SENZA shared memory.
    df2 = pd.read_json(sys.argv[1] + "/file_cuda_prod_ms.json") #df2 è relativo all'esecuzione CON shared memory.

    # Eliminazione delle entry di df2 relative alle esecuzioni che richiedono più shared memory del massimo consentito
    df2 = df2[df2["K"]*df2["YDB"]*8<49000]

    # Unione dei DataFrame in base alle colonne XDB, YBD, M, K
    merged_df = pd.merge(df1, df2, on=['XDB','YDB','M','K'], suffixes=('_df1', '_df2'))

    # Calcolo del valore massimo di GFLOPS per ogni quadrupla (XDB, YDB, M, K)
    merged_df['Max_GFLOPS'] = merged_df[['GFLOPS_df1', 'GFLOPS_df2']].max(axis=1)

    # Confronto tra i valori massimi per determinare quale DataFrame ha il valore maggiore
    merged_df['Max_GFLOPS_df'] = merged_df.apply(lambda row: 'df1' if row['Max_GFLOPS'] == row['GFLOPS_df1'] else 'df2', axis=1)

    # Calcolo della differenza tra il numero di GFLOPS nel caso senza shared memory e il numero di GFLOPS nel caso con shared memory
    merged_df["DIFF_df"]=abs(merged_df["GFLOPS_df2"]-merged_df["GFLOPS_df1"])

    """
    merged_df=merged_df[merged_df["DIFF_df"]<2]
    tabellaDF1 = merged_df[merged_df["Max_GFLOPS_df"]=="df1"]
    tabellaDF2 = merged_df[merged_df["Max_GFLOPS_df"]=="df2"]
    """

    sorted_df=merged_df.sort_values("Max_GFLOPS") 
    # Stampa del DataFrame risultante
    print(sorted_df)



if __name__ == "__main__":
    main()
