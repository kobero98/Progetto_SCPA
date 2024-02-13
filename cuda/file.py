import pandas as pd
from pandasql import sqldf

df = pd.read_json("file_cuda_prod.json")
df2 = pd.read_json("file_cuda_prod_ms.json")

df1 = df[df["GFLOPS"]<1000]
df2 = df2[df2["GFLOPS"]<1000]
print(df1.shape)
print(df2.shape)

# Unione dei DataFrame in base alle colonne XDB e YDB
merged_df = pd.merge(df1, df2, on=['XDB', 'YDB','M','K'], suffixes=('_df1', '_df2'))

# Calcolo del valore massimo di GFLOPS per ogni coppia (XDB, YDB)
merged_df['Max_GFLOPS'] = merged_df[['GFLOPS_df1', 'GFLOPS_df2']].max(axis=1)

# Confronto tra i valori massimi per determinare quale DataFrame ha il valore maggiore
merged_df['Max_GFLOPS_df'] = merged_df.apply(lambda row: 'df1' if row['Max_GFLOPS'] == row['GFLOPS_df1'] else 'df2', axis=1)

# Stampa del DataFrame risultante
print(merged_df)
