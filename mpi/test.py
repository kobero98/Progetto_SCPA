import json

# f=open("result/result-kobero-10-10")
# j=json.load(f)
# print(j)
dictResultTime={}
for i in [10,50,100,1000,2000]:
    f=open("result/result-kobero-"+str(i)+"-1")
    j=json.load(f)
    dictResultTime[i]=j["tempo_senza_creazione"];
    print(j)
#per ogni iterazione prendo il tempo massimo di esecuzione e inseguito faccio una media per le 20 iterazioni questo lo faccio per ogni matrice e per ogni numero di processi eseguitiù
dictionaryMediaValore={}
for i in [5,10,15,20]:
    for j in [10,50,100,1000,2000]:
        f=open("result/result-kobero-"+str(j)+"-"+str(i))
        js=json.load(f)
        somma=0
        for k in range(1,21):
            max= 0.0
            for v in js[k]:
                    print(i,j,k,js[k]," ciao ",v)
                    for x in v['vettore']:
                        if(max<x["tempo_senza_creazione"]):
                            max=x["tempo_senza_creazione"]
            somma=somma+max
        dictionaryMediaValore[(j,i)]=somma/20
print(dictionaryMediaValore)