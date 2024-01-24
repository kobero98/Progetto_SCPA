import json
import matplotlib.pyplot as plt

# f=open("result/result-10-10")
# j=json.load(f)
# print(j)
dictResultTime={}
for i in [10,50,100,1000,2000]:
    f=open("result/result-"+str(i)+"-1")
    j=json.load(f)
    dictResultTime[i]=j["tempo_senza_creazione"];
    print(j)
#per ogni iterazione prendo il tempo massimo di esecuzione e inseguito faccio una media per le 20 iterazioni questo lo faccio per ogni matrice e per ogni numero di processi eseguiti
dictionaryMediaValore={}
for i in [5,10,15,20]:
    for j in [10,50,100,1000,2000]:
        f=open("result/result-"+str(j)+"-"+str(i))
        js=json.load(f)
        somma=0
        for k in range(1,20):
            max= 0.0
            for v in js[k]["vettore"]:
                    if "tempo_senza_creazione" in v and max<v["tempo_senza_creazione"]:
                        max=v["tempo_senza_creazione"]
            somma=somma+max
        dictionaryMediaValore[(j,i)]=somma/20
print(dictionaryMediaValore)


def genera_dati_grafico(dizionario):
    dimensioni_matrice = []
    tempi_esecuzione = []

    for config, tempo in dizionario.items():
        num_processi, dim_matrice = config
        dimensioni_matrice.append(dim_matrice)
        tempi_esecuzione.append(tempo)

    return dimensioni_matrice, tempi_esecuzione

def disegna_grafici(dizionario):
    dimensioni_matrice, tempi_esecuzione = genera_dati_grafico(dizionario)

    # Grafico con l'aumentare delle dimensioni
    plt.figure(figsize=(10, 5))
    for num_processi in set([config[0] for config in dizionario.keys()]):
        tempi_processi = [tempi_esecuzione[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][0] == num_processi]
        dimensioni_processi = [dimensioni_matrice[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][0] == num_processi]
        plt.plot(dimensioni_processi, tempi_processi, label=f'{num_processi} processi')

    plt.xlabel('Dimensione Matrice')
    plt.ylabel('Tempo di Esecuzione')
    plt.title('Variazione del Tempo di Esecuzione con l\'Aumentare delle Dimensioni (per numero di processi)')
    plt.legend()
    plt.show()

    # Grafico con l'aumentare del numero di processi
    plt.figure(figsize=(10, 5))
    for dim_matrice in set([config[1] for config in dizionario.keys()]):
        tempi_dimensioni = [tempi_esecuzione[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][1] == dim_matrice]
        processi_dimensioni = [config[0] for config in dizionario.keys() if config[1] == dim_matrice]
        plt.plot(processi_dimensioni, tempi_dimensioni, label=f'Dimensione {dim_matrice}')

    plt.xlabel('Numero di Processi')
    plt.ylabel('Tempo di Esecuzione')
    plt.title('Variazione del Tempo di Esecuzione con l\'Aumentare del Numero di Processi (per dimensione)')
    plt.legend()
    plt.show()

# Esempio di utilizzo con un dizionario di prova
dizionario_prova = {
    (1, 100): 10,
    (1, 200): 15,
    (2, 100): 8,
    (2, 200): 12,
    (2, 300): 18,
    (4, 100): 5,
    (4, 200): 14
}
print(dictionaryMediaValore)
disegna_grafici(dictionaryMediaValore)