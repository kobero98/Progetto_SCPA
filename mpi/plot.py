import json
import matplotlib.pyplot as plt

def genera_dati_grafico(dizionario):
    numeri_processi = []
    dimensioni_matrice = []
    tempi_esecuzione = []

    for config, tempo in dizionario.items():
        num_processi, dim_matrice = config
        numeri_processi.append(num_processi)
        dimensioni_matrice.append(dim_matrice)
        tempi_esecuzione.append(tempo)

    return numeri_processi, dimensioni_matrice, tempi_esecuzione



def disegna_grafici(dizionario):
    numeri_processi, dimensioni_matrice, tempi_esecuzione = genera_dati_grafico(dizionario)

    # Grafico con l'aumentare delle dimensioni
    plt.figure(figsize=(10, 5))
    for num_processi in numeri_processi:
        tempo = [tempi_esecuzione[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][0] == num_processi]
        dimensione = [dimensioni_matrice[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][0] == num_processi]
        plt.plot(dimensione, tempo, label=f'{num_processi} processi')

    plt.xlabel('Dimensione Matrice')
    plt.ylabel('Tempo di Esecuzione')
    plt.title('Variazione del Tempo di Esecuzione con l\'Aumentare delle Dimensioni (per numero di processi)')
    plt.legend()
    plt.show()

    # Grafico con l'aumentare del numero di processi
    plt.figure(figsize=(10, 5))
    for dim_matrice in dimensioni_matrice:
        tempo = [tempi_esecuzione[i] for i in range(len(dizionario)) if list(dizionario.keys())[i][1] == dim_matrice]
        processi = [config[0] for config in dizionario.keys() if config[1] == dim_matrice]
        plt.plot(processi, tempo, label=f'Dimensione {dim_matrice}')

    plt.xlabel('Numero di Processi')
    plt.ylabel('Tempo di Esecuzione')
    plt.title('Variazione del Tempo di Esecuzione con l\'Aumentare del Numero di Processi (per dimensione)')
    plt.legend()
    plt.show()



def main():
    # f=open("result/result-10-10")
    # j=json.load(f)
    # print(j)
    dictResultTime={}
    for dim in [500,750,1000,2000,5000,10000]:
        f=open("result/result-"+str(dim)+"-1")
        js=json.load(f)
        dictResultTime[dim]=js["tempo_senza_creazione"]
        print(js)

    #per ogni iterazione prendo il tempo massimo di esecuzione e inseguito faccio una media per le 20 iterazioni
    #questo lo faccio per ogni matrice e per ogni numero di processi eseguiti
    dictionaryMediaValore={}
    for p in [4,8,12,16,20]:
        for dim in [500,750,1000,2000,5000,10000]:
            f=open("result/result-"+str(dim)+"-"+str(p))
            js=json.load(f)
            somma=0
            for k in range(1,20):
                max= 0.0
                for v in js[k]["vettore"]:
                        if "tempo_senza_creazione" in v and max<v["tempo_senza_creazione"]:
                            max=v["tempo_senza_creazione"]
                somma=somma+max
            dictionaryMediaValore[(p,dim)]=somma/20
    print(dictionaryMediaValore)

    disegna_grafici(dictionaryMediaValore)



if __name__ == "__main__":
    main()
