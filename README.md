# Pileup

## Autori

*856114 - Costantini Davide*  
*869877 - Giannaccari Mattia*

## Descrizione

Uno dei metodi di esplorazione dei file SAM/BAM consiste nell'iterare sul risultato della funzione pileup, ovvero una lista di oggetti (PileupColumn) che contengono le basi coperte dagli allineamenti. Il progetto consiste nel reimplementare questa funzione, aggiungendo un requisito sulla qualità delle singole basi della query.

## Requisiti per l'utilizzo

Per eseguire il progetto sono necessari:

* Python 3
* Conda
* Pysam

Le librerie necessarie sono incluse nei file environment.yml, che possono essere importati con il seguente comando:
`conda env create -f environment.yml`
