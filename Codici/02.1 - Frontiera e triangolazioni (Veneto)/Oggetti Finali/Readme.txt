Ci sono tre tipi di definizione della frontiera. In tutti i casi sono stati eliminati
i triangoli totalmente con vertici di frontiera nell'isola di Venezia, in quanto coprono una
grossa area raggruppata, e creano una grossa distanza con il comune.

Tutti gli RData contengono un vettore Codici, che in modo corrispondente a x e y ha il codice di comune se
il punto è di comune, 0 se di frontiera.

1) ConTriangoli
Non è stato eliminato nessun triangolo con vertici totalmente di frontiera

1) ConMenoTriangoli
Sono stati eliminati solo alcuni triangoli che causano un bordo troppo marcato di triangoli da eliminare

3) SenzaTriangoli
Sono stati eliminati tutti i triangoli con vertici di frontiera. Quelli che, se eliminati, creerebbero
un dominio diviso in due parti, sono stati ovviamente tenuti