# PROJET ALG - MAUFROY Benjamin & COUVRAT Maëlle

### Organisation du projet 
ALG/
├── index_colored_fm_index.py
├── query.py
└── tests.py

## MODE D'EMPLOI DU PROGRAMME

### Programme index colored fm index.py

python index_colored_fm_index.py -i fof_file -o serialized_file -m <generalized|bf>
[-h] [other params]

o Entrées :
∗ -i un fichier FOF2 (contenant dans notre cas N lignes, chacune d´ecrivant un g´enome Gi).
∗ -o le fichier de sortie contenant la structure d’indexation s´erialis´ee3
∗ -m la m´ethode d’indexation utilis´ee, choix entre generalized pour l’approche 1 (FM-index
g´en´eralis´e et color array) et bf pour l’approche 2 avec les filtres de bloom.
∗ -h affiche l’aide et ne fait rien.

O Sortie :
∗ la structure d’indexation s´erialis´ee, dans le fichier serialized file (voir Section 6.2),
∗ des informations de performances sorties dans la console comme indiqu´e ci-dessous.
OUT TIME_BUILD: temps pour cr´eer l'index avant s´erialisation
OUT TIME_SERIALISATION: temps pour s´erialiser l'index
D’autres informations pourront ˆetre ajout´ees tant que ces lignes apparaissent sous ce format.


### Programme query.py

python query.py -i serialized_file -q query.fa -o result.txt [-h]

o Entrées :
∗ -i le fichier contenant la structure d’indexation s´erialis´ee
∗ -q un fichier fasta contenant une ou plusieurs s´equences `a requˆeter : un ensemble de s´equences
génomiques sur l’alphabet Σ = {A, C, G, T }.
∗ -o le fichier de sortie contenant les r´esultats de pr´esence/absence entre chaque s´equence
requête et chaque g´enome Gi.
∗ -h affiche l’aide et ne fait rien.

o Sortie :
Un fichier de sortie au format txt, avec (1) un header pour les indexes des g´enomes, s´epar´es
par une tabulation, (2) une ligne de s´eparation, laiss´ee vide, (3) une ligne par s´equence
requˆete compos´ee du header de la s´equence et une suite de valeurs 0 ou 1, s´epar´ees par une
tabulation (la i`eme valeur vaut 1 si la s´equence requˆete est pr´esente dans le g´enome Gi, 0
sinon)


## Exemples pour reproduire le code du rapport