# PROJET ALG - MAUFROY Benjamin & COUVRAT Maëlle

### RESTE A FAIRE
- Dans model.py:
  - Commentaire
  - Docstring
  - Mise au propre
  - Ajouter la fonction d'indexation generalized
- Dans index_colored_fm_index.py
  - Compléter les commentaires/docstring
  - Mettre au propre
  - Ajouter la méthode d'indexation generalized
  - Tester si ça fonctionne
- Faire query.py
- Faire test.py
- Faire le rapport

### Organisation du projet 
ALG/  
├── make_fof.py  
├── index_colored_fm_index.py  
├── query.py  
└── tests.py  

## MODE D'EMPLOI DU PROGRAMME

### Programme index colored fm index.py

python index_colored_fm_index.py -i fof_file -o serialized_file -m <generalized|bf>
[-h] [other params]

o Entrées :
∗ -i un fichier FOF2 (contenant dans notre cas N lignes, chacune décrivant un génome Gi).  
∗ -o le fichier de sortie contenant la structure d’indexation s´erialis´ee3  
∗ -m la méthode d’indexation utilisée, choix entre generalized pour l’approche 1 
(FM-index généralisé et color array) et bf pour l’approche 2 avec les filtres de bloom.  
∗ -h affiche l’aide et ne fait rien.  

O Sortie :
∗ la structure d’indexation sérialisée, dans le fichier serialized file (voir Section 6.2),  
∗ des informations de performances sorties dans la console comme indiqué ci-dessous.  
OUT TIME_BUILD: temps pour créer l'index avant sérialisation  
OUT TIME_SERIALISATION: temps pour sérialiser l'index  
D’autres informations pourront ˆetre ajout´ees tant que ces lignes apparaissent sous ce format.  


### Programme query.py

python query.py -i serialized_file -q query.fa -o result.txt [-h]

o Entrées :
∗ -i le fichier contenant la structure d’indexation sérialisée  
∗ -q un fichier fasta contenant une ou plusieurs s´equences à requêter : un ensemble de séquences
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

### Commencer par la création du FOF (File of files)

Constituer un dossier "genomes" avec des génomes choisis.  
Par exemple avec les commandes suivantes qui vous permettent de récupérer des fragments de genomes humains :  

""
sudo install curl
curl -s "https://rest.ensembl.org/sequence/region/human/21:30000000..30001000:1?content-type=text/plain" -o fragment_chr21.txt
curl -s "https://rest.ensembl.org/sequence/region/human/1:1000000..1001000:1?content-type=text/x-fasta" -o chr1_fragment.fa
curl -s "https://rest.ensembl.org/sequence/region/human/7:117000000..117001000:1?content-type=text/x-fasta" -o chr7_fragment.fa
curl -s "https://rest.ensembl.org/sequence/region/human/X:50000000..50001000:1?content-type=text/x-fasta" -o chrX_fragment.fa
curl -s "https://rest.ensembl.org/sequence/region/human/MT:1..1000:1?content-type=text/x-fasta" -o chrMT_fragment.fa
""

Puis créer votre fof à l'aide du script make_fof.py joint au projet. La commande utilisée pour le projet :  

""
    python make_fof.py ./genomes                   
""
