# PROJET ALG - MAUFROY Benjamin & COUVRAT Maëlle


### Organisation du projet 
ALG/  
├── make_fof.py | Script pour construire un FOF (File of files) à partir d'un fichier de séquences  
├── index_colored_fm_index_V2.py  | Programme d'indexation : construit l'index et le sérialise  
├── model.py  | Bibliothèque : classes `FM_Index` (avec color array) et `BloomFilter`, lecture du FOF et des FASTA, fonction `bf` d'indexation    
└── query.py  | Programme de requête : interroge un index sérialisé   

----

# Indexation de collections de génomes : FM-index généralisé et filtres de Bloom

Ce projet vise à comparer deux approches pour indexer N > 1 génomes et déterminer, pour une séquence requête q, dans quel(s) génome(s) elle apparaît sans erreur :

- **Approche 1 (`generalized`)** : FM-index (BWT, Count, Rank) de la concaténation des N génomes, séparés par `$`, associé à un **color array**. Ce tableau est aligné sur le suffix array : au lieu de la position de chaque suffixe, il stocke l'identifiant (la couleur, un entier de 0 à N-1) du génome d'où il provient.
- **Approche 2 (`bf`)** : N FM-index (un par génome) et N filtres de Bloom contenant les k-mers de chaque génome. Une requête n'interroge que les FM-index dont le filtre contient le premier k-mer de q.

## Installation

Python 3 est requis. Bibliothèques externes :

```
pip install pydivsufsort mmh3 bitarray numpy
```

## Format du FOF (File of Files) - Création du FOF

Fichier texte contenant sur chaque ligne le chemin d'un fichier FASTA de génome. L'ordre des lignes définit l'identifiant de chaque génome (le premier est le génome 1 dans les résultats, la couleur 0 en interne).

```
/mon/chemin/genome1.fasta
/mon/chemin/genome2.fasta
/mon/autre/chemin/ecl.fa
```

Constituer un dossier avec des génomes choisis.  
Puis créer votre fof à l'aide du script make_fof.py joint au projet. La commande utilisée pour le projet :  

```
python make_fof.py ./dossier                   
```

## Programme d'indexation : `index_colored_fm_index.py`

```
python3 index_colored_fm_index.py -i fof_file -o serialized_file -m <generalized|bf> [-h]
```

**Entrées**

| Option | Description |
|---|---|
| `-i` | Fichier FOF contenant les N génomes |
| `-o` | Fichier de sortie contenant la structure d'indexation sérialisée (`pickle`). Il est créé s'il n'existe pas, ainsi que les dossiers parents. Sans cette option, `index.pkl` est créé dans le dossier du FOF |
| `-m` | Méthode d'indexation : `generalized` (approche 1) ou `bf` (approche 2) |
| `-h` | Affiche l'aide et ne fait rien |
| `--help_methods <generalized\|bf>` | Décrit une méthode d'indexation et quitte |

**Sorties**

- L'index sérialisé dans `serialized_file`.
- Deux lignes de performance dans la console :

```
OUT TIME_BUILD: temps pour créer l'index avant sérialisation
OUT TIME_SERIALISATION: temps pour sérialiser l'index
```

Les temps sont des temps CPU, en secondes.

**Contenu de l'index sérialisé**

- Méthode `generalized` : un dictionnaire de clés `FM_index` (`bwt`, `n`, `r`), `Color_array` (génome d'origine de chaque suffixe) et `paths` (chemins des génomes, l'indice dans la liste est la couleur).
- Méthode `bf` : une liste de N dictionnaires (`Filtre_Bloom`, `FM_index`), un par génome, dans l'ordre du FOF.

**Exemple**

```
python3 index_colored_fm_index.py -i genomes.fof -o sortie/index.pkl -m generalized
python3 index_colored_fm_index.py -i genomes.fof -o sortie/index_bf.pkl -m bf
```

## Programme de requête : `query.py`

```
python3 query.py -i serialized_file -q query.fa -o result.txt [-h]
```

**Entrées**

| Option | Description |
|---|---|
| `-i` | Fichier contenant la structure d'indexation sérialisée |
| `-q` | Fichier FASTA contenant une ou plusieurs séquences requêtes sur l'alphabet {A, C, G, T} |
| `-o` | Fichier de sortie des résultats de présence/absence entre chaque séquence requête et chaque génome |
| `-h` | Affiche l'aide et ne fait rien |

La méthode d'indexation n'est pas à préciser : elle est reconnue d'après le contenu de l'index sérialisé.

**Sortie**

Un fichier texte contenant :

1. un en-tête avec les indices des génomes, séparés par une tabulation ;
2. une ligne de séparation laissée vide ;
3. une ligne par séquence requête : son nom, puis une suite de 0 ou 1 séparés par une tabulation (la i-ème valeur vaut 1 si la requête est présente dans le génome Gi, 0 sinon).

L'ordre des colonnes correspond à l'ordre des fichiers du FOF.

```
i_gen   1   2   3
                       
query1  0   1   0
query2  1   1   0
```

Deux lignes de performance sont aussi affichées dans la console :

```
OUT TIME_DESERIALISATION: temps pour relire la structure de données sur disque
OUT TIME_QUERY: temps pour réaliser les requêtes
```

**Exemple**

```
python3 query.py -i sortie/index.pkl -q requetes.fa -o resultat.txt
```

## Paramètres conseillés

| Paramètre | Valeur | Justification |
|---|---|---|
| Taille des k-mers (`bf`) | 31 | Valeur fixée dans `model.bf` : assez longue pour être spécifique à un génome |
| Nombre de fonctions de hachage H (`bf`) | 7 | Valeur fixée dans `model.bf`, cohérente avec le taux de faux positifs visé |
| Taux de faux positifs toléré (`bf`) | 0,01 | Le filtre est dimensionné par `bf_size` pour atteindre ce taux |

Les requêtes `bf` doivent être au moins aussi longues que la taille des k-mers, puisque c'est le premier k-mer de la requête qui est testé dans le filtre.



_____________________


## Reproduire les résultats du rapport

Les commandes ci-dessous ce basent sur le jeux de données fournis sur Moodle dans le cadre de ce projet. À adapter aux jeux de données de chaque expérience.

**1. Construire les deux index et relever les temps**

```
python3 index_colored_fm_index.py -i genomes.fof -o sortie/index_generalized.pkl -m generalized
python3 index_colored_fm_index.py -i genomes.fof -o sortie/index_bf.pkl -m bf
```

**2. Comparer les temps de construction et de sérialisation (graphique)**

Version simple (modifier les trois variables en haut du fichier, puis lancer) :

```
python3 bench_simple.py
```

Version complète (plusieurs répétitions, plusieurs jeux de données, nombre de génomes croissant) :

```
python3 bench_index.py -f genomes.fof -r 5
python3 bench_index.py -f petit.fof moyen.fof grand.fof -r 5 --log
python3 bench_index.py -f genomes.fof -n 1 2 3 5 -r 5
```

Le graphique est enregistré dans `comparaison.png` et un tableau de résultats dans `comparaison.csv`.

**3. Comparer les tailles des structures**

La taille mémoire est estimée par la taille du fichier sérialisé :

```
ls -l sortie/index_generalized.pkl sortie/index_bf.pkl
```

**4. Comparer les résultats et les temps de requête**

```
python3 query.py -i sortie/index_generalized.pkl -q requetes.fa -o sortie/resultat_generalized.txt
python3 query.py -i sortie/index_bf.pkl -q requetes.fa -o sortie/resultat_bf.txt
diff sortie/resultat_generalized.txt sortie/resultat_bf.txt
```

Aucune différence attendue : les deux structures doivent renvoyer les mêmes résultats.


