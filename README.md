# PROJET ALG - MAUFROY Benjamin & COUVRAT Maëlle


### Organisation du projet 
ALG/  
├── make_fof.py | Script pour construire un FOF (File of files) à partir d'un fichier de séquences  
├── index_colored_fm_index.py  | Programme d'indexation : construit l'index et le sérialise  
├── model.py  | Bibliothèque : classes `FM_Index` (avec color array) et `BloomFilter`, lecture du FOF et des FASTA, fonction `bf` d'indexation    
├── query.py  | Programme de requête : interroge un index sérialisé   
└── Rapport.pdf | Rapport d'analyses sur le jeux de données fournis permettant d'évaluer le programme  

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


# Reproduire les résultats du rapport

Les commandes ci-dessous ce basent sur le jeux de données fournis sur Moodle dans le cadre de ce projet. Les commandes sont à adapter aux jeux de données de chaque expérience.
Le script de benchmark automatisé utilisé pour lancer les mesures est optionnel et n'est pas livré avec le projet. Les mesures peuvent être refaites à la main avec les commandes ci-dessous.

#### Environnement

- Python 3
- Bibliothèques : pydivsufsort, mmh3, bitarray, numpy (voir « Installation »).

#### Données

Les données viennent de Moodle, dans le dossier `data_test`. Il contient 4 jeux de 49 génomes chacun :

| Jeu | Description |
|---|---|
| `coli_10k` | Séquences d'E.coli |
| `coli_100k` | Séquences d'E.coli |
| `random_10k` | Séquences random de ACGT |
| `random_100k` |  Séquences random de ACGT |

Chaque jeu est interrogé avec deux fichiers de requêtes :

- `queries_<jeu>_len100bp.fa` : 1 960 requêtes de 100 pb ;
- `queries_<jeu>_len1000bp.fa` : 196 requêtes de 1 000 pb.

Toutes les requêtes sont plus longues que la taille des k-mers (31), condition nécessaire pour la méthode `bf`.

#### Mesures relevées

Les temps sont des temps CPU (`time.process_time`), en secondes, affichés par les scripts :

| Programme | Valeurs affichées |
|---|---|
| `index_colored_fm_index.py` | `OUT TIME_BUILD`, `OUT TIME_SERIALISATION` |
| `query.py` | `OUT TIME_DESERIALISATION`, `OUT TIME_QUERY` |

Le poids d'un index est la taille du fichier `.pkl` sur disque (`ls -l`). L'index est construit une seule fois par jeu et par configuration : ses temps et son poids sont répétés pour les deux fichiers de requêtes.

## 1. Créer les FoF

Un FoF par jeu, avec `make_fof.py` :

```bash
python3 make_fof.py ./data_test/coli_10k
python3 make_fof.py ./data_test/coli_100k
python3 make_fof.py ./data_test/random_10k
python3 make_fof.py ./data_test/random_100k
```

Les FoF utilisés dans les fichiers de résultats sont `resultats/fof/<jeu>.fof`.

## 2. Expérience 1 : configuration par défaut (`temps_et_poids.txt`)

Paramètres par défaut : fp_rate = 1 % (0.01 dans le code), H = 7, kmer = 31. Pour chaque jeu `<jeu>` et chaque méthode :

#### Indexation (relever `TIME_BUILD` et `TIME_SERIALISATION`)

```bash
python3 index_colored_fm_index.py -i resultats/fof/<jeu>.fof -o resultats/index/<jeu>_generalized.pkl -m generalized
python3 index_colored_fm_index.py -i resultats/fof/<jeu>.fof -o resultats/index/<jeu>_bf.pkl -m bf
```

#### Poids des index

On récupère le poid des index avec l'exemple de commande suivante :

```bash
ls -l resultats/index/<jeu>_generalized.pkl resultats/index/<jeu>_bf.pkl
```

#### Requêtage (relever `TIME_DESERIALISATION` et `TIME_QUERY`), pour chacun des deux fichiers de requêtes

```bash
python3 query.py -i resultats/index/<jeu>_generalized.pkl -q data_test/<jeu>/queries_<jeu>_len100bp.fa  -o resultats/sortie/<jeu>_generalized_100.txt
python3 query.py -i resultats/index/<jeu>_generalized.pkl -q data_test/<jeu>/queries_<jeu>_len1000bp.fa -o resultats/sortie/<jeu>_generalized_1000.txt
python3 query.py -i resultats/index/<jeu>_bf.pkl          -q data_test/<jeu>/queries_<jeu>_len100bp.fa  -o resultats/sortie/<jeu>_bf_100.txt
python3 query.py -i resultats/index/<jeu>_bf.pkl          -q data_test/<jeu>/queries_<jeu>_len1000bp.fa -o resultats/sortie/<jeu>_bf_1000.txt
```

Remarques :

- Les chemins `resultats/index/`, `resultats/sortie/` et `data_test/<jeu>/` sont des exemples à adapter.
- L'option `-o` est déclarée `required=True` dans `index_colored_fm_index.py`.

## 3. Expérience 2 : variation des paramètres du filtre de Bloom (`temps_et_poids_variations.txt`)

Seule la méthode `bf` est concernée pour les variations de paramètres. Un seul paramètre varie à la fois, à partir de la configuration par défaut (fp_rate = 1 %, H = 7, kmer = 31) :

| Série | Valeurs | Valeurs fixées |
|---|---|---|
| 1 : fp_rate | 0,01 ; 0,1 ; 1 ; 5 ; 10 % | H = 7, kmer = 31 |
| 2 : H | 3 ; 5 ; 7 ; 9 ; 11 | fp_rate = 1 %, kmer = 31 |
| 3 : kmer | 11 ; 15 ; 20 ; 25 ; 31 | fp_rate = 1 %, H = 7 |

La configuration par défaut est commune aux trois séries : elle n'est mesurée qu'une fois, soit 13 configurations distinctes. La méthode `generalized` est mesurée comme référence ; elle ne dépend pas de ces paramètres.

**Où modifier les paramètres.** fp_rate est exprimé en fraction dans le code (1 % = 0.01). Les valeurs doivent être identiques à l'indexation et au requêtage :

| Paramètre | `index_colored_fm_index.py` | `query.py` |
|---|---|---|
| H | `n_hashage` | `H` |
| fp_rate | `fp_rate` | (non utilisé) |
| kmer | `kmer_size` | `kmer` |

Exemple pour la série H, jeu `coli_10k` (à répéter avec chaque valeur, après avoir modifié `n_hashage` dans `index_colored_fm_index.py` et `H` dans `query.py`) :

```bash
python3 index_colored_fm_index.py -i resultats/fof/coli_10k.fof -o resultats/index/coli_10k_bf_H3.pkl -m bf
python3 query.py -i resultats/index/coli_10k_bf_H3.pkl -q data_test/coli_10k/queries_coli_10k_len100bp.fa -o resultats/sortie/coli_10k_bf_H3_100.txt
```

**Répétitions et médianes.** Chaque mesure (indexation complète + requêtage) est répétée 3 fois, les répétitions étant entrelacées entre configurations :
Les valeurs du fichier sont les médianes des 3 répétitions ; les essais bruts figurent en fin de fichier.

**Bits par k-mer.** Valeur théorique du filtre, donnée par `bf_size` :

```
bits/k-mer = -H / ln(1 - fp_rate^(1/H))
```

## 4. Contrôles

- Poids d'index identique d'une répétition à l'autre.
- Résultats de `query.py` identiques à une référence : 336 / 336 requêtages comparés (toutes configurations, toutes répétitions).
- Variabilité entre répétitions, (max − min) / médiane, médiane sur l'ensemble des mesures : `TIME_BUILD` 8 %, `TIME_SERIALISATION` 12 %, `TIME_DESERIALISATION` 11 %, `TIME_QUERY` 15 %. Les écarts de temps inférieurs à cet ordre de grandeur ne sont donc pas significatifs.

Pour vérifier que les deux méthodes donnent les mêmes résultats sur un jeu :

```bash
diff resultats/sortie/<jeu>_generalized_100.txt resultats/sortie/<jeu>_bf_100.txt
```

> Aucune différence n'est attendue. Les deux méthodes doivent fournir les mêmes résultats.

#### Format des fichiers de résultats

- `temps_et_poids.txt` : un bloc par jeu, méthode et fichier de requêtes, puis un récapitulatif tabulé (jeu, methode, requetes, n_requetes, TIME_BUILD, TIME_SERIALISATION, poids_index_octets, TIME_DESERIALISATION, TIME_QUERY).
- `temps_et_poids_variations.txt` : protocole, contrôles, un bloc par jeu et fichier de requêtes (référence generalized et séries fp_rate, H, kmer), un récapitulatif des médianes tabulé, puis les essais bruts (une ligne par répétition).

#### Paramètres par défaut (méthode `bf`)

| Paramètre | Valeur | Justification |
|---|---|---|
| Taille des k-mers | 31 | Assez longue pour être spécifique à un génome |
| Nombre de fonctions de hachage H | 7 | Cohérent avec le taux de faux positifs visé |
| Taux de faux positifs | 0,01 | Le filtre est dimensionné par `bf_size` pour atteindre ce taux |

Ces valeurs se règlent en tête de `index_colored_fm_index.py` (`kmer_size`, `n_hashage`, `fp_rate`) et dans `query.py` (`kmer`, `H`), et non dans `model.bf` : `kmer` et `H` doivent rester identiques dans les deux scripts. Les requêtes doivent être au moins aussi longues que la taille des k-mers, puisque c'est le premier k-mer de la requête qui est testé dans le filtre.


