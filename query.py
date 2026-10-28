#!/usr/bin/env python3
"""Requêtage d'un index de génomes à partir d'une ou plusieurs séquences.

Ce script désérialise (pickle) un index produit par le script d'indexation, lit un
fichier FASTA de séquences à requêter et indique, pour chaque séquence, dans quels
génomes elle est présente. Le résultat est écrit dans un fichier texte.

Deux types d'index sont reconnus automatiquement :

- index ``bf`` (liste de dictionnaires, un par génome) : pour chaque génome, le
  premier k-mer (k = 31) de la requête est d'abord testé dans le filtre de Bloom,
  puis la requête entière est recherchée dans le FM-index du génome ;
- index ``generalized`` (dictionnaire contenant ``FM_index``, ``Color_array`` et
  ``paths``) : la requête est recherchée une seule fois dans le FM-index des
  génomes concaténés, puis le color array donne les génomes concernés.

Utilisation :
    python3 query.py -i <index_serialise> -q <fichier_fasta> [-o <fichier_sortie>]

Arguments :
    -i, --index        Chemin vers le fichier d'index sérialisé (obligatoire).
    -q, --query        Chemin vers un fichier FASTA contenant une ou plusieurs
                       séquences à requêter, sur l'alphabet {A,C,G,T} (obligatoire).
    -o, --outer_file   Chemin vers le fichier de sortie. S'il n'existe pas, il est
                       créé (ainsi que les dossiers parents manquants). Si l'option
                       n'est pas renseignée, un fichier ``query_result.txt`` est
                       créé dans le dossier contenant le fichier de requêtes.

Entrée :
    Dans le fichier FASTA, chaque entrée est une requête, numérotée dans l'ordre du
    fichier (query1, query2, ...). Les séquences sont lues en majuscules. Une entrée
    sans séquence est conservée dans la numérotation et considérée comme absente de
    tous les génomes.

Sortie :
    Un fichier texte dont la première ligne contient ``i_gen`` puis le numéro de
    chaque génome (dans l'ordre du FoF), suivie d'une ligne vide, puis d'une ligne
    par requête avec 1 si la requête est présente dans le génome et 0 sinon.
    Les noms (``i_gen``, ``query1``, ...) sont complétés par des espaces jusqu'à une
    largeur fixe (multiple de 8 moins 1, au moins égale à la longueur du plus long
    nom) puis suivis de deux tabulations ; les cellules suivantes sont séparées par
    une tabulation. Ainsi la première cellule de chaque ligne démarre toujours à la
    même colonne, même avec plus de 99 requêtes.
    Les temps CPU de désérialisation et de requêtage sont affichés dans la console.

Dépendances :
    Le module local ``model``, qui doit fournir les fonctions ``query``,
    ``backward_search`` et ``backward_search_interval``.
"""

# import
import pickle
import time
import argparse
import model
import os
import pathlib

# Arguments
# Définition de l'interface en ligne de commande.
parser  = argparse.ArgumentParser(
    description = "Ce programme requête un index de génomes à partir d'une ou plusieurs séquences en entrée"
)
# Index sérialisé à interroger
parser.add_argument("-i", "--index", type = str, required = True, help = "Un fichier d'indexation sérialisé")
# Fichier FASTA des séquences à requêter
parser.add_argument("-q", "--query", type = str, required = True, help = "Un fichier fasta contenant une ou plusieurs séquences à requêter sur l'alphabet {A,C,G,T}")
# Fichier de sortie : la valeur par défaut "/query_result.txt" sert de marqueur pour
# créer le fichier dans le dossier du fichier de requêtes (voir la partie Main)
parser.add_argument("-o", "--outer_file", type = str, help = "Un chemin vers un fichier en sortie de type .txt", default = "/query_result.txt")
# Lecture et analyse des arguments passés en ligne de commande
args = parser.parse_args()

# Main
# Création d'un fichier de sortie
# Vérification de la présence de l'option -o. Si non, création d'un fichier par défaut dans le dossier contenant le fichier de requêtes.
if args.outer_file == "/query_result.txt":
    # Option -o non renseignée : on utilise le dossier contenant le fichier de requêtes
    default_path = pathlib.Path(args.query).parent
    # Chemin final : <dossier_du_fichier_de_requêtes>/query_result.txt
    storage_path = f"{default_path}{args.outer_file}"
else:
    # Option -o renseignée : vérification de l'existence du fichier de sortie
    if os.path.exists(args.outer_file):
        # Le fichier existe déjà : on le réutilise tel quel (son contenu sera écrasé)
        storage_path = args.outer_file
    else:
        # Le fichier n'existe pas : création des dossiers parents manquants avec os.makedirs.
        # On passe le dossier parent et non le chemin complet, sinon un dossier
        # portant le nom du fichier serait créé à la place du fichier.
        parent_dir = os.path.dirname(args.outer_file)
        # parent_dir est vide si seul un nom de fichier est donné (ex. : out.txt)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok = True)
        # Le fichier lui-même est créé à l'ouverture en écriture, à la fin du script
        storage_path = args.outer_file

# Désérialisation
# Chargement de l'index depuis le fichier pickle
start_deserialisation = time.process_time()
with open(args.index, "rb") as filin:
    index = pickle.load(filin)

# Durée CPU de la désérialisation
print(f"OUT TIME_DESERIALISATION: {time.process_time() - start_deserialisation}")

# Requêtage
# Ouverture du fichier de requêtage
seq_list = []   # Une séquence par entrée FASTA, dans l'ordre du fichier
seq = None      # Séquence en cours de lecture (None tant qu'aucune entrée n'a été rencontrée)
with open(args.query, "r") as filin:
    for line in filin:
        if line.startswith(">"):
            # Nouvelle entrée : on enregistre la précédente (si elle existe)
            if seq is not None:
                seq_list.append(seq)
            seq = ""
        else:
            # Ligne de séquence : espaces retirés et passage en majuscules ; les lignes vides sont ignorées
            residues = line.strip().upper()
            if residues:
                if seq is None:
                    # Séquence placée avant tout en-tête : elle forme la première entrée
                    seq = ""
                seq += residues
    # Dernière entrée du fichier
    if seq is not None:
        seq_list.append(seq)

# Requêtage de l'index

# Détection de la méthode d'indexation : pour generalized = dictionnaire
generalized_index = None
if isinstance(index, dict) and "Color_array" in index:
    generalized_index = index
    index = []  # la boucle des filtres de Bloom ci-dessous n'a alors rien à parcourir
    # Nombre de génomes indexés = nombre de colonnes du fichier de sortie
    nb_genomes = len(generalized_index["paths"])
else:
    nb_genomes = len(index)



kmer = 31  # Taille des k-mers du filtre de Bloom (identique à celle utilisée à l'indexation)
results = {} # Création d'un dictionnaire associant les séquences requêtes aux génomes où elles sont présentes

start_query = time.process_time()

# Requêtage pour index bf
for i, query in enumerate(seq_list):
    results[f"query{i+1}"] = []
    for j, genome in enumerate(index):
        # Par défaut, la requête est considérée absente du génome j
        results[f"query{i+1}"].append("0\t")
        # Pré-filtre : le premier k-mer de la requête doit être dans le filtre de Bloom.
        # Le filtre pouvant donner des faux positifs, la requête entière est ensuite
        # vérifiée dans le FM-index. Une requête vide est considérée comme absente.
        if query and model.query(genome["Filtre_Bloom"], query[:kmer]):
            if model.backward_search(genome["FM_index"], query):
                results[f"query{i+1}"][j] = "1\t"


# Requêtage pour index generalized
if generalized_index is not None:
    fm_index = generalized_index["FM_index"]
    color = generalized_index["Color_array"]

    for i, query in enumerate(seq_list):
        # Par défaut, la requête est considérée absente de tous les génomes
        results[f"query{i+1}"] = ["0\t"] * nb_genomes
        # Intervalle de lignes du FM-index où la requête est présente (None si absente).
        # Une requête vide est considérée comme absente.
        interval = model.backward_search_interval(fm_index, query) if query else None

        if interval is not None:
            start, stop = interval
            # Génomes d'origine des suffixes de l'intervalle, donnés par le color array
            for genome_index in set(int(c) for c in color[start:stop + 1] if c is not None):
                results[f"query{i+1}"][genome_index] = "1\t"


# Durée CPU du requêtage
print(f"OUT TIME_QUERY: {time.process_time() - start_query}")

# Création du fichier de sortie

# Nombre de colonnes de résultats (un par génome)
i_gen = nb_genomes

# Largeur de la colonne des noms : un multiple de 8 moins 1, au moins égale à la longueur du
# plus long nom (i_gen compris). Les noms sont complétés avec des espaces puis suivis de deux
# tabulations : la première tabulation atteint un multiple de 8 et la seconde fixe la colonne
# de départ des résultats, identique sur toutes les lignes (avec des tabulations seules, les
# lignes à partir de query100 étaient décalées).
longest_name = max([len("i_gen")] + [len(name) for name in results])
name_width = 8 * (longest_name // 8 + 1) - 1

with open(storage_path, "w") as querin:
    # En-tête : numéro de chaque génome
    querin.write("i_gen".ljust(name_width) + "\t\t")
    for i in range(i_gen):
        querin.write(f"{i+1}\t")
    querin.write("\n\n")
    # Une ligne par requête : 1 si présente dans le génome, 0 sinon
    for name, values in results.items():
        querin.write(name.ljust(name_width) + "\t\t")
        for result in values:
            querin.write(result)
        querin.write("\n")