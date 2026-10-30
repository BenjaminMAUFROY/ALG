"""Requêtage d'un index de génomes à partir d'une ou plusieurs séquences.

Ce script désérialise (pickle) un index produit par le script d'indexation,
lit un fichier FASTA de séquences à requêter et indique, pour chaque séquence,
dans quels génomes elle est présente. Le résultat est écrit dans un fichier
texte.

Deux types d'index sont reconnus :

- index ``bf`` (liste de dictionnaires, un par génome) : pour chaque génome,
  le premier kmer de la requête (de taille ``kmer``, par défaut 31) est
  d'abord testé dans le filtre de Bloom avec ``H``(par défaut 7) fonctions de
  hachage, puis la requête entière est recherchée dans le FM-index du génome ;
- index ``generalized`` (dictionnaire contenant ``FM_index``, ``Color_array``
  et ``paths``) : la requête est recherchée une seule fois dans le FM-index
  des génomes concaténés, puis le color array donne les génomes concernés.

Paramètres du filtre de Bloom (index ``bf``) :
    ``kmer`` (taille des kmers) et ``H`` (nombre de fonctions de hachage)
    sont définis dans la partie « Requêtage pour index bf » du script. Ils
    doivent être identiques à ``kmer_size`` et ``n_hashage`` du script
    d'indexation.

Utilisation :
    python3 query.py -i <index_serialise> -q <fichier_fasta>
        [-o <fichier_sortie>]

Arguments :
    -i, --index        Chemin vers le fichier d'index sérialisé (obligatoire).
    -q, --query        Chemin vers un fichier FASTA contenant une ou plusieurs
                       séquences à requêter, sur l'alphabet {A,C,G,T}
                       (obligatoire).
    -o, --outer_file   Chemin vers le fichier de sortie. S'il n'existe pas, il
                       est créé (ainsi que les dossiers parents manquants). Si
                       l'option n'est pas renseignée, un fichier
                       ``query_result.txt`` est créé dans le dossier contenant
                       le fichier de requêtes.

Entrée :
    Dans le fichier FASTA, chaque entrée est une requête, numérotée dans
    l'ordre du fichier (query1, query2, ...).

Sortie :
    Un fichier texte dont la première ligne contient ``i_gen`` puis le numéro
    de chaque génome (dans l'ordre du FoF), suivie d'une ligne vide, puis d'une
    ligne par requête avec 1 si la requête est présente dans le génome et 0
    sinon.
    Les temps CPU de désérialisation et de requêtage sont affichés dans la
    console.

Dépendances :
    Le module local ``model``, qui doit fournir les fonctions ``query``,
    ``backward_search`` et ``backward_search_interval``.
"""

# import
import argparse
import os
import pathlib
import pickle
import time

import model

# Arguments
# Définition de l'interface en ligne de commande.
parser = argparse.ArgumentParser(
    description=(
        "Ce programme requête un index de génomes à partir d'une ou "
        "plusieurs séquences en entrée"
    )
)
parser.add_argument(
    "-i", "--index",
    type=str,
    required=True,
    help="Un fichier d'indexation sérialisé"
)
parser.add_argument(
    "-q", "--query",
    type=str,
    required=True,
    help=("Un fichier fasta contenant une ou plusieurs séquences à requêter"
          "sur l'alphabet {A,C,G,T}")
)
parser.add_argument(
    "-o", "--outer_file",
    type=str,
    help="Un chemin vers un fichier en sortie de type .txt",
    default="/query_result.txt"
)
args = parser.parse_args()

# Création d'un fichier de sortie
# Vérification de la présence de l'option -o. Si non, création d'un fichier
# par défaut dans le dossier contenant le fichier de requêtes.
if args.outer_file == "/query_result.txt":
    default_path = pathlib.Path(args.query).parent
    storage_path = f"{default_path}{args.outer_file}"
else:
    # Option -o renseignée : vérification de l'existence du fichier de sortie
    if os.path.exists(args.outer_file):
        storage_path = args.outer_file
    else:
        parent_dir = os.path.dirname(args.outer_file)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        storage_path = args.outer_file

# Désérialisation
# Chargement de l'index depuis le fichier pickle
start_deserialisation = time.process_time()
with open(args.index, "rb") as filin:
    index = pickle.load(filin)

# Durée CPU de la désérialisation
print("OUT TIME_DESERIALISATION: "
      f"{time.process_time() - start_deserialisation}")

# Requêtage
# Ouverture du fichier de requêtage

seq_list = []

seq = None
with open(args.query, "r") as filin:
    for line in filin:
        if line.startswith(">"):
            # Nouvelle entrée : on enregistre la précédente (si elle existe)
            if seq is not None:
                seq_list.append(seq)
            seq = ""
        else:
            residues = line.strip().upper()
            if residues:
                if seq is None:
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
    
    index = []
    # Nombre de génomes indexés = nombre de colonnes du fichier de sortie
    nb_genomes = len(generalized_index["paths"])
else:
    nb_genomes = len(index)


# Requêtage pour index bf
# Taille des k-mers du filtre de Bloom (identique à celle utilisée à
# l'indexation)
kmer = 31
# Nombre de fonctions de hachage du filtre de Bloom (identique à celui
# utilisé à l'indexation)
H = 7

results = {}

start_query = time.process_time()

for i, query in enumerate(seq_list):
    results[f"query{i + 1}"] = []
    for j, genome in enumerate(index):
        results[f"query{i + 1}"].append("0\t")
        # Pré-filtre : le premier kmer de la requête doit être dans le
        # filtre de Bloom. Le filtre pouvant donner des faux positifs, la
        # requête entière est ensuite vérifiée dans le FM-index.
        if query and model.query(genome["Filtre_Bloom"], query[:kmer], h=H):
            if model.backward_search(genome["FM_index"], query):
                results[f"query{i + 1}"][j] = "1\t"


# Requêtage pour index generalized
if generalized_index is not None:
    fm_index = generalized_index["FM_index"]
    color = generalized_index["Color_array"]

    for i, query in enumerate(seq_list):
        # Par défaut, la requête est considérée absente de tous les génomes
        results[f"query{i + 1}"] = ["0\t"] * nb_genomes
        # Intervalle de lignes du FM-index où la requête est présente (None
        # si absente).
        interval = (model.backward_search_interval(fm_index, query)
                    if query else None)

        if interval is not None:
            start, stop = interval
            # Génomes d'origine des suffixes de l'intervalle, donnés par le
            # color array
            for genome_index in set(int(c) for c in color[start:stop + 1]
                                    if c is not None):
                results[f"query{i + 1}"][genome_index] = "1\t"


# Durée CPU du requêtage
print(f"OUT TIME_QUERY: {time.process_time() - start_query}")

# Création du fichier de sortie

i_gen = nb_genomes

longest_name = max([len("i_gen")] + [len(name) for name in results])
name_width = 8 * (longest_name // 8 + 1) - 1

with open(storage_path, "w") as querin:
    # En-tête : numéro de chaque génome
    querin.write("i_gen".ljust(name_width) + "\t\t")
    for i in range(i_gen):
        querin.write(f"{i + 1}\t")
    querin.write("\n\n")
    # Une ligne par requête : 1 si présente dans le génome, 0 sinon
    for name, values in results.items():
        querin.write(name.ljust(name_width) + "\t\t")
        for result in values:
            querin.write(result)
        querin.write("\n")
